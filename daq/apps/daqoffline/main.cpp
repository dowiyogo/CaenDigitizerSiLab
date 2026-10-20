// daqoffline: el monitor, los bits de datos, los paneles y el reloj del TTT de daqrun sobre datos
// ya escritos (fase 5, DAQ2-9). Solo lectura de los datos; resultados en /data/daq/sim/offline/.
//   daqoffline --config <cfg.json> --out <dir> ARCHIVO_seg*.root...     (runs reales, en orden)
//   daqoffline --config <cfg.json> --out <dir> --fake <escenario.json> --duration S
// - Archivos reales: árbol events (adc, ttt_raw, t_read_utc_ns, t_utc_ns) y árbol hv (muestras de
//   hvmon guardadas por daqrun). Bloques = eventos con la misma hora de lectura. Los períodos del
//   monitor se cierran con la hora de lectura (reloj de los datos). El estado de los bits sigue de
//   un run al siguiente (runs 7 y 8 se analizan como una sola toma); cada run tiene su TttClock.
// - Fake: FakeDigitizer con reloj manual (acelerado) y muestras HV sintéticas (todo encendido salvo
//   los hv_group_off del escenario, con Pw = 0 en su canal).
// Salidas: units.csv (por período y unidad), episodes.csv, mails.jsonl (correos que se habrían
// mandado), ttt_clock.csv (t_read − t_utc por hora: TttClock, tick fijo, archivo), panels.root,
// gains.json, summary.json.
#include <sys/stat.h>

#include <TFile.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/core/logger.h"
#include "daq/core/time.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/data_flags.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/monitor.h"
#include "daq/dig/fake_digitizer.h"
#include "daq/hvmon/mapping.h"
#include "daq/io/panels_root.h"

using namespace daq;
using namespace daq::daqrun;

namespace {

constexpr const char* kOutRoot = "/data/daq/sim/offline/";

struct HvRow {
  UtcNs t;
  HvSample s;
};

double quantile(std::vector<double> v, double q) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  return v[static_cast<std::size_t>(q * static_cast<double>(v.size() - 1))];
}
double median(std::vector<double> v) { return quantile(std::move(v), 0.5); }

struct Episode {
  std::string key, kind;
  std::uint64_t episode;
  UtcNs start = 0, end = 0;
  double maxRatio = 0, minP = 1;
};

class Offline {
 public:
  Offline(const DaqrunConfig& cfg, const hvmon::Mapping& m, const std::string& out, UtcNs origin)
      : cfg_(cfg), map_(m), out_(out), monitor_(cfg, origin, {}, &m), engine_(cfg.dataFlags, m, "daqoffline") {
    units_.open(out + "/units.csv");
    units_ << "t_utc_ns,t_utc,key,kind,rate,baseline,mu,n,ratio,p,decided,active,live\n";
    mails_.open(out + "/mails.jsonl");
  }

  void block(const dig::EventBlock& b, const HvSample* hv) {
    monitor_.process(b);
    events_ += b.events.size();
    const auto& bits = engine_.bits();
    for (std::size_t i = 0; i < b.events.size(); ++i) {
      bool b6 = false, b7 = false, b11 = false;
      for (auto x : bits.ch) {
        b6 = b6 || (x & chflag::DataBlackout);
        b7 = b7 || (x & chflag::LightLeak);
        b11 = b11 || (x & chflag::StatInsufficient);
      }
      evBit6_ += b6;
      evBit7_ += b7;
      evBit11_ += b11;
      evProvisional_ += bits.provisional;
    }
    MonitorSummary s;
    if (monitor_.takeSummary(b.tReadUtcNs, s, 0)) period(s, hv);
  }

  void finish(const nlohmann::json& extra) {
    for (auto& [k, e] : open_) {
      e.end = lastT_;
      episodes_.push_back(e);
    }
    std::ofstream ep(out_ + "/episodes.csv");
    ep << "key,kind,episode,start_utc,end_utc,duration_min,max_ratio,min_p\n";
    std::map<std::string, int> count;
    for (const auto& e : episodes_) {
      ep << e.key << "," << e.kind << "," << e.episode << "," << formatUtcIso(e.start) << "," << formatUtcIso(e.end) << ","
         << static_cast<double>(e.end - e.start) / 6e10 << "," << e.maxRatio << "," << e.minP << "\n";
      ++count[e.kind];
    }
    monitor_.publishView();
    const auto pv = monitor_.panels().first;
    io::PanelRoot pr;
    pr.update(pv, engine_.view(), cfg_.dataFlags.lightLeakRatio);
    std::unique_ptr<TFile> f(TFile::Open((out_ + "/panels.root").c_str(), "RECREATE"));
    pr.write(*f);
    f->Close();
    nlohmann::json gains = nlohmann::json::array();
    std::map<std::string, int> gainStatus;
    for (const auto& g : pv.gains) {
      gains.push_back({{"hour_start_utc", formatUtcIso(g.hourStartUtcNs)}, {"ch", g.ch}, {"value", g.value}, {"status", g.status}, {"peaks", g.peaks}});
      if (g.status != "dark_rate") ++gainStatus[g.status];
    }
    std::ofstream(out_ + "/gains.json") << gains.dump(1) << "\n";
    nlohmann::json sum = extra;
    sum["events"] = events_;
    sum["periods"] = periods_;
    sum["episodes_by_kind"] = count;
    sum["mails_would_send"] = mailCount_;
    sum["events_with_bit6"] = evBit6_;
    sum["events_with_bit7"] = evBit7_;
    sum["events_with_bit11"] = evBit11_;
    sum["events_provisional"] = evProvisional_;
    sum["gain_status_count"] = gainStatus;
    sum["data_flags"] = cfg_.dataFlags.toJson();
    sum["panels"] = cfg_.panels.toJson();
    std::ofstream(out_ + "/summary.json") << sum.dump(1) << "\n";
    std::printf("eventos %llu, períodos %llu, episodios:", static_cast<unsigned long long>(events_), static_cast<unsigned long long>(periods_));
    for (const auto& [k, n] : count) std::printf(" %s %d", k.c_str(), n);
    std::printf("; correos que se habrían mandado %d; ajustes de ganancia:", mailCount_);
    for (const auto& [k, n] : gainStatus) std::printf(" %s %d", k.c_str(), n);
    std::printf("\n");
    for (const auto& e : episodes_)
      std::printf("  episodio %-22s %-14s #%llu  %s -> %s  (%.0f min)  razón máx %.2f  p mín %.2g\n", e.key.c_str(), e.kind.c_str(),
                  static_cast<unsigned long long>(e.episode), formatUtcIso(e.start).c_str(), formatUtcIso(e.end).c_str(),
                  static_cast<double>(e.end - e.start) / 6e10, e.maxRatio, e.minP);
  }

 private:
  void period(const MonitorSummary& s, const HvSample* hv) {
    ++periods_;
    lastT_ = s.tEndUtcNs;
    for (auto& e : engine_.addPeriod(s, hv)) {
      mails_ << nlohmann::json{{"t_utc", formatUtcIso(e.tUtcNs)}, {"kind", e.kind}, {"ch", e.ch}, {"episode", e.episode},
                               {"reason", e.reason}, {"context", e.context}}.dump()
             << "\n";
      ++mailCount_;
    }
    for (const auto& u : engine_.view().units) {
      units_ << s.tEndUtcNs << "," << formatUtcIso(s.tEndUtcNs) << "," << u.key << "," << u.kind << "," << u.rateHz << ","
             << u.baselineHz << "," << u.mu << "," << u.n << "," << u.ratio << "," << u.p << "," << u.decided << "," << u.active << ","
             << u.live << "\n";
      const std::string kind = u.kind == "face" ? "bit7_light_leak" : u.kind == "group" ? "bit6_blackout" : "imon_warning";
      auto it = open_.find(u.key);
      if (u.active) {
        if (it == open_.end()) it = open_.emplace(u.key, Episode{u.key, kind, u.episode, u.since, 0, 0, 1}).first;
        it->second.maxRatio = std::max(it->second.maxRatio, u.ratio);
        it->second.minP = std::min(it->second.minP, u.p);
      } else if (it != open_.end()) {
        it->second.end = s.tEndUtcNs;
        episodes_.push_back(it->second);
        open_.erase(it);
      }
    }
  }

  const DaqrunConfig& cfg_;
  const hvmon::Mapping& map_;
  std::string out_;
  Monitor monitor_;
  DataFlagEngine engine_;
  std::ofstream units_, mails_;
  std::map<std::string, Episode> open_;
  std::vector<Episode> episodes_;
  std::uint64_t events_ = 0, periods_ = 0, evBit6_ = 0, evBit7_ = 0, evBit11_ = 0, evProvisional_ = 0;
  int mailCount_ = 0;
  UtcNs lastT_ = 0;
};

std::vector<HvRow> readHv(TFile& f) {
  std::vector<HvRow> out;
  auto* t = dynamic_cast<TTree*>(f.Get("hv"));
  if (!t) return out;
  Long64_t rec = 0, treq = 0, tresp = 0;
  ULong64_t rseq = 0;
  Bool_t alive = false, stale = true, conn = false, ok = false;
  Double_t v0[12], i0[12], vm[12], im[12];
  UInt_t st[12], pw[12];
  t->SetBranchAddress("received_utc_ns", &rec);
  t->SetBranchAddress("read_seq", &rseq);
  t->SetBranchAddress("t_req_utc_ns", &treq);
  t->SetBranchAddress("t_resp_utc_ns", &tresp);
  t->SetBranchAddress("alive", &alive);
  t->SetBranchAddress("stale", &stale);
  t->SetBranchAddress("connected", &conn);
  t->SetBranchAddress("file_ok", &ok);
  t->SetBranchAddress("v0set_V", v0);
  t->SetBranchAddress("i0set_A", i0);
  t->SetBranchAddress("vmon_V", vm);
  t->SetBranchAddress("imon_A", im);
  t->SetBranchAddress("status", st);
  t->SetBranchAddress("pw", pw);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    HvSample s;
    s.receivedUtcNs = rec;
    s.readSeq = rseq;
    s.tReqUtcNs = treq;
    s.tRespUtcNs = tresp;
    s.alive = alive;
    s.stale = stale;
    s.connected = conn;
    s.fileOk = ok;
    s.stateUtcNs = rec;
    for (int h = 0; h < 12; ++h) {
      auto& c = s.ch[static_cast<std::size_t>(h)];
      c.present = true;
      c.v0setV = v0[h];
      c.i0setA = i0[h];
      c.vmonV = vm[h];
      c.imonA = im[h];
      c.status = st[h];
      c.pw = pw[h];
    }
    out.push_back({rec, s});
  }
  t->ResetBranchAddresses();
  return out;
}

int usage() {
  std::fprintf(stderr,
               "Uso: daqoffline --config cfg.json --out DIR ARCHIVO_seg*.root...\n"
               "     daqoffline --config cfg.json --out DIR --fake escenario.json --duration S\n"
               "DIR debe estar bajo %s\n",
               kOutRoot);
  return 2;
}

}  // namespace

int main(int argc, char** argv) {
  std::string cfgPath, out, fake;
  double duration = 0;
  std::vector<std::string> files;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--config" && i + 1 < argc) cfgPath = argv[++i];
    else if (a == "--out" && i + 1 < argc) out = argv[++i];
    else if (a == "--fake" && i + 1 < argc) fake = argv[++i];
    else if (a == "--duration" && i + 1 < argc) duration = std::stod(argv[++i]);
    else if (!a.empty() && a[0] != '-') files.push_back(a);
    else return usage();
  }
  if (cfgPath.empty() || out.empty() || (files.empty() && fake.empty())) return usage();
  if (out.rfind(kOutRoot, 0) != 0 || out.find("..") != std::string::npos) {
    std::fprintf(stderr, "daqoffline: la salida debe estar bajo %s\n", kOutRoot);
    return 2;
  }
  ::mkdir("/data/daq/sim/offline", 0775);
  ::mkdir(out.c_str(), 0775);
  std::ifstream ci(cfgPath);
  nlohmann::json cj = nlohmann::json::parse(ci);
  // Lo que importa del archivo de configuración son los bits, los paneles, la ventana y el mapeo;
  // la campaña, las rutas y el backend no se usan (no se escribe ningún run).
  cj["digitizer"]["backend"] = "fake";
  cj["campaign"] = {{"name", "offline"}};
  auto cfgR = DaqrunConfig::fromJson(cj);
  if (!cfgR) {
    std::fprintf(stderr, "daqoffline: %s\n", cfgR.error().message.c_str());
    return 2;
  }
  const DaqrunConfig cfg = cfgR.value();
  auto mapR = hvmon::Mapping::fromFile(cfg.mappingFile);
  if (!mapR) {
    std::fprintf(stderr, "daqoffline: mapeo: %s\n", mapR.error().message.c_str());
    return 2;
  }
  const hvmon::Mapping& mapping = mapR.value();
  nlohmann::json extra{{"config", cfgPath}, {"mapping", cfg.mappingFile}};

  if (!fake.empty()) {
    std::ifstream si(fake);
    nlohmann::json sj = nlohmann::json::parse(si);
    auto sp = dig::SyntheticParams::fromJson(sj);
    if (!sp) {
      std::fprintf(stderr, "daqoffline: escenario: %s\n", sp.error().message.c_str());
      return 2;
    }
    ManualClock clock(1792500000LL * 1000000000LL);
    auto info = dig::defaultChannelInfo();
    for (const auto& e : mapping.channels) {
      if (e.digCh < 0 || e.digCh >= dig::kChannels) continue;
      auto& c = info[static_cast<std::size_t>(e.digCh)];
      c.station = e.station == "hod2019" ? 1 : 0;
      c.plane = e.plane == "Y" ? 'Y' : 'X';
      if (e.strip >= 0) c.strip = e.strip;
      c.hvGroup = e.hvGroup ? *e.hvGroup : 0;
      c.face = e.face;
    }
    dig::FakeDigitizer fd(sp.value(), clock, info);
    if (auto s = fd.configure(cfg.digitizer); !s) return 2;
    auto origin = fd.start();
    if (!origin) return 2;
    const UtcNs t0 = origin.value().utcBeforeStart;
    Offline off(cfg, mapping, out, t0);
    // HV sintética: todo encendido salvo los hv_group_off del escenario.
    std::vector<std::tuple<int, UtcNs, UtcNs>> offs;
    for (const auto& e : sj.value("events", nlohmann::json::array()))
      if (e.value("type", "") == "hv_group_off")
        offs.emplace_back(e.value("group", 0) - 1, t0 + static_cast<UtcNs>(e["t_s"].get<double>() * 1e9),
                          t0 + static_cast<UtcNs>((e["t_s"].get<double>() + e.value("duration_s", 300.0)) * 1e9));
    for (double t = 0; t < duration; t += 1.0) {
      clock.advance(std::chrono::seconds(1));
      dig::EventBlock b;
      if (!fd.readBlock(b, 100000)) return 3;
      HvSample hv;
      hv.alive = hv.connected = hv.fileOk = true;
      hv.stale = false;
      hv.receivedUtcNs = hv.stateUtcNs = clock.now_utc_ns();
      for (int h = 0; h < 12; ++h) {
        auto& c = hv.ch[static_cast<std::size_t>(h)];
        c.present = true;
        c.pw = h < 8 ? 1 : 0;
        c.i0setA = 15e-6;
        c.imonA = 5e-6;
        for (const auto& [g, a, z] : offs)
          if (g == h && hv.receivedUtcNs >= a && hv.receivedUtcNs < z) {
            c.pw = 0;
            c.imonA = 0;
          }
      }
      if (b.events.empty()) {
        b.tReadUtcNs = clock.now_utc_ns();
        b.recordLength = cfg.digitizer.recordLength;
      }
      off.block(b, &hv);
    }
    extra["fake_scenario"] = fake;
    extra["duration_s"] = duration;
    off.finish(extra);
    return 0;
  }

  // Archivos reales (solo lectura).
  std::sort(files.begin(), files.end());
  std::unique_ptr<Offline> off;
  std::ofstream tc(out + "/ttt_clock.csv");
  tc << "run,bin_start_utc,events,median_read_minus_utc_fit_ms,p05_fit_ms,median_read_minus_utc_fixed_ms,p05_fixed_ms,"
        "median_read_minus_utc_file_ms\n";
  int curRun = -1;
  std::unique_ptr<TttClock> fit, fixed;
  std::map<std::pair<int, UtcNs>, std::array<std::vector<double>, 3>> diffs;
  nlohmann::json filesJ = nlohmann::json::array();
  for (const auto& path : files) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::fprintf(stderr, "daqoffline: no se pudo abrir %s\n", path.c_str());
      return 3;
    }
    auto* ri = f->Get<TNamed>("run_info_json");
    const auto info = nlohmann::json::parse(ri ? ri->GetTitle() : "{}");
    const int run = info.value("run", 0);
    const UtcNs origin = info["ttt_origin"].value("utc_before_start_ns", UtcNs{0});
    if (!off) off = std::make_unique<Offline>(cfg, mapping, out, origin);
    if (run != curRun) {
      TttClock::Params p;
      p.initialTickNs = info.value("ttt_tick_ns", 8.0);
      p.wrapBits = info.value("ttt_wrap_bits", 31);
      p.windowS = cfg.tttClockWindowS;  // digitizer.ttt_clock de la configuración
      p.reanchorS = cfg.tttClockReanchorS;
      p.minSpanS = cfg.tttClockMinSpanS;
      p.minPoints = static_cast<std::size_t>(std::max(2, cfg.tttClockMinPoints));
      p.firstFitS = cfg.tttClockFirstFitS;
      fit = std::make_unique<TttClock>(origin, p);
      p.fit = false;
      p.initialTickNs = 8.0;
      fixed = std::make_unique<TttClock>(origin, p);
      curRun = run;
    } else {
      fit->requestReanchor("segmento");
    }
    const auto hvRows = readHv(*f);
    auto* ev = dynamic_cast<TTree*>(f->Get("events"));
    const int n = static_cast<int>(ev->GetLeaf("adc")->GetLenStatic() / 64);
    std::vector<UShort_t> adc(static_cast<std::size_t>(64 * n));
    UInt_t raw = 0, cnt = 0;
    Long64_t tread = 0, tutc = 0;
    UChar_t kind = 0;
    ev->SetBranchStatus("*", 0);
    for (const char* b : {"adc", "ttt_raw", "evt_counter_hw", "t_read_utc_ns", "t_utc_ns", "trigger_kind"}) ev->SetBranchStatus(b, 1);
    ev->SetBranchAddress("adc", adc.data());
    ev->SetBranchAddress("ttt_raw", &raw);
    ev->SetBranchAddress("evt_counter_hw", &cnt);
    ev->SetBranchAddress("t_read_utc_ns", &tread);
    ev->SetBranchAddress("t_utc_ns", &tutc);
    ev->SetBranchAddress("trigger_kind", &kind);
    std::size_t hvIdx = 0;
    dig::EventBlock blk;
    blk.recordLength = n;
    auto flush = [&] {
      if (blk.events.empty()) return;
      while (hvIdx + 1 < hvRows.size() && hvRows[hvIdx + 1].t <= blk.tReadUtcNs) ++hvIdx;
      const HvSample* hv = hvRows.empty() || hvRows[hvIdx].t > blk.tReadUtcNs ? nullptr : &hvRows[hvIdx].s;
      off->block(blk, hv);
      blk.events.clear();
    };
    const Long64_t ne = ev->GetEntries();
    for (Long64_t i = 0; i < ne; ++i) {
      ev->GetEntry(i);
      if (!blk.events.empty() && tread != blk.tReadUtcNs) flush();
      blk.tReadUtcNs = tread;
      dig::Event e;
      e.counterHw = cnt;
      e.tttRaw = raw;
      e.triggerKind = kind;
      e.adc.assign(adc.begin(), adc.end());
      blk.events.push_back(std::move(e));
      const auto a = fit->map(raw, tread);
      const auto b = fixed->map(raw, tread);
      const UtcNs hour = tread - (tread % 600000000000LL);  // bins de 10 min
      auto& d = diffs[{run, hour}];
      d[0].push_back(static_cast<double>(tread - a.second) * 1e-6);
      d[1].push_back(static_cast<double>(tread - b.second) * 1e-6);
      d[2].push_back(static_cast<double>(tread - tutc) * 1e-6);
    }
    flush();
    filesJ.push_back({{"file", path}, {"run", run}, {"events", ne}});
    std::printf("%s: run %d, %lld eventos, %zu muestras HV\n", path.c_str(), run, ne, hvRows.size());
  }
  std::printf("t_read − t_utc por 10 min (mediana y percentil 5, ms): TttClock | tick fijo 8 ns | archivo\n");
  for (const auto& [k, d] : diffs) {
    tc << k.first << "," << formatUtcIso(k.second) << "," << d[0].size() << "," << median(d[0]) << "," << quantile(d[0], 0.05) << ","
       << median(d[1]) << "," << quantile(d[1], 0.05) << "," << median(d[2]) << "\n";
    std::printf("  run %d %s  %5zu ev  %8.3f (p05 %6.3f) | %8.3f (p05 %8.3f) | %8.3f\n", k.first, formatUtcIso(k.second).substr(0, 16).c_str(),
                d[0].size(), median(d[0]), quantile(d[0], 0.05), median(d[1]), quantile(d[1], 0.05), median(d[2]));
  }
  if (fit) {
    std::int64_t maxJump = 0;
    for (std::size_t i = 2; i < fit->anchors().size(); ++i) maxJump = std::max<std::int64_t>(maxJump, static_cast<std::int64_t>(std::llabs(fit->anchors()[i].jumpNs)));
    std::printf("tick ajustado al final del último run: %.6f ns, %zu anclas, salto máximo después del primer ajuste %.3f ms\n",
                fit->tickNs(), fit->anchors().size(), static_cast<double>(maxJump) * 1e-6);
  }
  extra["files"] = filesJ;
  off->finish(extra);
  return 0;
}
