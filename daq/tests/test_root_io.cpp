// Tests que usan ROOT: escritor de runs (segmentos, índice, SHA256, permisos, notas, flags con la
// línea de retardo), replay de archivos con el formato de v8 (archivo pequeño generado aquí, nunca
// /data/daq/raw), hodo::resolveRawFile sobre una campaña simulada y lector de control lento.
#include <gtest/gtest.h>

#include <sys/stat.h>
#include <unistd.h>

#include <TFile.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

#include "daq/core/file_system.h"
#include "daq/core/gzip.h"
#include "daq/core/sha256.h"
#include "daq/daqrun/campaign.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/hv_client.h"
#include "daq/dig/fake_digitizer.h"
#include "daq/io/hist_util.h"
#include <TH1D.h>
#include <TH2F.h>
#include "daq/dig/replay_digitizer.h"
#include "daq/io/run_writer.h"
#include "daq/io/sc_reader.h"
#include "hodo_paths.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::daqrun;

// LeakSanitizer (build ASan): cling, el intérprete de ROOT, deja memoria sin liberar al salir (227
// asignaciones también en un programa mínimo con un TFile y un TTree, DAQ2-5). Se suprime solo
// libCling.so y solo en este binario; nuestro código no pasa por ahí.
extern "C" const char* __lsan_default_suppressions() { return "leak:libCling.so\n"; }

namespace {

constexpr UtcNs kS = 1000000000LL;

std::string tmpDir(const std::string& tag) {
  const std::string d = ::testing::TempDir() + "daqroot_" + tag + "_" + std::to_string(::getpid());
  PosixFileSystem().makeDirectories(d, 0755);
  return d;
}

// Estado de hvmon al instante t: CH00-07 encendidos; CH02 con trip desde tripAt.
nlohmann::json hvState(UtcNs t, std::uint64_t k, UtcNs tripAt, nlohmann::json changes = nlohmann::json::array()) {
  nlohmann::json chs = nlohmann::json::array();
  for (int h = 0; h < 12; ++h) {
    const bool on = h < 8 && !(h == 2 && t >= tripAt);
    chs.push_back({{"ch", h},
                   {"values",
                    {{"Pw", on ? 1 : 0},
                     {"Status", on ? 1 : (h == 2 ? 512 : 0)},
                     {"V0Set_V", h < 8 ? 65.75 : 0.0},
                     {"VMon_V", on ? 65.75 : 0.2},
                     {"IMon_A", on ? 2e-6 : -4e-8},
                     {"I0Set_A", 1.5e-5}}}});
  }
  return {{"schema", "daqv2.hvmon_state/1"},
          {"seq", k},
          {"t_utc_ns", t},
          {"pid", 1},
          {"final", false},
          {"connection", {{"connected", true}, {"stale", false}}},
          {"last_read", {{"read_seq", k}, {"t_req_utc_ns", t - 10000000}, {"t_resp_utc_ns", t}}},
          {"channels", chs},
          {"changes", changes}};
}

std::map<std::string, std::string> readSums(const std::string& dir) {
  std::map<std::string, std::string> m;
  std::ifstream in(dir + "/SHA256SUMS");
  std::string sha, file;
  while (in >> sha >> file) m[file] = sha;
  return m;
}

mode_t modeOf(const std::string& p) {
  struct stat st {};
  ::stat(p.c_str(), &st);
  return st.st_mode & 0777;
}

}  // namespace

TEST(RunWriter, SegmentsIndexShaPermissionsNotesAndDelayedFlags) {
  const std::string base = tmpDir("writer");
  const std::string notesDir = base + "/sc";
  PosixFileSystem().makeDirectories(notesDir + "/2026", 0755);
  ManualClock clock;  // 2026-10-18T00:00:00Z
  const UtcNs t0 = clock.now_utc_ns();
  {
    std::ofstream o(notesDir + "/2026/notes_20261018.jsonl");
    o << nlohmann::json{{"t_utc_ns", t0 + 150 * kS}, {"user", "rene"}, {"host", "h"}, {"text", "nota en el seg 2"}}.dump() << "\n";
  }
  auto cfgR = DaqrunConfig::fromJson({{"campaign", {{"name", "wtest"}, {"raw_root", base + "/raw"}, {"manifests_dir", base + "/man"}}},
                                      {"run", {{"segment_max_s", 120}, {"autosave_s", 30}, {"flag_delay_s", 5.0}}},
                                      {"notes", {{"slowcontrol_dir", notesDir}}},
                                      {"mapping", {{"file", DAQ_SOURCE_DIR "/config/mapping/hodo_mapping_2026_10.json"}}}});
  ASSERT_TRUE(cfgR.ok()) << cfgR.error().message;
  const DaqrunConfig cfg = cfgR.value();
  const CampaignLocation loc{cfg.rawRoot, cfg.manifestsDir, cfg.campaign};
  ASSERT_TRUE(openCampaign(loc, "test", true, t0).ok());
  auto run = RunReservation::reserve(loc, t0);
  ASSERT_TRUE(run.ok());
  auto m = hvmon::Mapping::fromFile(cfg.mappingFile);
  ASSERT_TRUE(m.ok());
  auto sp = dig::SyntheticParams::fromJson({{"seed", 5}, {"trigger_rate_Hz", 2.0}});
  dig::FakeDigitizer fd(sp.value(), clock);
  ASSERT_TRUE(fd.configure(cfg.digitizer).ok());
  auto origin = fd.start();
  io::RunStatic rs;
  rs.configJson = cfg.toJson();
  rs.configSha256 = sha256Hex(rs.configJson.dump());
  rs.mappingText = "{}";
  rs.digitizerInfo = fd.info();
  rs.origin = origin.value();
  rs.hvAtStart = hvState(t0, 0, t0 + 200 * kS);
  io::WriterStats ws;
  io::RunWriter w(cfg, run.value(), m.value(), rs, clock, ws);
  ASSERT_TRUE(w.open().ok());
  std::uint64_t delivered = 0, maxPendingSeen = 0;
  const UtcNs tripAt = t0 + 200 * kS;
  for (std::uint64_t k = 1; k <= 400; ++k) {
    clock.advance(1s);
    const UtcNs now = clock.now_utc_ns();
    nlohmann::json ch = nlohmann::json::array();
    if (now == tripAt)
      ch.push_back({{"id", 1}, {"t_utc_ns", now}, {"t_prev_read_utc_ns", now - kS}, {"ch", 2}, {"param", "Pw"}, {"old", 1}, {"new", 0}, {"origin", "trip"}});
    w.addHv(parseHvState(hvState(now, k, tripAt, ch), now + 1000000, 3.0, true));
    auto b = std::make_shared<dig::EventBlock>();
    ASSERT_TRUE(fd.readBlock(*b, 1000).ok());
    delivered += b->events.size();
    ASSERT_TRUE(w.addBlock(b).ok());
    ASSERT_TRUE(w.tick().ok());
    maxPendingSeen = std::max<std::uint64_t>(maxPendingSeen, ws.pendingDelay.load());
  }
  // La línea de retardo retiene ~5 s de eventos (2 Hz) antes de escribirlos.
  EXPECT_GE(maxPendingSeen, 5u);
  EXPECT_LT(ws.eventsWritten.load(), delivered);
  ASSERT_TRUE(w.finish(true).ok());
  EXPECT_EQ(ws.eventsWritten.load(), delivered);
  // Índice: 4 segmentos (rotación a 120, 240 y 360 s), contiguos, complete, sin escritura.
  const auto& idx = w.index();
  EXPECT_TRUE(idx["complete"].get<bool>());
  ASSERT_EQ(idx["segments"].size(), 4u);
  EXPECT_EQ(modeOf(w.indexPath()), 0444u);
  const auto sums = readSums(loc.dir());
  std::uint64_t expectSeq = 0, total = 0;
  for (const auto& s : idx["segments"]) {
    const std::string f = s["file"];
    const std::string p = loc.dir() + "/" + f;
    EXPECT_EQ(modeOf(p), 0444u) << f;
    auto sha = sha256File(p);
    ASSERT_TRUE(sha.ok());
    EXPECT_EQ(s["sha256"], sha.value());
    EXPECT_EQ(sums.at(f), sha.value());
    EXPECT_EQ(s["first_evt_seq"].get<std::uint64_t>(), expectSeq) << f;
    expectSeq = s["last_evt_seq"].get<std::uint64_t>() + 1;
    total += s["n_events"].get<std::uint64_t>();
    EXPECT_TRUE(s["closed_cleanly"].get<bool>());
    // Contenido del segmento.
    std::unique_ptr<TFile> in(TFile::Open(p.c_str(), "READ"));
    ASSERT_TRUE(in && !in->IsZombie());
    auto* ev = dynamic_cast<TTree*>(in->Get("events"));
    ASSERT_NE(ev, nullptr);
    EXPECT_EQ(static_cast<std::uint64_t>(ev->GetEntries()), s["n_events"].get<std::uint64_t>());
    EXPECT_GT(dynamic_cast<TTree*>(in->Get("hv"))->GetEntries(), 0);
    // Get() de un TNamed devuelve un objeto que es del llamador.
    std::unique_ptr<TNamed> ri(in->Get<TNamed>("run_info_json"));
    ASSERT_NE(ri, nullptr);
    const auto info = nlohmann::json::parse(ri->GetTitle());
    EXPECT_EQ(info["segment"], s["seg"]);
    EXPECT_EQ(info["calibration_state"], "none");
    EXPECT_EQ(info["adc_conversion"]["NBit"], 12);
    EXPECT_EQ(info["ttt_tick_source"], "nominal");  // de dónde sale ttt_tick_ns (DAQ2-7)
    std::unique_ptr<TNamed> cv(in->Get<TNamed>("calibration_version"));
    ASSERT_NE(cv, nullptr);
    EXPECT_EQ(std::string(cv->GetTitle()), "none");
    std::unique_ptr<TNamed> hs(in->Get<TNamed>("hvmon_state_at_start"));
    EXPECT_NE(hs, nullptr);
    auto* notes = dynamic_cast<TTree*>(in->Get("notes"));
    EXPECT_EQ(notes->GetEntries(), s["seg"] == 2 ? 1 : 0) << f;
    auto* chg = dynamic_cast<TTree*>(in->Get("hv_changes"));
    EXPECT_EQ(chg->GetEntries(), s["seg"] == 2 ? 1 : 0) << f;
    // Flags: CH02 (grupo 3 = ch 33, 38, 49) apagado desde 200 s; bit 5 alrededor del cambio.
    ULong64_t seq = 0;
    Long64_t tUtc = 0;
    UShort_t chf[64];
    UInt_t evf = 0;
    ev->SetBranchAddress("evt_seq", &seq);
    ev->SetBranchAddress("t_utc_ns", &tUtc);
    ev->SetBranchAddress("ch_flags", chf);
    ev->SetBranchAddress("evt_flags", &evf);
    for (Long64_t i = 0; i < ev->GetEntries(); ++i) {
      ev->GetEntry(i);
      const double ts = static_cast<double>(tUtc - t0) / 1e9;
      EXPECT_TRUE(evf & evtflag::NoCalibration);
      if (ts > 202 && ts < 399) {
        EXPECT_TRUE(chf[33] & chflag::HvOff) << ts;
        EXPECT_FALSE(chf[34] & chflag::HvOff) << ts;
        EXPECT_FALSE(chf[33] & chflag::HvAmbiguous) << ts;
      } else if (ts > 2 && ts < 190) {
        EXPECT_FALSE(chf[33] & (chflag::HvOff | chflag::HvAmbiguous | chflag::HvUnknown)) << ts;
      } else if (ts > 195.5 && ts <= 199.9) {
        EXPECT_TRUE(chf[33] & chflag::HvAmbiguous) << ts;  // backdate de trip_backdate_s = 5 s
      }
      EXPECT_TRUE(chf[5] & chflag::HvUnknown);  // ch5 de hod2018: grupo desconocido
    }
  }
  EXPECT_EQ(total, delivered);
  EXPECT_EQ(sums.count(run.value().indexFile()), 1u);
  // Manifiesto con el run.
  const auto st = campaignStatus(loc);
  ASSERT_EQ(st["manifest"]["runs"].size(), 1u);
  EXPECT_EQ(st["manifest"]["runs"][0]["n_events"].get<std::uint64_t>(), delivered);
}

TEST(ReplayDigitizer, CropsV8FormatAndRecoversAdcExactly) {
  const std::string raw = tmpDir("replay");
  PosixFileSystem().makeDirectories(raw + "/hodo_test", 0755);
  const std::string file = raw + "/hodo_test/lectura_test_v8fmt.root";
  // Archivo con el formato de v8: Ch0..Ch63 vector<float> en mV (adc2mV del DAQ antiguo, 0x7FFF),
  // time, event, timestamp (Long_t, µs).
  auto adcOf = [](int e, int c, int k) { return 1900 + ((e * 7 + c * 13 + k * 3) % 400); };
  {
    TFile out(file.c_str(), "RECREATE");
    TTree t("data_dig", "formato v8");
    std::vector<std::vector<float>> ch(64);
    std::vector<float> tm(624);
    int event = 0;
    Long_t ts = 0;
    for (int c = 0; c < 64; ++c) t.Branch(("Ch" + std::to_string(c)).c_str(), &ch[static_cast<std::size_t>(c)]);
    t.Branch("time", &tm);
    t.Branch("event", &event);
    t.Branch("timestamp", &ts);
    for (int e = 0; e < 20; ++e) {
      event = e;
      ts = 1091270226715L + e * 1000;
      for (int c = 0; c < 64; ++c) {
        ch[static_cast<std::size_t>(c)].resize(624);
        for (int k = 0; k < 624; ++k) {
          const int adc = adcOf(e, c, k);
          ch[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)] =
              static_cast<float>((adc / ((1 << 12) - 1.0) - 1.0 + static_cast<float>(0x7FFF) / static_cast<float>(0xffff)) * 2 * 1000);
        }
      }
      t.Fill();
    }
    t.Write();
  }
  ::setenv("HODO_RAW_DIR", raw.c_str(), 1);
  ManualClock clock;
  auto rp = dig::ReplayParams::fromJson({{"file", "lectura_test_v8fmt.root"}, {"speed", 10.0}, {"max_entries", 15}});
  ASSERT_TRUE(rp.ok());
  dig::ReplayDigitizer rd(rp.value(), clock);
  dig::DigitizerConfig dc;
  ASSERT_TRUE(rd.configure(dc).ok());
  EXPECT_EQ(rd.resolvedPath(), file);
  ASSERT_TRUE(rd.start().ok());
  std::vector<dig::Event> evs;
  std::vector<UtcNs> tRead;
  for (int i = 0; i < 1000 && !rd.exhausted(); ++i) {
    clock.advance(50us);  // speed 10: un evento cada 100 µs (1 ms del original)
    dig::EventBlock b;
    ASSERT_TRUE(rd.readBlock(b, 4).ok());
    for (auto& e : b.events) {
      evs.push_back(std::move(e));
      tRead.push_back(b.tReadUtcNs);
    }
  }
  ASSERT_EQ(evs.size(), 15u);
  for (std::size_t e = 0; e < evs.size(); ++e) {
    EXPECT_EQ(evs[e].counterHw, e);
    for (int c = 0; c < 64; c += 7)
      for (int k = 0; k < 144; ++k)
        ASSERT_EQ(evs[e].adc[static_cast<std::size_t>(c * 144 + k)], adcOf(static_cast<int>(e), c, 209 + k)) << e << " " << c << " " << k;
  }
  EXPECT_LT(rd.maxRoundingResidualAdc(), 0.01);
  // TTT al ritmo de entrega: 3 ms del original / speed 10 = 300 µs, en ticks de 8 ns; coherente con
  // la hora de lectura (el evento se entrega en cuanto vence).
  EXPECT_EQ(evs[3].tttRaw, 300u * 1000u / 8u);
  EXPECT_LE(static_cast<double>(tRead[3] - tRead[0]) / 1e3, 300.0 + 50.0);
  ::unsetenv("HODO_RAW_DIR");
}

TEST(ResolveRawFile, FindsSegmentsOfSimulatedCampaignWithHodoRawDir) {
  const std::string root = tmpDir("resolve");
  PosixFileSystem fs;
  fs.makeDirectories(root + "/simA", 0755);
  fs.makeDirectories(root + "/simB", 0755);
  { std::ofstream(root + "/simA/simA_20261019T000000Z_run001_seg001.root") << "x"; }
  ::setenv("HODO_RAW_DIR", root.c_str(), 1);
  EXPECT_EQ(hodo::resolveRawFile("simA_20261019T000000Z_run001_seg001.root"), root + "/simA/simA_20261019T000000Z_run001_seg001.root");
  // Un nombre repetido en dos campañas sería ambiguo (por eso el nombre lleva la campaña).
  { std::ofstream(root + "/simB/dup.root") << "x"; }
  { std::ofstream(root + "/simA/dup.root") << "x"; }
  EXPECT_EQ(hodo::resolveRawFile("dup.root"), "");
  EXPECT_TRUE(fileNameUsedInAnyCampaign(root, "dup.root"));
  ::unsetenv("HODO_RAW_DIR");
}

TEST(SlowControlReader, MergesPartsGzAndSchemas) {
  const std::string d = tmpDir("screader");
  PosixFileSystem fs;
  fs.makeDirectories(d + "/2026", 0755);
  auto schema = [](const std::vector<std::pair<std::string, std::string>>& cols) {
    nlohmann::json c = nlohmann::json::array();
    for (const auto& x : cols) c.push_back({{"name", x.first}, {"type", x.second}});
    return nlohmann::json{{"schema", "daqv2.hv_full/x"}, {"columns", c}}.dump();
  };
  {
    std::ofstream o(d + "/2026/hv_full_20261019.csv");
    o << "t_resp_utc_ns,ch,VMon_V\n1,0,65.5\n3,0,65.6\n5,0,65.7\n";
    std::ofstream(d + "/2026/hv_full_20261019.csv.schema.json")
        << schema({{"t_resp_utc_ns", "int64"}, {"ch", "int"}, {"VMon_V", "float64"}});
  }
  {
    std::ofstream o(d + "/2026/hv_full_20261019.2.csv");
    o << "t_resp_utc_ns,ch,VMon_V,bd_Temp_degC\n2,1,65.1,25\n4,1,65.2,26\n";
    std::ofstream(d + "/2026/hv_full_20261019.2.csv.schema.json")
        << schema({{"t_resp_utc_ns", "int64"}, {"ch", "int"}, {"VMon_V", "float64"}, {"bd_Temp_degC", "float64"}});
  }
  ASSERT_TRUE(gzipFileVerified(d + "/2026/hv_full_20261019.2.csv", d + "/2026/hv_full_20261019.2.csv.gz", true).ok());
  {
    std::ofstream o(d + "/2026/hv_full_20261020.csv");
    o << "t_resp_utc_ns,ch,VMon_V\n9,0,65.9\n";
  }
  std::string err;
  std::unique_ptr<TTree> t(io::readSlowControl(d, "hv_full", "20261019", "20261020", &err));
  ASSERT_TRUE(t) << err;
  ASSERT_EQ(t->GetEntries(), 6);
  Long64_t ts = 0;
  Int_t ch = 0, part = 0;
  Double_t temp = 0, vmon = 0;
  t->SetBranchAddress("t_resp_utc_ns", &ts);
  t->SetBranchAddress("ch", &ch);
  t->SetBranchAddress("part", &part);
  t->SetBranchAddress("bd_Temp_degC", &temp);
  t->SetBranchAddress("VMon_V", &vmon);
  const Long64_t expT[] = {1, 2, 3, 4, 5, 9};
  const Int_t expPart[] = {1, 2, 1, 2, 1, 1};
  for (Long64_t i = 0; i < 6; ++i) {
    t->GetEntry(i);
    EXPECT_EQ(ts, expT[i]);
    EXPECT_EQ(part, expPart[i]);
    if (part == 2) EXPECT_GE(temp, 25.0);
    else EXPECT_TRUE(std::isnan(temp));
  }
  EXPECT_EQ(std::string(t->GetBranch("ch")->GetLeaf("ch")->GetTypeName()), "Int_t");
  // Con dirección en una sola rama, GetEntry no escribe en las variables locales del lector (ASan).
  std::unique_ptr<TTree> t2(io::readSlowControl(d, "hv_full", "20261019", "20261020", &err));
  Long64_t onlyT = 0;
  t2->SetBranchAddress("t_resp_utc_ns", &onlyT);
  for (Long64_t i = 0; i < t2->GetEntries(); ++i) t2->GetEntry(i);
  EXPECT_EQ(onlyT, 9);
  EXPECT_EQ(t2->GetEntries("ch == 1"), 2);
  EXPECT_EQ(io::readSlowControl(d, "hv", "20261019", "20261019", &err), nullptr);
  EXPECT_NE(err.find("no hay archivos"), std::string::npos);
  {
    std::ofstream o(d + "/2026/notes_20261019.jsonl");
    o << R"({"t_utc_ns":7,"user":"rene","host":"h","text":"hola"})" << "\n";
  }
  std::unique_ptr<TTree> n(io::readSlowControl(d, "notes", "20261019", "20261019", &err));
  ASSERT_TRUE(n) << err;
  EXPECT_EQ(n->GetEntries(), 1);
  EXPECT_EQ(n->GetEntries("text == \"hola\""), 1);
}

// pulses/amp de la página de daqrun (DAQ2-8): SetBinContent repetido suma entradas en cada refresco;
// refreshEntries() las deja iguales a la suma de contenidos.
TEST(HistUtil, RefreshEntriesEqualsSumOfContentsAfterRepeatedRefills) {
  TH2F h("amp_test", "", 4, 0, 4, 3, 0, 3);
  h.SetDirectory(nullptr);
  for (int refresh = 0; refresh < 3; ++refresh)
    for (int x = 1; x <= 4; ++x)
      for (int y = 1; y <= 3; ++y) h.SetBinContent(x, y, x * y);
  EXPECT_EQ(h.GetEntries(), 36.0);  // 3 refrescos × 12 bins: el error visto en la página
  daq::io::refreshEntries(h);
  EXPECT_EQ(h.GetEntries(), 60.0);  // (1+2+3+4) × (1+2+3)
  TH1D g("dt_test", "", 3, 0, 3);
  g.SetDirectory(nullptr);
  g.SetBinContent(1, 7);
  g.SetBinContent(1, 7);
  daq::io::refreshEntries(g);
  EXPECT_EQ(g.GetEntries(), 7.0);
}
