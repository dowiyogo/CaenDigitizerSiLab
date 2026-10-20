// panels_root.cpp
#include "daq/io/panels_root.h"

#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THttpServer.h>

#include <algorithm>
#include <cstdio>

#include "daq/io/hist_util.h"

namespace daq::io {

namespace {
double hoursOf(UtcNs t, UtcNs t0) { return static_cast<double>(t - t0) / 3.6e12; }
}  // namespace

PanelRoot::PanelRoot() = default;
PanelRoot::~PanelRoot() = default;

TH2F* PanelRoot::h2(const std::string& folder, const std::string& name, const std::string& title, int nx, double x0, double x1,
                    int ny, double y0, double y1) {
  auto it = obj_.find(name);
  if (it != obj_.end()) {
    auto* h = static_cast<TH2F*>(it->second.second.get());
    if (h->GetNbinsX() == nx && h->GetXaxis()->GetXmax() == x1) return h;
    h->SetBins(nx, x0, x1, ny, y0, y1);  // eje de tiempo que crece
    return h;
  }
  auto h = std::make_unique<TH2F>(name.c_str(), title.c_str(), nx, x0, x1, ny, y0, y1);
  h->SetDirectory(nullptr);
  h->SetStats(false);
  TH2F* raw = h.get();
  obj_[name] = {folder, std::move(h)};
  return raw;
}

TH1* PanelRoot::h1(const std::string& folder, const std::string& name, const std::string& title, int nx, double x0, double x1) {
  auto it = obj_.find(name);
  if (it != obj_.end()) return static_cast<TH1*>(it->second.second.get());
  auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), nx, x0, x1);
  h->SetDirectory(nullptr);
  TH1* raw = h.get();
  obj_[name] = {folder, std::move(h)};
  return raw;
}

TGraph* PanelRoot::graph(const std::string& folder, const std::string& name, const std::string& title) {
  auto it = obj_.find(name);
  if (it != obj_.end()) return static_cast<TGraph*>(it->second.second.get());
  auto g = std::make_unique<TGraph>();
  g->SetName(name.c_str());
  g->SetTitle(title.c_str());
  TGraph* raw = g.get();
  obj_[name] = {folder, std::move(g)};
  return raw;
}

void PanelRoot::update(const daqrun::PanelView& p, const daqrun::DataFlagView& d, double leakRatio) {
  // Histogramas de eventos (acumulado y última hora).
  for (const auto& ph : p.hists) {
    std::string folder = "hits";
    if (ph.name.rfind("intime", 0) == 0 || ph.name == "min_sample") folder = "intime";
    for (int v = 0; v < 2; ++v) {
      const std::string name = ph.name + (v ? "_1h" : "");
      const std::string title = ph.title.substr(0, ph.title.find(';')) + (v ? " [última hora]" : " [run]") +
                                ph.title.substr(std::min(ph.title.size(), ph.title.find(';')));
      TH1* h = ph.ny > 1 ? static_cast<TH1*>(h2(folder, name, title, ph.nx, ph.x0, ph.x1, ph.ny, ph.y0, ph.y1))
                         : h1(folder, name, title, ph.nx, ph.x0, ph.x1);
      for (int iy = 0; iy < ph.ny; ++iy)
        for (int ix = 0; ix < ph.nx; ++ix) {
          const auto i = static_cast<std::size_t>(iy * ph.nx + ix);
          const double c = v ? static_cast<double>(ph.hour[i]) : static_cast<double>(ph.cum[i]);
          if (ph.ny > 1)
            h->SetBinContent(ix + 1, iy + 1, c);
          else
            h->SetBinContent(ix + 1, c);
        }
      refreshEntries(*h);
    }
  }
  // Persistencia y último evento.
  const int N = std::max(1, p.recordLength);
  for (std::size_t g = 0; g < 8; ++g) {
    char name[32], title[160];
    std::snprintf(name, sizeof name, "persistence_g%zu", g);
    std::snprintf(title, sizeof title, "Persistencia grupo %zu (canales %zu-%zu, últimos ~1000 eventos);muestra;ADC", g, 8 * g, 8 * g + 7);
    TH2F* h = h2("scope", name, title, N, -0.5, N - 0.5, p.persistAdcBins, 0, 16.0 * p.persistAdcBins);
    const auto& v = p.persistence[g];
    for (int k = 0; k < N; ++k)
      for (int a = 0; a < p.persistAdcBins; ++a) {
        const auto i = static_cast<std::size_t>(k * p.persistAdcBins + a);
        h->SetBinContent(k + 1, a + 1, i < v.size() ? v[i] : 0);
      }
    refreshEntries(*h);
  }
  {
    TH2F* h = h2("scope", "last_event", "Último evento (ADC);muestra;canal", N, -0.5, N - 0.5, 64, -0.5, 63.5);
    for (int c = 0; c < 64; ++c)
      for (int k = 0; k < N; ++k) {
        const auto i = static_cast<std::size_t>(c * N + k);
        h->SetBinContent(k + 1, c + 1, i < p.lastEvent.size() ? p.lastEvent[i] : 0);
      }
    refreshEntries(*h);
  }
  // Mapa canal × tiempo (períodos) y vivos por plano.
  const UtcNs t0 = p.rateMapT.empty() ? p.tUtcNs : p.rateMapT.front();
  const int nt = std::max(1, static_cast<int>(p.rateMapT.size()));
  {
    const double span = std::max(1.0 / 60, hoursOf(p.rateMapT.empty() ? t0 : p.rateMapT.back(), t0));
    TH2F* h = h2("overview", "rate_map", "Tasa fuera de tiempo por canal (umbral PROVISIONAL) [Hz];horas desde el primer período;canal",
                 nt, -span / std::max(1, nt - 1) / 2, span + span / std::max(1, nt - 1) / 2, 64, -0.5, 63.5);
    for (int i = 0; i < static_cast<int>(p.rateMap.size()); ++i)
      for (int c = 0; c < 64; ++c) h->SetBinContent(i + 1, c + 1, p.rateMap[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)]);
    refreshEntries(*h);
  }
  const char* planes[4] = {"hod2018_X", "hod2018_Y", "hod2019_X", "hod2019_Y"};
  for (int q = 0; q < 4; ++q) {
    TGraph* g = graph("overview", std::string("live_") + planes[q],
                      std::string("Canales vivos ") + planes[q] + " (> 2 LSB por período);horas;canales");
    g->Set(0);
    for (std::size_t i = 0; i < p.livePerPlane.size(); ++i) g->AddPoint(hoursOf(p.rateMapT[i], t0), p.livePerPlane[i][static_cast<std::size_t>(q)]);
  }
  // Ganancia por hora: G (hod2018) y tasa > 3 mV (hod2019) por canal y hora.
  {
    std::vector<UtcNs> hours;
    for (const auto& gp : p.gains)
      if (std::find(hours.begin(), hours.end(), gp.hourStartUtcNs) == hours.end()) hours.push_back(gp.hourStartUtcNs);
    const int nh = std::max(1, static_cast<int>(hours.size()));
    TH2F* gG = h2("gain", "gain_hod2018", "G por hora hod2018 [mV/p.e.] PROVISIONAL (0: sin ajuste aceptado);hora del run;canal", nh, -0.5,
                  nh - 0.5, 32, -0.5, 31.5);
    TH2F* gS = h2("gain", "gain_status_hod2018",
                  "Estado del ajuste (1 ok, 2 cambio > 20 %, 3 ambiguo G/2-2G, 4 rechazado por espaciado, 5 insuficiente);hora;canal",
                  nh, -0.5, nh - 0.5, 32, -0.5, 31.5);
    TH2F* gD = h2("gain", "dark_rate_hod2019", "Tasa > 3 mV por hora hod2019 [Hz] (indicador de ganancia);hora;canal", nh, -0.5, nh - 0.5,
                  32, 31.5, 63.5);
    gG->Reset();
    gS->Reset();
    gD->Reset();
    for (const auto& gp : p.gains) {
      const int ih = static_cast<int>(std::find(hours.begin(), hours.end(), gp.hourStartUtcNs) - hours.begin()) + 1;
      if (gp.status == "dark_rate") {
        gD->SetBinContent(ih, gp.ch - 31, gp.value);
        continue;
      }
      const bool acc = gp.status == "ok" || gp.status == "changed";
      gG->SetBinContent(ih, gp.ch + 1, acc ? gp.value : 0);
      const int code = gp.status == "ok" ? 1 : gp.status == "changed" ? 2 : gp.status == "ambiguous" ? 3 : gp.status == "rejected" ? 4 : 5;
      gS->SetBinContent(ih, gp.ch + 1, code);
    }
    for (TH1* h : {static_cast<TH1*>(gG), static_cast<TH1*>(gS), static_cast<TH1*>(gD)}) refreshEntries(*h);
  }
  // Bits de datos: tasa y base por unidad contra el tiempo.
  UtcNs td0 = 0;
  for (const auto& [k, h] : d.history)
    if (!h.empty() && (td0 == 0 || h.front().t < td0)) td0 = h.front().t;
  for (const auto& u : d.units) {
    auto it = d.history.find(u.key);
    if (it == d.history.end()) continue;
    std::string base = u.key;
    std::replace(base.begin(), base.end(), ':', '_');
    std::replace(base.begin(), base.end(), ' ', '_');
    const bool isImon = u.kind == "imon";
    const std::string unit = isImon ? "A" : "Hz por canal";
    const std::string folder = u.kind == "face" ? "overview" : "dataflags";
    TGraph* g = graph(folder, "rate_" + base, u.name + (isImon ? " IMon" : " tasa") + ";horas;" + unit);
    g->Set(0);
    auto bandIt = obj_.find("baseline_" + base);
    TGraphAsymmErrors* band = nullptr;
    if (bandIt == obj_.end()) {
      auto b = std::make_unique<TGraphAsymmErrors>();
      b->SetName(("baseline_" + base).c_str());
      b->SetTitle((u.name + " base (banda: hasta el umbral del aviso);horas;" + unit).c_str());
      band = b.get();
      obj_["baseline_" + base] = {folder, std::move(b)};
    } else {
      band = static_cast<TGraphAsymmErrors*>(bandIt->second.second.get());
    }
    band->Set(0);
    const double k = u.kind == "face" ? leakRatio : 1.0;
    for (const auto& pt : it->second) {
      const double x = hoursOf(pt.t, td0);
      g->AddPoint(x, pt.rate);
      if (pt.baseline > 0) {
        const int i = band->GetN();
        band->SetPoint(i, x, pt.baseline);
        band->SetPointError(i, 0, 0, 0, pt.baseline * (k - 1.0));
      }
    }
  }
}

void PanelRoot::registerIn(THttpServer& s) {
  for (auto& [name, fo] : obj_) {
    TObject* o = fo.second.get();
    if (std::find(registered_.begin(), registered_.end(), o) != registered_.end()) continue;
    s.Register(("/" + fo.first).c_str(), o);
    registered_.push_back(o);
  }
}

void PanelRoot::unregisterFrom(THttpServer& s) {
  for (TObject* o : registered_) s.Unregister(o);
  registered_.clear();
}

void PanelRoot::write(TDirectory& dir) const {
  for (const auto& [name, fo] : obj_) {
    TDirectory* sub = dir.GetDirectory(fo.first.c_str());
    if (!sub) sub = dir.mkdir(fo.first.c_str());
    sub->WriteTObject(fo.second.get(), name.c_str(), "Overwrite");
  }
}

std::vector<std::string> PanelRoot::defaultDrawItems() {
  return {"trigger/trigger_rate", "overview/rate_map", "overview/rate_face_hod2019_Trasera", "overview/live_hod2019_X",
          "intime/min_sample", "flags/ch_flags"};
}

}  // namespace daq::io
