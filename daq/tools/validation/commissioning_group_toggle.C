// commissioning_group_toggle.C — DAQ2-10, parte C (solo lectura). Prueba de grupos HV: un grupo
// apagado a la vez con SW trigger a 100 Hz y registro completo. Ventanas (inicio/fin en s UTC) del
// controlador (/data/daq/logs/daq2-10_toggle.jsonl, eventos "window"): ref_inicio, off (grupo g),
// ref_entre (tras volver el grupo g), ref_final.
//  - Pulsos oscuros por canal en las muestras [0,49) + [102,622) (todo el registro salvo la zona de
//    señal, sin muón con SW trigger), buscador del DAQ (pulse_finder.h), separación 4; umbrales:
//    hod2018 3 p.e. con G de v8 (solo referencia), hod2019 1.5 mV y 3 mV.
//  - Por paso "off" g y canal: n, μ = tasa de la referencia más cercana × tiempo expuesto,
//    razón n/μ y P(N ≤ n | μ). Asignación: el canal pertenece al grupo g si en el paso g la razón es
//    < 0.1 y p < 1e-6 y en ningún otro paso pasa lo mismo.
//  - Corriente: IMon de los canales HV al inicio y fin de cada ventana (del controlador).
// Salidas: CSV por paso y canal (out.csv) y tabla de asignación (out_map.csv).
// Uso: root -l -b -q 'daq/tools/validation/commissioning_group_toggle.C+("<run>_seg*.root", "toggle.jsonl", "out")'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>
#include <TLeaf.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/pulse_finder.h"
#include "daq/trigger_interval_stats.h"

namespace {
struct Win {
  std::string tag;
  int group = -1;
  double t0 = 0, t1 = 0;
  long events = 0;
  std::vector<std::array<long, 3>> n = std::vector<std::array<long, 3>>(64);
};
}  // namespace

void commissioning_group_toggle(const char* pattern, const char* jsonl, const char* out,
                                const char* gains = "outputs/pe_calibration_v8.txt") {
  std::vector<double> G(64, NAN);
  {
    std::ifstream in(gains);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int ch;
      std::string pl;
      double p0, g;
      if (ss >> ch >> pl >> p0 >> g && ch >= 0 && ch < 64) G[ch] = g;
    }
  }
  std::vector<Win> wins;
  {
    std::ifstream in(jsonl);
    std::string line;
    while (std::getline(in, line)) {
      auto j = nlohmann::json::parse(line);
      if (j.value("ev", "") != "window") continue;
      Win w;
      w.tag = j["tag"];
      w.group = j["group"].is_null() ? -1 : j["group"].get<int>();
      w.t0 = j["start"];
      w.t1 = j["end"];
      wins.push_back(w);
    }
  }
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  Long64_t tread = 0;
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchStatus("t_read_utc_ns", 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("t_read_utc_ns", &tread);
  const double lsb = daq::pulse::lsbMv();
  const std::vector<daq::pulse::Region> zones{{0, 49}, {102, std::min(622, n)}};
  const double exposure = (49 + std::min(622, n) - 102) * 16e-9;
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    const Long64_t le = t.LoadTree(e);  // índice local del archivo de la cadena
    t.GetTree()->GetBranch("t_read_utc_ns")->GetEntry(le);
    const double ts = static_cast<double>(tread) * 1e-9;
    Win* w = nullptr;
    for (auto& x : wins)
      if (ts >= x.t0 && ts < x.t1) w = &x;
    if (!w) continue;
    t.GetEntry(e);
    ++w->events;
    for (int c = 0; c < 64; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      std::vector<daq::pulse::Pulse> p;
      daq::pulse::findPulses(x, n, bl, lsb, 1.0, 4, zones, p);
      const double thr[3] = {c < 32 ? 3.0 * G[c] : 1.5, c < 32 ? NAN : 3.0, NAN};
      for (const auto& q : p)
        for (int k = 0; k < 2; ++k)
          if (std::isfinite(thr[k]) && q.amplitudeMv > thr[k]) ++w->n[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)];
    }
  }
  // Referencias: todas las ventanas con todo encendido; para cada paso off, la más cercana en tiempo.
  auto nearestRef = [&](const Win& off) -> const Win* {
    const Win* best = nullptr;
    double bd = 1e18;
    for (const auto& r : wins)
      if (r.tag != "off" && r.events > 0) {
        const double d = std::fabs(0.5 * (r.t0 + r.t1) - 0.5 * (off.t0 + off.t1));
        if (d < bd) {
          bd = d;
          best = &r;
        }
      }
    return best;
  };
  FILE* f = fopen((std::string(out) + ".csv").c_str(), "w");
  fprintf(f, "step,group,ch,threshold,ref_tag,events,n,mu,ratio,p_drop\n");
  std::map<int, std::vector<int>> dropIn;  // canal -> pasos donde cae
  printf("ventanas: ");
  for (const auto& w : wins) printf("%s%s(%ld) ", w.tag.c_str(), w.group >= 0 ? std::to_string(w.group).c_str() : "", w.events);
  printf("\n");
  for (const auto& w : wins) {
    if (w.tag != "off") continue;
    const Win* r = nearestRef(w);
    if (!r) continue;
    std::vector<int> drops;
    for (int c = 0; c < 64; ++c)
      for (int k = 0; k < 2; ++k) {
        const double rr = static_cast<double>(r->n[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)]) / (r->events * exposure);
        const double mu = rr * w.events * exposure;
        const long nn = w.n[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)];
        if (!(mu > 0)) continue;
        const double p = daq::tis::detail::gammaQ(static_cast<double>(nn) + 1.0, mu);
        const char* thrName = c < 32 ? (k == 0 ? "3pe_v8" : "") : (k == 0 ? "1.5mV" : "3mV");
        if (!*thrName) continue;
        fprintf(f, "%d,%d,%d,%s,%s%s,%ld,%ld,%.2f,%.4f,%.3g\n", w.group, w.group + 1, c, thrName, r->tag.c_str(),
                r->group >= 0 ? std::to_string(r->group).c_str() : "", w.events, nn, mu, nn / mu, p);
        if (k == 0 && nn / mu < 0.1 && p < 1e-6) drops.push_back(c);
      }
    for (int c : drops) dropIn[c].push_back(w.group);
    printf("grupo %d (CH%02d) apagado: %ld eventos, ref %s%s; canales que caen (razón < 0.1, p < 1e-6):", w.group + 1, w.group, w.events,
           r->tag.c_str(), r->group >= 0 ? std::to_string(r->group).c_str() : "");
    for (int c : drops) printf(" %d", c);
    printf("\n");
  }
  fclose(f);
  FILE* m = fopen((std::string(out) + "_map.csv").c_str(), "w");
  fprintf(m, "ch,measured_group,steps_dropped\n");
  for (int c = 0; c < 64; ++c) {
    const auto it = dropIn.find(c);
    std::string steps;
    if (it != dropIn.end())
      for (int g : it->second) steps += (steps.empty() ? "" : " ") + std::to_string(g + 1);
    const int g = it != dropIn.end() && it->second.size() == 1 ? it->second[0] + 1 : 0;
    fprintf(m, "%d,%d,%s\n", c, g, steps.c_str());
  }
  fclose(m);
}
