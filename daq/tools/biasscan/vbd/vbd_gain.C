// vbd_gain.C — copia de analysis/followup8_2026_10_31/gain_spectra.C (seguimiento 2026-10-31 (8)) para el análisis del barrido de V_bd
// (seguimiento 2026-10-31 (9), vbd_analysis.py). Solo lectura: variante de daq/tools/characterization/charact_gain.C (reglas de la fase 5,
// las de c1) que además guarda el espectro de cada canal (bins de 1 LSB = 0.4884 mV, hasta 100 mV) y ajusta G con gain_fit.h con dos juegos de
// opciones: "c1" (por omisión: picos entre 5 y 26 mV, al menos 3 picos de 20 cuentas) y "low" (picos desde 2.5 mV, para tensiones bajas donde
// la serie de crosstalk sobre 5 mV ya no tiene 3 picos con cuentas suficientes). Mismos pulsos que charact_gain.C: regiones fuera de tiempo
// [0,40) ∪ [40,49) ∪ [102,n), pedestal = media recortada de [0,40), umbral 0.9 mV, separación 4 muestras.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'vbd_gain.C+("<run>_seg*.root", "<semillas.json>", "<prefijo>")'   -> <prefijo>_fits.csv, <prefijo>_spec.csv
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/gain_fit.h"
#include "daq/pulse_finder.h"

void vbd_gain(const char* pattern, const char* calibJson, const char* prefix) {
  std::vector<double> G(32, 4.0);
  {
    std::ifstream in(calibJson);
    auto j = nlohmann::json::parse(in);
    for (const auto& e : j["channels"])
      if (e["ch"].get<int>() < 32 && e.contains("G_mV") && e["G_mV"].is_number()) G[e["ch"].get<int>()] = e["G_mV"].get<double>();
  }
  struct Win {
    std::string label;
    double t0, t1;
    long events = 0;
    std::vector<std::vector<unsigned>> spec;
  };
  std::vector<Win> wins{{"all", 0.0, 1e11}};
  const double lsb = daq::pulse::lsbMv();
  const int nb = static_cast<int>(100.0 / lsb);
  for (auto& w : wins) w.spec.assign(32, std::vector<unsigned>(static_cast<std::size_t>(nb) + 1, 0));
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  Long64_t tutc = 0;
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchStatus("t_utc_ns", 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("t_utc_ns", &tutc);
  const std::vector<daq::pulse::Region> oot{{0, 40}, {40, 49}, {102, n}};
  std::vector<daq::pulse::Pulse> p;
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double ts = static_cast<double>(tutc) * 1e-9;
    std::vector<Win*> in;
    for (auto& w : wins)
      if (ts >= w.t0 && ts < w.t1) in.push_back(&w);
    if (in.empty()) continue;
    for (auto* w : in) ++w->events;
    for (int c = 0; c < 32; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      p.clear();
      daq::pulse::findPulses(x, n, bl, lsb, 0.9, 4, oot, p);
      for (const auto& u : p) {
        const auto b = static_cast<std::size_t>(std::min(nb, static_cast<int>(u.amplitudeMv / lsb)));
        for (auto* w : in) ++w->spec[static_cast<std::size_t>(c)][b];
      }
    }
  }
  const auto& w = wins[0];
  FILE* f = fopen((std::string(prefix) + "_fits.csv").c_str(), "w");
  fprintf(f, "ch,events,options,seed_mV,G_mV,status,peaks,positions_mV\n");
  for (int c = 0; c < 32; ++c)
    for (int k = 0; k < 2; ++k) {
      daq::gain::Options o;
      if (k == 1) o.loMv = 2.5;
      const auto r = daq::gain::fitGain(w.spec[static_cast<std::size_t>(c)], lsb, G[static_cast<std::size_t>(c)], o);
      std::string pos;
      for (double v : r.positionsMv) pos += (pos.empty() ? "" : " ") + std::to_string(v);
      fprintf(f, "%d,%ld,%s,%.5f,%.5f,%s,%d,%s\n", c, w.events, k ? "low" : "c1", G[static_cast<std::size_t>(c)], r.gainMv, r.status.c_str(), r.peaks, pos.c_str());
    }
  fclose(f);
  FILE* fs = fopen((std::string(prefix) + "_spec.csv").c_str(), "w");
  fprintf(fs, "ch,bin_lo_mV,n\n");
  for (int c = 0; c < 32; ++c)
    for (int b = 0; b <= nb; ++b) fprintf(fs, "%d,%.4f,%u\n", c, b * lsb, w.spec[static_cast<std::size_t>(c)][static_cast<std::size_t>(b)]);
  fclose(fs);
  printf("eventos %lld\n", ne);
}
