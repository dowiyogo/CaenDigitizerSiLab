// selftrig_muon.C — DAQ2-13, parte A.3 (solo lectura). Amplitudes en tiempo de hod2018 por grupo del
// digitizador en runs de trigger externo, para la eficiencia a muones contra el umbral.
//  - Por evento y grupo 0-3 (ch 0-7, 8-15 = hod2018 X; 16-23, 24-31 = Y): máximo sobre los canales del
//    grupo (sin los canales de masked, p. ej. "2,9,15": fuera del auto-trigger) de (línea base recortada de [0,40) − mínimo) en la ventana signal [49,58) y en la de control
//    [40,49) (accidentales, mismo ancho), en mV y en p.e. (G de la calibración; canales sin G: mediana).
// Salida: CSV evento,s0_mV..s3_mV,c0_mV..c3_mV,s0_pe..s3_pe,c0_pe..c3_pe.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/validation/selftrig_muon.C+("<run>_seg*.root", "<calibración.json>", "<salida.csv>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/pulse_finder.h"

void selftrig_muon(const char* pattern, const char* calibJson, const char* out, const char* masked = "") {
  std::vector<bool> skip(64, false);
  {
    std::string m = masked ? masked : "";
    std::size_t p = 0;
    while (p < m.size()) {
      const std::size_t q = m.find(',', p);
      const std::string tok = m.substr(p, q == std::string::npos ? std::string::npos : q - p);
      if (!tok.empty()) skip[static_cast<std::size_t>(std::stoi(tok))] = true;
      if (q == std::string::npos) break;
      p = q + 1;
    }
  }
  std::vector<double> G(64, NAN);
  {
    std::ifstream in(calibJson);
    auto j = nlohmann::json::parse(in);
    for (const auto& e : j["channels"])
      if (e.contains("G_mV") && e["G_mV"].is_number()) G[e["ch"].get<int>()] = e["G_mV"].get<double>();
  }
  std::vector<double> ok;
  for (int c = 0; c < 32; ++c)
    if (std::isfinite(G[c])) ok.push_back(G[c]);
  std::sort(ok.begin(), ok.end());
  for (int c = 0; c < 32; ++c)
    if (!std::isfinite(G[c])) G[c] = ok[ok.size() / 2];
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchAddress("adc", adc.data());
  const double lsb = daq::pulse::lsbMv();
  FILE* f = fopen(out, "w");
  fprintf(f, "ev,s0_mV,s1_mV,s2_mV,s3_mV,c0_mV,c1_mV,c2_mV,c3_mV,s0_pe,s1_pe,s2_pe,s3_pe,c0_pe,c1_pe,c2_pe,c3_pe\n");
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    double sm[4] = {0, 0, 0, 0}, cm[4] = {0, 0, 0, 0}, sp[4] = {0, 0, 0, 0}, cp[4] = {0, 0, 0, 0};
    for (int c = 0; c < 32; ++c) {
      if (skip[static_cast<std::size_t>(c)]) continue;
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      UShort_t ms = x[49], mc = x[40];
      for (int k = 49; k < 58; ++k) ms = std::min(ms, x[k]);
      for (int k = 40; k < 49; ++k) mc = std::min(mc, x[k]);
      const double as = (bl - ms) * lsb, ac = (bl - mc) * lsb;
      const int g = c / 8;
      sm[g] = std::max(sm[g], as);
      cm[g] = std::max(cm[g], ac);
      sp[g] = std::max(sp[g], as / G[c]);
      cp[g] = std::max(cp[g], ac / G[c]);
    }
    fprintf(f, "%lld,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", e, sm[0], sm[1], sm[2], sm[3], cm[0], cm[1],
            cm[2], cm[3], sp[0], sp[1], sp[2], sp[3], cp[0], cp[1], cp[2], cp[3]);
  }
  fclose(f);
  printf("eventos %lld, muestras %d\n", ne, n);
}
