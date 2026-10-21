// commissioning_neighbours.C — DAQ2-10, parte E (solo lectura). Coincidencias entre tiras vecinas y
// f_w por canal de hod2018 con trigger externo (handoff 6), para comprobar el orden de coordenadas del
// mapeo de octubre y encontrar barras con poca luz.
//  - Impacto: amplitud en tiempo (línea base recortada de [0,40) − mínimo en [49,58)) > thr_pe × G del
//    canal (G de la calibración de octubre, JSON daqv2.pe_calibration/1; canales sin G: mediana).
//    Accidentales: la misma cuenta en una ventana de 9 muestras fuera de señal, [111,120).
//  - Matriz 16 × 16 por plano en el orden de tiras del mapeo (pares de tiras encendidas en el mismo
//    evento, i ≠ j); por plano: fracción de pares con |Δtira| = 1 y, por tira, si su compañera más
//    frecuente es vecina.
//  - f_w = C_w / S_c: S_c impactos del canal menos accidentales; C_w = Σ 1/n_plano sobre los eventos con
//    al menos un impacto en X y uno en Y de hod2018 (peso 1/(nX·nY) por par repartido entre los canales
//    del evento), menos lo mismo en la ventana de accidentales.
// Uso: root -l -b -q 'daq/tools/validation/commissioning_neighbours.C+("<run>_seg*.root", "<calibración.json>", "<mapeo.json>", 10)'
R__ADD_INCLUDE_PATH(daq/core/include)
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

void commissioning_neighbours(const char* pattern, const char* calibJson, const char* mappingJson, double thrPe = 10) {
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
  const double gMed = ok.empty() ? 4.0 : ok[ok.size() / 2];
  std::vector<int> plane(64, -1), strip(64, -1);
  {
    std::ifstream in(mappingJson);
    auto j = nlohmann::json::parse(in);
    for (const auto& e : j["channels"]) {
      const int c = e["dig_ch"];
      if (c >= 32) continue;
      plane[c] = e["plane"] == "Y" ? 1 : 0;
      strip[c] = e["strip"];
    }
  }
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchAddress("adc", adc.data());
  const double lsb = daq::pulse::lsbMv();
  std::vector<std::vector<double>> pairs(2, std::vector<double>(256, 0));
  std::vector<double> S(32, 0), Sacc(32, 0), Cw(32, 0), Cwacc(32, 0);
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    for (int w = 0; w < 2; ++w) {  // 0: señal [49,58), 1: accidentales [111,120)
      const int s0 = w == 0 ? 49 : 111, s1 = w == 0 ? 58 : 120;
      std::vector<int> hits[2];
      for (int c = 0; c < 32; ++c) {
        const UShort_t* x = adc.data() + c * n;
        const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
        UShort_t mn = x[s0];
        for (int k = s0; k < s1; ++k) mn = std::min(mn, x[k]);
        const double g = std::isfinite(G[c]) ? G[c] : gMed;
        if ((bl - mn) * lsb > thrPe * g && plane[c] >= 0) hits[plane[c]].push_back(c);
      }
      for (int p = 0; p < 2; ++p)
        for (int c : hits[p]) (w == 0 ? S : Sacc)[static_cast<std::size_t>(c)] += 1;
      if (w == 0)
        for (int p = 0; p < 2; ++p)
          for (int a : hits[p])
            for (int b : hits[p])
              if (a != b) pairs[static_cast<std::size_t>(p)][static_cast<std::size_t>(strip[a] * 16 + strip[b])] += 1;
      if (!hits[0].empty() && !hits[1].empty())
        for (int p = 0; p < 2; ++p)
          for (int c : hits[p]) (w == 0 ? Cw : Cwacc)[static_cast<std::size_t>(c)] += 1.0 / hits[p].size();
    }
  }
  printf("eventos %lld, umbral %.0f p.e. (G mediana %.3f mV para canales sin G)\n", ne, thrPe, gMed);
  const char* pn[2] = {"X", "Y"};
  for (int p = 0; p < 2; ++p) {
    double tot = 0, adj = 0;
    for (int i = 0; i < 16; ++i)
      for (int j = 0; j < 16; ++j) {
        const double v = pairs[static_cast<std::size_t>(p)][static_cast<std::size_t>(i * 16 + j)];
        tot += v;
        if (std::abs(i - j) == 1) adj += v;
      }
    printf("hod2018 %s: pares %.0f, con |Δtira| = 1: %.0f (%.3f); por tira (tira: compañera más frecuente, vecina?):\n  ", pn[p], tot / 2,
           adj / 2, tot > 0 ? adj / tot : 0);
    int bad = 0;
    for (int i = 0; i < 16; ++i) {
      int best = -1;
      double bv = -1;
      for (int j = 0; j < 16; ++j) {
        const double v = pairs[static_cast<std::size_t>(p)][static_cast<std::size_t>(i * 16 + j)];
        if (j != i && v > bv) {
          bv = v;
          best = j;
        }
      }
      const bool nb = std::abs(best - i) == 1;
      bad += !nb;
      printf("%d:%d%s ", i, best, nb ? "" : "*");
    }
    printf("\n  tiras cuya compañera más frecuente no es vecina: %d de 16\n", bad);
    printf("  matriz (fila = tira, columnas 0-15):\n");
    for (int i = 0; i < 16; ++i) {
      printf("   %2d:", i);
      for (int j = 0; j < 16; ++j) printf(" %4.0f", pairs[static_cast<std::size_t>(p)][static_cast<std::size_t>(i * 16 + j)]);
      printf("\n");
    }
  }
  printf("f_w por canal (S = impactos − accidentales, C_w − accidentales):\n");
  std::vector<double> fx[2];
  for (int c = 0; c < 32; ++c) {
    const double s = S[static_cast<std::size_t>(c)] - Sacc[static_cast<std::size_t>(c)];
    const double cw = Cw[static_cast<std::size_t>(c)] - Cwacc[static_cast<std::size_t>(c)];
    const double f = s > 0 ? cw / s : NAN;
    const double ef = s > 0 && cw > 0 ? f * std::sqrt(1.0 / cw + 1.0 / s) : NAN;
    if (plane[c] >= 0 && std::isfinite(f)) fx[plane[c]].push_back(f);
    printf("  ch%02d %s tira %2d: S %.0f (acc %.0f), C_w %.1f (acc %.1f), f_w %.3f ± %.3f\n", c, plane[c] >= 0 ? pn[plane[c]] : "?", strip[c],
           S[static_cast<std::size_t>(c)], Sacc[static_cast<std::size_t>(c)], Cw[static_cast<std::size_t>(c)], Cwacc[static_cast<std::size_t>(c)], f,
           ef);
  }
  for (int p = 0; p < 2; ++p) {
    auto v = fx[p];
    std::sort(v.begin(), v.end());
    if (!v.empty()) printf("mediana f_w %s: %.3f (%zu canales)\n", pn[p], v[v.size() / 2], v.size());
  }
}
