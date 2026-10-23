// charact_extract.C — DAQ2-15, parte B (solo lectura). Extracción compacta de runs con trigger externo para la
// caracterización del detector (análisis en charact_analyse.py).
//  - Por evento: t_utc (s, TTT ajustado por run) y, por canal 0-63, amplitud en la ventana signal [49,58) y en la de
//    control [40,49) (línea base recortada de [0,40) − mínimo, mV), integral en signal Σ(línea base − x)·16 ns
//    (mV·ns) y muestra del mínimo de todo el registro. Binario: registros {double t; float s[64], c[64], q[64];
//    uint8 k[64]} (840 bytes).
//  - Por hora UTC y canal: pulsos fuera de tiempo (buscador del DAQ, separación 4) en pedestal + control + post
//    (91 muestras = 1.456 µs por evento; sin la cola tail), sobre 0.5 / 1.5 / 3 / 5 p.e. (hod2018, T·G de la
//    calibración) o 1 / 3 mV (hod2019); exposición; nivel de pedestal (suma y suma de cuadrados) y RMS recortado
//    del pedestal (muestras a ±2 mV de la mediana) por evento.
//  - Por canal: espectro de amplitud fuera de tiempo (0.25 mV hasta 100 mV) e integral media de los pulsos fuera de
//    tiempo de 1 p.e. (0.5-1.5 G, ventana [k−2, k+7), 9 muestras como signal) para pasar la integral a p.e.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/characterization/charact_extract.C+("<run>_seg*.root", "<calibración.json>", "<prefijo de salida>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/pulse_finder.h"

namespace {
#pragma pack(push, 1)
struct Rec {
  double t;
  float s[64], c[64], q[64];
  uint8_t k[64];
};
#pragma pack(pop)
struct Hour {
  double exposure[64] = {0};
  long n[64] = {0};
  long cnt[64][4] = {{0}};
  double pedSum[64] = {0}, pedSq[64] = {0}, rmsSum[64] = {0};
};
}  // namespace

void charact_extract(const char* pattern, const char* calibJson, const char* prefix) {
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
    if (!std::isfinite(G[c])) G[c] = ok[ok.size() / 2];  // canales sin G (ch3 muerto, ch22): mediana
  double thr[64][4];
  for (int c = 0; c < 64; ++c) {
    const double h18[4] = {0.5, 1.5, 3, 5}, h19[4] = {1, 3, 1e9, 1e9};
    for (int i = 0; i < 4; ++i) thr[c][i] = c < 32 ? h18[i] * G[c] : h19[i];
  }
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
  const double lsb = daq::pulse::lsbMv();
  const std::vector<daq::pulse::Region> oot{{0, 40}, {40, 49}, {102, n}};
  const double ootS = (40 + 9 + (n - 102)) * 16e-9;
  const double binMv = 0.25;
  const int nb = static_cast<int>(100.0 / binMv);
  std::vector<std::vector<long>> spec(64, std::vector<long>(static_cast<std::size_t>(nb) + 1, 0));
  std::vector<double> q1Sum(64, 0);
  std::vector<long> q1N(64, 0);
  std::map<long, Hour> hours;
  FILE* fb = fopen((std::string(prefix) + "_events.bin").c_str(), "wb");
  std::vector<daq::pulse::Pulse> p;
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    Rec r{};
    r.t = static_cast<double>(tutc) * 1e-9;
    Hour& h = hours[static_cast<long>(std::floor(r.t / 3600.0))];
    for (int c = 0; c < 64; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      const double med = daq::pulse::median(x, daq::pulse::Region{0, 40});
      double s2 = 0;
      int ns = 0;
      for (int i = 0; i < 40; ++i)
        if (std::fabs(x[i] - med) * lsb <= 2.0) {
          s2 += (x[i] - bl) * (x[i] - bl);
          ++ns;
        }
      UShort_t ms = x[49], mc = x[40];
      double q = 0;
      for (int i = 49; i < 58; ++i) {
        ms = std::min(ms, x[i]);
        q += (bl - x[i]) * lsb * 16.0;
      }
      for (int i = 40; i < 49; ++i) mc = std::min(mc, x[i]);
      r.s[c] = static_cast<float>((bl - ms) * lsb);
      r.c[c] = static_cast<float>((bl - mc) * lsb);
      r.q[c] = static_cast<float>(q);
      r.k[c] = static_cast<uint8_t>(std::min_element(x, x + n) - x);
      p.clear();
      daq::pulse::findPulses(x, n, bl, lsb, 0.9, 4, oot, p);
      for (const auto& u : p) {
        for (int i = 0; i < 4; ++i)
          if (u.amplitudeMv > thr[c][i]) ++h.cnt[c][i];
        ++spec[c][static_cast<std::size_t>(std::min(nb, static_cast<int>(u.amplitudeMv / binMv)))];
        if (c < 32 && u.amplitudeMv > 0.5 * G[c] && u.amplitudeMv < 1.5 * G[c] && u.sample >= 2 && u.sample + 7 <= n) {
          double qq = 0;
          for (int i = u.sample - 2; i < u.sample + 7; ++i) qq += (bl - x[i]) * lsb * 16.0;
          q1Sum[c] += qq;
          ++q1N[c];
        }
      }
      h.exposure[c] += ootS;
      ++h.n[c];
      const double blMv = bl * lsb;
      h.pedSum[c] += blMv;
      h.pedSq[c] += blMv * blMv;
      h.rmsSum[c] += ns > 1 ? std::sqrt(s2 / ns) * lsb : 0;
    }
    fwrite(&r, sizeof(Rec), 1, fb);
  }
  fclose(fb);
  FILE* f = fopen((std::string(prefix) + "_hourly.csv").c_str(), "w");
  fprintf(f, "hour_start_utc_s,ch,events,exposure_s,cnt0,cnt1,cnt2,cnt3,thr0_mV,thr1_mV,thr2_mV,thr3_mV,ped_sum_mV,ped_sumsq_mV2,rms_sum_mV\n");
  for (const auto& [hh, h] : hours)
    for (int c = 0; c < 64; ++c)
      fprintf(f, "%ld,%d,%ld,%.6f,%ld,%ld,%ld,%ld,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f\n", hh * 3600, c, h.n[c], h.exposure[c], h.cnt[c][0], h.cnt[c][1],
              h.cnt[c][2], h.cnt[c][3], thr[c][0], thr[c][1], thr[c][2] > 1e8 ? -1.0 : thr[c][2], thr[c][3] > 1e8 ? -1.0 : thr[c][3], h.pedSum[c],
              h.pedSq[c], h.rmsSum[c]);
  fclose(f);
  f = fopen((std::string(prefix) + "_ootspec.csv").c_str(), "w");
  fprintf(f, "ch,bin_lo_mV,counts\n");
  for (int c = 0; c < 64; ++c)
    for (int b = 0; b <= nb; ++b)
      if (spec[c][static_cast<std::size_t>(b)]) fprintf(f, "%d,%.2f,%ld\n", c, b * binMv, spec[c][static_cast<std::size_t>(b)]);
  fclose(f);
  f = fopen((std::string(prefix) + "_q1.csv").c_str(), "w");
  fprintf(f, "ch,G_mV,q1_mean_mVns,n\n");
  for (int c = 0; c < 32; ++c) fprintf(f, "%d,%.5f,%.4f,%ld\n", c, G[c], q1N[c] ? q1Sum[c] / q1N[c] : NAN, q1N[c]);
  fclose(f);
  printf("eventos %lld, muestras %d, registro %zu bytes, horas %zu\n", ne, n, sizeof(Rec), hours.size());
}
