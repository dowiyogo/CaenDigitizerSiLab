// charge_extract.C — 2026-11-01 (solo lectura de los crudos): amplitud y carga integrada de hod2018 (ch 0-31) por evento, y calibración de
// 1 p.e. en carga con pulsos oscuros, para los espectros MIP etiquetados (mip_tagged.py).
// Mismo tramo, pedestal y buscador que hitdef_extract.C: tramo de 144 muestras (16 ns) alineado con la ventana signal (mayo s0 = 258,
// octubre s0 = 49); pedestal [0,40) = media recortada ±2 mV de la mediana; pulsos con pulse_finder.h (umbral 0.5 mV, separación 4).
//  - aS / aC: amplitud del pulso más alto con el mínimo en la ventana signal [49,58) / de control [99,108) (idéntica a hitdef_extract.C).
//  - QS / QC: carga del mismo pulso = Σ (p_loc − x_i) · 16 ns para i en [k−3, k+12] (k = muestra del mínimo; −48 a +192 ns), en mV·ns, con
//    pedestal LOCAL p_loc = media de [k−12, k−4] (la convención de aM del análisis de mayo). Es la ventana de Q de hod2018 del análisis de mayo
//    (ANALYSIS_CHANGES_2026-10-17, eedcebc). El pedestal del tramo no sirve para integrar: con la tasa de pulsos oscuros queda ~4 mV por debajo
//    de la línea base de los tramos sin pulsos (comprobado en mayo ch13).
//  - Pulsos oscuros para la calibración: mínimo en [16,37) o [102,132) (pedestal local y ventana dentro del tramo, fuera de la ventana signal),
//    sin otro pulso > 0.5 G en [k−16, k+16] (aislar contra todo pulso > 0.5 mV elige tramos con la línea base desplazada: el pulso medio
//    de 1 p.e. no volvía a cero); amplitud y carga con la misma ventana y pedestal. Se guardan los de amplitud en [0.5 G, 1.5 G] (1 p.e.).
// Salida: <prefijo>_charge.bin, registros {double t_utc_s; float aS[32], aC[32], QS[32], QC[32]} (520 bytes);
//         <prefijo>_dark.csv (ch, G_mV, n_1pe, Q1_median_mVns, Q1_mean_mVns, Q1_rms_mVns, A_mean_mV, Q_over_A_ns).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/followup10_2026_11_01/charge_extract.C+("<archivo(s)>", "may"|"daqv2", "<G: 32 valores con coma>", "<prefijo>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
constexpr int kSpan = 144;
constexpr double kDt = 16.0;  // ns
double charge(const std::vector<double>& x, int k) {
  double ped = 0;
  for (int i = k - 12; i <= k - 4; ++i) ped += x[static_cast<std::size_t>(i)];
  ped /= 9.0;
  double q = 0;
  for (int i = k - 3; i <= std::min(kSpan - 1, k + 12); ++i) q += (ped - x[static_cast<std::size_t>(i)]) * kDt;
  return q;
}
}  // namespace

void charge_extract(const char* pattern, const char* format, const char* gains, const char* prefix) {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  }
  const bool may = std::string(format) == "may";
  const int s0 = may ? 258 : 49;
  const int off = s0 - 49;
  const double lsb = daq::pulse::lsbMv();
  TChain t(may ? "data_dig" : "events");
  t.Add(pattern);
  t.SetBranchStatus("*", 0);
  std::vector<float>* cw[32] = {};
  std::vector<UShort_t> adc;
  Long_t ts = 0;
  Long64_t tutc = 0;
  int n = 0;
  if (may) {
    for (int c = 0; c < 32; ++c) {
      t.SetBranchStatus(Form("Ch%d", c), 1);
      t.SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    t.SetBranchStatus("timestamp", 1);
    t.SetBranchAddress("timestamp", &ts);
    t.SetCacheSize(256LL << 20);
  } else {
    n = t.GetLeaf("adc")->GetLenStatic() / 64;
    adc.resize(64 * n);
    t.SetBranchStatus("adc", 1);
    t.SetBranchStatus("t_utc_ns", 1);
    t.SetBranchAddress("adc", adc.data());
    t.SetBranchAddress("t_utc_ns", &tutc);
  }
  const double bootRefEpoch = 1778886208.8;  // como hitdef_extract.C
  const std::vector<daq::pulse::Region> sig{{49, 58}}, ctl{{99, 108}}, all{{0, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  std::vector<std::vector<float>> dq(32), da(32);
  FILE* fb = fopen((std::string(prefix) + "_charge.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  Long64_t nw = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    double tt = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    float aS[32] = {}, aC[32] = {}, QS[32] = {}, QC[32] = {};
    bool okEv = true;
    for (int c = 0; c < 32; ++c) {
      if (may) {
        if (!cw[c] || static_cast<int>(cw[c]->size()) < off + kSpan) {
          okEv = false;
          break;
        }
        for (int i = 0; i < kSpan; ++i) x[i] = (*cw[c])[static_cast<std::size_t>(off + i)];
      } else {
        for (int i = 0; i < kSpan; ++i) x[i] = adc[static_cast<std::size_t>(c * n + i)] * lsb;
      }
      const double bl = daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0);
      for (int w = 0; w < 2; ++w) {
        p.clear();
        daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, w == 0 ? sig : ctl, p);
        float best = 0;
        int kb = -1;
        for (const auto& u : p)
          if (u.amplitudeMv > best) best = static_cast<float>(u.amplitudeMv), kb = u.sample;
        (w == 0 ? aS : aC)[c] = best;
        (w == 0 ? QS : QC)[c] = kb >= 0 ? static_cast<float>(charge(x, kb)) : 0.f;
      }
      // pulsos oscuros aislados
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, all, p);
      for (std::size_t i = 0; i < p.size(); ++i) {
        const int k = p[i].sample;
        if (!((k >= 16 && k < 37) || (k >= 102 && k < 132))) continue;
        const double a = p[i].amplitudeMv, g = G[static_cast<std::size_t>(c)];
        bool iso = true;
        for (std::size_t j = 0; j < p.size(); ++j)
          if (j != i && p[j].amplitudeMv > 0.5 * g && p[j].sample >= k - 16 && p[j].sample <= k + 16) iso = false;
        if (iso && a >= 0.5 * g && a <= 1.5 * g) {
          dq[static_cast<std::size_t>(c)].push_back(static_cast<float>(charge(x, k)));
          da[static_cast<std::size_t>(c)].push_back(static_cast<float>(a));
        }
      }
    }
    if (!okEv) continue;
    fwrite(&tt, sizeof(double), 1, fb);
    fwrite(aS, sizeof(float), 32, fb);
    fwrite(aC, sizeof(float), 32, fb);
    fwrite(QS, sizeof(float), 32, fb);
    fwrite(QC, sizeof(float), 32, fb);
    ++nw;
  }
  fclose(fb);
  FILE* f = fopen((std::string(prefix) + "_dark.csv").c_str(), "w");
  fprintf(f, "ch,G_mV,n_1pe,Q1_median_mVns,Q1_mean_mVns,Q1_rms_mVns,A_mean_mV,Q_over_A_ns\n");
  for (int c = 0; c < 32; ++c) {
    auto& q = dq[static_cast<std::size_t>(c)];
    auto& a = da[static_cast<std::size_t>(c)];
    if (q.empty()) {
      fprintf(f, "%d,%.5f,0,,,,,\n", c, G[static_cast<std::size_t>(c)]);
      continue;
    }
    double sq = 0, sq2 = 0, sa = 0;
    for (std::size_t i = 0; i < q.size(); ++i) sq += q[i], sq2 += q[i] * q[i], sa += a[i];
    const double m = sq / q.size(), rms = std::sqrt(std::max(0.0, sq2 / q.size() - m * m)), am = sa / a.size();
    std::vector<float> s = q;
    std::nth_element(s.begin(), s.begin() + s.size() / 2, s.end());
    fprintf(f, "%d,%.5f,%zu,%.4f,%.4f,%.4f,%.4f,%.4f\n", c, G[static_cast<std::size_t>(c)], q.size(), s[s.size() / 2], m, rms, am, m / am);
  }
  fclose(f);
  printf("eventos %lld, escritos %lld\n", ne, nw);
}
