// time_extract.C — seguimiento 2026-11-03 (13), solo lectura del crudo de mayo v8 (data_dig, READ): forma del frente de subida de los pulsos de
// hod2018 (ch 0-31) para la resolución temporal (time_res.py) y amplitudes para la posición en el solapamiento (overlap_ratio.py).
// Mismo tramo, pedestal y buscador que charge_extract.C / hitdef_extract.C: tramo de 144 muestras (16 ns) alineado con la ventana signal
// (s0 = 258); pedestal del tramo bl = media recortada ±2 mV de la mediana de [0,40); pulsos con pulse_finder.h (umbral 0.5 mV, separación 4).
// Por canal, el pulso más alto con el mínimo en la ventana signal [49,58): aS (mV, idéntica a hitdef_extract.C), kS (muestra del mínimo en el
// tramo), bl, p_loc = media de [k−12, k−4] (pedestal local de charge_extract.C) y las 10 muestras x[k−7 … k+2] (mV, sin restar pedestal).
// Sin pulso: aS = 0, kS = 0xFFFF, muestras a 0.
// Salida: <prefijo>_time.bin, registros {double t_utc_s; float aS[32]; uint16 kS[32]; float bl[32]; float ploc[32]; float w[32][10]} (1736 bytes).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/followup13_2026_11_03/time_extract.C+("<v8.root>", "<prefijo>")'
#include <TChain.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

void time_extract(const char* pattern, const char* prefix) {
  constexpr int kSpan = 144, kNw = 10, kPre = 7;
  const int off = 258 - 49;
  TChain t("data_dig");
  t.Add(pattern);
  t.SetBranchStatus("*", 0);
  std::vector<float>* cw[32] = {};
  Long_t ts = 0;
  for (int c = 0; c < 32; ++c) {
    t.SetBranchStatus(Form("Ch%d", c), 1);
    t.SetBranchAddress(Form("Ch%d", c), &cw[c]);
  }
  t.SetBranchStatus("timestamp", 1);
  t.SetBranchAddress("timestamp", &ts);
  t.SetCacheSize(256LL << 20);
  const double bootRefEpoch = 1778886208.8;  // como hitdef_extract.C
  const std::vector<daq::pulse::Region> sig{{49, 58}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  FILE* fb = fopen((std::string(prefix) + "_time.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  Long64_t nw = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double tt = bootRefEpoch + static_cast<double>(ts) * 1e-6;
    float aS[32] = {}, bl[32] = {}, pl[32] = {}, w[32][kNw] = {};
    uint16_t kS[32];
    bool okEv = true;
    for (int c = 0; c < 32; ++c) {
      kS[c] = 0xFFFF;
      if (!cw[c] || static_cast<int>(cw[c]->size()) < off + kSpan) {
        okEv = false;
        break;
      }
      for (int i = 0; i < kSpan; ++i) x[i] = (*cw[c])[static_cast<std::size_t>(off + i)];
      bl[c] = static_cast<float>(daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0));
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl[c], 1.0, 0.5, 4, sig, p);
      int kb = -1;
      for (const auto& u : p)
        if (u.amplitudeMv > aS[c]) aS[c] = static_cast<float>(u.amplitudeMv), kb = u.sample;
      if (kb < 0) continue;
      kS[c] = static_cast<uint16_t>(kb);
      double s = 0;
      for (int i = kb - 12; i <= kb - 4; ++i) s += x[static_cast<std::size_t>(i)];
      pl[c] = static_cast<float>(s / 9.0);
      for (int i = 0; i < kNw; ++i) w[c][i] = static_cast<float>(x[static_cast<std::size_t>(kb - kPre + i)]);
    }
    if (!okEv) continue;
    fwrite(&tt, sizeof(double), 1, fb);
    fwrite(aS, sizeof(float), 32, fb);
    fwrite(kS, sizeof(uint16_t), 32, fb);
    fwrite(bl, sizeof(float), 32, fb);
    fwrite(pl, sizeof(float), 32, fb);
    fwrite(w, sizeof(float), 32 * kNw, fb);
    ++nw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  fclose(fb);
  printf("eventos %lld, escritos %lld\n", ne, nw);
}
