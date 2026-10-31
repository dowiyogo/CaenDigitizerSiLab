// h19_live_check.C — seguimiento 2026-10-31 (5), punto 1: comprobación rápida (solo lectura) de que la lectura de hod2019 está viva, para el
// checklist del 2026-11-02. Sobre los primeros nev eventos de un segmento cerrado o abierto (formato daqv2, o "may" para la referencia):
// por canal, la tasa de cruces descendentes de −thr mV respecto del pedestal del evento (media recortada ±2 mV de las 40 primeras muestras
// del tramo), con rearme al volver por encima de −thr/3, en las muestras [0,40) del tramo (pre-trigger), en kHz. No usa el buscador de pulsos.
// Salida: una línea por canal y la mediana de hod2018 (0-31) y de hod2019 (32-63).
// Uso: root -l -b -q 'h19_live_check.C+("<segmento.root>", "daqv2"|"may", nev, thr_mV)'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

void h19_live_check(const char* raw, const char* format = "daqv2", int nev = 500, double thr = 1.5) {
  const bool may = std::string(format) == "may";
  const int off = may ? 209 : 0;
  TChain ch(may ? "data_dig" : "events");
  ch.Add(raw);
  std::vector<std::vector<float>*> cw(64, nullptr);
  std::vector<UShort_t> adc;
  int len = 0;
  ch.SetBranchStatus("*", 0);
  if (may) {
    for (int c = 0; c < 64; ++c) {
      ch.SetBranchStatus(Form("Ch%d", c), 1);
      ch.SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
  } else {
    len = ch.GetLeaf("adc")->GetLenStatic() / 64;
    adc.resize(64 * len);
    ch.SetBranchStatus("adc", 1);
    ch.SetBranchAddress("adc", adc.data());
  }
  const double lsb = 2000.0 / 4095.0;
  std::vector<long> nx(64, 0);
  int used = 0;
  double x[40];
  for (Long64_t e = 0; e < ch.GetEntries() && used < nev; ++e, ++used) {
    ch.GetEntry(e);
    for (int c = 0; c < 64; ++c) {
      for (int s = 0; s < 40; ++s) x[s] = may ? (*cw[c])[off + s] : adc[c * len + s] * lsb;
      std::vector<double> v(x, x + 40);
      std::nth_element(v.begin(), v.begin() + 20, v.end());
      double sum = 0;
      int k = 0;
      for (double y : x)
        if (std::fabs(y - v[20]) <= 2.0) { sum += y; ++k; }
      const double bl = k ? sum / k : v[20];
      bool armed = true;
      for (int s = 0; s < 40; ++s) {
        const double d = x[s] - bl;
        if (armed && d < -thr) { ++nx[c]; armed = false; }
        if (!armed && d > -thr / 3) armed = true;
      }
    }
  }
  const double live_s = used * 40 * 16e-9;
  std::vector<double> r18, r19;
  printf("eventos %d, umbral %.2f mV, tiempo vivo %.3g s\nch,cruces,tasa_kHz\n", used, thr, live_s);
  for (int c = 0; c < 64; ++c) {
    const double r = live_s > 0 ? nx[c] / live_s / 1e3 : NAN;
    printf("%d,%ld,%.3f\n", c, nx[c], r);
    (c < 32 ? r18 : r19).push_back(r);
  }
  auto med = [](std::vector<double> v) { std::nth_element(v.begin(), v.begin() + 16, v.end()); return v[16]; };
  auto above = [](const std::vector<double>& v, double t) { return std::count_if(v.begin(), v.end(), [t](double r) { return r > t; }); };
  printf("RESUMEN hod2018 mediana %.3f kHz | hod2019 mediana %.3f kHz, canales de hod2019 > 50 kHz: %ld de 32\n", med(r18), med(r19), above(r19, 50.0));
}
