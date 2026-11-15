// time_extract2.C — seguimiento 2026-11-15 (ventana CFD), solo lectura de los crudos (READ): time_extract.C del seguimiento 13 extendido a la
// ventana fuera de tiempo y al formato DAQ v2.
// Formatos: "may" (árbol data_dig, Ch0..Ch31 vector<float> en mV, timestamp en µs; tramo desde la muestra 209, ventana signal en 258) y
// "daqv2" (árbol events, adc UShort_t[64·n] · lsb, t_utc_ns; tramo desde la muestra 0, ventana signal en 49), como hitdef_extract.C.
// Tramo de 144 muestras; pedestal del tramo bl = media recortada ±2 mV de la mediana de [0,40); pulse_finder.h (umbral 0.5 mV, separación 4),
// una llamada por ventana: signal [49,58) y fuera de tiempo [99,108) del tramo (50 muestras = 800 ns después), como hitdef_extract.C.
// Por canal de hod2018 (0-31) y ventana, el pulso más alto: amplitud (mV; aS idéntica a time_extract.C y a hitdef_extract.C, aC idéntica a
// hitdef_extract.C), k (muestra del mínimo en el tramo), p_loc = media de [k−12, k−4] y las 10 muestras x[k−7 … k+2] (mV, sin pedestal).
// Sin pulso: amplitud 0, k = 0xFFFF, muestras 0.
// Salida: <prefijo>_time2.bin, registros {double t_utc_s; float bl[32]; float aS[32]; uint16 kS[32]; float plS[32]; float wS[32][10];
//         float aC[32]; uint16 kC[32]; float plC[32]; float wC[32][10]} (3336 bytes).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/cfdwin_2026_11_15/time_extract2.C+("<archivo(s)>", "may"|"daqv2", "<prefijo>")'
#include <TChain.h>
#include <TLeaf.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
struct Win {
  float a[32];
  uint16_t k[32];
  float pl[32];
  float w[32][10];
};
void fillWin(const std::vector<double>& x, double bl, const daq::pulse::Region& r, int c, Win& o, std::vector<daq::pulse::Pulse>& p) {
  constexpr int kNw = 10, kPre = 7, kSpan = 144;
  o.a[c] = 0;
  o.k[c] = 0xFFFF;
  o.pl[c] = 0;
  for (int i = 0; i < kNw; ++i) o.w[c][i] = 0;
  p.clear();
  daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, std::vector<daq::pulse::Region>{r}, p);
  int kb = -1;
  for (const auto& u : p)
    if (u.amplitudeMv > o.a[c]) o.a[c] = static_cast<float>(u.amplitudeMv), kb = u.sample;
  if (kb < 0) return;
  o.k[c] = static_cast<uint16_t>(kb);
  double s = 0;
  for (int i = kb - 12; i <= kb - 4; ++i) s += x[static_cast<std::size_t>(i)];
  o.pl[c] = static_cast<float>(s / 9.0);
  for (int i = 0; i < kNw; ++i) {
    const int j = kb - kPre + i;
    o.w[c][i] = (j >= 0 && j < kSpan) ? static_cast<float>(x[static_cast<std::size_t>(j)]) : 0.f;
  }
}
}  // namespace

void time_extract2(const char* pattern, const char* format, const char* prefix) {
  constexpr int kSpan = 144;
  const bool may = std::string(format) == "may";
  const int off = may ? 258 - 49 : 0;
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
  } else {
    t.GetEntry(0);
    n = t.GetLeaf("adc")->GetLenStatic() / 64;
    adc.resize(static_cast<std::size_t>(64 * n));
    t.SetBranchStatus("adc", 1);
    t.SetBranchStatus("t_utc_ns", 1);
    t.SetBranchAddress("adc", adc.data());
    t.SetBranchAddress("t_utc_ns", &tutc);
  }
  t.SetCacheSize(256LL << 20);
  const double bootRefEpoch = 1778886208.8;  // como hitdef_extract.C
  const daq::pulse::Region sig{49, 58}, ctl{99, 108};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  FILE* fb = fopen((std::string(prefix) + "_time2.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  Long64_t nw = 0;
  Win S{}, C{};
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double tt = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    float bl[32] = {};
    bool okEv = true;
    for (int c = 0; c < 32; ++c) {
      if (may) {
        if (!cw[c] || static_cast<int>(cw[c]->size()) < off + kSpan) {
          okEv = false;
          break;
        }
        for (int i = 0; i < kSpan; ++i) x[static_cast<std::size_t>(i)] = (*cw[c])[static_cast<std::size_t>(off + i)];
      } else {
        if (n < kSpan) {
          okEv = false;
          break;
        }
        for (int i = 0; i < kSpan; ++i) x[static_cast<std::size_t>(i)] = adc[static_cast<std::size_t>(c * n + i)] * lsb;
      }
      bl[c] = static_cast<float>(daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0));
      fillWin(x, bl[c], sig, c, S, p);
      fillWin(x, bl[c], ctl, c, C, p);
    }
    if (!okEv) continue;
    fwrite(&tt, sizeof(double), 1, fb);
    fwrite(bl, sizeof(float), 32, fb);
    for (Win* wv : {&S, &C}) {
      fwrite(wv->a, sizeof(float), 32, fb);
      fwrite(wv->k, sizeof(uint16_t), 32, fb);
      fwrite(wv->pl, sizeof(float), 32, fb);
      fwrite(wv->w, sizeof(float), 32 * 10, fb);
    }
    ++nw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  fclose(fb);
  printf("eventos %lld, escritos %lld\n", ne, nw);
}
