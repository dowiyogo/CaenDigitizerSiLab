// dark2_extract.C — seguimiento 2026-11-15 (consistencia de la calibración de amplitud mayo / octubre), solo lectura de los crudos (READ).
// Pulsos oscuros aislados de hod2018 (ch 0-31) con amplitud y carga, para comparar amplitud/carga en 1 y 2 p.e. con la de los muones.
// Definiciones idénticas a charge_extract.C (seguimiento 10): tramo de 144 muestras (16 ns) alineado con la ventana signal (mayo s0 = 258,
// octubre s0 = 49); pedestal del tramo [0,40) = media recortada ±2 mV de la mediana; pulse_finder.h (umbral 0.5 mV, separación 4) sobre todo
// el tramo; pulso oscuro = mínimo en [16,37) o [102,132), sin otro pulso > 0.5 G en [k−16, k+16]; carga = Σ (p_loc − x_i) · 16 ns en
// [k−3, k+12], p_loc = media de [k−12, k−4]. Diferencia: se guardan todos los pulsos aislados con amplitud en [0.3 G, 5.5 G] (no solo
// [0.5 G, 1.5 G]) y uno por registro, para elegir 1 y 2 p.e. después con cualquier juego de G.
// Salida: <prefijo>_dark2.bin, registros {uint32 evento (índice entre los eventos escritos); uint8 ch; uint8 región (0: [16,37), 1: [102,132));
//         uint16 k; float a_mV; float q_mVns} (16 bytes); <prefijo>_dark2_t.bin: double t_utc_s por evento escrito.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/ampcal_2026_11_15/dark2_extract.C+("<archivo(s)>", "may"|"daqv2", "<G: 32 valores con coma>", "<prefijo>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
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
#pragma pack(push, 1)
struct Rec {
  uint32_t ev;
  uint8_t ch, region;
  uint16_t k;
  float a, q;
};
#pragma pack(pop)
static_assert(sizeof(Rec) == 16, "registro de 16 bytes");
}  // namespace

void dark2_extract(const char* pattern, const char* format, const char* gains, const char* prefix) {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  }
  if (G.size() != 32) {
    printf("se esperaban 32 G, hay %zu\n", G.size());
    return;
  }
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
  const std::vector<daq::pulse::Region> all{{0, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  FILE* fb = fopen((std::string(prefix) + "_dark2.bin").c_str(), "wb");
  FILE* ft = fopen((std::string(prefix) + "_dark2_t.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  uint32_t nw = 0;
  Long64_t nrec = 0;
  std::vector<Rec> buf;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double tt = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    bool okEv = true;
    buf.clear();
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
      const double bl = daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0);
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, all, p);
      const double g = G[static_cast<std::size_t>(c)];
      for (std::size_t i = 0; i < p.size(); ++i) {
        const int k = p[i].sample;
        const bool r0 = k >= 16 && k < 37, r1 = k >= 102 && k < 132;
        if (!r0 && !r1) continue;
        const double a = p[i].amplitudeMv;
        if (a < 0.3 * g || a > 5.5 * g) continue;
        bool iso = true;
        for (std::size_t j = 0; j < p.size(); ++j)
          if (j != i && p[j].amplitudeMv > 0.5 * g && p[j].sample >= k - 16 && p[j].sample <= k + 16) iso = false;
        if (!iso) continue;
        buf.push_back(Rec{nw, static_cast<uint8_t>(c), static_cast<uint8_t>(r0 ? 0 : 1), static_cast<uint16_t>(k), static_cast<float>(a),
                          static_cast<float>(charge(x, k))});
      }
    }
    if (!okEv) continue;
    fwrite(buf.data(), sizeof(Rec), buf.size(), fb);
    fwrite(&tt, sizeof(double), 1, ft);
    nrec += static_cast<Long64_t>(buf.size());
    ++nw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  fclose(fb);
  fclose(ft);
  printf("eventos %lld, escritos %u, pulsos %lld\n", ne, nw, nrec);
}
