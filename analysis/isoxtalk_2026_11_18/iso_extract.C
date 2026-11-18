// iso_extract.C — seguimiento 2026-11-18: pulsos oscuros de hod2018 para la diafonía con pulsos aislados. Solo lectura de los crudos (READ).
// Geometría común a mayo y octubre: tramo de 144 muestras (16 ns) alineado con la ventana signal [49,58) (mayo: muestras 209-352 del registro de
// 624; octubre: el registro completo de 144); pedestal del tramo = media recortada ±2 mV de la mediana de [0,40); pulse_finder.h (umbral 0.5 mV,
// separación 4) sobre todo el tramo, como dark2_extract.C y charge_extract.C.
// Se guardan todos los pulsos con amplitud ≥ 0.3 G y mínimo en la región pre-trigger k ∈ [15, 33] (≥ 200 ns antes de la ventana signal, con una
// ventana de ±15 muestras que no la toca) o en la región posterior k ∈ [73, 128] (≥ 200 ns después de la ventana signal y ventana de ±15 muestras
// dentro del tramo). Por pulso:
//   isod = distancia en muestras al pulso más cercano del mismo canal con amplitud > 0.5 G (255 si no hay);
//   flags: bit0 = alguna barra vecina del mismo plano (índice de tira ±1, ±2) con un pulso ≥ 3 G en |Δk| ≤ 2 (32 ns) en el mismo evento;
//          bit1 = alguna barra del plano ortogonal, mismo criterio, mismo evento; bit2, bit3 = lo mismo con el evento siguiente (mezclado, mismas
//          muestras); bit4 = el canal tiene un pulso > 0.5 G en la ventana signal [49,58) de este evento.
// Salida: <prefijo>_iso.bin, registros {uint32 ev; uint8 ch; uint8 region; uint16 k; float a_mV; uint8 isod; uint8 flags; uint16 pad} (16 bytes);
//         <prefijo>_iso_t.bin: double t_utc_s por evento escrito.
// Uso: root -l -b -q -e 'gSystem->SetBuildDir("<dir del repositorio>",true); gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/isoxtalk_2026_11_18/iso_extract.C+("<archivo(s)>", "may"|"daqv2", "<G 32 valores>", "<orden X 16>", "<orden Y 16>", "<prefijo>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
#pragma pack(push, 1)
struct Rec {
  uint32_t ev;
  uint8_t ch, region;
  uint16_t k;
  float a;
  uint8_t isod, flags;
  uint16_t pad;
};
#pragma pack(pop)
static_assert(sizeof(Rec) == 16, "16 bytes");
std::vector<double> split(const char* s) {
  std::vector<double> v;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ',')) v.push_back(std::stod(t));
  return v;
}
struct P {
  int k;
  double a;
};
}  // namespace

void iso_extract(const char* pattern, const char* format, const char* gains, const char* xorder, const char* yorder, const char* prefix) {
  constexpr int kSpan = 144;
  const auto G = split(gains);
  const auto xo = split(xorder), yo = split(yorder);
  int plane[32], idx[32];
  for (int c = 0; c < 32; ++c) plane[c] = -1, idx[c] = -1;
  for (int i = 0; i < 16; ++i) {
    plane[static_cast<int>(xo[i])] = 0, idx[static_cast<int>(xo[i])] = i;
    plane[static_cast<int>(yo[i])] = 1, idx[static_cast<int>(yo[i])] = i;
  }
  int byIdx[2][16];
  for (int i = 0; i < 16; ++i) byIdx[0][i] = static_cast<int>(xo[i]), byIdx[1][i] = static_cast<int>(yo[i]);
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
  const double bootRefEpoch = 1778886208.8;
  const std::vector<daq::pulse::Region> all{{0, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  FILE* fb = fopen((std::string(prefix) + "_iso.bin").c_str(), "wb");
  FILE* ft = fopen((std::string(prefix) + "_iso_t.bin").c_str(), "wb");
  std::vector<std::vector<P>> prev(32), cur(32);
  std::vector<Rec> pend;  // registros del evento anterior, a completar con el evento actual (mezclado)
  auto hasBig = [&](const std::vector<std::vector<P>>& ev, int ch, int k) {
    for (const auto& q : ev[static_cast<std::size_t>(ch)])
      if (std::abs(q.k - k) <= 2 && q.a >= 3.0 * G[static_cast<std::size_t>(ch)]) return true;
    return false;
  };
  auto nbFlags = [&](const std::vector<std::vector<P>>& ev, int c, int k) {
    uint8_t f = 0;
    const int pl = plane[c], i = idx[c];
    for (int d : {-2, -1, 1, 2}) {
      const int j = i + d;
      if (j < 0 || j > 15) continue;
      const int c2 = byIdx[pl][j];
      if (c2 != 3 && hasBig(ev, c2, k)) f |= 1;
    }
    for (int j = 0; j < 16; ++j) {
      const int c2 = byIdx[1 - pl][j];
      if (c2 != 3 && hasBig(ev, c2, k)) f |= 2;
    }
    return f;
  };
  const Long64_t ne = t.GetEntries();
  uint32_t nw = 0;
  Long64_t nrec = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double tt = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    bool okEv = true;
    for (int c = 0; c < 32; ++c) {
      cur[static_cast<std::size_t>(c)].clear();
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
      for (const auto& u : p) cur[static_cast<std::size_t>(c)].push_back({u.sample, u.amplitudeMv});
    }
    if (!okEv) continue;
    // completar los registros del evento anterior con este evento como mezclado
    for (auto& r : pend) {
      const uint8_t fm = nbFlags(cur, r.ch, r.k);
      r.flags |= static_cast<uint8_t>((fm & 1 ? 4 : 0) | (fm & 2 ? 8 : 0));
    }
    fwrite(pend.data(), sizeof(Rec), pend.size(), fb);
    nrec += static_cast<Long64_t>(pend.size());
    pend.clear();
    for (int c = 0; c < 32; ++c) {
      if (c == 3) continue;
      const double g = G[static_cast<std::size_t>(c)];
      const auto& v = cur[static_cast<std::size_t>(c)];
      bool sig = false;
      for (const auto& q : v)
        if (q.k >= 49 && q.k < 58 && q.a > 0.5 * g) sig = true;
      for (std::size_t i = 0; i < v.size(); ++i) {
        const int k = v[i].k;
        const int reg = (k >= 15 && k <= 33) ? 0 : ((k >= 73 && k <= 128) ? 1 : -1);
        if (reg < 0 || v[i].a < 0.3 * g) continue;
        int dmin = 255;
        for (std::size_t j = 0; j < v.size(); ++j)
          if (j != i && v[j].a > 0.5 * g) dmin = std::min(dmin, std::abs(v[j].k - k));
        Rec r{nw, static_cast<uint8_t>(c), static_cast<uint8_t>(reg), static_cast<uint16_t>(k), static_cast<float>(v[i].a),
              static_cast<uint8_t>(std::min(dmin, 255)), 0, 0};
        r.flags = nbFlags(cur, c, k);
        if (sig) r.flags |= 16;
        pend.push_back(r);
      }
    }
    fwrite(&tt, sizeof(double), 1, ft);
    ++nw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  // el último evento no tiene siguiente: bits de mezclado en 0
  fwrite(pend.data(), sizeof(Rec), pend.size(), fb);
  nrec += static_cast<Long64_t>(pend.size());
  fclose(fb);
  fclose(ft);
  printf("eventos %lld, escritos %u, pulsos %lld\n", ne, nw, nrec);
}
