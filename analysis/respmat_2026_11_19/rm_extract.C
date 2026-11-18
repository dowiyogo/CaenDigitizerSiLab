// rm_extract.C — seguimiento 2026-11-19 (matriz de respuesta), primera pasada sobre los crudos de hod2018 (solo lectura, READ).
// Geometría y estimador idénticos a hitdef_extract.C (DAQ2-15), que produce el aS de los espectros de muones etiquetados:
//   tramo de 144 muestras (16 ns) alineado con la ventana signal [49,58) (mayo: muestras 209-352 del registro de 624; octubre: el registro
//   completo de 144); pedestal bl = media recortada ±2 mV de la mediana de [0,40); pulse_finder.h (umbral 0.5 mV, separación 4) restringido a
//   la región [49,58); aS = mayor amplitud (bl − x en la muestra del mínimo local), 0 si no hay pulso. Sin filtro ni interpolación.
// Por evento escribe aS (para comprobar contra hd64_*_hits.bin) y el tiempo CFD del pulso de aS: 50 % de aS en el flanco de bajada,
// interpolación lineal entre las dos muestras que lo cruzan (en muestras del tramo; −1 si aS = 0).
// Por canal guarda además los pulsos oscuros de la región pre-trigger k ∈ [15,33] (búsqueda en todo el tramo, como iso_extract.C) con
// amplitud ≥ 0.3 G y sin otro pulso > 0.5 G del mismo canal a |Δk| ≤ 9 (±144 ns, la ventana de ±150 ns): amplitud, tiempo CFD y las 18
// muestras y = bl − x en k−8 … k+9.
// Salida: <prefijo>_ev.bin  registros {double t_utc_s; float aS[32]; float tS[32]} (264 bytes) por evento escrito;
//         <prefijo>_iso.bin registros {uint32 ev; uint8 ch; uint8 isod; uint16 k; float a_mV; float tcfd; float y[18]} (88 bytes).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/respmat_2026_11_19/rm_extract.C+("<archivo(s)>", "may"|"daqv2", "<G 32 valores>", "<prefijo>")'
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
#pragma pack(push, 1)
struct Ev {
  double t;
  float aS[32], tS[32];
};
struct Iso {
  uint32_t ev;
  uint8_t ch, isod;
  uint16_t k;
  float a, tcfd;
  float y[18];
};
#pragma pack(pop)
static_assert(sizeof(Ev) == 264, "264 bytes");
static_assert(sizeof(Iso) == 88, "88 bytes");
constexpr int kSpan = 144;
// Tiempo CFD (50 %) de un pulso con mínimo en k y amplitud a; y = bl − x.
double cfd(const std::vector<double>& y, int k, double a) {
  const double f = 0.5 * a;
  int j = k;
  while (j > 0 && k - j < 6 && y[static_cast<std::size_t>(j - 1)] >= f) --j;
  if (j == 0 || k - j >= 6) return j;
  const double y0 = y[static_cast<std::size_t>(j - 1)], y1 = y[static_cast<std::size_t>(j)];
  return (j - 1) + (y1 > y0 ? (f - y0) / (y1 - y0) : 1.0);
}
}  // namespace

void rm_extract(const char* pattern, const char* format, const char* gains, const char* prefix) {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
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
  const std::vector<daq::pulse::Region> sig{{49, 58}}, all{{0, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan), y(kSpan);
  FILE* fe = fopen((std::string(prefix) + "_ev.bin").c_str(), "wb");
  FILE* fi = fopen((std::string(prefix) + "_iso.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  uint32_t nw = 0;
  Long64_t niso = 0;
  std::vector<Iso> pend;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    Ev h{};
    h.t = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    bool okEv = true;
    pend.clear();
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
      for (int i = 0; i < kSpan; ++i) y[static_cast<std::size_t>(i)] = bl - x[static_cast<std::size_t>(i)];
      // estimador de los espectros (hitdef_extract.C): región signal sola
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, sig, p);
      int kS = -1;
      for (const auto& u : p)
        if (u.amplitudeMv > h.aS[c]) h.aS[c] = static_cast<float>(u.amplitudeMv), kS = u.sample;
      h.tS[c] = kS >= 0 ? static_cast<float>(cfd(y, kS, h.aS[c])) : -1.0f;
      // pulsos oscuros aislados de la región pre-trigger (búsqueda en todo el tramo, como iso_extract.C)
      if (c == 3) continue;
      const double g = G[static_cast<std::size_t>(c)];
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, all, p);
      for (std::size_t i = 0; i < p.size(); ++i) {
        const int k = p[i].sample;
        if (k < 15 || k > 33 || p[i].amplitudeMv < 0.3 * g) continue;
        int dmin = 255;
        for (std::size_t j = 0; j < p.size(); ++j)
          if (j != i && p[j].amplitudeMv > 0.5 * g) dmin = std::min(dmin, std::abs(p[j].sample - k));
        if (dmin <= 9) continue;
        Iso r{};
        r.ch = static_cast<uint8_t>(c);
        r.isod = static_cast<uint8_t>(std::min(dmin, 255));
        r.k = static_cast<uint16_t>(k);
        r.a = static_cast<float>(p[i].amplitudeMv);
        r.tcfd = static_cast<float>(cfd(y, k, p[i].amplitudeMv));
        for (int j = 0; j < 18; ++j) r.y[j] = static_cast<float>(y[static_cast<std::size_t>(k - 8 + j)]);
        pend.push_back(r);
      }
    }
    if (!okEv) continue;
    for (auto& r : pend) r.ev = nw;
    fwrite(pend.data(), sizeof(Iso), pend.size(), fi);
    niso += static_cast<Long64_t>(pend.size());
    fwrite(&h, sizeof(Ev), 1, fe);
    ++nw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  fclose(fe);
  fclose(fi);
  printf("eventos %lld, escritos %u, pulsos aislados %lld\n", ne, nw, niso);
}
