// hitdef_extract.C — DAQ2-15 (solo lectura). Impactos de hod2018 (ch 0-31) con una definición única para mayo
// (árbol data_dig: Ch0..Ch31 vector<float> en mV, time, timestamp) y octubre (DAQ v2: árbol events, adc).
// Todo se calcula en un tramo de 144 muestras (16 ns) alineado con el inicio de la ventana signal s0:
// octubre s0 = 49 (tramo = registro completo), mayo s0 = 258 (4128 ns, la ventana W2 del análisis de mayo; el
// registro de octubre empieza en 3344 ns del de mayo). En coordenadas del tramo:
//   pedestal [0,40) (media recortada ±2 mV de la mediana), signal [49,58), fuera de tiempo [99,108).
//  - Definición nueva (impacto = pulso): el buscador de pulsos del DAQ (pulse_finder.h: mínimo local, separación
//    4 muestras) con el mínimo dentro de la ventana; amplitud = pedestal − muestra del mínimo (mV). Se guarda la
//    mayor de la ventana signal (aS) y de la de control (aC).
//  - Definición de DAQ2-13 (selftrig_check.C): mínimo de todo el tramo; amplitud aG y muestra kG.
//  - Definición del análisis de mayo (hod2018_good_time_step2.C): mínimo de la ventana de 10 muestras [49,59) (4128-4272 ns,
//    ambos extremos) con pedestal local = media de las muestras [k−12, k−4]; amplitud aM; la misma en [99,109) (aMc).
//  - Comprobación de la calibración con este pedestal: espectro de los pulsos fuera de tiempo ([0,49) y [102,144))
//    por canal y ajuste de picos con gain_fit.h (semilla = G de la calibración): G' y posiciones de los picos.
// PROC-1: nch = 64 extrae también hod2019 (ch 32-63, misma ventana: en mayo v8 las dos estaciones tienen el máximo de
// tiempos en las muestras 256-263); los registros llevan entonces nch canales por campo (32 = formato anterior, idéntico).
// Salida: <prefijo>_hits.bin, registros {double t_utc_s; float aS[nch], aC[nch], aG[nch]; uint8 kG[nch]; float aM[nch], aMc[nch]}
// (680 bytes con nch = 32, 1352 con nch = 64),
//         <prefijo>_peaks.csv (ch,G_seed,G_fit,status,peaks,positions_mV).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/characterization/hitdef_extract.C+("<archivo(s)>", "may"|"daqv2", "<G semilla: 32 valores con coma>", "<prefijo>"[, nch])'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "daq/gain_fit.h"
#include "daq/pulse_finder.h"

namespace {
#pragma pack(push, 1)
struct Hit {
  double t;
  float aS[64], aC[64], aG[64];
  uint8_t kG[64];
  float aM[64], aMc[64];
};
#pragma pack(pop)
constexpr int kSpan = 144;
}  // namespace

void hitdef_extract(const char* pattern, const char* format, const char* seeds, const char* prefix, int nch = 32) {
  nch = nch == 64 ? 64 : 32;
  std::vector<double> seed;
  {
    std::stringstream ss(seeds);
    std::string tok;
    while (std::getline(ss, tok, ',')) seed.push_back(std::stod(tok));
  }
  const bool may = std::string(format) == "may";
  const int s0 = may ? 258 : 49;
  const int off = s0 - 49;  // primera muestra del tramo en el registro
  const double lsb = daq::pulse::lsbMv();
  TChain t(may ? "data_dig" : "events");
  t.Add(pattern);
  t.SetBranchStatus("*", 0);
  std::vector<float>* cw[64] = {};
  std::vector<UShort_t> adc;
  Long_t ts = 0;
  Long64_t tutc = 0;
  int n = 0;
  if (may) {
    for (int c = 0; c < nch; ++c) {
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
  const double bootRefEpoch = 1778886208.8;  // mayo: época UTC de 2026-05-15T19:03:28.8-04:00 (hod2018_good_time_step2.C)
  const double binMv = 0.25;
  const int nb = static_cast<int>(100.0 / binMv);
  std::vector<std::vector<unsigned>> spec(32, std::vector<unsigned>(static_cast<std::size_t>(nb) + 1, 0));
  const std::vector<daq::pulse::Region> sig{{49, 58}}, ctl{{99, 108}}, oot{{0, 49}, {102, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  FILE* fb = fopen((std::string(prefix) + "_hits.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  Long64_t nw = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    Hit h{};
    h.t = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    bool okEv = true;
    for (int c = 0; c < nch; ++c) {
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
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, sig, p);
      for (const auto& u : p) h.aS[c] = std::max(h.aS[c], static_cast<float>(u.amplitudeMv));
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, ctl, p);
      for (const auto& u : p) h.aC[c] = std::max(h.aC[c], static_cast<float>(u.amplitudeMv));
      const int k = static_cast<int>(std::min_element(x.begin(), x.end()) - x.begin());
      h.aG[c] = static_cast<float>(bl - x[k]);
      h.kG[c] = static_cast<uint8_t>(k);
      for (int w = 0; w < 2; ++w) {
        const int a = w == 0 ? 49 : 99;
        int km = a;
        for (int i = a; i < a + 10; ++i)
          if (x[i] < x[km]) km = i;
        double ped = 0;
        for (int j = -12; j <= -4; ++j) ped += x[km + j];
        (w == 0 ? h.aM : h.aMc)[c] = static_cast<float>(ped / 9.0 - x[km]);
      }
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.9, 4, oot, p);
      if (c < 32)  // espectro para la comprobación de G: solo hod2018
        for (const auto& u : p) ++spec[static_cast<std::size_t>(c)][static_cast<std::size_t>(std::min(nb, static_cast<int>(u.amplitudeMv / binMv)))];
    }
    if (!okEv) continue;
    // campo por campo con nch canales (con nch = 32, los mismos bytes que el registro empaquetado de antes)
    const auto n_ = static_cast<std::size_t>(nch);
    fwrite(&h.t, sizeof(double), 1, fb);
    fwrite(h.aS, sizeof(float), n_, fb);
    fwrite(h.aC, sizeof(float), n_, fb);
    fwrite(h.aG, sizeof(float), n_, fb);
    fwrite(h.kG, sizeof(uint8_t), n_, fb);
    fwrite(h.aM, sizeof(float), n_, fb);
    fwrite(h.aMc, sizeof(float), n_, fb);
    ++nw;
  }
  fclose(fb);
  FILE* f = fopen((std::string(prefix) + "_peaks.csv").c_str(), "w");
  fprintf(f, "ch,G_seed,G_fit,status,peaks,positions_mV\n");
  for (int c = 0; c < 32; ++c) {
    const double sd = c < static_cast<int>(seed.size()) && seed[static_cast<std::size_t>(c)] > 0 ? seed[static_cast<std::size_t>(c)] : 4.0;
    const auto r = daq::gain::fitGain(spec[static_cast<std::size_t>(c)], binMv, sd);
    std::string pos;
    for (double v : r.positionsMv) pos += (pos.empty() ? "" : " ") + std::to_string(v);
    fprintf(f, "%d,%.5f,%.5f,%s,%d,%s\n", c, sd, r.gainMv, r.status.c_str(), r.peaks, pos.c_str());
  }
  fclose(f);
  printf("eventos %lld, escritos %lld, registro %d bytes, canales %d\n", ne, nw, 8 + nch * 21, nch);
}
