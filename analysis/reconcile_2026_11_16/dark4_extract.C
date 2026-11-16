// dark4_extract.C — seguimiento 2026-11-16 (reconciliación de razones), solo lectura de los crudos (READ): cuentas de pulsos oscuros por evento y
// canal de hod2018 en la región pre-trigger [16,37) del tramo, con la definición de dark3_extract.C (tramo de 144 muestras alineado con la ventana
// signal, pedestal del tramo = media recortada ±2 mV de [0,40), pulse_finder.h con umbral 0.5 mV y separación 4), sin aislamiento: pulsos con el
// mínimo en [16,37) y amplitud > 0.5 G y > 1.5 G (G = la del análisis, pasada como argumento).
// Salida: <prefijo>_dark4.bin, registros {double t_utc_s; uint8 n05[32]; uint8 n15[32]} (72 bytes) por evento.
// Uso: root -l -b -q -e 'gSystem->SetBuildDir("<dir del repositorio>",true); gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/reconcile_2026_11_16/dark4_extract.C+("<archivo(s)>", "may"|"daqv2", "<G: 32 valores con coma>", "<prefijo>")'
#include <TChain.h>
#include <TLeaf.h>

#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

void dark4_extract(const char* pattern, const char* format, const char* gains, const char* prefix) {
  constexpr int kSpan = 144;
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
  const std::vector<daq::pulse::Region> all{{0, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan);
  FILE* fb = fopen((std::string(prefix) + "_dark4.bin").c_str(), "wb");
  const Long64_t ne = t.GetEntries();
  Long64_t nw = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double tt = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    uint8_t n05[32] = {}, n15[32] = {};
    bool okEv = true;
    for (int c = 0; c < 32; ++c) {
      if (may) {
        if (!cw[c] || static_cast<int>(cw[c]->size()) < off + kSpan) {
          okEv = false;
          break;
        }
        for (int i = 0; i < kSpan; ++i) x[static_cast<std::size_t>(i)] = (*cw[c])[static_cast<std::size_t>(off + i)];
      } else {
        for (int i = 0; i < kSpan; ++i) x[static_cast<std::size_t>(i)] = adc[static_cast<std::size_t>(c * n + i)] * lsb;
      }
      const double bl = daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0);
      p.clear();
      daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, all, p);
      const double g = G[static_cast<std::size_t>(c)];
      for (const auto& u : p) {
        if (u.sample < 16 || u.sample >= 37) continue;
        if (u.amplitudeMv > 0.5 * g && n05[c] < 255) ++n05[c];
        if (u.amplitudeMv > 1.5 * g && n15[c] < 255) ++n15[c];
      }
    }
    if (!okEv) continue;
    fwrite(&tt, sizeof(double), 1, fb);
    fwrite(n05, 1, 32, fb);
    fwrite(n15, 1, 32, fb);
    ++nw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  fclose(fb);
  printf("eventos %lld, escritos %lld\n", ne, nw);
}
