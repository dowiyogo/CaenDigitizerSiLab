// ev_dump.C — seguimiento 2026-11-03 (12), solo lectura de las ntuplas XPLUS-4b (ntuple_extract.C, /data/daq/derived/night_2026-10-25/ntuples/).
// Una pasada por el árbol "ev":
//  (1) <prefijo>_ev.bin: por evento, t_utc_s (f8), aS[64] (f4, mV), kS[64] (u2, muestra del pulso de aS en coordenadas del registro), aC[64] (f4);
//      aS/aC = mayor amplitud en la ventana signal / fuera de tiempo (definición de hitdef_extract.C).
//  (2) <prefijo>_dark.bin: por hora UTC (inicio = floor(t/3600)·3600, como el árbol "spec") y canal, el espectro de los pulsos oscuros de la región
//      pre-trigger [lo, hi) (buscador del DAQ, separación 3, piso 1 mV, pedestal local; los pulsos ya guardados en pch/psmp/pamp) en bins de
//      0.25 mV de 0 a 100 mV (bin 400 = desborde) y la exposición en muestras Σ (hi − lo − 3·pulsos del canal en la región), como dark_stage2.C.
//      Registro: hour i8, ch i4, pad i4, events i8, expo_samples f8, counts u4[401].
// Uso: root -l -b -q 'analysis/followup12_2026_11_03/ev_dump.C+("<ntupla.root>", lo, hi, "<prefijo>")'   (mayo [20, 258), octubre [20, 49))
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

void ev_dump(const char* ntuple, int lo, int hi, const char* prefix) {
  constexpr int kNch = 64, kNb = 400;
  constexpr double kBin = 0.25;
  TFile f(ntuple, "READ");
  auto* t = static_cast<TTree*>(f.Get("ev"));
  double tu;
  float aS[kNch], aC[kNch];
  UShort_t kS[kNch];
  std::vector<UChar_t>* pch = nullptr;
  std::vector<UShort_t>* psmp = nullptr;
  std::vector<float>* pamp = nullptr;
  t->SetBranchStatus("*", 0);
  for (const char* b : {"t_utc_s", "aS", "kS", "aC", "pch", "psmp", "pamp"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("t_utc_s", &tu);
  t->SetBranchAddress("aS", aS);
  t->SetBranchAddress("kS", kS);
  t->SetBranchAddress("aC", aC);
  t->SetBranchAddress("pch", &pch);
  t->SetBranchAddress("psmp", &psmp);
  t->SetBranchAddress("pamp", &pamp);
  struct Acc {
    long events = 0;
    double expo[kNch] = {};
    std::vector<std::vector<uint32_t>> h = std::vector<std::vector<uint32_t>>(kNch, std::vector<uint32_t>(kNb + 1, 0));
  };
  std::map<int64_t, Acc> hours;
  FILE* fe = fopen((std::string(prefix) + "_ev.bin").c_str(), "wb");
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    fwrite(&tu, 8, 1, fe);
    fwrite(aS, 4, kNch, fe);
    fwrite(kS, 2, kNch, fe);
    fwrite(aC, 4, kNch, fe);
    const int64_t hr = static_cast<int64_t>(std::floor(tu / 3600.0)) * 3600;
    auto& a = hours[hr];
    ++a.events;
    int np[kNch] = {};
    for (size_t q = 0; q < pch->size(); ++q) {
      const int s = (*psmp)[q];
      if (s < lo || s >= hi) continue;
      const int c = (*pch)[q];
      ++np[c];
      const int b = std::min(kNb, static_cast<int>((*pamp)[q] / kBin));
      if (b >= 0) ++a.h[c][b];
    }
    for (int c = 0; c < kNch; ++c) a.expo[c] += (hi - lo) - 3.0 * np[c];
  }
  fclose(fe);
  FILE* fd = fopen((std::string(prefix) + "_dark.bin").c_str(), "wb");
  for (auto& [hr, a] : hours)
    for (int c = 0; c < kNch; ++c) {
      const int32_t pad = 0;
      const int32_t ch = c;
      fwrite(&hr, 8, 1, fd);
      fwrite(&ch, 4, 1, fd);
      fwrite(&pad, 4, 1, fd);
      int64_t ev = a.events;
      fwrite(&ev, 8, 1, fd);
      fwrite(&a.expo[c], 8, 1, fd);
      fwrite(a.h[c].data(), 4, kNb + 1, fd);
    }
  fclose(fd);
  printf("%s: eventos %lld, horas %zu, región [%d, %d)\n", ntuple, n, hours.size(), lo, hi);
}
