// dark_spectrum_long.C — seguimiento 2026-10-31 (6), punto 3 (solo lectura): espectro de amplitud de los pulsos oscuros de los 64 canales en
// la región pre-trigger larga de mayo v8, para el ajuste de G de hod2019 (y de hod2018 como control contra el reajuste de PROC).
// Pulsos: listas de la ntupla XPLUS-4b (ntuple_extract.C: buscador del DAQ, separación 3 muestras, piso 1 mV, pedestal local = media
// recortada ±2 mV de la propia región), solo los de la región [lo, hi) en muestras del registro (mayo: [20, 258), 238 muestras).
// Eventos: los de tiempo bueno (archivo de intervalos UTC excluidos, "t0 t1" por línea; vacío = todos).
// Salida: <prefijo>_spec.csv (ch, bin_lo_mV, n) en bins de 0.25 mV de 0 a 100 mV (el último bin acumula el desborde) y
//         <prefijo>_spec_meta.json (eventos, región, muestras por evento).
// Uso: root -l -b -q 'dark_spectrum_long.C+("<ntupla.root>", "<excluidos.txt o vacío>", lo, hi, "<prefijo>")'
#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

void dark_spectrum_long(const char* ntuple, const char* excl, int lo, int hi, const char* prefix) {
  std::vector<std::pair<double, double>> ex;
  if (excl && *excl) {
    std::ifstream f(excl);
    double a, b;
    while (f >> a >> b) ex.push_back({a, b});
  }
  constexpr int kNb = 400;
  constexpr double kBin = 0.25;
  std::vector<std::vector<long>> h(64, std::vector<long>(kNb + 1, 0));
  TFile f(ntuple, "READ");
  TTree* t = (TTree*)f.Get("ev");
  double tu;
  std::vector<UChar_t>* pch = nullptr;
  std::vector<UShort_t>* psmp = nullptr;
  std::vector<float>* pamp = nullptr;
  t->SetBranchStatus("*", 0);
  for (const char* b : {"t_utc_s", "pch", "psmp", "pamp"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("t_utc_s", &tu);
  t->SetBranchAddress("pch", &pch);
  t->SetBranchAddress("psmp", &psmp);
  t->SetBranchAddress("pamp", &pamp);
  long used = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetBranch("t_utc_s")->GetEntry(i);
    bool good = true;
    for (auto& [a, b] : ex)
      if (tu >= a && tu < b) good = false;
    if (!good) continue;
    t->GetEntry(i);
    ++used;
    for (size_t q = 0; q < pch->size(); ++q) {
      const int s = (*psmp)[q];
      if (s < lo || s >= hi) continue;
      const int b = std::min(kNb, static_cast<int>((*pamp)[q] / kBin));
      if (b >= 0) ++h[(*pch)[q]][b];
    }
  }
  FILE* fo = fopen((std::string(prefix) + "_spec.csv").c_str(), "w");
  fprintf(fo, "ch,bin_lo_mV,n\n");
  for (int c = 0; c < 64; ++c)
    for (int b = 0; b <= kNb; ++b) fprintf(fo, "%d,%.2f,%ld\n", c, b * kBin, h[c][b]);
  fclose(fo);
  FILE* fm = fopen((std::string(prefix) + "_spec_meta.json").c_str(), "w");
  fprintf(fm, "{\"ntuple\": \"%s\", \"excluded_file\": \"%s\", \"events\": %ld, \"region_samples\": [%d, %d], \"samples_per_event\": %d, \"bin_mV\": %.2f}\n",
          ntuple, excl, used, lo, hi, hi - lo, kBin);
  fclose(fm);
  printf("eventos usados %ld, región [%d, %d)\n", used, lo, hi);
}
