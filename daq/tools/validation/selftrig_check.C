// selftrig_check.C — DAQ2-13, parte B.3 (solo lectura). Pulsos de hod2018 en runs de auto-trigger (o de
// trigger externo, para comparar): por evento y canal 0-31 con amplitud > minPe, la amplitud máxima del
// registro (línea base recortada de [0,40) − mínimo) en mV y p.e. (G de la calibración; sin G: mediana) y la
// muestra del mínimo. El análisis (ventana en tiempo alrededor de la moda, X e Y en tiempo, mapa X × Y,
// espectros) se hace sobre este CSV.
// Salida: CSV ev,ch,amp_mV,amp_pe,kmin.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/validation/selftrig_check.C+("<run>_seg*.root", "<calibración.json>", "<salida.csv>", 3.0)'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/pulse_finder.h"

void selftrig_check(const char* pattern, const char* calibJson, const char* out, double minPe = 3.0) {
  std::vector<double> G(64, NAN);
  {
    std::ifstream in(calibJson);
    auto j = nlohmann::json::parse(in);
    for (const auto& e : j["channels"])
      if (e.contains("G_mV") && e["G_mV"].is_number()) G[e["ch"].get<int>()] = e["G_mV"].get<double>();
  }
  std::vector<double> ok;
  for (int c = 0; c < 32; ++c)
    if (std::isfinite(G[c])) ok.push_back(G[c]);
  std::sort(ok.begin(), ok.end());
  for (int c = 0; c < 32; ++c)
    if (!std::isfinite(G[c])) G[c] = ok[ok.size() / 2];
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchAddress("adc", adc.data());
  const double lsb = daq::pulse::lsbMv();
  FILE* f = fopen(out, "w");
  fprintf(f, "ev,ch,amp_mV,amp_pe,kmin\n");
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    for (int c = 0; c < 32; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      int k = 0;
      for (int i = 1; i < n; ++i)
        if (x[i] < x[k]) k = i;
      const double a = (bl - x[k]) * lsb;
      if (a / G[c] > minPe) fprintf(f, "%lld,%d,%.2f,%.3f,%d\n", e, c, a, a / G[c], k);
    }
  }
  fclose(f);
  printf("eventos %lld, muestras %d\n", ne, n);
}
