// record_info.C — seguimiento 2026-11-18: longitud de registro de los crudos (solo lectura). Mayo: tamaño de Ch0 en las primeras entradas;
// DAQ v2: longitud estática de la hoja adc / 64 canales. Uso: root -l -b -q 'analysis/isoxtalk_2026_11_18/record_info.C("<archivo>", "may"|"daqv2")'
#include <TChain.h>
#include <TLeaf.h>
#include <cstdio>
#include <string>
#include <vector>
void record_info(const char* f, const char* fmt) {
  const bool may = std::string(fmt) == "may";
  TChain t(may ? "data_dig" : "events");
  t.Add(f);
  if (may) {
    std::vector<float>* c = nullptr;
    t.SetBranchStatus("*", 0);
    t.SetBranchStatus("Ch0", 1);
    t.SetBranchAddress("Ch0", &c);
    for (Long64_t e = 0; e < 3; ++e) {
      t.GetEntry(e);
      printf("entrada %lld: muestras %zu\n", e, c->size());
    }
  } else {
    t.GetEntry(0);
    printf("adc: %d valores = %d muestras x 64 canales\n", t.GetLeaf("adc")->GetLenStatic(), t.GetLeaf("adc")->GetLenStatic() / 64);
  }
  printf("entradas %lld\n", t.GetEntries());
}
