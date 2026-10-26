// gain_blocks.C — XPLUS-4b etapas 3 y 5: G por canal y bloque de tiempo con gain_fit.h, sobre el árbol "spec" de la ntupla
// (espectro oscuro por hora y canal con el método del reajuste de hitdef_extract.C: misma G que PROC cuando el bloque es el run entero).
// Bloques: archivo con "etiqueta t0 t1" por línea (UTC s; una hora entra en el bloque si su inicio está en [t0, t1)); varias líneas con la
// misma etiqueta forman un solo bloque (así un bloque de horas de tiempo bueno no incluye las horas excluidas entre ellas).
// Semillas: G por canal (32 valores, coma). Además, por bloque y canal, la tasa de pulsos oscuros sobre 0.5·G (con la G del bloque si
// el ajuste se aceptó, si no con la semilla) por evento, como control de DCR (91 muestras por evento en el método de hitdef).
// Salida CSV: block,ch,events,G_mV,status,peaks,seed_mV,dark_gt0p5G_per_event
// Uso: root -l -b -q 'analysis/night_2026-10-25/gain_blocks.C+("<ntupla.root>", "<bloques.txt>", "<semillas>", "<salida.csv>")'
#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "daq/gain_fit.h"

void gain_blocks(const char* ntuple, const char* blocks, const char* seeds, const char* outCsv) {
  std::vector<double> seed;
  {
    std::stringstream ss(seeds);
    std::string tok;
    while (std::getline(ss, tok, ',')) seed.push_back(std::stod(tok));
  }
  std::vector<std::tuple<std::string, double, double>> bl;
  {
    std::ifstream in(blocks);
    std::string lab;
    double a, b;
    while (in >> lab >> a >> b) bl.emplace_back(lab, a, b);
  }
  TFile f(ntuple);
  auto* t = static_cast<TTree*>(f.Get("spec"));
  Long64_t hr = 0, nev = 0;
  int ch = 0;
  unsigned counts[401];
  t->SetBranchAddress("hour_utc_s", &hr);
  t->SetBranchAddress("ch", &ch);
  t->SetBranchAddress("events", &nev);
  t->SetBranchAddress("counts", counts);
  std::map<std::pair<std::string, int>, std::pair<long, std::vector<double>>> acc;
  for (Long64_t e = 0; e < t->GetEntries(); ++e) {
    t->GetEntry(e);
    if (ch >= 32) continue;
    for (std::size_t k = 0; k < bl.size(); ++k) {
      if (hr < std::get<1>(bl[k]) || hr >= std::get<2>(bl[k])) continue;
      auto& a = acc[{std::get<0>(bl[k]), ch}];
      if (a.second.empty()) a.second.assign(401, 0.0);
      a.first += nev;
      for (int b = 0; b < 401; ++b) a.second[static_cast<std::size_t>(b)] += counts[b];
    }
  }
  FILE* fo = fopen(outCsv, "w");
  fprintf(fo, "block,ch,events,G_mV,status,peaks,seed_mV,dark_gt0p5G_per_event\n");
  for (auto& [key, a] : acc) {
    const int c = key.second;
    const double sd = c < static_cast<int>(seed.size()) && seed[static_cast<std::size_t>(c)] > 0 ? seed[static_cast<std::size_t>(c)] : 4.0;
    const auto r = daq::gain::fitGain(a.second, 0.25, sd);
    const double g = r.accepted() ? r.gainMv : sd;
    double dk = 0;
    for (int b = static_cast<int>(0.5 * g / 0.25); b < 401; ++b) dk += a.second[static_cast<std::size_t>(b)];
    fprintf(fo, "%s,%d,%ld,%.5f,%s,%d,%.5f,%.6f\n", key.first.c_str(), c, a.first, r.gainMv, r.status.c_str(), r.peaks, sd,
            a.first ? dk / a.first : 0.0);
  }
  fclose(fo);
  printf("bloques %zu\n", bl.size());
}
