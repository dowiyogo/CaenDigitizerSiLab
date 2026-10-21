// selftrig_noise.C — DAQ2-13, parte A.2 (solo lectura). Tasa de pulsos oscuros de hod2018 contra el umbral,
// desde runs de SW trigger con registro completo (el SW trigger muestrea al azar: todo el registro sirve).
//  - Por canal 0-31 (o 0-63 con nch = 64, DAQ2-14): espectro de amplitud (línea base recortada de [0,40) − mínimo del pulso, buscador del
//    DAQ, separación 4, piso 1.0 mV) en todo el registro, bins de 0.25 mV hasta 200 mV (último = desborde).
//  - Exposición por canal = eventos usados × muestras × 16 ns.
//  - Con jsonl (controlador de la prueba de grupos): solo eventos dentro de ventanas; en las ventanas "off"
//    del grupo g se excluyen los canales cuyo grupo HV medido (mapeo) es g + 1.
// Salida: CSV ch,bin_lo_mV,counts (bins con cuentas) y una línea "#exposure ch s".
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/validation/selftrig_noise.C+("<run>_seg*.root", "<toggle.jsonl o vacío>", "<mapeo.json>", "<salida.csv>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/pulse_finder.h"

namespace {
struct Win {
  double t0 = 0, t1 = 0;
  int offGroup = -1;  // grupo HV apagado (1-8) o −1
};
}  // namespace

void selftrig_noise(const char* pattern, const char* jsonl, const char* mappingJson, const char* out, int nch = 32) {
  std::vector<int> hvGroup(64, 0);
  {
    std::ifstream in(mappingJson);
    auto m = nlohmann::json::parse(in);
    for (const auto& c : m["channels"])
      if (c["hv_group"].is_number()) hvGroup[c["dig_ch"].get<int>()] = c["hv_group"].get<int>();
  }
  std::vector<Win> wins;
  const bool useWins = jsonl && *jsonl;
  if (useWins) {
    std::ifstream in(jsonl);
    std::string line;
    while (std::getline(in, line)) {
      auto j = nlohmann::json::parse(line);
      if (j.value("ev", "") != "window") continue;
      Win w;
      w.t0 = j["start"];
      w.t1 = j["end"];
      if (j.value("tag", "") == "off") w.offGroup = j["group"].get<int>() + 1;
      wins.push_back(w);
    }
  }
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  Long64_t tread = 0;
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchStatus("t_read_utc_ns", 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("t_read_utc_ns", &tread);
  const double lsb = daq::pulse::lsbMv();
  const double binMv = 0.25;
  const int nb = static_cast<int>(200.0 / binMv);
  std::vector<std::vector<long>> spec(static_cast<std::size_t>(nch), std::vector<long>(static_cast<std::size_t>(nb) + 1, 0));
  std::vector<long> used(static_cast<std::size_t>(nch), 0);
  const std::vector<daq::pulse::Region> zones{{0, n}};
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    const Long64_t le = t.LoadTree(e);
    t.GetTree()->GetBranch("t_read_utc_ns")->GetEntry(le);
    int off = -1;
    if (useWins) {
      const double ts = static_cast<double>(tread) * 1e-9;
      const Win* w = nullptr;
      for (const auto& x : wins)
        if (ts >= x.t0 && ts < x.t1) w = &x;
      if (!w) continue;
      off = w->offGroup;
    }
    t.GetEntry(e);
    for (int c = 0; c < nch; ++c) {
      if (off > 0 && hvGroup[static_cast<std::size_t>(c)] == off) continue;
      ++used[static_cast<std::size_t>(c)];
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      std::vector<daq::pulse::Pulse> p;
      daq::pulse::findPulses(x, n, bl, lsb, 1.0, 4, zones, p);
      for (const auto& q : p) {
        const int b = std::min(nb, static_cast<int>(q.amplitudeMv / binMv));
        ++spec[static_cast<std::size_t>(c)][static_cast<std::size_t>(b)];
      }
    }
  }
  FILE* f = fopen(out, "w");
  fprintf(f, "# selftrig_noise %s; muestras por evento %d\n", pattern, n);
  for (int c = 0; c < nch; ++c) fprintf(f, "#exposure %d %.6f\n", c, used[static_cast<std::size_t>(c)] * n * 16e-9);
  fprintf(f, "ch,bin_lo_mV,counts\n");
  for (int c = 0; c < nch; ++c)
    for (int b = 0; b <= nb; ++b)
      if (spec[static_cast<std::size_t>(c)][static_cast<std::size_t>(b)]) fprintf(f, "%d,%.2f,%ld\n", c, b * binMv, spec[static_cast<std::size_t>(c)][static_cast<std::size_t>(b)]);
  fclose(f);
  printf("eventos %lld, muestras %d, exposición ch0 %.4f s\n", ne, n, used[0] * n * 16e-9);
}
