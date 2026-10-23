// charact_gain.C — DAQ2-15, parte B (solo lectura). G de hod2018 por ventana de tiempo con las reglas de la fase 5
// (daq/core/include/daq/gain_fit.h, igual que el monitor de daqrun): espectro de amplitud de los pulsos fuera de
// tiempo en pedestal + control + post, bins de 1 LSB (0.4884 mV), semilla = G de la calibración. Sirve para ventanas
// largas con pocos eventos por hora (p. ej. tarde de los runs 4 y 7), donde la hora del monitor queda "insufficient".
// Ventanas: CSV label,t0_utc_s,t1_utc_s. Salida: CSV label,ch,events,G_mV,status,peaks.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/characterization/charact_gain.C+("<run>_seg*.root", "<calibración.json>", "<ventanas.csv>", "<salida.csv>")'
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/gain_fit.h"
#include "daq/pulse_finder.h"

void charact_gain(const char* pattern, const char* calibJson, const char* windowsCsv, const char* out) {
  std::vector<double> G(32, 4.0);
  {
    std::ifstream in(calibJson);
    auto j = nlohmann::json::parse(in);
    for (const auto& e : j["channels"])
      if (e["ch"].get<int>() < 32 && e.contains("G_mV") && e["G_mV"].is_number()) G[e["ch"].get<int>()] = e["G_mV"].get<double>();
  }
  struct Win {
    std::string label;
    double t0, t1;
    long events = 0;
    std::vector<std::vector<unsigned>> spec;
  };
  std::vector<Win> wins;
  {
    std::ifstream in(windowsCsv);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#' || line.rfind("label", 0) == 0) continue;
      std::stringstream ss(line);
      Win w;
      std::string a, b;
      std::getline(ss, w.label, ',');
      std::getline(ss, a, ',');
      std::getline(ss, b, ',');
      w.t0 = std::stod(a);
      w.t1 = std::stod(b);
      wins.push_back(w);
    }
  }
  const double lsb = daq::pulse::lsbMv();
  const int nb = static_cast<int>(100.0 / lsb);
  for (auto& w : wins) w.spec.assign(32, std::vector<unsigned>(static_cast<std::size_t>(nb) + 1, 0));
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  Long64_t tutc = 0;
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchStatus("t_utc_ns", 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("t_utc_ns", &tutc);
  const std::vector<daq::pulse::Region> oot{{0, 40}, {40, 49}, {102, n}};
  std::vector<daq::pulse::Pulse> p;
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const double ts = static_cast<double>(tutc) * 1e-9;
    std::vector<Win*> in;
    for (auto& w : wins)
      if (ts >= w.t0 && ts < w.t1) in.push_back(&w);
    if (in.empty()) continue;
    for (auto* w : in) ++w->events;
    for (int c = 0; c < 32; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      p.clear();
      daq::pulse::findPulses(x, n, bl, lsb, 0.9, 4, oot, p);
      for (const auto& u : p) {
        const auto b = static_cast<std::size_t>(std::min(nb, static_cast<int>(u.amplitudeMv / lsb)));
        for (auto* w : in) ++w->spec[static_cast<std::size_t>(c)][b];
      }
    }
  }
  FILE* f = fopen(out, "w");
  fprintf(f, "label,ch,events,G_mV,status,peaks\n");
  for (const auto& w : wins)
    for (int c = 0; c < 32; ++c) {
      const auto r = daq::gain::fitGain(w.spec[static_cast<std::size_t>(c)], lsb, G[static_cast<std::size_t>(c)]);
      fprintf(f, "%s,%d,%ld,%.5f,%s,%d\n", w.label.c_str(), c, w.events, r.gainMv, r.status.c_str(), r.peaks);
    }
  fclose(f);
  printf("eventos %lld, ventanas %zu\n", ne, wins.size());
}
