// timing_extract.C — PROC-1 etapa 4 (solo lectura): distribución de tiempos de pulso en todo el registro de mayo (árbol
// data_dig, Ch0..Ch63 vector<float> en mV, 624 muestras de 16 ns). hitdef_extract.C solo guarda un tramo de 144 muestras, así
// que el tiempo de los pulsos fuera de ese tramo se saca aquí, con el mismo buscador de pulsos del DAQ (pulse_finder.h).
//  - Definición única extendida al registro: pedestal = media recortada (±2 mV de la mediana) de las muestras [0, 150)
//    (0-2384 ns); todos los pulsos del buscador (mínimo local, separación 4 muestras) en [150, 624); por cada pulso con
//    amplitud > umbral se llena la muestra de su mínimo. Umbrales: hod2018 5 p.e. (T·G) y 50 mV; hod2019 3 mV y 35 mV.
//  - Definición del manuscrito (gen_beamer_figs.C, figura timesync): pedestal = media de t < 2500 ns, mínimo global en
//    [3000, 7000] ns (ch0) o [3000, 5000] ns (hod2019), un tiempo por evento si ped − mín > umbral; ch0 a 5 y 50 mV, ch50,
//    ch53, ch34, ch46 a 35 mV; histograma de 100 bins en [0, 8000] ns como el original.
// Salida: <prefijo>_pulsetimes.csv (ch,threshold_mV_or_label,counts por muestra 0..623) y <prefijo>_manuscript.csv
// (ch,threshold_mV,counts de los 100 bins de 80 ns).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/proceeding/timing_extract.C+("<archivo>", "<G de hod2018: 32 valores con coma>", "<prefijo>")'
#include <TChain.h>

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

void timing_extract(const char* file, const char* gains, const char* prefix) {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  }
  constexpr int kN = 624, kPed = 150, kNb = 100;
  TChain t("data_dig");
  t.Add(file);
  t.SetBranchStatus("*", 0);
  std::vector<float>* cw[64] = {};
  for (int c = 0; c < 64; ++c) {
    t.SetBranchStatus(Form("Ch%d", c), 1);
    t.SetBranchAddress(Form("Ch%d", c), &cw[c]);
  }
  t.SetCacheSize(256LL << 20);
  // umbrales por canal: [0] = bajo (5 p.e. o 3 mV), [1] = alto (50 o 35 mV)
  std::vector<std::vector<unsigned>> h(128, std::vector<unsigned>(kN, 0));
  double thr[64][2];
  for (int c = 0; c < 64; ++c) {
    thr[c][0] = c < 32 ? 5.0 * (c < static_cast<int>(G.size()) ? G[static_cast<std::size_t>(c)] : 4.5) : 3.0;
    thr[c][1] = c < 32 ? 50.0 : 35.0;
  }
  // definición del manuscrito
  struct Ms {
    int ch;
    double thr, lo, hi;
  };
  const std::vector<Ms> ms{{0, 5, 3000, 7000}, {0, 50, 3000, 7000}, {50, 35, 3000, 5000}, {53, 35, 3000, 5000}, {34, 35, 3000, 5000}, {46, 35, 3000, 5000}};
  std::vector<std::vector<unsigned>> hm(ms.size(), std::vector<unsigned>(kNb, 0));
  const std::vector<daq::pulse::Region> all{{kPed, kN}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kN);
  const Long64_t ne = t.GetEntries();
  Long64_t nok = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    bool ok = true;
    for (int c = 0; c < 64; ++c)
      if (!cw[c] || static_cast<int>(cw[c]->size()) < kN) ok = false;
    if (!ok) continue;
    ++nok;
    for (int c = 0; c < 64; ++c) {
      for (int i = 0; i < kN; ++i) x[i] = (*cw[c])[static_cast<std::size_t>(i)];
      const double bl = daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, kPed}, 2.0);
      p.clear();
      daq::pulse::findPulses(x.data(), kN, bl, 1.0, 0.5, 4, all, p);
      for (const auto& u : p)
        for (int k = 0; k < 2; ++k)
          if (u.amplitudeMv > thr[c][k]) ++h[static_cast<std::size_t>(2 * c + k)][static_cast<std::size_t>(u.sample)];
    }
    for (std::size_t m = 0; m < ms.size(); ++m) {
      const auto& w = *cw[ms[m].ch];
      double s = 0;
      int nb = 0;
      for (int i = 0; i < kN; ++i)
        if (i * 16.0 < 2500.0) {
          s += w[static_cast<std::size_t>(i)];
          ++nb;
        }
      const double ped = s / nb;
      double vmin = ped, tmin = -1;
      for (int i = 0; i < kN; ++i) {
        const double ti = i * 16.0;
        if (ti >= ms[m].lo && ti <= ms[m].hi && w[static_cast<std::size_t>(i)] < vmin) {
          vmin = w[static_cast<std::size_t>(i)];
          tmin = ti;
        }
      }
      if (tmin > 0 && ped - vmin > ms[m].thr) ++hm[m][static_cast<std::size_t>(std::min(kNb - 1, static_cast<int>(tmin / 80.0)))];
    }
  }
  FILE* f = fopen((std::string(prefix) + "_pulsetimes.csv").c_str(), "w");
  fprintf(f, "# eventos %lld; pulsos por muestra (16 ns) del buscador del DAQ, pedestal [0,150)\nch,threshold_mV", nok);
  for (int i = 0; i < kN; ++i) fprintf(f, ",s%d", i);
  fprintf(f, "\n");
  for (int c = 0; c < 64; ++c)
    for (int k = 0; k < 2; ++k) {
      fprintf(f, "%d,%.4f", c, thr[c][k]);
      for (int i = 0; i < kN; ++i) fprintf(f, ",%u", h[static_cast<std::size_t>(2 * c + k)][static_cast<std::size_t>(i)]);
      fprintf(f, "\n");
    }
  fclose(f);
  f = fopen((std::string(prefix) + "_manuscript.csv").c_str(), "w");
  fprintf(f, "# eventos %lld; definición de gen_beamer_figs.C, 100 bins de 80 ns en [0, 8000) ns\nch,threshold_mV,lo_ns,hi_ns", nok);
  for (int i = 0; i < kNb; ++i) fprintf(f, ",b%d", i);
  fprintf(f, "\n");
  for (std::size_t m = 0; m < ms.size(); ++m) {
    fprintf(f, "%d,%.1f,%.0f,%.0f", ms[m].ch, ms[m].thr, ms[m].lo, ms[m].hi);
    for (int i = 0; i < kNb; ++i) fprintf(f, ",%u", hm[m][static_cast<std::size_t>(i)]);
    fprintf(f, "\n");
  }
  fclose(f);
  printf("eventos %lld, usados %lld\n", ne, nok);
}
