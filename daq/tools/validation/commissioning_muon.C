// commissioning_muon.C — DAQ2-6, paso 4.1 (solo lectura). Run con trigger externo y HV encendida.
//  1. Posición del muón: por canal, histograma de la muestra del mínimo del registro en los eventos
//     con amplitud (línea base − mínimo) > amp_mV (por omisión 20 mV en hod2018, 10 mV en hod2019).
//     Línea base = media recortada del pedestal [0,40). Se imprime la moda por estación y la
//     fracción en [sig_lo, sig_hi) (49-58 esperado, DAQ2-1).
//  2. Tasa fuera de tiempo por canal (pedestal+control+post = [0,49)+[102,144), 91 muestras),
//     hod2018 > 3 p.e. con las ganancias de v8 (solo referencia), hod2019 > 3 mV; separación 4.
//  3. Espectros de amplitud fuera de tiempo (mV, paso 1 LSB) y los histogramas de posición, en un
//     archivo ROOT.
// Uso: root -l -b -q 'daq/tools/validation/commissioning_muon.C("<dir>/<base>_seg*.root", "out.root")'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#include "daq/pulse_finder.h"

void commissioning_muon(const char* pattern, const char* outRoot, int sigLo = 49, int sigHi = 58, double amp18 = 20.0,
                        double amp19 = 10.0, const char* gains = "outputs/pe_calibration_v8.txt") {
  std::vector<double> G(64, NAN);
  {
    std::ifstream in(gains);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int ch;
      std::string pl;
      double p0, g;
      if (ss >> ch >> pl >> p0 >> g && ch >= 0 && ch < 64) G[ch] = g;
    }
  }
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchAddress("adc", adc.data());
  const double lsb = daq::pulse::lsbMv();
  TFile out(outRoot, "RECREATE");
  TH2D hPos("min_sample", "Muestra del mínimo en eventos con pulso grande;canal;muestra", 64, -0.5, 63.5, n, -0.5, n - 0.5);
  TH1D hPos18("min_sample_hod2018", "hod2018: muestra del mínimo (pulso > umbral);muestra;canales×eventos", n, -0.5, n - 0.5);
  TH1D hPos19("min_sample_hod2019", "hod2019: muestra del mínimo (pulso > umbral);muestra;canales×eventos", n, -0.5, n - 0.5);
  const int nb = static_cast<int>(100.0 / lsb);
  TH2D hAmp("oot_amp", "Amplitud de pulsos fuera de tiempo;canal;amplitud [mV]", 64, -0.5, 63.5, nb, 0, nb * lsb);
  std::vector<long> oot(64, 0);
  const std::vector<daq::pulse::Region> zones{{0, 49}, {102, std::min(144, n)}};
  const Long64_t ne = t.GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    for (int c = 0; c < 64; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      int kmin = 0;
      for (int k = 1; k < n; ++k)
        if (x[k] < x[kmin]) kmin = k;
      const double amp = (bl - x[kmin]) * lsb;
      if (amp > (c < 32 ? amp18 : amp19)) {
        hPos.Fill(c, kmin);
        (c < 32 ? hPos18 : hPos19).Fill(kmin);
      }
      std::vector<daq::pulse::Pulse> p;
      daq::pulse::findPulses(x, n, bl, lsb, 2 * lsb, 4, zones, p);
      const double thr = c < 32 ? 3.0 * G[c] : 3.0;
      for (const auto& q : p) {
        hAmp.Fill(c, q.amplitudeMv);
        if (std::isfinite(thr) && q.amplitudeMv > thr) ++oot[c];
      }
    }
  }
  hPos.Write();
  hPos18.Write();
  hPos19.Write();
  hAmp.Write();
  const double T = ne * 91 * 16e-9;
  auto report = [&](TH1D& h, const char* name) {
    const int mode = h.GetMaximumBin() - 1;
    const double inWin = h.Integral(sigLo + 1, sigHi), tot = h.Integral();
    printf("%s: %0.f canales×eventos con pulso grande; moda muestra %d; fracción en [%d,%d): %.3f; percentiles 10/50/90: ", name, tot,
           mode, sigLo, sigHi, tot > 0 ? inWin / tot : 0.0);
    double q[3], pr[3] = {0.1, 0.5, 0.9};
    h.GetQuantiles(3, q, pr);
    printf("%.1f %.1f %.1f\n", q[0], q[1], q[2]);
  };
  printf("eventos %lld, N = %d\n", ne, n);
  report(hPos18, "hod2018 (> 20 mV)");
  report(hPos19, "hod2019 (> 10 mV)");
  std::vector<double> r18, r19;
  FILE* csv = fopen((std::string(outRoot) + ".rates.csv").c_str(), "w");
  fprintf(csv, "ch,threshold_mV,oot_count,oot_rate_kHz,muon_mode_sample,muon_events\n");
  for (int c = 0; c < 64; ++c) {
    const double thr = c < 32 ? 3.0 * G[c] : 3.0;
    if (!std::isfinite(thr)) continue;
    TH1D* proj = hPos.ProjectionY(Form("p%d", c), c + 1, c + 1);
    const double rate = oot[c] / T / 1e3;
    if (c != 3) (c < 32 ? r18 : r19).push_back(rate);
    fprintf(csv, "%d,%.3f,%ld,%.3f,%d,%.0f\n", c, thr, oot[c], rate, proj->GetMaximumBin() - 1, proj->Integral());
    delete proj;
  }
  fclose(csv);
  auto med = [](std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return std::vector<double>{v[v.size() / 2], v.front(), v.back()};
  };
  const auto a = med(r18), b = med(r19);
  printf("tasa fuera de tiempo (91 muestras, %.4f s por canal): hod2018 > 3 p.e. (G de v8) mediana %.2f kHz (mín %.2f, máx %.2f, %zu "
         "canales); hod2019 > 3 mV mediana %.2f kHz (mín %.2f, máx %.2f, %zu canales)\n",
         T, a[0], a[1], a[2], r18.size(), b[0], b[1], b[2], r19.size());
  out.Close();
}
