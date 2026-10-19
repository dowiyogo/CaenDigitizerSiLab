// commissioning_pe_preview.C — DAQ2-6, paso 4.2 (solo lectura). Run de SW trigger con registro
// completo y HV encendida: espectro de amplitud de pulsos oscuros por canal de hod2018 y un primer
// ajuste multigaussiano (vista previa de la calibración de octubre; NO se guarda como calibración).
//  - Pulsos en todo el registro salvo la zona del muón del registro antiguo [250, 280) (con SW
//    trigger no hay muón, se excluye por si acaso), línea base = media recortada de [20, 200).
//  - Ajuste: 4 picos, posiciones p0 + n·G (n = 1..4), anchos sqrt(s0² + n·s1²), alturas libres, en
//    [lo, hi] = [5, 26] mV (sobre el bulto de 1-4 mV de mínimos de ruido y apilamiento, que tapa el
//    pico de 1 p.e.), p0 en ±3 mV. Semilla de G: máximo de la autocorrelación del espectro entre 2 y
//    8 mV; se repite con semillas G/2 y 2G y se marca "ambiguo" si una alternativa converge a una G
//    distinta en más de 20 % con χ²/ndf a menos de 1.2×.
// Uso: root -l -b -q 'daq/tools/validation/commissioning_pe_preview.C("<dir>/<base>_seg*.root", "out.root")'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <cstdio>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
double multiGauss(double* x, double* p) {
  // p[0] p0, p[1] G, p[2] s0, p[3] s1, p[4..7] alturas de 1..4 p.e.
  double s = 0;
  for (int k = 0; k < 4; ++k) {
    const double mu = p[0] + (k + 1) * p[1];
    const double sg = std::sqrt(p[2] * p[2] + (k + 1) * p[3] * p[3]);
    s += p[4 + k] * std::exp(-0.5 * (x[0] - mu) * (x[0] - mu) / (sg * sg));
  }
  return s;
}
}  // namespace

void commissioning_pe_preview(const char* pattern, const char* outRoot, double seedG = 4.8, double lo = 5.0, double hi = 26.0) {
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchAddress("adc", adc.data());
  const double lsb = daq::pulse::lsbMv();
  const int nb = static_cast<int>(40.0 / lsb);
  std::vector<TH1D*> h(32);
  for (int c = 0; c < 32; ++c) h[c] = new TH1D(Form("amp_ch%d", c), Form("ch%d: amplitud de pulsos oscuros;mV;pulsos", c), nb, 0, nb * lsb);
  const std::vector<daq::pulse::Region> zones{{20, 250}, {280, n - 2}};
  const Long64_t ne = t.GetEntries();
  long pulses = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    for (int c = 0; c < 32; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{20, 200}, 2.0 / lsb);
      std::vector<daq::pulse::Pulse> p;
      daq::pulse::findPulses(x, n, bl, lsb, 1.5, 4, zones, p);
      for (const auto& q : p) h[c]->Fill(q.amplitudeMv);
      pulses += static_cast<long>(p.size());
    }
  }
  TFile out(outRoot, "RECREATE");
  printf("eventos %lld, N = %d, pulsos > 1.5 mV en hod2018: %ld\n", ne, n, pulses);
  std::vector<double> Gv8(32, NAN);
  {
    FILE* g = fopen("outputs/pe_calibration_v8.txt", "r");
    char line[512];
    while (g && fgets(line, sizeof line, g)) {
      int ch;
      char pl[16];
      double p0, gg;
      if (line[0] != '#' && sscanf(line, "%d %15s %lf %lf", &ch, pl, &p0, &gg) == 4 && ch >= 0 && ch < 32) Gv8[ch] = gg;
    }
    if (g) fclose(g);
  }
  printf("ch  G_mV   G_err  p0_mV  chi2/ndf  semillas alternativas       G_v8 (mayo, referencia)  estado\n");
  std::vector<double> gs;
  for (int c = 0; c < 32; ++c) {
    TH1D* hh = h[c];
    // Semilla de G: máximo de la autocorrelación del espectro entre 2 y 8 mV, usando solo la parte
    // sobre fit_lo (el bulto de 1-4 mV de mínimos de ruido y apilamiento tapa el pico de 1 p.e.).
    const int b0 = hh->FindBin(lo), b1 = hh->FindBin(hi);
    double bestAc = -1, seed = seedG;
    for (int lag = static_cast<int>(2.0 / lsb); lag <= static_cast<int>(8.0 / lsb); ++lag) {
      double ac = 0, m = 0;
      int nn = 0;
      for (int b = b0; b + lag <= b1; ++b) {
        m += hh->GetBinContent(b);
        ++nn;
      }
      m /= std::max(1, nn);
      for (int b = b0; b + lag <= b1; ++b) ac += (hh->GetBinContent(b) - m) * (hh->GetBinContent(b + lag) - m);
      if (ac > bestAc) {
        bestAc = ac;
        seed = lag * lsb;
      }
    }
    auto fit = [&](double g0, double& chi) {
      TF1 f(Form("f%d_%g", c, g0), multiGauss, lo, hi, 8);
      const double a = hh->GetBinContent(hh->FindBin(-0.5 + 2 * g0));
      f.SetParameters(-0.5, g0, 0.8, 0.4, 4 * a, 2 * a, a, a / 2);
      f.SetParLimits(1, 1.5, 12);
      f.SetParLimits(2, 0.1, 4);
      f.SetParLimits(3, 0.0, 2);
      for (int k = 4; k < 8; ++k) f.SetParLimits(k, 0, 1e9);
      f.SetParLimits(0, -3, 3);
      const int st = hh->Fit(&f, "QRN0");
      chi = f.GetNDF() > 0 ? f.GetChisquare() / f.GetNDF() : 1e9;
      return std::vector<double>{f.GetParameter(1), f.GetParError(1), f.GetParameter(0), static_cast<double>(st), seed};
    };
    double c0 = 0, ch = 0, c2 = 0;
    const auto r = fit(seed, c0);
    const auto rh = fit(seed / 2, ch);
    const auto r2 = fit(seed * 2, c2);
    // Ambiguo: una semilla alternativa converge a otra G (> 20 % distinta) con χ²/ndf parecido.
    auto alt = [&](const std::vector<double>& x, double chi) { return std::fabs(x[0] - r[0]) / r[0] > 0.2 && chi < 1.2 * c0; };
    const bool ambiguous = alt(rh, ch) || alt(r2, c2);
    const bool atLimit = std::fabs(r[2]) > 2.95;
    const char* state = hh->GetEntries() < 200 ? "pocos pulsos"
                        : (r[3] != 0 ? "fallo" : (ambiguous ? "ambiguo" : (atLimit ? "p0 en el límite" : "ok")));
    printf("%2d %6.3f %6.3f %6.3f %8.2f  semilla %5.2f  G(G/2)=%6.3f G(2G)=%6.3f  G_v8=%6.3f  %s (%.0f pulsos)\n", c, r[0], r[1], r[2],
           c0, seed, rh[0], r2[0], Gv8[c], state, hh->GetEntries());
    if (std::string(state) == "ok") gs.push_back(r[0]);
    hh->Write();
  }
  std::sort(gs.begin(), gs.end());
  if (!gs.empty())
    printf("G (ajustes ok, %zu canales): mediana %.3f mV, mín %.3f, máx %.3f — VISTA PREVIA, no es calibración\n", gs.size(), gs[gs.size() / 2],
           gs.front(), gs.back());
  out.Close();
}
