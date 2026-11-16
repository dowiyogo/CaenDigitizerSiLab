// fit_dark2.C — seguimiento 2026-11-15: ajuste de los picos de amplitud de los pulsos oscuros aislados de dark2_extract.C (misma definición de
// amplitud que aS: pedestal del tramo, pulse_finder.h). Por canal de hod2018: histograma de a de los pulsos de una región (0 = pre-trigger
// [16,37), 1 = [102,132), 2 = las dos) en bins de 0.2 mV en [lo·G, hi·G] (G = la del análisis, solo para la ventana y la semilla) y ajuste χ² de
// gaussianas en p0 + n·Gf (n = nlo … nhi), alturas libres, σ_n = √(s0² + n·s1²). Gf = separación de picos con la definición de amplitud de la
// señal; p0 = desplazamiento común.
// Salida: <out>.csv (ch, G_used_mV, n_pulses, Gf_mV, Gf_err_mV, p0_mV, p0_err_mV, s0_mV, s1_mV, chi2, ndf, status).
// Uso: root -l -b -q 'analysis/ampcal_2026_11_15/fit_dark2.C+("<prefijo>_dark2.bin", "<G: 32 valores con coma>", "<out>", región, nlo, nhi,
//        lo, hi)'
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

namespace {
#pragma pack(push, 1)
struct Rec {
  uint32_t ev;
  uint8_t ch, region;
  uint16_t k;
  float a, q;
};
#pragma pack(pop)
int gNlo = 1, gNhi = 3;
double peaks(double* x, double* p) {
  // p: 0 p0, 1 G, 2 s0, 3 s1, 4… alturas de n = gNlo … gNhi
  double y = 0;
  for (int n = gNlo; n <= gNhi; ++n) {
    const double s = std::sqrt(p[2] * p[2] + n * p[3] * p[3]);
    const double m = p[0] + n * p[1];
    y += p[4 + n - gNlo] * std::exp(-0.5 * (x[0] - m) * (x[0] - m) / (s * s));
  }
  return y;
}
}  // namespace

void fit_dark2(const char* bin, const char* gains, const char* out, int region, int nlo, int nhi, double flo, double fhi) {
  gNlo = nlo;
  gNhi = nhi;
  const int np = 4 + nhi - nlo + 1;
  std::vector<double> G;
  std::stringstream ss(gains);
  std::string tok;
  while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  std::vector<std::vector<float>> a(32);
  FILE* f = fopen(bin, "rb");
  Rec r;
  while (fread(&r, sizeof r, 1, f) == 1)
    if (region == 2 || r.region == region) a[r.ch].push_back(r.a);
  fclose(f);
  FILE* fo = fopen((std::string(out) + ".csv").c_str(), "w");
  fprintf(fo, "ch,G_used_mV,n_pulses,Gf_mV,Gf_err_mV,p0_mV,p0_err_mV,s0_mV,s1_mV,chi2,ndf,status\n");
  for (int c = 0; c < 32; ++c) {
    const double g = G[c], lo = flo * g, hi = fhi * g;
    const int nb = static_cast<int>(std::floor((hi - lo) / 0.2));
    TH1D h(Form("h%d", c), "", nb, lo, lo + 0.2 * nb);
    for (float v : a[c]) h.Fill(v);
    if (h.GetEntries() < 500) {
      fprintf(fo, "%d,%.5f,%zu,,,,,,,,,insufficient\n", c, g, a[c].size());
      continue;
    }
    TF1 fn(Form("f%d", c), peaks, lo, lo + 0.2 * nb, np);
    const double hmax = h.GetMaximum();
    fn.SetParameter(0, -0.25 * g);
    fn.SetParameter(1, g);
    fn.SetParameter(2, 0.4);
    fn.SetParameter(3, 0.3);
    for (int n = nlo; n <= nhi; ++n) fn.SetParameter(4 + n - nlo, h.GetBinContent(h.FindBin(n * g - 0.25 * g)) + 1);
    fn.SetParLimits(1, 0.6 * g, 1.4 * g);
    fn.SetParLimits(2, 0.02, 3.0);
    fn.SetParLimits(3, 0.0, 3.0);
    for (int i = 4; i < np; ++i) fn.SetParLimits(i, 0, 10 * hmax);
    TFitResultPtr fr = h.Fit(&fn, "QSRN");
    const int st = fr;
    fprintf(fo, "%d,%.5f,%zu,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%d,%d\n", c, g, a[c].size(), fn.GetParameter(1), fn.GetParError(1), fn.GetParameter(0),
            fn.GetParError(0), std::fabs(fn.GetParameter(2)), std::fabs(fn.GetParameter(3)), fn.GetChisquare(), fn.GetNDF(), st);
  }
  fclose(fo);
}
