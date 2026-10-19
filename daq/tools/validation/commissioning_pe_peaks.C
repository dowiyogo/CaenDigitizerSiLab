// commissioning_pe_peaks.C — DAQ2-6, paso 4.2 (solo lectura). Vista previa de la ganancia de hod2018
// desde los espectros de pulsos oscuros que guarda commissioning_pe_preview.C (amp_chN):
//  1. espectro suavizado (media de 3 bins), máximos locales sobre lo_mV (5 mV: el bulto de 1-4 mV
//     de mínimos de ruido y apilamiento tapa el pico de 1 p.e.) hasta hi_mV;
//  2. una gaussiana por pico en ±1 mV alrededor del máximo;
//  3. G = pendiente de la recta posición(pico) contra su índice consecutivo (separación media).
// No determina a qué número de p.e. corresponde el primer pico visible (ni p0): VISTA PREVIA, no es
// calibración. Compara con G de v8 (mayo, otros voltajes) solo como referencia.
// Uso: root -l -b -q 'daq/tools/validation/commissioning_pe_peaks.C("outputs/daq2-6/s4b_pe_preview.root")'
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

void commissioning_pe_peaks(const char* in, double lo = 5.0, double hi = 26.0, const char* gains = "outputs/pe_calibration_v8.txt") {
  std::vector<double> Gv8(32, NAN);
  if (FILE* g = fopen(gains, "r")) {
    char line[512];
    while (fgets(line, sizeof line, g)) {
      int ch;
      char pl[16];
      double p0, gg;
      if (line[0] != '#' && sscanf(line, "%d %15s %lf %lf", &ch, pl, &p0, &gg) == 4 && ch >= 0 && ch < 32) Gv8[ch] = gg;
    }
    fclose(g);
  }
  TFile f(in, "READ");
  printf("ch  picos  posiciones_mV                         G_mV   G_err  residuo_mV  G_v8(mayo)\n");
  std::vector<double> gs;
  for (int c = 0; c < 32; ++c) {
    TH1D* h = dynamic_cast<TH1D*>(f.Get(Form("amp_ch%d", c)));
    if (!h || h->GetEntries() < 1000) {
      printf("%2d  sin datos suficientes\n", c);
      continue;
    }
    const int nb = h->GetNbinsX();
    std::vector<double> s(nb + 2, 0);
    for (int b = 2; b < nb; ++b) s[b] = (h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1)) / 3;
    std::vector<double> pos;
    for (int b = h->FindBin(lo); b <= h->FindBin(hi) && b + 2 < nb; ++b) {
      bool peak = s[b] > 0;
      for (int d = -2; d <= 2 && peak; ++d)
        if (d != 0 && s[b + d] > s[b]) peak = false;
      if (!peak || s[b] <= s[b - 3] * 1.02 || s[b] <= s[b + 3] * 1.02) continue;  // pico de verdad, no meseta
      const double x0 = h->GetBinCenter(b);
      TF1 g("g", "gaus", x0 - 1.0, x0 + 1.0);
      g.SetParameters(h->GetBinContent(b), x0, 0.6);
      h->Fit(&g, "QRN0");
      const double m = g.GetParameter(1);
      if (std::fabs(m - x0) < 1.0 && (pos.empty() || m - pos.back() > 1.5)) pos.push_back(m);
    }
    if (pos.size() < 3) {
      printf("%2d  %zu picos: insuficiente\n", c, pos.size());
      continue;
    }
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double n = static_cast<double>(pos.size());
    for (size_t i = 0; i < pos.size(); ++i) {
      sx += i;
      sy += pos[i];
      sxx += double(i) * i;
      sxy += i * pos[i];
    }
    const double G = (n * sxy - sx * sy) / (n * sxx - sx * sx), a = (sy - G * sx) / n;
    double r2 = 0;
    for (size_t i = 0; i < pos.size(); ++i) r2 += std::pow(pos[i] - a - G * i, 2);
    const double res = std::sqrt(r2 / n);
    const double gErr = pos.size() > 2 ? std::sqrt(r2 / (n - 2) / (sxx - sx * sx / n)) : NAN;
    printf("%2d  %5zu  ", c, pos.size());
    for (double p : pos) printf("%5.2f ", p);
    for (size_t k = pos.size(); k < 6; ++k) printf("      ");
    printf(" %6.3f %6.3f  %6.3f     %6.3f\n", G, gErr, res, Gv8[c]);
    if (res < 0.4 * G / 4) gs.push_back(G);
  }
  std::sort(gs.begin(), gs.end());
  if (!gs.empty())
    printf("G con picos alineados (residuo < G/10), %zu canales: mediana %.3f mV, mín %.3f, máx %.3f — VISTA PREVIA, no es calibración\n",
           gs.size(), gs[gs.size() / 2], gs.front(), gs.back());
}
