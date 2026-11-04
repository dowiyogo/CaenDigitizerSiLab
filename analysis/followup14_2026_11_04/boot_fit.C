// boot_fit.C — 2026-11-04: repite el ajuste Landau ⊗ Gauss de layer_fit.C (seguimiento 11: mismos rango de 4 a 105 p.e., semilla, límites, opciones y langaufun)
// sobre las réplicas del bootstrap de eventos de layer_targets.py. Entrada: <out>/boot_rows.csv y <out>/boot_spectra.bin (double, (fila, réplica, {neto, error}, bin));
// réplica 0 = espectro nominal. Salida: <out>/boot_fit_<irow0>_<irow1>.csv (layer,pos_bin,rep,status,mode,f5,f10,chi2,ndf). Filas [irow0, irow1) para repartir procesos.
// Uso: root -l -b -q 'analysis/followup14_2026_11_04/boot_fit.C+("<out>", irow0, irow1)'
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {
// copia de langaufun de mip_fit.C (langaus del tutorial de ROOT)
double langaufun(double* x, double* par) {
  const double invsq2pi = 0.3989422804014, mpshift = -0.22278298;
  const double np = 100.0, sc = 5.0;
  const double mpc = par[1] - mpshift * par[0];
  const double xlow = x[0] - sc * par[3], xupp = x[0] + sc * par[3], step = (xupp - xlow) / np;
  double sum = 0.0;
  for (double i = 1.0; i <= np / 2; i++) {
    double xx = xlow + (i - .5) * step;
    double fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }
  return par[2] * step * sum * invsq2pi / par[3];
}
}  // namespace

void boot_fit(const char* out, int irow0, int irow1, int nrep = 201, int nbin = 300) {
  std::vector<std::string> keys;
  {
    std::ifstream in(std::string(out) + "/boot_rows.csv");
    std::string l;
    std::getline(in, l);
    while (std::getline(in, l)) keys.push_back(l.substr(0, l.rfind(',')));
  }
  FILE* fb = fopen((std::string(out) + "/boot_spectra.bin").c_str(), "rb");
  FILE* fo = fopen(Form("%s/boot_fit_%d_%d.csv", out, irow0, irow1), "w");
  const double lo = 4.0, hi = 105.0, top = 262.5;
  const double plo[4] = {0.01, 0.1, 0, 0.01};
  std::vector<double> buf(2 * nbin);
  TH1D h("h", "", nbin, 0, nbin);
  h.Sumw2();
  TF1 f("fl", langaufun, lo, hi, 4);
  f.SetNpx(1500);
  TF1 g("gl", langaufun, 0, top, 4);
  g.SetNpx(3000);
  for (int ir = irow0; ir < irow1 && ir < static_cast<int>(keys.size()); ++ir) {
    for (int rep = 0; rep < nrep; ++rep) {
      fseek(fb, (static_cast<long>(ir) * nrep + rep) * 2 * nbin * sizeof(double), SEEK_SET);
      if (fread(buf.data(), sizeof(double), buf.size(), fb) != buf.size()) { printf("lectura corta\n"); return; }
      for (int b = 0; b < nbin; ++b) {
        h.SetBinContent(b + 1, buf[b]);
        h.SetBinError(b + 1, buf[nbin + b]);
      }
      const double integ = h.Integral(h.FindBin(lo + 0.5), h.FindBin(hi - 0.5));
      int bmax = h.FindBin(lo + 0.5);
      double best = -1e9;
      for (int b = h.FindBin(lo + 0.5) + 1; b < h.FindBin(hi - 0.5); b++) {
        const double m3 = h.GetBinContent(b - 1) + h.GetBinContent(b) + h.GetBinContent(b + 1);
        if (m3 > best) best = m3, bmax = b;
      }
      const double mp0 = h.GetBinCenter(bmax);
      f.SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
      const double phi[4] = {50, 105, 10 * std::max(integ, 1.0), 50};
      for (int i = 0; i < 4; i++) f.SetParLimits(i, plo[i], phi[i]);
      TFitResultPtr r = h.Fit(&f, "RBQ0S");
      const int st = r;
      g.SetParameters(f.GetParameters());
      const double md = g.GetMaximumX(lo, hi), tot = g.Integral(0, top, 1e-5);
      fprintf(fo, "%s,%d,%d,%.4f,%.6f,%.6f,%.3f,%d\n", keys[ir].c_str(), rep, st, md, g.Integral(5, top, 1e-5) / tot, g.Integral(10, top, 1e-5) / tot, f.GetChisquare(), f.GetNDF());
    }
    fflush(fo);
    printf("fila %d %s lista\n", ir, keys[ir].c_str());
    fflush(stdout);
  }
  fclose(fo);
  fclose(fb);
}
