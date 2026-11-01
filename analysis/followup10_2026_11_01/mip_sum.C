// mip_sum.C — 2026-11-01: suma de las cuatro capas de mayo (XA, XB, YA, YB) de los espectros etiquetados de mip_tagged.py, por intervalo
// (21-29, 29-54 mm) y variable (amplitud, carga), en p.e. primarios; mismo ajuste Landau ⊗ Gauss que mip_fit.C (langaufun, desde 2 hasta
// 60 p.e. primarios, es decir sobre el umbral de 3 p.e. medidos = 1.71 primarios). Neto y error por bin: suma de los netos, errores en
// cuadratura. Fracción de la distribución AJUSTADA por encima de 5 y 10 p.e. medidos (5/1.75 y 10/1.75 primarios): ∫_T^150 f / ∫_0^150 f, con
// la función extrapolada por debajo del inicio del ajuste.
// Salida: <out>/fits_may_sum.csv. Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/followup10_2026_11_01/mip_sum.C("<out>")'
#include "mip_fit.C"

void mip_sum(const char* out) {
  std::map<std::string, TH1D*> H;
  std::ifstream in(std::string(out) + "/spectra.csv");
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("campaign", 0) == 0) continue;
    std::stringstream ss(line);
    std::vector<std::string> v;
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(tok);
    if (v[0] != "may") continue;
    const std::string key = v[2] + "_" + v[3];
    if (!H.count(key)) H[key] = new TH1D(("hs_" + key).c_str(), "", 150, 0, 150);
    TH1D* h = H[key];
    const int b = h->FindBin(std::stod(v[4]) + 0.5);
    const double e = std::stod(v[9]);
    h->SetBinContent(b, h->GetBinContent(b) + std::stod(v[8]));
    h->SetBinError(b, std::hypot(h->GetBinError(b), e));
  }
  FILE* fo = fopen((std::string(out) + "/fits_may_sum.csv").c_str(), "w");
  fprintf(fo, "pos_bin,variable,n_net_total,n_net_in_fit_range,status,width,width_err,mpv,mpv_err,gsigma,gsigma_err,chi2,ndf,peak_pe,frac_above_5meas,frac_above_10meas\n");
  for (auto& [key, h] : H) {
    const double lo = 2.0, hi = 60.0;
    const double integ = h->Integral(h->FindBin(lo + 0.5), h->FindBin(hi - 0.5));
    int bmax = h->FindBin(lo + 0.5);
    double best = -1e9;
    for (int b = h->FindBin(lo + 0.5) + 1; b < h->FindBin(hi - 0.5); b++) {
      const double m3 = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
      if (m3 > best) best = m3, bmax = b;
    }
    const double mp0 = h->GetBinCenter(bmax);
    auto* f = new TF1(("fs_" + key).c_str(), langaufun, lo, hi, 4);
    f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
    f->SetParLimits(0, 0.01, 50); f->SetParLimits(1, 0.1, 60); f->SetParLimits(2, 0, 10 * std::max(integ, 1.0)); f->SetParLimits(3, 0.01, 50);
    f->SetNpx(1500);
    const int st = h->Fit(f, "RBQ0");
    f->SetRange(0, 150);
    const double tot = f->Integral(0, 150, 1e-6);
    const double f5 = f->Integral(5 / 1.75, 150, 1e-6) / tot, f10 = f->Integral(10 / 1.75, 150, 1e-6) / tot;
    const double peak = f->GetMaximumX(lo, hi);
    std::string kcsv = key;
    kcsv[kcsv.find('_')] = ',';
    fprintf(fo, "%s,%.1f,%.1f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%d,%.3f,%.4f,%.4f\n", kcsv.c_str(), h->Integral(), integ, st,
            f->GetParameter(0), f->GetParError(0), f->GetParameter(1), f->GetParError(1), f->GetParameter(3), f->GetParError(3), f->GetChisquare(), f->GetNDF(),
            peak, f5, f10);
    printf("%-12s n_net %.0f (en el rango %.0f) status %d MPV %.2f ± %.2f anchura %.2f ± %.2f σ %.2f ± %.2f χ²/ndf %.1f/%d pico %.2f | f(>5 med) %.4f f(>10 med) %.4f\n",
           key.c_str(), h->Integral(), integ, st, f->GetParameter(1), f->GetParError(1), f->GetParameter(0), f->GetParError(0), f->GetParameter(3),
           f->GetParError(3), f->GetChisquare(), f->GetNDF(), peak, f5, f10);
  }
  fclose(fo);
}
