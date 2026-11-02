// layer_fit.C — 2026-11-02: ajuste Landau ⊗ Gauss de los espectros por capa e intervalo de may_layer_targets.py (<out>/spectra_layers.csv,
// p.e. medidos). Ajuste idéntico al de mip_paper.C (langaufun de mip_fit.C, de 4 a 105 p.e., misma semilla y límites); moda = máximo de la
// función, f(>T) = fracción de la función sobre T en [0, 262.5].
// Errores: NSAMP vectores de parámetros de la gaussiana multivariada con la covarianza del ajuste (descomposición de Cholesky); se descartan los
// que salen de los límites del ajuste o dan parámetros no positivos; error = desviación estándar de moda, f(>5) y f(>10) sobre las muestras.
// Salida: <out>/layer_fit.csv. Uso (desde la raíz): root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/followup11_2026_11_02/layer_fit.C("<out>")'
#include "../followup10_2026_11_01/mip_fit.C"

#include <TDecompChol.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>

void layer_fit(const char* out, int NSAMP = 400) {
  std::map<std::string, TH1D*> H;
  std::ifstream in(std::string(out) + "/spectra_layers.csv");
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("layer", 0) == 0) continue;
    std::stringstream ss(line);
    std::vector<std::string> v;
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(tok);
    const std::string key = v[0] + "," + v[1];
    if (!H.count(key)) H[key] = new TH1D(Form("hl_%zu", H.size()), "", 300, 0, 300);
    const int b = H[key]->FindBin(std::stod(v[2]) + 0.5);
    H[key]->SetBinContent(b, std::stod(v[6]));
    H[key]->SetBinError(b, std::stod(v[7]));
  }
  TRandom3 rnd(1102);
  FILE* fo = fopen((std::string(out) + "/layer_fit.csv").c_str(), "w");
  fprintf(fo, "# layer_fit.C: Landau x Gauss de 4 a 105 p.e. medidos; errores con %d muestras de la covarianza\n", NSAMP);
  fprintf(fo, "layer,pos_bin,status,n_net_fit,mode_pe,mode_pe_err,mpv_param,chi2,ndf,f5,f5_err,f10,f10_err,n_samples_used\n");
  const double lo = 4.0, hi = 105.0, top = 262.5;
  const double plo[4] = {0.01, 0.1, 0, 0.01};
  for (auto& [key, h] : H) {
    const double integ = h->Integral(h->FindBin(lo + 0.5), h->FindBin(hi - 0.5));
    int bmax = h->FindBin(lo + 0.5);
    double best = -1e9;
    for (int b = h->FindBin(lo + 0.5) + 1; b < h->FindBin(hi - 0.5); b++) {
      const double m3 = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
      if (m3 > best) best = m3, bmax = b;
    }
    const double mp0 = h->GetBinCenter(bmax);
    auto* f = new TF1(Form("fl_%s", h->GetName()), langaufun, lo, hi, 4);
    f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
    const double phi[4] = {50, 105, 10 * std::max(integ, 1.0), 50};
    for (int i = 0; i < 4; i++) f->SetParLimits(i, plo[i], phi[i]);
    f->SetNpx(1500);
    TFitResultPtr r = h->Fit(f, "RBQ0S");
    const int st = r;
    auto* g = new TF1(Form("gl_%s", h->GetName()), langaufun, 0, top, 4);
    g->SetNpx(3000);
    auto stats = [&](const double* p, double& md, double& f5, double& f10) {
      g->SetParameters(p);
      md = g->GetMaximumX(lo, hi);
      const double tot = g->Integral(0, top, 1e-5);
      f5 = g->Integral(5, top, 1e-5) / tot;
      f10 = g->Integral(10, top, 1e-5) / tot;
    };
    double md, f5, f10;
    stats(f->GetParameters(), md, f5, f10);
    TMatrixDSym cov = r->GetCovarianceMatrix();
    TDecompChol ch(cov);
    double sm = 0, sm2 = 0, s5 = 0, s52 = 0, s10 = 0, s102 = 0;
    int used = 0;
    if (ch.Decompose()) {
      TMatrixD U = ch.GetU();  // cov = Uᵀ U
      for (int k = 0; k < NSAMP; k++) {
        double z[4], p[4];
        for (int i = 0; i < 4; i++) z[i] = rnd.Gaus();
        bool ok = true;
        for (int i = 0; i < 4; i++) {
          p[i] = f->GetParameter(i);
          for (int j = 0; j <= i; j++) p[i] += U(j, i) * z[j];
          if (p[i] <= plo[i] || p[i] >= phi[i] || p[i] <= 0) ok = false;
        }
        if (!ok) continue;
        double a, b5, b10;
        stats(p, a, b5, b10);
        sm += a, sm2 += a * a, s5 += b5, s52 += b5 * b5, s10 += b10, s102 += b10 * b10;
        used++;
      }
    }
    auto sd = [&](double s, double s2) { return used > 1 ? std::sqrt(std::max(0.0, s2 / used - (s / used) * (s / used)) * used / (used - 1)) : NAN; };
    fprintf(fo, "%s,%d,%.2f,%.4f,%.4f,%.4f,%.3f,%d,%.5f,%.5f,%.5f,%.5f,%d\n", key.c_str(), st, integ, md, sd(sm, sm2), f->GetParameter(1), f->GetChisquare(),
            f->GetNDF(), f5, sd(s5, s52), f10, sd(s10, s102), used);
    printf("%-10s status %d n %.0f moda %.3f ± %.3f χ² %.1f/%d f5 %.4f ± %.4f f10 %.4f ± %.4f (%d muestras)\n", key.c_str(), st, integ, md, sd(sm, sm2),
           f->GetChisquare(), f->GetNDF(), f5, sd(s5, s52), f10, sd(s10, s102), used);
    fflush(stdout);
  }
  fclose(fo);
}
