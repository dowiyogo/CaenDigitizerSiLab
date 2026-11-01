// mip_fit.C — 2026-11-01: ajuste Landau ⊗ Gauss (función langaus del tutorial de ROOT, convolución numérica) de los espectros MIP etiquetados de
// mip_tagged.py (spectra.csv: neto = mismo evento − mezclados, error por bin), en p.e. primarios, por capa, campaña, intervalo y variable
// (amplitud, carga). Rango del ajuste: desde el primer bin entero por encima del umbral de la sonda (3 p.e. medidos = 3/F primarios:
// mayo desde 2, octubre desde 3) hasta 60 p.e. primarios; χ² con los errores del neto (TH1::Fit, bins vacíos fuera).
// Parámetros: anchura de Landau, MPV (parámetro de posición de la Landau, con el corrimiento de langaus), área, σ de la Gauss. Se da también el
// máximo de la función convolucionada.
// Salida: <out>/fits.csv; figura <out>/mip_tagged_spectra.{png,pdf} (2 × 2 capas, mayo, intervalo cerca, amplitud y carga superpuestas, con los
// ajustes) y <out>/mip_tagged_spectra_all.{png,pdf} (los 16 espectros, para inspección).
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/followup10_2026_11_01/mip_fit.C("<out>")'
#include <TF1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// langaus del tutorial de ROOT (langaus.C): par[0] anchura de Landau, par[1] MPV, par[2] área, par[3] σ de la Gauss
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

struct Spec {
  TH1D* h;
  TF1* f = nullptr;
};

void mip_fit(const char* out) {
  std::map<std::string, Spec> S;
  std::ifstream in(std::string(out) + "/spectra.csv");
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("campaign", 0) == 0) continue;
    std::stringstream ss(line);
    std::vector<std::string> v;
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(tok);
    const std::string key = v[0] + "_" + v[1] + "_" + v[2] + "_" + v[3];
    if (!S.count(key)) S[key].h = new TH1D(("h_" + key).c_str(), "", 150, 0, 150);
    const int b = S[key].h->FindBin(std::stod(v[4]) + 0.5);
    S[key].h->SetBinContent(b, std::stod(v[8]));
    S[key].h->SetBinError(b, std::stod(v[9]));
  }
  FILE* fo = fopen((std::string(out) + "/fits.csv").c_str(), "w");
  fprintf(fo, "campaign,layer,pos_bin,variable,fit_lo_pe,fit_hi_pe,n_net_in_range,status,width,width_err,mpv,mpv_err,area,area_err,gsigma,gsigma_err,chi2,ndf,peak_pe\n");
  for (auto& [key, sp] : S) {
    std::stringstream ks(key);
    std::string camp, layer, pos, var;
    std::getline(ks, camp, '_'); std::getline(ks, layer, '_'); std::getline(ks, pos, '_'); std::getline(ks, var, '_');
    const double lo = camp == "may" ? 2.0 : 3.0, hi = 60.0;
    TH1D* h = sp.h;
    double integ = h->Integral(h->FindBin(lo + 0.5), h->FindBin(hi - 0.5));
    // semilla: máximo de la suma móvil de 3 bins en el rango
    int bmax = h->FindBin(lo + 0.5);
    double best = -1e9;
    for (int b = h->FindBin(lo + 0.5) + 1; b < h->FindBin(hi - 0.5); b++) {
      const double m3 = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
      if (m3 > best) best = m3, bmax = b;
    }
    const double mp0 = h->GetBinCenter(bmax);
    auto* f = new TF1(("f_" + key).c_str(), langaufun, lo, hi, 4);
    f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
    f->SetParLimits(0, 0.01, 50); f->SetParLimits(1, 0.1, 60); f->SetParLimits(2, 0, 10 * std::max(integ, 1.0)); f->SetParLimits(3, 0.01, 50);
    f->SetNpx(600);
    const int st = h->Fit(f, "RBQ0");
    sp.f = f;
    const double peak = f->GetMaximumX(lo, hi);
    fprintf(fo, "%s,%s,%s,%s,%.0f,%.0f,%.1f,%d,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.4f,%.4f,%.2f,%d,%.3f\n", camp.c_str(), layer.c_str(), pos.c_str(), var.c_str(),
            lo, hi, integ, st, f->GetParameter(0), f->GetParError(0), f->GetParameter(1), f->GetParError(1), f->GetParameter(2), f->GetParError(2),
            f->GetParameter(3), f->GetParError(3), f->GetChisquare(), f->GetNDF(), peak);
    printf("%-22s status %d MPV %6.2f ± %5.2f width %5.2f ± %5.2f σ %5.2f ± %5.2f χ²/ndf %6.1f/%d peak %6.2f n %.0f\n", key.c_str(), st, f->GetParameter(1),
           f->GetParError(1), f->GetParameter(0), f->GetParError(0), f->GetParameter(3), f->GetParError(3), f->GetChisquare(), f->GetNDF(), peak, integ);
  }
  fclose(fo);

  auto draw = [&](TVirtualPad* pad, const std::string& camp, const std::string& layer, const std::string& pos, bool legend, double xmax) {
    pad->cd();
    Spec& a = S[camp + "_" + layer + "_" + pos + "_amp"];
    Spec& q = S[camp + "_" + layer + "_" + pos + "_charge"];
    const int ca = proc::color(5), cq = proc::color(6);
    for (auto* sp : {&a, &q}) {
      const int col = sp == &a ? ca : cq;
      sp->h->SetLineColor(col); sp->h->SetMarkerColor(col); sp->h->SetMarkerStyle(sp == &a ? 20 : 24); sp->h->SetMarkerSize(0.9);
      sp->f->SetLineColor(col); sp->f->SetLineWidth(3); sp->f->SetLineStyle(sp == &a ? 1 : 2);
    }
    a.h->GetXaxis()->SetRangeUser(0, xmax);
    const double ymax = 1.25 * std::max(a.h->GetMaximum(), q.h->GetMaximum());
    a.h->SetMaximum(ymax); a.h->SetMinimum(std::min(0.0, std::min(a.h->GetMinimum(), q.h->GetMinimum())) * 1.1);
    a.h->GetXaxis()->SetTitle("Primary photoelectrons");
    a.h->GetYaxis()->SetTitle("Entries / p.e. (tagged, net)");
    a.h->Draw("E1");
    q.h->Draw("E1 SAME");
    a.f->Draw("SAME"); q.f->Draw("SAME");
    TLatex t; t.SetNDC(); t.SetTextFont(43); t.SetTextSize(proc::kFontPx);
    t.DrawLatex(0.62, 0.86, (layer + ", " + (camp == "may" ? "May" : "Oct") + ", " + (pos == "near" ? "21-29 mm" : "29-54 mm")).c_str());
    if (legend) {
      auto* lg = new TLegend(0.55, 0.56, 0.95, 0.82);
      lg->AddEntry(a.h, "amplitude", "pe"); lg->AddEntry(q.h, "charge", "pe"); lg->AddEntry(a.f, "Landau #otimes Gauss", "l");
      lg->Draw();
    }
  };
  auto* c = proc_canvas("c_mip", 2, 0.8);
  c->Divide(2, 2, 0.001, 0.001);
  const char* L[4] = {"XA", "XB", "YA", "YB"};
  for (int i = 0; i < 4; i++) {
    c->cd(i + 1);
    gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04);
    draw(gPad, "may", L[i], "near", i == 0, 50);
  }
  proc_save(c, std::string(out) + "/mip_tagged_spectra");
  auto* c2 = proc_canvas("c_all", 2, 1.6);
  c2->Divide(4, 4, 0.001, 0.001);
  int k = 1;
  for (const char* camp : {"may", "oct"})
    for (const char* pos : {"near", "far"})
      for (int i = 0; i < 4; i++) {
        c2->cd(k++);
        gPad->SetLeftMargin(0.2); gPad->SetBottomMargin(0.2); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.04);
        draw(gPad, camp, L[i], pos, false, 50);
      }
  proc_save(c2, std::string(out) + "/mip_tagged_spectra_all");
}
