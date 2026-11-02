// mip_paper.C — 2026-11-02: figura del artículo de los espectros MIP etiquetados: mayo, 29-54 mm, cuatro capas (XA, XB, YA, YB), AMPLITUD
// en p.e. MEDIDOS (solo G del canal, sin dividir por F), desde <out>/spectra_meas.csv (mip_paper.py).
//  - Ajuste: el de mip_talk.C (langaufun de mip_fit.C, misma semilla y límites), con el rango pasado a p.e. medidos: de 4 a 105 (= 60 · 1.75).
//    El ajuste de p.e. primarios empezaba en 2 = 3.5 p.e. medidos; con bins de 1 p.e. medido el borde siguiente es 4. Empezar en 3 mete el bin
//    [3, 4) completo, que en YB lejos tiene un exceso en el umbral (109 ± 43) y lleva el ajuste al límite (MPV 0.10, σ 12).
//  - Tabla: moda de la función ajustada (no el parámetro MPV), f(>5) y f(>10) = fracción de la función ajustada sobre 5 y 10 p.e. medidos,
//    normalizada a su integral en [0, 262.5] (= [0, 150] primarios; extrapolada por debajo de 4), y f².
//  - Figura: estilo del artículo (proc_style.C, letra de 8 pt como eff_vs_threshold), dos columnas, 2 × 2 paneles; eje x de 0 a 60 p.e.
//    medidos, eje y desde 0; ajuste dibujado en su rango; líneas a trazos en 5 y 10 p.e.; en cada panel la moda y n_net (en el rango del ajuste).
// Salida: <out>/paper_table.csv, <out>/mip_tagged_spectra_paper.{png,pdf}.
// Uso (desde la raíz): root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/followup11_2026_11_02/mip_paper.C("<out>")'
#include "../followup10_2026_11_01/mip_fit.C"

#include <TLine.h>

void mip_paper(const char* out) {
  proc_style();
  std::map<std::string, TH1D*> H;
  std::ifstream in(std::string(out) + "/spectra_meas.csv");
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("layer", 0) == 0) continue;
    std::stringstream ss(line);
    std::vector<std::string> v;
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(tok);
    const std::string key = v[0] + "_" + v[1];
    if (!H.count(key)) H[key] = new TH1D(("hp_" + key).c_str(), "", 300, 0, 300);
    const int b = H[key]->FindBin(std::stod(v[2]) + 0.5);
    H[key]->SetBinContent(b, std::stod(v[6]));
    H[key]->SetBinError(b, std::stod(v[7]));
  }
  std::map<std::string, TF1*> Fn;
  std::map<std::string, double> mode, nnet;
  FILE* fo = fopen((std::string(out) + "/paper_table.csv").c_str(), "w");
  fprintf(fo, "layer,pos_bin,n_net_in_fit_range,mode_pe_meas,mpv_param,chi2,ndf,f_above_5,f_above_10,f2_5,f2_10\n");
  for (auto& [key, h] : H) {
    const double lo = 4.0, hi = 105.0, top = 262.5;
    const double integ = h->Integral(h->FindBin(lo + 0.5), h->FindBin(hi - 0.5));
    int bmax = h->FindBin(lo + 0.5);
    double best = -1e9;
    for (int b = h->FindBin(lo + 0.5) + 1; b < h->FindBin(hi - 0.5); b++) {
      const double m3 = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
      if (m3 > best) best = m3, bmax = b;
    }
    const double mp0 = h->GetBinCenter(bmax);
    auto* f = new TF1(("fp_" + key).c_str(), langaufun, lo, hi, 4);
    f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
    f->SetParLimits(0, 0.01, 50); f->SetParLimits(1, 0.1, 105); f->SetParLimits(2, 0, 10 * std::max(integ, 1.0)); f->SetParLimits(3, 0.01, 50);
    f->SetNpx(1500);
    const int st = h->Fit(f, "RBQ0");
    const double md = f->GetMaximumX(lo, hi);
    auto* g = new TF1(("fpg_" + key).c_str(), langaufun, 0, top, 4);
    g->SetParameters(f->GetParameters());
    g->SetNpx(3000);
    const double tot = g->Integral(0, top, 1e-6), f5 = g->Integral(5, top, 1e-6) / tot, f10 = g->Integral(10, top, 1e-6) / tot;
    std::string k = key;
    for (auto& ch : k) if (ch == '_') ch = ',';
    fprintf(fo, "%s,%.1f,%.3f,%.3f,%.2f,%d,%.4f,%.4f,%.4f,%.4f\n", k.c_str(), integ, md, f->GetParameter(1), f->GetChisquare(), f->GetNDF(), f5, f10, f5 * f5, f10 * f10);
    printf("%-9s status %d n %.0f moda %.2f (MPV %.2f, anchura %.2f, σ %.2f) χ² %.1f/%d f5 %.3f f10 %.3f f5² %.3f f10² %.3f\n", key.c_str(), st, integ, md,
           f->GetParameter(1), f->GetParameter(0), f->GetParameter(3), f->GetChisquare(), f->GetNDF(), f5, f10, f5 * f5, f10 * f10);
    Fn[key] = f; mode[key] = md; nnet[key] = integ;
  }
  fclose(fo);

  auto* c = proc_canvas("c_paper", 2, 0.62);
  c->Divide(2, 2, 0.0, 0.0);
  const char* L[4] = {"XA", "XB", "YA", "YB"};
  for (int i = 0; i < 4; i++) {
    c->cd(i + 1);
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.04);
    const std::string key = std::string(L[i]) + "_far";
    TH1D* h = H[key];
    h->SetLineColor(proc::color(0)); h->SetMarkerColor(proc::color(0)); h->SetMarkerStyle(20); h->SetMarkerSize(0.9); h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("Amplitude (p.e.)");
    h->GetYaxis()->SetTitle("Muons / p.e.");
    h->GetYaxis()->SetTitleOffset(1.75);
    h->GetXaxis()->SetRangeUser(0, 60);
    double ymax = 0;
    for (int b = h->FindBin(0.5); b <= h->FindBin(59.5); b++) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    const double ytop = 1.25 * ymax;
    h->SetMinimum(0); h->SetMaximum(ytop);
    h->Draw("E1");
    TF1* f = Fn[key];
    f->SetLineColor(proc::color(5)); f->SetLineWidth(3); f->SetRange(4.0, 60);
    f->Draw("SAME");
    auto* l5 = new TLine(5, 0, 5, ytop), *l10 = new TLine(10, 0, 10, ytop);
    for (auto* l : {l5, l10}) { l->SetLineStyle(9); l->SetLineWidth(3); l->SetLineColor(proc::color(6)); l->Draw(); }
    h->Draw("E1 SAME");
    gPad->RedrawAxis();
    TLatex t; t.SetNDC(); t.SetTextAlign(33);
    t.DrawLatex(0.95, 0.93, Form("(%c) %s, 29-54 mm", 'a' + i, L[i]));
    t.DrawLatex(0.95, 0.80, Form("mode %.1f p.e.", mode[key]));
    t.DrawLatex(0.95, 0.67, Form("#it{n}_{net} = %.0f", nnet[key]));
    if (i == 0) {
      auto* lg = new TLegend(0.55, 0.30, 0.96, 0.60);
      lg->AddEntry(h, "tagged, net", "pe");
      lg->AddEntry(f, "Landau#otimesGauss", "l");
      lg->AddEntry(l5, "#it{T} = 5, 10 p.e.", "l");
      lg->Draw();
    }
  }
  proc_save(c, std::string(out) + "/mip_tagged_spectra_paper");
}
