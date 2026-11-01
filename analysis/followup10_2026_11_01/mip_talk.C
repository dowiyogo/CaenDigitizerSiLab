// mip_talk.C — 2026-11-01: versión para diapositiva de los espectros MIP etiquetados (mismos datos y mismo ajuste que mip_fit.C).
//  - Tabla: mayo, por capa e intervalo, amplitud y carga: moda de la función ajustada (no el parámetro MPV), f(>5) y f(>10) = fracción de la
//    función ajustada sobre 5/1.75 y 10/1.75 p.e. primarios, normalizada a su integral en [0, 150] (extrapolada por debajo de 2 p.e.), y f².
//  - Figura: mayo, 29-54 mm, carga, 4 paneles (XA, XB, YA, YB), p.e. primarios; eje x de 2.5 a 40 (desde 2.5 y no 3 para que se vea la línea
//    de T = 5 p.e. medidos, en 2.86), eje y desde 0; ajuste dibujado en su rango; líneas a 2.86 y 5.71 con su leyenda; moda (en YB «peak below
//    3 p.e.») y n_net (en el rango del ajuste). Lienzo de 16 × 10 cm a 300 ppp; textos de 0.06 de la altura del panel, márgenes 0.16 izq. y 0.14 abajo.
// Salida: <out>/talk_table.csv, <out>/mip_tagged_spectra_talk.{png,pdf}.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/followup10_2026_11_01/mip_talk.C("<out>")' (desde analysis/followup10_2026_11_01)
#include "mip_fit.C"

#include <TLine.h>

void mip_talk(const char* out) {
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
    const std::string key = v[1] + "_" + v[2] + "_" + v[3];
    if (!H.count(key)) H[key] = new TH1D(("ht_" + key).c_str(), "", 150, 0, 150);
    const int b = H[key]->FindBin(std::stod(v[4]) + 0.5);
    H[key]->SetBinContent(b, std::stod(v[8]));
    H[key]->SetBinError(b, std::stod(v[9]));
  }
  std::map<std::string, TF1*> Fn;
  std::map<std::string, double> mode, nnet;
  FILE* fo = fopen((std::string(out) + "/talk_table.csv").c_str(), "w");
  fprintf(fo, "layer,pos_bin,variable,n_net_in_fit_range,mode_pe,mpv_param,chi2,ndf,f_above_5meas,f_above_10meas,f2_5,f2_10\n");
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
    auto* f = new TF1(("ft_" + key).c_str(), langaufun, lo, hi, 4);
    f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
    f->SetParLimits(0, 0.01, 50); f->SetParLimits(1, 0.1, 60); f->SetParLimits(2, 0, 10 * std::max(integ, 1.0)); f->SetParLimits(3, 0.01, 50);
    f->SetNpx(1500);
    h->Fit(f, "RBQ0");
    const double md = f->GetMaximumX(lo, hi);
    auto* g = new TF1(("fg_" + key).c_str(), langaufun, 0, 150, 4);
    g->SetParameters(f->GetParameters());
    g->SetNpx(1500);
    const double tot = g->Integral(0, 150, 1e-6), f5 = g->Integral(5 / 1.75, 150, 1e-6) / tot, f10 = g->Integral(10 / 1.75, 150, 1e-6) / tot;
    std::string k = key;
    for (auto& ch : k) if (ch == '_') ch = ',';
    fprintf(fo, "%s,%.1f,%.3f,%.3f,%.2f,%d,%.4f,%.4f,%.4f,%.4f\n", k.c_str(), integ, md, f->GetParameter(1), f->GetChisquare(), f->GetNDF(), f5, f10, f5 * f5, f10 * f10);
    printf("%-18s n %.0f moda %.2f (MPV %.2f) χ² %.1f/%d f5 %.3f f10 %.3f f5² %.3f f10² %.3f\n", key.c_str(), integ, md, f->GetParameter(1), f->GetChisquare(),
           f->GetNDF(), f5, f10, f5 * f5, f10 * f10);
    Fn[key] = f; mode[key] = md; nnet[key] = integ;
  }
  fclose(fo);

  // figura: 16 × 10 cm a 300 ppp
  const int W = static_cast<int>(16.0 * proc::kPxPerCm + 0.5), Hh = static_cast<int>(10.0 * proc::kPxPerCm + 0.5);
  const int fAx = 75, fTx = 66;  // 18 pt y 16 pt a la escala de impresión
  auto* c = new TCanvas("c_talk", "", W, Hh);
  c->SetWindowSize(W + (W - c->GetWw()), Hh + (Hh - c->GetWh()));
  c->Divide(2, 2, 0.001, 0.001);
  const double tsz = 0.06;  // títulos, números y textos: 0.06 de la altura del panel
  const char* L[4] = {"XA", "XB", "YA", "YB"};
  const int col = proc::color(6);
  for (int i = 0; i < 4; i++) {
    c->cd(i + 1);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04);
    gPad->SetTicks(1, 1);
    const std::string key = std::string(L[i]) + "_far_charge";
    TH1D* h = H[key];
    h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(20); h->SetMarkerSize(1.2); h->SetLineWidth(3);
    for (auto* ax : {h->GetXaxis(), h->GetYaxis()}) {
      ax->SetLabelFont(42); ax->SetTitleFont(42); ax->SetLabelSize(tsz); ax->SetTitleSize(tsz); ax->SetNdivisions(505);
    }
    h->GetXaxis()->SetTitle("Primary p.e. (charge)");
    h->GetYaxis()->SetTitle("Muons / p.e.");
    h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(1.25);
    h->GetXaxis()->SetRangeUser(2.5, 40);
    double ymax = 0;
    for (int b = h->FindBin(2.6); b <= h->FindBin(39.9); b++) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    const double ytop = 1.15 * ymax;
    h->SetMinimum(0); h->SetMaximum(ytop);
    h->Draw("E1");
    TF1* f = Fn[key];
    f->SetLineColor(proc::color(0)); f->SetLineWidth(4); f->SetRange(2.0, 40);
    f->Draw("SAME");
    auto* l5 = new TLine(5 / 1.75, 0, 5 / 1.75, ytop);
    auto* l10 = new TLine(10 / 1.75, 0, 10 / 1.75, ytop);
    l5->SetLineStyle(9); l5->SetLineWidth(4); l5->SetLineColor(proc::color(5));
    l10->SetLineStyle(2); l10->SetLineWidth(4); l10->SetLineColor(proc::color(3));
    l5->Draw(); l10->Draw();
    gPad->RedrawAxis();
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(tsz); t.SetTextAlign(33);
    t.DrawLatex(0.94, 0.93, (std::string(L[i]) + ", May, 29-54 mm").c_str());
    t.DrawLatex(0.94, 0.83, L[i] == std::string("YB") ? "peak below 3 p.e." : Form("mode %.1f p.e.", mode[key]));
    t.DrawLatex(0.94, 0.73, Form("n_{net} = %.0f", nnet[key]));
    auto* lg = new TLegend(0.60, 0.42, 0.94, 0.62);
    lg->SetTextFont(42); lg->SetTextSize(tsz); lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(l5, "T = 5 p.e.", "l");
    lg->AddEntry(l10, "T = 10 p.e.", "l");
    lg->Draw();
  }
  proc_save(c, std::string(out) + "/mip_tagged_spectra_talk");
}
