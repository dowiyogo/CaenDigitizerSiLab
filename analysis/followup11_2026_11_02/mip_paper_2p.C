// mip_paper_2p.C — 2026-11-02: figura del artículo mip_spectra_paper_2p, una columna, dos paneles apilados (estilo proc_style.C).
//  Márgenes izquierdos iguales (0.17) para alinear los marcos; en (b) techo del eje y a 1.6 × el máximo para la leyenda.
//  (a) espectro de amplitud en p.e. de un canal de XB con G cerca de la mediana (fig_amp_xb.json de amp_spectrum_xb.py), dibujado como
//      amp_spectra_hod2018_pe (proc_plot_stage4.C): ventana signal en línea, fuera de tiempo en gris, log y, líneas a trazos en 5 y 10 p.e.
//  (b) el espectro etiquetado de XB, 29-54 mm, del panel (b) de mip_tagged_spectra_paper: mismos datos (spectra_meas.csv), mismo ajuste
//      (Landau ⊗ Gauss de 4 a 105 p.e. medidos, semilla y límites de mip_paper.C) y mismo dibujo.
// Uso (desde la raíz): root -l -b -q daq/tools/proceeding/proc_style.C \
//   'analysis/followup11_2026_11_02/mip_paper_2p.C("<fig_amp_xb.json>", "<spectra_meas.csv>", "<base>")'
#include "../followup10_2026_11_01/mip_fit.C"

#include <TLine.h>
#include <TPad.h>

#include <nlohmann/json.hpp>

void mip_paper_2p(const char* json, const char* spectra, const char* base) {
  proc_style();
  auto* c = proc_canvas("c_2p", 1, 1.25);
  c->Divide(1, 2, 0.0, 0.0);

  // (a) como proc_plot_stage4.C con un canal
  {
    std::ifstream in(json);
    auto j = nlohmann::json::parse(in);
    const auto edges = j["edges"].get<std::vector<double>>();
    const int nb = static_cast<int>(edges.size()) - 1;
    const auto& ch = j["channels"][0];
    c->cd(1);
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.04);
    gPad->SetLogy();
    auto* hi = new TH1F("hi_a", ";Amplitude (p.e.);Triggers / 0.5 p.e.", nb, edges.front(), edges.back());
    auto* ho = new TH1F("ho_a", "", nb, edges.front(), edges.back());
    for (int i = 0; i < nb; ++i) {
      hi->SetBinContent(i + 1, ch["in_time"][i].get<double>());
      ho->SetBinContent(i + 1, ch["out_of_time"][i].get<double>());
    }
    hi->GetYaxis()->SetTitleOffset(1.75);
    hi->SetMinimum(0.5);
    hi->SetMaximum(3.0 * std::max(hi->GetMaximum(), ho->GetMaximum()));
    hi->SetLineColor(proc::color(5));
    hi->Draw("HIST");
    ho->SetLineColor(proc::grey());
    ho->SetFillColorAlpha(proc::grey(), 0.35);
    ho->Draw("HIST SAME");
    hi->Draw("HIST SAME");
    for (double x : {5.0, 10.0}) {
      auto* l = new TLine(x, 0.5, x, hi->GetMaximum());
      l->SetLineStyle(2);
      l->SetLineColor(proc::color(6));
      l->Draw();
    }
    auto* lg = new TLegend(0.52, 0.58, 0.96, 0.78);
    lg->AddEntry(hi, "signal window", "l");
    lg->AddEntry(ho, "out-of-time window", "f");
    lg->Draw();
    auto* t = new TLatex(0.52, 0.93, Form("(a) ch%d, layer XB", ch["ch"].get<int>()));
    t->SetNDC(); t->SetTextAlign(13); t->Draw();
    auto* t2 = new TLatex(0.52, 0.85, Form("#it{G} = %.2f mV", ch["G_mV"].get<double>()));
    t2->SetNDC(); t2->SetTextAlign(13); t2->Draw();
    printf("(a) ch%d G %.5f mV\n", ch["ch"].get<int>(), ch["G_mV"].get<double>());
  }

  // (b) como el panel (b) de mip_paper.C
  {
    auto* h = new TH1D("hp_XB_far", "", 300, 0, 300);
    std::ifstream in(spectra);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#' || line.rfind("layer", 0) == 0) continue;
      std::stringstream ss(line);
      std::vector<std::string> v;
      std::string tok;
      while (std::getline(ss, tok, ',')) v.push_back(tok);
      if (v[0] != "XB" || v[1] != "far") continue;
      const int b = h->FindBin(std::stod(v[2]) + 0.5);
      h->SetBinContent(b, std::stod(v[6]));
      h->SetBinError(b, std::stod(v[7]));
    }
    const double lo = 4.0, hi = 105.0;
    const double integ = h->Integral(h->FindBin(lo + 0.5), h->FindBin(hi - 0.5));
    int bmax = h->FindBin(lo + 0.5);
    double best = -1e9;
    for (int b = h->FindBin(lo + 0.5) + 1; b < h->FindBin(hi - 0.5); b++) {
      const double m3 = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
      if (m3 > best) best = m3, bmax = b;
    }
    const double mp0 = h->GetBinCenter(bmax);
    auto* f = new TF1("fp_XB_far", langaufun, lo, hi, 4);
    f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
    f->SetParLimits(0, 0.01, 50); f->SetParLimits(1, 0.1, 105); f->SetParLimits(2, 0, 10 * std::max(integ, 1.0)); f->SetParLimits(3, 0.01, 50);
    f->SetNpx(1500);
    const int st = h->Fit(f, "RBQ0");
    const double md = f->GetMaximumX(lo, hi);
    printf("(b) XB far: status %d n %.0f moda %.2f (MPV %.2f) χ² %.1f/%d\n", st, integ, md, f->GetParameter(1), f->GetChisquare(), f->GetNDF());
    c->cd(2);
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.17); gPad->SetTopMargin(0.04);
    h->SetLineColor(proc::color(0)); h->SetMarkerColor(proc::color(0)); h->SetMarkerStyle(20); h->SetMarkerSize(0.9); h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("Amplitude (p.e.)");
    h->GetYaxis()->SetTitle("Muons / p.e.");
    h->GetYaxis()->SetTitleOffset(1.75);
    h->GetXaxis()->SetRangeUser(0, 60);
    double ymax = 0;
    for (int b = h->FindBin(0.5); b <= h->FindBin(59.5); b++) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    const double ytop = 1.6 * ymax;
    h->SetMinimum(0); h->SetMaximum(ytop);
    h->Draw("E1");
    f->SetLineColor(proc::color(5)); f->SetLineWidth(3); f->SetRange(4.0, 60);
    f->Draw("SAME");
    auto* l5 = new TLine(5, 0, 5, ytop), *l10 = new TLine(10, 0, 10, ytop);
    for (auto* l : {l5, l10}) { l->SetLineStyle(9); l->SetLineWidth(3); l->SetLineColor(proc::color(6)); l->Draw(); }
    h->Draw("E1 SAME");
    gPad->RedrawAxis();
    TLatex t; t.SetNDC(); t.SetTextAlign(33);
    t.DrawLatex(0.95, 0.93, "(b) XB, 29-54 mm");
    t.DrawLatex(0.95, 0.85, Form("mode %.1f p.e.", md));
    t.DrawLatex(0.95, 0.77, Form("#it{n}_{net} = %.0f", integ));
    auto* lg = new TLegend(0.50, 0.48, 0.96, 0.71);
    lg->AddEntry(h, "tagged, net", "pe");
    lg->AddEntry(f, "Landau#otimesGauss", "l");
    lg->AddEntry(l5, "#it{T} = 5, 10 p.e.", "l");
    lg->Draw();
  }
  proc_save(c, base);
}
