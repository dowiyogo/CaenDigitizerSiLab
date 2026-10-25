// proc_plot_along_bar.C — PROC-2 etapa 2: figura eff_along_bar (eficiencia de cada capa de hod2018 contra la posición a lo largo
// de sus barras, dada por la tira de referencia de la coordenada ortogonal; dividida por la probabilidad geométrica del MC de barras),
// cuatro paneles (barras X y barras Y, 5 y 10 p.e.), capas A y B, recta (continua) y exponencial (discontinua) ajustadas. Datos:
// work/fig_along_bar.json (along_bar.py). Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_along_bar.C("<fig_along_bar.json>", "<base>")'
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>

#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

void proc_plot_along_bar(const char* json, const char* base) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  auto* c = proc_canvas("c_bar", 2, 0.62);
  c->Divide(2, 2, 0.0, 0.0);
  int p = 0;
  for (const auto& pan : j["panels"]) {
    c->cd(++p);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.17);
    gPad->SetTopMargin(0.04);
    auto* fr = new TH1F(Form("fr%d", p), Form(";%s;#varepsilon / #varepsilon_{geom. MC}", pan["xlabel"].get<std::string>().c_str()), 10, -17, 17);
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.7);
    fr->Draw("AXIS");
    auto* lg = new TLegend(0.18, 0.78, 0.70, 0.95);
    lg->SetNColumns(2);
    int k = 0;
    for (const auto& s : pan["series"]) {
      const int n = static_cast<int>(s["x"].size());
      auto* g = new TGraphErrors(n);
      for (int i = 0; i < n; ++i) {
        g->SetPoint(i, s["x"][i].get<double>() + (k == 0 ? -0.2 : 0.2), s["y"][i].get<double>());
        g->SetPointError(i, 0, s["e"][i].get<double>());
      }
      const int col = proc::color(k == 0 ? 5 : 6);
      g->SetMarkerStyle(k == 0 ? 20 : 21);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->Draw("P SAME");
      auto* fl = new TF1(Form("fl%d%d", p, k), "[0]+[1]*x", -16.5, 16.5);
      fl->SetParameters(s["line"][0].get<double>(), s["line"][1].get<double>());
      fl->SetLineColor(col);
      fl->SetLineWidth(2);
      fl->Draw("SAME");
      auto* fe = new TF1(Form("fe%d%d", p, k), "[0]*exp([1]*x)", -16.5, 16.5);
      fe->SetParameters(s["exp"][0].get<double>(), s["exp"][1].get<double>());
      fe->SetLineColor(col);
      fe->SetLineStyle(2);
      fe->SetLineWidth(2);
      fe->Draw("SAME");
      lg->AddEntry(g, s["name"].get<std::string>().c_str(), "pe");
      ++k;
    }
    lg->Draw();
    auto* t = new TLatex(0.18, 0.24, Form("(%c) %s", 'a' + p - 1, pan["title"].get<std::string>().c_str()));
    t->SetNDC();
    t->Draw();
    gPad->RedrawAxis();
  }
  proc_save(c, base);
}
