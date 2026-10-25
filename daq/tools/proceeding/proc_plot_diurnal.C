// proc_plot_diurnal.C — PROC-2 etapa 5: figura trigger_rate_vs_time (una columna, dos paneles). (a) Tasa del trigger por hora
// (sin interrupciones > 20 s) relativa a la mediana de las tasas por hora de cada run, contra las horas desde el inicio del run; (b) la misma plegada
// en la hora local de Chile (UTC−4): media por hora del día ± error de la media, con el coseno ajustado (v8 y v10; en v1 el ajuste no
// describe los datos y no se dibuja). Runs v1, v8 y v10.
// Datos: work/fig_trigger_diurnal.json (trigger_diurnal.py). Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_diurnal.C("<json>", "<base>")'
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>

#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

void proc_plot_diurnal(const char* json, const char* base) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  auto* c = proc_canvas("c_diu", 1, 1.25);
  c->Divide(1, 2, 0.0, 0.0);
  const int col[3] = {proc::grey(), proc::color(5), proc::color(6)};
  const int mk[3] = {24, 20, 21};
  for (int p = 0; p < 2; ++p) {
    c->cd(p + 1);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.03);
    auto* fr = new TH1F(Form("fr%d", p), p == 0 ? ";Hours since start of run;Rate / run median" : ";Local hour (UTC#minus4);Rate / run median", 10, p == 0 ? -2 : 0,
                        p == 0 ? 122 : 24);
    fr->SetMinimum(p == 0 ? 0.0 : 0.6);
    fr->SetMaximum(p == 0 ? 1.4 : 1.6);
    fr->Draw("AXIS");
    auto* lg = new TLegend(0.20, 0.80, 0.95, 0.95);
    lg->SetNColumns(3);
    int k = 0;
    for (const auto& r : j["runs"]) {
      const double a = r["median_Hz"].get<double>();
      auto* g = new TGraphErrors();
      const char* xs = p == 0 ? "t_h" : "fold_h";
      const char* ys = p == 0 ? "rate" : "fold_mean";
      const char* es = p == 0 ? "err" : "fold_err";
      for (std::size_t i = 0; i < r[xs].size(); ++i) {
        g->SetPoint(static_cast<int>(i), r[xs][i].get<double>(), r[ys][i].get<double>() / a);
        g->SetPointError(static_cast<int>(i), 0, r[es][i].get<double>() / a);
      }
      g->SetMarkerStyle(mk[k]);
      g->SetMarkerColor(col[k]);
      g->SetLineColor(col[k]);
      g->SetMarkerSize(0.8);
      g->Draw("P SAME");
      if (p == 1 && r["draw_cos"].get<bool>()) {
        auto* f = new TF1(Form("cos%d", k), "[2] + [0]*cos(2*TMath::Pi()*x/24) + [1]*sin(2*TMath::Pi()*x/24)", 0, 24);
        f->SetParameters(r["cos"][1].get<double>() / a, r["cos"][2].get<double>() / a, r["cos"][0].get<double>() / a);
        f->SetLineColor(col[k]);
        f->SetLineWidth(2);
        f->Draw("SAME");
      }
      lg->AddEntry(g, r["name"].get<std::string>().c_str(), "pe");
      ++k;
    }
    lg->Draw();
    auto* t = new TLatex(0.20, 0.22, p == 0 ? "(a)" : "(b)");
    t->SetNDC();
    t->Draw();
    gPad->RedrawAxis();
  }
  proc_save(c, base);
}
