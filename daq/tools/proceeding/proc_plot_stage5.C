// proc_plot_stage5.C — PROC-1 etapa 5: figura double_hit_fraction (pares de tiras vecinas verdaderos por trigger con impacto en
// el plano, contra el umbral; bandas = MC de barras cos³θ con detección perfecta y corte de camino 0, 0.2 y 0.5 mm, ancho = las
// dos opciones de bloque arriba; línea = incidencia normal 15/33) desde work/fig_stage5.json (stage5.py). Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage5.C("<fig_stage5.json>", "<base>")'
#include <TBox.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>

#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

void proc_plot_stage5(const char* json, const char* base) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  auto* c = proc_canvas("c_f2", 1, 0.85);
  auto* fr = new TH1F("fr", ";Threshold #it{T} (p.e.);Neighbour pairs per hit plane", 10, 2.5, 12.5);
  fr->SetMinimum(0.0);
  fr->SetMaximum(0.75);
  fr->Draw("AXIS");
  int k = 0;
  const int mcCol[3] = {proc::color(2), proc::color(3), proc::color(1)};
  std::vector<TBox*> bands;
  for (const auto& m : j["mc"]) {
    auto* b = new TBox(2.5, m["min"].get<double>(), 12.5, m["max"].get<double>());
    b->SetFillColorAlpha(mcCol[k], 0.45);
    b->SetLineColor(mcCol[k]);
    b->Draw();
    bands.push_back(b);
    ++k;
  }
  auto* ln = new TLine(2.5, j["normal"].get<double>(), 12.5, j["normal"].get<double>());
  ln->SetLineStyle(2);
  ln->SetLineColor(proc::grey());
  ln->Draw();
  auto* lg = new TLegend(0.18, 0.74, 0.96, 0.96);
  lg->SetNColumns(2);
  const std::vector<int> mk{20, 21, 24, 25}, col{proc::color(5), proc::color(6), proc::color(5), proc::color(6)};
  const std::vector<double> dx{-0.15, -0.05, 0.05, 0.15};
  k = 0;
  for (const auto& s : j["series"]) {
    const int n = static_cast<int>(s["f2"].size());
    auto* g = new TGraphErrors(n), *gs = new TGraphErrors(n);
    for (int i = 0; i < n; ++i) {
      const double x = j["T"][i].get<double>() + dx[k];
      g->SetPoint(i, x, s["f2"][i].get<double>());
      g->SetPointError(i, 0, s["stat"][i].get<double>());
      gs->SetPoint(i, x, s["f2"][i].get<double>());
      gs->SetPointError(i, 0.05, s["syst"][i].get<double>());
    }
    gs->SetFillStyle(0);
    gs->SetLineColor(col[k]);
    gs->SetLineWidth(1);
    gs->Draw("2 SAME");
    g->SetMarkerStyle(mk[k]);
    g->SetMarkerColor(col[k]);
    g->SetLineColor(col[k]);
    g->Draw("P SAME");
    lg->AddEntry(g, s["name"].get<std::string>().c_str(), "pe");
    ++k;
  }
  k = 0;
  for (const auto& m : j["mc"]) lg->AddEntry(bands[k++], Form("MC, path > %.1f mm", m["cut_mm"].get<double>()), "f");
  lg->AddEntry(ln, "normal incidence 15/33", "l");
  lg->Draw();
  gPad->RedrawAxis();
  proc_save(c, base);
}
