// proc_plot_stage7.C — PROC-1 etapa 7: figura hod2019_rates. (a) Mayo: mediana por canal de la tasa fuera de tiempo contra el
// índice del run (sesgo desconocido): hod2019 > 1 y > 3 mV, hod2018 > 0.5 y > 5 p.e. (G de v8); los valores 0 no se dibujan
// (escala logarítmica). (b) Octubre: tasa > 1 mV de hod2019 por hora, suma de los 32 canales; banda gris = intervalo sin datos
// en que cambió el nivel (causa no registrada). Datos: work/fig_stage7.json (stage7.py). Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage7.C("<fig_stage7.json>", "<base>")'
#include <TAxis.h>
#include <TBox.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>

#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

void proc_plot_stage7(const char* json, const char* base) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  auto* c = proc_canvas("c_h19", 2, 0.42);
  c->Divide(2, 1, 0.0, 0.0);
  // (a) mayo
  c->cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.05);
  gPad->SetLogy();
  auto* fa = new TH1F("fa", ";May run index (bias not recorded);Median out-of-time rate per channel (Hz)", 10, 0.5, 10.5);
  fa->SetMinimum(50);
  fa->SetMaximum(2e7);
  fa->Draw("AXIS");
  const auto& m = j["may"];
  const std::vector<std::pair<std::string, std::string>> ser{{"h19_1mV", "hod2019 > 1 mV"}, {"h19_3mV", "hod2019 > 3 mV"}, {"h18_0p5pe", "hod2018 > 0.5 p.e."},
                                                            {"h18_5pe", "hod2018 > 5 p.e."}};
  const int col[4] = {proc::color(6), proc::color(1), proc::color(5), proc::color(2)};
  const int mk[4] = {20, 24, 21, 25};
  auto* lg = new TLegend(0.50, 0.19, 0.96, 0.45);
  lg->SetFillStyle(1001);
  lg->SetFillColor(0);
  for (int s = 0; s < 4; ++s) {
    auto* g = new TGraphErrors();
    for (int i = 0; i < 10; ++i) {
      const double y = m[ser[s].first][i].get<double>();
      if (y > 0) g->SetPoint(g->GetN(), m["index"][i].get<double>(), y);
    }
    g->SetMarkerStyle(mk[s]);
    g->SetMarkerColor(col[s]);
    g->SetLineColor(col[s]);
    g->SetLineWidth(1);
    g->Draw("PL SAME");
    lg->AddEntry(g, ser[s].second.c_str(), "pl");
  }
  lg->Draw();
  auto* ta = new TLatex(0.19, 0.88, "(a) May v1-v10");
  ta->SetNDC();
  ta->Draw();
  gPad->RedrawAxis();
  // (b) octubre por hora
  c->cd(2);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.05);
  gPad->SetLogy();
  const double t0 = j["oct"]["t_s"][0].get<double>();
  const int n = static_cast<int>(j["oct"]["t_s"].size());
  const double tmax = (j["oct"]["t_s"][n - 1].get<double>() - t0) / 3600.0 + 2;
  auto* fb = new TH1F("fb", ";Hours since 2026-10-19 21:00 UTC;hod2019 rate > 1 mV, sum of 32 channels (Hz)", 10, -2, tmax);
  fb->SetMinimum(2e4);
  fb->SetMaximum(5e7);
  fb->Draw("AXIS");
  auto* band = new TBox((j["oct"]["change_lo_s"].get<double>() - t0) / 3600.0, 2e4, (j["oct"]["change_hi_s"].get<double>() - t0) / 3600.0, 5e7);
  band->SetFillColorAlpha(proc::grey(), 0.35);
  band->SetLineColor(0);
  band->Draw();
  auto* g = new TGraphErrors();
  for (int i = 0; i < n; ++i) {
    const double y = j["oct"]["rate"][i].get<double>();
    if (y <= 0) continue;
    const int k = g->GetN();
    g->SetPoint(k, (j["oct"]["t_s"][i].get<double>() - t0) / 3600.0 + 0.5, y);
    g->SetPointError(k, 0, j["oct"]["err"][i].get<double>());
  }
  g->SetMarkerStyle(20);
  g->SetMarkerColor(proc::color(6));
  g->SetLineColor(proc::color(6));
  g->Draw("P SAME");
  auto* tb = new TLatex(0.18, 0.88, "(b) October, per hour");
  tb->SetNDC();
  tb->Draw();
  auto* tc = new TLatex(0.5 * ((j["oct"]["change_lo_s"].get<double>() + j["oct"]["change_hi_s"].get<double>()) / 2 - t0) / 3600.0 * 2, 3e4, "no data");
  tc->SetTextAlign(21);
  tc->SetTextColor(proc::grey());
  tc->Draw();
  gPad->RedrawAxis();
  proc_save(c, base);
}
