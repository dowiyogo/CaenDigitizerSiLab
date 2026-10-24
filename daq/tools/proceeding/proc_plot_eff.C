// proc_plot_eff.C — PROC-1 etapa 2: figura eff_vs_threshold (q_mix / f_2018 contra el umbral) desde work/fig_eff.json
// (stage2.py). Barras = error estadístico; cajas = syst del estimador (|q_mix − q_fuera|) / f_2018; banda gris = syst de
// aceptancia (error del MC ⊕ cota de posición) sobre la serie de mayo; estrella = convención publicada (p0 + 3.5 G) en su T
// equivalente (mediana por canal); eje derecho = q. Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_eff.C("<fig_eff.json>", "<base de salida>")'
#include <TAxis.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>

#include <fstream>
#include <vector>

#include <nlohmann/json.hpp>

void proc_plot_eff(const char* json, const char* base) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  const double f18 = j["f18"], rel = j["f18_rel_syst"];
  auto* c = proc_canvas("c_eff", 1, 0.78);
  c->SetRightMargin(0.15);
  auto* fr = new TH1F("fr", ";Threshold #it{T} (p.e.);#it{q} / #it{f}_{2018}", 10, 2.5, 12.5);
  fr->SetMinimum(0.0);
  fr->SetMaximum(1.3);
  fr->Draw("AXIS");
  const std::vector<int> mk{20, 21, 22};
  const std::vector<double> dx{-0.15, 0.0, 0.15};
  auto* lg = new TLegend(0.30, 0.72, 0.74, 0.95);
  int k = 0;
  for (const auto& s : j["series"]) {
    const int n = static_cast<int>(s["eps"].size());
    auto* g = new TGraphErrors(n), *gs = new TGraphErrors(n);
    for (int i = 0; i < n; ++i) {
      const double x = j["T"][i].get<double>() + dx[k];
      g->SetPoint(i, x, s["eps"][i]);
      g->SetPointError(i, 0, s["stat"][i]);
      gs->SetPoint(i, x, s["eps"][i]);
      gs->SetPointError(i, 0.08, s["syst"][i]);
    }
    const int col = proc::color(k == 0 ? 0 : (k == 1 ? 6 : 5));
    if (k == 0) {  // banda de aceptancia sobre mayo
      auto* band = new TGraphErrors(n);
      for (int i = 0; i < n; ++i) {
        band->SetPoint(i, j["T"][i].get<double>(), s["eps"][i]);
        band->SetPointError(i, 0.5, s["eps"][i].get<double>() * rel);
      }
      band->SetFillColorAlpha(proc::grey(), 0.35);
      band->SetLineColor(0);
      band->Draw("E3 SAME");
      lg->AddEntry(band, Form("MC acceptance #pm%.1f%%", 100 * rel), "f");
    }
    gs->SetFillStyle(0);
    gs->SetLineColor(col);
    gs->SetLineWidth(1);
    gs->Draw("2 SAME");
    g->SetMarkerStyle(mk[k]);
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->Draw("P SAME");
    lg->AddEntry(g, s["name"].get<std::string>().c_str(), "pe");
    ++k;
  }
  auto* p = new TGraphErrors(1);
  p->SetPoint(0, j["published"]["T"], j["published"]["eps"]);
  p->SetPointError(0, 0, j["published"]["err"]);
  p->SetMarkerStyle(29);
  p->SetMarkerSize(2.0);
  p->SetMarkerColor(proc::color(7));
  p->SetLineColor(proc::color(7));
  p->Draw("P SAME");
  lg->AddEntry(p, "May v8, workshop convention", "pe");
  lg->Draw();
  auto* ax = new TGaxis(12.5, 0.0, 12.5, 1.3, 0.0, 1.3 * f18, 505, "+L");
  ax->SetLabelFont(43);
  ax->SetLabelSize(proc::kFontPx);
  ax->SetTitleFont(43);
  ax->SetTitleSize(proc::kFontPx);
  ax->SetTitle("#it{q} (fraction of triggers)");
  ax->SetTitleOffset(1.9);
  ax->SetMaxDigits(2);
  ax->Draw();
  proc_save(c, base);
}
