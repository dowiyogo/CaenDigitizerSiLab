// proc_plot_sensitive_length.C — PROC-3 etapa 2: figura sensitive_length_hod2019 (una columna). Modelo de longitud sensible (estimado):
// P(X), P(Y) y P(X ∧ Y) por trigger contra L (la misma L para las dos coordenadas), y los valores medidos del run B como líneas
// horizontales (rX, rY y N_XY/N a 35 mV y a 3 mV; a 3 mV N_XY incluye el azar). Datos: work/fig_sensitive_length.json
// (sensitive_length.py). Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_sensitive_length.C("<json>", "<base>")'
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>

#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

void proc_plot_sensitive_length(const char* json, const char* base) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  auto* c = proc_canvas("c_sl", 1, 0.85);
  gPad->SetLogy();
  auto* fr = new TH1F("fr", ";Sensitive length #it{L} from the MPPC end (mm);Probability per trigger", 10, 0, 76);
  fr->SetMinimum(1e-6);
  fr->SetMaximum(2.0);
  fr->Draw("AXIS");
  const char* keys[3] = {"PX", "PY", "PXY"};
  const char* names[3] = {"model #it{P}(X)", "model #it{P}(Y)", "model #it{P}(X#wedgeY)"};
  const int col[3] = {proc::color(5), proc::color(6), proc::color(0)};
  const int sty[3] = {1, 2, 1};
  auto* lg = new TLegend(0.42, 0.17, 0.96, 0.50);
  lg->SetFillStyle(0);
  for (int k = 0; k < 3; ++k) {
    auto* g = new TGraph();
    for (std::size_t i = 0; i < j["L"].size(); ++i)
      if (j[keys[k]][i].get<double>() > 0) g->SetPoint(g->GetN(), j["L"][i].get<double>(), j[keys[k]][i].get<double>());
    g->SetLineColor(col[k]);
    g->SetLineStyle(sty[k]);
    g->SetLineWidth(3);
    g->Draw("L SAME");
    lg->AddEntry(g, names[k], "l");
  }
  const std::pair<const char*, int> thr[2] = {{"35mV", 2}, {"3mV", 3}};
  for (const auto& [t, ls] : thr) {
    const auto& m = j["measured"][t];
    for (const auto& [key, cc] : {std::pair<const char*, int>{"rX", col[0]}, {"rY", col[1]}, {"pXY", col[2]}}) {
      auto* l = new TLine(0, m[key].get<double>(), 76, m[key].get<double>());
      l->SetLineColor(cc);
      l->SetLineStyle(ls);
      l->SetLineWidth(1);
      l->Draw();
    }
    auto* tx = new TLatex(4, m["rX"].get<double>() * 1.3, t[0] == '3' && t[1] == '5' ? "35 mV" : "3 mV");
    tx->SetTextAlign(11);
    tx->Draw();
  }
  lg->AddEntry((TObject*)nullptr, "lines: run B", "");
  lg->AddEntry((TObject*)nullptr, "35 mV dashed, 3 mV dotted", "");
  lg->Draw();
  gPad->RedrawAxis();
  proc_save(c, base);
}
