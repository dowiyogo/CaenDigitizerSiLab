// proc_plot_stage4.C — PROC-1 etapa 4: figura amp_spectra_hod2018_pe (espectros de amplitud en p.e. de 4 canales de hod2018
// elegidos por regla, ventana signal en color y ventana fuera de tiempo en gris; líneas a 5 y 10 p.e.) desde
// work/fig_stage4.json (stage4.py). Estilo: proc_style.C.
// PROC-4: chans = "" (los cuatro canales, doble columna) o una lista de canales, p. ej. "13,2" (una columna).
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage4.C("<fig_stage4.json>", "<base>"[, "13,2"])'
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>

#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

void proc_plot_stage4(const char* json, const char* base, const char* chans = "") {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  const auto edges = j["edges"].get<std::vector<double>>();
  const int nb = static_cast<int>(edges.size()) - 1;
  const std::map<std::string, std::string> rule{{"G mediana", "median gain"}, {"G mínima", "lowest gain"},
                                                {"más ruidoso a 5 p.e.", "noisiest at 5 p.e."}, {"Y de ocupación mediana", "Y, median occupancy"}};
  std::vector<int> sel;
  {
    std::string s(chans), tok;
    std::stringstream ss(s);
    while (std::getline(ss, tok, ','))
      if (!tok.empty()) sel.push_back(std::stoi(tok));
  }
  const bool two = !sel.empty();
  auto* c = two ? proc_canvas("c_spec", 1, 1.25) : proc_canvas("c_spec", 2, 0.62);
  if (two) c->Divide(1, static_cast<int>(sel.size()), 0.0, 0.0);
  else c->Divide(2, 2, 0.0, 0.0);
  int p = 0;
  std::vector<nlohmann::json> order;
  if (two) {
    for (int want : sel)
      for (const auto& ch : j["channels"])
        if (ch["ch"].get<int>() == want) order.push_back(ch);
  } else {
    for (const auto& ch : j["channels"]) order.push_back(ch);
  }
  for (const auto& ch : order) {
    c->cd(++p);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.17);
    gPad->SetTopMargin(0.04);
    gPad->SetLogy();
    auto* hi = new TH1F(Form("hi%d", p), ";Amplitude (p.e.);Triggers / 0.5 p.e.", nb, edges.front(), edges.back());
    auto* ho = new TH1F(Form("ho%d", p), "", nb, edges.front(), edges.back());
    for (int i = 0; i < nb; ++i) {
      hi->SetBinContent(i + 1, ch["in_time"][i].get<double>());
      ho->SetBinContent(i + 1, ch["out_of_time"][i].get<double>());
    }
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
    if (p == 1) {
      auto* lg = new TLegend(0.52, 0.58, 0.96, 0.78);
      lg->AddEntry(hi, "signal window", "l");
      lg->AddEntry(ho, "out-of-time window", "f");
      lg->Draw();
    }
    auto* t = new TLatex(0.52, 0.93, Form("(%c) ch%d, %s", 'a' + p - 1, ch["ch"].get<int>(), rule.at(ch["rule"].get<std::string>()).c_str()));
    t->SetNDC();
    t->SetTextAlign(13);
    t->Draw();
    auto* t2 = new TLatex(0.52, 0.85, Form("#it{G} = %.2f mV", ch["G_mV"].get<double>()));
    t2->SetNDC();
    t2->SetTextAlign(13);
    t2->Draw();
  }
  proc_save(c, base);
}
