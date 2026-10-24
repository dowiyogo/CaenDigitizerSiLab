// proc_plot_stage3.C — PROC-1 etapa 3: figuras map2d_hod2018_xy_sub (mapa X × Y de hod2018 en tiempo − mezclados con el término
// rX·rY devuelto, y mapa de mezclados) y occupancy_2018_2019 (impactos por tira en tiempo − fuera de tiempo; banda = MC de barras
// normalizado a la suma medida) desde work/fig_stage3.json (stage3.py). Tiras cuya posición no está validada: cajas rayadas grises
// (fuera del mapa, junto al eje; en la ocupación, detrás de la tira). La tira muerta (ch3) va rayada y con "dead". Estilo: proc_style.C.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C \
//        'daq/tools/proceeding/proc_plot_stage3.C("<fig_stage3.json>", "<dir figs>")'
#include <TAxis.h>
#include <TBox.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>

#include <algorithm>
#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

namespace {
TBox* hatch(double x1, double y1, double x2, double y2) {
  auto* b = new TBox(x1, y1, x2, y2);
  b->SetFillStyle(3354);
  b->SetFillColor(proc::grey());
  b->SetLineColor(proc::grey());
  b->Draw();
  return b;
}

void draw_map(const nlohmann::json& m, const char* key, const char* name, const char* ztitle) {
  const std::string tt = std::string(";hod2018 X strip coordinate;hod2018 Y strip coordinate;") + ztitle;
  auto* h = new TH2F(name, tt.c_str(), 16, -0.5, 15.5, 16, -0.5, 15.5);
  for (int i = 0; i < 16; ++i)
    for (int j = 0; j < 16; ++j) h->SetBinContent(i + 1, j + 1, m[key][i][j].get<double>());
  auto* fr = new TH2F(Form("%s_fr", name), tt.c_str(), 1, -1.6, 15.5, 1, -1.6, 15.5);
  double lo = 0, hi = 0;
  for (int i = 1; i <= 16; ++i)
    for (int k = 1; k <= 16; ++k) {
      lo = std::min(lo, h->GetBinContent(i, k));
      hi = std::max(hi, h->GetBinContent(i, k));
    }
  for (auto* x : {fr, h}) {
    x->SetMinimum(lo);
    x->SetMaximum(hi);
  }
  fr->Draw("AXIS");
  h->GetZaxis()->SetTitleOffset(2.1);
  h->Draw("COLZ SAME");
  const int dead = m["dead_X_coord"];
  hatch(dead - 0.5, -0.5, dead + 0.5, 15.5);
  for (int k = 0; k < 16; ++k) {
    if (!m["valid_X"][k].get<bool>()) hatch(k - 0.5, -1.6, k + 0.5, -0.6);
    if (!m["valid_Y"][k].get<bool>()) hatch(-1.6, k - 0.5, -0.6, k + 0.5);
  }
  gPad->RedrawAxis();
}
}  // namespace

void proc_plot_stage3(const char* json, const char* dir) {
  proc_style();
  std::ifstream in(json);
  auto j = nlohmann::json::parse(in);
  // ---------------------------------------------------------------- mapas (mayo y octubre, dos paneles cada uno)
  for (const auto& [key, out] : {std::pair<std::string, std::string>{"may8_good", "map2d_hod2018_xy_sub"},
                                 std::pair<std::string, std::string>{"oct8", "map2d_hod2018_xy_sub_oct"}}) {
    const auto& m = j["maps"][key];
    auto* c = proc_canvas(("c_" + out).c_str(), 2, 0.46);
    c->Divide(2, 1, 0.0, 0.0);
    const char* zt[2] = {"In time minus mixed (weighted triggers)", "Mixed events (weighted triggers)"};
    const char* kk[2] = {"sub", "mixed"};
    for (int p = 0; p < 2; ++p) {
      c->cd(p + 1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.25);
      gPad->SetBottomMargin(0.15);
      gPad->SetTopMargin(0.09);
      draw_map(m, kk[p], Form("h_%s_%d", key.c_str(), p), zt[p]);
      auto* t = new TLatex(0.14, 0.985, Form("(%c) %s", 'a' + p, m["label"].get<std::string>() == "mayo v8, tiempo bueno" ? "May v8, 5 p.e." : "Oct run 8, 5 p.e."));
      t->SetNDC();
      t->SetTextAlign(13);
      t->SetTextSize(proc::kFontPx);
      t->Draw();
    }
    proc_save(c, std::string(dir) + "/" + out);
  }
  // ---------------------------------------------------------------- ocupación (cuatro paneles)
  auto* c = proc_canvas("c_occ", 2, 0.62);
  c->Divide(2, 2, 0.0, 0.0);
  int p = 0;
  for (const auto& panel : j["occ"]) {
    c->cd(++p);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.17);
    gPad->SetTopMargin(0.04);
    double ymax = 0, ymin = 0;
    for (const auto& s : panel["series"])
      for (int i = 0; i < 16; ++i) {
        ymax = std::max(ymax, s["y"][i].get<double>() + s["e"][i].get<double>());
        ymin = std::min(ymin, s["y"][i].get<double>() - s["e"][i].get<double>());
      }
    auto* fr = new TH1F(Form("fr_occ%d", p), Form(";%s strip coordinate;%s", panel["panel"].get<std::string>().c_str(),
                                                   panel["ylabel"].get<std::string>().c_str()), 16, -0.5, 15.5);
    fr->SetMinimum(ymin < 0 ? 1.15 * ymin : 0);
    fr->SetMaximum(1.35 * ymax);
    fr->Draw("AXIS");
    for (int k = 0; k < 16; ++k)
      if (!panel["valid"][k].get<bool>() || panel["dead"][k].get<bool>()) hatch(k - 0.5, fr->GetMinimum(), k + 0.5, fr->GetMaximum());
    for (int k = 0; k < 16; ++k)
      if (panel["dead"][k].get<bool>()) {
        auto* t = new TLatex(k, 0.5 * ymax, "dead");
        t->SetTextAngle(90);
        t->SetTextAlign(22);
        t->Draw();
      }
    auto* lg = new TLegend(0.18, 0.80, 0.95, 0.95);
    lg->SetNColumns(3);
    lg->SetFillStyle(1001);
    lg->SetFillColor(0);
    if (!panel["mc"].is_null()) {
      auto* band = new TGraphErrors(16);
      for (int k = 0; k < 16; ++k) {
        band->SetPoint(k, k, panel["mc"][k].get<double>());
        band->SetPointError(k, 0.45, panel["mc_err"][k].get<double>());
      }
      band->SetFillColorAlpha(proc::color(2), 0.45);
      band->SetLineColor(proc::color(2));
      band->Draw("2 SAME");
      lg->AddEntry(band, "MC (norm.)", "f");
    }
    int k = 0;
    for (const auto& s : panel["series"]) {
      auto* g = new TGraphErrors(16);
      for (int i = 0; i < 16; ++i) {
        g->SetPoint(i, i + (k == 0 ? -0.12 : 0.12), s["y"][i].get<double>());
        g->SetPointError(i, 0, s["e"][i].get<double>());
      }
      const int col = proc::color(k == 0 ? 5 : 6);
      g->SetMarkerStyle(k == 0 ? 20 : 24);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->Draw("P SAME");
      lg->AddEntry(g, s["name"].get<std::string>().c_str(), "pe");
      ++k;
    }
    lg->Draw();
    auto* t = new TLatex(0.17, 0.78, Form("(%c)", 'a' + p - 1));
    t->SetNDC();
    t->SetTextAlign(13);
    t->Draw();
    gPad->RedrawAxis();
  }
  proc_save(c, std::string(dir) + "/occupancy_2018_2019");
}
