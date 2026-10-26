// plot_day.C — XPLUS-4b etapa 6: figuras de las etapas 2, 3, 4b y 5 (estilo PROC, proc_style.C) desde los CSV de export_plot_data.py.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/night_2026-10-25/plot_day.C("<dir plots>", "<dir salida>")'
#include <TAxis.h>
#include <TBox.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<std::vector<double>> readRows(const std::string& f) {
  std::vector<std::vector<double>> out;
  std::ifstream in(f);
  for (std::string l; std::getline(in, l);) {
    if (l.empty() || l[0] == '#') continue;
    std::istringstream s(l);
    std::vector<double> r;
    std::string tok;
    while (s >> tok) {
      try {
        r.push_back(std::stod(tok));
      } catch (...) {
        r.push_back(tok == "A" ? 0 : 1);
      }
    }
    out.push_back(r);
  }
  return out;
}
}  // namespace

void plot_day(const char* pdir, const char* odir) {
  proc_style();
  const std::string P(pdir), O(odir);
  // etapa 2: escalera por canal, una figura por condición, X e Y en dos paneles
  for (const char* tag : {"sw0V_all", "swm0p5V_all", "v8good_pre", "run8_post"}) {
    auto rows = readRows(P + "/stair_" + tag + ".csv");
    auto* c = proc_canvas(Form("c_st_%s", tag), 2, 0.45);
    c->Divide(2, 1);
    for (int pl = 0; pl < 2; ++pl) {
      c->cd(pl + 1);
      gPad->SetLogy();
      gPad->SetLeftMargin(0.16);
      auto* fr = new TH1F(Form("fr_%s_%d", tag, pl), Form(";umbral (p.e.);tasa oscura, hod2018 %s (s^{-1})", pl ? "Y" : "X"), 10, 0.3, 6.3);
      fr->SetMinimum(1e2);
      fr->SetMaximum(2e7);
      fr->Draw("AXIS");
      std::map<int, TGraph*> g;
      for (auto& r : rows) {
        const int ch = static_cast<int>(r[0]);
        if ((ch < 16) != (pl == 0)) continue;
        auto*& gg = g[ch];
        if (!gg) gg = new TGraph();
        gg->SetPoint(gg->GetN(), r[1], r[2]);
      }
      int k = 0;
      for (auto& [ch, gg] : g) {
        gg->SetLineColor(proc::color(k % 7 + 1));
        gg->SetMarkerColor(proc::color(k % 7 + 1));
        gg->SetMarkerStyle(20);
        gg->SetMarkerSize(0.5);
        gg->Draw("PL SAME");
        ++k;
      }
    }
    proc_save(c, O + "/stage2_staircase_" + tag);
  }
  // etapa 2: intervalos
  for (const char* tag : {"sw0V_all", "v8good_pre"}) {
    auto rows = readRows(P + "/dt_" + tag + ".csv");
    auto* go = new TGraph();
    auto* gp = new TGraph();
    for (auto& r : rows) {
      if (r[1] > 0) go->SetPoint(go->GetN(), r[0], r[1]);
      if (r[2] > 0) gp->SetPoint(gp->GetN(), r[0], r[2]);
    }
    auto* c = proc_canvas(Form("c_dt_%s", tag), 1, 0.8);
    c->SetLogy();
    go->SetMarkerStyle(20);
    go->SetMarkerSize(0.5);
    gp->SetLineColor(proc::color(6));
    gp->SetLineWidth(2);
    auto* mg = new TMultiGraph();
    mg->Add(go, "P");
    mg->Add(gp, "L");
    mg->SetTitle(";intervalo entre pulsos consecutivos > 0.5 p.e. (ns);pares (hod2018, suma de canales)");
    mg->Draw("A");
    auto* lg = new TLegend(0.45, 0.75, 0.93, 0.93);
    lg->AddEntry(go, "medido", "p");
    lg->AddEntry(gp, "exponencial ajustada en 0.64-2.4 #mus", "l");
    lg->Draw();
    proc_save(c, O + "/stage2_intervals_" + tag);
  }
  {
    auto rows = readRows(P + "/dcr_vs_dV.csv");
    auto* g1 = new TGraph();
    auto* g2 = new TGraph();
    for (auto& r : rows) {
      if (r[1] > 0) g1->SetPoint(g1->GetN(), r[0], r[1]);
      if (r[2] > 0) g2->SetPoint(g2->GetN(), r[0], r[2]);
    }
    auto* c = proc_canvas("c_dv", 1, 0.8);
    c->SetLogy();
    g1->SetMarkerStyle(20);
    g2->SetMarkerStyle(24);
    g2->SetMarkerColor(proc::color(6));
    auto* mg = new TMultiGraph();
    mg->Add(g1, "PL");
    mg->Add(g2, "P");
    mg->SetTitle(";#DeltaV respecto del nominal (V);DCR mediana de hod2018 (s^{-1})");
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(-1.7, 0.2);
    auto* lg = new TLegend(0.18, 0.75, 0.75, 0.93);
    lg->AddEntry(g1, "> 0.5 p.e. de G(0 V), mV fijos", "p");
    lg->AddEntry(g2, "> 0.5 p.e. del escal#acute{o}n", "p");
    lg->Draw();
    proc_save(c, O + "/stage2_dcr_vs_dV");
  }
  // etapa 4b: perfiles
  for (const char* tag : {"v8_T5", "run8_T5", "sw0V_T5"}) {
    auto rows = readRows(P + "/profile_" + std::string(tag) + ".csv");
    auto* gx = new TGraph();
    auto* gy = new TGraph();
    for (auto& r : rows) (r[0] == 0 ? gx : gy)->SetPoint((r[0] == 0 ? gx : gy)->GetN(), r[1] + 72, r[2]);
    auto* c = proc_canvas(Form("c_pr_%s", tag), 1, 0.65);
    c->SetLeftMargin(0.21);
    gx->SetLineColor(proc::color(5));
    gy->SetLineColor(proc::color(6));
    gx->SetLineWidth(2);
    gy->SetLineWidth(2);
    auto* mg = new TMultiGraph();
    mg->Add(gx, "L");
    mg->Add(gy, "L");
    mg->SetTitle(";centro de la ventana de 9 muestras (ns);fracci#acute{o}n con pulso > 5 p.e.");
    mg->Draw("A");
    mg->GetYaxis()->SetTitleOffset(2.0);
    c->Modified();
    const std::string t(tag);
    const double sig = t == "v8_T5" ? 258 * 16 : (t == "run8_T5" ? 49 * 16 : -1);
    if (sig > 0) {
      const double ymin = mg->GetYaxis()->GetXmin(), ymax = mg->GetYaxis()->GetXmax();
      auto* bs = new TBox(sig, ymin, sig + 144, ymax);
      bs->SetFillColorAlpha(proc::grey(), 0.4);
      bs->Draw();
      auto* bo = new TBox(sig + 800, ymin, sig + 944, ymax);
      bo->SetFillColorAlpha(proc::color(1), 0.3);
      bo->Draw();
    }
    auto* lg = new TLegend(0.74, 0.74, 0.95, 0.93);
    lg->AddEntry(gx, "hod2018 X", "l");
    lg->AddEntry(gy, "hod2018 Y", "l");
    lg->Draw();
    proc_save(c, O + "/stage4b_profile_" + t);
  }
  // etapa 3: q contra g_b
  for (const char* run : {"v8", "v10", "run8"}) {
    auto rows = readRows(P + "/stage3_" + std::string(run) + ".csv");
    auto fits = readRows(P + "/stage3_" + std::string(run) + "_fit.csv");
    auto* c = proc_canvas(Form("c_s3_%s", run), 1, 0.8);
    c->SetLeftMargin(0.20);
    auto* ga = new TGraphErrors();
    auto* gb = new TGraphErrors();
    for (auto& r : rows) {
      ga->SetPoint(ga->GetN(), r[0] - 0.001, r[2]);
      ga->SetPointError(ga->GetN() - 1, r[1], r[3]);
      gb->SetPoint(gb->GetN(), r[0] + 0.001, r[4]);
      gb->SetPointError(gb->GetN() - 1, r[1], r[5]);
    }
    ga->SetMarkerStyle(24);
    gb->SetMarkerStyle(20);
    gb->SetMarkerColor(proc::color(6));
    gb->SetLineColor(proc::color(6));
    auto* mg = new TMultiGraph();
    mg->Add(ga, "P");
    mg->Add(gb, "P");
    mg->SetTitle(";g_{b} = G_{bloque}/G_{run};q (5 p.e.)");
    mg->Draw("A");
    mg->GetYaxis()->SetTitleOffset(1.9);
    c->Modified();
    for (auto& f : fits) {
      auto* fl = new TF1(Form("f_%s_%d", run, static_cast<int>(f[0])), "[0]*(1+[1]*(x-1))", 0.9, 1.1);
      fl->SetParameters(f[1], f[2]);
      fl->SetLineColor(f[0] == 0 ? proc::color(0) : proc::color(6));
      fl->SetLineStyle(f[0] == 0 ? 2 : 1);
      fl->Draw("SAME");
    }
    auto* lg = new TLegend(0.18, 0.78, 0.80, 0.93);
    lg->AddEntry(ga, "modo A (umbral en mV fijo)", "p");
    lg->AddEntry(gb, "modo B (umbral 5 G_{bloque})", "p");
    lg->Draw();
    proc_save(c, O + "/stage3_q_vs_g_" + std::string(run));
  }
  // etapa 5: ε por run con la banda del combinado; G mediana contra fecha
  {
    auto rows = readRows(P + "/stage5_eps.csv");
    auto comb = readRows(P + "/stage5_comb.csv");
    auto* gv = new TGraphErrors();
    auto* gn = new TGraphErrors();
    for (auto& r : rows) {
      auto* g = r[3] > 0 ? gv : gn;
      g->SetPoint(g->GetN(), r[0], r[1]);
      g->SetPointError(g->GetN() - 1, 0, r[2]);
    }
    auto* c = proc_canvas("c_s5", 1, 0.75);
    auto* fr = new TH1F("fr_s5", ";run de mayo (v1 ... v10);#varepsilon_{XY} (5 p.e.)", 10, 0.5, 10.5);
    fr->SetMinimum(0);
    fr->SetMaximum(1.3);
    fr->Draw("AXIS");
    auto* bx = new TBox(0.5, comb[0][0] - comb[0][1], 10.5, comb[0][0] + comb[0][1]);
    bx->SetFillColorAlpha(proc::grey(), 0.4);
    bx->Draw();
    gv->SetMarkerStyle(20);
    gn->SetMarkerStyle(24);
    gn->Draw("P SAME");
    gv->Draw("P SAME");
    auto* lg = new TLegend(0.18, 0.75, 0.80, 0.93);
    lg->AddEntry(gv, "orden de tiras validado (v8, v10)", "p");
    lg->AddEntry(gn, "no validado (pocos eventos o condiciones distintas)", "p");
    lg->AddEntry(bx, "combinado v8 + v10", "f");
    lg->Draw();
    proc_save(c, O + "/stage5_eps_by_run");
    auto* gg = new TGraph();
    auto* gf = new TGraph();
    for (auto& r : rows) (r[6] > 0 ? gg : gf)->SetPoint((r[6] > 0 ? gg : gf)->GetN(), r[4], r[5]);
    auto* c2 = proc_canvas("c_s5g", 1, 0.75);
    gg->SetMarkerStyle(20);
    gf->SetMarkerStyle(24);
    auto* mg = new TMultiGraph();
    mg->Add(gg, "P");
    mg->Add(gf, "P");
    mg->SetTitle(";fecha (UTC);G mediana de hod2018 (mV/p.e.)");
    mg->Draw("A");
    mg->GetXaxis()->SetTimeDisplay(1);
    mg->GetXaxis()->SetTimeFormat("%d/%m%F1970-01-01 00:00:00");
    auto* lg2 = new TLegend(0.18, 0.20, 0.80, 0.35);
    lg2->AddEntry(gg, "G reajustada en el run", "p");
    lg2->AddEntry(gf, "sin ajuste: G de v8", "p");
    lg2->Draw();
    proc_save(c2, O + "/stage5_G_vs_date");
  }
}
