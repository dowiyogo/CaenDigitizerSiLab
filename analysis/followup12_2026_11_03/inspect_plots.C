// inspect_plots.C — seguimiento 2026-11-03 (12): figuras de inspección (no para el artículo) desde las tablas de h19_hvgroup.py y diurnal.py.
//  A: (1) mayor pulso de hod2018 en triggers con hod2019 en tiempo (> 35 y > 3 mV) contra todos, normalizado, mayo (mapa de la planilla); (2) mayor
//     pulso dentro del mismo grupo HV, pares con hod2019 contra triggers sin hod2019; (3) Δk mismo grupo / otros grupos; (4) DR por grupo HV.
//  B: medianas de la estación por bin de 1 h (G/G_ref, DCR/DCR_ref, p_ct/p_ct,ref) contra la hora local, mayo y octubre; bins buenos llenos.
// Uso (desde la raíz): root -l -b -q 'analysis/followup12_2026_11_03/inspect_plots.C("<dir A>", "<dir B>")'
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TStyle.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<std::map<std::string, std::string>> readCsv(const std::string& p) {
  std::ifstream in(p);
  std::string line;
  std::vector<std::string> hdr;
  std::vector<std::map<std::string, std::string>> rows;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    std::vector<std::string> v;
    std::string tok;
    while (std::getline(ss, tok, ',')) v.push_back(tok);
    if (hdr.empty()) {
      hdr = v;
      continue;
    }
    std::map<std::string, std::string> r;
    for (size_t i = 0; i < hdr.size() && i < v.size(); ++i) r[hdr[i]] = v[i];
    rows.push_back(r);
  }
  return rows;
}
double num(const std::string& x) {
  try {
    return std::stod(x);
  } catch (...) {
    printf("valor no numérico: '%s'\n", x.c_str());
    return NAN;
  }
}
}  // namespace

void inspect_plots(const char* dirA, const char* dirB) {
  gStyle->SetOptStat(0);
  const std::string A(dirA), B(dirB);
  auto mx = readCsv(A + "/h19_hvgroup_maxamp.csv");
  auto dk = readCsv(A + "/h19_hvgroup_dk.csv");
  auto bg = readCsv(A + "/h19_hvgroup_by_group.csv");
  auto* c = new TCanvas("cA", "", 1400, 1000);
  c->Divide(2, 2);
  int col[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
  // (1)
  c->cd(1);
  gPad->SetLogy();
  {
    auto* lg = new TLegend(0.45, 0.65, 0.88, 0.88);
    int k = 0;
    for (const char* T : {"35", "3"}) {
      auto* ha = new TH1D(Form("ha%s", T), "May run B: largest hod2018 pulse;max_{j} a_{j}/G_{j} (p.e.);fraction / 0.5 p.e.", 200, 0, 100);
      auto* hh = new TH1D(Form("hh%s", T), "", 200, 0, 100);
      for (auto& r : mx)
        if (r["campaign"] == "may" && r["map"] == "sheet" && r["T_mV"] == T) {
          const double x = num(r["bin_lo_pe"]) + 0.25;
          ha->Fill(x, num(r["all_triggers"]));
          hh->Fill(x, num(r["triggers_with_h19"]));
        }
      ha->Scale(1 / ha->Integral());
      hh->Scale(1 / hh->Integral());
      ha->GetXaxis()->SetRangeUser(0, 60);
      ha->SetLineColor(kGray + 2);
      hh->SetLineColor(std::string(T) == "35" ? kRed + 1 : kBlue + 1);
      hh->SetLineWidth(2);
      if (k == 0) {
        ha->Draw("HIST");
        lg->AddEntry(ha, "all triggers", "l");
      }
      hh->Draw("HIST SAME");
      lg->AddEntry(hh, Form("hod2019 in time > %s mV", T), "l");
      ++k;
    }
    lg->Draw();
  }
  // (2)
  c->cd(2);
  gPad->SetLogy();
  {
    auto* lg = new TLegend(0.40, 0.65, 0.88, 0.88);
    int k = 0;
    for (const char* T : {"35", "3"}) {
      auto* ha = new TH1D(Form("sa%s", T), "May run B, sheet map: largest hod2018 pulse in the SAME HV group;p.e.;fraction / 0.5 p.e.", 200, 0, 100);
      auto* hb = new TH1D(Form("sb%s", T), "", 200, 0, 100);
      for (auto& r : mx)
        if (r["campaign"] == "may" && r["map"] == "sheet" && r["T_mV"] == T) {
          const double x = num(r["bin_lo_pe"]) + 0.25;
          ha->Fill(x, num(r["same_group_A_pairs"]));
          hb->Fill(x, num(r["same_group_B_pairs"]));
        }
      ha->Scale(1 / ha->Integral());
      hb->Scale(1 / hb->Integral());
      ha->GetXaxis()->SetRangeUser(0, 40);
      ha->SetLineColor(std::string(T) == "35" ? kRed + 1 : kBlue + 1);
      ha->SetLineWidth(2);
      hb->SetLineColor(std::string(T) == "35" ? kRed + 1 : kBlue + 1);
      hb->SetLineWidth(2);
      hb->SetLineStyle(2);
      ha->Draw(k == 0 ? "HIST" : "HIST SAME");
      hb->Draw("HIST SAME");
      lg->AddEntry(ha, Form("with hod2019 > %s mV", T), "l");
      lg->AddEntry(hb, Form("no hod2019 > %s mV", T), "l");
      ++k;
    }
    lg->Draw();
  }
  // (3)
  c->cd(3);
  {
    auto* lg = new TLegend(0.12, 0.68, 0.50, 0.88);
    int k = 0;
    for (const char* T : {"35", "3"}) {
      auto* hs = new TH1D(Form("ds%s", T), "May run B, sheet map: #Deltak = k_{hod2019} - k_{hod2018} (> 5 p.e.);#Deltak (samples of 16 ns);fraction",
                          25, -12.5, 12.5);
      auto* ho = new TH1D(Form("do%s", T), "", 25, -12.5, 12.5);
      for (auto& r : dk)
        if (r["campaign"] == "may" && r["map"] == "sheet" && r["T_mV"] == T) {
          hs->Fill(num(r["dk_samples"]), num(r["pairs_same_group"]));
          ho->Fill(num(r["dk_samples"]), num(r["pairs_other_groups"]));
        }
      hs->Scale(1 / hs->Integral());
      ho->Scale(1 / ho->Integral());
      hs->SetMaximum(0.3);
      hs->SetLineColor(std::string(T) == "35" ? kRed + 1 : kBlue + 1);
      hs->SetLineWidth(2);
      ho->SetLineColor(std::string(T) == "35" ? kRed + 1 : kBlue + 1);
      ho->SetLineWidth(2);
      ho->SetLineStyle(2);
      hs->Draw(k == 0 ? "HIST" : "HIST SAME");
      ho->Draw("HIST SAME");
      lg->AddEntry(hs, Form("same group, > %s mV", T), "l");
      lg->AddEntry(ho, Form("other groups, > %s mV", T), "l");
      ++k;
    }
    lg->Draw();
  }
  // (4)
  c->cd(4);
  {
    auto* mg = new TMultiGraph();
    mg->SetTitle("May run B, in time: DR = R_{same}/R_{other} by HV group;HV group (+ offset per map / threshold);DR");
    auto* lg = new TLegend(0.12, 0.70, 0.60, 0.88);
    int k = 0;
    for (const char* m : {"sheet", "oct_measured", "may_measured_g12"})
      for (const char* T : {"35", "3"}) {
        auto* g = new TGraphErrors();
        for (auto& r : bg)
          if (r["campaign"] == "may" && r["map"] == m && r["T_mV"] == T && r["window"] == "S" && r["DR"] != "nan") {
            const int n = g->GetN();
            g->SetPoint(n, num(r["hv_group"]) + 0.12 * (k - 2.5), num(r["DR"]));
            g->SetPointError(n, 0, num(r["DR_err"]));
          }
        g->SetMarkerStyle(std::string(T) == "35" ? 20 : 24);
        g->SetMarkerColor(col[1 + k / 2]);
        g->SetLineColor(col[1 + k / 2]);
        mg->Add(g, "P");
        lg->AddEntry(g, Form("%s, > %s mV", m, T), "p");
        ++k;
      }
    mg->Draw("A");
    mg->SetMinimum(0);
    mg->SetMaximum(3);
    auto* l = new TLine(0.5, 1, 8.5, 1);
    l->SetLineStyle(2);
    l->Draw();
    lg->SetNColumns(2);
    lg->Draw();
  }
  c->SaveAs((A + "/h19_hvgroup_inspect.png").c_str());
  c->SaveAs((A + "/h19_hvgroup_inspect.pdf").c_str());

  auto bins = readCsv(B + "/diurnal_bins.csv");
  auto* d = new TCanvas("cB", "", 1400, 1000);
  d->Divide(2, 3);
  int pad = 1;
  for (const char* q : {"G_rel", "DCR_rel", "pct_rel"})
    for (const char* camp : {"may", "oct"}) {
      d->cd(pad++);
      gPad->SetGridx();
      auto* gg = new TGraphErrors();
      auto* gb = new TGraphErrors();
      double t0 = -1;
      for (auto& r : bins) {
        if (r["campaign"] != camp || r[std::string(q) + "_median"] == "nan") continue;
        const double tl = num(r["t_since_first_trigger_h"]);
        const double loc = num(r["local_h_mid"]);
        if (t0 < 0) t0 = loc - tl;  // hora local del primer trigger
        auto* g = r["good_bin"] == "1" ? gg : gb;
        const int n = g->GetN();
        g->SetPoint(n, tl + t0, num(r[std::string(q) + "_median"]));
        g->SetPointError(n, 0, num(r[std::string(q) + "_median_err"]));
      }
      auto* mg = new TMultiGraph();
      mg->SetTitle(Form("%s %s (station median per 1 h bin);local time since 00:00 of the first day (h);%s", camp, q, q));
      gg->SetMarkerStyle(20);
      gg->SetMarkerSize(0.7);
      gb->SetMarkerStyle(24);
      gb->SetMarkerSize(0.7);
      gb->SetMarkerColor(kGray + 2);
      gb->SetLineColor(kGray + 2);
      if (gg->GetN()) mg->Add(gg, "P");
      if (gb->GetN()) mg->Add(gb, "P");
      mg->Draw("A");
      for (double x = 0; x < 100; x += 24) {
        auto* l = new TLine(x, mg->GetYaxis()->GetXmin(), x, mg->GetYaxis()->GetXmax());
        l->SetLineStyle(3);
        if (x > mg->GetXaxis()->GetXmin() && x < mg->GetXaxis()->GetXmax()) l->Draw();
      }
    }
  d->SaveAs((B + "/diurnal_inspect.png").c_str());
  d->SaveAs((B + "/diurnal_inspect.pdf").c_str());
}
