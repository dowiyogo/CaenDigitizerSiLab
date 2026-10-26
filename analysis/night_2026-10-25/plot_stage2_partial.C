// plot_stage2_partial.C — XPLUS-4 etapa 2 (parcial): figuras de la tasa oscura de hod2018 desde los CSV de stage2_partial.py.
//  (a) escalera: tasa mediana sobre umbral (en mV fijos = T·G a 0 V) para ΔV = 0, −0.5, −1.0, −1.5 V;
//  (b) R(> 0.5 p.e. a 0 V) y tasa sobre el piso de 1 mV contra ΔV (mediana de los canales);
//  (c) DCR mediana de hod2018 por hora en el run 8 de octubre.
// Uso: root -l -b -q daq/tools/proceeding/proc_style.C 'analysis/night_2026-10-25/plot_stage2_partial.C("<dir>")'
#include <TAxis.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

void plot_stage2_partial(const char* dir) {
  proc_style();
  const std::string d(dir);
  // (a) escalera
  std::map<double, TGraph*> st;
  std::ifstream a(d + "/staircase_median.csv");
  for (std::string l; std::getline(a, l);) {
    if (l.empty() || l[0] == '#') continue;
    std::istringstream s(l);
    double dv, t, r;
    s >> dv >> t >> r;
    if (r <= 0) continue;
    auto*& g = st[dv];
    if (!g) g = new TGraph();
    g->SetPoint(g->GetN(), t, r);
  }
  auto* c1 = proc_canvas("c_st", 1, 0.8);
  c1->SetLogy();
  auto* mg = new TMultiGraph();
  auto* lg = new TLegend(0.18, 0.20, 0.50, 0.40);
  int k = 0;
  for (auto it = st.rbegin(); it != st.rend(); ++it, ++k) {
    it->second->SetMarkerStyle(20 + k);
    it->second->SetMarkerColor(proc::color(k == 0 ? 0 : k + 4));
    it->second->SetLineColor(proc::color(k == 0 ? 0 : k + 4));
    mg->Add(it->second, "PL");
    lg->AddEntry(it->second, Form("#DeltaV = %+.1f V", it->first), "pl");
  }
  mg->SetTitle(";umbral (p.e. de G a 0 V, en mV fijos);tasa mediana de hod2018 (s^{-1})");
  mg->Draw("A");
  lg->Draw();
  proc_save(c1, d + "/staircase_scan_h18");
  // (b) contra ΔV
  std::map<double, std::vector<double>> r05, rfl;
  std::ifstream b(d + "/dcr_vs_dV.csv");
  for (std::string l; std::getline(b, l);) {
    if (l.empty() || l[0] == '#') continue;
    std::istringstream s(l);
    int ch;
    double dv, x, y;
    s >> ch >> dv >> x >> y;
    r05[dv].push_back(x);
    rfl[dv].push_back(y);
  }
  auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v.empty() ? 0.0 : v[v.size() / 2]; };
  auto* g05 = new TGraph();
  auto* gfl = new TGraph();
  for (auto& [dv, v] : r05) {
    if (med(v) > 0) g05->SetPoint(g05->GetN(), dv, med(v));
    gfl->SetPoint(gfl->GetN(), dv, med(rfl[dv]));
  }
  auto* c2 = proc_canvas("c_dv", 1, 0.8);
  c2->SetLogy();
  g05->SetMarkerStyle(20);
  gfl->SetMarkerStyle(24);
  gfl->SetMarkerColor(proc::color(6));
  gfl->SetLineColor(proc::color(6));
  auto* mg2 = new TMultiGraph();
  mg2->Add(g05, "PL");
  mg2->Add(gfl, "PL");
  mg2->SetTitle(";#DeltaV respecto del nominal (V);tasa mediana de hod2018 (s^{-1})");
  mg2->Draw("A");
  auto* lg2 = new TLegend(0.42, 0.20, 0.93, 0.38);
  lg2->AddEntry(g05, "> 0.5 p.e. (G a 0 V)", "pl");
  lg2->AddEntry(gfl, "> 1 mV (piso del buscador)", "pl");
  lg2->Draw();
  proc_save(c2, d + "/dcr_vs_dV_h18");
  // (c) por hora, run 8
  auto* gh = new TGraph();
  std::ifstream h(d + "/run8_hourly_dcr.csv");
  for (std::string l; std::getline(h, l);) {
    if (l.empty() || l[0] == '#') continue;
    std::istringstream s(l);
    double t, r;
    s >> t >> r;
    gh->SetPoint(gh->GetN(), t, r / 1e6);
  }
  auto* c3 = proc_canvas("c_h", 1, 0.7);
  gh->SetMarkerStyle(20);
  gh->SetTitle(";hora (UTC);DCR mediana de hod2018 (10^{6} s^{-1})");
  gh->Draw("APL");
  gh->GetXaxis()->SetTimeDisplay(1);
  gh->GetXaxis()->SetTimeFormat("%H:%M%F1970-01-01 00:00:00");
  proc_save(c3, d + "/run8_hourly_dcr_h18");
}
