// plot_imon_24h.C — seguimiento 2026-10-31, 1b: IMon de CH00-07 (medianas por minuto, imon_24h_1min.csv de las últimas 24 h) contra la
// hora UTC, con CH03 marcado (línea gruesa roja), la línea de I0Set = 15 µA y el instante del trip de CH03 (02:08:30 UTC del 2026-10-30).
// Segundo panel: CH03 a 1 s alrededor del barrido (ch03_scan_1s.csv), IMon y VMon.
// Uso: root -l -b -q 'plot_imon_24h.C("<dir>")'
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void plot_imon_24h(const char* dir) {
  const std::string D(dir);
  gStyle->SetOptStat(0);
  const int col[8] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1, kMagenta + 1, kOrange + 7, kCyan + 2, kGray + 2};
  std::vector<TGraph*> g(8);
  for (int c = 0; c < 8; ++c) g[c] = new TGraph();
  std::ifstream f(D + "/imon_24h_1min.csv");
  std::string l;
  while (std::getline(f, l)) {
    if (l.empty() || l[0] == '#' || l[0] == 'm') continue;
    std::stringstream s(l);
    std::string a, b, c;
    std::getline(s, a, ',');
    std::getline(s, b, ',');
    std::getline(s, c, ',');
    g[std::stoi(b)]->AddPoint(std::stod(a), std::stod(c));
  }
  TCanvas cv("cv", "", 1000, 900);
  cv.Divide(1, 2);
  cv.cd(1)->SetGrid();
  g[0]->SetTitle("IMon de CH00-07, #acute{u}ltimas 24 h (mediana por minuto);hora UTC;IMon (#muA)");
  g[0]->GetXaxis()->SetTimeDisplay(1);
  g[0]->GetXaxis()->SetTimeFormat("%d %H:%M");
  g[0]->GetXaxis()->SetTimeOffset(0, "gmt");
  g[0]->SetMinimum(-1);
  g[0]->SetMaximum(17);
  TLegend lg(0.12, 0.62, 0.32, 0.88);
  lg.SetNColumns(2);
  for (int c = 0; c < 8; ++c) {
    g[c]->SetLineColor(col[c]);
    g[c]->SetLineWidth(c == 3 ? 3 : 1);
    g[c]->Draw(c == 0 ? "AL" : "L same");
    lg.AddEntry(g[c], Form(c == 3 ? "CH%02d (marcado)" : "CH%02d", c), "l");
  }
  const double x0 = g[0]->GetXaxis()->GetXmin(), x1 = g[0]->GetXaxis()->GetXmax();
  TLine li(x0, 15, x1, 15);
  li.SetLineStyle(2);
  li.SetLineColor(kRed);
  li.Draw();
  TLine lt(1793326110, -1, 1793326110, 17);
  lt.SetLineStyle(3);
  lt.Draw();
  lg.Draw();
  cv.cd(2)->SetGrid();
  auto* gi = new TGraph();
  auto* gv = new TGraph();
  std::ifstream f2(D + "/ch03_scan_1s.csv");
  std::getline(f2, l);
  while (std::getline(f2, l)) {
    std::stringstream s(l);
    std::string t, u, i, v;
    std::getline(s, t, ',');
    std::getline(s, u, ',');
    std::getline(s, i, ',');
    std::getline(s, v, ',');
    gi->AddPoint(std::stod(t), std::stod(i));
    gv->AddPoint(std::stod(t), (std::stod(v) - 60.0) * 2.0);
  }
  gi->SetTitle("CH03 alrededor del barrido, 1 s: IMon (rojo, #muA) y (VMon #minus 60 V)#times2 (azul);hora UTC del 2026-10-30;IMon (#muA)");
  gi->GetXaxis()->SetTimeDisplay(1);
  gi->GetXaxis()->SetTimeFormat("%H:%M");
  gi->GetXaxis()->SetTimeOffset(0, "gmt");
  gi->GetXaxis()->SetLimits(1793325900, 1793326320);
  gi->SetMinimum(-1);
  gi->SetMaximum(17);
  gi->SetLineColor(kRed + 1);
  gi->SetLineWidth(2);
  gi->Draw("AL");
  gv->SetLineColor(kBlue + 1);
  gv->Draw("L same");
  TLine li2(1793325900, 15, 1793326320, 15);
  li2.SetLineStyle(2);
  li2.SetLineColor(kRed);
  li2.Draw();
  cv.SaveAs((D + "/imon_24h.png").c_str());
}
