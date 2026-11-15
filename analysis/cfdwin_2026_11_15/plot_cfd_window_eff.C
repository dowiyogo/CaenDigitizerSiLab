// plot_cfd_window_eff.C — seguimiento 2026-11-15: ε_XY de hod2018 contra el umbral T para la ventana de 144 ns (stage2) y las ventanas CFD de
// ±8 y ±16 ns (y, si están, las coincidencias por pares |t_X − t_Y| < 8, 16 ns, marcadores abiertos), una columna por campaña presente en el CSV
// (formato de cfd_window_eff.py, 16 columnas, o de cfd_window_eff2.py / eps_xy_cfd_window.csv, 17 columnas con noise_times). Barras internas: estadística; externas: estadística ⊕ estimador
// (la aceptancia, común a las tres ventanas, va en el texto). Figura de inspección.
// Uso: root -l -b -q 'analysis/cfdwin_2026_11_15/plot_cfd_window_eff.C("<csv>", "<prefijo de salida>")'
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

void plot_cfd_window_eff(const char* csv, const char* prefix) {
  gStyle->SetOptStat(0);
  struct P {
    double T, e, s, x;
  };
  std::map<std::string, std::map<std::string, std::vector<P>>> d;  // campaña -> ventana -> puntos
  std::vector<std::string> camps;
  std::ifstream f(csv);
  std::string l;
  while (std::getline(f, l)) {
    if (l.empty() || l[0] == '#' || l.rfind("campaign", 0) == 0) continue;
    std::vector<std::string> c;
    std::stringstream ss(l);
    std::string x;
    while (std::getline(ss, x, ',')) c.push_back(x);
    if (c.size() < 16) continue;
    const int o = c.size() >= 17 ? 1 : 0;  // columna noise_times en el formato de la versión 2
    if (o && c[3] == "uniform") continue;
    std::string win = c[2];
    if (win == "median8") win = "cfd8";
    if (win == "median16") win = "cfd16";
    if (!d.count(c[0])) camps.push_back(c[0]);
    d[c[0]][win].push_back({std::stod(c[1]), std::stod(c[12 + o]), std::stod(c[13 + o]), std::stod(c[14 + o])});
  }
  const std::vector<std::string> wins{"144ns", "cfd16", "cfd8", "pair16", "pair8"};
  const std::vector<std::string> labs{"144 ns signal window", "CFD window #pm16 ns around median t_{ref}", "CFD window #pm8 ns around median t_{ref}",
                                      "X-Y pair, #Deltat < 16 ns", "X-Y pair, #Deltat < 8 ns"};
  const std::vector<int> col{TColor::GetColor("#2a78d6"), TColor::GetColor("#eb6834"), TColor::GetColor("#1baf7a"), TColor::GetColor("#eb6834"),
                              TColor::GetColor("#1baf7a")};
  const std::vector<int> mk{20, 21, 33, 25, 27};
  const std::vector<double> dx{-0.3, -0.15, 0.0, 0.15, 0.3};
  const int nc = static_cast<int>(camps.size());
  TCanvas cv("cv", "", 700 * nc, 560);
  cv.Divide(nc, 1);
  std::vector<TObject*> keep;
  for (int ic = 0; ic < nc; ++ic) {
    cv.cd(ic + 1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gPad->SetGridy();
    auto* fr = new TH1F(Form("fr%d", ic), "", 10, 2, 11);
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.5);
    fr->GetXaxis()->SetTitle("threshold T (measured p.e.)");
    fr->GetYaxis()->SetTitle("#varepsilon_{XY}");
    fr->Draw();
    keep.push_back(fr);
    auto* lg = new TLegend(0.50, 0.56, 0.89, 0.83);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->SetHeader(camps[ic].c_str());
    for (size_t iw = 0; iw < wins.size(); ++iw) {
      const auto& v = d[camps[ic]][wins[iw]];
      if (v.empty()) continue;
      auto* go = new TGraphErrors(static_cast<int>(v.size()));
      auto* gi = new TGraphErrors(static_cast<int>(v.size()));
      for (size_t i = 0; i < v.size(); ++i) {
        go->SetPoint(i, v[i].T + dx[iw], v[i].e);
        go->SetPointError(i, 0, std::hypot(v[i].s, v[i].x));
        gi->SetPoint(i, v[i].T + dx[iw], v[i].e);
        gi->SetPointError(i, 0, v[i].s);
      }
      for (auto* g : {go, gi}) {
        g->SetLineColor(col[iw]);
        g->SetMarkerColor(col[iw]);
        g->SetMarkerStyle(mk[iw]);
        g->SetMarkerSize(1.2);
        g->SetLineWidth(2);
      }
      go->SetLineStyle(1);
      go->SetLineWidth(1);
      go->Draw("P Z same");
      gi->Draw("P same");
      lg->AddEntry(gi, labs[iw].c_str(), "lp");
      keep.push_back(go);
      keep.push_back(gi);
    }
    lg->Draw();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.03);
    t.DrawLatex(0.15, 0.86, "inner bar: stat; outer: stat #oplus estimator; acceptance 4.65 % common");
  }
  cv.SaveAs((std::string(prefix) + ".png").c_str());
  cv.SaveAs((std::string(prefix) + ".pdf").c_str());
}
