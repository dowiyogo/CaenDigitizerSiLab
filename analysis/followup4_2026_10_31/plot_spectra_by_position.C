// plot_spectra_by_position.C — seguimiento 2026-10-31 (4), punto 2: espectros de coincidencias X–Y verdaderas (mismo evento − mezclados) por capa
// y por distancia al extremo de lectura, desde spectrum_<capa>_pos<k>.csv (bins de 1 p.e.); se muestran reagrupados a 4 p.e. hasta 120 p.e.
// Un PNG por capa: spectrum_<capa>.png. Uso: root -l -b -q 'plot_spectra_by_position.C("<dir>")'
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

void plot_spectra_by_position(const char* dir) {
  gStyle->SetOptStat(0);
  const char* lab[4] = {"21-29 mm", "29-37 mm", "37-45 mm", "45-54 mm"};
  const int col[4] = {kBlue + 1, kGreen + 2, kOrange + 7, kRed + 1};
  for (const char* L : {"XA", "XB", "YA", "YB"}) {
    TCanvas c("c", "", 900, 650);
    c.SetGrid();
    TLegend lg(0.60, 0.62, 0.89, 0.88);
    lg.SetHeader(Form("capa %s, distancia al extremo de lectura", L));
    TH1D* hs[4];
    double ymax = 0;
    for (int k = 0; k < 4; ++k) {
      hs[k] = new TH1D(Form("h%s%d", L, k), "", 30, 0, 120);
      std::ifstream f(std::string(dir) + "/spectrum_" + L + "_pos" + std::to_string(k) + ".csv");
      std::string l;
      while (std::getline(f, l)) {
        if (l.empty() || l[0] == '#' || l[0] == 'p') continue;
        std::stringstream s(l);
        std::string a, b, ns, nm, nn, ne;
        std::getline(s, a, ','); std::getline(s, b, ','); std::getline(s, ns, ','); std::getline(s, nm, ','); std::getline(s, nn, ','); std::getline(s, ne, ',');
        const double x = std::stod(a) + 0.5;
        if (x >= 120) continue;
        const int bin = hs[k]->FindBin(x);
        hs[k]->SetBinContent(bin, hs[k]->GetBinContent(bin) + std::stod(nn));
        hs[k]->SetBinError(bin, std::hypot(hs[k]->GetBinError(bin), std::stod(ne)));
      }
      hs[k]->SetLineColor(col[k]);
      hs[k]->SetMarkerColor(col[k]);
      hs[k]->SetMarkerStyle(20 + k);
      hs[k]->SetMarkerSize(0.7);
      ymax = std::max(ymax, hs[k]->GetMaximum());
      lg.AddEntry(hs[k], lab[k], "lp");
    }
    hs[0]->SetTitle(Form("hod2018 octubre run 8, coincidencias X-Y verdaderas, capa %s;amplitud (p.e., bins de 4);impactos netos (mismo #minus mezclados)", L));
    hs[0]->SetMaximum(1.25 * ymax);
    hs[0]->SetMinimum(-10);
    hs[0]->Draw("E1");
    for (int k = 1; k < 4; ++k) hs[k]->Draw("E1 same");
    lg.Draw();
    c.SaveAs((std::string(dir) + "/spectrum_" + L + ".png").c_str());
  }
}
