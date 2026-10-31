// plot_amp_profiles.C — seguimiento 2026-10-31 (3), punto 2: perfiles de amplitud (media de la resta mismo − mezclados, p.e.) contra la posición
// a lo largo de la barra, desde los CSV de amp_vs_position.py; curvas = ajustes exponencial (sólida) y lineal (trazos) de fits.csv.
// Salidas en <dir>: profiles_layers.png (XA, XB, YA, YB) y profiles_channels_<capa>.png (un panel por canal).
// Uso: root -l -b -q 'plot_amp_profiles.C("<dir del run>", "<etiqueta>")'
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<std::string> split(const std::string& l) {
  std::vector<std::string> v;
  std::stringstream s(l);
  std::string t;
  while (std::getline(s, t, ',')) v.push_back(t);
  return v;
}
double num(const std::string& s) { return s.empty() || s == "None" || s == "nan" ? NAN : std::stod(s); }
std::map<std::string, std::vector<double>> readFits(const std::string& fn) {
  std::map<std::string, std::vector<double>> m;
  std::ifstream f(fn);
  std::string l;
  std::getline(f, l);
  while (std::getline(f, l)) {
    auto v = split(l);
    std::vector<double> p;
    for (size_t i = 5; i < v.size(); ++i) p.push_back(num(v[i]));
    m[v[0]] = p;
  }
  return m;
}
void panel(const std::string& csv, const std::string& title, const std::vector<double>& p) {
  auto* g = new TGraphErrors();
  std::ifstream f(csv);
  std::string l;
  std::getline(f, l);
  while (std::getline(f, l)) {
    auto v = split(l);
    if (v.size() < 12) continue;
    const double y = num(v[10]), e = num(v[11]);
    if (!std::isfinite(y)) continue;
    const int n = g->GetN();
    g->SetPoint(n, num(v[1]), y);
    g->SetPointError(n, 0, e);
  }
  g->SetTitle((title + ";d desde el extremo del #acute{i}ndice 0 (mm);amplitud media (p.e.)").c_str());
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.7);
  g->SetMinimum(0);
  g->SetMaximum(45);
  g->GetXaxis()->SetLimits(19, 56);
  g->Draw("AP");
  if (p.size() >= 11 && std::isfinite(p[0])) {
    auto* fe = new TF1("", "[0]*exp([1]*(x-37.5))", 21, 54);
    fe->SetParameters(p[0], p[1]);
    fe->SetLineColor(kRed + 1);
    fe->Draw("same");
    auto* fl = new TF1("", "[0]+[1]*(x-37.5)", 21, 54);
    fl->SetParameters(p[9], p[10]);
    fl->SetLineColor(kBlue + 1);
    fl->SetLineStyle(2);
    fl->Draw("same");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.05);
    t.DrawLatex(0.15, 0.82, Form("L_{eff} = %.0f #pm %.0f mm (signo %+.0f), raz#acute{o}n %.2f #pm %.2f", p[2], p[3], p[4], p[5], p[6]));
  }
}
}  // namespace

void plot_amp_profiles(const char* dir, const char* label) {
  gStyle->SetOptStat(0);
  const std::string D(dir);
  auto fits = readFits(D + "/fits.csv");
  TCanvas c("c", "", 1200, 900);
  c.Divide(2, 2);
  int k = 1;
  for (const char* L : {"XA", "XB", "YA", "YB"}) {
    c.cd(k++);
    gPad->SetGrid();
    panel(D + "/profile_" + L + ".csv", std::string(label) + ", capa " + L + " (suma de canales)", fits[L]);
  }
  c.SaveAs((D + "/profiles_layers.png").c_str());
  for (const char* L : {"XA", "XB", "YA", "YB"}) {
    std::vector<std::pair<std::string, std::string>> chans;
    for (auto& [name, p] : fits)
      if (name.size() > 3 && name[0] == L[0] && name.find("_ch") != std::string::npos) {
        // capa del canal: índice par = A, impar = B
        const int idx = std::stoi(name.substr(1, 2));
        if ((idx % 2 == 0 ? 'A' : 'B') == L[1]) chans.push_back({name, name.substr(name.find("_ch") + 1)});
      }
    TCanvas cc("cc", "", 1600, 800);
    cc.Divide(4, 2);
    int i = 1;
    for (auto& [name, ch] : chans) {
      cc.cd(i++);
      gPad->SetGrid();
      panel(D + "/profile_" + std::string(1, L[0]) + "_" + ch + ".csv", std::string(label) + " " + name, fits[name]);
    }
    cc.SaveAs((D + "/profiles_channels_" + L + ".png").c_str());
  }
}
