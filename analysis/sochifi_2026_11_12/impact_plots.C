// impact_plots.C — SOCHIFI 2026 (póster), 2026-11-12: histogramas ROOT y figuras de las proyecciones de impacto de impact_map.py.
// Lee <dir>/<segmento>/impact_map_runB_T{5,10}.csv; escribe <dir>/<segmento>/impact_map_runB.root (por umbral y eje: en tiempo, accidentales,
// restadas, por trigger, por mm; bordes variables de las 31 celdas) y <dir>/<segmento>/impacto_{x,y}.{pdf,png}: cuentas restadas por mm contra
// la posición, 5 y 10 p.e., celdas de ch3 sombreadas en x. Sin corrección de eficiencia ni aceptancia.
// Uso (desde la raíz): root -l -b -q 'analysis/sochifi_2026_11_12/impact_plots.C("<dir de salida>", "<segmento>")'
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
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
}  // namespace

void impact_plots(const char* dir, const char* seg) {
  gStyle->SetOptStat(0);
  const std::string base = std::string(dir) + "/" + seg;
  std::vector<double> edges{-16.5};
  for (int k = 0; k < 30; ++k) edges.push_back(-14.5 + k);
  edges.push_back(16.5);
  TFile fo((base + "/impact_map_runB.root").c_str(), "RECREATE");
  std::map<std::string, TH1D*> perMm;
  long ntrig = 0;
  for (int T : {5, 10}) {
    auto rows = readCsv(base + Form("/impact_map_runB_T%d.csv", T));
    for (const char* ax : {"x", "y"}) {
      auto mk = [&](const char* what, const char* title) {
        auto* h = new TH1D(Form("h_%s_%s_T%d", what, ax, T), Form("%s, eje %s, %d p.e.;posici#acute{o}n %s (mm);%s", seg, ax, T, ax, title), 31, edges.data());
        h->SetDirectory(&fo);
        return h;
      };
      auto* hOn = mk("on_time", "cuentas en tiempo");
      auto* hAcc = mk("accidental", "accidentales estimados");
      auto* hRaw = mk("raw", "cuentas restadas");
      auto* hTrig = mk("per_trigger", "cuentas restadas por trigger");
      auto* hMm = mk("per_mm", "cuentas restadas por mm");
      for (auto& r : rows) {
        if (r["axis"] != ax) continue;
        const int b = std::stoi(r["cell_index"]) + 1;
        hOn->SetBinContent(b, std::stod(r["on_time_count"]));
        hOn->SetBinError(b, std::sqrt(std::stod(r["on_time_count"])));
        hAcc->SetBinContent(b, std::stod(r["accidental_estimate"]));
        hRaw->SetBinContent(b, std::stod(r["raw_count"]));
        hRaw->SetBinError(b, std::stod(r["stat_error_count"]));
        hTrig->SetBinContent(b, std::stod(r["counts_per_trigger"]));
        hTrig->SetBinError(b, std::stod(r["stat_error_per_trigger"]));
        hMm->SetBinContent(b, std::stod(r["counts_per_mm"]));
        hMm->SetBinError(b, std::stod(r["counts_per_mm_err"]));
        ntrig = std::stol(r["triggers_denominator"]);
      }
      perMm[Form("%s_%d", ax, T)] = hMm;
    }
  }
  fo.Write();
  for (const char* ax : {"x", "y"}) {
    auto* c = new TCanvas(Form("c_%s", ax), "", 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);
    auto* h5 = perMm[Form("%s_5", ax)];
    auto* h10 = perMm[Form("%s_10", ax)];
    double ymax = 0;
    for (auto* h : {h5, h10})
      for (int b = 1; b <= 31; ++b) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    h5->SetTitle(Form(";posici#acute{o}n %s (mm, centrada);cuentas restadas por mm", ax));
    h5->SetMinimum(0);
    h5->SetMaximum(1.45 * ymax);
    h5->SetLineColor(kBlue + 1);
    h5->SetMarkerColor(kBlue + 1);
    h5->SetMarkerStyle(20);
    h5->SetLineWidth(2);
    h10->SetLineColor(kRed + 1);
    h10->SetMarkerColor(kRed + 1);
    h10->SetMarkerStyle(24);
    h10->SetLineWidth(2);
    h5->GetXaxis()->SetTitleSize(0.045);
    h5->GetYaxis()->SetTitleSize(0.045);
    h5->Draw("E1");
    if (std::string(ax) == "x") {  // ch3: celdas 23-25 (7.5 a 10.5 mm) sin cuentas posibles; 22 y 26 (6.5-7.5, 10.5-11.5) vecinas afectadas
      auto* bx = new TBox(7.5, 0, 10.5, 1.3 * ymax);
      bx->SetFillColorAlpha(kGray + 1, 0.5);
      bx->Draw();
      for (double lo : {6.5, 10.5}) {
        auto* bn = new TBox(lo, 0, lo + 1.0, 1.3 * ymax);
        bn->SetFillColorAlpha(kGray, 0.3);
        bn->Draw();
      }
    }
    h5->Draw("E1 SAME");
    h10->Draw("E1 SAME");
    auto* lg = new TLegend(0.15, 0.72, 0.70, 0.89);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.033);
    lg->SetHeader(Form("hod2018, run B, %ld triggers, accidentales restados", ntrig));
    lg->AddEntry(h5, "umbral 5 p.e.", "lp");
    lg->AddEntry(h10, "umbral 10 p.e.", "lp");
    if (std::string(ax) == "x") lg->AddEntry((TObject*)nullptr, "gris: canal 3 muerto (claro: celdas vecinas)", "");
    lg->Draw();
    c->SaveAs((base + Form("/impacto_%s.pdf", ax)).c_str());
    c->SaveAs((base + Form("/impacto_%s.png", ax)).c_str());
  }
}
