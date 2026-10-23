// charact_plots.C — DAQ2-15, parte B: figuras PNG de la caracterización desde work/figs.json (charact_analyse.py).
// Tipos: "points" / "lines" (series con x numérica o etiquetas de texto, yerr opcional, logy), "map" (matriz z con
// etiquetas de ejes, COLZ con el valor en cada celda) y "grid" (16 paneles de histogramas, uno o dos por panel).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'daq/tools/characterization/charact_plots.C+("<work>/figs.json", "<dir de salida>")'
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace {
const int kColors[] = {1, 2, 4, 8, 6, 7, 9, 28, 30, 38, 46, 41};

double num(const nlohmann::json& v) { return v.is_number() ? v.get<double>() : NAN; }

void drawSeries(const nlohmann::json& f, const std::string& out) {
  TCanvas c("c", "", 1100, 650);
  c.SetGrid();
  if (f.value("logy", false)) c.SetLogy();
  const auto& S = f["series"];
  // eje x con etiquetas si la primera serie tiene x de texto
  bool text = !S.empty() && !S[0]["x"].empty() && S[0]["x"][0].is_string();
  std::vector<std::string> labels;
  if (text)
    for (const auto& x : S[0]["x"]) labels.push_back(x.get<std::string>());
  TMultiGraph mg;
  TLegend lg(0.80, 0.55, 0.995, 0.995);
  lg.SetTextSize(0.022);
  int k = 0;
  for (const auto& s : S) {
    auto* g = new TGraphErrors();
    for (std::size_t i = 0; i < s["x"].size(); ++i) {
      double x = NAN;
      if (text) {
        const std::string lab = s["x"][i].is_string() ? s["x"][i].get<std::string>() : "";
        for (std::size_t j = 0; j < labels.size(); ++j)
          if (labels[j] == lab) x = static_cast<double>(j);
        if (std::isnan(x)) {
          labels.push_back(lab);
          x = static_cast<double>(labels.size() - 1);
        }
      } else {
        x = num(s["x"][i]);
      }
      const double y = num(s["y"][i]);
      if (std::isnan(x) || std::isnan(y) || (f.value("logy", false) && y <= 0)) continue;
      const int n = g->GetN();
      g->SetPoint(n, x, y);
      if (s.contains("yerr") && s["yerr"].is_array() && i < s["yerr"].size()) g->SetPointError(n, 0, num(s["yerr"][i]) == num(s["yerr"][i]) ? num(s["yerr"][i]) : 0);
    }
    const int col = kColors[k % 12];
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->SetMarkerStyle(20 + (k / 12) % 4);
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2);
    if (g->GetN() > 0) {
      mg.Add(g, f["kind"] == "lines" ? "L" : "PE");
      lg.AddEntry(g, s["name"].get<std::string>().c_str(), f["kind"] == "lines" ? "l" : "pe");
    }
    ++k;
  }
  const std::string title = f["title"].get<std::string>() + ";" + f["xlabel"].get<std::string>() + ";" + f["ylabel"].get<std::string>();
  TH1D* frame = nullptr;
  if (text && mg.GetListOfGraphs()) {
    // marco con etiquetas de texto en x; el rango y sale de los puntos
    double lo = 1e300, hi = -1e300;
    for (auto* o : *mg.GetListOfGraphs()) {
      auto* g = static_cast<TGraphErrors*>(o);
      for (int i = 0; i < g->GetN(); ++i) {
        lo = std::min(lo, g->GetY()[i] - g->GetEY()[i]);
        hi = std::max(hi, g->GetY()[i] + g->GetEY()[i]);
      }
    }
    const bool logy = f.value("logy", false);
    frame = new TH1D("frame", title.c_str(), static_cast<int>(labels.size()), -0.5, static_cast<double>(labels.size()) - 0.5);
    const int step = std::max<int>(1, static_cast<int>(labels.size()) / 30);
    for (std::size_t j = 0; j < labels.size(); j += static_cast<std::size_t>(step)) frame->GetXaxis()->SetBinLabel(static_cast<int>(j) + 1, labels[j].c_str());
    frame->GetXaxis()->LabelsOption("v");
    frame->SetMinimum(logy ? std::max(lo * 0.5, 1e-6) : lo - 0.05 * (hi - lo));
    frame->SetMaximum(logy ? hi * 2 : hi + 0.05 * (hi - lo));
    frame->SetStats(false);
    frame->Draw("AXIS");
    mg.Draw("");
  } else {
    mg.SetTitle(title.c_str());
    mg.Draw("A");
  }
  if (S.size() <= 14) lg.Draw();
  c.SaveAs(out.c_str());
  delete frame;
}

void drawMap(const nlohmann::json& f, const std::string& out) {
  const auto& z = f["z"];
  const int nx = static_cast<int>(z.size()), ny = static_cast<int>(z[0].size());
  TH2D h("h", (f["title"].get<std::string>() + ";" + f["xlabel"].get<std::string>() + ";" + f["ylabel"].get<std::string>()).c_str(), nx, -0.5, nx - 0.5, ny,
         -0.5, ny - 0.5);
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) h.SetBinContent(i + 1, j + 1, num(z[i][j]));
  for (int i = 0; i < nx && i < static_cast<int>(f["xlabels"].size()); ++i) h.GetXaxis()->SetBinLabel(i + 1, (std::to_string(i) + ":" + f["xlabels"][i].get<std::string>()).c_str());
  for (int j = 0; j < ny && j < static_cast<int>(f["ylabels"].size()); ++j) h.GetYaxis()->SetBinLabel(j + 1, (std::to_string(j) + ":" + f["ylabels"][j].get<std::string>()).c_str());
  TCanvas c("c", "", 1000, 850);
  c.SetRightMargin(0.14);
  gStyle->SetPaintTextFormat(".0f");
  h.SetMarkerSize(0.7);
  h.GetXaxis()->LabelsOption("v");
  h.Draw("COLZ TEXT");
  c.SaveAs(out.c_str());
}

void drawGrid(const nlohmann::json& f, const std::string& out) {
  TCanvas c("c", "", 1600, 1200);
  c.Divide(4, 4);
  std::vector<TH1D*> keep;
  int p = 0;
  for (const auto& pad : f["pads"]) {
    c.cd(++p);
    if (f.value("logy", false)) gPad->SetLogy();
    const auto& x = pad["x"];
    const int n = static_cast<int>(x.size());
    const double w = n > 1 ? num(x[1]) - num(x[0]) : 1;
    int k = 0;
    for (const auto& ys : pad["ys"]) {
      auto* h = new TH1D(Form("h%d_%d", p, k), (pad["title"].get<std::string>() + ";" + f["xlabel"].get<std::string>() + ";" + f["ylabel"].get<std::string>()).c_str(), n,
                         num(x[0]) - w / 2, num(x[n - 1]) + w / 2);
      for (int i = 0; i < n; ++i) h->SetBinContent(i + 1, num(ys[i]));
      h->SetLineColor(k == 0 ? 1 : 2);
      h->SetStats(false);
      h->SetMinimum(0.5);
      h->Draw(k == 0 ? "HIST" : "HIST SAME");
      keep.push_back(h);
      ++k;
    }
  }
  c.SaveAs(out.c_str());
  for (auto* h : keep) delete h;
}
}  // namespace

void charact_plots(const char* figsJson, const char* outDir) {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.03);
  std::ifstream in(figsJson);
  auto F = nlohmann::json::parse(in);
  for (const auto& f : F) {
    const std::string out = std::string(outDir) + "/" + f["file"].get<std::string>();
    if (f["kind"] == "map")
      drawMap(f, out);
    else if (f["kind"] == "grid")
      drawGrid(f, out);
    else
      drawSeries(f, out);
  }
}
