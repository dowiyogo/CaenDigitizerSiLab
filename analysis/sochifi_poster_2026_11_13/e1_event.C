// e1_event.C — SOCHIFI 2026 (póster), revisión 2 (2026-11-13): figura E1, un muón real en hod2018 (event display), lienzo 2.2:1 de 3300 × 1500 px
// (37 cm de ancho impreso: 5 mm = 0.030 de la altura; textos ≥ 0.05 de la altura del pad).
// Entrada: e1_selection.json de e1_select.py (entrada elegida, G por canal, orden de tiras de mayo). Lee SOLO esa entrada del crudo v8 en modo READ
// (árbol data_dig, Ch0..Ch31 en mV) y comprueba que su timestamp es el de hits.bin.
//  - izquierda: vista desde arriba del área activa de 33 × 33 mm, 16 tiras X (verticales, miden x) y 16 tiras Y (horizontales) de 3 mm en
//    sus posiciones (centro de la tira i = −15 + 2i mm, convención de bench_mc: dos subcapas por eje desplazadas 2 mm, solapamiento de 1 mm),
//    semitransparentes; las tiras encendidas con la escala de color en fotoelectrones = aS/G (la amplitud de hitdef_extract.C); celda del
//    cruce (bordes de impact_map.py) recuadrada. El sentido absoluto del orden de tiras no está fijado (como en P6).
//  - derecha: formas de onda de los canales encendidos, (pedestal − muestra)/G en fotoelectrones, pedestal = media recortada (±2 mV de la
//    mediana) de las 40 muestras del tramo de hitdef_extract.C (registro 209-248); t = 0 en el inicio de la ventana signal (muestra 258 del
//    registro, 4128 ns); ventana signal [0, 144) ns sombreada.
// Uso (desde la raíz): root -l -b -q analysis/sochifi_poster_2026_11_13/poster_style.C \
//   -e 'gInterpreter->AddIncludePath("daq/core/include")' 'analysis/sochifi_poster_2026_11_13/e1_event.C("<dir>/e1_selection.json", "<dir>")'
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/pulse_finder.h"

// poster::color, poster::grey, poster_save: de poster_style.C (cargado antes, como en poster_figs.C)

void e1_event(const char* selJson, const char* outDir) {
  std::ifstream in(selJson);
  const auto sel = nlohmann::json::parse(in);
  const auto& ch = sel["chosen"];
  const Long64_t entry = ch["entry"].get<Long64_t>();
  const Long64_t tsWant = ch["timestamp_us"].get<Long64_t>();
  std::vector<double> G = sel["G_per_channel_mV"].get<std::vector<double>>();
  const std::vector<int> XO = sel["strip_order_may"]["X"].get<std::vector<int>>(), YO = sel["strip_order_may"]["Y"].get<std::vector<int>>();

  // --- una entrada del crudo, en modo READ
  TFile f(sel["raw_file"].get<std::string>().c_str(), "READ");
  auto* t = static_cast<TTree*>(f.Get("data_dig"));
  if (!t) throw std::runtime_error("sin data_dig");
  t->SetBranchStatus("*", 0);
  std::vector<float>* w[32] = {};
  Long_t ts = 0;
  for (int c = 0; c < 32; ++c) {
    t->SetBranchStatus(Form("Ch%d", c), 1);
    t->SetBranchAddress(Form("Ch%d", c), &w[c]);
  }
  t->SetBranchStatus("timestamp", 1);
  t->SetBranchAddress("timestamp", &ts);
  t->GetEntry(entry);
  printf("E1: entrada %lld, timestamp %ld us (hits.bin: %lld) -> %s\n", entry, ts, tsWant, ts == tsWant ? "IGUAL" : "DISTINTO");
  if (ts != tsWant) throw std::runtime_error("timestamp distinto");

  constexpr int kOff = 209, kS0 = 258;  // primera muestra del tramo de hitdef_extract.C; inicio de la ventana signal
  struct Bar {
    int channel, pos;
    std::string axis, name;
    double pe_sel;
  };
  std::vector<Bar> bars;
  for (const auto& b : ch["bars"])
    bars.push_back({b["channel"].get<int>(), b["position"].get<int>(), b["axis"].get<std::string>(), b["bar"].get<std::string>(), b["pe"].get<double>()});

  // formas de onda en p.e.
  std::vector<TGraph*> gw;
  double wmax = 0;
  for (const auto& b : bars) {
    const auto& x = *w[b.channel];
    const double bl = daq::pulse::trimmedMean(x.data() + kOff, daq::pulse::Region{0, 40}, 2.0);
    auto* g = new TGraph();
    double amin = 0;
    for (int s = kS0 - 14; s <= kS0 + 26; ++s) {
      const double a = (bl - x[static_cast<size_t>(s)]) / G[static_cast<size_t>(b.channel)];
      g->SetPoint(g->GetN(), 16.0 * (s - kS0), a);
      wmax = std::max(wmax, a);
      if (s >= kS0 && s < kS0 + 9) amin = std::max(amin, a);
    }
    printf("E1: %s ch%d pedestal %.2f mV, máx en la ventana signal %.2f p.e. (aS/G de la selección %.2f p.e.)\n", b.name.c_str(), b.channel, bl, amin,
           b.pe_sel);
    gw.push_back(g);
  }

  auto* c = new TCanvas("c_e1", "", 3300, 1500);
  c->SetWindowSize(3300 + (3300 - c->GetWw()), 1500 + (1500 - c->GetWh()));
  const double xs = 0.47;
  auto* pl = new TPad("e1l", "", 0.0, 0.0, xs, 1.0);
  auto* pr = new TPad("e1r", "", xs, 0.0, 1.0, 1.0);
  pl->SetLeftMargin(0.155); pl->SetRightMargin(0.215); pl->SetTopMargin(0.05); pl->SetBottomMargin(0.16);
  pr->SetLeftMargin(0.15); pr->SetRightMargin(0.05); pr->SetTopMargin(0.05); pr->SetBottomMargin(0.16);
  pl->Draw();
  pr->Draw();

  // --- izquierda: vista desde arriba
  pl->cd();
  const double zmax = 30.0;
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(255);
  auto* fr = new TH2D("fr_e1", ";x (mm);y (mm);fotoelectrones", 36, -18, 18, 36, -18, 18);
  fr->SetBinContent(1, 1, 0.0);
  fr->SetEntries(1);
  fr->SetMinimum(0);
  fr->SetMaximum(zmax);
  for (TAxis* a : {fr->GetXaxis(), fr->GetYaxis(), fr->GetZaxis()}) {
    a->SetLabelSize(0.055);
    a->SetTitleSize(0.060);
  }
  fr->GetXaxis()->SetTitleOffset(1.05);
  fr->GetYaxis()->SetTitleOffset(1.05);
  fr->GetZaxis()->SetTitleOffset(1.15);
  fr->GetXaxis()->SetNdivisions(506);
  fr->GetYaxis()->SetNdivisions(506);
  fr->Draw("COLZ");
  gPad->Update();
  if (auto* pal = static_cast<TPaletteAxis*>(fr->GetListOfFunctions()->FindObject("palette"))) {
    pal->SetX1NDC(0.80); pal->SetX2NDC(0.84); pal->SetY1NDC(0.16); pal->SetY2NDC(0.95);
  }
  auto palColor = [&](double pe) {
    const int n = gStyle->GetNumberOfColors();
    const int k = std::clamp(static_cast<int>(pe / zmax * (n - 1) + 0.5), 0, n - 1);
    return gStyle->GetColorPalette(k);
  };
  auto firedPe = [&](const std::string& ax, int pos) {
    for (const auto& b : bars)
      if (b.axis == ax && b.pos == pos) return b.pe_sel;
    return -1.0;
  };
  // tiras apagadas primero, luego las encendidas (encima); X verticales, Y horizontales
  for (int pass = 0; pass < 2; ++pass)
    for (const char* ax : {"X", "Y"})
      for (int i = 0; i < 16; ++i) {
        const double pe = firedPe(ax, i);
        if ((pass == 0) != (pe < 0)) continue;
        const double c0 = -15 + 2 * i;
        auto* b = std::string(ax) == "X" ? new TBox(c0 - 1.5, -16.5, c0 + 1.5, 16.5) : new TBox(-16.5, c0 - 1.5, 16.5, c0 + 1.5);
        if (pe < 0) {
          b->SetFillColorAlpha(poster::grey(), 0.12);
          b->SetLineColorAlpha(poster::grey(), 0.6);
        } else {
          b->SetFillColorAlpha(palColor(pe), 0.80);
          b->SetLineColor(kBlack);
        }
        b->SetLineWidth(2);
        b->Draw("L");
        if (std::string(ax) == "X" && XO[static_cast<size_t>(i)] == 3) {  // ch3 muerto
          auto* d = new TBox(c0 - 1.5, -16.5, c0 + 1.5, 16.5);
          d->SetFillStyle(3354); d->SetFillColor(poster::grey()); d->SetLineWidth(0);
          d->Draw();
        }
      }
  // nombres de las tiras encendidas, a lo largo de la tira, del lado opuesto al cruce (centro del cruce: medias de los centros por eje)
  double mx = 0, my = 0;
  int nx = 0, ny = 0;
  for (const auto& b : bars) (b.axis == "X" ? mx : my) += (-15 + 2 * b.pos), (b.axis == "X" ? nx : ny)++;
  mx /= std::max(nx, 1);
  my /= std::max(ny, 1);
  nx = ny = 0;
  for (const auto& b : bars) {
    const double c0 = -15 + 2 * b.pos;
    TLatex* l;
    if (b.axis == "X") {
      l = new TLatex(c0, (my > 0 ? -1 : 1) * (12.5 - 6.0 * nx++), b.name.c_str());
      l->SetTextAngle(90);
    } else {
      l = new TLatex((mx > 0 ? -1 : 1) * (12.5 - 6.0 * ny++), c0, b.name.c_str());
    }
    l->SetTextAlign(22); l->SetTextSize(0.050); l->SetTextColor(kBlack);
    l->Draw();
  }
  // celda del cruce (bordes de impact_map.py: −16.5, −14.5, −13.5, …, 14.5, 16.5)
  auto edge = [](int k) { return k == 0 ? -16.5 : (k == 31 ? 16.5 : -15.5 + k); };
  const int cx = ch["cell_x"].get<int>(), cy = ch["cell_y"].get<int>();
  auto* cell = new TBox(edge(cx), edge(cy), edge(cx + 1), edge(cy + 1));
  cell->SetFillStyle(0); cell->SetLineColor(poster::color(6)); cell->SetLineWidth(7);
  cell->Draw("L");
  printf("E1: celda x %d [%.1f, %.1f] mm, celda y %d [%.1f, %.1f] mm\n", cx, edge(cx), edge(cx + 1), cy, edge(cy), edge(cy + 1));
  gPad->RedrawAxis();

  // --- derecha: formas de onda
  pr->cd();
  const double ytop = 1.30 * wmax;
  auto* fw = new TH1F("fw_e1", ";tiempo (ns);amplitud (fotoelectrones)", 10, -200, 400);
  fw->SetMinimum(-0.08 * wmax);
  fw->SetMaximum(ytop);
  for (TAxis* a : {fw->GetXaxis(), fw->GetYaxis()}) {
    a->SetLabelSize(0.055);
    a->SetTitleSize(0.060);
  }
  fw->GetXaxis()->SetTitleOffset(1.05);
  fw->GetYaxis()->SetTitleOffset(1.05);
  fw->Draw("AXIS");
  auto* win = new TBox(0, fw->GetMinimum(), 144, ytop);
  win->SetFillColorAlpha(poster::grey(), 0.30);
  win->SetLineWidth(0);
  win->Draw();
  const int cols[4] = {poster::color(1), poster::color(6), poster::color(2), poster::color(5)};
  auto* lg = new TLegend(0.55, 0.50, 0.94, 0.93);
  lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.055);
  for (size_t i = 0; i < gw.size(); ++i) {
    gw[i]->SetLineColor(cols[i % 4]); gw[i]->SetLineWidth(5);
    gw[i]->SetMarkerColor(cols[i % 4]); gw[i]->SetMarkerStyle(20); gw[i]->SetMarkerSize(1.6);
    gw[i]->Draw("LP SAME");
    lg->AddEntry(gw[i], bars[i].name.c_str(), "lp");
  }
  lg->AddEntry(win, "ventana de se\361al", "f");
  lg->Draw();
  gPad->RedrawAxis();
  c->cd();
  poster_save(c, std::string(outDir) + "/E1_muon_real");
}
