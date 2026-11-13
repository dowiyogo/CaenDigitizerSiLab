// poster_style.C — SOCHIFI 2026 (póster A0), 2026-11-13: estilo único de las figuras del póster, derivado de daq/tools/proceeding/proc_style.C.
//  - lienzo 4:3 de 2400 × 1800 px (20.3 × 15.2 cm a 300 ppp); PDF vectorial y PNG;
//  - letra Helvetica (fuente 42) con tamaño relativo a la altura del pad: rótulos y títulos ≥ 0.05, leyenda 0.05;
//  - líneas de ancho ≥ 3 (histogramas, gráficos, marco, funciones);
//  - sin títulos dentro de la figura ni cajas de estadística; paleta Okabe-Ito de proc_style.C; gris para fondos y bandas.
// Uso: gROOT->ProcessLine(".L analysis/sochifi_poster_2026_11_13/poster_style.C"); poster_style(); auto c = poster_canvas("c"); …;
//      poster_save(c, "<ruta sin extensión>");
#include <TCanvas.h>
#include <TColor.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>

namespace poster {
inline const char* kHex[8] = {"#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"};
inline int color(int i) { return TColor::GetColor(kHex[i % 8]); }
inline int grey() { return TColor::GetColor("#999999"); }
constexpr int kW = 2400, kH = 1800;
constexpr double kLabel = 0.052, kTitle = 0.058, kLegend = 0.050, kText = 0.050;
constexpr int kLine = 3;
}  // namespace poster

void poster_style() {
  auto* s = new TStyle("poster", "SOCHIFI 2026");
  s->SetCanvasColor(0);
  s->SetPadColor(0);
  s->SetFrameFillColor(0);
  s->SetCanvasBorderMode(0);
  s->SetPadBorderMode(0);
  s->SetFrameBorderMode(0);
  s->SetOptTitle(0);
  s->SetOptStat(0);
  s->SetOptFit(0);
  for (const char* ax : {"x", "y", "z"}) {
    s->SetLabelFont(42, ax);
    s->SetTitleFont(42, ax);
    s->SetLabelSize(poster::kLabel, ax);
    s->SetTitleSize(poster::kTitle, ax);
  }
  s->SetTextFont(42);
  s->SetTextSize(poster::kText);
  s->SetLegendFont(42);
  s->SetLegendTextSize(poster::kLegend);
  s->SetLegendBorderSize(0);
  s->SetLegendFillColor(0);
  s->SetTitleOffset(1.15, "x");
  s->SetTitleOffset(1.30, "y");
  s->SetPadLeftMargin(0.16);
  s->SetPadRightMargin(0.04);
  s->SetPadBottomMargin(0.16);
  s->SetPadTopMargin(0.05);
  s->SetPadTickX(1);
  s->SetPadTickY(1);
  s->SetEndErrorSize(0);
  s->SetMarkerSize(2.0);
  s->SetLineWidth(poster::kLine);
  s->SetHistLineWidth(poster::kLine);
  s->SetFuncWidth(poster::kLine + 1);
  s->SetFrameLineWidth(poster::kLine);
  s->SetGridWidth(1);
  s->SetNdivisions(506, "x");
  s->SetNdivisions(506, "y");
  gROOT->SetStyle("poster");
  gROOT->ForceStyle();
}

TCanvas* poster_canvas(const char* name) {
  auto* c = new TCanvas(name, "", poster::kW, poster::kH);
  c->SetWindowSize(poster::kW + (poster::kW - c->GetWw()), poster::kH + (poster::kH - c->GetWh()));
  return c;
}

void poster_save(TCanvas* c, const std::string& base) {
  c->SaveAs((base + ".pdf").c_str());
  c->SaveAs((base + ".png").c_str());
}
