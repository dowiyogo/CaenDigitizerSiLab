// proc_style.C — PROC-1: estilo único de las figuras del proceeding (A2).
//  - lienzos para una columna (8.6 cm) y dos columnas (17.8 cm), alto proporcional; PNG a 300 ppp;
//  - letra Helvetica con tamaño absoluto (precisión 3) equivalente a 8 pt a la escala de impresión;
//  - sin títulos dentro de la figura; los ejes llevan unidades (las pone cada macro);
//  - paleta apta para daltonismo (Okabe-Ito) y gris para las series fuera de tiempo o mezcladas.
// Uso: gROOT->ProcessLine(".L daq/tools/proceeding/proc_style.C"); proc_style(); auto c = proc_canvas("c", 1, 0.75); …;
//      proc_save(c, "<ruta sin extensión>");
#include <TCanvas.h>
#include <TColor.h>
#include <TImage.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>

namespace proc {
// Okabe-Ito: negro, naranja, celeste, verde azulado, amarillo, azul, bermellón, púrpura rojizo
inline const char* kHex[8] = {"#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"};
inline int color(int i) { return TColor::GetColor(kHex[i % 8]); }
inline int grey() { return TColor::GetColor("#999999"); }
// 8 pt a la escala de impresión: el lienzo de una columna tiene 8.6 cm = 1016 px a 300 ppp; 8 pt = 2.822 mm = 33 px.
constexpr double kPxPerCm = 300.0 / 2.54;
constexpr int kFontPx = 33;
}  // namespace proc

void proc_style() {
  auto* s = new TStyle("proc", "PROC-1");
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
    s->SetLabelFont(43, ax);
    s->SetTitleFont(43, ax);
    s->SetLabelSize(proc::kFontPx, ax);
    s->SetTitleSize(proc::kFontPx, ax);
  }
  s->SetTextFont(43);
  s->SetTextSize(proc::kFontPx);
  s->SetLegendFont(43);
  s->SetLegendTextSize(proc::kFontPx);
  s->SetLegendBorderSize(0);
  s->SetLegendFillColor(0);
  s->SetTitleOffset(1.25, "x");
  s->SetTitleOffset(1.55, "y");
  s->SetPadLeftMargin(0.15);
  s->SetPadRightMargin(0.04);
  s->SetPadBottomMargin(0.15);
  s->SetPadTopMargin(0.04);
  s->SetPadTickX(1);
  s->SetPadTickY(1);
  s->SetEndErrorSize(0);
  s->SetMarkerSize(1.2);
  s->SetLineWidth(2);
  s->SetHistLineWidth(2);
  s->SetFrameLineWidth(2);
  s->SetPalette(kCividis);  // mapas: paleta secuencial apta para daltonismo
  s->SetNumberContours(64);
  gROOT->SetStyle("proc");
  gROOT->ForceStyle();
}

// columns = 1 (8.6 cm) o 2 (17.8 cm); aspect = alto / ancho
TCanvas* proc_canvas(const char* name, int columns, double aspect) {
  const double wcm = columns == 2 ? 17.8 : 8.6;
  const int w = static_cast<int>(wcm * proc::kPxPerCm + 0.5), h = static_cast<int>(wcm * aspect * proc::kPxPerCm + 0.5);
  auto* c = new TCanvas(name, "", w, h);
  c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
  return c;
}

void proc_save(TCanvas* c, const std::string& base) {
  c->SaveAs((base + ".pdf").c_str());
  c->SaveAs((base + ".png").c_str());
}
