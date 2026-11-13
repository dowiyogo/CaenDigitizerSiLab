// poster_figs.C — SOCHIFI 2026 (póster A0), 2026-11-13: figuras P1-P8 con el estilo de poster_style.C, a partir de salidas existentes de run B
// (sin selecciones, ajustes ni estimadores nuevos; ver el informe y README de outputs/sochifi_poster/2026-11-13_figures/).
//  P1 intervalos entre triggers, t ≥ 5.25 h: mismo histograma que analysis/followup11_2026_11_02/interarrival_stationary.C (100 bines de ΔT/⟨ΔT⟩
//     en [0, 10), densidad normalizada) releyendo los tiempos del crudo v8 en modo READ; λ y p KS copiados del log de esa macro (sin reajustar).
//  P2 espectro etiquetado de XB, 29-54 mm (spectra_meas.csv): misma curva Landau ⊗ Gauss que mip_paper_2p.C (mismo código, semilla y límites;
//     se comprueba que la moda, n y χ² salen iguales al log); moda marcada.
//  P3 ε_XY contra el umbral (fig_eff_g09551_bandrunB.json): puntos, barras estadísticas, cajas del estimador, banda de aceptancia.
//  P4 fracción de dos barras contra el umbral (doublehit_confirmed.json): X e Y medidos (todos los pares) y MC con detección perfecta.
//  P5 σ(t_X − t_Y) por bin de la menor amplitud (time_res.csv: CFD 0.5, con desplazamientos y propagación).
//  P6 perfiles de impacto x e y a 10 p.e., stationary_goodtime (impact_map_runB_T10.csv): canal muerto sombreado, uniones del cableado marcadas.
//  P7 G/G_ref por bin de 1 h contra la hora local de Chile (diurnal_bins.csv, bins buenos), mayo y octubre; pico a pico de poster_numbers.json.
//  P8 ocupación por tira de hod2018 (5 p.e.) y hod2019 (3 mV), dos paneles 2.2:1 (fig_stage3.json, la entrada de occupancy_2018_2019).
// Revisión 1 (2026-11-13): tamaños de texto por tamaño impreso (TextSizes), «disparo(s)» en lugar de «trigger(s)», leyendas de P1, P4, P6.
// Revisión 2 (2026-11-13): P3poster y P6xposter (claves 'a', 'b'): rótulos para público general, mismos datos; E1 en e1_event.C.
// Uso (desde la raíz): root -l -b -q analysis/sochifi_poster_2026_11_13/poster_style.C 'analysis/sochifi_poster_2026_11_13/poster_figs.C("<out>")'
#include <TBox.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TMath.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace {
const std::string R = "outputs/";
std::vector<std::map<std::string, std::string>> readCsv(const std::string& p) {
  std::ifstream in(p);
  std::string line;
  std::vector<std::string> hdr;
  std::vector<std::map<std::string, std::string>> rows;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("\"#", 0) == 0) continue;
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
nlohmann::json readJson(const std::string& p) {
  std::ifstream in(p);
  return nlohmann::json::parse(in);
}
TLegend* legend(double x1, double y1, double x2, double y2) {
  auto* l = new TLegend(x1, y1, x2, y2);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(poster::kLegend);
  return l;
}
// copia de langaufun de analysis/followup10_2026_11_01/mip_fit.C (langaus del tutorial de ROOT)
double langaufun(double* x, double* par) {
  const double invsq2pi = 0.3989422804014, mpshift = -0.22278298;
  const double np = 100.0, sc = 5.0;
  const double mpc = par[1] - mpshift * par[0];
  const double xlow = x[0] - sc * par[3], xupp = x[0] + sc * par[3], step = (xupp - xlow) / np;
  double sum = 0.0;
  for (double i = 1.0; i <= np / 2; i++) {
    double xx = xlow + (i - .5) * step;
    double fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }
  return par[2] * step * sum * invsq2pi / par[3];
}
// tamaños de texto por figura, en fracción de la altura del pad (revisión 1, 2026-11-13): ≥ 5 mm impresos en A0.
//  P1 columna completa (≈ 18 × 13.5 cm): 5 mm = 0.037 de la altura; P4, P5, P6 media columna (≈ 12.5 × 9.4 cm): 5 mm = 0.053;
//  P8 columna completa 2.2:1 (≈ 25 × 11.4 cm): 5 mm = 0.044.
struct TextSizes {
  double label, title, legend, xoff, yoff;
};
const TextSizes kFull{0.052, 0.058, 0.050, 1.15, 1.30};
const TextSizes kHalf{0.060, 0.065, 0.058, 1.12, 1.30};
const TextSizes kWide{0.055, 0.060, 0.052, 1.05, 1.05};
void applySizes(TH1* h, const TextSizes& t) {
  for (TAxis* a : {h->GetXaxis(), h->GetYaxis()}) {
    a->SetLabelSize(t.label);
    a->SetTitleSize(t.title);
  }
  h->GetXaxis()->SetTitleOffset(t.xoff);
  h->GetYaxis()->SetTitleOffset(t.yoff);
}
void halfMargins() {
  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.19); gPad->SetTopMargin(0.05);
}
TLegend* legendT(double x1, double y1, double x2, double y2, double ts) {
  auto* l = legend(x1, y1, x2, y2);
  l->SetTextSize(ts);
  return l;
}
}  // namespace

// ------------------------------------------------------------------------------------------------ P1
void P1(const std::string& out) {
  const char* fname = "/data/daq/raw/hodoscope_extnim_2026/lectura_vx1740_ext_nim_1M_v8.root";
  TFile* f = TFile::Open(fname, "READ");
  auto* t = static_cast<TTree*>(f->Get("data_dig"));
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("timestamp", 1);
  Long_t ts_val = 0;
  t->SetBranchAddress("timestamp", &ts_val);
  const long long N0 = t->GetEntries();
  std::vector<double> ts(N0);
  for (long long i = 0; i < N0; i++) { t->GetEntry(i); ts[i] = static_cast<double>(ts_val); }
  f->Close();
  std::sort(ts.begin(), ts.end());
  const double t0 = ts[0];
  std::vector<double> keep;
  for (double v : ts) if (v - t0 >= 5.25 * 3.6e9) keep.push_back(v - t0);
  std::vector<double> dt;
  for (size_t i = 0; i + 1 < keep.size(); i++) dt.push_back(keep[i + 1] - keep[i]);
  double s = 0;
  for (double v : dt) s += v;
  const double dt_mean = s / dt.size();
  auto* h = new TH1D("h_dt", ";intervalo entre disparos / media;densidad", 100, 0, 10);
  long long nout = 0;
  for (double v : dt) { const double x = v / dt_mean; if (x < 10) h->Fill(x); else ++nout; }
  printf("P1: triggers %zu, intervalos %zu, <dT> = %.1f us, en [0,10) %.0f, sobre 10 %lld\n", keep.size(), dt.size(), dt_mean, h->Integral(), nout);
  h->Scale(1.0 / h->Integral() / 0.1);
  auto* c = poster_canvas("c_p1");
  c->SetLogy();
  applySizes(h, kFull);
  h->SetLineColor(poster::color(5));
  h->SetLineWidth(4);
  h->SetMinimum(0.5 * h->GetMinimum(0));
  h->SetMaximum(3.0);
  h->Draw("HIST");
  auto* fth = new TF1("fth", "exp(-x)", 0, 10);
  fth->SetLineColor(poster::color(6));
  fth->SetLineWidth(5);
  fth->SetLineStyle(9);
  fth->Draw("SAME");
  // abajo a la izquierda: para x < 4.6 la densidad es > 0.01, la curva pasa por encima de la leyenda
  auto* lg = legendT(0.19, 0.19, 0.66, 0.47, kFull.legend);
  lg->AddEntry(h, "run B, desde 5.25 h", "l");
  lg->AddEntry(fth, "exponencial", "l");
  lg->AddEntry((TObject*)nullptr, "#lambda = 0.5947#pm0.0017 Hz", "");
  lg->AddEntry((TObject*)nullptr, "Kolmogorov-Smirnov: p = 0.27", "");
  lg->Draw();
  poster_save(c, out + "/P1_intervalos_trigger");
}

// ------------------------------------------------------------------------------------------------ P2
void P2(const std::string& out) {
  auto* h = new TH1D("hp_XB_far", ";amplitud (p.e.);muones / p.e.", 300, 0, 300);
  for (auto& r : readCsv(R + "proceeding_revision/2026-11-02_revision/spectra_meas.csv")) {
    if (r["layer"] != "XB" || r["pos_bin"] != "far") continue;
    const int b = h->FindBin(std::stod(r["pe_lo"]) + 0.5);
    h->SetBinContent(b, std::stod(r["n_net"]));
    h->SetBinError(b, std::stod(r["err"]));
  }
  const double lo = 4.0, hi = 105.0;
  const double integ = h->Integral(h->FindBin(lo + 0.5), h->FindBin(hi - 0.5));
  int bmax = h->FindBin(lo + 0.5);
  double best = -1e9;
  for (int b = h->FindBin(lo + 0.5) + 1; b < h->FindBin(hi - 0.5); b++) {
    const double m3 = h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1);
    if (m3 > best) best = m3, bmax = b;
  }
  const double mp0 = h->GetBinCenter(bmax);
  auto* f = new TF1("fp_XB_far", langaufun, lo, hi, 4);
  f->SetParameters(0.15 * mp0 + 0.3, mp0, integ, 0.3 * mp0 + 0.5);
  f->SetParLimits(0, 0.01, 50); f->SetParLimits(1, 0.1, 105); f->SetParLimits(2, 0, 10 * std::max(integ, 1.0)); f->SetParLimits(3, 0.01, 50);
  f->SetNpx(1500);
  const int st = h->Fit(f, "RBQ0");
  const double md = f->GetMaximumX(lo, hi);
  printf("P2: XB far status %d n %.0f moda %.2f (MPV %.2f) chi2 %.1f/%d (log de mip_paper_2p.C: status 0 n 865 moda 13.23 (MPV 12.81) chi2 60.2/80)\n", st, integ, md,
         f->GetParameter(1), f->GetChisquare(), f->GetNDF());
  auto* c = poster_canvas("c_p2");
  h->SetLineColor(poster::color(0)); h->SetMarkerColor(poster::color(0)); h->SetMarkerStyle(20); h->SetMarkerSize(1.6);
  h->GetXaxis()->SetRangeUser(0, 60);
  double ymax = 0;
  for (int b = h->FindBin(0.5); b <= h->FindBin(59.5); b++) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
  const double ytop = 1.55 * ymax;
  h->SetMinimum(0); h->SetMaximum(ytop);
  h->Draw("E1");
  f->SetLineColor(poster::color(5)); f->SetLineWidth(5); f->SetRange(4.0, 60);
  f->Draw("SAME");
  auto* lm = new TGraph(2);
  lm->SetPoint(0, md, 0); lm->SetPoint(1, md, f->Eval(md));
  lm->SetLineColor(poster::color(6)); lm->SetLineWidth(5); lm->SetLineStyle(9);
  lm->Draw("L SAME");
  h->Draw("E1 SAME");
  gPad->RedrawAxis();
  auto* lg = legend(0.48, 0.62, 0.95, 0.93);
  lg->AddEntry(h, "capa XB, 29-54 mm", "pe");
  lg->AddEntry(f, "Landau #otimes Gauss", "l");
  lg->AddEntry(lm, Form("moda = %.1f p.e.", md), "l");
  lg->AddEntry((TObject*)nullptr, Form("n_{neto} = %.0f muones", integ), "");
  lg->Draw();
  poster_save(c, out + "/P2_espectro_XB");
}

// ------------------------------------------------------------------------------------------------ P3
void P3(const std::string& out) {
  auto j = readJson(R + "proceeding_revision/2026-11-03_band_legend/fig_eff_g09551_bandrunB.json");
  const auto& s = j["series"][0];
  const double rel = j["f18_rel_syst"];
  const int n = static_cast<int>(s["eps"].size());
  auto* c = poster_canvas("c_p3");
  auto* fr = new TH1F("fr3", ";umbral T (p.e.);#varepsilon_{XY}", 10, 2.5, 12.5);
  fr->SetMinimum(0.0);
  fr->SetMaximum(1.4);
  fr->Draw("AXIS");
  auto* band = new TGraphErrors(n), *g = new TGraphErrors(n), *gs = new TGraphErrors(n);
  for (int i = 0; i < n; ++i) {
    const double T = j["T"][i], e = s["eps"][i];
    band->SetPoint(i, T, e); band->SetPointError(i, 0.5, e * rel);
    g->SetPoint(i, T, e); g->SetPointError(i, 0, s["stat"][i].get<double>());
    gs->SetPoint(i, T, e); gs->SetPointError(i, 0.12, s["syst"][i].get<double>());
    printf("P3: T %d eps %.4f stat %.4f estimador %.4f aceptancia %.4f\n", j["T"][i].get<int>(), e, s["stat"][i].get<double>(), s["syst"][i].get<double>(), e * rel);
  }
  band->SetFillColorAlpha(poster::grey(), 0.40); band->SetLineColor(0);
  band->Draw("E3 SAME");
  gs->SetFillStyle(0); gs->SetLineColor(poster::color(0)); gs->SetLineWidth(3);
  gs->Draw("2 SAME");
  g->SetMarkerStyle(20); g->SetMarkerColor(poster::color(0)); g->SetLineColor(poster::color(0)); g->SetLineWidth(3);
  g->Draw("P SAME");
  gPad->RedrawAxis();
  auto* lg = legend(0.38, 0.68, 0.95, 0.93);
  lg->AddEntry(g, "run B, error estad\355stico", "pe");
  lg->AddEntry(gs, "estimador de accidentales", "f");
  lg->AddEntry(band, Form("aceptancia (%.2f %%)", 100 * rel), "f");
  lg->Draw();
  poster_save(c, out + "/P3_eficiencia_umbral");
}

// ------------------------------------------------------------------------------------------------ P4
void P4(const std::string& out) {
  auto j = readJson(R + "proceeding_hvdiag5/2026-10-31_full/work/doublehit_confirmed.json");
  auto* c = poster_canvas("c_p4");
  halfMargins();
  auto* fr = new TH1F("fr4", ";umbral T (p.e.);fracci\363n de dos barras f_{2}", 10, 2.5, 12.5);
  applySizes(fr, kHalf);
  fr->SetMinimum(0.0);
  fr->SetMaximum(0.82);
  fr->Draw("AXIS");
  auto* lg = legendT(0.21, 0.72, 0.95, 0.93, kHalf.legend);
  lg->SetNColumns(2);
  int k = 0;
  for (const char* pl : {"X", "Y"}) {
    const auto& rows = j["data"]["may8_good"][pl]["rows"];
    auto* g = new TGraphErrors();
    double mc = 0;
    for (const auto& r : rows) {
      const int i = g->GetN();
      g->SetPoint(i, r["T"].get<double>() + (k == 0 ? -0.1 : 0.1), r["all"].get<double>());
      g->SetPointError(i, 0, r["all_err"].get<double>());
      mc = r["mc_all_cut0.0"].get<double>();
      if (r["T"] == 5 || r["T"] == 10) printf("P4: %s T %d f2 %.4f +- %.4f MC %.4f\n", pl, r["T"].get<int>(), r["all"].get<double>(), r["all_err"].get<double>(), mc);
    }
    const int col = poster::color(k == 0 ? 5 : 6);
    g->SetMarkerStyle(k == 0 ? 20 : 21); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(3);
    auto* l = new TGraph(2);
    l->SetPoint(0, 2.5, mc); l->SetPoint(1, 12.5, mc);
    l->SetLineColor(col); l->SetLineStyle(9); l->SetLineWidth(5);
    l->Draw("L SAME");
    g->Draw("P SAME");
    lg->AddEntry(g, Form("%s medido", pl), "pe");
    lg->AddEntry(l, Form("%s MC ideal", pl), "l");
    ++k;
  }
  lg->Draw();
  gPad->RedrawAxis();
  poster_save(c, out + "/P4_fraccion_dos_barras");
}

// ------------------------------------------------------------------------------------------------ P5
void P5(const std::string& out) {
  auto rows = readCsv(R + "hod2018_timing/2026-11-03_followup13/A_time/time_res.csv");
  const std::vector<std::string> bins{"10-15", "15-20", "20-30", ">30"};
  auto* c = poster_canvas("c_p5");
  halfMargins();
  auto* fr = new TH1F("fr5", ";menor amplitud de X e Y (p.e.);#sigma(t_{X}#minust_{Y}) (ns)", 4, 0.5, 4.5);
  applySizes(fr, kHalf);
  for (int i = 0; i < 4; ++i) fr->GetXaxis()->SetBinLabel(i + 1, i == 3 ? "> 30" : bins[i].c_str());
  fr->GetXaxis()->SetLabelSize(0.075);  // con «más de 30» ROOT achicaba las etiquetas para que cupieran en el bin
  fr->SetMinimum(0);
  fr->SetMaximum(3.5);
  fr->Draw("AXIS");
  auto* g = new TGraphErrors();
  for (size_t i = 0; i < bins.size(); ++i)
    for (auto& r : rows)
      if (r["method"] == "cfd50" && r["correction"] == "offsets+prop" && r["min_amp_bin_pe"] == bins[i]) {
        g->SetPoint(i, i + 1, std::stod(r["sigma_ns"]));
        g->SetPointError(i, 0, std::stod(r["sigma_err_ns"]));
        printf("P5: %s sigma %s +- %s ns (n %s)\n", bins[i].c_str(), r["sigma_ns"].c_str(), r["sigma_err_ns"].c_str(), r["n_same"].c_str());
      }
  g->SetMarkerStyle(20); g->SetMarkerColor(poster::color(5)); g->SetLineColor(poster::color(5)); g->SetLineWidth(3);
  g->Draw("P SAME");
  auto* tx = new TLatex(0.93, 0.90, "fracci\363n constante 0.5");
  tx->SetNDC(); tx->SetTextAlign(33); tx->SetTextSize(kHalf.legend);
  tx->Draw();
  auto* tx2 = new TLatex(0.93, 0.82, "muestreo de 16 ns");
  tx2->SetNDC(); tx2->SetTextAlign(33); tx2->SetTextSize(kHalf.legend);
  tx2->Draw();
  poster_save(c, out + "/P5_resolucion_temporal");
}

// ------------------------------------------------------------------------------------------------ P6
void P6(const std::string& out) {
  auto rows = readCsv(R + "hod2018_impact_map/2026-11-12_impact_map_runB/stationary_goodtime/impact_map_runB_T10.csv");
  std::vector<double> edges{-16.5};
  for (int k = 0; k < 30; ++k) edges.push_back(-14.5 + k);
  edges.push_back(16.5);
  for (const char* ax : {"x", "y"}) {
    auto* h = new TH1D(Form("h6%s", ax), Form(";posici\363n %s (mm);cuentas restadas por mm", ax), 31, edges.data());
    std::vector<int> junction, dead;
    long ntrig = 0;
    for (auto& r : rows) {
      if (r["axis"] != ax) continue;
      const int b = std::stoi(r["cell_index"]) + 1;
      h->SetBinContent(b, std::stod(r["counts_per_mm"]));
      h->SetBinError(b, std::stod(r["counts_per_mm_err"]));
      if (r["flags"].find("junction_from_cabling") != std::string::npos) junction.push_back(b);
      if (r["dead_channel_flag"] == "1") dead.push_back(b);
      ntrig = std::stol(r["triggers_denominator"]);
    }
    double ymax = 0;
    for (int b = 1; b <= 31; ++b) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    const double ytop = 1.75 * ymax;
    auto* c = poster_canvas(Form("c_p6%s", ax));
    halfMargins();
    applySizes(h, kHalf);
    h->SetMinimum(0); h->SetMaximum(ytop);
    h->SetLineColor(poster::color(5)); h->SetMarkerColor(poster::color(5)); h->SetMarkerStyle(20); h->SetLineWidth(3);
    h->Draw("E1");
    TBox* bdead = nullptr;
    if (!dead.empty()) {  // celdas 23-25 sin cuentas posibles (7.5-10.5 mm) y vecinas 22, 26 en gris claro
      for (int b : dead) {
        const bool core = (b != dead.front() && b != dead.back());
        auto* bx = new TBox(h->GetBinLowEdge(b), 0, h->GetBinLowEdge(b + 1), ytop);
        bx->SetFillColorAlpha(poster::grey(), core ? 0.55 : 0.25);
        bx->SetLineWidth(0);
        bx->Draw();
        if (core && !bdead) bdead = bx;
      }
      auto* t = new TLatex(0.5 * (h->GetBinLowEdge(dead.front() + 1) + h->GetBinLowEdge(dead.back())), 0.5 * ymax, "canal muerto");
      t->SetTextAngle(90); t->SetTextAlign(22); t->SetTextSize(kHalf.legend);
      t->Draw();
    }
    TBox* bj = nullptr;
    for (int b : junction) {
      bj = new TBox(h->GetBinLowEdge(b), 0, h->GetBinLowEdge(b + 1), ytop);
      bj->SetFillStyle(3354); bj->SetFillColor(poster::color(1)); bj->SetLineWidth(0);
      bj->Draw();
    }
    h->Draw("E1 SAME");
    gPad->RedrawAxis();
    auto* lg = legendT(0.21, 0.76, 0.95, 0.93, kHalf.legend);
    lg->SetNColumns(1);
    lg->AddEntry(h, "run B, 10 p.e.", "pe");
    if (bj) lg->AddEntry(bj, "uni\363n seg\372n cableado", "f");
    lg->Draw();
    printf("P6 %s: %ld disparos, celdas de union %zu, celdas del canal muerto %zu\n", ax, ntrig, junction.size(), dead.size());
    poster_save(c, out + Form("/P6_perfil_impacto_%s", ax));
  }
}

// ------------------------------------------------------------------------------------------------ P7
void P7(const std::string& out) {
  auto rows = readCsv(R + "hod_stability/2026-11-03_followup12/B_diurnal/diurnal_bins.csv");
  // diurnal_summary.json tiene NaN (no es JSON estándar): la modulación se toma de poster_numbers.json (make_numbers.py la copia de ahí)
  auto nm = readJson(out + "/poster_numbers.json")["numbers"];
  auto* c = poster_canvas("c_p7");
  auto* fr = new TH1F("fr7", ";hora local de Chile (h, UTC#minus4 mayo, UTC#minus3 octubre);G / G_{ref}", 24, 0, 24);
  fr->SetMinimum(0.93);
  fr->SetMaximum(1.09);
  fr->GetXaxis()->SetNdivisions(406);
  fr->Draw("AXIS");
  auto* lg = legend(0.20, 0.74, 0.95, 0.93);
  int k = 0;
  for (const char* camp : {"may", "oct"}) {
    auto* g = new TGraphErrors();
    for (auto& r : rows) {
      if (r["campaign"] != camp || r["good_bin"] != "1" || r["G_rel_median"] == "nan") continue;
      const int i = g->GetN();
      g->SetPoint(i, std::stod(r["local_h_mid"]), std::stod(r["G_rel_median"]));
      g->SetPointError(i, 0, std::stod(r["G_rel_median_err"]));
    }
    const auto& fit = nm[std::string("G_diurnal_peak_to_peak_") + camp];
    const int col = poster::color(k == 0 ? 5 : 6);
    g->SetMarkerStyle(k == 0 ? 20 : 25); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(3);
    g->Draw("P SAME");
    lg->AddEntry(g, Form("%s: pico a pico %.1f#pm%.1f %%", std::string(camp) == "may" ? "mayo" : "octubre",
                         fit["value"].get<double>(), fit["error"].get<double>()), "pe");
    printf("P7: %s %d bins buenos, pico a pico %.4f +- %.4f %%\n", camp, g->GetN(), fit["value"].get<double>(), fit["error"].get<double>());
    ++k;
  }
  auto* l1 = new TLine(0, 1, 24, 1);
  l1->SetLineStyle(3); l1->SetLineWidth(3);
  l1->Draw();
  lg->Draw();
  gPad->RedrawAxis();
  poster_save(c, out + "/P7_estabilidad_ganancia");
}

// ------------------------------------------------------------------------------------------------ P8
// revisión 1: dos paneles lado a lado (2.2:1), eje y común «impactos por 1000 disparos»: izquierda hod2018 a 5 p.e. (X e Y, bandas MC,
// canal muerto), derecha hod2019 a 3 mV (X e Y). Mismos datos que occupancy_2018_2019 (fig_stage3.json). Bandas MC normalizadas a la suma medida
// (como en occupancy_2018_2019).
void P8(const std::string& out) {
  auto j = readJson(R + "proceeding_hvdiag5/2026-10-31_full/work/fig_stage3.json");
  std::map<std::string, const nlohmann::json*> P;
  for (const auto& panel : j["occ"]) P[panel["panel"].get<std::string>()] = &panel;
  auto series = [&](const std::string& panel, const std::string& name) -> const nlohmann::json& {
    for (const auto& s : (*P.at(panel))["series"])
      if (s["name"] == name) return s;
    throw std::runtime_error("serie no encontrada: " + panel + " " + name);
  };
  const std::vector<std::pair<std::string, std::string>> left{{"hod2018 X", "5 p.e."}, {"hod2018 Y", "5 p.e."}},
      right{{"hod2019 X", "3 mV"}, {"hod2019 Y", "3 mV"}};
  double ymax = 0;
  for (const auto* side : {&left, &right})
    for (const auto& [pn, sn] : *side) {
      const auto& s = series(pn, sn);
      for (int i = 0; i < 16; ++i) ymax = std::max(ymax, s["y"][i].get<double>() + s["e"][i].get<double>());
    }
  for (const char* pn : {"hod2018 X", "hod2018 Y"})
    for (int i = 0; i < 16; ++i) ymax = std::max(ymax, (*P.at(pn))["mc"][i].get<double>() + (*P.at(pn))["mc_err"][i].get<double>());
  const double ytop = 1.85 * ymax;
  auto* c = new TCanvas("c_p8", "", 3300, 1500);
  c->SetWindowSize(3300 + (3300 - c->GetWw()), 1500 + (1500 - c->GetWh()));
  const double xs = 0.535;
  auto* pl = new TPad("p8l", "", 0.0, 0.0, xs, 1.0);
  auto* pr = new TPad("p8r", "", xs, 0.0, 1.0, 1.0);
  for (auto* pd : {pl, pr}) { pd->SetTopMargin(0.05); pd->SetBottomMargin(0.18); pd->Draw(); }
  pl->SetLeftMargin(0.17); pl->SetRightMargin(0.015);
  pr->SetLeftMargin(0.045); pr->SetRightMargin(0.03);
  int ip = 0;
  for (auto* pd : {pl, pr}) {
    pd->cd();
    const bool isL = pd == pl;
    auto* fr = new TH1F(Form("fr8_%d", ip++), isL ? ";coordenada de la tira;impactos por 1000 disparos" : ";coordenada de la tira;", 16, -0.5, 15.5);
    applySizes(fr, kWide);
    fr->GetYaxis()->SetTitleOffset(1.15);
    if (!isL) fr->GetYaxis()->SetLabelSize(0);  // eje y común: rótulos solo a la izquierda
    fr->SetMinimum(0);
    fr->SetMaximum(ytop);
    fr->Draw("AXIS");
    // la banda del canal muerto (X12) empieza en x ≈ 0.78 del pad izquierdo: la leyenda termina antes
    auto* lg = legendT(isL ? 0.20 : 0.08, isL ? 0.62 : 0.72, isL ? 0.76 : 0.80, 0.94, kWide.legend);
    lg->SetNColumns(2);
    lg->SetHeader(isL ? "hod2018, 5 p.e." : "hod2019, 3 mV");
    if (isL) {
      const auto& px = *P.at("hod2018 X");
      for (int k = 0; k < 16; ++k)
        if (px["dead"][k].get<bool>()) {
          auto* bx = new TBox(k - 0.5, 0, k + 0.5, ytop);
          bx->SetFillColorAlpha(poster::grey(), 0.55); bx->SetLineWidth(0);
          bx->Draw();
          auto* t = new TLatex(k, 0.30 * ytop, "canal muerto (X)");
          t->SetTextAngle(90); t->SetTextAlign(22); t->SetTextSize(kWide.legend);
          t->Draw();
        }
      int kk = 0;
      for (const char* pn : {"hod2018 X", "hod2018 Y"}) {
        const auto& pp = *P.at(pn);
        auto* band = new TGraphErrors(16);
        for (int k = 0; k < 16; ++k) { band->SetPoint(k, k, pp["mc"][k].get<double>()); band->SetPointError(k, 0.45, pp["mc_err"][k].get<double>()); }
        band->SetFillColorAlpha(poster::color(kk == 0 ? 5 : 6), 0.25); band->SetLineColor(0);
        band->Draw("2 SAME");
        lg->AddEntry(band, kk == 0 ? "MC X" : "MC Y", "f");
        ++kk;
      }
    }
    int k = 0;
    for (const auto& [pn, sn] : isL ? left : right) {
      const auto& s = series(pn, sn);
      auto* g = new TGraphErrors(16);
      for (int i = 0; i < 16; ++i) { g->SetPoint(i, i + (k == 0 ? -0.12 : 0.12), s["y"][i].get<double>()); g->SetPointError(i, 0, s["e"][i].get<double>()); }
      const int col = poster::color(k == 0 ? 5 : 6);
      g->SetMarkerStyle(k == 0 ? 20 : 21); g->SetMarkerSize(1.7); g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(3);
      g->Draw("P SAME");
      lg->AddEntry(g, k == 0 ? "plano X" : "plano Y", "pe");
      printf("P8: %s %s, suma de la serie %.3f\n", pn.c_str(), sn.c_str(), [&] { double t = 0; for (int i = 0; i < 16; ++i) t += s["y"][i].get<double>(); return t; }());
      ++k;
    }
    lg->Draw();
    gPad->RedrawAxis();
  }
  c->cd();
  poster_save(c, out + "/P8_ocupacion_hod2018_hod2019");
}

// ------------------------------------------------------------------------------------------------ revisión 2: P3 y P6x para público general
// Mismos datos y números que P3 y P6 (x); solo rótulos, leyendas y tamaños (media columna, kHalf). Claves 'a' (P3) y 'b' (P6x).
void P3poster(const std::string& out) {
  auto j = readJson(R + "proceeding_revision/2026-11-03_band_legend/fig_eff_g09551_bandrunB.json");
  const auto& s = j["series"][0];
  const double rel = j["f18_rel_syst"];
  const int n = static_cast<int>(s["eps"].size());
  auto* c = poster_canvas("c_p3p");
  halfMargins();
  auto* fr = new TH1F("fr3p", ";umbral (fotoelectrones);eficiencia X\267Y", 10, 2.5, 12.5);
  applySizes(fr, kHalf);
  fr->SetMinimum(0.0);
  fr->SetMaximum(1.2);
  fr->Draw("AXIS");
  auto* band = new TGraphErrors(n), *g = new TGraphErrors(n), *gs = new TGraphErrors(n);
  for (int i = 0; i < n; ++i) {
    const double T = j["T"][i], e = s["eps"][i];
    band->SetPoint(i, T, e); band->SetPointError(i, 0.5, e * rel);
    g->SetPoint(i, T, e); g->SetPointError(i, 0, s["stat"][i].get<double>());
    gs->SetPoint(i, T, e); gs->SetPointError(i, 0.12, s["syst"][i].get<double>());
    printf("P3p: T %d eps %.4f stat %.4f estimador %.4f aceptancia %.4f\n", j["T"][i].get<int>(), e, s["stat"][i].get<double>(), s["syst"][i].get<double>(), e * rel);
  }
  band->SetFillColorAlpha(poster::grey(), 0.40); band->SetLineColor(0);
  band->Draw("E3 SAME");
  gs->SetFillStyle(0); gs->SetLineColor(poster::color(0)); gs->SetLineWidth(3);
  gs->Draw("2 SAME");
  g->SetMarkerStyle(20); g->SetMarkerColor(poster::color(0)); g->SetLineColor(poster::color(0)); g->SetLineWidth(3);
  g->Draw("P SAME");
  gPad->RedrawAxis();
  auto* lg = legendT(0.22, 0.22, 0.95, 0.45, kHalf.legend);
  lg->AddEntry(g, "medido (run B)", "pe");
  lg->AddEntry(gs, "incertidumbre del m\351todo", "f");
  lg->AddEntry(band, "incertidumbre de aceptancia", "f");
  lg->Draw();
  poster_save(c, out + "/P3_eficiencia_umbral");
}

void P6xposter(const std::string& out) {
  auto rows = readCsv(R + "hod2018_impact_map/2026-11-12_impact_map_runB/stationary_goodtime/impact_map_runB_T10.csv");
  std::vector<double> edges{-16.5};
  for (int k = 0; k < 30; ++k) edges.push_back(-14.5 + k);
  edges.push_back(16.5);
  auto* h = new TH1D("h6xp", ";posici\363n x (mm);cuentas restadas por mm", 31, edges.data());
  std::vector<int> junction, dead;
  long ntrig = 0;
  for (auto& r : rows) {
    if (r["axis"] != "x") continue;
    const int b = std::stoi(r["cell_index"]) + 1;
    h->SetBinContent(b, std::stod(r["counts_per_mm"]));
    h->SetBinError(b, std::stod(r["counts_per_mm_err"]));
    if (r["flags"].find("junction_from_cabling") != std::string::npos) junction.push_back(b);
    if (r["dead_channel_flag"] == "1") dead.push_back(b);
    ntrig = std::stol(r["triggers_denominator"]);
  }
  double ymax = 0;
  for (int b = 1; b <= 31; ++b) ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
  const double ytop = 1.75 * ymax;
  auto* c = poster_canvas("c_p6xp");
  halfMargins();
  applySizes(h, kHalf);
  h->SetMinimum(0); h->SetMaximum(ytop);
  h->SetLineColor(poster::color(5)); h->SetMarkerColor(poster::color(5)); h->SetMarkerStyle(20); h->SetLineWidth(3);
  h->Draw("E1");
  for (int b : dead) {  // celdas 23-25 sin cuentas posibles y vecinas 22, 26 en gris claro (como P6)
    const bool core = (b != dead.front() && b != dead.back());
    auto* bx = new TBox(h->GetBinLowEdge(b), 0, h->GetBinLowEdge(b + 1), ytop);
    bx->SetFillColorAlpha(poster::grey(), core ? 0.55 : 0.25);
    bx->SetLineWidth(0);
    bx->Draw();
  }
  if (!dead.empty()) {
    auto* t = new TLatex(0.5 * (h->GetBinLowEdge(dead.front() + 1) + h->GetBinLowEdge(dead.back())), 0.5 * ymax, "canal sin se\361al");
    t->SetTextAngle(90); t->SetTextAlign(22); t->SetTextSize(kHalf.legend);
    t->Draw();
  }
  double jlo = 1e9, jhi = -1e9;
  for (int b : junction) {
    auto* bj = new TBox(h->GetBinLowEdge(b), 0, h->GetBinLowEdge(b + 1), 1.12 * ymax);  // termina bajo el rótulo
    bj->SetFillStyle(3354); bj->SetFillColor(poster::color(1)); bj->SetLineWidth(0);
    bj->Draw();
    jlo = std::min(jlo, h->GetBinLowEdge(b));
    jhi = std::max(jhi, h->GetBinLowEdge(b + 1));
  }
  if (!junction.empty()) {
    auto* t = new TLatex(0.5 * (jlo + jhi), 1.16 * ymax, "uni\363n incierta");
    t->SetTextAlign(21); t->SetTextSize(kHalf.legend); t->SetTextColor(poster::color(1));
    t->Draw();
  }
  h->Draw("E1 SAME");
  gPad->RedrawAxis();
  auto* lg = legendT(0.21, 0.82, 0.60, 0.93, kHalf.legend);
  lg->AddEntry(h, "medido, 10 p.e.", "pe");
  lg->Draw();
  printf("P6xp: %ld disparos, celdas de union %zu, celdas del canal sin senal %zu\n", ntrig, junction.size(), dead.size());
  poster_save(c, out + "/P6_perfil_impacto_x");
}

void poster_figs(const char* out, const char* which = "12345678") {
  poster_style();
  setvbuf(stdout, nullptr, _IONBF, 0);
  const std::string o(out), w(which);
  if (w.find('1') != std::string::npos) P1(o);
  if (w.find('2') != std::string::npos) P2(o);
  if (w.find('3') != std::string::npos) P3(o);
  if (w.find('4') != std::string::npos) P4(o);
  if (w.find('5') != std::string::npos) P5(o);
  if (w.find('6') != std::string::npos) P6(o);
  if (w.find('7') != std::string::npos) P7(o);
  if (w.find('8') != std::string::npos) P8(o);
  if (w.find('a') != std::string::npos) P3poster(o);
  if (w.find('b') != std::string::npos) P6xposter(o);
}
