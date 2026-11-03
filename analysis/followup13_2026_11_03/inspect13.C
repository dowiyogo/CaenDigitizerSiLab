// inspect13.C — seguimiento 2026-11-03 (13): figuras de inspección (no para el artículo) desde las tablas de time_res.py y overlap_ratio.py.
//  A: Δt = t_X − t_Y neto (mismo − mezclados) por bin de la menor amplitud, cfd50 con desplazamientos y propagación, con la gaussiana de
//     time_res.csv; kmin (todos); σ por tira contra el bin, cfd50 / cfd50_bl / kmin (RMS).
//  B: R = a1/(a1+a2) neto en X e Y, R dentro del par, MC R_path, MC fotoestadística (pois) y MC con espectro medido (spec), normalizados;
//     u_mid del MC.
// Uso (desde la raíz): root -l -b -q 'analysis/followup13_2026_11_03/inspect13.C("<dir A>", "<dir B>")'
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>

#include <cmath>
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
    if (line.empty() || line[0] == '#' || line.rfind("\"#", 0) == 0) continue;  // comentario (entre comillas si tiene comas)
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
double num(const std::string& x) {
  try {
    return std::stod(x);
  } catch (...) {
    return NAN;
  }
}
}  // namespace

void inspect13(const char* dirA, const char* dirB) {
  gStyle->SetOptStat(0);
  const std::string A(dirA), B(dirB);
  auto res = readCsv(A + "/time_res.csv");
  auto hist = readCsv(A + "/time_res_hist.csv");
  auto* c = new TCanvas("cA", "", 1500, 1000);
  c->Divide(3, 2);
  const std::vector<std::string> bins{"10-15", "15-20", "20-30", ">30"};
  int pad = 1;
  for (const auto& b : bins) {
    c->cd(pad++);
    auto* h = new TH1D(Form("h%d", pad), Form("May run B: #Deltat = t_{X} - t_{Y}, CFD 0.5, offsets + prop., min amp %s p.e.;#Deltat (ns);net / 1 ns", b.c_str()),
                       160, -80, 80);
    for (auto& r : hist)
      if (r["method"] == "cfd50" && r["correction"] == "offsets+prop" && r["min_amp_bin_pe"] == b) h->Fill(0.5 * (num(r["lo_ns"]) + num(r["hi_ns"])), num(r["n_net"]));
    h->GetXaxis()->SetRangeUser(-20, 20);
    h->SetLineColor(kBlack);
    h->Draw("HIST");
    for (auto& r : res)
      if (r["method"] == "cfd50" && r["correction"] == "offsets+prop" && r["min_amp_bin_pe"] == b) {
        auto* f = new TF1(Form("g%d", pad), "gaus", -20, 20);
        const double s = num(r["sigma_ns"]), mu = num(r["mu_ns"]);
        double amp = 0;
        for (int i = 1; i <= h->GetNbinsX(); ++i)
          if (std::fabs(h->GetBinCenter(i) - mu) < 2.5 * s) amp = std::max(amp, h->GetBinContent(i));
        f->SetParameters(amp, mu, s);
        f->SetLineColor(kRed + 1);
        f->Draw("SAME");
        auto* lg = new TLegend(0.55, 0.75, 0.88, 0.88);
        lg->AddEntry(f, Form("#sigma = %.2f ns", s), "l");
        lg->Draw();
      }
  }
  c->cd(5);
  {
    auto* h = new TH1D("hk", "May run B: #Deltat from sample of minimum (all amplitudes);#Deltat (ns);net / 16 ns", 10, -88, 72);
    for (auto& r : hist)
      if (r["method"] == "kmin" && r["correction"] == "none" && r["min_amp_bin_pe"] == "all") h->Fill(0.5 * (num(r["lo_ns"]) + num(r["hi_ns"])), num(r["n_net"]));
    h->SetLineColor(kBlue + 1);
    h->Draw("HIST");
  }
  c->cd(6);
  {
    auto* mg = new TMultiGraph();
    mg->SetTitle("#sigma per strip (= #sigma/#sqrt{2}) vs smaller amplitude;bin (1: 10-15, 2: 15-20, 3: 20-30, 4: > 30 p.e.);#sigma_{t} per strip (ns)");
    auto* lg = new TLegend(0.40, 0.65, 0.88, 0.88);
    int k = 0;
    for (const char* m : {"cfd50", "cfd50_bl", "kmin"}) {
      auto* g = new TGraphErrors();
      for (size_t i = 0; i < bins.size(); ++i)
        for (auto& r : res)
          if (r["method"] == m && r["correction"] == (std::string(m) == "kmin" ? "none" : "offsets+prop") && r["min_amp_bin_pe"] == bins[i]) {
            const int n = g->GetN();
            const bool km = std::string(m) == "kmin";
            g->SetPoint(n, i + 1 + 0.08 * (k - 1), num(r[km ? "rms_strip_ns" : "sigma_strip_ns"]));
            g->SetPointError(n, 0, km ? 0 : num(r["sigma_strip_err_ns"]));
          }
      const int col[3] = {kRed + 1, kGreen + 2, kBlue + 1};
      g->SetMarkerStyle(20 + k);
      g->SetMarkerColor(col[k]);
      g->SetLineColor(col[k]);
      mg->Add(g, "P");
      lg->AddEntry(g, std::string(m) == "kmin" ? "sample of minimum (RMS/#sqrt{2})" : Form("%s, Gaussian #sigma/#sqrt{2}", m), "p");
      ++k;
    }
    mg->Draw("A");
    mg->SetMinimum(0);
    mg->SetMaximum(8);
    lg->Draw();
  }
  c->SaveAs((A + "/time_res_inspect.png").c_str());
  c->SaveAs((A + "/time_res_inspect.pdf").c_str());

  auto ov = readCsv(B + "/overlap_ratio.csv");
  auto* d = new TCanvas("cB", "", 1500, 500);
  d->Divide(3, 1);
  int q = 1;
  for (const char* co : {"X", "Y"}) {
    d->cd(q++);
    auto* lg = new TLegend(0.12, 0.62, 0.50, 0.88);
    struct S {
      const char* series;
      const char* top;
      const char* col;
      int color;
      const char* label;
    };
    const std::vector<S> ss{{"data", "", "n_net", kBlack, "data, net"},
                            {"mc_path", "X", "n_same", kGray + 2, "MC L_{1}/(L_{1}+L_{2}) (peak at 0.5 off scale)"},
                            {"mc_pois", "X", "n_same", kRed + 1, "MC + photostatistics (13 p.e.)"},
                            {"mc_spec", "X", "n_same", kBlue + 1, "MC + measured spectrum"}};
    int k = 0;
    for (const auto& s : ss) {
      auto* h = new TH1D(Form("r%s%d", co, k), Form("May run B, %s: a_{1}/(a_{1}+a_{2}), adjacent strips > 5 p.e., tag > 10 p.e.;ratio;normalised / 0.025", co), 40, 0, 1);
      for (auto& r : ov)
        if (r["coord"] == co && r["series"] == s.series && (std::string(s.top).empty() || r["mc_top"] == s.top))
          h->Fill(0.5 * (num(r["lo"]) + num(r["hi"])), num(r[s.col]));
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(s.color);
      h->SetLineWidth(2);
      h->SetMaximum(0.30);
      h->Draw(k == 0 ? "HIST" : "HIST SAME");
      lg->AddEntry(h, s.label, "l");
      ++k;
    }
    lg->Draw();
  }
  d->cd(3);
  {
    auto* h = new TH1D("u", "MC: crossing at the mid-plane between the layers, relative to the overlap centre;u_{mid} (mm);fraction / 0.1 mm", 60, -3, 3);
    for (auto& r : ov)
      if (r["coord"] == "X" && r["series"] == "mc_u_mid" && r["mc_top"] == "X") h->Fill(num(r["lo"]) + 0.05, num(r["n_same"]));
    h->Scale(1.0 / h->Integral());
    h->SetLineColor(kBlack);
    h->Draw("HIST");
  }
  d->SaveAs((B + "/overlap_ratio_inspect.png").c_str());
  d->SaveAs((B + "/overlap_ratio_inspect.pdf").c_str());
}
