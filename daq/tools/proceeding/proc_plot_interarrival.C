// proc_plot_interarrival.C — figura «trigger_interarrival» del proceeding Muographers 2026 (figures/data/trigger_interarrival.{png,pdf}).
// Lee la rama timestamp (µs) del árbol data_dig de run B, lectura_vx1740_ext_nim_1M_v8.root (143 814 triggers), solo en modo READ.
// Escribe <base>.png (300 ppp) y <base>.pdf (vectorial): ΔT/⟨ΔT⟩ sin los 15 huecos de DAQ (> 50 T/(N−1)), 100 bines en [0, 10], contra exp(−x).
// Imprime λ y p de KS (como fit_exp_ks() de analysis/src/poisson_analysis.cxx), ⟨ΔT⟩ usada y la pendiente de ln(densidad) en x ∈ [0.5, 6].
// Necesita proc_style.C: root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_interarrival.C("<v8.root>","<base>")'
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

void proc_plot_interarrival(const char* fname, const char* base) {
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) { std::printf("ERROR abriendo %s\n", fname); return; }
  TTree* t = static_cast<TTree*>(f->Get("data_dig"));
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("timestamp", 1);
  Long_t ts_val = 0;
  t->SetBranchAddress("timestamp", &ts_val);
  const long long N = t->GetEntries();
  std::vector<double> ts(N);
  for (long long i = 0; i < N; i++) { t->GetEntry(i); ts[i] = static_cast<double>(ts_val); }
  f->Close();

  std::sort(ts.begin(), ts.end());
  const double t0 = ts[0];
  for (auto& v : ts) v -= t0;
  const double T_total = ts.back();            // µs
  const double lam_hat = (N - 1.0) / T_total;  // µs⁻¹
  const double t_mean = 1.0 / lam_hat;         // µs
  std::vector<double> dt_clean;
  int n_gaps = 0;
  for (long long i = 0; i < N - 1; i++) {
    const double dt = ts[i + 1] - ts[i];
    if (dt < 50 * t_mean) dt_clean.push_back(dt); else ++n_gaps;
  }
  const double tau_dead = *std::min_element(dt_clean.begin(), dt_clean.end());

  // fit_exp_ks() de poisson_analysis.cxx
  const int nbins = 200;
  const double range = 10.0 * t_mean, bin_w = range / nbins;
  const double fit_lo = std::max(tau_dead * 1.05, t_mean * 0.02), fit_hi = 5.0 * t_mean;
  TH1D hfit("_hfit_", "", nbins, 0, range);
  for (double v : dt_clean) if (v < range) hfit.Fill(v);
  TF1 fE("_fE_", "[0]*exp(-[1]*x)", fit_lo, fit_hi);
  fE.SetParameters(dt_clean.size() * bin_w * lam_hat, lam_hat);
  fE.SetParLimits(1, lam_hat * 0.05, lam_hat * 20.0);
  const int st = hfit.Fit(&fE, "RQN");
  const double lam_fit = fE.GetParameter(1), lam_fit_err = fE.GetParError(1);
  std::vector<double> dts = dt_clean;
  std::sort(dts.begin(), dts.end());
  const int n = static_cast<int>(dts.size());
  double D = 0;
  for (int i = 0; i < n; i++) {
    const double F = 1.0 - std::exp(-lam_fit * dts[i]);
    D = std::max(D, std::max(std::abs((i + 1.0) / n - F), std::abs(double(i) / n - F)));
  }
  const double pKS = TMath::KolmogorovProb(D * std::sqrt(n));
  std::printf("file %s\n", fname);
  std::printf("N triggers %lld, intervals %lld, DAQ gaps (dT > 50 t_mean) %d, intervals used %zu\n", N, N - 1, n_gaps, dt_clean.size());
  std::printf("duration %.2f h, lambda_hat %.4f Hz, t_mean %.1f us, dead time %.0f us\n", T_total / 3.6e9, lam_hat * 1e6, t_mean, tau_dead);
  std::printf("fit status %d: lambda_fit = %.4f +- %.4f Hz (%.3f +- %.3f Hz)\n", st, lam_fit * 1e6, lam_fit_err * 1e6, lam_fit * 1e6, lam_fit_err * 1e6);
  std::printf("KS vs exp(lambda_fit): D = %.5f, p = %.4f (%.3f)\n", D, pKS, pKS);

  // panel (gen_beamer_figs.C, figura C1, panel 1), pero en unidades de la media de los intervalos histogramados (sin huecos)
  double dt_sum = 0;
  for (double v : dt_clean) dt_sum += v;
  const double dt_mean = dt_sum / dt_clean.size();  // µs
  std::printf("mean used: <dT> of the %zu intervals without gaps = %.1f us = %.6f s (1/<dT> = %.4f Hz); mean with gaps T/(N-1) = %.1f us\n",
              dt_clean.size(), dt_mean, dt_mean * 1e-6, 1e6 / dt_mean, t_mean);
  TH1D h("h_dt", "", 100, 0, 10);
  long long n_out = 0;
  for (double v : dt_clean) { const double x = v / dt_mean; if (x < 10) h.Fill(x); else ++n_out; }
  std::printf("histogram: %.0f entries in [0,10), %lld above 10\n", h.Integral(), n_out);
  // recta a ln(densidad) en los bines con centro en [0.5, 6]; error de cada bin 1/sqrt(cuentas)
  const double norm = 1.0 / h.Integral() / (10.0 / 100);
  TGraphErrors g;
  for (int b = 1; b <= h.GetNbinsX(); b++) {
    const double xc = h.GetBinCenter(b), cnt = h.GetBinContent(b);
    if (xc < 0.5 || xc > 6.0 || cnt <= 0) continue;
    const int k = g.GetN();
    g.SetPoint(k, xc, std::log(cnt * norm));
    g.SetPointError(k, 0, 1.0 / std::sqrt(cnt));
  }
  TF1 fl("fl", "pol1", 0.5, 6.0);
  g.Fit(&fl, "QN");
  std::printf("line fit to ln(density), x in [0.5, 6] (%d bins): slope = %.4f +- %.4f, intercept = %.4f +- %.4f, chi2/ndf = %.1f/%d\n", g.GetN(),
              fl.GetParameter(1), fl.GetParError(1), fl.GetParameter(0), fl.GetParError(0), fl.GetChisquare(), fl.GetNDF());
  h.Scale(norm);
  auto* c = proc_canvas("c_ia", 1, 0.75);
  c->SetLogy();
  c->SetLeftMargin(0.16);
  c->SetBottomMargin(0.17);
  h.SetLineColor(proc::color(5));
  h.GetXaxis()->SetTitle("Time between triggers / mean");
  h.GetYaxis()->SetTitle("Density");
  h.SetMinimum(0.5 * h.GetMinimum(0));
  h.SetMaximum(2.0);
  h.Draw("HIST");
  TF1 fth("fth", "exp(-x)", 0, 10);
  fth.SetNpx(40);
  fth.SetLineColor(proc::color(6));
  fth.SetLineWidth(3);
  fth.SetLineStyle(9);
  fth.Draw("SAME");
  TLegend lg(0.56, 0.76, 0.94, 0.94);
  lg.AddEntry(&h, "data (run B)", "l");
  lg.AddEntry(&fth, "exp(-x)", "l");
  lg.Draw();
  proc_save(c, base);
}
