// analysis/macros/proceeding_figs_v2.C
// Generates proceeding figures v2 for run v8:
// 3a Hit times (16 ns bins in [3900, 4500] ns)
// 3b Hits per X strip in W2 [4128, 4272] ns
// 3c hod2018 X-Y map in W2 (same - mixed)
// 3d Q spectra of hod2018 channels 0, 1, 2, 4, 5, 6, 7, 8 in W2 (signal - control)
// Outputs PNG, PDF, and figdata.txt with every plotted number.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBox.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGraphErrors.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <filesystem>
#include <ctime>

namespace fs = std::filesystem;

#include "../src/hodo_consts.h"

constexpr double kPrePulse = 2500.0;

// Authoritative mapping for hod2019
static const std::array<int,16> kH19X = {53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y = {37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};

void proceeding_figs_v2(const char* input_file = "lectura_vx1740_ext_nim_1M_v8.root",
                        const char* custom_outdir = "") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    std::string outdir_str;
    if (custom_outdir && std::string(custom_outdir).length() > 0) {
        outdir_str = custom_outdir;
    } else {
        std::time_t now = std::time(nullptr);
        char ds[32];
        std::strftime(ds, sizeof(ds), "%Y%m%d_%H%M%S", std::localtime(&now));
        outdir_str = "outputs/proceeding_figs_v2_" + std::string(ds);
    }
    fs::path outdir(outdir_str);
    fs::create_directories(outdir);
    std::cout << "Writing figures to: " << outdir << "\n";

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open input file " << input_file << "\n";
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(f->Get("data_dig"));
    if (!tree) {
        std::cerr << "ERROR: data_dig tree not found in " << input_file << "\n";
        return;
    }

    std::array<std::vector<float>*, 64> cw{};
    std::vector<float>* tw = nullptr;
    tree->SetBranchStatus("*", 0);
    tree->SetCacheSize(256LL << 20);
    for (int c = 0; c < 64; ++c) {
        char bn[16]; std::snprintf(bn, sizeof(bn), "Ch%d", c);
        if (tree->GetBranch(bn)) {
            tree->SetBranchStatus(bn, 1);
            tree->SetBranchAddress(bn, &cw[c]);
        }
    }
    tree->SetBranchStatus("time", 1);
    tree->SetBranchAddress("time", &tw);

    const Long64_t nEnt = tree->GetEntries();
    std::cout << "Processing " << nEnt << " events...\n";

    // -------------------------------------------------------------
    // Data structures for 3a: Hit times [3900, 4500] ns
    // 38 bins of 16 ns in [3896, 4504] -> bin centers 3904..4496
    // -------------------------------------------------------------
    constexpr int kNBins3a = 38;
    constexpr double kTlo3a = 3896.0;
    constexpr double kThi3a = 4504.0;
    TH1D h_18x_sig("h_18x_sig", ";Hit time [ns];Hits / 16 ns", kNBins3a, kTlo3a, kThi3a);
    TH1D h_18x_ctl("h_18x_ctl", ";Hit time [ns];Hits / 16 ns", kNBins3a, kTlo3a, kThi3a);
    TH1D h_19_sig("h_19_sig", ";Hit time [ns];Hits / 16 ns", kNBins3a, kTlo3a, kThi3a);

    // -------------------------------------------------------------
    // Data structures for 3b: Hits per X strip in W2 [4128, 4272]
    // -------------------------------------------------------------
    long long n18x_s[16] = {};
    long long n18x_c[16] = {};
    long long n19x_s[16] = {};

    // -------------------------------------------------------------
    // Data structures for 3c: hod2018 X-Y map in W2
    // -------------------------------------------------------------
    TH2D h_xy_same("h_xy_same", ";hod2018 X strip;hod2018 Y strip", 16, -0.5, 15.5, 16, -0.5, 15.5);
    TH2D h_xy_mixed("h_xy_mixed", ";hod2018 X strip;hod2018 Y strip", 16, -0.5, 15.5, 16, -0.5, 15.5);

    // Per-event best strips in W2 for 3c
    struct EvBest3c {
        int best_x = -1;
        double best_x_ph = 0;
        int best_y = -1;
        double best_y_ph = 0;
    };
    std::vector<EvBest3c> ev_bests(nEnt);

    // -------------------------------------------------------------
    // Data structures for 3d: Q spectra of hod2018 channels 0,1,2,4,5,6,7,8
    // -------------------------------------------------------------
    const std::vector<int> chs3d = {0, 1, 2, 4, 5, 6, 7, 8};
    constexpr int kNBins3d = 200;
    constexpr double kQlo3d = 0.0;
    constexpr double kQhi3d = 40000.0;
    std::map<int, TH1D*> h_q_sig;
    std::map<int, TH1D*> h_q_ctl;
    for (int ch : chs3d) {
        h_q_sig[ch] = new TH1D(Form("h_q_sig_ch%d", ch), Form("Ch %d;Q [mV#upointns];Hits / 200 mV#upointns", ch), kNBins3d, kQlo3d, kQhi3d);
        h_q_ctl[ch] = new TH1D(Form("h_q_ctl_ch%d", ch), Form("Ch %d;Q [mV#upointns];Hits / 200 mV#upointns", ch), kNBins3d, kQlo3d, kQhi3d);
    }

    // Inverse map for hod2019 X strips
    std::array<int, 64> inv19x; inv19x.fill(-1);
    for (int k = 0; k < 16; ++k) inv19x[kH19X[k]] = k;

    bool checked_step = false;

    // Event loop
    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 30000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int n_samp = (int)tv.size();
        if (!checked_step && n_samp > 1) {
            checked_step = true;
            if (std::fabs((tv[1] - tv[0]) - kDt_ns) > 1e-3) {
                std::fprintf(stderr, "Sampling step %.4f ns != kDt_ns %.4f ns; aborting\n", tv[1]-tv[0], kDt_ns);
                std::exit(1);
            }
        }

        EvBest3c ev_b;

        for (int c = 0; c < 64; ++c) {
            if (!cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // Pedestal pre-pulse (t < 2500 ns)
            double s_pre = 0; int n_pre = 0;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] < kPrePulse) { s_pre += wf[i]; n_pre++; }
            }
            if (n_pre < 2) continue;
            const double ped_glob = s_pre / n_pre;
            const double thr = (c < 32) ? kPh18_mV : kPh19_mV;

            // ---------------------------------------------------------
            // 3a: Hit times [3900, 4500] ns
            // ---------------------------------------------------------
            // Signal in [3900, 4500]
            double vmin_3a = 1e30; int imin_3a = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= 3900.0 && tv[i] <= 4500.0 && wf[i] < vmin_3a) {
                    vmin_3a = wf[i]; imin_3a = i;
                }
            }
            if (imin_3a >= 0 && (ped_glob - vmin_3a) > thr) {
                if (c >= 0 && c <= 15) {
                    h_18x_sig.Fill(tv[imin_3a]);
                } else if (c >= 32 && c <= 63) {
                    h_19_sig.Fill(tv[imin_3a]);
                }
            }

            // Control in [5900, 6500] (for 18X, shifted by -2000 ns)
            if (c >= 0 && c <= 15) {
                double vmin_ctl = 1e30; int imin_ctl = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= 5900.0 && tv[i] <= 6500.0 && wf[i] < vmin_ctl) {
                        vmin_ctl = wf[i]; imin_ctl = i;
                    }
                }
                if (imin_ctl >= 0 && (ped_glob - vmin_ctl) > thr) {
                    h_18x_ctl.Fill(tv[imin_ctl] - 2000.0);
                }
            }

            // ---------------------------------------------------------
            // 3b, 3c, 3d: Window W2 [kSigTlo, kSigThi] = [4128, 4272] ns
            // Control W2 [6128, 6272] ns
            // ---------------------------------------------------------
            double vmin_w2 = 1e30; int imin_w2 = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_w2) {
                    vmin_w2 = wf[i]; imin_w2 = i;
                }
            }
            const bool hit_w2_sig = (imin_w2 >= 0 && (ped_glob - vmin_w2) > thr);

            double vmin_w2_c = 1e30; int imin_w2_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo + kCtlShift && tv[i] <= kSigThi + kCtlShift && wf[i] < vmin_w2_c) {
                    vmin_w2_c = wf[i]; imin_w2_c = i;
                }
            }
            const bool hit_w2_ctl = (imin_w2_c >= 0 && (ped_glob - vmin_w2_c) > thr);

            // 3b: Hits per X strip
            if (c >= 0 && c <= 15) {
                int k = c;
                if (hit_w2_sig) n18x_s[k]++;
                if (hit_w2_ctl) n18x_c[k]++;
            } else if (c >= 32 && c <= 63) {
                int k = inv19x[c];
                if (k >= 0 && hit_w2_sig) n19x_s[k]++;
            }

            // 3c: Best strip in 18X (0..15) and 18Y (16..31) in W2
            if (hit_w2_sig) {
                double ph = ped_glob - vmin_w2;
                if (c >= 0 && c <= 15) {
                    if (ph > ev_b.best_x_ph) {
                        ev_b.best_x_ph = ph;
                        ev_b.best_x = c;
                    }
                } else if (c >= 16 && c <= 31) {
                    if (ph > ev_b.best_y_ph) {
                        ev_b.best_y_ph = ph;
                        ev_b.best_y = c - 16;
                    }
                }
            }

            // 3d: Q spectra for hod2018 channels 0, 1, 2, 4, 5, 6, 7, 8 in W2
            if (h_q_sig.find(c) != h_q_sig.end()) {
                // Signal Q
                if (hit_w2_sig) {
                    double ped_q = ped_glob;
                    if (imin_w2 >= 12 && imin_w2 - 4 < n_samp) {
                        double s_loc = 0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin_w2 + j];
                        ped_q = s_loc / 9.0;
                    }
                    double Q = 0;
                    double tilo = tv[imin_w2] - kNpre * kDt_ns;
                    double tihi = tv[imin_w2] + kNpost18 * kDt_ns;
                    for (int i = 0; i < n_samp; ++i) {
                        if (tv[i] >= tilo && tv[i] <= tihi) Q += (ped_q - wf[i]) * kDt_ns;
                    }
                    h_q_sig[c]->Fill(Q);
                }

                // Control Q
                if (hit_w2_ctl) {
                    double ped_q = ped_glob;
                    if (imin_w2_c >= 12 && imin_w2_c - 4 < n_samp) {
                        double s_loc = 0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin_w2_c + j];
                        ped_q = s_loc / 9.0;
                    }
                    double Q = 0;
                    double tilo = tv[imin_w2_c] - kNpre * kDt_ns;
                    double tihi = tv[imin_w2_c] + kNpost18 * kDt_ns;
                    for (int i = 0; i < n_samp; ++i) {
                        if (tv[i] >= tilo && tv[i] <= tihi) Q += (ped_q - wf[i]) * kDt_ns;
                    }
                    h_q_ctl[c]->Fill(Q);
                }
            }
        }

        ev_bests[iev] = ev_b;
    }
    std::cout << "  " << nEnt << "/" << nEnt << " ✓\n";

    // -------------------------------------------------------------
    // Compute 3c: hod2018 X-Y map in W2 (same and mixed)
    // -------------------------------------------------------------
    const size_t n_acc_ev = (nEnt > 1) ? (nEnt - 1) : 0;
    for (size_t i = 0; i < n_acc_ev; ++i) {
        size_t next_i = i + 1;
        if (ev_bests[i].best_x >= 0 && ev_bests[i].best_y >= 0) {
            h_xy_same.Fill(ev_bests[i].best_x, ev_bests[i].best_y);
        }
        if (ev_bests[i].best_x >= 0 && ev_bests[next_i].best_y >= 0) {
            h_xy_mixed.Fill(ev_bests[i].best_x, ev_bests[next_i].best_y);
        }
    }
    TH2D h_xy_diff("h_xy_diff", ";hod2018 X strip;hod2018 Y strip;Hits (same - mixed)", 16, -0.5, 15.5, 16, -0.5, 15.5);
    h_xy_diff.Add(&h_xy_same, 1.0);
    h_xy_diff.Add(&h_xy_mixed, -1.0);

    // =============================================================
    // GENERATE FIGURES
    // =============================================================

    // -------------------------------------------------------------
    // Figure 3a: Hit times
    // -------------------------------------------------------------
    TCanvas c3a("c3a", "Hit times", 900, 600);
    c3a.SetLeftMargin(0.12);
    c3a.SetRightMargin(0.08);
    c3a.SetBottomMargin(0.12);
    c3a.SetTopMargin(0.08);

    h_18x_sig.SetLineColor(kBlue + 1);
    h_18x_sig.SetLineWidth(2);
    h_18x_ctl.SetLineColor(kGray + 1);
    h_18x_ctl.SetLineWidth(2);
    h_18x_ctl.SetLineStyle(2);
    h_19_sig.SetLineColor(kRed + 1);
    h_19_sig.SetLineWidth(2);

    double ymax3a = std::max({h_18x_sig.GetMaximum(), h_18x_ctl.GetMaximum(), h_19_sig.GetMaximum()}) * 1.25;
    h_18x_sig.SetMaximum(ymax3a);
    h_18x_sig.Draw("HIST");

    // Shading for W2 [4128, 4272]
    TBox box_w2(kSigTlo, 0, kSigThi, ymax3a * 0.95);
    box_w2.SetFillColorAlpha(kYellow - 9, 0.45);
    box_w2.SetLineColor(kOrange + 2);
    box_w2.SetLineStyle(3);
    box_w2.Draw("SAME");

    h_18x_sig.Draw("HIST SAME");
    h_18x_ctl.Draw("HIST SAME");
    h_19_sig.Draw("HIST SAME");

    TLegend leg3a(0.42, 0.65, 0.88, 0.88);
    leg3a.SetBorderSize(0);
    leg3a.SetFillStyle(0);
    leg3a.AddEntry(&h_18x_sig, "hod2018 X (ph > 50 mV)", "l");
    leg3a.AddEntry(&h_18x_ctl, "hod2018 X control (-2000 ns, ph > 50 mV)", "l");
    leg3a.AddEntry(&h_19_sig, "hod2019 all ch (ph > 35 mV)", "l");
    leg3a.AddEntry(&box_w2, "W2 signal window [4128, 4272] ns", "f");
    leg3a.Draw();

    c3a.SaveAs((outdir / "fig_3a_hit_times.png").string().c_str());
    c3a.SaveAs((outdir / "fig_3a_hit_times.pdf").string().c_str());

    // -------------------------------------------------------------
    // Figure 3b: Hits per X strip in W2
    // -------------------------------------------------------------
    TCanvas c3b("c3b", "Hits per X strip in W2", 900, 700);
    c3b.Divide(1, 2, 0.001, 0.001);

    // Pad 1: hod2018
    c3b.cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.10);

    std::vector<double> x_18(16), y_18(16), ex_18(16, 0.0), ey_18(16);
    for (int k = 0; k < 16; ++k) {
        x_18[k] = k;
        y_18[k] = n18x_s[k] - n18x_c[k];
        ey_18[k] = std::sqrt(n18x_s[k] + n18x_c[k]);
    }
    TGraphErrors g18(16, x_18.data(), y_18.data(), ex_18.data(), ey_18.data());
    g18.SetTitle("hod2018 X (W2);Strip k;N_{S} - N_{C} [hits]");
    g18.SetMarkerStyle(20);
    g18.SetMarkerColor(kBlue + 1);
    g18.SetLineColor(kBlue + 1);
    g18.SetLineWidth(2);
    g18.GetXaxis()->SetLimits(-0.5, 15.5);
    g18.Draw("AP");

    // Pad 2: hod2019
    c3b.cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.10);

    TH1D h19_g6("h19_g6", "hod2019 X (W2);Strip k;N_{S} [hits]", 16, -0.5, 15.5);
    TH1D h19_g7("h19_g7", "hod2019 X (W2);Strip k;N_{S} [hits]", 16, -0.5, 15.5);
    for (int k = 0; k < 16; ++k) {
        int ch = kH19X[k];
        if (ch >= 48 && ch <= 55) {
            h19_g6.SetBinContent(k + 1, n19x_s[k]);
        } else {
            h19_g7.SetBinContent(k + 1, n19x_s[k]);
        }
    }
    h19_g6.SetFillColor(kAzure + 1);
    h19_g6.SetBarWidth(0.8);
    h19_g7.SetFillColor(kOrange + 7);
    h19_g7.SetBarWidth(0.8);

    double max_19 = 0;
    for (int k = 0; k < 16; ++k) if (n19x_s[k] > max_19) max_19 = n19x_s[k];
    h19_g6.SetMaximum(max_19 * 1.3);
    h19_g6.Draw("BAR");
    h19_g7.Draw("BAR SAME");

    TLegend leg3b(0.60, 0.68, 0.88, 0.88);
    leg3b.SetBorderSize(0);
    leg3b.SetFillStyle(0);
    leg3b.AddEntry(&h19_g6, "hod2019 Group 6 (Ch 48-55)", "f");
    leg3b.AddEntry(&h19_g7, "hod2019 Group 7 (Ch 56-63)", "f");
    leg3b.Draw();

    c3b.SaveAs((outdir / "fig_3b_hits_per_strip.png").string().c_str());
    c3b.SaveAs((outdir / "fig_3b_hits_per_strip.pdf").string().c_str());

    // -------------------------------------------------------------
    // Figure 3c: hod2018 X-Y map in W2
    // -------------------------------------------------------------
    TCanvas c3c("c3c", "hod2018 X-Y map in W2", 800, 700);
    c3c.SetLeftMargin(0.12);
    c3c.SetRightMargin(0.14);
    c3c.SetBottomMargin(0.12);
    c3c.SetTopMargin(0.08);

    h_xy_diff.Draw("COLZ TEXT");
    c3c.SaveAs((outdir / "fig_3c_hod2018_xy_map.png").string().c_str());
    c3c.SaveAs((outdir / "fig_3c_hod2018_xy_map.pdf").string().c_str());

    // -------------------------------------------------------------
    // Figure 3d: Q spectra of hod2018 channels 0,1,2,4,5,6,7,8 in W2
    // -------------------------------------------------------------
    TCanvas c3d("c3d", "Q spectra in W2", 1200, 700);
    c3d.Divide(4, 2, 0.001, 0.001);

    std::map<int, TH1D*> h_q_diff;
    for (size_t idx = 0; idx < chs3d.size(); ++idx) {
        int ch = chs3d[idx];
        c3d.cd((int)idx + 1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.10);

        TH1D* h_diff = new TH1D(Form("h_q_diff_ch%d", ch),
                                Form("hod2018 Ch %d (W2);Q [mV#upointns];Hits / 200 mV#upointns", ch),
                                kNBins3d, kQlo3d, kQhi3d);
        h_diff->Add(h_q_sig[ch], 1.0);
        h_diff->Add(h_q_ctl[ch], -1.0);
        h_diff->SetLineColor(kBlue + 1);
        h_diff->SetLineWidth(2);
        h_diff->Draw("HIST");
        h_q_diff[ch] = h_diff;
    }

    c3d.SaveAs((outdir / "fig_3d_q_spectra.png").string().c_str());
    c3d.SaveAs((outdir / "fig_3d_q_spectra.pdf").string().c_str());

    // =============================================================
    // WRITE figdata.txt
    // =============================================================
    std::ofstream df(outdir / "figdata.txt");
    df << "=== FIGDATA: Proceeding Figures v2 (v8) ===\n\n";

    // 3a Data
    df << "--- 3a. Hit Times (16 ns bins in [3900, 4500] ns) ---\n";
    df << std::setw(12) << "t_center_ns"
       << std::setw(14) << "18X_sig"
       << std::setw(14) << "18X_ctl"
       << std::setw(14) << "19_sig" << "\n";
    for (int b = 1; b <= kNBins3a; ++b) {
        df << std::setw(12) << std::fixed << std::setprecision(1) << h_18x_sig.GetBinCenter(b)
           << std::setw(14) << (long long)h_18x_sig.GetBinContent(b)
           << std::setw(14) << (long long)h_18x_ctl.GetBinContent(b)
           << std::setw(14) << (long long)h_19_sig.GetBinContent(b) << "\n";
    }

    // 3b Data
    df << "\n--- 3b. Hits per X strip in W2 [4128, 4272] ns ---\n";
    df << "hod2018 X:\n";
    df << std::setw(8) << "strip"
       << std::setw(8) << "ch"
       << std::setw(12) << "N_S"
       << std::setw(12) << "N_C"
       << std::setw(14) << "N_S - N_C"
       << std::setw(14) << "err" << "\n";
    for (int k = 0; k < 16; ++k) {
        df << std::setw(8) << k
           << std::setw(8) << k
           << std::setw(12) << n18x_s[k]
           << std::setw(12) << n18x_c[k]
           << std::setw(14) << (n18x_s[k] - n18x_c[k])
           << std::setw(14) << std::fixed << std::setprecision(2) << std::sqrt(n18x_s[k] + n18x_c[k]) << "\n";
    }
    df << "\nhod2019 X:\n";
    df << std::setw(8) << "strip"
       << std::setw(8) << "ch"
       << std::setw(8) << "grp"
       << std::setw(12) << "N_S" << "\n";
    for (int k = 0; k < 16; ++k) {
        int ch = kH19X[k];
        int grp = ch / 8;
        df << std::setw(8) << k
           << std::setw(8) << ch
           << std::setw(8) << grp
           << std::setw(12) << n19x_s[k] << "\n";
    }

    // 3c Data
    double tot_same = h_xy_same.Integral();
    double tot_mixed = h_xy_mixed.Integral();
    df << "\n--- 3c. hod2018 X-Y map in W2 ---\n";
    df << "Total entries same:  " << std::fixed << std::setprecision(0) << tot_same << "\n";
    df << "Total entries mixed: " << tot_mixed << "\n";
    df << "Total same - mixed:  " << (tot_same - tot_mixed) << "\n\n";

    df << "Same - Mixed matrix (X columns 0..15, Y rows 0..15):\n";
    df << "Y\\X";
    for (int kx = 0; kx < 16; ++kx) df << std::setw(7) << kx;
    df << "\n";
    for (int ky = 0; ky < 16; ++ky) {
        df << std::setw(3) << ky;
        for (int kx = 0; kx < 16; ++kx) {
            double v = h_xy_diff.GetBinContent(kx + 1, ky + 1);
            df << std::setw(7) << (long long)v;
        }
        df << "\n";
    }

    // 3d Data
    df << "\n--- 3d. Q spectra of hod2018 channels in W2 (Signal - Control) [mV·ns] ---\n";
    df << std::setw(12) << "Q_center";
    for (int ch : chs3d) df << std::setw(10) << Form("ch%d", ch);
    df << "\n";
    for (int b = 1; b <= kNBins3d; ++b) {
        df << std::setw(12) << std::fixed << std::setprecision(1) << h_q_sig[chs3d[0]]->GetBinCenter(b);
        for (int ch : chs3d) {
            double v = h_q_diff[ch]->GetBinContent(b);
            df << std::setw(10) << (long long)v;
        }
        df << "\n";
    }

    df.close();
    std::cout << "Done! Summary written to: " << (outdir / "figdata.txt") << "\n";
    if (f) f->Close();
}

