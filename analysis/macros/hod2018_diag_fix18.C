// analysis/macros/hod2018_diag_fix18.C
// TASK: ANALYSIS-FIX-18
// Origin of extra hits in the same plane: dstrip, dch, and amplitude correlation.

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>

#include "../src/hodo_mapping.h"

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::string channels;
    std::string reason;
};

struct EventV8 {
    float ph_s[32]; // channels 0..31
    float ph_c[32];
};

struct ChannelCalib {
    int ch;
    std::string plane;
    double p0;
    double G;
    int status;
};

struct PairStats {
    long long n = 0;
    double sum_r = 0;
    double sum_r2 = 0;
    double sum_x = 0;
    double sum_x2 = 0;
    double sum_y = 0;
    double sum_y2 = 0;
    double sum_xy = 0;

    void add(double x, double y) {
        double r = y / x;
        n++;
        sum_r += r;
        sum_r2 += r * r;
        sum_x += x;
        sum_x2 += x * x;
        sum_y += y;
        sum_y2 += y * y;
        sum_xy += x * y;
    }

    double mean_r() const {
        return (n > 0) ? (sum_r / n) : 0.0;
    }

    double rms_r() const {
        if (n <= 1) return 0.0;
        double m = mean_r();
        double var = (sum_r2 / n) - (m * m);
        return (var > 0.0) ? std::sqrt(var) : 0.0;
    }

    double pearson_corr() const {
        if (n <= 1) return 0.0;
        double mx = sum_x / n;
        double my = sum_y / n;
        double sxx = sum_x2 - n * mx * mx;
        double syy = sum_y2 - n * my * my;
        double sxy = sum_xy - n * mx * my;
        if (sxx <= 0.0 || syy <= 0.0) return 0.0;
        return sxy / std::sqrt(sxx * syy);
    }
};

void hod2018_diag_fix18(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // 1. Read calibration parameters from outputs/pe_calibration_v8.txt
    ChannelCalib calib[32] = {};
    for (int c = 0; c < 32; ++c) {
        calib[c].ch = c;
        calib[c].plane = (c < 16) ? "18X" : "18Y";
        calib[c].p0 = 0.0;
        calib[c].G = 0.0;
        calib[c].status = -1;
    }

    std::ifstream pf("outputs/pe_calibration_v8.txt");
    if (!pf.is_open()) {
        std::cerr << "Error: outputs/pe_calibration_v8.txt not found!\n";
        return;
    }
    std::string line;
    while (std::getline(pf, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int ch; std::string pl; double p0, G, s0, s1, chi2; int ndf, status;
        if (iss >> ch >> pl >> p0 >> G >> s0 >> s1 >> chi2 >> ndf >> status) {
            if (ch >= 0 && ch < 32) {
                calib[ch].ch = ch;
                calib[ch].plane = pl;
                calib[ch].p0 = p0;
                calib[ch].G  = G;
                calib[ch].status = status;
            }
        }
    }
    pf.close();

    auto get_thr = [&](int c, double n_pe) -> double {
        if (c == 3) return 1e9; // dead channel
        return calib[c].p0 + (n_pe - 1.5) * calib[c].G;
    };

    // 2. Read intervals from runs/good_intervals_2026.csv
    std::vector<Interval> intervals;
    std::ifstream fcsv("runs/good_intervals_2026.csv");
    if (!fcsv.is_open()) {
        std::cerr << "Error: runs/good_intervals_2026.csv not found!\n";
        return;
    }
    std::getline(fcsv, line); // header
    while (std::getline(fcsv, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Interval iv;
        if (std::getline(iss, iv.run, ',') &&
            std::getline(iss, iv.start_local, ',') &&
            std::getline(iss, iv.end_local, ',') &&
            std::getline(iss, iv.status, ',') &&
            std::getline(iss, iv.channels, ',') &&
            std::getline(iss, iv.reason)) {
            intervals.push_back(iv);
        }
    }
    fcsv.close();

    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // -04:00 offset

    // Nominal signal window W2 and control window
    constexpr double kNomSigTlo = 4128.0, kNomSigThi = 4272.0;
    constexpr double kNomCtlTlo = 6128.0, kNomCtlThi = 6272.0;

    TFile fin(fname, "READ");
    if (fin.IsZombie()) {
        std::cerr << "Error opening " << fname << "\n";
        return;
    }
    TTree* tree = (TTree*)fin.Get("data_dig");
    if (!tree) {
        std::cerr << "Error finding tree data_dig in " << fname << "\n";
        return;
    }

    Long64_t nEnt = tree->GetEntries();
    std::vector<float>* tw = nullptr;
    Long_t timestamp = 0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("time", 1);
    tree->SetBranchAddress("time", &tw);
    tree->SetBranchStatus("timestamp", 1);
    tree->SetBranchAddress("timestamp", &timestamp);

    std::vector<float>* cw[32] = {};
    for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    std::vector<EventV8> kept_events;
    kept_events.reserve(120000);

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (!tw || tw->empty()) continue;

        double t_sec = boot_ref_epoch + ((double)timestamp) * 1e-6 + tz_offset;
        time_t t_local = (time_t)t_sec;
        struct tm tm_val; gmtime_r(&t_local, &tm_val);
        char buf_min[32];
        std::snprintf(buf_min, sizeof(buf_min), "%04d-%02d-%02d %02d:%02d",
                      tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday, tm_val.tm_hour, tm_val.tm_min);
        std::string s_min(buf_min);

        bool is_exclude = false;
        bool is_hv_off = false;
        for (const auto& iv : intervals) {
            if (iv.run == "v8") {
                if (s_min >= iv.start_local && s_min <= iv.end_local) {
                    if (iv.status == "EXCLUDE") is_exclude = true;
                    if (iv.status == "HV_OFF")  is_hv_off = true;
                }
            }
        }

        // Good time only
        if (is_exclude || is_hv_off) continue;

        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        EventV8 ev{};

        for (int c = 0; c < 32; ++c) {
            if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // Signal window W2 [4128, 4272] ns
            double vmin_s = 1e30; int imin_s = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kNomSigTlo && tv[i] <= kNomSigThi && wf[i] < vmin_s) {
                    vmin_s = wf[i]; imin_s = i;
                }
            }
            double ph_s = -1e9;
            if (imin_s >= 0) {
                double ped = 0;
                if (imin_s >= 12 && imin_s - 4 < n_samp) {
                    double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_s + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                ph_s = ped - vmin_s;
            }

            // Control window [6128, 6272] ns
            double vmin_c = 1e30; int imin_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kNomCtlTlo && tv[i] <= kNomCtlThi && wf[i] < vmin_c) {
                    vmin_c = wf[i]; imin_c = i;
                }
            }
            double ph_c = -1e9;
            if (imin_c >= 0) {
                double ped = 0;
                if (imin_c >= 12 && imin_c - 4 < n_samp) {
                    double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_c + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                ph_c = ped - vmin_c;
            }

            ev.ph_s[c] = (float)ph_s;
            ev.ph_c[c] = (float)ph_c;
        }

        kept_events.push_back(ev);
    }

    Long64_t n_v8 = (Long64_t)kept_events.size();
    std::cout << "Loaded v8: " << n_v8 << " good-time events.\n\n";

    // 2D Histograms for PNG: A_e vs A_m (0-60 p.e., 1 p.e. bins) at n = 5 p.e.
    TH2D h2_18x_ds1("h2_18x_ds1", "18X, |#Delta strip| = 1 (n = 5 p.e.);A_{m} [p.e.];A_{e} [p.e.]", 60, 0, 60, 60, 0, 60);
    TH2D h2_18x_ds3("h2_18x_ds3", "18X, |#Delta strip| #geq 3 (n = 5 p.e.);A_{m} [p.e.];A_{e} [p.e.]", 60, 0, 60, 60, 0, 60);
    TH2D h2_18y_ds1("h2_18y_ds1", "18Y, |#Delta strip| = 1 (n = 5 p.e.);A_{m} [p.e.];A_{e} [p.e.]", 60, 0, 60, 60, 0, 60);
    TH2D h2_18y_ds3("h2_18y_ds3", "18Y, |#Delta strip| #geq 3 (n = 5 p.e.);A_{m} [p.e.];A_{e} [p.e.]", 60, 0, 60, 60, 0, 60);

    const std::vector<double> n_pe_list = {5.0, 10.0};
    const std::vector<std::string> plane_names = {"18X", "18Y"};

    for (double n_pe : n_pe_list) {
        for (const std::string& pl_name : plane_names) {
            bool is_x = (pl_name == "18X");
            int ch_start = is_x ? 0 : 16;
            int ch_end   = is_x ? 16 : 32;

            long long cnt_sig_dstrip[31] = {};
            long long cnt_ctl_dstrip[31] = {};
            long long cnt_sig_dch[31]    = {};
            long long cnt_ctl_dch[31]    = {};

            PairStats stats_ds1;
            PairStats stats_ds3;
            PairStats stats_dch1_ds3;

            long long n_sel_sig = 0;
            long long n_sel_ctl = 0;
            long long n_mult_sig = 0;
            long long n_mult_ctl = 0;
            long long n_pairs_sig = 0;
            long long n_pairs_ctl = 0;

            for (Long64_t i = 0; i < n_v8; ++i) {
                // Determine hits across 18X and 18Y
                std::vector<int> hits_sig_x, hits_sig_y;
                std::vector<int> hits_ctl_x, hits_ctl_y;

                for (int c = 0; c < 32; ++c) {
                    if (c == 3) continue;
                    double thr = get_thr(c, n_pe);
                    if (kept_events[i].ph_s[c] > thr) {
                        if (c < 16) hits_sig_x.push_back(c);
                        else        hits_sig_y.push_back(c);
                    }
                    if (kept_events[i].ph_c[c] > thr) {
                        if (c < 16) hits_ctl_x.push_back(c);
                        else        hits_ctl_y.push_back(c);
                    }
                }

                // Signal selection: nX_sig >= 1 AND nY_sig >= 1
                if (!hits_sig_x.empty() && !hits_sig_y.empty()) {
                    n_sel_sig++;
                    const auto& plane_hits = is_x ? hits_sig_x : hits_sig_y;
                    if (plane_hits.size() >= 2) n_mult_sig++;

                    // Find m: channel with largest amplitude in p.e. in this plane
                    int m = -1;
                    double max_pe = -1e9;
                    for (int c : plane_hits) {
                        double pe = (kept_events[i].ph_s[c] - calib[c].p0) / calib[c].G;
                        if (pe > max_pe) {
                            max_pe = pe;
                            m = c;
                        }
                    }

                    if (m >= 0) {
                        double A_m = max_pe;
                        int strip_m = ch_to_strip(m);

                        for (int e : plane_hits) {
                            if (e == m) continue;
                            n_pairs_sig++;
                            double A_e = (kept_events[i].ph_s[e] - calib[e].p0) / calib[e].G;
                            int strip_e = ch_to_strip(e);
                            int dstrip = strip_e - strip_m;
                            int dch = e - m;

                            if (dstrip >= -15 && dstrip <= 15) {
                                cnt_sig_dstrip[dstrip + 15]++;
                            }
                            if (dch >= -15 && dch <= 15) {
                                cnt_sig_dch[dch + 15]++;
                            }

                            // Correlations for signal window
                            int abs_ds = std::abs(dstrip);
                            int abs_dch = std::abs(dch);

                            if (abs_ds == 1) {
                                stats_ds1.add(A_m, A_e);
                                if (n_pe == 5.0) {
                                    if (is_x) h2_18x_ds1.Fill(A_m, A_e);
                                    else      h2_18y_ds1.Fill(A_m, A_e);
                                }
                            }
                            if (abs_ds >= 3) {
                                stats_ds3.add(A_m, A_e);
                                if (n_pe == 5.0) {
                                    if (is_x) h2_18x_ds3.Fill(A_m, A_e);
                                    else      h2_18y_ds3.Fill(A_m, A_e);
                                }
                            }
                            if (abs_dch == 1 && abs_ds >= 3) {
                                stats_dch1_ds3.add(A_m, A_e);
                            }
                        }
                    }
                }

                // Control selection: nX_ctl >= 1 AND nY_ctl >= 1
                if (!hits_ctl_x.empty() && !hits_ctl_y.empty()) {
                    n_sel_ctl++;
                    const auto& plane_hits = is_x ? hits_ctl_x : hits_ctl_y;
                    if (plane_hits.size() >= 2) n_mult_ctl++;

                    // Find m: channel with largest amplitude in p.e. in this plane
                    int m = -1;
                    double max_pe = -1e9;
                    for (int c : plane_hits) {
                        double pe = (kept_events[i].ph_c[c] - calib[c].p0) / calib[c].G;
                        if (pe > max_pe) {
                            max_pe = pe;
                            m = c;
                        }
                    }

                    if (m >= 0) {
                        int strip_m = ch_to_strip(m);

                        for (int e : plane_hits) {
                            if (e == m) continue;
                            n_pairs_ctl++;
                            int strip_e = ch_to_strip(e);
                            int dstrip = strip_e - strip_m;
                            int dch = e - m;

                            if (dstrip >= -15 && dstrip <= 15) {
                                cnt_ctl_dstrip[dstrip + 15]++;
                            }
                            if (dch >= -15 && dch <= 15) {
                                cnt_ctl_dch[dch + 15]++;
                            }
                        }
                    }
                }
            }

            // Report output for this (n_pe, plane)
            std::cout << "====================================================================================================\n";
            std::cout << "THRESHOLD n = " << (int)n_pe << " p.e., PLANE " << pl_name << "\n";
            std::cout << "====================================================================================================\n";

            // (e) Number of selected events
            std::cout << "(e) Number of selected events:\n";
            std::cout << "  Selected events (nX >= 1 AND nY >= 1): count_sig = " << n_sel_sig
                      << ", count_ctl = " << n_sel_ctl << "\n";
            std::cout << "  Events with >= 2 hits in " << pl_name << ": count_sig = " << n_mult_sig
                      << ", count_ctl = " << n_mult_ctl << "\n";
            std::cout << "  Total extra-hit entries in " << pl_name << ": count_sig = " << n_pairs_sig
                      << ", count_ctl = " << n_pairs_ctl << "\n\n";

            // (a) Table of dstrip (-15..15)
            std::cout << "(a) Table of dstrip (-15..15):\n";
            std::cout << std::left << std::setw(8) << "dstrip"
                      << std::right << std::setw(12) << "count_sig"
                      << std::setw(12) << "count_ctl" << "\n";
            for (int d = -15; d <= 15; ++d) {
                std::cout << std::left << std::setw(8) << d
                          << std::right << std::setw(12) << cnt_sig_dstrip[d + 15]
                          << std::setw(12) << cnt_ctl_dstrip[d + 15] << "\n";
            }
            std::cout << "\n";

            // (b) Table of dch (-15..15)
            std::cout << "(b) Table of dch (-15..15):\n";
            std::cout << std::left << std::setw(8) << "dch"
                      << std::right << std::setw(12) << "count_sig"
                      << std::setw(12) << "count_ctl" << "\n";
            for (int d = -15; d <= 15; ++d) {
                std::cout << std::left << std::setw(8) << d
                          << std::right << std::setw(12) << cnt_sig_dch[d + 15]
                          << std::setw(12) << cnt_ctl_dch[d + 15] << "\n";
            }
            std::cout << "\n";

            // (c) Amplitude correlation by dstrip (signal window)
            std::cout << "(c) Signal window amplitude correlation by dstrip:\n";
            std::cout << "|dstrip| = 1:  entries = " << stats_ds1.n
                      << ", mean(A_e/A_m) = " << std::fixed << std::setprecision(4) << stats_ds1.mean_r()
                      << ", RMS(A_e/A_m) = " << stats_ds1.rms_r()
                      << ", Pearson r = " << stats_ds1.pearson_corr() << "\n";
            std::cout << "|dstrip| >= 3: entries = " << stats_ds3.n
                      << ", mean(A_e/A_m) = " << std::fixed << std::setprecision(4) << stats_ds3.mean_r()
                      << ", RMS(A_e/A_m) = " << stats_ds3.rms_r()
                      << ", Pearson r = " << stats_ds3.pearson_corr() << "\n\n";

            // (d) Amplitude correlation by dch (|dch| = 1 restricted to |dstrip| >= 3)
            std::cout << "(d) Signal window amplitude correlation for |dch| = 1 restricted to |dstrip| >= 3:\n";
            std::cout << "|dch| = 1 (|dstrip| >= 3): entries = " << stats_dch1_ds3.n
                      << ", mean(A_e/A_m) = " << std::fixed << std::setprecision(4) << stats_dch1_ds3.mean_r()
                      << ", RMS(A_e/A_m) = " << stats_dch1_ds3.rms_r()
                      << ", Pearson r = " << stats_dch1_ds3.pearson_corr() << "\n\n";
        }
    }

    // Save PNG outputs/fix18_Ae_vs_Am.png
    TCanvas c1("c1", "A_e vs A_m (n = 5 p.e.)", 1100, 950);
    c1.Divide(2, 2);

    c1.cd(1);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h2_18x_ds1.Draw("COLZ");

    c1.cd(2);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h2_18x_ds3.Draw("COLZ");

    c1.cd(3);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h2_18y_ds1.Draw("COLZ");

    c1.cd(4);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h2_18y_ds3.Draw("COLZ");

    c1.SaveAs("outputs/fix18_Ae_vs_Am.png");
    std::cout << "Saved figure: outputs/fix18_Ae_vs_Am.png\n";
}
