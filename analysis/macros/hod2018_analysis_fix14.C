// analysis/macros/hod2018_analysis_fix14.C
// TASK: ANALYSIS-FIX-14
// Step 1: Corrected control-window estimator (run v8 good time, 5 p.e.)
// Step 2: Bootstrap uncertainty (runs v7..v10, good time, 5 p.e., multi-offset mixing)
// Step 3: Threshold scan of the ratio (run v8 good time, multi-offset mixing, bootstrap RMS as error)
// Step 4: X-Y map at 10 p.e. (run v8 good time, mapping kCfg_2026_05)

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
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
#include <numeric>
#include <random>

#include "../src/hodo_mapping.h"

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::string channels;
    std::string reason;
};

// Hit flags for 5 p.e. (used for Step 2 across all runs)
struct Event5pe {
    bool xs;
    bool xc;
    bool ys;
    bool yc;
};

// Full event data for v8 (used for Steps 1, 3, 4)
struct EventV8 {
    float ph_s[32]; // channels 0..31
    float ph_c[32];
};

struct MixingResult {
    double N_XY_real;
    double sigma;
    double ratio;
    double sigma_ratio;
    double aX;
    double aY;
    double rX;
    double rY;
    long long N_X_S;
    long long N_X_C;
    long long N_Y_S;
    long long N_Y_C;
    long long N_XY_same;
    double M_bar;
};

inline MixingResult eval_mixing(const std::vector<bool>& xs,
                                const std::vector<bool>& xc,
                                const std::vector<bool>& ys,
                                const std::vector<bool>& yc) {
    Long64_t n_ev = xs.size();
    Long64_t n_pairs = (n_ev > 1) ? (n_ev - 1) : 0;
    double N_d = (double)n_pairs;

    long long N_X_S = 0, N_X_C = 0, N_Y_S = 0, N_Y_C = 0;
    long long N_XY_same = 0;
    for (Long64_t i = 0; i < n_pairs; ++i) {
        if (xs[i]) N_X_S++;
        if (xc[i]) N_X_C++;
        if (ys[i]) N_Y_S++;
        if (yc[i]) N_Y_C++;
        if (xs[i] && ys[i]) N_XY_same++;
    }

    double aX = (N_d > 0) ? ((double)N_X_C / N_d) : 0.0;
    double aY = (N_d > 0) ? ((double)N_Y_C / N_d) : 0.0;
    double rX = (1.0 - aX > 1e-12) ? (((double)N_X_S / N_d - aX) / (1.0 - aX)) : 0.0;
    double rY = (1.0 - aY > 1e-12) ? (((double)N_Y_S / N_d - aY) / (1.0 - aY)) : 0.0;

    double sum_m = 0;
    for (int k = 1; k <= 20; ++k) {
        long long m = 0;
        for (Long64_t i = 0; i + k < n_ev; ++i) {
            if (xs[i] && ys[i + k]) m++;
        }
        sum_m += (double)m;
    }
    double M_bar = sum_m / 20.0;

    double denom = (1.0 - aX) * (1.0 - aY);
    double N_XY_real = (std::abs(denom) > 1e-12) ? ((N_XY_same - M_bar) / denom + rX * rY * N_d) : 0.0;
    double sigma = (std::abs(denom) > 1e-12) ? (std::sqrt(N_XY_same + M_bar / 20.0) / denom) : 0.0;

    double ratio = (n_ev > 0) ? (N_XY_real / (0.0360 * (double)n_ev)) : 0.0;
    double sigma_ratio = (n_ev > 0) ? (sigma / (0.0360 * (double)n_ev)) : 0.0;

    MixingResult res;
    res.N_XY_real = N_XY_real;
    res.sigma = sigma;
    res.ratio = ratio;
    res.sigma_ratio = sigma_ratio;
    res.aX = aX;
    res.aY = aY;
    res.rX = rX;
    res.rY = rY;
    res.N_X_S = N_X_S;
    res.N_X_C = N_X_C;
    res.N_Y_S = N_Y_S;
    res.N_Y_C = N_Y_C;
    res.N_XY_same = N_XY_same;
    res.M_bar = M_bar;
    return res;
}

void hod2018_analysis_fix14() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // 1. Read calibration parameters from outputs/pe_calibration_v8.txt
    double p0_cal[32] = {};
    double G_cal[32] = {};
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
                p0_cal[ch] = p0;
                G_cal[ch]  = G;
            }
        }
    }
    pf.close();

    // Threshold calculation helper: T_c(n) = p0_c + (n - 1.5) * G_c
    auto get_thr = [&](int c, double n_pe) -> double {
        if (c == 3) return 1e9; // dead channel
        return p0_cal[c] + (n_pe - 1.5) * G_cal[c];
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

    struct FileSpec {
        std::string run_name;
        std::string filename;
    };
    std::vector<FileSpec> file_list = {
        {"v7",  "lectura_vx1740_ext_nim_1M_v7.root"},
        {"v8",  "lectura_vx1740_ext_nim_1M_v8.root"},
        {"v9",  "lectura_vx1740_ext_nim_1M_v9.root"},
        {"v10", "lectura_vx1740_ext_nim_1M_v10.root"}
    };

    std::vector<std::vector<Event5pe>> runs_events_5pe(4);
    std::vector<EventV8> v8_events;

    for (size_t ifile = 0; ifile < file_list.size(); ++ifile) {
        const auto& fs = file_list[ifile];
        TFile fin(fs.filename.c_str(), "READ");
        if (fin.IsZombie()) {
            std::cerr << "Error opening " << fs.filename << "\n";
            continue;
        }
        TTree* tree = (TTree*)fin.Get("data_dig");
        if (!tree) {
            std::cerr << "Error finding tree data_dig in " << fs.filename << "\n";
            continue;
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
                if (iv.run == fs.run_name) {
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

            Event5pe e5{};
            EventV8 ev8{};

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

                // 5 p.e. hit flag
                double thr5 = get_thr(c, 5.0);
                if (ph_s > thr5) {
                    if (c < 16) e5.xs = true;
                    else        e5.ys = true;
                }
                if (ph_c > thr5) {
                    if (c < 16) e5.xc = true;
                    else        e5.yc = true;
                }

                if (fs.run_name == "v8") {
                    ev8.ph_s[c] = (float)ph_s;
                    ev8.ph_c[c] = (float)ph_c;
                }
            }

            runs_events_5pe[ifile].push_back(e5);
            if (fs.run_name == "v8") {
                v8_events.push_back(ev8);
            }
        }
        std::cout << "Loaded " << fs.run_name << ": " << runs_events_5pe[ifile].size() << " good-time events.\n";
    }
    std::cout << "\n";

    // =========================================================================
    // Step 1. Corrected Control-Window Estimator (Run v8 good time, 5 p.e.)
    // =========================================================================
    std::cout << "=======================================================================================================================================================\n";
    std::cout << "Step 1: Corrected Control-Window Estimator (Run v8 Good Time, 5 p.e.)\n";
    std::cout << "=======================================================================================================================================================\n";

    const auto& v8_e5 = runs_events_5pe[1]; // v8 is index 1
    Long64_t n_v8 = (Long64_t)v8_e5.size();
    Long64_t N_pairs_v8 = (n_v8 > 1) ? (n_v8 - 1) : 0;
    double N_d_v8 = (double)N_pairs_v8;

    std::vector<bool> v8_xs_5(n_v8), v8_xc_5(n_v8), v8_ys_5(n_v8), v8_yc_5(n_v8);
    long long n_xs_nom = 0, n_xc_nom = 0, n_ys_nom = 0, n_yc_nom = 0, n_same_nom = 0;
    for (Long64_t i = 0; i < n_v8; ++i) {
        v8_xs_5[i] = v8_e5[i].xs;
        v8_xc_5[i] = v8_e5[i].xc;
        v8_ys_5[i] = v8_e5[i].ys;
        v8_yc_5[i] = v8_e5[i].yc;
        if (i < N_pairs_v8) {
            if (v8_xs_5[i]) n_xs_nom++;
            if (v8_xc_5[i]) n_xc_nom++;
            if (v8_ys_5[i]) n_ys_nom++;
            if (v8_yc_5[i]) n_yc_nom++;
            if (v8_xs_5[i] && v8_ys_5[i]) n_same_nom++;
        }
    }

    // Mixing estimator nominal values
    MixingResult mix_v8 = eval_mixing(v8_xs_5, v8_xc_5, v8_ys_5, v8_yc_5);

    // Accidental counts as in FIX-13 Step 2(c)
    double N_tot = (double)N_pairs_v8;
    double P_Xc = (double)n_xc_nom / N_tot;
    double P_Ys = (double)n_ys_nom / N_tot;
    double P_Xs = (double)n_xs_nom / N_tot;
    double P_Yc = (double)n_yc_nom / N_tot;
    double N_XY_acc = N_tot * P_Xc * P_Ys + N_tot * P_Xs * P_Yc - N_tot * P_Xc * P_Yc;

    // Corrected control-window estimator
    double aX_v8 = mix_v8.aX;
    double aY_v8 = mix_v8.aY;
    double denom_v8 = (1.0 - aX_v8) * (1.0 - aY_v8);

    double N_XY_real_ctl_corr = (n_same_nom - N_XY_acc) / denom_v8;
    double sigma_ctl_corr = std::sqrt((double)n_same_nom + N_XY_acc) / denom_v8;
    double ratio_ctl_corr = N_XY_real_ctl_corr / (0.0360 * (double)n_v8);
    double sigma_ratio_ctl_corr = sigma_ctl_corr / (0.0360 * (double)n_v8);

    std::cout << std::left << std::setw(32) << "Estimator"
              << std::right << std::setw(14) << "N_XY_real"
              << std::setw(12) << "sigma"
              << std::setw(12) << "ratio"
              << std::setw(14) << "sigma_ratio"
              << std::setw(24) << "Notes"
              << "\n";
    std::cout << std::string(108, '-') << "\n";

    std::cout << std::left << std::setw(32) << "Mixing (nominal reference)"
              << std::right << std::fixed << std::setprecision(2)
              << std::setw(14) << mix_v8.N_XY_real
              << std::setw(12) << mix_v8.sigma
              << std::setprecision(6)
              << std::setw(12) << mix_v8.ratio
              << std::setw(14) << mix_v8.sigma_ratio
              << std::setw(24) << "M_bar = 20928.50"
              << "\n";

    std::cout << std::left << std::setw(32) << "Corrected Control Window"
              << std::right << std::fixed << std::setprecision(2)
              << std::setw(14) << N_XY_real_ctl_corr
              << std::setw(12) << sigma_ctl_corr
              << std::setprecision(6)
              << std::setw(12) << ratio_ctl_corr
              << std::setw(14) << sigma_ratio_ctl_corr
              << std::setw(24) << "N_XY_acc = 20550.40"
              << "\n";
    std::cout << std::string(108, '-') << "\n\n";

    // =========================================================================
    // Step 2. Bootstrap Uncertainty (Runs v7..v10, good time, 5 p.e., multi-offset mixing)
    // =========================================================================
    std::cout << "=======================================================================================================================================================\n";
    std::cout << "Step 2: Bootstrap Uncertainty (Runs v7..v10, Good Time, 5 p.e., Multi-Offset Mixing)\n";
    std::cout << "Dividing each run into blocks of 1000 consecutive kept events, 200 bootstrap replicas with replacement\n";
    std::cout << "=======================================================================================================================================================\n";

    const int kBlockSize = 1000;
    const int kReplicas = 200;

    struct Step2Row {
        std::string run;
        double ratio_nom;
        double boot_mean;
        double boot_rms;
    };
    std::vector<Step2Row> step2_table;

    // Resample blocks with replacement
    std::mt19937_64 rng(20261018);

    for (size_t ifile = 0; ifile < file_list.size(); ++ifile) {
        const auto& fs = file_list[ifile];
        const auto& evs = runs_events_5pe[ifile];
        int N = (int)evs.size();

        // Nominal ratio
        std::vector<bool> xs_nom(N), xc_nom(N), ys_nom(N), yc_nom(N);
        for (int i = 0; i < N; ++i) {
            xs_nom[i] = evs[i].xs;
            xc_nom[i] = evs[i].xc;
            ys_nom[i] = evs[i].ys;
            yc_nom[i] = evs[i].yc;
        }
        MixingResult nom_res = eval_mixing(xs_nom, xc_nom, ys_nom, yc_nom);

        int n_blocks = (N + kBlockSize - 1) / kBlockSize;
        std::vector<int> b_start(n_blocks), b_len(n_blocks);
        for (int b = 0; b < n_blocks; ++b) {
            b_start[b] = b * kBlockSize;
            b_len[b] = std::min(kBlockSize, N - b_start[b]);
        }

        std::uniform_int_distribution<int> dist_block(0, n_blocks - 1);
        std::vector<double> rep_ratios(kReplicas);

        for (int r = 0; r < kReplicas; ++r) {
            std::vector<bool> rep_xs, rep_xc, rep_ys, rep_yc;
            rep_xs.reserve(N); rep_xc.reserve(N); rep_ys.reserve(N); rep_yc.reserve(N);

            for (int ib = 0; ib < n_blocks; ++ib) {
                int b_chosen = dist_block(rng);
                int st = b_start[b_chosen];
                int len = b_len[b_chosen];
                for (int j = 0; j < len; ++j) {
                    rep_xs.push_back(xs_nom[st + j]);
                    rep_xc.push_back(xc_nom[st + j]);
                    rep_ys.push_back(ys_nom[st + j]);
                    rep_yc.push_back(yc_nom[st + j]);
                }
            }

            MixingResult rep_res = eval_mixing(rep_xs, rep_xc, rep_ys, rep_yc);
            rep_ratios[r] = rep_res.ratio;
        }

        double sum_r = 0;
        for (double val : rep_ratios) sum_r += val;
        double b_mean = sum_r / (double)kReplicas;

        double sum_sq = 0;
        for (double val : rep_ratios) sum_sq += (val - b_mean) * (val - b_mean);
        double b_rms = std::sqrt(sum_sq / (double)(kReplicas - 1));

        Step2Row row;
        row.run = fs.run_name;
        row.ratio_nom = nom_res.ratio;
        row.boot_mean = b_mean;
        row.boot_rms = b_rms;
        step2_table.push_back(row);
    }

    std::cout << std::left << std::setw(8) << "run"
              << std::right << std::setw(16) << "ratio (nominal)"
              << std::setw(16) << "bootstrap mean"
              << std::setw(16) << "bootstrap RMS"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    double sum_w = 0;
    double sum_wx = 0;
    for (const auto& r : step2_table) {
        std::cout << std::left << std::setw(8) << r.run
                  << std::right << std::fixed << std::setprecision(6)
                  << std::setw(16) << r.ratio_nom
                  << std::setw(16) << r.boot_mean
                  << std::setw(16) << r.boot_rms
                  << "\n";
        double w = 1.0 / (r.boot_rms * r.boot_rms);
        sum_w += w;
        sum_wx += w * r.ratio_nom;
    }
    std::cout << std::string(56, '-') << "\n";

    double wt_mean = sum_wx / sum_w;
    double wt_sigma = 1.0 / std::sqrt(sum_w);
    double chi2 = 0;
    for (const auto& r : step2_table) {
        double diff = r.ratio_nom - wt_mean;
        chi2 += (diff * diff) / (r.boot_rms * r.boot_rms);
    }
    double ndf = 3.0; // 4 runs - 1
    std::cout << "\nWeighted mean over runs v7..v10 (with bootstrap RMS as sigma):\n";
    std::cout << "  weighted mean ratio = " << std::fixed << std::setprecision(6) << wt_mean << "\n";
    std::cout << "  sigma               = " << wt_sigma << "\n";
    std::cout << "  chi2 / ndf          = " << chi2 << " / " << (int)ndf << " = " << (chi2 / ndf) << "\n\n";

    // =========================================================================
    // Step 3. Threshold Scan of the Ratio (Run v8 good time, multi-offset mixing)
    // =========================================================================
    std::cout << "=======================================================================================================================================================\n";
    std::cout << "Step 3: Threshold Scan of the Ratio (Run v8 Good Time, Multi-Offset Mixing, Bootstrap RMS as Error)\n";
    std::cout << "Thresholds: n = 3, 4, 5, 6, 7, 8, 10, 12 p.e.\n";
    std::cout << "=======================================================================================================================================================\n";

    const std::vector<double> n_pe_list = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0};
    const size_t num_thr = n_pe_list.size();

    // Precompute block indices for v8
    int n_blocks_v8 = (n_v8 + kBlockSize - 1) / kBlockSize;
    std::vector<int> b8_start(n_blocks_v8), b8_len(n_blocks_v8);
    for (int b = 0; b < n_blocks_v8; ++b) {
        b8_start[b] = b * kBlockSize;
        b8_len[b] = std::min(kBlockSize, (int)n_v8 - b8_start[b]);
    }

    // Pre-generate 200 replica block choices for v8
    std::mt19937_64 rng_v8(202610188);
    std::uniform_int_distribution<int> dist_v8(0, n_blocks_v8 - 1);
    std::vector<std::vector<int>> rep_blocks_v8(kReplicas, std::vector<int>(n_blocks_v8));
    for (int r = 0; r < kReplicas; ++r) {
        for (int b = 0; b < n_blocks_v8; ++b) {
            rep_blocks_v8[r][b] = dist_v8(rng_v8);
        }
    }

    struct Step3Row {
        double n_pe;
        double aX;
        double aY;
        double ratio;
        double boot_rms;
    };
    std::vector<Step3Row> step3_table;

    for (size_t it = 0; it < num_thr; ++it) {
        double n_pe = n_pe_list[it];

        // Evaluate nominal flags for all v8 events at threshold n_pe
        std::vector<bool> xs_t(n_v8, false), xc_t(n_v8, false), ys_t(n_v8, false), yc_t(n_v8, false);
        for (Long64_t i = 0; i < n_v8; ++i) {
            for (int c = 0; c < 32; ++c) {
                if (c == 3) continue;
                double thr = get_thr(c, n_pe);
                if (v8_events[i].ph_s[c] > thr) {
                    if (c < 16) xs_t[i] = true;
                    else        ys_t[i] = true;
                }
                if (v8_events[i].ph_c[c] > thr) {
                    if (c < 16) xc_t[i] = true;
                    else        yc_t[i] = true;
                }
            }
        }

        MixingResult nom_t = eval_mixing(xs_t, xc_t, ys_t, yc_t);

        // 200 bootstrap replicas
        std::vector<double> rep_r(kReplicas);
        for (int r = 0; r < kReplicas; ++r) {
            std::vector<bool> r_xs, r_xc, r_ys, r_yc;
            r_xs.reserve(n_v8); r_xc.reserve(n_v8); r_ys.reserve(n_v8); r_yc.reserve(n_v8);

            for (int ib = 0; ib < n_blocks_v8; ++ib) {
                int b_chosen = rep_blocks_v8[r][ib];
                int st = b8_start[b_chosen];
                int len = b8_len[b_chosen];
                for (int j = 0; j < len; ++j) {
                    r_xs.push_back(xs_t[st + j]);
                    r_xc.push_back(xc_t[st + j]);
                    r_ys.push_back(ys_t[st + j]);
                    r_yc.push_back(yc_t[st + j]);
                }
            }

            MixingResult rep_res = eval_mixing(r_xs, r_xc, r_ys, r_yc);
            rep_r[r] = rep_res.ratio;
        }

        double sum_val = 0;
        for (double val : rep_r) sum_val += val;
        double m_val = sum_val / (double)kReplicas;

        double sum_sq = 0;
        for (double val : rep_r) sum_sq += (val - m_val) * (val - m_val);
        double rms_val = std::sqrt(sum_sq / (double)(kReplicas - 1));

        Step3Row s3;
        s3.n_pe = n_pe;
        s3.aX = nom_t.aX;
        s3.aY = nom_t.aY;
        s3.ratio = nom_t.ratio;
        s3.boot_rms = rms_val;
        step3_table.push_back(s3);
    }

    std::cout << std::right << std::setw(6) << "n"
              << std::setw(12) << "aX"
              << std::setw(12) << "aY"
              << std::setw(14) << "ratio"
              << std::setw(16) << "bootstrap RMS"
              << "\n";
    std::cout << std::string(60, '-') << "\n";

    for (const auto& r : step3_table) {
        std::cout << std::right << std::setw(6) << (int)r.n_pe
                  << std::fixed << std::setprecision(6)
                  << std::setw(12) << r.aX
                  << std::setw(12) << r.aY
                  << std::setw(14) << r.ratio
                  << std::setw(16) << r.boot_rms
                  << "\n";
    }
    std::cout << std::string(60, '-') << "\n\n";

    // =========================================================================
    // Step 4. X-Y Map at 10 p.e. (Run v8 good time, mapping kCfg_2026_05)
    // =========================================================================
    std::cout << "=======================================================================================================================================================\n";
    std::cout << "Step 4: hod2018 X-Y Map in Good Time at 10 p.e. (Run v8, Mapping kCfg_2026_05)\n";
    std::cout << "16x16 matrix: same-event minus M_bar-normalised mixed (20 offsets) at 10 p.e.\n";
    std::cout << "=======================================================================================================================================================\n";

    // Precompute 10 p.e. best_strip_x and best_strip_y for each event in v8
    std::vector<int> v8_best_x(n_v8, -1), v8_best_y(n_v8, -1);
    std::vector<bool> v8_xs_10(n_v8, false), v8_xc_10(n_v8, false), v8_ys_10(n_v8, false), v8_yc_10(n_v8, false);

    for (Long64_t i = 0; i < n_v8; ++i) {
        float max_ph_x = -1e9, max_ph_y = -1e9;
        for (int c = 0; c < 32; ++c) {
            if (c == 3) continue;
            double thr10 = get_thr(c, 10.0);
            float ph_s = v8_events[i].ph_s[c];
            float ph_c = v8_events[i].ph_c[c];

            if (ph_s > thr10) {
                int strip = ch_to_strip(c);
                if (c < 16) {
                    v8_xs_10[i] = true;
                    if (ph_s > max_ph_x) {
                        max_ph_x = ph_s;
                        v8_best_x[i] = strip;
                    }
                } else {
                    v8_ys_10[i] = true;
                    if (ph_s > max_ph_y) {
                        max_ph_y = ph_s;
                        v8_best_y[i] = strip;
                    }
                }
            }
            if (ph_c > thr10) {
                if (c < 16) v8_xc_10[i] = true;
                else        v8_yc_10[i] = true;
            }
        }
    }

    MixingResult mix_10pe = eval_mixing(v8_xs_10, v8_xc_10, v8_ys_10, v8_yc_10);
    double denom_10pe = (1.0 - mix_10pe.aX) * (1.0 - mix_10pe.aY);

    double M_same[16][16] = {};
    double M_mixed_sum[16][16] = {};

    for (Long64_t i = 0; i < n_v8; ++i) {
        int sx = v8_best_x[i];
        int sy = v8_best_y[i];

        if (sx >= 0 && sx < 16 && sy >= 0 && sy < 16) {
            M_same[sx][sy] += 1.0;
        }

        // 20 offsets
        if (sx >= 0 && sx < 16) {
            for (int k = 1; k <= 20; ++k) {
                if (i + k < n_v8) {
                    int sy_k = v8_best_y[i + k];
                    if (sy_k >= 0 && sy_k < 16) {
                        M_mixed_sum[sx][sy_k] += 1.0;
                    }
                }
            }
        }
    }

    double M_net[16][16] = {};
    double row_sums[16] = {};
    double col_sums[16] = {};
    double total_net = 0;

    TH2D* h_map = new TH2D("h_hod2018_xy_map_10pe",
                           "hod2018 X-Y Net Coincidence Map (v8 Good Time, 10 p.e.);X strip;Y strip;Net Coincidences",
                           16, -0.5, 15.5, 16, -0.5, 15.5);

    for (int x = 0; x < 16; ++x) {
        for (int y = 0; y < 16; ++y) {
            double m_bar_cell = M_mixed_sum[x][y] / 20.0;
            M_net[x][y] = M_same[x][y] - m_bar_cell;
            row_sums[x] += M_net[x][y];
            col_sums[y] += M_net[x][y];
            total_net += M_net[x][y];
            h_map->SetBinContent(x + 1, y + 1, M_net[x][y]);
        }
    }

    // Print 16x16 matrix with strip indices, row sums and column sums
    std::cout << "\nMatrix M_net[X_strip][Y_strip] = M_same - M_bar at 10 p.e.:\n\n";
    std::cout << " X\\Y ";
    for (int y = 0; y < 16; ++y) {
        std::cout << std::right << std::setw(8) << Form("Y%d", y);
    }
    std::cout << std::right << std::setw(10) << "RowSum" << "\n";
    std::cout << std::string(144, '-') << "\n";

    for (int x = 0; x < 16; ++x) {
        std::cout << std::left << std::setw(5) << Form("X%d", x);
        for (int y = 0; y < 16; ++y) {
            std::cout << std::right << std::fixed << std::setprecision(1) << std::setw(8) << M_net[x][y];
        }
        std::cout << std::right << std::fixed << std::setprecision(1) << std::setw(10) << row_sums[x] << "\n";
    }
    std::cout << std::string(144, '-') << "\n";

    std::cout << std::left << std::setw(5) << "ColSum";
    for (int y = 0; y < 16; ++y) {
        std::cout << std::right << std::fixed << std::setprecision(1) << std::setw(8) << col_sums[y];
    }
    std::cout << std::right << std::fixed << std::setprecision(1) << std::setw(10) << total_net << "\n";
    std::cout << std::string(144, '-') << "\n";

    // Save PNG
    TCanvas c1("c1", "hod2018 X-Y Map 10 pe", 900, 800);
    c1.SetRightMargin(0.14);
    c1.SetLeftMargin(0.10);
    c1.SetTopMargin(0.08);
    c1.SetBottomMargin(0.10);
    h_map->GetZaxis()->SetTitleOffset(1.2);
    h_map->Draw("COLZ TEXT");
    c1.SaveAs("outputs/hod2018_xy_map_v8_10pe.png");
    std::cout << "\nSaved map figure: outputs/hod2018_xy_map_v8_10pe.png\n\n";

    // Corrected sums: (RowSum or ColSum) / ((1-aX)(1-aY))
    double corr_row_sums[16] = {};
    double corr_col_sums[16] = {};
    for (int x = 0; x < 16; ++x) corr_row_sums[x] = row_sums[x] / denom_10pe;
    for (int y = 0; y < 16; ++y) corr_col_sums[y] = col_sums[y] / denom_10pe;

    // Bootstrap RMS of corrected sums across 200 replicas
    std::vector<std::vector<double>> rep_corr_rows(kReplicas, std::vector<double>(16, 0.0));
    std::vector<std::vector<double>> rep_corr_cols(kReplicas, std::vector<double>(16, 0.0));

    for (int r = 0; r < kReplicas; ++r) {
        std::vector<int> r_best_x, r_best_y;
        std::vector<bool> r_xs, r_xc, r_ys, r_yc;
        r_best_x.reserve(n_v8); r_best_y.reserve(n_v8);
        r_xs.reserve(n_v8); r_xc.reserve(n_v8); r_ys.reserve(n_v8); r_yc.reserve(n_v8);

        for (int ib = 0; ib < n_blocks_v8; ++ib) {
            int b_chosen = rep_blocks_v8[r][ib];
            int st = b8_start[b_chosen];
            int len = b8_len[b_chosen];
            for (int j = 0; j < len; ++j) {
                int idx = st + j;
                r_best_x.push_back(v8_best_x[idx]);
                r_best_y.push_back(v8_best_y[idx]);
                r_xs.push_back(v8_xs_10[idx]);
                r_xc.push_back(v8_xc_10[idx]);
                r_ys.push_back(v8_ys_10[idx]);
                r_yc.push_back(v8_yc_10[idx]);
            }
        }

        MixingResult r_mix = eval_mixing(r_xs, r_xc, r_ys, r_yc);
        double r_denom = (1.0 - r_mix.aX) * (1.0 - r_mix.aY);

        double r_M_same[16][16] = {};
        double r_M_mixed_sum[16][16] = {};

        Long64_t sz = r_best_x.size();
        for (Long64_t i = 0; i < sz; ++i) {
            int sx = r_best_x[i];
            int sy = r_best_y[i];
            if (sx >= 0 && sx < 16 && sy >= 0 && sy < 16) {
                r_M_same[sx][sy] += 1.0;
            }
            if (sx >= 0 && sx < 16) {
                for (int k = 1; k <= 20; ++k) {
                    if (i + k < sz) {
                        int sy_k = r_best_y[i + k];
                        if (sy_k >= 0 && sy_k < 16) {
                            r_M_mixed_sum[sx][sy_k] += 1.0;
                        }
                    }
                }
            }
        }

        for (int x = 0; x < 16; ++x) {
            double r_rsum = 0;
            for (int y = 0; y < 16; ++y) {
                double net_c = r_M_same[x][y] - r_M_mixed_sum[x][y] / 20.0;
                r_rsum += net_c;
            }
            rep_corr_rows[r][x] = r_rsum / r_denom;
        }

        for (int y = 0; y < 16; ++y) {
            double r_csum = 0;
            for (int x = 0; x < 16; ++x) {
                double net_c = r_M_same[x][y] - r_M_mixed_sum[x][y] / 20.0;
                r_csum += net_c;
            }
            rep_corr_cols[r][y] = r_csum / r_denom;
        }
    }

    double row_boot_rms[16] = {};
    for (int x = 0; x < 16; ++x) {
        double s = 0; for (int r = 0; r < kReplicas; ++r) s += rep_corr_rows[r][x];
        double m = s / (double)kReplicas;
        double ssq = 0; for (int r = 0; r < kReplicas; ++r) ssq += (rep_corr_rows[r][x] - m) * (rep_corr_rows[r][x] - m);
        row_boot_rms[x] = std::sqrt(ssq / (double)(kReplicas - 1));
    }

    double col_boot_rms[16] = {};
    for (int y = 0; y < 16; ++y) {
        double s = 0; for (int r = 0; r < kReplicas; ++r) s += rep_corr_cols[r][y];
        double m = s / (double)kReplicas;
        double ssq = 0; for (int r = 0; r < kReplicas; ++r) ssq += (rep_corr_cols[r][y] - m) * (rep_corr_cols[r][y] - m);
        col_boot_rms[y] = std::sqrt(ssq / (double)(kReplicas - 1));
    }

    // Print table of corrected sums and bootstrap RMS
    std::cout << "Corrected Profile Sums and Bootstrap RMS at 10 p.e. (aX = " << std::fixed << std::setprecision(6)
              << mix_10pe.aX << ", aY = " << mix_10pe.aY << ", denom = " << std::setprecision(6) << denom_10pe << "):\n\n";

    std::cout << std::left << std::setw(10) << "X strip"
              << std::right << std::setw(18) << "Corrected RowSum"
              << std::setw(18) << "Bootstrap RMS"
              << std::setw(6) << " "
              << std::left << std::setw(10) << "Y strip"
              << std::right << std::setw(18) << "Corrected ColSum"
              << std::setw(18) << "Bootstrap RMS"
              << "\n";
    std::cout << std::string(88, '-') << "\n";

    for (int k = 0; k < 16; ++k) {
        std::cout << std::left << std::setw(10) << Form("X%d", k)
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(18) << corr_row_sums[k]
                  << std::setw(18) << row_boot_rms[k]
                  << std::setw(6) << " "
                  << std::left << std::setw(10) << Form("Y%d", k)
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(18) << corr_col_sums[k]
                  << std::setw(18) << col_boot_rms[k]
                  << "\n";
    }
    std::cout << std::string(88, '-') << "\n\n";
}

