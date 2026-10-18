#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
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

#include "../src/hodo_mapping.h"

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::string channels;
    std::string reason;
};

// Hit info for an event
struct EventHits {
    Long_t timestamp;
    std::string dh_min;
    // For nominal window [4128, 4272] ns:
    bool xs[3]; // thr 4, 5, 6 pe
    bool xc[3];
    bool ys[3];
    bool yc[3];
    // For narrow window [4144, 4256] ns at 5 pe:
    bool xs_narrow;
    bool xc_narrow;
    bool ys_narrow;
    bool yc_narrow;
    // For 2D map (nominal window, 5 pe):
    int best_strip_x; // -1 if no hit
    int best_strip_y; // -1 if no hit
    double max_ph_x;
    double max_ph_y;
};

struct Step1Result {
    std::string run;
    long long n_used;
    long long N_XY_same;
    double M_bar;
    long long N_XY_mixed_1;
    double aX;
    double aY;
    double rX;
    double rY;
    double N_XY_real;
    double sigma;
    double ratio;
    double sigma_ratio;
};

void hod2018_analysis_step1_3() {
    gStyle->SetOptStat(0);

    // 1. Calibration parameters from outputs/pe_calibration_v8.txt
    double p0_cal[32] = {};
    double G_cal[32] = {};
    double T_pe[32][3] = {}; // 0: 4 pe, 1: 5 pe, 2: 6 pe
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
                T_pe[ch][0] = p0 + (4.0 - 1.5) * G; // 4 pe
                T_pe[ch][1] = p0 + (5.0 - 1.5) * G; // 5 pe
                T_pe[ch][2] = p0 + (6.0 - 1.5) * G; // 6 pe
            }
        }
    }
    pf.close();
    for (int it = 0; it < 3; ++it) T_pe[3][it] = 20.0; // dead channel

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

    // Nominal windows
    constexpr double kNomSigTlo = 4128.0, kNomSigThi = 4272.0;
    constexpr double kNomCtlTlo = 6128.0, kNomCtlThi = 6272.0;

    // Narrow windows
    constexpr double kNarSigTlo = 4144.0, kNarSigThi = 4256.0;
    constexpr double kNarCtlTlo = 6144.0, kNarCtlThi = 6256.0;

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

    std::vector<Step1Result> step1_results;
    std::vector<EventHits> v8_good_events; // save for step 2 and 3

    for (const auto& fs : file_list) {
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

        std::vector<EventHits> kept;

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

            EventHits eh{};
            eh.timestamp = timestamp;
            eh.dh_min = s_min;
            eh.best_strip_x = -1;
            eh.best_strip_y = -1;
            eh.max_ph_x = -1e9;
            eh.max_ph_y = -1e9;

            for (int c = 0; c < 32; ++c) {
                if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
                const auto& wf = *cw[c];

                // 1. Nominal signal window [4128, 4272] ns
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

                // 2. Nominal control window [6128, 6272] ns
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

                // Checks for 4, 5, 6 pe thresholds (indices 0, 1, 2)
                for (int it = 0; it < 3; ++it) {
                    if (ph_s > T_pe[c][it]) {
                        if (c < 16) eh.xs[it] = true;
                        else        eh.ys[it] = true;
                    }
                    if (ph_c > T_pe[c][it]) {
                        if (c < 16) eh.xc[it] = true;
                        else        eh.yc[it] = true;
                    }
                }

                // Check for highest-ph hit at 5 pe (it = 1) in nominal signal window
                if (ph_s > T_pe[c][1]) {
                    int strip = ch_to_strip(c);
                    if (c < 16) {
                        if (ph_s > eh.max_ph_x) {
                            eh.max_ph_x = ph_s;
                            eh.best_strip_x = strip;
                        }
                    } else {
                        if (ph_s > eh.max_ph_y) {
                            eh.max_ph_y = ph_s;
                            eh.best_strip_y = strip;
                        }
                    }
                }

                // 3. Narrow window [4144, 4256] ns at 5 pe
                if (fs.run_name == "v8") {
                    double vmin_sn = 1e30; int imin_sn = -1;
                    for (int i = 0; i < n_samp; ++i) {
                        if (tv[i] >= kNarSigTlo && tv[i] <= kNarSigThi && wf[i] < vmin_sn) {
                            vmin_sn = wf[i]; imin_sn = i;
                        }
                    }
                    if (imin_sn >= 0) {
                        double ped = 0;
                        if (imin_sn >= 12 && imin_sn - 4 < n_samp) {
                            double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_sn + j];
                            ped = s / 9.0;
                        } else {
                            double s = 0; int n = 0;
                            for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                            ped = (n > 0) ? s / n : 0;
                        }
                        if (ped - vmin_sn > T_pe[c][1]) {
                            if (c < 16) eh.xs_narrow = true;
                            else        eh.ys_narrow = true;
                        }
                    }

                    // Narrow control [6144, 6256] ns
                    double vmin_cn = 1e30; int imin_cn = -1;
                    for (int i = 0; i < n_samp; ++i) {
                        if (tv[i] >= kNarCtlTlo && tv[i] <= kNarCtlThi && wf[i] < vmin_cn) {
                            vmin_cn = wf[i]; imin_cn = i;
                        }
                    }
                    if (imin_cn >= 0) {
                        double ped = 0;
                        if (imin_cn >= 12 && imin_cn - 4 < n_samp) {
                            double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_cn + j];
                            ped = s / 9.0;
                        } else {
                            double s = 0; int n = 0;
                            for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                            ped = (n > 0) ? s / n : 0;
                        }
                        if (ped - vmin_cn > T_pe[c][1]) {
                            if (c < 16) eh.xc_narrow = true;
                            else        eh.yc_narrow = true;
                        }
                    }
                }
            }

            kept.push_back(eh);
        }

        if (fs.run_name == "v8") {
            v8_good_events = kept;
        }

        // Compute multi-offset mixing for Step 1
        Long64_t n_used = (Long64_t)kept.size();
        Long64_t N_pairs = (n_used > 1) ? (n_used - 1) : 0;
        double N_d = (double)N_pairs;

        long long N_X_S = 0, N_X_C = 0, N_Y_S = 0, N_Y_C = 0;
        long long N_XY_same = 0;
        for (Long64_t i = 0; i < N_pairs; ++i) {
            if (kept[i].xs[1]) N_X_S++;
            if (kept[i].xc[1]) N_X_C++;
            if (kept[i].ys[1]) N_Y_S++;
            if (kept[i].yc[1]) N_Y_C++;
            if (kept[i].xs[1] && kept[i].ys[1]) N_XY_same++;
        }

        double aX = (N_d > 0) ? ((double)N_X_C / N_d) : 0.0;
        double aY = (N_d > 0) ? ((double)N_Y_C / N_d) : 0.0;
        double rX = (1.0 - aX > 1e-12) ? (((double)N_X_S / N_d - aX) / (1.0 - aX)) : 0.0;
        double rY = (1.0 - aY > 1e-12) ? (((double)N_Y_S / N_d - aY) / (1.0 - aY)) : 0.0;

        // Multi-offset mixing k = 1..20
        double sum_m = 0;
        long long m_k[21] = {};
        for (int k = 1; k <= 20; ++k) {
            long long m = 0;
            for (Long64_t i = 0; i + k < n_used; ++i) {
                if (kept[i].xs[1] && kept[i + k].ys[1]) m++;
            }
            m_k[k] = m;
            sum_m += (double)m;
        }
        double M_bar = sum_m / 20.0;
        long long N_XY_mixed_1 = m_k[1];

        double denom = (1.0 - aX) * (1.0 - aY);
        double N_XY_real = (std::abs(denom) > 1e-12) ? ((N_XY_same - M_bar) / denom + rX * rY * N_d) : 0.0;
        double sigma = (std::abs(denom) > 1e-12) ? (std::sqrt(N_XY_same + M_bar / 20.0) / denom) : 0.0;

        double ratio = (n_used > 0) ? (N_XY_real / (0.0360 * (double)n_used)) : 0.0;
        double sigma_ratio = (n_used > 0) ? (sigma / (0.0360 * (double)n_used)) : 0.0;

        Step1Result s1;
        s1.run = fs.run_name;
        s1.n_used = n_used;
        s1.N_XY_same = N_XY_same;
        s1.M_bar = M_bar;
        s1.N_XY_mixed_1 = N_XY_mixed_1;
        s1.aX = aX;
        s1.aY = aY;
        s1.rX = rX;
        s1.rY = rY;
        s1.N_XY_real = N_XY_real;
        s1.sigma = sigma;
        s1.ratio = ratio;
        s1.sigma_ratio = sigma_ratio;
        step1_results.push_back(s1);
    }

    // =========================================================================
    // Print Step 1 Table
    // =========================================================================
    std::cout << "\n=======================================================================================================================================================\n";
    std::cout << "Step 1: Statistical Uncertainty with Multi-Offset Mixing (Good Time Only, Runs v7..v10, 5 p.e.)\n";
    std::cout << "=======================================================================================================================================================\n";
    std::cout << std::left << std::setw(8) << "run"
              << std::right << std::setw(15) << "N_events_used"
              << std::setw(12) << "N_XY_same"
              << std::setw(12) << "M_bar"
              << std::setw(16) << "N_XY_mixed(1)"
              << std::setw(10) << "aX"
              << std::setw(10) << "aY"
              << std::setw(10) << "rX"
              << std::setw(10) << "rY"
              << std::setw(14) << "N_XY_real"
              << std::setw(12) << "sigma"
              << std::setw(12) << "ratio"
              << std::setw(14) << "sigma_ratio"
              << "\n";
    std::cout << std::string(155, '-') << "\n";

    double sum_w = 0;
    double sum_wx = 0;
    for (const auto& r : step1_results) {
        std::cout << std::left << std::setw(8) << r.run
                  << std::right << std::setw(15) << r.n_used
                  << std::setw(12) << r.N_XY_same
                  << std::fixed << std::setprecision(2)
                  << std::setw(12) << r.M_bar
                  << std::setw(16) << r.N_XY_mixed_1
                  << std::setprecision(6)
                  << std::setw(10) << r.aX
                  << std::setw(10) << r.aY
                  << std::setw(10) << r.rX
                  << std::setw(10) << r.rY
                  << std::setprecision(2)
                  << std::setw(14) << r.N_XY_real
                  << std::setw(12) << r.sigma
                  << std::setprecision(6)
                  << std::setw(12) << r.ratio
                  << std::setw(14) << r.sigma_ratio
                  << "\n";

        if (r.sigma_ratio > 0) {
            double w = 1.0 / (r.sigma_ratio * r.sigma_ratio);
            sum_w += w;
            sum_wx += w * r.ratio;
        }
    }
    std::cout << std::string(155, '-') << "\n";

    double mean_ratio = (sum_w > 0) ? (sum_wx / sum_w) : 0.0;
    double sigma_mean = (sum_w > 0) ? (1.0 / std::sqrt(sum_w)) : 0.0;

    double chi2 = 0;
    for (const auto& r : step1_results) {
        if (r.sigma_ratio > 0) {
            double d = (r.ratio - mean_ratio) / r.sigma_ratio;
            chi2 += d * d;
        }
    }
    int ndf = (int)step1_results.size() - 1;
    double chi2_ndf = (ndf > 0) ? (chi2 / (double)ndf) : 0.0;

    std::cout << "\nInverse-variance weighted mean over runs v7..v10:\n";
    std::cout << "  weighted mean ratio = " << std::fixed << std::setprecision(6) << mean_ratio << "\n";
    std::cout << "  sigma               = " << sigma_mean << "\n";
    std::cout << "  chi2 / ndf          = " << std::setprecision(4) << chi2 << " / " << ndf
              << " = " << chi2_ndf << "\n";

    // =========================================================================
    // Step 2: Systematic Variations (run v8 good time, multi-offset mixing)
    // =========================================================================
    std::cout << "\n=======================================================================================================================================================\n";
    std::cout << "Step 2: Systematic Variations (Run v8 Good Time, Multi-Offset Mixing)\n";
    std::cout << "=======================================================================================================================================================\n";
    std::cout << std::left << std::setw(42) << "Variation"
              << std::right << std::setw(14) << "N_XY_real"
              << std::setw(12) << "sigma"
              << std::setw(12) << "ratio"
              << std::setw(14) << "sigma_ratio"
              << std::setw(30) << "Notes"
              << "\n";
    std::cout << std::string(112, '-') << "\n";

    Long64_t n_v8 = (Long64_t)v8_good_events.size();
    Long64_t N_pairs_v8 = (n_v8 > 1) ? (n_v8 - 1) : 0;
    double N_d_v8 = (double)N_pairs_v8;

    auto eval_variation = [&](const std::string& name,
                              const std::vector<bool>& xs_vec,
                              const std::vector<bool>& xc_vec,
                              const std::vector<bool>& ys_vec,
                              const std::vector<bool>& yc_vec,
                              const std::string& note = "") {
        long long N_X_S = 0, N_X_C = 0, N_Y_S = 0, N_Y_C = 0;
        long long N_XY_same = 0;
        for (Long64_t i = 0; i < N_pairs_v8; ++i) {
            if (xs_vec[i]) N_X_S++;
            if (xc_vec[i]) N_X_C++;
            if (ys_vec[i]) N_Y_S++;
            if (yc_vec[i]) N_Y_C++;
            if (xs_vec[i] && ys_vec[i]) N_XY_same++;
        }

        double aX = (N_d_v8 > 0) ? ((double)N_X_C / N_d_v8) : 0.0;
        double aY = (N_d_v8 > 0) ? ((double)N_Y_C / N_d_v8) : 0.0;
        double rX = (1.0 - aX > 1e-12) ? (((double)N_X_S / N_d_v8 - aX) / (1.0 - aX)) : 0.0;
        double rY = (1.0 - aY > 1e-12) ? (((double)N_Y_S / N_d_v8 - aY) / (1.0 - aY)) : 0.0;

        double sum_m = 0;
        for (int k = 1; k <= 20; ++k) {
            long long m = 0;
            for (Long64_t i = 0; i + k < n_v8; ++i) {
                if (xs_vec[i] && ys_vec[i + k]) m++;
            }
            sum_m += (double)m;
        }
        double M_bar = sum_m / 20.0;

        double denom = (1.0 - aX) * (1.0 - aY);
        double N_XY_real = (std::abs(denom) > 1e-12) ? ((N_XY_same - M_bar) / denom + rX * rY * N_d_v8) : 0.0;
        double sigma = (std::abs(denom) > 1e-12) ? (std::sqrt(N_XY_same + M_bar / 20.0) / denom) : 0.0;
        double ratio = (n_v8 > 0) ? (N_XY_real / (0.0360 * (double)n_v8)) : 0.0;
        double sigma_ratio = (n_v8 > 0) ? (sigma / (0.0360 * (double)n_v8)) : 0.0;

        std::cout << std::left << std::setw(42) << name
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(14) << N_XY_real
                  << std::setw(12) << sigma
                  << std::setprecision(6)
                  << std::setw(12) << ratio
                  << std::setw(14) << sigma_ratio
                  << std::setw(30) << note
                  << "\n";
    };

    // Prepare vectors
    std::vector<bool> xs_4(n_v8), xc_4(n_v8), ys_4(n_v8), yc_4(n_v8);
    std::vector<bool> xs_5(n_v8), xc_5(n_v8), ys_5(n_v8), yc_5(n_v8);
    std::vector<bool> xs_6(n_v8), xc_6(n_v8), ys_6(n_v8), yc_6(n_v8);
    std::vector<bool> xs_nar(n_v8), xc_nar(n_v8), ys_nar(n_v8), yc_nar(n_v8);

    long long n_xs_nom = 0, n_xc_nom = 0, n_ys_nom = 0, n_yc_nom = 0, n_same_nom = 0;

    for (Long64_t i = 0; i < n_v8; ++i) {
        xs_4[i] = v8_good_events[i].xs[0]; xc_4[i] = v8_good_events[i].xc[0];
        ys_4[i] = v8_good_events[i].ys[0]; yc_4[i] = v8_good_events[i].yc[0];

        xs_5[i] = v8_good_events[i].xs[1]; xc_5[i] = v8_good_events[i].xc[1];
        ys_5[i] = v8_good_events[i].ys[1]; yc_5[i] = v8_good_events[i].yc[1];

        xs_6[i] = v8_good_events[i].xs[2]; xc_6[i] = v8_good_events[i].xc[2];
        ys_6[i] = v8_good_events[i].ys[2]; yc_6[i] = v8_good_events[i].yc[2];

        xs_nar[i] = v8_good_events[i].xs_narrow; xc_nar[i] = v8_good_events[i].xc_narrow;
        ys_nar[i] = v8_good_events[i].ys_narrow; yc_nar[i] = v8_good_events[i].yc_narrow;

        if (i < N_pairs_v8) {
            if (xs_5[i]) n_xs_nom++;
            if (xc_5[i]) n_xc_nom++;
            if (ys_5[i]) n_ys_nom++;
            if (yc_5[i]) n_yc_nom++;
            if (xs_5[i] && ys_5[i]) n_same_nom++;
        }
    }

    // (a) Threshold variations
    eval_variation("(a) threshold 4 p.e. (nominal window)", xs_4, xc_4, ys_4, yc_4, "");
    eval_variation("(a) threshold 5 p.e. (nominal, reference)", xs_5, xc_5, ys_5, yc_5, "nominal selection");
    eval_variation("(a) threshold 6 p.e. (nominal window)", xs_6, xc_6, ys_6, yc_6, "");

    // (b) Window variations
    eval_variation("(b) window [4128, 4272] ns (nominal 144 ns)", xs_5, xc_5, ys_5, yc_5, "nominal window");
    eval_variation("(b) window [4144, 4256] ns (narrow 112 ns)", xs_nar, xc_nar, ys_nar, yc_nar, "narrow window");

    // (c) Accidental estimate from control window instead of mixing
    double N_tot = (double)N_pairs_v8; // or n_v8
    double P_Xc = (double)n_xc_nom / N_tot;
    double P_Ys = (double)n_ys_nom / N_tot;
    double P_Xs = (double)n_xs_nom / N_tot;
    double P_Yc = (double)n_yc_nom / N_tot;

    double N_XY_acc = N_tot * P_Xc * P_Ys + N_tot * P_Xs * P_Yc - N_tot * P_Xc * P_Yc;
    double N_XY_real_ctl = (double)n_same_nom - N_XY_acc;
    double sigma_ctl = std::sqrt((double)n_same_nom + N_XY_acc);
    double ratio_ctl = N_XY_real_ctl / (0.0360 * (double)n_v8);
    double sigma_ratio_ctl = sigma_ctl / (0.0360 * (double)n_v8);

    std::cout << std::left << std::setw(42) << "(c) accidental from control window (nominal)"
              << std::right << std::fixed << std::setprecision(2)
              << std::setw(14) << N_XY_real_ctl
              << std::setw(12) << sigma_ctl
              << std::setprecision(6)
              << std::setw(12) << ratio_ctl
              << std::setw(14) << sigma_ratio_ctl
              << std::setw(30) << "ctl-window formula"
              << "\n";
    std::cout << std::string(112, '-') << "\n";

    std::cout << "\nComparison for nominal selection (run v8 good time, 5 p.e.):\n";
    std::cout << "  Mixing (M_bar = " << std::fixed << std::setprecision(2) << step1_results[1].M_bar << "): ";
    std::cout << "N_XY_real = " << step1_results[1].N_XY_real
              << " ± " << step1_results[1].sigma
              << ", ratio = " << std::setprecision(6) << step1_results[1].ratio
              << " ± " << step1_results[1].sigma_ratio << "\n";
    std::cout << "  Control window (N_XY_acc = " << std::setprecision(2) << N_XY_acc << "): ";
    std::cout << "N_XY_real_ctl = " << N_XY_real_ctl
              << " ± " << sigma_ctl
              << ", ratio = " << std::setprecision(6) << ratio_ctl
              << " ± " << sigma_ratio_ctl << "\n";

    // =========================================================================
    // Step 3: hod2018 X-Y map in good time (run v8, 5 p.e., mapping kCfg_2026_05)
    // =========================================================================
    std::cout << "\n=======================================================================================================================================================\n";
    std::cout << "Step 3: hod2018 X-Y Map in Good Time (Run v8, 5 p.e., Mapping kCfg_2026_05)\n";
    std::cout << "16x16 matrix: same-event minus M_bar-normalised mixed (20 offsets)\n";
    std::cout << "=======================================================================================================================================================\n";

    double M_same[16][16] = {};
    double M_mixed_sum[16][16] = {};

    for (Long64_t i = 0; i < n_v8; ++i) {
        int sx = v8_good_events[i].best_strip_x;
        int sy = v8_good_events[i].best_strip_y;

        if (sx >= 0 && sx < 16 && sy >= 0 && sy < 16) {
            M_same[sx][sy] += 1.0;
        }

        // 20 offsets
        if (sx >= 0 && sx < 16) {
            for (int k = 1; k <= 20; ++k) {
                if (i + k < n_v8) {
                    int sy_k = v8_good_events[i + k].best_strip_y;
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

    TH2D* h_map = new TH2D("h_hod2018_xy_map",
                           "hod2018 X-Y Net Coincidence Map (v8 Good Time, 5 p.e.);X strip;Y strip;Net Coincidences",
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
    std::cout << "\nMatrix M_net[X_strip][Y_strip] = M_same - M_bar:\n\n";
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
    TCanvas c1("c1", "hod2018 X-Y Map", 900, 800);
    c1.SetRightMargin(0.14);
    c1.SetLeftMargin(0.10);
    c1.SetTopMargin(0.08);
    c1.SetBottomMargin(0.10);
    h_map->GetZaxis()->SetTitleOffset(1.2);
    h_map->Draw("COLZ TEXT");
    c1.SaveAs("outputs/hod2018_xy_map_v8_goodtime.png");
    std::cout << "\nSaved map figure: outputs/hod2018_xy_map_v8_goodtime.png\n";
}
