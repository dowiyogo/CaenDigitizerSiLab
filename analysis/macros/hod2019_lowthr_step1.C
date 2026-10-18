// analysis/macros/hod2019_lowthr_step1.C
// TASK: HOD2019-LOWTHR-1
// Measures hod2019 X-Y coincidences at low thresholds (3, 5, 8, 12, 20, 35 mV)
// on run v8 good time only, using multi-offset mixing (20 offsets).

#include <TFile.h>
#include <TTree.h>
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

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::string channels;
    std::string reason;
};

struct EventData {
    Long_t timestamp;
    std::string dh_min;
    float ph_S[32]; // channels 32..63: ph_S[c - 32]
    float ph_C[32]; // channels 32..63: ph_C[c - 32]
};

void hod2019_lowthr_step1(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root") {
    gROOT->SetBatch(kTRUE);

    // 1. Read good intervals
    std::vector<Interval> intervals;
    std::ifstream fcsv("runs/good_intervals_2026.csv");
    if (!fcsv.is_open()) {
        std::cerr << "Error: runs/good_intervals_2026.csv not found!\n";
        return;
    }
    std::string line;
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

    // Signal window W2 and control window
    constexpr double kSigTlo = 4128.0, kSigThi = 4272.0;
    constexpr double kCtlTlo = 6128.0, kCtlThi = 6272.0;

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
    for (int c = 32; c <= 63; ++c) {
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c - 32]);
    }
    tree->SetCacheSize(256LL << 20);

    std::vector<EventData> kept;
    kept.reserve(120000);

    std::cout << "Reading entries from " << fname << " (total " << nEnt << ")...\n";

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

        // Good time only: exclude EXCLUDE and HV_OFF
        if (is_exclude || is_hv_off) continue;

        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        EventData ed{};
        ed.timestamp = timestamp;
        ed.dh_min = s_min;

        for (int c = 32; c <= 63; ++c) {
            int idx = c - 32;
            if (!cw[idx] || (int)cw[idx]->size() != n_samp) continue;
            const auto& wf = *cw[idx];

            // Global pedestal: samples with t < 2500 ns
            double s_ped = 0;
            int n_ped = 0;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] < 2500.0) {
                    s_ped += wf[i];
                    n_ped++;
                }
            }
            double ped = (n_ped > 0) ? (s_ped / n_ped) : 0.0;

            // Signal window W2 = [4128, 4272] ns
            double min_S = 1e30;
            int imin_S = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < min_S) {
                    min_S = wf[i];
                    imin_S = i;
                }
            }
            ed.ph_S[idx] = (imin_S >= 0) ? (float)(ped - min_S) : 0.0f;

            // Control window = [6128, 6272] ns
            double min_C = 1e30;
            int imin_C = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < min_C) {
                    min_C = wf[i];
                    imin_C = i;
                }
            }
            ed.ph_C[idx] = (imin_C >= 0) ? (float)(ped - min_C) : 0.0f;
        }

        kept.push_back(ed);
    }

    Long64_t n_used = (Long64_t)kept.size();
    std::cout << "Good time events kept: " << n_used << "\n\n";

    Long64_t N_pairs = (n_used > 1) ? (n_used - 1) : 0;
    double N_d = (double)N_pairs;
    double N_tot = (double)n_used;

    const std::vector<double> thresholds = {3.0, 5.0, 8.0, 12.0, 20.0, 35.0};

    struct RowOutput {
        double T;
        long long N_X_S;
        long long N_X_C;
        long long N_Y_S;
        long long N_Y_C;
        long long N_XY_same;
        double M_bar;
        double aX;
        double aY;
        double rX;
        double rY;
        double N_XY_real;
        double sigma;
        double ratio;
        double sigma_ratio;
        long long n_ev_2ch_2grp;
    };

    std::vector<RowOutput> table1_rows;

    for (double T : thresholds) {
        // Precompute per-event flags for threshold T
        // hod2019 Y: ch 32..47 (indices 0..15)
        // hod2019 X: ch 48..63 (indices 16..31)
        // Groups:
        // G4: ch 32..39 (indices 0..7)
        // G5: ch 40..47 (indices 8..15)
        // G6: ch 48..55 (indices 16..23)
        // G7: ch 56..63 (indices 24..31)
        std::vector<bool> xs(n_used, false);
        std::vector<bool> xc(n_used, false);
        std::vector<bool> ys(n_used, false);
        std::vector<bool> yc(n_used, false);
        std::vector<bool> ev_2ch_2grp(n_used, false);

        for (Long64_t i = 0; i < n_used; ++i) {
            const auto& ed = kept[i];
            int n_hit_ch_S = 0;
            bool grp_hit_S[8] = {};

            for (int c = 32; c <= 63; ++c) {
                int idx = c - 32;
                int grp = c / 8;

                if (ed.ph_S[idx] > T) {
                    n_hit_ch_S++;
                    grp_hit_S[grp] = true;
                    if (c >= 48) xs[i] = true;
                    else         ys[i] = true;
                }
                if (ed.ph_C[idx] > T) {
                    if (c >= 48) xc[i] = true;
                    else         yc[i] = true;
                }
            }

            int n_grp_hit_S = 0;
            for (int g = 4; g <= 7; ++g) {
                if (grp_hit_S[g]) n_grp_hit_S++;
            }

            if (n_hit_ch_S >= 2 && n_grp_hit_S >= 2) {
                ev_2ch_2grp[i] = true;
            }
        }

        // Accumulate counts over N_pairs
        long long N_X_S = 0, N_X_C = 0, N_Y_S = 0, N_Y_C = 0;
        long long N_XY_same = 0;
        long long count_2ch_2grp = 0;

        for (Long64_t i = 0; i < N_pairs; ++i) {
            if (xs[i]) N_X_S++;
            if (xc[i]) N_X_C++;
            if (ys[i]) N_Y_S++;
            if (yc[i]) N_Y_C++;
            if (xs[i] && ys[i]) N_XY_same++;
            if (ev_2ch_2grp[i]) count_2ch_2grp++;
        }

        double aX = (N_d > 0) ? ((double)N_X_C / N_d) : 0.0;
        double aY = (N_d > 0) ? ((double)N_Y_C / N_d) : 0.0;
        double rX = (1.0 - aX > 1e-12) ? (((double)N_X_S / N_d - aX) / (1.0 - aX)) : 0.0;
        double rY = (1.0 - aY > 1e-12) ? (((double)N_Y_S / N_d - aY) / (1.0 - aY)) : 0.0;

        // Multi-offset mixing k = 1..20
        double sum_m = 0;
        for (int k = 1; k <= 20; ++k) {
            long long m = 0;
            for (Long64_t i = 0; i + k < n_used; ++i) {
                if (xs[i] && ys[i + k]) m++;
            }
            sum_m += (double)m;
        }
        double M_bar = sum_m / 20.0;

        double denom = (1.0 - aX) * (1.0 - aY);
        double N_XY_real = (std::abs(denom) > 1e-12) ? ((N_XY_same - M_bar) / denom + rX * rY * N_d) : 0.0;
        double sigma = (std::abs(denom) > 1e-12) ? (std::sqrt(N_XY_same + M_bar / 20.0) / denom) : 0.0;

        double ratio = (N_tot > 0) ? (N_XY_real / (0.198 * N_tot)) : 0.0;
        double sigma_ratio = (N_tot > 0) ? (sigma / (0.198 * N_tot)) : 0.0;

        RowOutput ro;
        ro.T = T;
        ro.N_X_S = N_X_S;
        ro.N_X_C = N_X_C;
        ro.N_Y_S = N_Y_S;
        ro.N_Y_C = N_Y_C;
        ro.N_XY_same = N_XY_same;
        ro.M_bar = M_bar;
        ro.aX = aX;
        ro.aY = aY;
        ro.rX = rX;
        ro.rY = rY;
        ro.N_XY_real = N_XY_real;
        ro.sigma = sigma;
        ro.ratio = ratio;
        ro.sigma_ratio = sigma_ratio;
        ro.n_ev_2ch_2grp = count_2ch_2grp;

        table1_rows.push_back(ro);
    }

    // =========================================================================
    // Table 1: Low-threshold scan
    // =========================================================================
    std::cout << "====================================================================================================================================================================================\n";
    std::cout << "Table 1: hod2019 X-Y Coincidences and Low-Threshold Scan (Run v8 Good Time Only, N = " << n_used << ", N_pairs = " << N_pairs << ")\n";
    std::cout << "====================================================================================================================================================================================\n";
    std::cout << std::right << std::setw(4) << "T"
              << std::setw(10) << "N_X_S"
              << std::setw(10) << "N_X_C"
              << std::setw(10) << "N_Y_S"
              << std::setw(10) << "N_Y_C"
              << std::setw(12) << "N_XY_same"
              << std::setw(12) << "M_bar"
              << std::setw(11) << "aX"
              << std::setw(11) << "aY"
              << std::setw(11) << "rX"
              << std::setw(11) << "rY"
              << std::setw(14) << "N_XY_real"
              << std::setw(12) << "sigma"
              << std::setw(24) << "N_XY_real/(0.198*N)"
              << std::setw(20) << "sigma/(0.198*N)"
              << std::setw(18) << "N_ev_>=2ch_>=2grp"
              << "\n";
    std::cout << std::string(180, '-') << "\n";

    for (const auto& r : table1_rows) {
        std::cout << std::right << std::setw(4) << (int)r.T
                  << std::setw(10) << r.N_X_S
                  << std::setw(10) << r.N_X_C
                  << std::setw(10) << r.N_Y_S
                  << std::setw(10) << r.N_Y_C
                  << std::setw(12) << r.N_XY_same
                  << std::fixed << std::setprecision(2)
                  << std::setw(12) << r.M_bar
                  << std::setprecision(6)
                  << std::setw(11) << r.aX
                  << std::setw(11) << r.aY
                  << std::setw(11) << r.rX
                  << std::setw(11) << r.rY
                  << std::setprecision(2)
                  << std::setw(14) << r.N_XY_real
                  << std::setw(12) << r.sigma
                  << std::setprecision(6)
                  << std::setw(24) << r.ratio
                  << std::setw(20) << r.sigma_ratio
                  << std::setw(18) << r.n_ev_2ch_2grp
                  << "\n";
    }
    std::cout << std::string(180, '-') << "\n\n";

    // Also CSV format for exact parsing
    std::cout << "CSV Table 1:\n";
    std::cout << "T,N_X_S,N_X_C,N_Y_S,N_Y_C,N_XY_same,M_bar,aX,aY,rX,rY,N_XY_real,sigma,N_XY_real/(0.198*N),sigma/(0.198*N),N_ev_2ch_2grp\n";
    for (const auto& r : table1_rows) {
        std::cout << (int)r.T << ","
                  << r.N_X_S << ","
                  << r.N_X_C << ","
                  << r.N_Y_S << ","
                  << r.N_Y_C << ","
                  << r.N_XY_same << ","
                  << std::fixed << std::setprecision(2) << r.M_bar << ","
                  << std::setprecision(6) << r.aX << ","
                  << r.aY << ","
                  << r.rX << ","
                  << r.rY << ","
                  << std::setprecision(2) << r.N_XY_real << ","
                  << r.sigma << ","
                  << std::setprecision(6) << r.ratio << ","
                  << r.sigma_ratio << ","
                  << r.n_ev_2ch_2grp << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // Table 2: Per-channel counts at T = 5 mV (channels 32..63)
    // =========================================================================
    std::cout << "========================================================================\n";
    std::cout << "Table 2: hod2019 Per-Channel Counts at T = 5 mV (Run v8 Good Time Only)\n";
    std::cout << "========================================================================\n";
    std::cout << std::right << std::setw(10) << "channel"
              << std::setw(12) << "N_S"
              << std::setw(12) << "N_C"
              << std::setw(14) << "N_S - N_C"
              << "\n";
    std::cout << std::string(48, '-') << "\n";

    struct ChCount {
        int ch;
        long long N_S;
        long long N_C;
    };
    std::vector<ChCount> ch_counts;

    for (int c = 32; c <= 63; ++c) {
        int idx = c - 32;
        long long n_s = 0;
        long long n_c = 0;
        for (Long64_t i = 0; i < n_used; ++i) {
            if (kept[i].ph_S[idx] > 5.0) n_s++;
            if (kept[i].ph_C[idx] > 5.0) n_c++;
        }
        ChCount cc{c, n_s, n_c};
        ch_counts.push_back(cc);
        std::cout << std::right << std::setw(10) << c
                  << std::setw(12) << n_s
                  << std::setw(12) << n_c
                  << std::setw(14) << (n_s - n_c)
                  << "\n";
    }
    std::cout << std::string(48, '-') << "\n\n";

    // Also CSV format
    std::cout << "CSV Table 2:\n";
    std::cout << "channel,N_S,N_C,N_S - N_C\n";
    for (const auto& cc : ch_counts) {
        std::cout << cc.ch << "," << cc.N_S << "," << cc.N_C << "," << (cc.N_S - cc.N_C) << "\n";
    }
}

