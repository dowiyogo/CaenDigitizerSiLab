// analysis/macros/hod2018_channel_coinc_fraction_v2.C
// TASK: ANALYSIS-FIX-16
// Measures per-channel coincidence fractions with all-pairs (C_all), weighted (C_w),
// and best-p.e. (C_best) coincidence estimators on run v8 good time.

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

struct ChannelRow {
    int ch;
    std::string plane;
    int strip;
    double G;
    double S_c;
    double C_all;
    double f_all;
    double boot_all;
    double C_w;
    double f_w;
    double boot_w;
    double C_best;
    double f_best;
    double boot_best;
};

double get_median(std::vector<double> vals) {
    if (vals.empty()) return 0.0;
    std::sort(vals.begin(), vals.end());
    size_t n = vals.size();
    if (n % 2 == 1) return vals[n / 2];
    else return 0.5 * (vals[n / 2 - 1] + vals[n / 2]);
}

void hod2018_channel_coinc_fraction_v2(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root") {
    gROOT->SetBatch(kTRUE);

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

    // Bootstrap setup (blocks of 1000 events, 200 replicas)
    const int kBlockSize = 1000;
    const int kReplicas = 200;
    int n_blocks = (n_v8 + kBlockSize - 1) / kBlockSize;
    std::vector<int> b_start(n_blocks), b_len(n_blocks);
    for (int b = 0; b < n_blocks; ++b) {
        b_start[b] = b * kBlockSize;
        b_len[b] = std::min(kBlockSize, (int)n_v8 - b_start[b]);
    }

    std::mt19937_64 rng(202610189);
    std::uniform_int_distribution<int> dist_block(0, n_blocks - 1);
    std::vector<std::vector<int>> rep_blocks(kReplicas, std::vector<int>(n_blocks));
    for (int r = 0; r < kReplicas; ++r) {
        for (int b = 0; b < n_blocks; ++b) {
            rep_blocks[r][b] = dist_block(rng);
        }
    }

    const std::vector<double> n_pe_list = {5.0, 10.0};

    // Helper struct for event-level hit information
    struct EventHitInfo {
        std::vector<int> fx; // strips in 18X over thr
        std::vector<int> fy; // strips in 18Y over thr
        int best_x;          // strip with max (ph - p0)/G in 18X (-1 if none)
        int best_y;          // strip with max (ph - p0)/G in 18Y (-1 if none)
        bool pl_xc;
        bool pl_yc;
        bool ch_s[32];
        bool ch_c[32];
    };

    for (double n_pe : n_pe_list) {
        std::cout << "=======================================================================================================================================================\n";
        std::cout << "Threshold n = " << (int)n_pe << " p.e. (Run v8 Good Time Only, N = " << n_v8 << ")\n";
        std::cout << "=======================================================================================================================================================\n";

        std::vector<EventHitInfo> ev_hits(n_v8);

        long long n_ev_nx_ge1 = 0;
        long long sum_nx = 0;
        long long n_ev_ny_ge1 = 0;
        long long sum_ny = 0;

        for (Long64_t i = 0; i < n_v8; ++i) {
            auto& eh = ev_hits[i];
            eh.best_x = -1;
            eh.best_y = -1;
            eh.pl_xc = false;
            eh.pl_yc = false;
            for (int c = 0; c < 32; ++c) {
                eh.ch_s[c] = false;
                eh.ch_c[c] = false;
            }

            float max_pe_x = -1e9f, max_pe_y = -1e9f;

            for (int c = 0; c < 32; ++c) {
                if (c == 3) continue;
                double thr = get_thr(c, n_pe);
                float ph_s = kept_events[i].ph_s[c];
                float ph_c = kept_events[i].ph_c[c];

                if (ph_s > thr) {
                    eh.ch_s[c] = true;
                    int strip = ch_to_strip(c);
                    float pe_s = (ph_s - (float)calib[c].p0) / (float)calib[c].G;

                    if (c < 16) {
                        eh.fx.push_back(strip);
                        if (pe_s > max_pe_x) {
                            max_pe_x = pe_s;
                            eh.best_x = strip;
                        }
                    } else {
                        eh.fy.push_back(strip);
                        if (pe_s > max_pe_y) {
                            max_pe_y = pe_s;
                            eh.best_y = strip;
                        }
                    }
                }

                if (ph_c > thr) {
                    eh.ch_c[c] = true;
                    if (c < 16) eh.pl_xc = true;
                    else        eh.pl_yc = true;
                }
            }

            int nx = (int)eh.fx.size();
            int ny = (int)eh.fy.size();
            if (nx >= 1) {
                n_ev_nx_ge1++;
                sum_nx += nx;
            }
            if (ny >= 1) {
                n_ev_ny_ge1++;
                sum_ny += ny;
            }
        }

        double mean_mult_X = (n_ev_nx_ge1 > 0) ? ((double)sum_nx / (double)n_ev_nx_ge1) : 0.0;
        double mean_mult_Y = (n_ev_ny_ge1 > 0) ? ((double)sum_ny / (double)n_ev_ny_ge1) : 0.0;

        // Plane accidentals over N_pairs
        Long64_t N_pairs = n_v8 - 1;
        long long n_xc_tot = 0, n_yc_tot = 0;
        for (Long64_t i = 0; i < N_pairs; ++i) {
            if (ev_hits[i].pl_xc) n_xc_tot++;
            if (ev_hits[i].pl_yc) n_yc_tot++;
        }
        double aX_nom = (double)n_xc_tot / (double)N_pairs;
        double aY_nom = (double)n_yc_tot / (double)N_pairs;
        double denom_nom = (1.0 - aX_nom) * (1.0 - aY_nom);

        // Coincidence matrices for nominal dataset
        double M_same_all[16][16] = {}, M_mixed_all[16][16] = {};
        double M_same_w[16][16] = {},   M_mixed_w[16][16] = {};
        double M_same_best[16][16] = {}, M_mixed_best[16][16] = {};

        for (Long64_t i = 0; i < n_v8; ++i) {
            const auto& eh = ev_hits[i];
            int nx = (int)eh.fx.size();
            int ny = (int)eh.fy.size();

            // Same event
            if (nx >= 1 && ny >= 1) {
                double w_s = 1.0 / (double)(nx * ny);
                for (int sx : eh.fx) {
                    for (int sy : eh.fy) {
                        M_same_all[sx][sy] += 1.0;
                        M_same_w[sx][sy] += w_s;
                    }
                }
            }
            if (eh.best_x >= 0 && eh.best_y >= 0) {
                M_same_best[eh.best_x][eh.best_y] += 1.0;
            }

            // Mixed events k = 1..20
            for (int k = 1; k <= 20; ++k) {
                if (i + k < n_v8) {
                    const auto& eh_k = ev_hits[i + k];
                    int ny_k = (int)eh_k.fy.size();

                    if (nx >= 1 && ny_k >= 1) {
                        double w_m = 1.0 / (double)(nx * ny_k);
                        for (int sx : eh.fx) {
                            for (int sy : eh_k.fy) {
                                M_mixed_all[sx][sy] += 1.0;
                                M_mixed_w[sx][sy] += w_m;
                            }
                        }
                    }
                    if (eh.best_x >= 0 && eh_k.best_y >= 0) {
                        M_mixed_best[eh.best_x][eh_k.best_y] += 1.0;
                    }
                }
            }
        }

        // Net matrices and profile sums
        double row_all[16] = {}, col_all[16] = {};
        double row_w[16] = {},   col_w[16] = {};
        double row_best[16] = {}, col_best[16] = {};

        for (int x = 0; x < 16; ++x) {
            for (int y = 0; y < 16; ++y) {
                double net_all = M_same_all[x][y] - M_mixed_all[x][y] / 20.0;
                row_all[x] += net_all;
                col_all[y] += net_all;

                double net_w = M_same_w[x][y] - M_mixed_w[x][y] / 20.0;
                row_w[x] += net_w;
                col_w[y] += net_w;

                double net_best = M_same_best[x][y] - M_mixed_best[x][y] / 20.0;
                row_best[x] += net_best;
                col_best[y] += net_best;
            }
        }

        // Divide by denom
        for (int k = 0; k < 16; ++k) {
            row_all[k] /= denom_nom;  col_all[k] /= denom_nom;
            row_w[k] /= denom_nom;    col_w[k] /= denom_nom;
            row_best[k] /= denom_nom; col_best[k] /= denom_nom;
        }

        // Per-channel single hits S_c and nominal rows
        std::vector<ChannelRow> table_rows;
        for (int c = 0; c < 32; ++c) {
            if (c == 3) continue;
            int strip = ch_to_strip(c);
            long long ns = 0, nc = 0;
            for (Long64_t i = 0; i < n_v8; ++i) {
                if (ev_hits[i].ch_s[c]) ns++;
                if (ev_hits[i].ch_c[c]) nc++;
            }
            double a_c = (double)nc / (double)n_v8;
            double S_c = (1.0 - a_c > 1e-12) ? ((double)(ns - nc) / (1.0 - a_c)) : 0.0;

            double C_all_val = (c < 16) ? row_all[strip] : col_all[strip];
            double C_w_val   = (c < 16) ? row_w[strip]   : col_w[strip];
            double C_best_val = (c < 16) ? row_best[strip] : col_best[strip];

            double f_all_val = (S_c > 0) ? (C_all_val / S_c) : 0.0;
            double f_w_val   = (S_c > 0) ? (C_w_val / S_c) : 0.0;
            double f_best_val = (S_c > 0) ? (C_best_val / S_c) : 0.0;

            ChannelRow cr{};
            cr.ch = c;
            cr.plane = calib[c].plane;
            cr.strip = strip;
            cr.G = calib[c].G;
            cr.S_c = S_c;
            cr.C_all = C_all_val;
            cr.f_all = f_all_val;
            cr.boot_all = 0.0;
            cr.C_w = C_w_val;
            cr.f_w = f_w_val;
            cr.boot_w = 0.0;
            cr.C_best = C_best_val;
            cr.f_best = f_best_val;
            cr.boot_best = 0.0;
            table_rows.push_back(cr);
        }

        // 2. Bootstrap replicas (200 replicas)
        size_t n_chans = table_rows.size();
        std::vector<std::vector<double>> rep_f_all(n_chans, std::vector<double>(kReplicas, 0.0));
        std::vector<std::vector<double>> rep_f_w(n_chans, std::vector<double>(kReplicas, 0.0));
        std::vector<std::vector<double>> rep_f_best(n_chans, std::vector<double>(kReplicas, 0.0));

        for (int r = 0; r < kReplicas; ++r) {
            std::vector<int> r_idx_list;
            r_idx_list.reserve(n_v8);
            for (int ib = 0; ib < n_blocks; ++ib) {
                int b_chosen = rep_blocks[r][ib];
                int st = b_start[b_chosen];
                int len = b_len[b_chosen];
                for (int j = 0; j < len; ++j) {
                    r_idx_list.push_back(st + j);
                }
            }

            Long64_t sz = r_idx_list.size();
            Long64_t r_npairs = sz - 1;

            long long r_n_xc = 0, r_n_yc = 0;
            for (Long64_t i = 0; i < r_npairs; ++i) {
                int idx = r_idx_list[i];
                if (ev_hits[idx].pl_xc) r_n_xc++;
                if (ev_hits[idx].pl_yc) r_n_yc++;
            }
            double r_aX = (double)r_n_xc / (double)r_npairs;
            double r_aY = (double)r_n_yc / (double)r_npairs;
            double r_denom = (1.0 - r_aX) * (1.0 - r_aY);

            double r_M_same_all[16][16] = {}, r_M_mixed_all[16][16] = {};
            double r_M_same_w[16][16] = {},   r_M_mixed_w[16][16] = {};
            double r_M_same_best[16][16] = {}, r_M_mixed_best[16][16] = {};

            for (Long64_t i = 0; i < sz; ++i) {
                int idx = r_idx_list[i];
                const auto& eh = ev_hits[idx];
                int nx = (int)eh.fx.size();
                int ny = (int)eh.fy.size();

                if (nx >= 1 && ny >= 1) {
                    double w_s = 1.0 / (double)(nx * ny);
                    for (int sx : eh.fx) {
                        for (int sy : eh.fy) {
                            r_M_same_all[sx][sy] += 1.0;
                            r_M_same_w[sx][sy] += w_s;
                        }
                    }
                }
                if (eh.best_x >= 0 && eh.best_y >= 0) {
                    r_M_same_best[eh.best_x][eh.best_y] += 1.0;
                }

                for (int k = 1; k <= 20; ++k) {
                    if (i + k < sz) {
                        int idx_k = r_idx_list[i + k];
                        const auto& eh_k = ev_hits[idx_k];
                        int ny_k = (int)eh_k.fy.size();

                        if (nx >= 1 && ny_k >= 1) {
                            double w_m = 1.0 / (double)(nx * ny_k);
                            for (int sx : eh.fx) {
                                for (int sy : eh_k.fy) {
                                    r_M_mixed_all[sx][sy] += 1.0;
                                    r_M_mixed_w[sx][sy] += w_m;
                                }
                            }
                        }
                        if (eh.best_x >= 0 && eh_k.best_y >= 0) {
                            r_M_mixed_best[eh.best_x][eh_k.best_y] += 1.0;
                        }
                    }
                }
            }

            double r_row_all[16] = {}, r_col_all[16] = {};
            double r_row_w[16] = {},   r_col_w[16] = {};
            double r_row_best[16] = {}, r_col_best[16] = {};

            for (int x = 0; x < 16; ++x) {
                for (int y = 0; y < 16; ++y) {
                    r_row_all[x] += (r_M_same_all[x][y] - r_M_mixed_all[x][y] / 20.0);
                    r_col_all[y] += (r_M_same_all[x][y] - r_M_mixed_all[x][y] / 20.0);

                    r_row_w[x] += (r_M_same_w[x][y] - r_M_mixed_w[x][y] / 20.0);
                    r_col_w[y] += (r_M_same_w[x][y] - r_M_mixed_w[x][y] / 20.0);

                    r_row_best[x] += (r_M_same_best[x][y] - r_M_mixed_best[x][y] / 20.0);
                    r_col_best[y] += (r_M_same_best[x][y] - r_M_mixed_best[x][y] / 20.0);
                }
            }

            for (int k = 0; k < 16; ++k) {
                r_row_all[k] /= r_denom;  r_col_all[k] /= r_denom;
                r_row_w[k] /= r_denom;    r_col_w[k] /= r_denom;
                r_row_best[k] /= r_denom; r_col_best[k] /= r_denom;
            }

            for (size_t ic = 0; ic < n_chans; ++ic) {
                int c = table_rows[ic].ch;
                int strip = table_rows[ic].strip;

                long long r_ns = 0, r_nc = 0;
                for (Long64_t i = 0; i < sz; ++i) {
                    int idx = r_idx_list[i];
                    if (ev_hits[idx].ch_s[c]) r_ns++;
                    if (ev_hits[idx].ch_c[c]) r_nc++;
                }
                double r_ac = (double)r_nc / (double)sz;
                double r_Sc = (1.0 - r_ac > 1e-12) ? ((double)(r_ns - r_nc) / (1.0 - r_ac)) : 0.0;

                double r_Call = (c < 16) ? r_row_all[strip] : r_col_all[strip];
                double r_Cw   = (c < 16) ? r_row_w[strip]   : r_col_w[strip];
                double r_Cbest = (c < 16) ? r_row_best[strip] : r_col_best[strip];

                rep_f_all[ic][r]  = (r_Sc > 0) ? (r_Call / r_Sc) : 0.0;
                rep_f_w[ic][r]    = (r_Sc > 0) ? (r_Cw / r_Sc) : 0.0;
                rep_f_best[ic][r] = (r_Sc > 0) ? (r_Cbest / r_Sc) : 0.0;
            }
        }

        // Bootstrap RMS
        for (size_t ic = 0; ic < n_chans; ++ic) {
            auto calc_rms = [&](const std::vector<double>& vec) -> double {
                double s = 0; for (double v : vec) s += v;
                double m = s / (double)vec.size();
                double sq = 0; for (double v : vec) sq += (v - m) * (v - m);
                return std::sqrt(sq / (double)(vec.size() - 1));
            };
            table_rows[ic].boot_all  = calc_rms(rep_f_all[ic]);
            table_rows[ic].boot_w    = calc_rms(rep_f_w[ic]);
            table_rows[ic].boot_best = calc_rms(rep_f_best[ic]);
        }

        // Print Main Table
        // Columns EXACTLY:
        // c plane strip G_c S_c C_all f_all boot_all C_w f_w boot_w C_best f_best boot_best
        std::cout << std::left << std::setw(3) << "c"
                  << std::setw(6) << "plane"
                  << std::right << std::setw(6) << "strip"
                  << std::setw(8) << "G_c"
                  << std::setw(9) << "S_c"
                  << std::setw(10) << "C_all"
                  << std::setw(9) << "f_all"
                  << std::setw(10) << "boot_all"
                  << std::setw(10) << "C_w"
                  << std::setw(9) << "f_w"
                  << std::setw(10) << "boot_w"
                  << std::setw(10) << "C_best"
                  << std::setw(9) << "f_best"
                  << std::setw(11) << "boot_best"
                  << "\n";
        std::cout << std::string(120, '-') << "\n";

        for (const auto& r : table_rows) {
            std::cout << std::left << std::setw(3) << r.ch
                      << std::setw(6) << r.plane
                      << std::right << std::setw(6) << r.strip
                      << std::fixed << std::setprecision(3)
                      << std::setw(8) << r.G
                      << std::setprecision(1)
                      << std::setw(9) << r.S_c
                      << std::setw(10) << r.C_all
                      << std::setprecision(4)
                      << std::setw(9) << r.f_all
                      << std::setw(10) << r.boot_all
                      << std::setprecision(1)
                      << std::setw(10) << r.C_w
                      << std::setprecision(4)
                      << std::setw(9) << r.f_w
                      << std::setw(10) << r.boot_w
                      << std::setprecision(1)
                      << std::setw(10) << r.C_best
                      << std::setprecision(4)
                      << std::setw(9) << r.f_best
                      << std::setw(11) << r.boot_best
                      << "\n";
        }
        std::cout << std::string(120, '-') << "\n\n";

        // Plane summaries
        for (const std::string& pl : {"18X", "18Y"}) {
            std::vector<double> v_fall, v_fw, v_fbest;
            std::vector<double> v_boot_w, v_boot_best;
            double sum_Cw = 0.0, sum_Sc = 0.0;

            for (const auto& r : table_rows) {
                if (r.plane == pl) {
                    v_fall.push_back(r.f_all);
                    v_fw.push_back(r.f_w);
                    v_fbest.push_back(r.f_best);
                    v_boot_w.push_back(r.boot_w);
                    v_boot_best.push_back(r.boot_best);
                    sum_Cw += r.C_w;
                    sum_Sc += r.S_c;
                }
            }

            double med_all  = get_median(v_fall);
            double med_w    = get_median(v_fw);
            double med_best = get_median(v_fbest);

            // Weighted mean and chi2/ndf for f_w
            double sum_w_w = 0, sum_wx_w = 0;
            for (size_t i = 0; i < v_fw.size(); ++i) {
                double w = 1.0 / (v_boot_w[i] * v_boot_w[i]);
                sum_w_w += w;
                sum_wx_w += w * v_fw[i];
            }
            double mean_w = sum_wx_w / sum_w_w;
            double chi2_w = 0;
            for (size_t i = 0; i < v_fw.size(); ++i) {
                double diff = v_fw[i] - mean_w;
                chi2_w += (diff * diff) / (v_boot_w[i] * v_boot_w[i]);
            }
            double ndf = (double)(v_fw.size() - 1);
            double chi2_ndf_w = chi2_w / ndf;

            // Weighted mean and chi2/ndf for f_best
            double sum_w_best = 0, sum_wx_best = 0;
            for (size_t i = 0; i < v_fbest.size(); ++i) {
                double w = 1.0 / (v_boot_best[i] * v_boot_best[i]);
                sum_w_best += w;
                sum_wx_best += w * v_fbest[i];
            }
            double mean_best = sum_wx_best / sum_w_best;
            double chi2_best = 0;
            for (size_t i = 0; i < v_fbest.size(); ++i) {
                double diff = v_fbest[i] - mean_best;
                chi2_best += (diff * diff) / (v_boot_best[i] * v_boot_best[i]);
            }
            double chi2_ndf_best = chi2_best / ndf;

            double mean_mult = (pl == "18X") ? mean_mult_X : mean_mult_Y;

            std::cout << "Summary for Plane " << pl << " at n = " << (int)n_pe << " p.e. (" << v_fw.size() << " live channels):\n";
            std::cout << "  median f_all  = " << std::fixed << std::setprecision(4) << med_all << "\n";
            std::cout << "  median f_w    = " << med_w << "\n";
            std::cout << "  median f_best = " << med_best << "\n";
            std::cout << "  chi2/ndf (f_w against weighted mean " << std::setprecision(4) << mean_w << "): "
                      << std::setprecision(2) << chi2_w << " / " << (int)ndf << " = " << std::setprecision(2) << chi2_ndf_w << "\n";
            std::cout << "  chi2/ndf (f_best against weighted mean " << std::setprecision(4) << mean_best << "): "
                      << std::setprecision(2) << chi2_best << " / " << (int)ndf << " = " << std::setprecision(2) << chi2_ndf_best << "\n";
            std::cout << "  sum_c C_w     = " << std::fixed << std::setprecision(1) << sum_Cw << "\n";
            std::cout << "  sum_c S_c     = " << sum_Sc << "\n";
            if (pl == "18X") {
                std::cout << "  mean multiplicity <nX> (nX >= 1) = " << std::fixed << std::setprecision(3) << mean_mult << "\n";
            } else {
                std::cout << "  mean multiplicity <nY> (nY >= 1) = " << std::fixed << std::setprecision(3) << mean_mult << "\n";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

