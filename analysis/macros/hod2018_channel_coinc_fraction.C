// analysis/macros/hod2018_channel_coinc_fraction.C
// TASK: ANALYSIS-FIX-15
// Measures per-channel single hits S_c, coincidence sum C_c, coincidence fraction f_c = C_c / S_c,
// bootstrap RMS of f_c, and f_c / median(plane) at n = 5 and n = 10 p.e. on run v8 good time.

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

struct ChannelResult {
    int ch;
    std::string plane;
    int strip;
    double G;
    int status;
    long long N_S;
    long long N_C;
    double a_c;
    double S_c;
    double C_c;
    double f_c;
    double boot_rms_f;
    double f_ratio_to_median;
};

double get_median(std::vector<double> vals) {
    if (vals.empty()) return 0.0;
    std::sort(vals.begin(), vals.end());
    size_t n = vals.size();
    if (n % 2 == 1) return vals[n / 2];
    else return 0.5 * (vals[n / 2 - 1] + vals[n / 2]);
}

void hod2018_channel_coinc_fraction(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root") {
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

    for (double n_pe : n_pe_list) {
        std::cout << "=======================================================================================================================================================\n";
        std::cout << "Threshold n = " << (int)n_pe << " p.e. (Run v8 Good Time Only, N = " << n_v8 << ")\n";
        std::cout << "=======================================================================================================================================================\n";

        // 1. Single channel hit flags & plane hit flags for nominal v8
        std::vector<bool> ch_hit_s[32], ch_hit_c[32];
        for (int c = 0; c < 32; ++c) {
            ch_hit_s[c].assign(n_v8, false);
            ch_hit_c[c].assign(n_v8, false);
        }
        std::vector<int> best_x(n_v8, -1), best_y(n_v8, -1);
        std::vector<bool> pl_xs(n_v8, false), pl_xc(n_v8, false), pl_ys(n_v8, false), pl_yc(n_v8, false);

        for (Long64_t i = 0; i < n_v8; ++i) {
            float max_ph_x = -1e9, max_ph_y = -1e9;
            for (int c = 0; c < 32; ++c) {
                if (c == 3) continue;
                double thr = get_thr(c, n_pe);
                float ph_s = kept_events[i].ph_s[c];
                float ph_c = kept_events[i].ph_c[c];

                if (ph_s > thr) {
                    ch_hit_s[c][i] = true;
                    int strip = ch_to_strip(c);
                    if (c < 16) {
                        pl_xs[i] = true;
                        if (ph_s > max_ph_x) {
                            max_ph_x = ph_s;
                            best_x[i] = strip;
                        }
                    } else {
                        pl_ys[i] = true;
                        if (ph_s > max_ph_y) {
                            max_ph_y = ph_s;
                            best_y[i] = strip;
                        }
                    }
                }
                if (ph_c > thr) {
                    ch_hit_c[c][i] = true;
                    if (c < 16) pl_xc[i] = true;
                    else        pl_yc[i] = true;
                }
            }
        }

        // Nominal plane accidental probabilities over N_pairs
        Long64_t N_pairs = n_v8 - 1;
        long long n_xc_tot = 0, n_yc_tot = 0;
        for (Long64_t i = 0; i < N_pairs; ++i) {
            if (pl_xc[i]) n_xc_tot++;
            if (pl_yc[i]) n_yc_tot++;
        }
        double aX_nom = (double)n_xc_tot / (double)N_pairs;
        double aY_nom = (double)n_yc_tot / (double)N_pairs;
        double denom_nom = (1.0 - aX_nom) * (1.0 - aY_nom);

        // Nominal coincidence matrix
        double M_same[16][16] = {};
        double M_mixed_sum[16][16] = {};
        for (Long64_t i = 0; i < n_v8; ++i) {
            int sx = best_x[i];
            int sy = best_y[i];
            if (sx >= 0 && sx < 16 && sy >= 0 && sy < 16) {
                M_same[sx][sy] += 1.0;
            }
            if (sx >= 0 && sx < 16) {
                for (int k = 1; k <= 20; ++k) {
                    if (i + k < n_v8) {
                        int sy_k = best_y[i + k];
                        if (sy_k >= 0 && sy_k < 16) {
                            M_mixed_sum[sx][sy_k] += 1.0;
                        }
                    }
                }
            }
        }

        double row_sums[16] = {};
        double col_sums[16] = {};
        for (int x = 0; x < 16; ++x) {
            for (int y = 0; y < 16; ++y) {
                double net = M_same[x][y] - M_mixed_sum[x][y] / 20.0;
                row_sums[x] += net;
                col_sums[y] += net;
            }
        }

        double corr_row_sums[16] = {};
        double corr_col_sums[16] = {};
        for (int x = 0; x < 16; ++x) corr_row_sums[x] = row_sums[x] / denom_nom;
        for (int y = 0; y < 16; ++y) corr_col_sums[y] = col_sums[y] / denom_nom;

        // Nominal channel results
        std::vector<ChannelResult> chan_results;
        for (int c = 0; c < 32; ++c) {
            if (c == 3) continue; // dead channel
            int strip = ch_to_strip(c);
            long long ns = 0, nc = 0;
            for (Long64_t i = 0; i < n_v8; ++i) {
                if (ch_hit_s[c][i]) ns++;
                if (ch_hit_c[c][i]) nc++;
            }
            double a_c = (double)nc / (double)n_v8;
            double S_c = (1.0 - a_c > 1e-12) ? ((double)(ns - nc) / (1.0 - a_c)) : 0.0;
            double C_c = (c < 16) ? corr_row_sums[strip] : corr_col_sums[strip];
            double f_c = (S_c > 0) ? (C_c / S_c) : 0.0;

            ChannelResult cr{};
            cr.ch = c;
            cr.plane = calib[c].plane;
            cr.strip = strip;
            cr.G = calib[c].G;
            cr.status = calib[c].status;
            cr.N_S = ns;
            cr.N_C = nc;
            cr.a_c = a_c;
            cr.S_c = S_c;
            cr.C_c = C_c;
            cr.f_c = f_c;
            cr.boot_rms_f = 0.0;
            cr.f_ratio_to_median = 0.0;
            chan_results.push_back(cr);
        }

        // 2. Bootstrap replicas (200 replicas)
        std::vector<std::vector<double>> rep_fc(chan_results.size(), std::vector<double>(kReplicas, 0.0));

        for (int r = 0; r < kReplicas; ++r) {
            std::vector<int> r_best_x, r_best_y;
            std::vector<bool> r_pl_xc, r_pl_yc;
            std::vector<std::vector<bool>> r_ch_s(chan_results.size()), r_ch_c(chan_results.size());

            r_best_x.reserve(n_v8); r_best_y.reserve(n_v8);
            r_pl_xc.reserve(n_v8); r_pl_yc.reserve(n_v8);
            for (size_t ic = 0; ic < chan_results.size(); ++ic) {
                r_ch_s[ic].reserve(n_v8);
                r_ch_c[ic].reserve(n_v8);
            }

            for (int ib = 0; ib < n_blocks; ++ib) {
                int b_chosen = rep_blocks[r][ib];
                int st = b_start[b_chosen];
                int len = b_len[b_chosen];
                for (int j = 0; j < len; ++j) {
                    int idx = st + j;
                    r_best_x.push_back(best_x[idx]);
                    r_best_y.push_back(best_y[idx]);
                    r_pl_xc.push_back(pl_xc[idx]);
                    r_pl_yc.push_back(pl_yc[idx]);
                    for (size_t ic = 0; ic < chan_results.size(); ++ic) {
                        int c = chan_results[ic].ch;
                        r_ch_s[ic].push_back(ch_hit_s[c][idx]);
                        r_ch_c[ic].push_back(ch_hit_c[c][idx]);
                    }
                }
            }

            Long64_t sz = r_best_x.size();
            Long64_t r_npairs = sz - 1;

            long long r_n_xc = 0, r_n_yc = 0;
            for (Long64_t i = 0; i < r_npairs; ++i) {
                if (r_pl_xc[i]) r_n_xc++;
                if (r_pl_yc[i]) r_n_yc++;
            }
            double r_aX = (double)r_n_xc / (double)r_npairs;
            double r_aY = (double)r_n_yc / (double)r_npairs;
            double r_denom = (1.0 - r_aX) * (1.0 - r_aY);

            double r_M_same[16][16] = {};
            double r_M_mixed_sum[16][16] = {};
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

            double r_row_sums[16] = {};
            double r_col_sums[16] = {};
            for (int x = 0; x < 16; ++x) {
                for (int y = 0; y < 16; ++y) {
                    double net = r_M_same[x][y] - r_M_mixed_sum[x][y] / 20.0;
                    r_row_sums[x] += net;
                    r_col_sums[y] += net;
                }
            }

            for (size_t ic = 0; ic < chan_results.size(); ++ic) {
                int c = chan_results[ic].ch;
                int strip = chan_results[ic].strip;
                long long r_ns = 0, r_nc = 0;
                for (Long64_t i = 0; i < sz; ++i) {
                    if (r_ch_s[ic][i]) r_ns++;
                    if (r_ch_c[ic][i]) r_nc++;
                }
                double r_ac = (double)r_nc / (double)sz;
                double r_Sc = (1.0 - r_ac > 1e-12) ? ((double)(r_ns - r_nc) / (1.0 - r_ac)) : 0.0;
                double r_Cc = (c < 16) ? (r_row_sums[strip] / r_denom) : (r_col_sums[strip] / r_denom);
                rep_fc[ic][r] = (r_Sc > 0) ? (r_Cc / r_Sc) : 0.0;
            }
        }

        // Calculate bootstrap RMS of f_c
        for (size_t ic = 0; ic < chan_results.size(); ++ic) {
            double sum_f = 0;
            for (int r = 0; r < kReplicas; ++r) sum_f += rep_fc[ic][r];
            double m_f = sum_f / (double)kReplicas;

            double sum_sq = 0;
            for (int r = 0; r < kReplicas; ++r) sum_sq += (rep_fc[ic][r] - m_f) * (rep_fc[ic][r] - m_f);
            chan_results[ic].boot_rms_f = std::sqrt(sum_sq / (double)(kReplicas - 1));
        }

        // 3. Medians per plane
        std::vector<double> f_18x, f_18y;
        for (const auto& cr : chan_results) {
            if (cr.ch < 16) f_18x.push_back(cr.f_c);
            else            f_18y.push_back(cr.f_c);
        }
        double med_18x = get_median(f_18x);
        double med_18y = get_median(f_18y);

        for (auto& cr : chan_results) {
            double med = (cr.ch < 16) ? med_18x : med_18y;
            cr.f_ratio_to_median = (med > 0) ? (cr.f_c / med) : 0.0;
        }

        // 4. Print Table
        std::cout << std::left << std::setw(4) << "c"
                  << std::setw(7) << "plane"
                  << std::right << std::setw(6) << "strip"
                  << std::setw(8) << "G_c"
                  << std::setw(8) << "status"
                  << std::setw(12) << "S_c"
                  << std::setw(12) << "C_c"
                  << std::setw(10) << "f_c"
                  << std::setw(14) << "boot_RMS"
                  << std::setw(16) << "f_c/median"
                  << "\n";
        std::cout << std::string(97, '-') << "\n";

        for (const auto& cr : chan_results) {
            std::cout << std::left << std::setw(4) << cr.ch
                      << std::setw(7) << cr.plane
                      << std::right << std::setw(6) << cr.strip
                      << std::fixed << std::setprecision(3)
                      << std::setw(8) << cr.G
                      << std::setw(8) << cr.status
                      << std::setprecision(1)
                      << std::setw(12) << cr.S_c
                      << std::setw(12) << cr.C_c
                      << std::setprecision(4)
                      << std::setw(10) << cr.f_c
                      << std::setw(14) << cr.boot_rms_f
                      << std::setw(16) << cr.f_ratio_to_median
                      << "\n";
        }
        std::cout << std::string(97, '-') << "\n";

        std::cout << "\nPlane Summary at n = " << (int)n_pe << " p.e.:\n";
        std::cout << "  Plane 18X (15 live channels): median f_c = " << std::fixed << std::setprecision(4) << med_18x << "\n";
        std::cout << "  Plane 18Y (16 live channels): median f_c = " << std::fixed << std::setprecision(4) << med_18y << "\n\n";

        // Print per plane list of f_c / median explicitly
        std::cout << "Channels in 18X (c, strip, f_c, f_c / median):\n";
        for (const auto& cr : chan_results) {
            if (cr.ch < 16) {
                std::cout << "  ch " << std::setw(2) << cr.ch
                          << " (strip " << std::setw(2) << cr.strip << "): f_c = "
                          << std::fixed << std::setprecision(4) << cr.f_c
                          << " ± " << cr.boot_rms_f
                          << ", f_c / median = " << std::setprecision(4) << cr.f_ratio_to_median << "\n";
            }
        }
        std::cout << "\nChannels in 18Y (c, strip, f_c, f_c / median):\n";
        for (const auto& cr : chan_results) {
            if (cr.ch >= 16) {
                std::cout << "  ch " << std::setw(2) << cr.ch
                          << " (strip " << std::setw(2) << cr.strip << "): f_c = "
                          << std::fixed << std::setprecision(4) << cr.f_c
                          << " ± " << cr.boot_rms_f
                          << ", f_c / median = " << std::setprecision(4) << cr.f_ratio_to_median << "\n";
            }
        }
        std::cout << "\n\n";
    }
}

