#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>

void baseline_scan() {
    gROOT->SetBatch(kTRUE);

    const char* fname = "lectura_vx1740_ext_nim_1M_v8.root";
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening " << fname << "\n";
        return;
    }
    TTree* t = (TTree*)f->Get("data_dig");
    if (!t) {
        std::cerr << "Error reading data_dig\n";
        f->Close();
        return;
    }

    std::vector<float>* tw = nullptr;
    t->SetBranchAddress("time", &tw);
    std::vector<float>* cw[64] = {nullptr};
    for (int c = 0; c < 64; ++c) {
        char bn[16];
        std::snprintf(bn, sizeof(bn), "Ch%d", c);
        if (t->GetBranch(bn)) t->SetBranchAddress(bn, &cw[c]);
    }

    const int kH19X[16] = {53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56};
    const int kH19Y[16] = {37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40};

    std::vector<int> planes[4];
    for (int i = 0; i < 16; ++i) {
        planes[0].push_back(i);        // hod2018X
        planes[1].push_back(i + 16);   // hod2018Y
        planes[2].push_back(kH19X[i]); // hod2019X
        planes[3].push_back(kH19Y[i]); // hod2019Y
    }
    const char* plane_names[4] = {"hod2018X", "hod2018Y", "hod2019X", "hod2019Y"};

    // Time ranges for A1
    struct RangeDef {
        const char* name;
        double tlo, thi;
        bool hi_inclusive;
    };
    const RangeDef ranges[9] = {
        {"[0, 2500)",     0.0,    2500.0, false},
        {"[2500, 3400)",  2500.0, 3400.0, false},
        {"[3400, 3700)",  3400.0, 3700.0, false},
        {"[3700, 3950)",  3700.0, 3950.0, false},
        {"[3950, 4000)",  3950.0, 4000.0, false},
        {"[4000, 4400]",  4000.0, 4400.0, true},
        {"[4400, 5000)",  4400.0, 5000.0, false},
        {"[5000, 7000)",  5000.0, 7000.0, false},
        {"[7000, 9968]",  7000.0, 9968.0, true}
    };

    // Accumulators for A1, A2 (all triggers)
    // Per plane: sum of range averages, count of channel-events
    std::vector<double> sum_range_all[4];
    std::vector<long long> count_range_all[4];
    for (int pl = 0; pl < 4; ++pl) {
        sum_range_all[pl].assign(9, 0.0);
        count_range_all[pl].assign(9, 0);
    }
    std::vector<double> b_shift_ch_all[4];       // per channel-event
    std::vector<double> b_shift_plane_all[4];    // per event (plane-averaged)

    // Accumulators for A3 (restricted to ph_glob < 5 mV in [2500, 8000] ns)
    std::vector<double> sum_range_nopulse[4];
    std::vector<long long> count_range_nopulse[4];
    for (int pl = 0; pl < 4; ++pl) {
        sum_range_nopulse[pl].assign(9, 0.0);
        count_range_nopulse[pl].assign(9, 0);
    }
    std::vector<double> b_shift_ch_nopulse[4];
    std::vector<double> b_shift_plane_nopulse[4];

    // Accumulators for Step C: Q pedestal noise
    // For hod2018 (ch 0..31): evaluate at kNpost=8 and kNpost=12
    // For hod2019 (ch 32..63): evaluate at kNpost=22
    std::vector<double> Q_glob_nopulse[4];
    std::vector<double> Q_loc_nopulse[4];

    // Also with kNpost18 = 12 for hod2018 to compare
    std::vector<double> Q_glob_nopulse_n12[2];
    std::vector<double> Q_loc_nopulse_n12[2];

    const Long64_t nEntries = std::min((Long64_t)20000, t->GetEntries());
    std::cout << "=== Baseline Diagnosis (v8, first " << nEntries << " entries) ===\n\n";

    for (Long64_t iev = 0; iev < nEntries; ++iev) {
        t->GetEntry(iev);
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int nSamp = (int)tv.size();

        // Sample at 4200 ns: index with minimum |tv[i] - 4200.0|
        // For vx1740, 16.0 ns step: sample 262 = 4192 ns, sample 263 = 4208 ns.
        // Let's use 262 (closest round(4200/16))
        const int i0 = 262;

        for (int pl = 0; pl < 4; ++pl) {
            double ev_shift_sum_all = 0.0;
            int ev_shift_cnt_all = 0;
            double ev_shift_sum_nopulse = 0.0;
            int ev_shift_cnt_nopulse = 0;

            for (int c : planes[pl]) {
                if (!cw[c] || cw[c]->empty()) continue;
                const auto& wf = *cw[c];
                if ((int)wf.size() != nSamp) continue;

                // ped_glob: samples with t < 2500 ns
                double s_ped = 0.0; int n_ped = 0;
                for (int i = 0; i < nSamp; ++i) {
                    if (tv[i] < 2500.0) { s_ped += wf[i]; n_ped++; }
                }
                if (n_ped < 2) continue;
                double ped_glob = s_ped / n_ped;

                // Range averages for A1
                double r_avgs[9] = {0.0};
                for (int r = 0; r < 9; ++r) {
                    double s_r = 0.0; int n_r = 0;
                    for (int i = 0; i < nSamp; ++i) {
                        bool in_r = false;
                        if (ranges[r].hi_inclusive) in_r = (tv[i] >= ranges[r].tlo && tv[i] <= ranges[r].thi);
                        else {
                            if (r == 6) in_r = (tv[i] > ranges[r].tlo && tv[i] < ranges[r].thi);
                            else in_r = (tv[i] >= ranges[r].tlo && tv[i] < ranges[r].thi);
                        }
                        if (in_r) { s_r += (wf[i] - ped_glob); n_r++; }
                    }
                    if (n_r > 0) {
                        r_avgs[r] = s_r / n_r;
                        sum_range_all[pl][r] += r_avgs[r];
                        count_range_all[pl][r]++;
                    }
                }

                // b_shift for A2
                double b_shift = r_avgs[3]; // range 3 is [3700, 3950)
                b_shift_ch_all[pl].push_back(b_shift);
                ev_shift_sum_all += b_shift;
                ev_shift_cnt_all++;

                // Pulse check for A3: vmin in [2500, 8000] ns
                double vmin_pulse = 1e30;
                for (int i = 0; i < nSamp; ++i) {
                    if (tv[i] >= 2500.0 && tv[i] <= 8000.0) {
                        if (wf[i] < vmin_pulse) vmin_pulse = wf[i];
                    }
                }
                double ph_glob = ped_glob - vmin_pulse;

                if (ph_glob < 5.0) {
                    // No pulse in [2500, 8000] ns
                    for (int r = 0; r < 9; ++r) {
                        sum_range_nopulse[pl][r] += r_avgs[r];
                        count_range_nopulse[pl][r]++;
                    }
                    b_shift_ch_nopulse[pl].push_back(b_shift);
                    ev_shift_sum_nopulse += b_shift;
                    ev_shift_cnt_nopulse++;

                    // Step C: Q pedestal noise
                    // ped_loc: offsets -12..-4 from i0 (4200 ns)
                    double s_loc = 0.0; int n_loc = 0;
                    for (int k = -12; k <= -4; ++k) {
                        int idx = i0 + k;
                        if (idx >= 0 && idx < nSamp) { s_loc += wf[idx]; n_loc++; }
                    }
                    if (n_loc == 9) {
                        double ped_loc = s_loc / n_loc;
                        int kNpost = (pl < 2) ? 8 : 22; // default kNpost18=8, kNpost19=22

                        double Q_g = 0.0, Q_l = 0.0;
                        for (int k = -3; k <= kNpost; ++k) {
                            int idx = i0 + k;
                            if (idx >= 0 && idx < nSamp) {
                                Q_g += (ped_glob - wf[idx]) * 16.0;
                                Q_l += (ped_loc - wf[idx]) * 16.0;
                            }
                        }
                        Q_glob_nopulse[pl].push_back(Q_g);
                        Q_loc_nopulse[pl].push_back(Q_l);

                        // If hod2018, also evaluate with kNpost=12
                        if (pl < 2) {
                            double Q_g12 = 0.0, Q_l12 = 0.0;
                            for (int k = -3; k <= 12; ++k) {
                                int idx = i0 + k;
                                if (idx >= 0 && idx < nSamp) {
                                    Q_g12 += (ped_glob - wf[idx]) * 16.0;
                                    Q_l12 += (ped_loc - wf[idx]) * 16.0;
                                }
                            }
                            Q_glob_nopulse_n12[pl].push_back(Q_g12);
                            Q_loc_nopulse_n12[pl].push_back(Q_l12);
                        }
                    }
                }
            }

            if (ev_shift_cnt_all > 0) {
                b_shift_plane_all[pl].push_back(ev_shift_sum_all / ev_shift_cnt_all);
            }
            if (ev_shift_cnt_nopulse > 0) {
                b_shift_plane_nopulse[pl].push_back(ev_shift_sum_nopulse / ev_shift_cnt_nopulse);
            }
        }
    }

    f->Close();

    // Helper for stats
    auto get_stats = [](const std::vector<double>& v, double& mean, double& rms) {
        if (v.empty()) { mean = 0; rms = 0; return; }
        double sum = 0.0;
        for (double x : v) sum += x;
        mean = sum / v.size();
        double sum2 = 0.0;
        for (double x : v) sum2 += (x - mean) * (x - mean);
        rms = std::sqrt(sum2 / v.size());
    };

    // Print A1
    std::cout << "## Step A1: Mean of (wf - ped_glob) [mV] per plane across time ranges (All Triggers)\n";
    std::cout << std::setw(15) << std::left << "Range [ns]";
    for (int pl = 0; pl < 4; ++pl) std::cout << std::setw(14) << plane_names[pl];
    std::cout << "\n" << std::string(71, '-') << "\n";
    for (int r = 0; r < 9; ++r) {
        std::cout << std::setw(15) << std::left << ranges[r].name;
        for (int pl = 0; pl < 4; ++pl) {
            double avg = (count_range_all[pl][r] > 0) ? (sum_range_all[pl][r] / count_range_all[pl][r]) : 0.0;
            std::cout << std::setw(14) << std::fixed << std::setprecision(4) << avg;
        }
        std::cout << "\n";
    }

    // Print A2
    std::cout << "\n## Step A2: b_shift = mean(wf in [3700, 3950]) - ped_glob [mV] (All Triggers)\n";
    std::cout << std::setw(14) << std::left << "Plane"
              << std::setw(14) << "Mean [mV]"
              << std::setw(18) << "RMS (events)"
              << std::setw(18) << "RMS (ch-events)"
              << std::setw(12) << "N_events"
              << "\n";
    std::cout << std::string(76, '-') << "\n";
    for (int pl = 0; pl < 4; ++pl) {
        double m_ev = 0, rms_ev = 0;
        get_stats(b_shift_plane_all[pl], m_ev, rms_ev);
        double m_ch = 0, rms_ch = 0;
        get_stats(b_shift_ch_all[pl], m_ch, rms_ch);
        std::cout << std::setw(14) << std::left << plane_names[pl]
                  << std::setw(14) << std::fixed << std::setprecision(4) << m_ev
                  << std::setw(18) << std::fixed << std::setprecision(4) << rms_ev
                  << std::setw(18) << std::fixed << std::setprecision(4) << rms_ch
                  << std::setw(12) << b_shift_plane_all[pl].size()
                  << "\n";
    }

    // Print A3
    std::cout << "\n## Step A3: Restricted to events with NO PULSE (ph_glob < 5 mV in [2500, 8000] ns)\n\n";
    std::cout << "### A3.1: Mean of (wf - ped_glob) [mV] across time ranges (No Pulse)\n";
    std::cout << std::setw(15) << std::left << "Range [ns]";
    for (int pl = 0; pl < 4; ++pl) std::cout << std::setw(14) << plane_names[pl];
    std::cout << "\n" << std::string(71, '-') << "\n";
    for (int r = 0; r < 9; ++r) {
        std::cout << std::setw(15) << std::left << ranges[r].name;
        for (int pl = 0; pl < 4; ++pl) {
            double avg = (count_range_nopulse[pl][r] > 0) ? (sum_range_nopulse[pl][r] / count_range_nopulse[pl][r]) : 0.0;
            std::cout << std::setw(14) << std::fixed << std::setprecision(4) << avg;
        }
        std::cout << "\n";
    }

    std::cout << "\n### A3.2: b_shift = mean(wf in [3700, 3950]) - ped_glob [mV] (No Pulse)\n";
    std::cout << std::setw(14) << std::left << "Plane"
              << std::setw(14) << "Mean [mV]"
              << std::setw(18) << "RMS (events)"
              << std::setw(18) << "RMS (ch-events)"
              << std::setw(12) << "N_ch-events"
              << "\n";
    std::cout << std::string(76, '-') << "\n";
    for (int pl = 0; pl < 4; ++pl) {
        double m_ev = 0, rms_ev = 0;
        get_stats(b_shift_plane_nopulse[pl], m_ev, rms_ev);
        double m_ch = 0, rms_ch = 0;
        get_stats(b_shift_ch_nopulse[pl], m_ch, rms_ch);
        std::cout << std::setw(14) << std::left << plane_names[pl]
                  << std::setw(14) << std::fixed << std::setprecision(4) << m_ev
                  << std::setw(18) << std::fixed << std::setprecision(4) << rms_ev
                  << std::setw(18) << std::fixed << std::setprecision(4) << rms_ch
                  << std::setw(12) << b_shift_ch_nopulse[pl].size()
                  << "\n";
    }

    // Print Step C: Q pedestal noise
    std::cout << "\n## Step C: Q Pedestal Noise (ph_glob < 5 mV in [2500, 8000] ns)\n";
    std::cout << "Window: [t0 - 3, t0 + kNpost] samples at t0 = 4200 ns (sample 262)\n";
    std::cout << "ped_glob: t < 2500 ns | ped_loc: offsets -12..-4 from t0\n\n";

    // Per plane
    std::cout << "### Per Plane:\n";
    std::cout << std::setw(12) << std::left << "Plane"
              << std::setw(8) << "kNpost"
              << std::setw(10) << "N"
              << std::setw(14) << "Mean_Q(glob)"
              << std::setw(14) << "RMS_Q(glob)"
              << std::setw(14) << "Mean_Q(loc)"
              << std::setw(14) << "RMS_Q(loc)"
              << std::setw(14) << "RMS_loc/glob"
              << "\n";
    std::cout << std::string(100, '-') << "\n";

    for (int pl = 0; pl < 4; ++pl) {
        double mg = 0, rmsg = 0, ml = 0, rmsl = 0;
        get_stats(Q_glob_nopulse[pl], mg, rmsg);
        get_stats(Q_loc_nopulse[pl], ml, rmsl);
        int kNpost = (pl < 2) ? 8 : 22;
        double ratio = (rmsg > 0) ? (rmsl / rmsg) : 0.0;
        std::cout << std::setw(12) << std::left << plane_names[pl]
                  << std::setw(8) << kNpost
                  << std::setw(10) << Q_glob_nopulse[pl].size()
                  << std::setw(14) << std::fixed << std::setprecision(2) << mg
                  << std::setw(14) << std::fixed << std::setprecision(2) << rmsg
                  << std::setw(14) << std::fixed << std::setprecision(2) << ml
                  << std::setw(14) << std::fixed << std::setprecision(2) << rmsl
                  << std::setw(14) << std::fixed << std::setprecision(4) << ratio
                  << "\n";
    }

    // Also for hod2018 at kNpost = 12
    for (int pl = 0; pl < 2; ++pl) {
        double mg = 0, rmsg = 0, ml = 0, rmsl = 0;
        get_stats(Q_glob_nopulse_n12[pl], mg, rmsg);
        get_stats(Q_loc_nopulse_n12[pl], ml, rmsl);
        double ratio = (rmsg > 0) ? (rmsl / rmsg) : 0.0;
        std::cout << std::setw(12) << std::left << (std::string(plane_names[pl]) + " (n12)")
                  << std::setw(8) << 12
                  << std::setw(10) << Q_glob_nopulse_n12[pl].size()
                  << std::setw(14) << std::fixed << std::setprecision(2) << mg
                  << std::setw(14) << std::fixed << std::setprecision(2) << rmsg
                  << std::setw(14) << std::fixed << std::setprecision(2) << ml
                  << std::setw(14) << std::fixed << std::setprecision(2) << rmsl
                  << std::setw(14) << std::fixed << std::setprecision(4) << ratio
                  << "\n";
    }

    // Per station summary
    std::cout << "\n### Per Station (Combined):\n";
    std::cout << std::setw(12) << std::left << "Station"
              << std::setw(8) << "kNpost"
              << std::setw(10) << "N"
              << std::setw(14) << "Mean_Q(glob)"
              << std::setw(14) << "RMS_Q(glob)"
              << std::setw(14) << "Mean_Q(loc)"
              << std::setw(14) << "RMS_Q(loc)"
              << std::setw(14) << "RMS_loc/glob"
              << "\n";
    std::cout << std::string(100, '-') << "\n";

    // hod2018 combined (kNpost=8)
    std::vector<double> Q18_g, Q18_l;
    Q18_g.insert(Q18_g.end(), Q_glob_nopulse[0].begin(), Q_glob_nopulse[0].end());
    Q18_g.insert(Q18_g.end(), Q_glob_nopulse[1].begin(), Q_glob_nopulse[1].end());
    Q18_l.insert(Q18_l.end(), Q_loc_nopulse[0].begin(), Q_loc_nopulse[0].end());
    Q18_l.insert(Q18_l.end(), Q_loc_nopulse[1].begin(), Q_loc_nopulse[1].end());
    double m18_g = 0, rms18_g = 0, m18_l = 0, rms18_l = 0;
    get_stats(Q18_g, m18_g, rms18_g);
    get_stats(Q18_l, m18_l, rms18_l);
    std::cout << std::setw(12) << std::left << "hod2018"
              << std::setw(8) << 8
              << std::setw(10) << Q18_g.size()
              << std::setw(14) << std::fixed << std::setprecision(2) << m18_g
              << std::setw(14) << std::fixed << std::setprecision(2) << rms18_g
              << std::setw(14) << std::fixed << std::setprecision(2) << m18_l
              << std::setw(14) << std::fixed << std::setprecision(2) << rms18_l
              << std::setw(14) << std::fixed << std::setprecision(4) << (rms18_l / rms18_g)
              << "\n";

    // hod2018 combined (kNpost=12)
    std::vector<double> Q18_g12, Q18_l12;
    Q18_g12.insert(Q18_g12.end(), Q_glob_nopulse_n12[0].begin(), Q_glob_nopulse_n12[0].end());
    Q18_g12.insert(Q18_g12.end(), Q_glob_nopulse_n12[1].begin(), Q_glob_nopulse_n12[1].end());
    Q18_l12.insert(Q18_l12.end(), Q_loc_nopulse_n12[0].begin(), Q_loc_nopulse_n12[0].end());
    Q18_l12.insert(Q18_l12.end(), Q_loc_nopulse_n12[1].begin(), Q_loc_nopulse_n12[1].end());
    double m18_g12 = 0, rms18_g12 = 0, m18_l12 = 0, rms18_l12 = 0;
    get_stats(Q18_g12, m18_g12, rms18_g12);
    get_stats(Q18_l12, m18_l12, rms18_l12);
    std::cout << std::setw(12) << std::left << "hod2018"
              << std::setw(8) << 12
              << std::setw(10) << Q18_g12.size()
              << std::setw(14) << std::fixed << std::setprecision(2) << m18_g12
              << std::setw(14) << std::fixed << std::setprecision(2) << rms18_g12
              << std::setw(14) << std::fixed << std::setprecision(2) << m18_l12
              << std::setw(14) << std::fixed << std::setprecision(2) << rms18_l12
              << std::setw(14) << std::fixed << std::setprecision(4) << (rms18_l12 / rms18_g12)
              << "\n";

    // hod2019 combined (kNpost=22)
    std::vector<double> Q19_g, Q19_l;
    Q19_g.insert(Q19_g.end(), Q_glob_nopulse[2].begin(), Q_glob_nopulse[2].end());
    Q19_g.insert(Q19_g.end(), Q_glob_nopulse[3].begin(), Q_glob_nopulse[3].end());
    Q19_l.insert(Q19_l.end(), Q_loc_nopulse[2].begin(), Q_loc_nopulse[2].end());
    Q19_l.insert(Q19_l.end(), Q_loc_nopulse[3].begin(), Q_loc_nopulse[3].end());
    double m19_g = 0, rms19_g = 0, m19_l = 0, rms19_l = 0;
    get_stats(Q19_g, m19_g, rms19_g);
    get_stats(Q19_l, m19_l, rms19_l);
    std::cout << std::setw(12) << std::left << "hod2019"
              << std::setw(8) << 22
              << std::setw(10) << Q19_g.size()
              << std::setw(14) << std::fixed << std::setprecision(2) << m19_g
              << std::setw(14) << std::fixed << std::setprecision(2) << rms19_g
              << std::setw(14) << std::fixed << std::setprecision(2) << m19_l
              << std::setw(14) << std::fixed << std::setprecision(2) << rms19_l
              << std::setw(14) << std::fixed << std::setprecision(4) << (rms19_l / rms19_g)
              << "\n";
}

