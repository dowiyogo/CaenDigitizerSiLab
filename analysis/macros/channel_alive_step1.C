#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <ctime>
#include <algorithm>

struct HourChannelData {
    std::string dh;
    int n_events = 0;
    int hits_ctl[64] = {};
    int hits_w2[64] = {};
    double rate_ctl[64] = {};
    double rate_w2[64] = {};
};

struct DeadInterval {
    int channel;
    std::string first_hour;
    std::string last_hour;
    int num_hours;
    double median_rate;
};

void run_channel_alive_for_run(const std::string& run_name, const std::string& filename,
                               const double* T_thr, bool is_v10,
                               std::vector<HourChannelData>& out_hourly) {
    TFile fin(filename.c_str(), "READ");
    if (fin.IsZombie()) {
        std::cerr << "Error opening " << filename << "\n";
        return;
    }
    TTree* tree = (TTree*)fin.Get("data_dig");
    if (!tree) {
        std::cerr << "Error finding tree data_dig in " << filename << "\n";
        return;
    }

    Long64_t nEnt = tree->GetEntries();
    std::cout << "\n========================================================================================\n";
    std::cout << "Step 1: Processing Run " << run_name << " (" << filename << ", " << nEnt << " entries)\n";
    std::cout << "========================================================================================\n";

    std::vector<float>* tw = nullptr;
    Long_t timestamp = 0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("time", 1);
    tree->SetBranchAddress("time", &tw);
    tree->SetBranchStatus("timestamp", 1);
    tree->SetBranchAddress("timestamp", &timestamp);

    std::vector<float>* cw[64] = {};
    for (int c = 0; c < 64; ++c) {
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    constexpr double kCtlTlo = 6128.0, kCtlThi = 6272.0;
    constexpr double kSigTlo = 4128.0, kSigThi = 4272.0;
    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // -04:00 offset

    std::map<std::string, HourChannelData> hr_map;
    std::vector<std::string> hr_order;

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 50000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;

        double t_sec = boot_ref_epoch + ((double)timestamp) * 1e-6 + tz_offset;
        time_t t_local = (time_t)t_sec;
        struct tm tm_val; gmtime_r(&t_local, &tm_val);
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:00",
                      tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday, tm_val.tm_hour);
        std::string key(buf);

        if (hr_map.find(key) == hr_map.end()) {
            HourChannelData hd;
            hd.dh = key;
            hr_map[key] = hd;
            hr_order.push_back(key);
        }

        HourChannelData& hd = hr_map[key];
        hd.n_events++;

        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        for (int c = 0; c < 64; ++c) {
            if (!cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // In-time W2 [4128, 4272] ns
            double vmin_s = 1e30; int imin_s = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_s) {
                    vmin_s = wf[i]; imin_s = i;
                }
            }
            if (imin_s >= 0) {
                double ped = 0;
                if (c < 32 && imin_s >= 12 && imin_s - 4 < n_samp) {
                    double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_s + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                if (ped - vmin_s > T_thr[c]) hd.hits_w2[c]++;
            }

            // Out-of-time control window [6128, 6272] ns
            double vmin_c = 1e30; int imin_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < vmin_c) {
                    vmin_c = wf[i]; imin_c = i;
                }
            }
            if (imin_c >= 0) {
                double ped = 0;
                if (c < 32 && imin_c >= 12 && imin_c - 4 < n_samp) {
                    double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_c + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                if (ped - vmin_c > T_thr[c]) hd.hits_ctl[c]++;
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n";

    // Compute rates per hour
    out_hourly.clear();
    for (const auto& k : hr_order) {
        HourChannelData& hd = hr_map[k];
        double norm_ms = (double)hd.n_events * 144.0e-6; // 144 ns in ms
        for (int c = 0; c < 64; ++c) {
            hd.rate_ctl[c] = (norm_ms > 0) ? ((double)hd.hits_ctl[c] / norm_ms) : 0.0;
            hd.rate_w2[c]  = (norm_ms > 0) ? ((double)hd.hits_w2[c] / norm_ms) : 0.0;
        }
        out_hourly.push_back(hd);
    }

    // 1a: Median hourly rate per channel and run
    // For hod2018 (ch 0..31): out-of-time pulse rate in control window
    // For hod2019 (ch 32..63): in-time W2 hit rate/count (and control window rate)
    double med_rate[64] = {};
    for (int c = 0; c < 64; ++c) {
        std::vector<double> vals;
        for (const auto& hd : out_hourly) {
            // For hod2018, primary rate is control-window out-of-time rate; for hod2019, in-time W2 rate
            double v = (c < 32) ? hd.rate_ctl[c] : hd.rate_w2[c];
            vals.push_back(v);
        }
        std::sort(vals.begin(), vals.end());
        med_rate[c] = vals[vals.size() / 2];
        if (med_rate[c] == 0.0) {
            // If zero for > 50% of the run, check median among active hours
            std::vector<double> vpos;
            for (double v : vals) if (v > 0) vpos.push_back(v);
            if (!vpos.empty()) med_rate[c] = vpos[vpos.size() / 2];
        }
    }

    // Identify DEAD intervals for 1a:
    // A channel-hour is DEAD if its rate is < 10% of that median.
    std::vector<DeadInterval> dead_intervals;
    int n_hr = (int)out_hourly.size();

    for (int c = 0; c < 64; ++c) {
        double med = med_rate[c];
        if (med <= 0.0) {
            // Channel never had any activity in this run
            DeadInterval di;
            di.channel = c;
            di.first_hour = out_hourly[0].dh;
            di.last_hour = out_hourly.back().dh;
            di.num_hours = n_hr;
            di.median_rate = 0.0;
            dead_intervals.push_back(di);
            continue;
        }
        double thr = 0.10 * med;
        int start_idx = -1;
        for (int i = 0; i <= n_hr; ++i) {
            bool is_dead = false;
            if (i < n_hr) {
                double r_val = (c < 32) ? out_hourly[i].rate_ctl[c] : out_hourly[i].rate_w2[c];
                if (r_val < thr) is_dead = true;
            }
            if (is_dead) {
                if (start_idx < 0) start_idx = i;
            } else {
                if (start_idx >= 0) {
                    int end_idx = i - 1;
                    DeadInterval di;
                    di.channel = c;
                    di.first_hour = out_hourly[start_idx].dh;
                    di.last_hour = out_hourly[end_idx].dh;
                    di.num_hours = end_idx - start_idx + 1;
                    di.median_rate = med;
                    dead_intervals.push_back(di);
                    start_idx = -1;
                }
            }
        }
    }

    // Print 1a table of DEAD intervals
    std::cout << "\nTable 1a: DEAD Intervals for Run " << run_name << " (" << n_hr << " hours)\n";
    std::cout << std::string(88, '-') << "\n";
    std::cout << std::left << std::setw(8) << "Channel"
              << std::setw(8) << "Plane"
              << std::setw(20) << "First Hour"
              << std::setw(20) << "Last Hour"
              << std::right << std::setw(12) << "Num Hours"
              << std::setw(18) << "Median Rate[kHz]"
              << "\n";
    std::cout << std::string(88, '-') << "\n";
    for (const auto& di : dead_intervals) {
        const char* pl = (di.channel < 16) ? "18X" : (di.channel < 32) ? "18Y" : (di.channel < 48) ? "19Y" : "19X";
        std::cout << std::left << std::setw(8) << di.channel
                  << std::setw(8) << pl
                  << std::setw(20) << di.first_hour
                  << std::setw(20) << di.last_hour
                  << std::right << std::setw(12) << di.num_hours
                  << std::fixed << std::setprecision(2) << std::setw(18) << di.median_rate
                  << "\n";
    }
    std::cout << std::string(88, '-') << "\n";
}

void channel_alive_step1() {
    double T_thr[64] = {};
    // Load hod2018 calibration (above 5 p.e.)
    std::ifstream pf("outputs/pe_calibration_v8.txt");
    std::string line;
    while (std::getline(pf, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int ch; std::string pl; double p0, G, s0, s1, chi2; int ndf, status;
        if (iss >> ch >> pl >> p0 >> G >> s0 >> s1 >> chi2 >> ndf >> status) {
            if (ch >= 0 && ch < 32) T_thr[ch] = p0 + 3.5 * G;
        }
    }
    pf.close();
    T_thr[3] = 20.0; // dead channel
    for (int c = 32; c < 64; ++c) T_thr[c] = 5.0; // hod2019 threshold: above 5 mV

    std::vector<HourChannelData> hourly_v7, hourly_v8, hourly_v9, hourly_v10;
    run_channel_alive_for_run("v7", "lectura_vx1740_ext_nim_1M_v7.root", T_thr, false, hourly_v7);
    run_channel_alive_for_run("v8", "lectura_vx1740_ext_nim_1M_v8.root", T_thr, false, hourly_v8);
    run_channel_alive_for_run("v9", "lectura_vx1740_ext_nim_1M_v9.root", T_thr, false, hourly_v9);
    run_channel_alive_for_run("v10", "lectura_vx1740_ext_nim_1M_v10.root", T_thr, true, hourly_v10);

    // 1b: Print the hourly table of the number of live channels per plane (18X, 18Y, 19Y, 19X) for v10
    // Compute medians for v10
    double med_v10[64] = {};
    for (int c = 0; c < 64; ++c) {
        std::vector<double> vals;
        for (const auto& hd : hourly_v10) {
            double v = (c < 32) ? hd.rate_ctl[c] : hd.rate_w2[c];
            vals.push_back(v);
        }
        std::sort(vals.begin(), vals.end());
        med_v10[c] = vals[vals.size() / 2];
        if (med_v10[c] == 0.0) {
            std::vector<double> vpos;
            for (double v : vals) if (v > 0) vpos.push_back(v);
            if (!vpos.empty()) med_v10[c] = vpos[vpos.size() / 2];
        }
    }

    std::cout << "\n========================================================================================\n";
    std::cout << "Table 1b: Number of Live Channels per Plane versus Time for Run v10 (120 hours)\n";
    std::cout << "========================================================================================\n";
    std::cout << std::left << std::setw(18) << "date-hour"
              << std::right << std::setw(10) << "N_events"
              << std::setw(8) << "18X"
              << std::setw(8) << "18Y"
              << std::setw(8) << "19Y"
              << std::setw(8) << "19X"
              << std::setw(10) << "Total"
              << "\n";
    std::cout << std::string(70, '-') << "\n";

    for (const auto& hd : hourly_v10) {
        int live_cnt[4] = {}; // 0:18X, 1:18Y, 2:19Y, 3:19X
        for (int c = 0; c < 64; ++c) {
            if (med_v10[c] <= 0.0) continue;
            double thr = 0.10 * med_v10[c];
            double r_val = (c < 32) ? hd.rate_ctl[c] : hd.rate_w2[c];
            if (r_val >= thr) {
                int pl = (c < 16) ? 0 : (c < 32) ? 1 : (c < 48) ? 2 : 3;
                live_cnt[pl]++;
            }
        }
        int tot = live_cnt[0] + live_cnt[1] + live_cnt[2] + live_cnt[3];
        std::cout << std::left << std::setw(18) << hd.dh
                  << std::right << std::setw(10) << hd.n_events
                  << std::setw(8) << live_cnt[0]
                  << std::setw(8) << live_cnt[1]
                  << std::setw(8) << live_cnt[2]
                  << std::setw(8) << live_cnt[3]
                  << std::setw(10) << tot
                  << "\n";
    }
    std::cout << std::string(70, '-') << "\n";
}

