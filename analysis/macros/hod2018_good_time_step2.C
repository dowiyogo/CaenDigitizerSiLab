#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <set>

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::string channels;
    std::string reason;
};

struct EventSummary {
    Long_t timestamp;
    std::string dh_min;  // YYYY-MM-DD HH:MM
    std::string dh_hour; // YYYY-MM-DD HH:00
    bool xs;
    bool xc;
    bool ys;
    bool yc;
};

struct RowResult {
    std::string label;
    long long n_used;
    int live_hours;
    long long N_X_S;
    long long N_X_C;
    long long N_Y_S;
    long long N_Y_C;
    long long N_XY_same;
    long long N_XY_mixed;
    double aX;
    double aY;
    double rX;
    double rY;
    double N_XY_real;
    double eff;
};

RowResult compute_row(const std::string& label,
                      const std::vector<EventSummary>& kept_events,
                      const std::set<std::string>& active_hour_bins) {
    RowResult res;
    res.label = label;
    res.n_used = (long long)kept_events.size();
    res.live_hours = (int)active_hour_bins.size();

    Long64_t N_pairs = (res.n_used > 1) ? (res.n_used - 1) : 0;
    double N_d = (double)N_pairs;

    res.N_X_S = 0; res.N_X_C = 0; res.N_Y_S = 0; res.N_Y_C = 0;
    res.N_XY_same = 0; res.N_XY_mixed = 0;

    for (Long64_t k = 0; k < N_pairs; ++k) {
        const auto& cur = kept_events[k];
        const auto& nxt = kept_events[k + 1];

        if (cur.xs) res.N_X_S++;
        if (cur.xc) res.N_X_C++;
        if (cur.ys) res.N_Y_S++;
        if (cur.yc) res.N_Y_C++;
        if (cur.xs && cur.ys) res.N_XY_same++;
        if (cur.xs && nxt.ys) res.N_XY_mixed++;
    }

    res.aX = (N_d > 0) ? ((double)res.N_X_C / N_d) : 0.0;
    res.aY = (N_d > 0) ? ((double)res.N_Y_C / N_d) : 0.0;
    res.rX = (1.0 - res.aX > 1e-12) ? (((double)res.N_X_S / N_d - res.aX) / (1.0 - res.aX)) : 0.0;
    res.rY = (1.0 - res.aY > 1e-12) ? (((double)res.N_Y_S / N_d - res.aY) / (1.0 - res.aY)) : 0.0;
    double denom = (1.0 - res.aX) * (1.0 - res.aY);
    res.N_XY_real = (std::abs(denom) > 1e-12) ? (((double)(res.N_XY_same - res.N_XY_mixed)) / denom + res.rX * res.rY * N_d) : 0.0;
    res.eff = (res.n_used > 0) ? (res.N_XY_real / (0.0360 * (double)res.n_used)) : 0.0;

    return res;
}

void print_row(const RowResult& r) {
    std::cout << std::left << std::setw(14) << r.label
              << std::right << std::setw(15) << r.n_used
              << std::setw(12) << r.live_hours
              << std::setw(10) << r.N_X_S
              << std::setw(10) << r.N_X_C
              << std::setw(10) << r.N_Y_S
              << std::setw(10) << r.N_Y_C
              << std::setw(12) << r.N_XY_same
              << std::setw(12) << r.N_XY_mixed
              << std::fixed << std::setprecision(6)
              << std::setw(10) << r.aX
              << std::setw(10) << r.aY
              << std::setw(10) << r.rX
              << std::setw(10) << r.rY
              << std::setprecision(2)
              << std::setw(14) << r.N_XY_real
              << std::setprecision(6)
              << std::setw(34) << r.eff
              << "\n";
}

void hod2018_good_time_step2() {
    // 1. Read calibration parameters T_c(5) from outputs/pe_calibration_v8.txt
    double T5[32] = {};
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
                T5[ch] = p0 + 3.5 * G;
            }
        }
    }
    pf.close();
    T5[3] = 20.0; // dead channel

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

    constexpr double kSigTlo = 4128.0, kSigThi = 4272.0;
    constexpr double kCtlTlo = 6128.0, kCtlThi = 6272.0;

    struct FileInfo {
        std::string run_name;
        std::string filename;
    };
    std::vector<FileInfo> files = {
        {"v7",  "lectura_vx1740_ext_nim_1M_v7.root"},
        {"v8",  "lectura_vx1740_ext_nim_1M_v8.root"},
        {"v9",  "lectura_vx1740_ext_nim_1M_v9.root"},
        {"v10", "lectura_vx1740_ext_nim_1M_v10.root"}
    };

    std::vector<RowResult> results;
    RowResult row_v8_hv_off;

    for (const auto& fi : files) {
        TFile fin(fi.filename.c_str(), "READ");
        if (fin.IsZombie()) {
            std::cerr << "Error opening " << fi.filename << "\n";
            continue;
        }
        TTree* tree = (TTree*)fin.Get("data_dig");
        if (!tree) {
            std::cerr << "Error finding tree data_dig in " << fi.filename << "\n";
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

        std::vector<EventSummary> events_good;
        std::vector<EventSummary> events_hv_off;
        std::set<std::string> hours_good;
        std::set<std::string> hours_hv_off;

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

            char buf_hr[32];
            std::snprintf(buf_hr, sizeof(buf_hr), "%04d-%02d-%02d %02d:00",
                          tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday, tm_val.tm_hour);
            std::string s_hr(buf_hr);

            bool is_exclude = false;
            bool is_hv_off = false;
            for (const auto& iv : intervals) {
                if (iv.run == fi.run_name) {
                    if (s_min >= iv.start_local && s_min <= iv.end_local) {
                        if (iv.status == "EXCLUDE") is_exclude = true;
                        if (iv.status == "HV_OFF")  is_hv_off = true;
                    }
                }
            }

            // If neither EXCLUDE nor HV_OFF -> good time
            bool keep_good = (!is_exclude && !is_hv_off);
            // If HV_OFF (and not EXCLUDE) -> hv_off interval
            bool keep_hv_off = (is_hv_off && !is_exclude);

            if (!keep_good && !keep_hv_off) continue;

            const auto& tv = *tw;
            const int n_samp = (int)tv.size();

            bool xs = false, xc = false, ys = false, yc = false;

            for (int c = 0; c < 32; ++c) {
                if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
                const auto& wf = *cw[c];

                // Signal window W2 [4128, 4272] ns
                double vmin_s = 1e30; int imin_s = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_s) {
                        vmin_s = wf[i]; imin_s = i;
                    }
                }
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
                    if (ped - vmin_s > T5[c]) {
                        if (c < 16) xs = true;
                        else        ys = true;
                    }
                }

                // Control window [6128, 6272] ns
                double vmin_c = 1e30; int imin_c = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < vmin_c) {
                        vmin_c = wf[i]; imin_c = i;
                    }
                }
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
                    if (ped - vmin_c > T5[c]) {
                        if (c < 16) xc = true;
                        else        yc = true;
                    }
                }
            }

            EventSummary es;
            es.timestamp = timestamp;
            es.dh_min = s_min;
            es.dh_hour = s_hr;
            es.xs = xs;
            es.xc = xc;
            es.ys = ys;
            es.yc = yc;

            if (keep_good) {
                events_good.push_back(es);
                hours_good.insert(s_hr);
            }
            if (keep_hv_off) {
                events_hv_off.push_back(es);
                hours_hv_off.insert(s_hr);
            }
        }

        if (fi.run_name == "v7") {
            results.push_back(compute_row("v7", events_good, hours_good));
        } else if (fi.run_name == "v8") {
            results.push_back(compute_row("v8", events_good, hours_good));
            row_v8_hv_off = compute_row("v8 (HV_OFF)", events_hv_off, hours_hv_off);
        } else if (fi.run_name == "v9") {
            results.push_back(compute_row("v9", events_good, hours_good));
        } else if (fi.run_name == "v10") {
            results.push_back(compute_row("v10", events_good, hours_good));
            results.push_back(compute_row("v10 (HV_OFF)", events_hv_off, hours_hv_off));
        }
    }

    results.push_back(row_v8_hv_off);

    std::cout << "\n==========================================================================================================================================\n";
    std::cout << "Step 2: hod2018 X-Y at 5 p.e. in Good Time Only (and HV_OFF intervals)\n";
    std::cout << "Threshold source: outputs/pe_calibration_v8.txt, T_c(5) = p0_c + 3.5 * G_c\n";
    std::cout << "Mixed partner = next event that is also kept (consecutive kept events k and k+1)\n";
    std::cout << "==========================================================================================================================================\n";
    std::cout << std::left << std::setw(14) << "run"
              << std::right << std::setw(15) << "N_events_used"
              << std::setw(12) << "live_hours"
              << std::setw(10) << "N_X_S"
              << std::setw(10) << "N_X_C"
              << std::setw(10) << "N_Y_S"
              << std::setw(10) << "N_Y_C"
              << std::setw(12) << "N_XY_same"
              << std::setw(12) << "N_XY_mixed"
              << std::setw(10) << "aX"
              << std::setw(10) << "aY"
              << std::setw(10) << "rX"
              << std::setw(10) << "rY"
              << std::setw(14) << "N_XY_real"
              << std::setw(34) << "N_XY_real / (0.0360 * N_events_used)"
              << "\n";
    std::cout << std::string(178, '-') << "\n";

    for (const auto& r : results) {
        print_row(r);
    }
    std::cout << std::string(178, '-') << "\n";
}

