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

struct HourBin {
    std::string date_hour;
    Long64_t n_events = 0;
    Long_t first_ts = 0;
    Long_t last_ts = 0;
    long long hits[32] = {};
};

void run_time_dependence_for_file(const char* input_file, const char* run_name, const double* T5) {
    TFile fin(input_file, "READ");
    if (fin.IsZombie()) {
        std::cerr << "Error opening " << input_file << "\n";
        return;
    }
    TTree* tree = (TTree*)fin.Get("data_dig");
    if (!tree) {
        std::cerr << "Error finding tree data_dig in " << input_file << "\n";
        return;
    }

    Long64_t nEnt = tree->GetEntries();
    std::cout << "\n========================================================================================\n";
    std::cout << "Step 4: Time Dependence for Run " << run_name << " (" << nEnt << " entries)\n";
    std::cout << "========================================================================================\n";

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

    constexpr double kCtlTlo = 6128.0;
    constexpr double kCtlThi = 6272.0;
    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // -04:00 offset

    std::map<std::string, HourBin> bins;
    std::vector<std::string> bin_order;

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 50000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;

        double t_sec = boot_ref_epoch + ((double)timestamp) * 1e-6;
        double t_local = t_sec + tz_offset;
        time_t t_local_epoch = (time_t)t_local;
        struct tm tm_val;
        gmtime_r(&t_local_epoch, &tm_val);

        char dh_buf[32];
        std::snprintf(dh_buf, sizeof(dh_buf), "%04d-%02d-%02d %02d:00",
                      tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday, tm_val.tm_hour);
        std::string key(dh_buf);

        if (bins.find(key) == bins.end()) {
            HourBin hb;
            hb.date_hour = key;
            hb.first_ts = timestamp;
            bins[key] = hb;
            bin_order.push_back(key);
        }

        HourBin& b = bins[key];
        b.n_events++;
        b.last_ts = timestamp;

        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        for (int c = 0; c < 32; ++c) {
            if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            double vmin = 1e30;
            int imin = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < vmin) {
                    vmin = wf[i];
                    imin = i;
                }
            }

            if (imin >= 0) {
                double ped_loc = 0;
                if (imin >= 12 && imin - 4 < n_samp) {
                    double s = 0;
                    for (int j = -12; j <= -4; ++j) s += wf[imin + j];
                    ped_loc = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped_loc = (n > 0) ? s / n : 0;
                }
                double ph = ped_loc - vmin;
                if (ph > T5[c]) {
                    b.hits[c]++;
                }
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n";

    // Print table
    std::cout << "\nTable: 1-Hour Binned Rates for Run " << run_name << "\n";
    std::cout << std::string(115, '-') << "\n";
    std::cout << std::left << std::setw(18) << "date-hour"
              << std::right << std::setw(10) << "N_events"
              << std::setw(14) << "trig_rate[Hz]"
              << std::setw(11) << "18X_mean"
              << std::setw(11) << "18Y_mean"
              << std::setw(9) << "ch2"
              << std::setw(9) << "ch9"
              << std::setw(9) << "ch15"
              << std::setw(9) << "ch16"
              << std::setw(9) << "ch31"
              << std::setw(9) << "ch12"
              << "\n";
    std::cout << std::left << std::setw(18) << ""
              << std::right << std::setw(10) << ""
              << std::setw(14) << ""
              << std::setw(11) << "[kHz/ch]"
              << std::setw(11) << "[kHz/ch]"
              << std::setw(9) << "[kHz]"
              << std::setw(9) << "[kHz]"
              << std::setw(9) << "[kHz]"
              << std::setw(9) << "[kHz]"
              << std::setw(9) << "[kHz]"
              << std::setw(9) << "[kHz]"
              << "\n";
    std::cout << std::string(115, '-') << "\n";

    for (const auto& k : bin_order) {
        const auto& b = bins[k];
        double dt_s = ((double)(b.last_ts - b.first_ts)) * 1e-6;
        double trig_rate = (dt_s > 0) ? ((double)b.n_events / dt_s) : 0.0;
        double norm_ch_ms = (double)b.n_events * 144.0e-6; // ms for 1 channel

        long long hits_18X = 0;
        for (int c = 0; c < 16; ++c) if (c != 3) hits_18X += b.hits[c];
        double rate_18X = (norm_ch_ms > 0) ? ((double)hits_18X / (15.0 * norm_ch_ms)) : 0.0;

        long long hits_18Y = 0;
        for (int c = 16; c < 32; ++c) hits_18Y += b.hits[c];
        double rate_18Y = (norm_ch_ms > 0) ? ((double)hits_18Y / (16.0 * norm_ch_ms)) : 0.0;

        auto ch_rate = [&](int c) {
            return (norm_ch_ms > 0) ? ((double)b.hits[c] / norm_ch_ms) : 0.0;
        };

        std::cout << std::left << std::setw(18) << b.date_hour
                  << std::right << std::setw(10) << b.n_events
                  << std::fixed << std::setprecision(4)
                  << std::setw(14) << trig_rate
                  << std::setprecision(3)
                  << std::setw(11) << rate_18X
                  << std::setw(11) << rate_18Y
                  << std::setw(9) << ch_rate(2)
                  << std::setw(9) << ch_rate(9)
                  << std::setw(9) << ch_rate(15)
                  << std::setw(9) << ch_rate(16)
                  << std::setw(9) << ch_rate(31)
                  << std::setw(9) << ch_rate(12)
                  << "\n";
    }
    std::cout << std::string(115, '-') << "\n";
}

void time_dependence_step4() {
    // Read calibration parameters from outputs/pe_calibration_v8.txt
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

    run_time_dependence_for_file("lectura_vx1740_ext_nim_1M_v8.root", "v8", T5);
    run_time_dependence_for_file("lectura_vx1740_ext_nim_1M_v10.root", "v10", T5);
}
