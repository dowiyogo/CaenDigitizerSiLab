#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

struct RunInfo {
    std::string name;
    std::string filename;
};

void dark_pulse_rates_step5() {
    std::vector<RunInfo> runs = {
        {"v1",  "lectura_vx1740_ext_nim_1M.root"},
        {"v2",  "lectura_vx1740_ext_nim_1M_v2.root"},
        {"v3",  "lectura_vx1740_ext_nim_1M_v3.root"},
        {"v4",  "lectura_vx1740_ext_nim_1M_v4.root"},
        {"v5",  "lectura_vx1740_ext_nim_1M_v5.root"},
        {"v6",  "lectura_vx1740_ext_nim_1M_v6.root"},
        {"v7",  "lectura_vx1740_ext_nim_1M_v7.root"},
        {"v8",  "lectura_vx1740_ext_nim_1M_v8.root"},
        {"v9",  "lectura_vx1740_ext_nim_1M_v9.root"},
        {"v10", "lectura_vx1740_ext_nim_1M_v10.root"}
    };

    const double thrs[5] = {5.0, 10.0, 20.0, 35.0, 50.0};
    constexpr double kCtlTlo = 6128.0;
    constexpr double kCtlThi = 6272.0;

    std::cout << "========================================================================================\n";
    std::cout << "Step 5: Per-Run Dark-Pulse Rates (Runs v1..v10)\n";
    std::cout << "Control window: [6128, 6272] ns (144 ns), global pedestal (t < 2500 ns)\n";
    std::cout << "========================================================================================\n\n";

    std::cout << std::left << std::setw(8) << "run"
              << std::setw(8) << "plane"
              << std::right << std::setw(12) << "rate(5)"
              << std::setw(12) << "rate(10)"
              << std::setw(12) << "rate(20)"
              << std::setw(12) << "rate(35)"
              << std::setw(12) << "rate(50)"
              << "\n";
    std::cout << std::left << std::setw(8) << ""
              << std::setw(8) << ""
              << std::right << std::setw(12) << "[kHz]"
              << std::setw(12) << "[kHz]"
              << std::setw(12) << "[kHz]"
              << std::setw(12) << "[kHz]"
              << std::setw(12) << "[kHz]"
              << "\n";
    std::cout << std::string(76, '-') << "\n";

    for (const auto& r : runs) {
        TFile fin(r.filename.c_str(), "READ");
        if (fin.IsZombie()) {
            std::cerr << "Error opening " << r.filename << "\n";
            continue;
        }
        TTree* tree = (TTree*)fin.Get("data_dig");
        if (!tree) {
            std::cerr << "Error finding tree data_dig in " << r.filename << "\n";
            continue;
        }

        Long64_t nEnt = tree->GetEntries();

        std::vector<float>* tw = nullptr;
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("time", 1);
        tree->SetBranchAddress("time", &tw);

        std::vector<float>* cw[64] = {};
        for (int c = 0; c < 64; ++c) {
            if (c == 3) continue; // skip dead ch3
            tree->SetBranchStatus(Form("Ch%d", c), 1);
            tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
        }
        tree->SetCacheSize(256LL << 20);

        // hits[plane][thr_idx]
        // Plane 0: 18X (ch 0..15, 15 live channels)
        // Plane 1: 18Y (ch 16..31, 16 live channels)
        // Plane 2: 19Y (ch 32..47, 16 live channels)
        // Plane 3: 19X (ch 48..63, 16 live channels)
        long long hits[4][5] = {};

        // Find sample indices for t < 2500 and [6128, 6272]
        tree->GetEntry(0);
        if (!tw || tw->empty()) {
            std::cerr << "Error: empty time vector in " << r.filename << "\n";
            continue;
        }
        const auto& tv0 = *tw;
        int n_samp = (int)tv0.size();

        int i_ped_end = 0;
        while (i_ped_end < n_samp && tv0[i_ped_end] < 2500.0) i_ped_end++;

        int i_ctl_start = -1, i_ctl_end = -1;
        for (int i = 0; i < n_samp; ++i) {
            if (tv0[i] >= kCtlTlo && tv0[i] <= kCtlThi) {
                if (i_ctl_start < 0) i_ctl_start = i;
                i_ctl_end = i;
            }
        }

        for (Long64_t iev = 0; iev < nEnt; ++iev) {
            tree->GetEntry(iev);

            for (int c = 0; c < 64; ++c) {
                if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
                const auto& wf = *cw[c];

                // Global pedestal
                double sum_ped = 0;
                for (int i = 0; i < i_ped_end; ++i) sum_ped += wf[i];
                double ped_glob = (i_ped_end > 0) ? (sum_ped / i_ped_end) : 0.0;

                // Min in control window
                double vmin = 1e30;
                for (int i = i_ctl_start; i <= i_ctl_end; ++i) {
                    if (wf[i] < vmin) vmin = wf[i];
                }
                double ph = ped_glob - vmin;

                int pl = -1;
                if (c >= 0 && c <= 15) pl = 0;
                else if (c >= 16 && c <= 31) pl = 1;
                else if (c >= 32 && c <= 47) pl = 2;
                else if (c >= 48 && c <= 63) pl = 3;

                for (int t = 0; t < 5; ++t) {
                    if (ph > thrs[t]) {
                        hits[pl][t]++;
                    }
                }
            }
        }

        // Output rows for this run
        const char* pl_names[4] = {"18X", "18Y", "19Y", "19X"};
        const int n_live_ch[4] = {15, 16, 16, 16};
        double norm_ms = (double)nEnt * 144.0e-6; // ms for 1 channel

        for (int pl = 0; pl < 4; ++pl) {
            double ch_norm_ms = norm_ms * (double)n_live_ch[pl];
            std::cout << std::left << std::setw(8) << r.name
                      << std::setw(8) << pl_names[pl];
            for (int t = 0; t < 5; ++t) {
                double rate_khz = (ch_norm_ms > 0) ? ((double)hits[pl][t] / ch_norm_ms) : 0.0;
                std::cout << std::right << std::fixed << std::setprecision(3) << std::setw(12) << rate_khz;
            }
            std::cout << "\n";
        }
    }
    std::cout << std::string(76, '-') << "\n";
}

