#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

void hod2018_threshold_scan(const char* input_file = "lectura_vx1740_ext_nim_1M_v8.root") {
    TFile fin(input_file, "READ");
    if (fin.IsZombie()) {
        std::cerr << "Error opening " << input_file << "\n";
        return;
    }
    TTree* tree = (TTree*)fin.Get("data_dig");
    if (!tree) {
        std::cerr << "Error finding tree data_dig\n";
        return;
    }

    Long64_t nEnt = tree->GetEntries();
    std::cout << "Reading " << nEnt << " entries from " << input_file << "...\n";

    std::vector<float>* tw = nullptr;
    tree->SetBranchAddress("time", &tw);

    std::vector<float>* cw[32] = {};
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("time", 1);
    for (int c = 0; c < 32; ++c) {
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    // Max ph arrays per event for Method A (global ped) and Method B (local ped)
    // 0: 18X sig, 1: 18X ctl, 2: 18Y sig, 3: 18Y ctl
    std::vector<float> max_ph_A[4];
    std::vector<float> max_ph_B[4];
    for (int p = 0; p < 4; ++p) {
        max_ph_A[p].resize(nEnt, -1e9f);
        max_ph_B[p].resize(nEnt, -1e9f);
    }

    constexpr double kSigTlo = 4128.0;
    constexpr double kSigThi = 4272.0;
    constexpr double kCtlTlo = 6128.0;
    constexpr double kCtlThi = 6272.0;

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 30000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        for (int c = 0; c < 32; ++c) {
            if (!cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // Global pedestal (t < 2500 ns)
            double s_glob = 0.0;
            int n_glob = 0;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] < 2500.0) { s_glob += wf[i]; n_glob++; }
            }
            if (n_glob < 2) continue;
            const double ped_glob = s_glob / n_glob;

            // Signal window [4128, 4272]
            double vmin_s = 1e30; int imin_s = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_s) {
                    vmin_s = wf[i]; imin_s = i;
                }
            }

            // Control window [6128, 6272]
            double vmin_c = 1e30; int imin_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < vmin_c) {
                    vmin_c = wf[i]; imin_c = i;
                }
            }

            // Method A: ph with global pedestal
            if (imin_s >= 0) {
                float ph_s_A = ped_glob - vmin_s;
                int pl = (c < 16) ? 0 : 2; // 0 for 18X, 2 for 18Y
                if (ph_s_A > max_ph_A[pl][iev]) max_ph_A[pl][iev] = ph_s_A;
            }
            if (imin_c >= 0) {
                float ph_c_A = ped_glob - vmin_c;
                int pl = (c < 16) ? 1 : 3; // 1 for 18X ctl, 3 for 18Y ctl
                if (ph_c_A > max_ph_A[pl][iev]) max_ph_A[pl][iev] = ph_c_A;
            }

            // Method B: ph with local pedestal (mean of imin-12 .. imin-4)
            if (imin_s >= 0) {
                double ped_loc_s = ped_glob;
                if (imin_s >= 12 && imin_s - 4 < n_samp) {
                    double s_loc = 0.0;
                    for (int j = -12; j <= -4; ++j) s_loc += wf[imin_s + j];
                    ped_loc_s = s_loc / 9.0;
                }
                float ph_s_B = ped_loc_s - vmin_s;
                int pl = (c < 16) ? 0 : 2;
                if (ph_s_B > max_ph_B[pl][iev]) max_ph_B[pl][iev] = ph_s_B;
            }
            if (imin_c >= 0) {
                double ped_loc_c = ped_glob;
                if (imin_c >= 12 && imin_c - 4 < n_samp) {
                    double s_loc = 0.0;
                    for (int j = -12; j <= -4; ++j) s_loc += wf[imin_c + j];
                    ped_loc_c = s_loc / 9.0;
                }
                float ph_c_B = ped_loc_c - vmin_c;
                int pl = (c < 16) ? 1 : 3;
                if (ph_c_B > max_ph_B[pl][iev]) max_ph_B[pl][iev] = ph_c_B;
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n\n";

    const Long64_t N = (nEnt > 1) ? (nEnt - 1) : 0;
    const double N_d = (double)N;

    const std::vector<double> thresholds = {20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0};

    auto run_scan = [&](const std::string& title, const std::vector<float> max_ph[4]) {
        std::cout << "========================================================================================================================\n";
        std::cout << title << "\n";
        std::cout << "========================================================================================================================\n";
        std::cout << std::setw(6)  << "T"
                  << std::setw(10) << "N_X_S"
                  << std::setw(10) << "N_X_C"
                  << std::setw(10) << "N_Y_S"
                  << std::setw(10) << "N_Y_C"
                  << std::setw(12) << "N_XY_same"
                  << std::setw(12) << "N_XY_mixed"
                  << std::setw(12) << "aX"
                  << std::setw(12) << "aY"
                  << std::setw(12) << "rX"
                  << std::setw(12) << "rY"
                  << std::setw(14) << "N_XY_real"
                  << std::setw(26) << "N_XY_real / (0.0360 * N)"
                  << "\n";
        std::cout << "------------------------------------------------------------------------------------------------------------------------\n";

        for (double T : thresholds) {
            long long N_X_S = 0, N_X_C = 0, N_Y_S = 0, N_Y_C = 0;
            long long N_XY_same = 0, N_XY_mixed = 0;

            for (Long64_t i = 0; i < N; ++i) {
                Long64_t next_i = i + 1;
                bool xs = (max_ph[0][i] > T);
                bool xc = (max_ph[1][i] > T);
                bool ys = (max_ph[2][i] > T);
                bool yc = (max_ph[3][i] > T);
                bool ys_next = (max_ph[2][next_i] > T);

                if (xs) N_X_S++;
                if (xc) N_X_C++;
                if (ys) N_Y_S++;
                if (yc) N_Y_C++;
                if (xs && ys) N_XY_same++;
                if (xs && ys_next) N_XY_mixed++;
            }

            double aX = (double)N_X_C / N_d;
            double aY = (double)N_Y_C / N_d;
            double rX = (1.0 - aX > 1e-12) ? (((double)N_X_S / N_d - aX) / (1.0 - aX)) : 0.0;
            double rY = (1.0 - aY > 1e-12) ? (((double)N_Y_S / N_d - aY) / (1.0 - aY)) : 0.0;
            double denom = (1.0 - aX) * (1.0 - aY);
            double N_XY_real = (std::abs(denom) > 1e-12) ? (((double)(N_XY_same - N_XY_mixed)) / denom + rX * rY * N_d) : 0.0;
            double eff = N_XY_real / (0.0360 * N_d);

            std::cout << std::fixed
                      << std::setw(6)  << std::setprecision(0) << T
                      << std::setw(10) << N_X_S
                      << std::setw(10) << N_X_C
                      << std::setw(10) << N_Y_S
                      << std::setw(10) << N_Y_C
                      << std::setw(12) << N_XY_same
                      << std::setw(12) << N_XY_mixed
                      << std::setw(12) << std::setprecision(6) << aX
                      << std::setw(12) << std::setprecision(6) << aY
                      << std::setw(12) << std::setprecision(6) << rX
                      << std::setw(12) << std::setprecision(6) << rY
                      << std::setw(14) << std::setprecision(2) << N_XY_real
                      << std::setw(26) << std::setprecision(6) << eff
                      << "\n";
        }
        std::cout << "========================================================================================================================\n\n";
    };

    std::cout << "Formulas of reconstruct_corrected:\n";
    std::cout << "  N = N_events - 1 = " << N << "\n";
    std::cout << "  aX = N_X_C / N\n";
    std::cout << "  aY = N_Y_C / N\n";
    std::cout << "  rX = (N_X_S / N - aX) / (1 - aX)\n";
    std::cout << "  rY = (N_Y_S / N - aY) / (1 - aY)\n";
    std::cout << "  N_XY_real = (N_XY_same - N_XY_mixed) / ((1 - aX) * (1 - aY)) + rX * rY * N\n";
    std::cout << "  N_XY_real / (0.0360 * N)\n\n";

    run_scan("Method (A): Global pedestal (t < 2500 ns)", max_ph_A);
    run_scan("Method (B): Local pedestal (mean of samples imin-12 .. imin-4)", max_ph_B);
}

