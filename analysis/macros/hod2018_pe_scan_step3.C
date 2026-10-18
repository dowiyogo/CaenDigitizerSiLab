#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

void hod2018_pe_scan_step3(const char* input_file = "lectura_vx1740_ext_nim_1M_v8.root") {
    // 1. Read calibration parameters from outputs/pe_calibration_v8.txt
    double p0_cal[32] = {};
    double G_cal[32] = {};
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
            }
        }
    }
    pf.close();

    // 2. Open input ROOT file
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
        if (c == 3) continue;
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    constexpr double kSigTlo = 4128.0;
    constexpr double kSigThi = 4272.0;
    constexpr double kCtlTlo = 6128.0;
    constexpr double kCtlThi = 6272.0;

    const std::vector<int> n_list = {3, 4, 5, 6, 8, 10, 12};
    const size_t num_n = n_list.size();

    // Thresholds per channel: T_c(n) = p0 + (n - 1.5) * G
    std::vector<double> thrs[32];
    for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        thrs[c].resize(num_n);
        for (size_t in = 0; in < num_n; ++in) {
            thrs[c][in] = p0_cal[c] + (n_list[in] - 1.5) * G_cal[c];
        }
    }

    // Per-event hit flags for planes:
    // 0: 18X sig, 1: 18X ctl, 2: 18Y sig, 3: 18Y ctl
    // Store as bitmasks for the 7 thresholds (fits in uint8_t)
    std::vector<uint8_t> ev_hits[4];
    for (int p = 0; p < 4; ++p) ev_hits[p].assign(nEnt, 0);

    // Also per-channel counts for n = 5 (index 2 in n_list)
    long long n5_ch_ns[32] = {};
    long long n5_ch_nc[32] = {};

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 30000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        for (int c = 0; c < 32; ++c) {
            if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // Signal window W2
            double vmin_s = 1e30; int imin_s = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_s) {
                    vmin_s = wf[i]; imin_s = i;
                }
            }
            double ph_s = -1e9;
            if (imin_s >= 0) {
                double ped_loc_s = 0;
                if (imin_s >= 12 && imin_s - 4 < n_samp) {
                    double s = 0;
                    for (int j = -12; j <= -4; ++j) s += wf[imin_s + j];
                    ped_loc_s = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped_loc_s = (n > 0) ? s / n : 0;
                }
                ph_s = ped_loc_s - vmin_s;
            }

            // Control window W2+2000 ns
            double vmin_c = 1e30; int imin_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < vmin_c) {
                    vmin_c = wf[i]; imin_c = i;
                }
            }
            double ph_c = -1e9;
            if (imin_c >= 0) {
                double ped_loc_c = 0;
                if (imin_c >= 12 && imin_c - 4 < n_samp) {
                    double s = 0;
                    for (int j = -12; j <= -4; ++j) s += wf[imin_c + j];
                    ped_loc_c = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped_loc_c = (n > 0) ? s / n : 0;
                }
                ph_c = ped_loc_c - vmin_c;
            }

            int p_sig = (c < 16) ? 0 : 2;
            int p_ctl = (c < 16) ? 1 : 3;

            for (size_t in = 0; in < num_n; ++in) {
                double T_cn = thrs[c][in];
                if (ph_s > T_cn) {
                    ev_hits[p_sig][iev] |= (1 << in);
                    if (in == 2) n5_ch_ns[c]++;
                }
                if (ph_c > T_cn) {
                    ev_hits[p_ctl][iev] |= (1 << in);
                    if (in == 2) n5_ch_nc[c]++;
                }
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n\n";

    const Long64_t N = (nEnt > 1) ? (nEnt - 1) : 0;
    const double N_d = (double)N;

    std::cout << "Formulas of reconstruct_corrected:\n";
    std::cout << "  N = N_events - 1 = " << N << "\n";
    std::cout << "  aX = N_X_C / N\n";
    std::cout << "  aY = N_Y_C / N\n";
    std::cout << "  rX = (N_X_S / N - aX) / (1 - aX)\n";
    std::cout << "  rY = (N_Y_S / N - aY) / (1 - aY)\n";
    std::cout << "  N_XY_real = (N_XY_same - N_XY_mixed) / ((1 - aX) * (1 - aY)) + rX * rY * N\n";
    std::cout << "  N_XY_real / (0.0360 * N)\n\n";

    std::cout << "========================================================================================================================\n";
    std::cout << "Step 3. hod2018 Threshold Scan in Photoelectrons (T_c(n) = p0_c + (n - 1.5) * G_c)\n";
    std::cout << "========================================================================================================================\n";
    std::cout << std::setw(6)  << "n"
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

    for (size_t in = 0; in < num_n; ++in) {
        int n_pe = n_list[in];
        uint8_t mask = (1 << in);

        long long N_X_S = 0, N_X_C = 0, N_Y_S = 0, N_Y_C = 0;
        long long N_XY_same = 0, N_XY_mixed = 0;

        for (Long64_t i = 0; i < N; ++i) {
            Long64_t next_i = i + 1;
            bool xs = (ev_hits[0][i] & mask);
            bool xc = (ev_hits[1][i] & mask);
            bool ys = (ev_hits[2][i] & mask);
            bool yc = (ev_hits[3][i] & mask);
            bool ys_next = (ev_hits[2][next_i] & mask);

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
                  << std::setw(6)  << n_pe
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

    // -------------------------------------------------------------------------
    // Per-channel (N_S - N_C) and N_C for n = 5
    // -------------------------------------------------------------------------
    std::cout << "================================================================================\n";
    std::cout << "hod2018 Per-Channel Counts at n = 5 p.e. in W2 and Control\n";
    std::cout << "================================================================================\n";
    std::cout << std::setw(8) << "Channel"
              << std::setw(8) << "Plane"
              << std::setw(12) << "T_c(5)[mV]"
              << std::setw(12) << "N_S"
              << std::setw(12) << "N_C"
              << std::setw(14) << "N_S - N_C"
              << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        std::string pl = (c < 16) ? "18X" : "18Y";
        long long ns = n5_ch_ns[c];
        long long nc = n5_ch_nc[c];
        std::cout << std::setw(8) << c
                  << std::setw(8) << pl
                  << std::fixed << std::setprecision(2)
                  << std::setw(12) << thrs[c][2] // index 2 is n = 5
                  << std::setw(12) << ns
                  << std::setw(12) << nc
                  << std::setw(14) << (ns - nc)
                  << "\n";
    }
    std::cout << "================================================================================\n\n";
}

