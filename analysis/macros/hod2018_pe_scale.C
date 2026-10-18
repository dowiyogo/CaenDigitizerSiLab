#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

void hod2018_pe_scale(const char* input_file = "lectura_vx1740_ext_nim_1M_v8.root") {
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
        if (c == 3) continue; // dead channel
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    // 1 mV histograms [0, 150] mV in Control window for 31 live channels
    TH1D* h_ctl_1mv[32] = {};
    for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        h_ctl_1mv[c] = new TH1D(Form("h_ctl_1mv_ch%d", c), Form("Ch %d Control;ph [mV];Events", c), 150, 0.0, 150.0);
    }

    // 2 mV histograms [0, 200] mV in W2 and Control for channels 12, 14, 24
    const std::vector<int> target_chs = {12, 14, 24};
    TH1D* h_w2_2mv[32] = {};
    TH1D* h_ctl_2mv[32] = {};
    for (int c : target_chs) {
        h_w2_2mv[c] = new TH1D(Form("h_w2_2mv_ch%d", c), Form("Ch %d W2;ph [mV];Events", c), 100, 0.0, 200.0);
        h_ctl_2mv[c] = new TH1D(Form("h_ctl_2mv_ch%d", c), Form("Ch %d Control;ph [mV];Events", c), 100, 0.0, 200.0);
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
            if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // Control window [6128, 6272] ns
            double vmin_c = 1e30; int imin_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kCtlTlo && tv[i] <= kCtlThi && wf[i] < vmin_c) {
                    vmin_c = wf[i]; imin_c = i;
                }
            }
            if (imin_c >= 0) {
                double ped_loc_c = 0.0;
                if (imin_c >= 12 && imin_c - 4 < n_samp) {
                    double s_loc = 0.0;
                    for (int j = -12; j <= -4; ++j) s_loc += wf[imin_c + j];
                    ped_loc_c = s_loc / 9.0;
                } else {
                    double s_pre = 0; int n_pre = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s_pre += wf[i]; n_pre++; }
                    ped_loc_c = (n_pre > 0) ? (s_pre / n_pre) : 0.0;
                }
                double ph_c = ped_loc_c - vmin_c;
                h_ctl_1mv[c]->Fill(ph_c);

                if (h_ctl_2mv[c]) {
                    h_ctl_2mv[c]->Fill(ph_c);
                }
            }

            // Signal window [4128, 4272] ns (for channels 12, 14, 24)
            if (h_w2_2mv[c]) {
                double vmin_s = 1e30; int imin_s = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_s) {
                        vmin_s = wf[i]; imin_s = i;
                    }
                }
                if (imin_s >= 0) {
                    double ped_loc_s = 0.0;
                    if (imin_s >= 12 && imin_s - 4 < n_samp) {
                        double s_loc = 0.0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin_s + j];
                        ped_loc_s = s_loc / 9.0;
                    } else {
                        double s_pre = 0; int n_pre = 0;
                        for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s_pre += wf[i]; n_pre++; }
                        ped_loc_s = (n_pre > 0) ? (s_pre / n_pre) : 0.0;
                    }
                    double ph_s = ped_loc_s - vmin_s;
                    h_w2_2mv[c]->Fill(ph_s);
                }
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n\n";

    // -------------------------------------------------------------------------
    // 3a. TSpectrum::Search on each of the 31 live channels
    // -------------------------------------------------------------------------
    std::cout << "================================================================================\n";
    std::cout << "Photoelectron Peaks in Control Window [6128, 6272] ns (Local Pedestal, 1 mV bins)\n";
    std::cout << "TSpectrum::Search(h, sigma = 2 bins, threshold = 0.01)\n";
    std::cout << "================================================================================\n";
    std::cout << std::setw(8) << "Channel"
              << std::setw(10) << "N_peaks"
              << "   Peak positions [mV]\n";
    std::cout << "--------------------------------------------------------------------------------\n";

    TSpectrum spec(50);
    for (int c = 0; c < 32; ++c) {
        if (c == 3 || !h_ctl_1mv[c]) continue;
        int nfound = spec.Search(h_ctl_1mv[c], 2.0, "goff", 0.01);
        Double_t* xpeaks = spec.GetPositionX();
        std::vector<double> pks;
        if (xpeaks && nfound > 0) {
            pks.assign(xpeaks, xpeaks + nfound);
            std::sort(pks.begin(), pks.end());
        }
        std::cout << std::setw(8) << c
                  << std::setw(10) << nfound << "   ";
        for (size_t ip = 0; ip < pks.size(); ++ip) {
            std::cout << std::fixed << std::setprecision(1) << pks[ip] << (ip + 1 < pks.size() ? ", " : "");
        }
        std::cout << "\n";
    }
    std::cout << "================================================================================\n\n";

    // -------------------------------------------------------------------------
    // 3b. Full histogram (bin center, content) for channels 12, 14, 24 (1 mV bins)
    // -------------------------------------------------------------------------
    for (int c : target_chs) {
        std::cout << "================================================================================\n";
        std::cout << "Full Control Histogram for Channel " << c << " (Local Pedestal, 1 mV bins [0, 150] mV)\n";
        std::cout << "================================================================================\n";
        std::cout << std::setw(6) << "bin"
                  << std::setw(14) << "bin_center_mV"
                  << std::setw(14) << "content" << "\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        for (int b = 1; b <= 150; ++b) {
            std::cout << std::setw(6) << b
                      << std::fixed << std::setprecision(1)
                      << std::setw(14) << h_ctl_1mv[c]->GetBinCenter(b)
                      << std::setw(14) << (long long)h_ctl_1mv[c]->GetBinContent(b)
                      << "\n";
        }
        std::cout << "================================================================================\n\n";
    }

    // -------------------------------------------------------------------------
    // 3c. W2 minus Control histogram for channels 12, 14, 24 (2 mV bins [0, 200] mV)
    // -------------------------------------------------------------------------
    for (int c : target_chs) {
        std::cout << "================================================================================\n";
        std::cout << "W2 Minus Control Histogram for Channel " << c << " (2 mV bins [0, 200] mV)\n";
        std::cout << "================================================================================\n";
        std::cout << std::setw(6) << "bin"
                  << std::setw(14) << "bin_center_mV"
                  << std::setw(12) << "W2"
                  << std::setw(12) << "Control"
                  << std::setw(14) << "W2 - Control" << "\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        for (int b = 1; b <= 100; ++b) {
            long long w2_cnt = (long long)h_w2_2mv[c]->GetBinContent(b);
            long long ctl_cnt = (long long)h_ctl_2mv[c]->GetBinContent(b);
            std::cout << std::setw(6) << b
                      << std::fixed << std::setprecision(1)
                      << std::setw(14) << h_w2_2mv[c]->GetBinCenter(b)
                      << std::setw(12) << w2_cnt
                      << std::setw(12) << ctl_cnt
                      << std::setw(14) << (w2_cnt - ctl_cnt)
                      << "\n";
        }
        std::cout << "================================================================================\n\n";
    }
}

