#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

void muon_timing_scan() {
    std::vector<std::pair<std::string, std::string>> runs = {
        {"v1", "lectura_vx1740_ext_nim_1M.root"},
        {"v2", "lectura_vx1740_ext_nim_1M_v2.root"},
        {"v3", "lectura_vx1740_ext_nim_1M_v3.root"},
        {"v4", "lectura_vx1740_ext_nim_1M_v4.root"},
        {"v5", "lectura_vx1740_ext_nim_1M_v5.root"},
        {"v6", "lectura_vx1740_ext_nim_1M_v6.root"},
        {"v7", "lectura_vx1740_ext_nim_1M_v7.root"},
        {"v8", "lectura_vx1740_ext_nim_1M_v8.root"},
        {"v9", "lectura_vx1740_ext_nim_1M_v9.root"},
        {"v10", "lectura_vx1740_ext_nim_1M_v10.root"}
    };

    std::cout << std::left 
              << std::setw(6) << "Run" 
              << std::setw(17) << "Station"
              << std::setw(12) << "Mode [ns]"
              << std::setw(18) << "N [4000,4400]"
              << std::setw(18) << "N [3000,3400]"
              << std::setw(12) << "Total Hits"
              << "\n";
    std::cout << std::string(83, '-') << "\n";

    for (const auto& [lbl, fname] : runs) {
        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cout << std::setw(6) << lbl << " Failed to open\n";
            continue;
        }
        TTree* t = (TTree*)f->Get("data_dig");
        if (!t) {
            std::cout << std::setw(6) << lbl << " No data_dig\n";
            f->Close();
            continue;
        }

        std::vector<float>* tw = nullptr;
        t->SetBranchAddress("time", &tw);
        std::vector<float>* cw[64] = {nullptr};
        for (int c = 0; c < 64; ++c) {
            char bn[16];
            std::snprintf(bn, sizeof(bn), "Ch%d", c);
            if (t->GetBranch(bn)) t->SetBranchAddress(bn, &cw[c]);
        }

        // Bins of 16 ns from 2500 to 8000 ns: (8000-2500)/16 = 343.75 -> 344 bins
        TH1D h18("h18", "", 344, 2500, 2500 + 344 * 16.0);
        TH1D h19("h19", "", 344, 2500, 2500 + 344 * 16.0);

        Long64_t nEnt = std::min((Long64_t)3000, t->GetEntries());
        for (Long64_t iev = 0; iev < nEnt; ++iev) {
            t->GetEntry(iev);
            if (!tw || tw->empty()) continue;
            const auto& tv = *tw;
            const int nSamp = (int)tv.size();

            for (int c = 0; c < 64; ++c) {
                if (!cw[c] || cw[c]->empty()) continue;
                const auto& wf = *cw[c];
                if ((int)wf.size() != nSamp) continue;

                // Baseline
                double s = 0; int nb = 0;
                for (int i = 0; i < nSamp; ++i) {
                    if (tv[i] < 2500.0) { s += wf[i]; nb++; }
                }
                if (nb < 2) continue;
                double ped = s / nb;

                // Minimum in [2500, 8000] ns
                double vmin = 1e30; int imin = -1;
                for (int i = 0; i < nSamp; ++i) {
                    if (tv[i] >= 2500.0 && tv[i] <= 8000.0) {
                        if (wf[i] < vmin) { vmin = wf[i]; imin = i; }
                    }
                }
                if (imin < 0) continue;
                double ph = ped - vmin;
                double pt = tv[imin];

                if (c < 32 && ph > 50.0) {
                    h18.Fill(pt);
                } else if (c >= 32 && ph > 35.0) {
                    h19.Fill(pt);
                }
            }
        }

        auto analyze_h = [](TH1D& h, const std::string& run, const std::string& st) {
            double mode = h.GetEntries() > 0 ? h.GetBinCenter(h.GetMaximumBin()) : 0.0;
            int b4000 = h.FindBin(4000.0);
            int b4400 = h.FindBin(4399.9);
            double n_4000_4400 = h.Integral(b4000, b4400);

            int b3000 = h.FindBin(3000.0);
            int b3400 = h.FindBin(3399.9);
            double n_3000_3400 = h.Integral(b3000, b3400);

            std::cout << std::left 
                      << std::setw(6) << run
                      << std::setw(17) << st
                      << std::setw(12) << std::fixed << std::setprecision(1) << mode
                      << std::setw(18) << (long long)n_4000_4400
                      << std::setw(18) << (long long)n_3000_3400
                      << std::setw(12) << (long long)h.GetEntries()
                      << "\n";
        };

        analyze_h(h18, lbl, "hod2018 (ph>50)");
        analyze_h(h19, lbl, "hod2019 (ph>35)");

        f->Close();
    }
}
