#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <vector>

void compare_ch0_q() {
    TFile f("lectura_vx1740_ext_nim_1M_v8.root");
    TTree* t = (TTree*)f.Get("data_dig");
    std::vector<float>* tw = nullptr; t->SetBranchAddress("time", &tw);
    std::vector<float>* ch0 = nullptr; t->SetBranchAddress("Ch0", &ch0);
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("time", 1);
    t->SetBranchStatus("Ch0", 1);
    t->SetCacheSize(256LL << 20);

    TH1D h_w2("h_w2", "", 80, 0, 80000);
    TH1D h_ctl("h_ctl", "", 80, 0, 80000);

    for (Long64_t iev = 0; iev < t->GetEntries(); ++iev) {
        t->GetEntry(iev);
        if (!ch0 || ch0->empty()) continue;
        const auto& wf = *ch0;
        const auto& tv = *tw;
        int n_samp = wf.size();
        double s = 0; int np = 0;
        for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; np++; }
        if (np < 2) continue;
        double ped = s / np;

        // W2
        double vmin = 1e30; int imin = -1;
        for (int i = 0; i < n_samp; ++i) if (tv[i] >= 4128.0 && tv[i] <= 4272.0 && wf[i] < vmin) { vmin = wf[i]; imin = i; }
        if (imin >= 0 && (ped - vmin) > 50.0) {
            double ped_q = ped;
            if (imin >= 12 && imin - 4 < n_samp) {
                double sl = 0;
                for (int j = -12; j <= -4; ++j) sl += wf[imin + j];
                ped_q = sl / 9.0;
            }
            double Q = 0;
            double tilo = tv[imin] - 3 * 16.0, tihi = tv[imin] + 12 * 16.0;
            for (int i = 0; i < n_samp; ++i) if (tv[i] >= tilo && tv[i] <= tihi) Q += (ped_q - wf[i]) * 16.0;
            h_w2.Fill(Q);
        }

        // Ctl
        double vmin_c = 1e30; int imin_c = -1;
        for (int i = 0; i < n_samp; ++i) if (tv[i] >= 6128.0 && tv[i] <= 6272.0 && wf[i] < vmin_c) { vmin_c = wf[i]; imin_c = i; }
        if (imin_c >= 0 && (ped - vmin_c) > 50.0) {
            double ped_q_c = ped;
            if (imin_c >= 12 && imin_c - 4 < n_samp) {
                double sl = 0;
                for (int j = -12; j <= -4; ++j) sl += wf[imin_c + j];
                ped_q_c = sl / 9.0;
            }
            double Q_c = 0;
            double tilo_c = tv[imin_c] - 3 * 16.0, tihi_c = tv[imin_c] + 12 * 16.0;
            for (int i = 0; i < n_samp; ++i) if (tv[i] >= tilo_c && tv[i] <= tihi_c) Q_c += (ped_q_c - wf[i]) * 16.0;
            h_ctl.Fill(Q_c);
        }
    }
    std::cout << "bin Q_center W2 Ctl Diff\n";
    for (int b = 1; b <= 80; ++b) {
        if (h_w2.GetBinContent(b) > 0 || h_ctl.GetBinContent(b) > 0) {
            std::cout << b << " " << h_w2.GetBinCenter(b) << " " << (long long)h_w2.GetBinContent(b) << " " << (long long)h_ctl.GetBinContent(b) << " " << (long long)(h_w2.GetBinContent(b) - h_ctl.GetBinContent(b)) << "\n";
        }
    }
}

