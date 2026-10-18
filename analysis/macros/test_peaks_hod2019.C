#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <iostream>
#include <vector>
#include <algorithm>

void test_peaks_hod2019(const char* fn = "lectura_vx1740_ext_nim_1M_v8.root") {
    TFile f(fn); TTree* t = (TTree*)f.Get("data_dig");
    std::vector<float>* tw = nullptr; t->SetBranchAddress("time", &tw);
    std::vector<float>* cw[64] = {};
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("time", 1);
    for (int c = 32; c < 64; ++c) {
        t->SetBranchStatus(Form("Ch%d", c), 1);
        t->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    t->SetCacheSize(256LL << 20);

    TH1D* h[64] = {};
    for (int c = 32; c < 64; ++c) h[c] = new TH1D(Form("h_%d", c), "", 240, 0, 60);

    for (Long64_t iev = 0; iev < 20000; ++iev) { // first 20k to inspect
        t->GetEntry(iev);
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        int n_samp = tv.size();
        for (int c = 32; c < 64; ++c) {
            if (!cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];
            double vm = 1e30; int im = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= 6128.0 && tv[i] <= 6272.0 && wf[i] < vm) { vm = wf[i]; im = i; }
            }
            if (im >= 0) {
                double ped = 0;
                if (im >= 12 && im - 4 < n_samp) {
                    double s = 0;
                    for (int j = -12; j <= -4; ++j) s += wf[im + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                h[c]->Fill(ped - vm);
            }
        }
    }

    std::cout << "ch  max1  max2  G_init\n";
    for (int c = 32; c < 64; ++c) {
        // Find local maxima above 1 mV
        // smooth a bit or simple peak search
        std::vector<double> pks;
        for (int b = h[c]->FindBin(1.0); b < h[c]->GetNbinsX() - 1; ++b) {
            double v = h[c]->GetBinContent(b);
            if (v > h[c]->GetBinContent(b-1) && v > h[c]->GetBinContent(b+1) && v > 5) {
                pks.push_back(h[c]->GetBinCenter(b));
            }
        }
        double m1 = pks.size() > 0 ? pks[0] : 0;
        double m2 = pks.size() > 1 ? pks[1] : 0;
        std::cout << c << "  " << m1 << "  " << m2 << "  " << (m2 - m1) << "\n";
    }
}
