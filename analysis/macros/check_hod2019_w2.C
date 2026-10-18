#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>

void check_hod2019_w2(const char* fn = "lectura_vx1740_ext_nim_1M_v8.root") {
    TFile f(fn); TTree* t = (TTree*)f.Get("data_dig");
    std::vector<float>* tw = nullptr; t->SetBranchAddress("time", &tw);
    std::vector<float>* cw[64] = {};
    for (int c = 32; c < 64; ++c) t->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    long long cnt[64] = {};
    for (Long64_t iev = 0; iev < t->GetEntries(); ++iev) {
        t->GetEntry(iev);
        for (int c = 32; c < 64; ++c) {
            if (!cw[c] || cw[c]->empty()) continue;
            double s = 0; int n = 0;
            for (size_t i = 0; i < tw->size(); ++i) if ((*tw)[i] < 2500.0) { s += (*cw[c])[i]; n++; }
            if (n < 2) continue;
            double ped = s / n, vm = 1e30; int im = -1;
            for (size_t i = 0; i < tw->size(); ++i)
                if ((*tw)[i] >= 4128.0 && (*tw)[i] <= 4272.0 && (*cw[c])[i] < vm) { vm = (*cw[c])[i]; im = i; }
            if (im >= 0 && (ped - vm) > 35.0) cnt[c]++;
        }
    }
    std::cout << "channel count\n";
    for (int c = 32; c < 64; ++c) std::cout << c << " " << cnt[c] << "\n";
}

