#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>

namespace fs = std::filesystem;

void pulse_shape() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const char* fname = "lectura_vx1740_ext_nim_1M_v8.root";
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening " << fname << "\n";
        return;
    }
    TTree* t = (TTree*)f->Get("data_dig");
    if (!t) {
        std::cerr << "Error reading data_dig\n";
        f->Close();
        return;
    }

    std::vector<float>* tw = nullptr;
    t->SetBranchAddress("time", &tw);
    std::vector<float>* cw[64] = {nullptr};
    for (int c = 0; c < 64; ++c) {
        char bn[16];
        std::snprintf(bn, sizeof(bn), "Ch%d", c);
        if (t->GetBranch(bn)) t->SetBranchAddress(bn, &cw[c]);
    }

    const int kH19X[16] = {53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56};
    std::vector<int> ch18x, ch19x;
    for (int i = 0; i < 16; ++i) {
        ch18x.push_back(i);
        ch19x.push_back(kH19X[i]);
    }

    const int kMinOff = -10;
    const int kMaxOff = 40;
    const int nOff = kMaxOff - kMinOff + 1; // 51 points

    std::vector<double> sum18(nOff, 0.0);
    std::vector<long long> count18(nOff, 0);
    long long n_pulses18 = 0;

    std::vector<double> sum19(nOff, 0.0);
    std::vector<long long> count19(nOff, 0);
    long long n_pulses19 = 0;

    const Long64_t nEntries = std::min((Long64_t)20000, t->GetEntries());
    std::cout << "Processing " << nEntries << " entries for pulse shape analysis...\n";

    for (Long64_t iev = 0; iev < nEntries; ++iev) {
        t->GetEntry(iev);
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int nSamp = (int)tv.size();

        // Process hod2018X
        for (int c : ch18x) {
            if (!cw[c] || cw[c]->empty()) continue;
            const auto& wf = *cw[c];
            if ((int)wf.size() != nSamp) continue;

            double s = 0; int nb = 0;
            for (int i = 0; i < nSamp; ++i) {
                if (tv[i] < 2500.0) { s += wf[i]; nb++; }
            }
            if (nb < 2) continue;
            double ped = s / nb;

            double vmin = 1e30; int imin = -1;
            for (int i = 0; i < nSamp; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0) {
                    if (wf[i] < vmin) { vmin = wf[i]; imin = i; }
                }
            }
            if (imin < 0) continue;
            double ph = ped - vmin;
            if (ph <= 50.0) continue;

            n_pulses18++;
            for (int k = kMinOff; k <= kMaxOff; ++k) {
                int idx = imin + k;
                if (idx >= 0 && idx < nSamp) {
                    int pt_idx = k - kMinOff;
                    sum18[pt_idx] += (ped - wf[idx]);
                    count18[pt_idx]++;
                }
            }
        }

        // Process hod2019X
        for (int c : ch19x) {
            if (!cw[c] || cw[c]->empty()) continue;
            const auto& wf = *cw[c];
            if ((int)wf.size() != nSamp) continue;

            double s = 0; int nb = 0;
            for (int i = 0; i < nSamp; ++i) {
                if (tv[i] < 2500.0) { s += wf[i]; nb++; }
            }
            if (nb < 2) continue;
            double ped = s / nb;

            double vmin = 1e30; int imin = -1;
            for (int i = 0; i < nSamp; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0) {
                    if (wf[i] < vmin) { vmin = wf[i]; imin = i; }
                }
            }
            if (imin < 0) continue;
            double ph = ped - vmin;
            if (ph <= 35.0) continue;

            n_pulses19++;
            for (int k = kMinOff; k <= kMaxOff; ++k) {
                int idx = imin + k;
                if (idx >= 0 && idx < nSamp) {
                    int pt_idx = k - kMinOff;
                    sum19[pt_idx] += (ped - wf[idx]);
                    count19[pt_idx]++;
                }
            }
        }
    }

    f->Close();

    std::vector<double> mean18(nOff, 0.0);
    std::vector<double> mean19(nOff, 0.0);
    for (int i = 0; i < nOff; ++i) {
        if (count18[i] > 0) mean18[i] = sum18[i] / count18[i];
        if (count19[i] > 0) mean19[i] = sum19[i] / count19[i];
    }

    std::cout << "\nPulses analyzed:\n";
    std::cout << "  hod2018X (ph > 50 mV, t in [4000, 4400] ns): " << n_pulses18 << "\n";
    std::cout << "  hod2019X (ph > 35 mV, t in [4000, 4400] ns): " << n_pulses19 << "\n\n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(8) << "offset"
              << std::setw(12) << "t_ns"
              << std::setw(18) << "mean_hod2018X"
              << std::setw(18) << "mean_hod2019X"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    for (int k = kMinOff; k <= kMaxOff; ++k) {
        int idx = k - kMinOff;
        double t_ns = k * 16.0;
        std::cout << std::setw(8) << k
                  << std::setw(12) << std::setprecision(1) << t_ns
                  << std::setw(18) << std::setprecision(3) << mean18[idx]
                  << std::setw(18) << std::setprecision(3) << mean19[idx]
                  << "\n";
    }

    // Peak values
    double peak18 = *std::max_element(mean18.begin(), mean18.end());
    double peak19 = *std::max_element(mean19.begin(), mean19.end());

    // Offsets where trailing edge falls below 10% and 1% of peak
    int off10_18 = kMaxOff, off1_18 = kMaxOff;
    int off10_19 = kMaxOff, off1_19 = kMaxOff;

    int zero_idx = -kMinOff;
    for (int i = zero_idx; i < nOff; ++i) {
        int k = i + kMinOff;
        if (off10_18 == kMaxOff && mean18[i] < 0.10 * peak18) off10_18 = k;
        if (off1_18 == kMaxOff && mean18[i] < 0.01 * peak18) off1_18 = k;
        if (off10_19 == kMaxOff && mean19[i] < 0.10 * peak19) off10_19 = k;
        if (off1_19 == kMaxOff && mean19[i] < 0.01 * peak19) off1_19 = k;
    }

    std::cout << "\n=== Fall Time Thresholds ===\n";
    std::cout << "hod2018X: Peak = " << peak18 << " mV\n";
    std::cout << "  Falls below 10% of peak (" << 0.10 * peak18 << " mV) at offset = " << off10_18
              << " (" << off10_18 * 16.0 << " ns)\n";
    std::cout << "  Falls below  1% of peak (" << 0.01 * peak18 << " mV) at offset = " << off1_18
              << " (" << off1_18 * 16.0 << " ns)\n";

    std::cout << "hod2019X: Peak = " << peak19 << " mV\n";
    std::cout << "  Falls below 10% of peak (" << 0.10 * peak19 << " mV) at offset = " << off10_19
              << " (" << off10_19 * 16.0 << " ns)\n";
    std::cout << "  Falls below  1% of peak (" << 0.01 * peak19 << " mV) at offset = " << off1_19
              << " (" << off1_19 * 16.0 << " ns)\n";

    // Integral fractions: window [-3, +8] vs [-10, +40]
    double tot_int18 = 0.0, win_int18 = 0.0;
    double tot_int19 = 0.0, win_int19 = 0.0;

    for (int k = kMinOff; k <= kMaxOff; ++k) {
        int idx = k - kMinOff;
        tot_int18 += mean18[idx] * 16.0;
        tot_int19 += mean19[idx] * 16.0;
        if (k >= -3 && k <= 8) {
            win_int18 += mean18[idx] * 16.0;
            win_int19 += mean19[idx] * 16.0;
        }
    }

    double frac18 = (tot_int18 > 0) ? (win_int18 / tot_int18) : 0.0;
    double frac19 = (tot_int19 > 0) ? (win_int19 / tot_int19) : 0.0;

    std::cout << "\n=== Integral Fractions (current window [-3, +8] samples vs [-10, +40]) ===\n";
    std::cout << "hod2018X: Integral [-10, +40] = " << tot_int18 << " mV*ns, "
              << "[-3, +8] = " << win_int18 << " mV*ns, "
              << "Fraction = " << frac18 * 100.0 << " %\n";
    std::cout << "hod2019X: Integral [-10, +40] = " << tot_int19 << " mV*ns, "
              << "[-3, +8] = " << win_int19 << " mV*ns, "
              << "Fraction = " << frac19 * 100.0 << " %\n";

    // Create output directory and plot
    std::time_t now = std::time(nullptr);
    char ds[32];
    std::strftime(ds, sizeof(ds), "%Y%m%d_%H%M%S", std::localtime(&now));
    fs::path outdir = fs::path("outputs") / ("pulse_shape_" + std::string(ds));
    fs::create_directories(outdir);

    std::vector<double> t_arr(nOff);
    for (int k = kMinOff; k <= kMaxOff; ++k) {
        t_arr[k - kMinOff] = k * 16.0;
    }

    TCanvas c("c_ps", "Average Pulse Shape", 900, 600);
    c.SetGrid();
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);

    TGraph g18(nOff, t_arr.data(), mean18.data());
    g18.SetTitle("Average Cosmic Muon Pulse Shape (v8, t_{min} aligned);t - t_{min} [ns];Average (pedestal - waveform) [mV]");
    g18.SetLineColor(kBlue + 1);
    g18.SetLineWidth(2);
    g18.SetMarkerStyle(20);
    g18.SetMarkerSize(0.7);
    g18.SetMarkerColor(kBlue + 1);

    TGraph g19(nOff, t_arr.data(), mean19.data());
    g19.SetLineColor(kRed + 1);
    g19.SetLineWidth(2);
    g19.SetMarkerStyle(21);
    g19.SetMarkerSize(0.7);
    g19.SetMarkerColor(kRed + 1);

    g18.Draw("ALP");
    g19.Draw("LP SAME");

    TLegend leg(0.60, 0.72, 0.88, 0.88);
    leg.AddEntry(&g18, "hod2018X (ph > 50 mV)", "lp");
    leg.AddEntry(&g19, "hod2019X (ph > 35 mV)", "lp");
    leg.Draw();

    fs::path outpng = outdir / "pulse_shape_v8.png";
    c.SaveAs(outpng.string().c_str());
    std::cout << "\nPlot saved to: " << outpng << "\n";
}
