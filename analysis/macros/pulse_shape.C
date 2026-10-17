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

    const int kMinOff = -12;
    const int kMaxOff = 40;
    const int nOff = kMaxOff - kMinOff + 1; // 53 points

    std::vector<double> sum18_glob(nOff, 0.0), sum18_loc(nOff, 0.0);
    std::vector<long long> count18(nOff, 0);
    long long n_pulses18 = 0;

    std::vector<double> sum19_glob(nOff, 0.0), sum19_loc(nOff, 0.0);
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
            double ped_glob = s / nb;

            double vmin = 1e30; int imin = -1;
            for (int i = 0; i < nSamp; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0) {
                    if (wf[i] < vmin) { vmin = wf[i]; imin = i; }
                }
            }
            if (imin < 0) continue;
            double ph = ped_glob - vmin;
            if (ph <= 50.0) continue;

            // Compute ped_loc (offsets -12..-4 from imin)
            double s_loc = 0; int nb_loc = 0;
            for (int k = -12; k <= -4; ++k) {
                int idx = imin + k;
                if (idx >= 0 && idx < nSamp) { s_loc += wf[idx]; nb_loc++; }
            }
            if (nb_loc < 9) continue;
            double ped_loc = s_loc / nb_loc;

            n_pulses18++;
            for (int k = kMinOff; k <= kMaxOff; ++k) {
                int idx = imin + k;
                if (idx >= 0 && idx < nSamp) {
                    int pt_idx = k - kMinOff;
                    sum18_glob[pt_idx] += (ped_glob - wf[idx]);
                    sum18_loc[pt_idx]  += (ped_loc - wf[idx]);
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
            double ped_glob = s / nb;

            double vmin = 1e30; int imin = -1;
            for (int i = 0; i < nSamp; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0) {
                    if (wf[i] < vmin) { vmin = wf[i]; imin = i; }
                }
            }
            if (imin < 0) continue;
            double ph = ped_glob - vmin;
            if (ph <= 35.0) continue;

            // Compute ped_loc (offsets -12..-4 from imin)
            double s_loc = 0; int nb_loc = 0;
            for (int k = -12; k <= -4; ++k) {
                int idx = imin + k;
                if (idx >= 0 && idx < nSamp) { s_loc += wf[idx]; nb_loc++; }
            }
            if (nb_loc < 9) continue;
            double ped_loc = s_loc / nb_loc;

            n_pulses19++;
            for (int k = kMinOff; k <= kMaxOff; ++k) {
                int idx = imin + k;
                if (idx >= 0 && idx < nSamp) {
                    int pt_idx = k - kMinOff;
                    sum19_glob[pt_idx] += (ped_glob - wf[idx]);
                    sum19_loc[pt_idx]  += (ped_loc - wf[idx]);
                    count19[pt_idx]++;
                }
            }
        }
    }

    f->Close();

    std::vector<double> mean18_glob(nOff, 0.0), mean18_loc(nOff, 0.0);
    std::vector<double> mean19_glob(nOff, 0.0), mean19_loc(nOff, 0.0);
    for (int i = 0; i < nOff; ++i) {
        if (count18[i] > 0) {
            mean18_glob[i] = sum18_glob[i] / count18[i];
            mean18_loc[i]  = sum18_loc[i] / count18[i];
        }
        if (count19[i] > 0) {
            mean19_glob[i] = sum19_glob[i] / count19[i];
            mean19_loc[i]  = sum19_loc[i] / count19[i];
        }
    }

    std::cout << "\nPulses analyzed:\n";
    std::cout << "  hod2018X (ph > 50 mV, t in [4000, 4400] ns): " << n_pulses18 << "\n";
    std::cout << "  hod2019X (ph > 35 mV, t in [4000, 4400] ns): " << n_pulses19 << "\n\n";

    std::cout << "=== Average Pulse Shapes (offsets -12..+40) ===\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(8) << "offset"
              << std::setw(12) << "t_ns"
              << std::setw(16) << "18X_glob"
              << std::setw(16) << "18X_loc"
              << std::setw(16) << "19X_glob"
              << std::setw(16) << "19X_loc"
              << "\n";
    std::cout << std::string(84, '-') << "\n";

    for (int k = kMinOff; k <= kMaxOff; ++k) {
        int idx = k - kMinOff;
        double t_ns = k * 16.0;
        std::cout << std::setw(8) << k
                  << std::setw(12) << std::setprecision(1) << t_ns
                  << std::setw(16) << std::setprecision(3) << mean18_glob[idx]
                  << std::setw(16) << std::setprecision(3) << mean18_loc[idx]
                  << std::setw(16) << std::setprecision(3) << mean19_glob[idx]
                  << std::setw(16) << std::setprecision(3) << mean19_loc[idx]
                  << "\n";
    }

    // Function to calculate n0, total, cumulative fractions, and n95
    auto eval_containment = [&](const std::string& label, const std::vector<double>& mean, int& out_n0, int& out_n95, double& out_tot) {
        int n0 = kMaxOff;
        for (int k = 1; k <= kMaxOff; ++k) {
            int idx = k - kMinOff;
            if (mean[idx] <= 0.0) {
                n0 = k;
                break;
            }
        }
        out_n0 = n0;

        double tot = 0.0;
        for (int k = -3; k <= n0; ++k) {
            int idx = k - kMinOff;
            tot += mean[idx] * 16.0;
        }
        out_tot = tot;

        int n95 = -1;
        for (int n = 8; n <= n0; ++n) {
            double cum = 0.0;
            for (int k = -3; k <= n; ++k) {
                int idx = k - kMinOff;
                cum += mean[idx] * 16.0;
            }
            double frac = (tot != 0.0) ? (cum / tot) : 0.0;
            if (n95 < 0 && frac >= 0.95) n95 = n;
        }
        out_n95 = n95;
    };

    int n0_18_glob = 40, n95_18_glob = -1; double tot_18_glob = 0;
    int n0_18_loc  = 40, n95_18_loc  = -1; double tot_18_loc  = 0;
    int n0_19_glob = 40, n95_19_glob = -1; double tot_19_glob = 0;
    int n0_19_loc  = 40, n95_19_loc  = -1; double tot_19_loc  = 0;

    eval_containment("hod2018X ped_glob", mean18_glob, n0_18_glob, n95_18_glob, tot_18_glob);
    eval_containment("hod2018X ped_loc",  mean18_loc,  n0_18_loc,  n95_18_loc,  tot_18_loc);
    eval_containment("hod2019X ped_glob", mean19_glob, n0_19_glob, n95_19_glob, tot_19_glob);
    eval_containment("hod2019X ped_loc",  mean19_loc,  n0_19_loc,  n95_19_loc,  tot_19_loc);

    std::cout << "\n=== Containment and n95 Summary ===\n";
    std::cout << "hod2018X ped_glob: n0 = " << n0_18_glob << ", total [-3, n0] = " << tot_18_glob << " mV*ns, n95 = " << n95_18_glob << "\n";
    std::cout << "hod2018X ped_loc:  n0 = " << n0_18_loc  << ", total [-3, n0] = " << tot_18_loc  << " mV*ns, n95 = " << n95_18_loc  << "\n";
    std::cout << "hod2019X ped_glob: n0 = " << n0_19_glob << ", total [-3, n0] = " << tot_19_glob << " mV*ns, n95 = " << n95_19_glob << "\n";
    std::cout << "hod2019X ped_loc:  n0 = " << n0_19_loc  << ", total [-3, n0] = " << tot_19_loc  << " mV*ns, n95 = " << n95_19_loc  << "\n";

    std::cout << "\n=== Cumulative fractions over [-3, n] for n = 8..n0 ===\n";

    auto print_series = [&](const std::string& title, const std::vector<double>& mean, int n0, double tot, int n95) {
        std::cout << "\n-- " << title << " (n0=" << n0 << ", total=" << tot << " mV*ns, n95=" << n95 << ") --\n";
        std::cout << std::setw(6) << "n"
                  << std::setw(12) << "t_n [ns]"
                  << std::setw(18) << "integral [-3,+n]"
                  << std::setw(14) << "fraction"
                  << "\n";
        std::cout << std::string(50, '-') << "\n";
        for (int n = 8; n <= n0; ++n) {
            double cum = 0.0;
            for (int k = -3; k <= n; ++k) {
                int idx = k - kMinOff;
                cum += mean[idx] * 16.0;
            }
            double frac = (tot != 0.0) ? (cum / tot) : 0.0;
            std::cout << std::setw(6) << n
                      << std::setw(12) << std::fixed << std::setprecision(1) << (n * 16.0)
                      << std::setw(18) << std::fixed << std::setprecision(3) << cum
                      << std::setw(14) << std::fixed << std::setprecision(4) << frac
                      << "\n";
        }
    };

    print_series("hod2018X ped_glob", mean18_glob, n0_18_glob, tot_18_glob, n95_18_glob);
    print_series("hod2018X ped_loc",  mean18_loc,  n0_18_loc,  tot_18_loc,  n95_18_loc);
    print_series("hod2019X ped_glob", mean19_glob, n0_19_glob, tot_19_glob, n95_19_glob);
    print_series("hod2019X ped_loc",  mean19_loc,  n0_19_loc,  tot_19_loc,  n95_19_loc);

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

    TCanvas c("c_ps", "Average Pulse Shape", 1000, 650);
    c.SetGrid();
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);

    TGraph g18_glob(nOff, t_arr.data(), mean18_glob.data());
    g18_glob.SetTitle("Average Pulse Shape with ped_glob vs ped_loc (v8);t - t_{min} [ns];Average Pulse [mV]");
    g18_glob.SetLineColor(kBlue + 1); g18_glob.SetLineWidth(2);
    g18_glob.SetMarkerStyle(20); g18_glob.SetMarkerSize(0.6); g18_glob.SetMarkerColor(kBlue + 1);

    TGraph g18_loc(nOff, t_arr.data(), mean18_loc.data());
    g18_loc.SetLineColor(kCyan + 2); g18_loc.SetLineWidth(2); g18_loc.SetLineStyle(2);
    g18_loc.SetMarkerStyle(24); g18_loc.SetMarkerSize(0.6); g18_loc.SetMarkerColor(kCyan + 2);

    TGraph g19_glob(nOff, t_arr.data(), mean19_glob.data());
    g19_glob.SetLineColor(kRed + 1); g19_glob.SetLineWidth(2);
    g19_glob.SetMarkerStyle(21); g19_glob.SetMarkerSize(0.6); g19_glob.SetMarkerColor(kRed + 1);

    TGraph g19_loc(nOff, t_arr.data(), mean19_loc.data());
    g19_loc.SetLineColor(kOrange + 2); g19_loc.SetLineWidth(2); g19_loc.SetLineStyle(2);
    g19_loc.SetMarkerStyle(25); g19_loc.SetMarkerSize(0.6); g19_loc.SetMarkerColor(kOrange + 2);

    g18_glob.Draw("ALP");
    g18_loc.Draw("LP SAME");
    g19_glob.Draw("LP SAME");
    g19_loc.Draw("LP SAME");

    TLegend leg(0.55, 0.70, 0.88, 0.88);
    leg.AddEntry(&g18_glob, Form("hod2018X ped_glob (n95=%d)", n95_18_glob), "lp");
    leg.AddEntry(&g18_loc,  Form("hod2018X ped_loc  (n95=%d)", n95_18_loc),  "lp");
    leg.AddEntry(&g19_glob, Form("hod2019X ped_glob (n95=%d)", n95_19_glob), "lp");
    leg.AddEntry(&g19_loc,  Form("hod2019X ped_loc  (n95=%d)", n95_19_loc),  "lp");
    leg.Draw();

    fs::path outpng = outdir / "pulse_shape_ped_comparison_v8.png";
    c.SaveAs(outpng.string().c_str());
    std::cout << "\nPlot saved to: " << outpng << "\n";
}
