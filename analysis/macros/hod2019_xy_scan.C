// analysis/macros/hod2019_xy_scan.C
// Test whether hod2019 X and Y planes fire in the same event, and at what Q.
// Run v8, all entries, READ only.
//
// Invocation: root -l -b -q analysis/macros/hod2019_xy_scan.C

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// Mapeo autoritativo (de reconstruct_corrected.cxx líneas 58-59)
static const std::array<int,16> kH19X={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};

struct EvHod19 {
    double Qmax_X = -1e9;
    double phmax_X = -1e9;
    double t_X = 0;
    double Qmax_Y = -1e9;
    double phmax_Y = -1e9;
    double t_Y = 0;
};

void hod2019_xy_scan(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    std::time_t now = std::time(nullptr);
    char ds[32];
    std::strftime(ds, sizeof(ds), "%Y%m%d_%H%M%S", std::localtime(&now));
    const fs::path outdir = fs::path("outputs") / ("hod2019_xy_scan_" + std::string(ds));
    fs::create_directories(outdir);

    TFile* fin = TFile::Open(fname, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open " << fname << "\n";
        return;
    }
    TTree* t = static_cast<TTree*>(fin->Get("data_dig"));
    if (!t) {
        std::cerr << "ERROR: no data_dig in " << fname << "\n";
        fin->Close();
        return;
    }

    t->SetBranchStatus("*", 0);
    t->SetCacheSize(128LL << 20);
    std::vector<float>* tw = nullptr;
    t->SetBranchStatus("time", 1);
    t->SetBranchAddress("time", &tw);

    std::array<std::vector<float>*, 64> cw{};
    for (int c = 32; c <= 63; ++c) {
        char bn[16];
        std::snprintf(bn, sizeof(bn), "Ch%d", c);
        if (t->GetBranch(bn)) {
            t->SetBranchStatus(bn, 1);
            t->SetBranchAddress(bn, &cw[c]);
        }
    }

    const Long64_t nEnt = t->GetEntries();
    std::cout << "Loading " << nEnt << " entries from " << fname << "...\n";

    // 2a: Qmax histograms, [-1000, 40000] mV·ns, 100 mV·ns bins (410 bins)
    TH1D h_Qmax_X("h_Qmax_X", "hod2019 Q_{max};Q_{max} [mV#upointns];Entries", 410, -1000, 40000);
    TH1D h_Qmax_Y("h_Qmax_Y", "hod2019 Q_{max};Q_{max} [mV#upointns];Entries", 410, -1000, 40000);

    // 2b: phmax histograms, 5 mV bins in [0, 100] mV (20 bins)
    TH1D h_phmax_X("h_phmax_X", "hod2019 ph_{max};ph_{max} [mV];Entries", 20, 0, 100);
    TH1D h_phmax_Y("h_phmax_Y", "hod2019 ph_{max};ph_{max} [mV];Entries", 20, 0, 100);

    std::vector<EvHod19> events;
    events.reserve(nEnt);

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        t->GetEntry(iev);
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int n_samples = (int)tv.size();

        EvHod19 ev;

        // Process plane X
        for (int c : kH19X) {
            if (!cw[c] || cw[c]->empty()) continue;
            const auto& wf = *cw[c];
            double s = 0; int nb = 0;
            for (int i = 0; i < n_samples; ++i) {
                if (tv[i] < 2500.0) { s += wf[i]; nb++; }
            }
            double ped = (nb > 0) ? s / nb : 0.0;
            double vmin = 1e30; int imin = -1;
            for (int i = 0; i < n_samples; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0 && wf[i] < vmin) {
                    vmin = wf[i];
                    imin = i;
                }
            }
            double ph = (imin >= 0) ? (ped - vmin) : 0.0;
            double Q = 0.0;
            if (imin >= 0) {
                for (int k = imin - 3; k <= imin + 22; ++k) {
                    if (k >= 0 && k < n_samples) {
                        Q += (ped - wf[k]) * 16.0;
                    }
                }
            }
            double pt = (imin >= 0) ? tv[imin] : 0.0;

            if (ph > ev.phmax_X) ev.phmax_X = ph;
            if (Q > ev.Qmax_X) {
                ev.Qmax_X = Q;
                ev.t_X = pt;
            }
        }

        // Process plane Y
        for (int c : kH19Y) {
            if (!cw[c] || cw[c]->empty()) continue;
            const auto& wf = *cw[c];
            double s = 0; int nb = 0;
            for (int i = 0; i < n_samples; ++i) {
                if (tv[i] < 2500.0) { s += wf[i]; nb++; }
            }
            double ped = (nb > 0) ? s / nb : 0.0;
            double vmin = 1e30; int imin = -1;
            for (int i = 0; i < n_samples; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0 && wf[i] < vmin) {
                    vmin = wf[i];
                    imin = i;
                }
            }
            double ph = (imin >= 0) ? (ped - vmin) : 0.0;
            double Q = 0.0;
            if (imin >= 0) {
                for (int k = imin - 3; k <= imin + 22; ++k) {
                    if (k >= 0 && k < n_samples) {
                        Q += (ped - wf[k]) * 16.0;
                    }
                }
            }
            double pt = (imin >= 0) ? tv[imin] : 0.0;

            if (ph > ev.phmax_Y) ev.phmax_Y = ph;
            if (Q > ev.Qmax_Y) {
                ev.Qmax_Y = Q;
                ev.t_Y = pt;
            }
        }

        h_Qmax_X.Fill(ev.Qmax_X);
        h_Qmax_Y.Fill(ev.Qmax_Y);
        h_phmax_X.Fill(ev.phmax_X);
        h_phmax_Y.Fill(ev.phmax_Y);

        events.push_back(ev);
    }
    fin->Close();
    const size_t Ntot = events.size();
    std::cout << "Processed " << Ntot << " events.\n";

    // 2d: Delta t histogram: 16 ns bins in [-400, 400] ns (50 bins)
    TH1D h_dt_same("h_dt_same", "#Delta t = t_{Y} - t_{X} (Tag ph_{X}>35, Probe Q_{Y}>687);#Delta t [ns];Entries", 50, -400, 400);
    TH1D h_dt_mixed("h_dt_mixed", "#Delta t mixed;#Delta t [ns];Entries", 50, -400, 400);

    for (size_t i = 0; i < Ntot; ++i) {
        if (events[i].phmax_X > 35.0) {
            if (events[i].Qmax_Y > 687.0) {
                h_dt_same.Fill(events[i].t_Y - events[i].t_X);
            }
            size_t inext = (i + 1 < Ntot) ? (i + 1) : 0;
            if (events[inext].Qmax_Y > 687.0) {
                h_dt_mixed.Fill(events[inext].t_Y - events[i].t_X);
            }
        }
    }

    // Plots
    {
        TCanvas c1("c1_qmax", "Qmax", 900, 650);
        c1.SetLogy();
        h_Qmax_X.SetLineColor(kBlue + 1);
        h_Qmax_X.SetLineWidth(2);
        h_Qmax_Y.SetLineColor(kRed + 1);
        h_Qmax_Y.SetLineWidth(2);
        h_Qmax_X.Draw("HIST");
        h_Qmax_Y.Draw("HIST SAME");
        TLegend leg(0.65, 0.75, 0.88, 0.88);
        leg.AddEntry(&h_Qmax_X, "hod2019 X", "l");
        leg.AddEntry(&h_Qmax_Y, "hod2019 Y", "l");
        leg.Draw();
        c1.SaveAs((outdir / "qmax_distributions_hod2019.png").string().c_str());
    }
    {
        TCanvas c2("c2_phmax", "phmax", 900, 650);
        c2.SetLogy();
        h_phmax_X.SetLineColor(kBlue + 1);
        h_phmax_X.SetLineWidth(2);
        h_phmax_Y.SetLineColor(kRed + 1);
        h_phmax_Y.SetLineWidth(2);
        h_phmax_X.Draw("HIST");
        h_phmax_Y.Draw("HIST SAME");
        TLegend leg(0.65, 0.75, 0.88, 0.88);
        leg.AddEntry(&h_phmax_X, "hod2019 X", "l");
        leg.AddEntry(&h_phmax_Y, "hod2019 Y", "l");
        leg.Draw();
        c2.SaveAs((outdir / "phmax_distributions_hod2019.png").string().c_str());
    }
    {
        TCanvas c3("c3_dt", "deltat", 900, 650);
        h_dt_same.SetLineColor(kBlue + 1);
        h_dt_same.SetLineWidth(2);
        h_dt_mixed.SetLineColor(kGray + 2);
        h_dt_mixed.SetLineStyle(2);
        h_dt_mixed.SetLineWidth(2);
        h_dt_same.Draw("HIST");
        h_dt_mixed.Draw("HIST SAME");
        TLegend leg(0.65, 0.75, 0.88, 0.88);
        leg.AddEntry(&h_dt_same, "Same event", "l");
        leg.AddEntry(&h_dt_mixed, "Mixed (event i+1)", "l");
        leg.Draw();
        c3.SaveAs((outdir / "deltat_distribution_hod2019.png").string().c_str());
    }

    // Prepare text summary
    auto print_all = [&](std::ostream& os) {
        os << "=== hod2019 X-Y Coincidence Scan ===\n\n";
        os << "Run: " << fname << "\n";
        os << "Total entries: " << Ntot << "\n";
        os << "Output dir: " << outdir << "\n\n";

        // 2a: Qmax rebinned to 500 mV·ns in [0, 12000] mV·ns
        os << "--- 2a. Qmax distribution (rebinned to 500 mV·ns, [0, 12000] mV·ns) ---\n";
        os << std::setw(12) << "center" << std::setw(15) << "content_X" << std::setw(15) << "content_Y\n";
        os << "-------------------------------------------\n";
        TH1D* h_Qmax_X_500 = static_cast<TH1D*>(h_Qmax_X.Rebin(5, "h_Qmax_X_500"));
        TH1D* h_Qmax_Y_500 = static_cast<TH1D*>(h_Qmax_Y.Rebin(5, "h_Qmax_Y_500"));
        for (int b = 1; b <= h_Qmax_X_500->GetNbinsX(); ++b) {
            double c = h_Qmax_X_500->GetBinCenter(b);
            if (c >= 0.0 && c <= 12000.0) {
                os << std::setw(12) << std::fixed << std::setprecision(1) << c
                   << std::setw(15) << (long long)h_Qmax_X_500->GetBinContent(b)
                   << std::setw(15) << (long long)h_Qmax_Y_500->GetBinContent(b) << "\n";
            }
        }
        os << "\n";

        // 2b: phmax distribution, 5 mV bins in [0, 100] mV
        os << "--- 2b. phmax distribution (5 mV bins in [0, 100] mV) ---\n";
        os << std::setw(12) << "center" << std::setw(15) << "content_X" << std::setw(15) << "content_Y\n";
        os << "-------------------------------------------\n";
        for (int b = 1; b <= h_phmax_X.GetNbinsX(); ++b) {
            double c = h_phmax_X.GetBinCenter(b);
            os << std::setw(12) << std::fixed << std::setprecision(1) << c
               << std::setw(15) << (long long)h_phmax_X.GetBinContent(b)
               << std::setw(15) << (long long)h_phmax_Y.GetBinContent(b) << "\n";
        }
        os << "\n";

        // 2c: Internal coincidence
        const std::vector<double> c_cuts = {687, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000};
        auto do_coin = [&](bool tag_is_x) {
            os << "--- 2c. Internal Coincidence: Tag " << (tag_is_x ? "phmax_X > 35 mV, Probe Qmax_Y > c" : "phmax_Y > 35 mV, Probe Qmax_X > c") << " ---\n";
            os << std::setw(8) << "c" << std::setw(10) << "N_tag" << std::setw(10) << "N_same" << std::setw(10) << "N_mixed"
               << std::setw(14) << "eps" << std::setw(14) << "sigma\n";
            os << "------------------------------------------------------------------\n";
            for (double c : c_cuts) {
                long long N_tag = 0, N_same = 0, N_mixed = 0;
                for (size_t i = 0; i < Ntot; ++i) {
                    double tag_val = tag_is_x ? events[i].phmax_X : events[i].phmax_Y;
                    if (tag_val > 35.0) {
                        ++N_tag;
                        double probe_same = tag_is_x ? events[i].Qmax_Y : events[i].Qmax_X;
                        if (probe_same > c) ++N_same;
                        size_t inext = (i + 1 < Ntot) ? (i + 1) : 0;
                        double probe_mixed = tag_is_x ? events[inext].Qmax_Y : events[inext].Qmax_X;
                        if (probe_mixed > c) ++N_mixed;
                    }
                }
                double eps = (N_tag > 0) ? (double)(N_same - N_mixed) / N_tag : 0.0;
                double sig = (N_tag > 0) ? std::sqrt((double)(N_same + N_mixed)) / N_tag : 0.0;
                os << std::setw(8) << std::fixed << std::setprecision(0) << c
                   << std::setw(10) << N_tag
                   << std::setw(10) << N_same
                   << std::setw(10) << N_mixed
                   << std::setw(14) << std::setprecision(6) << eps
                   << std::setw(14) << std::setprecision(6) << sig << "\n";
            }
            os << "\n";
        };
        do_coin(true);   // Tag X, Probe Y
        do_coin(false);  // Tag Y, Probe X

        // 2d: Timing difference
        os << "--- 2d. Delta t = t_Y - t_X (Tag phmax_X > 35 mV, Probe Qmax_Y > 687 mV·ns) ---\n";
        os << std::setw(12) << "center" << std::setw(12) << "same" << std::setw(12) << "mixed\n";
        os << "------------------------------------\n";
        for (int b = 1; b <= h_dt_same.GetNbinsX(); ++b) {
            double c = h_dt_same.GetBinCenter(b);
            os << std::setw(12) << std::fixed << std::setprecision(1) << c
               << std::setw(12) << (long long)h_dt_same.GetBinContent(b)
               << std::setw(12) << (long long)h_dt_mixed.GetBinContent(b) << "\n";
        }
        os << "\n";
    };

    print_all(std::cout);
    {
        std::ofstream sf(outdir / "summary.txt");
        print_all(sf);
    }
    std::cout << "Done. Results written to " << outdir << "\n";
}

