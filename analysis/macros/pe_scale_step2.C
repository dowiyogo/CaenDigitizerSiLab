#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

// 6-peak model
// par[0] = p0
// par[1] = G
// par[2] = s0
// par[3] = s1
// par[4..9] = A1..A6
Double_t fit_6peaks(Double_t* x, Double_t* par) {
    Double_t v = x[0];
    Double_t p0 = par[0];
    Double_t G  = par[1];
    Double_t s0 = par[2];
    Double_t s1 = par[3];
    Double_t sum = 0.0;
    for (int n = 1; n <= 6; ++n) {
        Double_t mean_n = p0 + (n - 1) * G;
        Double_t s2 = s0 * s0 + (n - 1) * s1 * s1;
        if (s2 <= 1e-6) continue;
        Double_t sigma_n = std::sqrt(s2);
        Double_t diff = (v - mean_n) / sigma_n;
        Double_t A_n = par[3 + n];
        sum += A_n * std::exp(-0.5 * diff * diff);
    }
    return sum;
}

void pe_scale_step2(const char* input_file = "lectura_vx1740_ext_nim_1M_v8.root") {
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

    std::vector<float>* cw[64] = {};
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("time", 1);
    for (int c = 0; c < 64; ++c) {
        if (c == 3) continue; // dead channel
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    // 0.25 mV bins from 0 to 60 mV (240 bins)
    constexpr int kNBins = 240;
    constexpr double kMaxPh = 60.0;
    TH1D* h_ph[64] = {};
    for (int c = 0; c < 64; ++c) {
        if (c == 3) continue;
        h_ph[c] = new TH1D(Form("h_ph_ch%d", c), Form("Ch %d Control ph;ph [mV];Events / 0.25 mV", c), kNBins, 0.0, kMaxPh);
    }

    constexpr double kCtlTlo = 6128.0;
    constexpr double kCtlThi = 6272.0;

    // Also store ph per event for channels to compute exact rates
    // To save memory, we can compute hits above thresholds directly from histograms or array
    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 30000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        for (int c = 0; c < 64; ++c) {
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
                h_ph[c]->Fill(ph_c);
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n\n";

    // Plane helper
    auto get_plane_str = [](int c) -> std::string {
        if (c >= 0 && c <= 15) return "18X";
        if (c >= 16 && c <= 31) return "18Y";
        if (c >= 32 && c <= 47) return "19Y";
        if (c >= 48 && c <= 63) return "19X";
        return "UNKNOWN";
    };

    struct FitOutput {
        int ch = -1;
        std::string plane;
        double p0_start = 0, G_start = 0;
        double p0 = 0, G = 0, s0 = 0, s1 = 0;
        double chi2 = 0;
        int ndf = 0;
        double chi2_ndf = 0;
        int status = -1;
    };
    std::vector<FitOutput> fit_results(64);

    std::cout << "========================================================================================================\n";
    std::cout << "Step 2b. Photoelectron Scale Fit (Sum of 6 Gaussians) per Channel\n";
    std::cout << "mean_n = p0 + (n-1)*G, sigma_n = sqrt(s0^2 + (n-1)*s1^2)\n";
    std::cout << "========================================================================================================\n";
    std::cout << std::setw(8) << "Channel"
              << std::setw(8) << "Plane"
              << std::setw(12) << "p0_start"
              << std::setw(12) << "G_start"
              << std::setw(10) << "p0 [mV]"
              << std::setw(10) << "G [mV]"
              << std::setw(10) << "s0 [mV]"
              << std::setw(10) << "s1 [mV]"
              << std::setw(12) << "chi2/ndf"
              << std::setw(12) << "fit status"
              << "\n";
    std::cout << "--------------------------------------------------------------------------------------------------------\n";

    for (int c = 0; c < 64; ++c) {
        if (c == 3 || !h_ph[c]) continue;
        std::string pl = get_plane_str(c);

        double p0_start = 4.3;
        double G_start  = 4.0;
        double s0_start = 0.8;
        double s1_start = 0.3;

        if (c >= 32) {
            // Find first two maxima above 1 mV
            TH1D* hs = (TH1D*)h_ph[c]->Clone("hs");
            hs->Smooth(1);
            std::vector<double> pks;
            int b_start = hs->FindBin(1.0);
            int b_end   = hs->FindBin(25.0);
            for (int b = b_start; b <= b_end; ++b) {
                double val = hs->GetBinContent(b);
                if (val > hs->GetBinContent(b-1) && val >= hs->GetBinContent(b+1) && val > 2.0) {
                    double pos = hs->GetBinCenter(b);
                    if (pks.empty() || pos > pks.back() + 0.6) {
                        pks.push_back(pos);
                    }
                }
            }
            delete hs;

            if (pks.size() >= 2) {
                p0_start = pks[0];
                G_start  = pks[1] - pks[0];
                if (G_start < 0.5) G_start = 1.5;
            } else if (pks.size() == 1) {
                p0_start = pks[0];
                G_start  = 2.0;
            } else {
                p0_start = 2.0;
                G_start  = 2.0;
            }
        }

        double fit_lo = p0_start - G_start / 2.0;
        double fit_hi = p0_start + 5.5 * G_start;
        if (fit_lo < 0.25) fit_lo = 0.25;
        if (fit_hi > 59.0) fit_hi = 59.0;

        TF1 f_fit(Form("f_fit_%d", c), fit_6peaks, fit_lo, fit_hi, 10);
        f_fit.SetParameter(0, p0_start);
        f_fit.SetParameter(1, G_start);
        f_fit.SetParameter(2, s0_start);
        f_fit.SetParameter(3, s1_start);

        // Parameter limits to prevent unphysical divergence
        f_fit.SetParLimits(0, 0.5, 12.0);
        f_fit.SetParLimits(1, 0.5, 12.0);
        f_fit.SetParLimits(2, 0.1, 3.5);
        f_fit.SetParLimits(3, 0.0, 3.5);

        for (int n = 1; n <= 6; ++n) {
            double mu = p0_start + (n - 1) * G_start;
            double A_init = h_ph[c]->GetBinContent(h_ph[c]->FindBin(mu));
            if (A_init <= 0) A_init = 1.0;
            f_fit.SetParameter(3 + n, A_init);
            f_fit.SetParLimits(3 + n, 0.0, h_ph[c]->GetMaximum() * 2.0);
        }

        TFitResultPtr res = h_ph[c]->Fit(&f_fit, "QRN0S", "", fit_lo, fit_hi);
        int status = (int)res;

        double p0 = f_fit.GetParameter(0);
        double G  = f_fit.GetParameter(1);
        double s0 = f_fit.GetParameter(2);
        double s1 = f_fit.GetParameter(3);
        double chi2 = f_fit.GetChisquare();
        int ndf = f_fit.GetNDF();
        double chi2_ndf = (ndf > 0) ? (chi2 / ndf) : 0.0;

        fit_results[c] = {c, pl, p0_start, G_start, p0, G, s0, s1, chi2, ndf, chi2_ndf, status};

        std::cout << std::fixed
                  << std::setw(8)  << c
                  << std::setw(8)  << pl
                  << std::setw(12) << std::setprecision(2) << p0_start
                  << std::setw(12) << std::setprecision(2) << G_start
                  << std::setw(10) << std::setprecision(2) << p0
                  << std::setw(10) << std::setprecision(2) << G
                  << std::setw(10) << std::setprecision(2) << s0
                  << std::setw(10) << std::setprecision(2) << s1
                  << std::setw(12) << std::setprecision(2) << chi2_ndf
                  << std::setw(12) << status
                  << "\n";
    }
    std::cout << "========================================================================================================\n\n";

    // -------------------------------------------------------------------------
    // 2b (continued): Histograms for channels 12, 14, 24, 48, 50, 36 (0 to 30 mV, 0.25 mV bins)
    // -------------------------------------------------------------------------
    const std::vector<int> target_chs = {12, 14, 24, 48, 50, 36};
    for (int c : target_chs) {
        std::cout << "================================================================================\n";
        std::cout << "0.25 mV Histogram [0, 30] mV for Channel " << c << " (" << get_plane_str(c) << ")\n";
        std::cout << "================================================================================\n";
        std::cout << std::setw(6) << "bin"
                  << std::setw(16) << "bin_center_mV"
                  << std::setw(14) << "content" << "\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        int max_b = h_ph[c]->FindBin(30.0 - 1e-4);
        for (int b = 1; b <= max_b; ++b) {
            std::cout << std::setw(6) << b
                      << std::fixed << std::setprecision(3)
                      << std::setw(16) << h_ph[c]->GetBinCenter(b)
                      << std::setw(14) << (long long)h_ph[c]->GetBinContent(b)
                      << "\n";
        }
        std::cout << "================================================================================\n\n";
    }

    // -------------------------------------------------------------------------
    // 2c. Out-of-time pulse rate in control window above n p.e. (n in {3, 5, 8, 12})
    // T_c(n) = p0_c + (n - 1.5) * G_c
    // rate = hits / (N_events * 144 ns) in kHz
    // -------------------------------------------------------------------------
    std::cout << "========================================================================================================\n";
    std::cout << "Step 2c. Out-of-Time Pulse Rates in Control Window [6128, 6272] ns (144 ns) above n p.e. [kHz]\n";
    std::cout << "T_c(n) = p0_c + (n - 1.5) * G_c,  rate = hits / (N_events * 144 ns) [kHz]\n";
    std::cout << "========================================================================================================\n";
    std::cout << std::setw(8) << "Channel"
              << std::setw(8) << "Plane"
              << std::setw(14) << "rate_3pe [kHz]"
              << std::setw(14) << "rate_5pe [kHz]"
              << std::setw(14) << "rate_8pe [kHz]"
              << std::setw(14) << "rate_12pe [kHz]"
              << std::setw(10) << "T(3)[mV]"
              << std::setw(10) << "T(5)[mV]"
              << std::setw(10) << "T(8)[mV]"
              << std::setw(10) << "T(12)[mV]"
              << "\n";
    std::cout << "--------------------------------------------------------------------------------------------------------\n";

    const double time_total_s = (double)nEnt * 144.0e-9; // total time in seconds
    const double factor_khz = 1.0 / (time_total_s * 1000.0);

    const std::vector<int> n_pe_list = {3, 5, 8, 12};

    for (int c = 0; c < 64; ++c) {
        if (c == 3 || !h_ph[c]) continue;
        const auto& fr = fit_results[c];
        double p0 = fr.p0;
        double G  = fr.G;
        // If fit was invalid, use start values
        if (fr.status != 0 || G <= 0.2 || G > 15.0) {
            p0 = fr.p0_start;
            G  = fr.G_start;
        }

        std::vector<double> rates(4, 0.0);
        std::vector<double> thrs(4, 0.0);

        for (size_t in = 0; in < n_pe_list.size(); ++in) {
            int n_pe = n_pe_list[in];
            double T_c = p0 + (n_pe - 1.5) * G;
            thrs[in] = T_c;
            int bin_thr = h_ph[c]->FindBin(T_c);
            double hits = h_ph[c]->Integral(bin_thr, kNBins + 1); // include overflow
            rates[in] = hits * factor_khz;
        }

        std::cout << std::fixed
                  << std::setw(8)  << c
                  << std::setw(8)  << fr.plane
                  << std::setw(14) << std::setprecision(3) << rates[0]
                  << std::setw(14) << std::setprecision(3) << rates[1]
                  << std::setw(14) << std::setprecision(3) << rates[2]
                  << std::setw(14) << std::setprecision(3) << rates[3]
                  << std::setw(10) << std::setprecision(1) << thrs[0]
                  << std::setw(10) << std::setprecision(1) << thrs[1]
                  << std::setw(10) << std::setprecision(1) << thrs[2]
                  << std::setw(10) << std::setprecision(1) << thrs[3]
                  << "\n";
    }
    std::cout << "========================================================================================================\n\n";

    // Write fit parameters to text file for use by Step 3 and Step 4
    std::ofstream pf("outputs/pe_calibration_v8.txt");
    pf << "# channel plane p0 G s0 s1 chi2 ndf status\n";
    for (int c = 0; c < 64; ++c) {
        if (c == 3 || !h_ph[c]) continue;
        const auto& fr = fit_results[c];
        double p0 = fr.p0;
        double G  = fr.G;
        if (fr.status != 0 || G <= 0.2 || G > 15.0) {
            p0 = fr.p0_start;
            G  = fr.G_start;
        }
        pf << c << " " << fr.plane << " " << p0 << " " << G << " " << fr.s0 << " " << fr.s1 << " "
           << fr.chi2 << " " << fr.ndf << " " << fr.status << "\n";
    }
    pf.close();
    std::cout << "Fit calibration saved to outputs/pe_calibration_v8.txt\n";
}

