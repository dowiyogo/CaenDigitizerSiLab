// analysis/macros/hod2018_diag_fix17.C
// TASK: ANALYSIS-FIX-17
// Multiplicity distributions (same vs control window) at n = 5 and 10 pe,
// and ch29 amplitude spectrum diagnostics at n = 10 pe on run v8 good time.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::string channels;
    std::string reason;
};

struct EventV8 {
    float ph_s[32]; // channels 0..31
    float ph_c[32];
};

struct ChannelCalib {
    int ch;
    std::string plane;
    double p0;
    double G;
    int status;
};

void hod2018_diag_fix17(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // 1. Read calibration parameters from outputs/pe_calibration_v8.txt
    ChannelCalib calib[32] = {};
    for (int c = 0; c < 32; ++c) {
        calib[c].ch = c;
        calib[c].plane = (c < 16) ? "18X" : "18Y";
        calib[c].p0 = 0.0;
        calib[c].G = 0.0;
        calib[c].status = -1;
    }

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
                calib[ch].ch = ch;
                calib[ch].plane = pl;
                calib[ch].p0 = p0;
                calib[ch].G  = G;
                calib[ch].status = status;
            }
        }
    }
    pf.close();

    auto get_thr = [&](int c, double n_pe) -> double {
        if (c == 3) return 1e9; // dead channel
        return calib[c].p0 + (n_pe - 1.5) * calib[c].G;
    };

    // 2. Read intervals from runs/good_intervals_2026.csv
    std::vector<Interval> intervals;
    std::ifstream fcsv("runs/good_intervals_2026.csv");
    if (!fcsv.is_open()) {
        std::cerr << "Error: runs/good_intervals_2026.csv not found!\n";
        return;
    }
    std::getline(fcsv, line); // header
    while (std::getline(fcsv, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Interval iv;
        if (std::getline(iss, iv.run, ',') &&
            std::getline(iss, iv.start_local, ',') &&
            std::getline(iss, iv.end_local, ',') &&
            std::getline(iss, iv.status, ',') &&
            std::getline(iss, iv.channels, ',') &&
            std::getline(iss, iv.reason)) {
            intervals.push_back(iv);
        }
    }
    fcsv.close();

    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // -04:00 offset

    // Nominal signal window W2 and control window
    constexpr double kNomSigTlo = 4128.0, kNomSigThi = 4272.0;
    constexpr double kNomCtlTlo = 6128.0, kNomCtlThi = 6272.0;

    TFile fin(fname, "READ");
    if (fin.IsZombie()) {
        std::cerr << "Error opening " << fname << "\n";
        return;
    }
    TTree* tree = (TTree*)fin.Get("data_dig");
    if (!tree) {
        std::cerr << "Error finding tree data_dig in " << fname << "\n";
        return;
    }

    Long64_t nEnt = tree->GetEntries();
    std::vector<float>* tw = nullptr;
    Long_t timestamp = 0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("time", 1);
    tree->SetBranchAddress("time", &tw);
    tree->SetBranchStatus("timestamp", 1);
    tree->SetBranchAddress("timestamp", &timestamp);

    std::vector<float>* cw[32] = {};
    for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        tree->SetBranchStatus(Form("Ch%d", c), 1);
        tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    tree->SetCacheSize(256LL << 20);

    std::vector<EventV8> kept_events;
    kept_events.reserve(120000);

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (!tw || tw->empty()) continue;

        double t_sec = boot_ref_epoch + ((double)timestamp) * 1e-6 + tz_offset;
        time_t t_local = (time_t)t_sec;
        struct tm tm_val; gmtime_r(&t_local, &tm_val);
        char buf_min[32];
        std::snprintf(buf_min, sizeof(buf_min), "%04d-%02d-%02d %02d:%02d",
                      tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday, tm_val.tm_hour, tm_val.tm_min);
        std::string s_min(buf_min);

        bool is_exclude = false;
        bool is_hv_off = false;
        for (const auto& iv : intervals) {
            if (iv.run == "v8") {
                if (s_min >= iv.start_local && s_min <= iv.end_local) {
                    if (iv.status == "EXCLUDE") is_exclude = true;
                    if (iv.status == "HV_OFF")  is_hv_off = true;
                }
            }
        }

        // Good time only
        if (is_exclude || is_hv_off) continue;

        const auto& tv = *tw;
        const int n_samp = (int)tv.size();

        EventV8 ev{};

        for (int c = 0; c < 32; ++c) {
            if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
            const auto& wf = *cw[c];

            // Signal window W2 [4128, 4272] ns
            double vmin_s = 1e30; int imin_s = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kNomSigTlo && tv[i] <= kNomSigThi && wf[i] < vmin_s) {
                    vmin_s = wf[i]; imin_s = i;
                }
            }
            double ph_s = -1e9;
            if (imin_s >= 0) {
                double ped = 0;
                if (imin_s >= 12 && imin_s - 4 < n_samp) {
                    double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_s + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                ph_s = ped - vmin_s;
            }

            // Control window [6128, 6272] ns
            double vmin_c = 1e30; int imin_c = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kNomCtlTlo && tv[i] <= kNomCtlThi && wf[i] < vmin_c) {
                    vmin_c = wf[i]; imin_c = i;
                }
            }
            double ph_c = -1e9;
            if (imin_c >= 0) {
                double ped = 0;
                if (imin_c >= 12 && imin_c - 4 < n_samp) {
                    double s = 0; for (int j = -12; j <= -4; ++j) s += wf[imin_c + j];
                    ped = s / 9.0;
                } else {
                    double s = 0; int n = 0;
                    for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; n++; }
                    ped = (n > 0) ? s / n : 0;
                }
                ph_c = ped - vmin_c;
            }

            ev.ph_s[c] = (float)ph_s;
            ev.ph_c[c] = (float)ph_c;
        }

        kept_events.push_back(ev);
    }

    Long64_t n_v8 = (Long64_t)kept_events.size();
    std::cout << "Loaded v8: " << n_v8 << " good-time events.\n\n";

    // =========================================================================
    // PART B: Multiplicity distributions
    // =========================================================================
    const std::vector<double> n_pe_list = {5.0, 10.0};

    for (double n_pe : n_pe_list) {
        std::cout << "=======================================================================================================================================================\n";
        std::cout << "PART B: Multiplicity Distributions at n = " << (int)n_pe << " p.e. (Run v8 Good Time Only, N = " << n_v8 << ")\n";
        std::cout << "=======================================================================================================================================================\n";

        // Counts of n = 0..8, >=9 (10 bins)
        long long all_sig_x[10] = {}, all_ctl_x[10] = {};
        long long all_sig_y[10] = {}, all_ctl_y[10] = {};

        long long rest_sig_x[10] = {}, rest_ctl_x[10] = {};
        long long rest_sig_y[10] = {}, rest_ctl_y[10] = {};

        long long n_high_mult_ev = 0;
        long long ch_high_mult_cnt[32] = {};

        for (Long64_t i = 0; i < n_v8; ++i) {
            int nx_s = 0, nx_c = 0;
            int ny_s = 0, ny_c = 0;
            bool hit_s[32] = {};

            for (int c = 0; c < 32; ++c) {
                if (c == 3) continue;
                double thr = get_thr(c, n_pe);
                if (kept_events[i].ph_s[c] > thr) {
                    hit_s[c] = true;
                    if (c < 16) nx_s++;
                    else        ny_s++;
                }
                if (kept_events[i].ph_c[c] > thr) {
                    if (c < 16) nx_c++;
                    else        ny_c++;
                }
            }

            int bin_nx_s = (nx_s >= 9) ? 9 : nx_s;
            int bin_nx_c = (nx_c >= 9) ? 9 : nx_c;
            int bin_ny_s = (ny_s >= 9) ? 9 : ny_s;
            int bin_ny_c = (ny_c >= 9) ? 9 : ny_c;

            // All events
            all_sig_x[bin_nx_s]++;
            all_ctl_x[bin_nx_c]++;
            all_sig_y[bin_ny_s]++;
            all_ctl_y[bin_ny_c]++;

            // Restricted events
            if (nx_s >= 1 && ny_s >= 1) {
                rest_sig_x[bin_nx_s]++;
                rest_sig_y[bin_ny_s]++;
            }
            if (nx_c >= 1 && ny_c >= 1) {
                rest_ctl_x[bin_nx_c]++;
                rest_ctl_y[bin_ny_c]++;
            }

            // High multiplicity check: nX_sig >= 4 OR nY_sig >= 4
            if (nx_s >= 4 || ny_s >= 4) {
                n_high_mult_ev++;
                for (int c = 0; c < 32; ++c) {
                    if (hit_s[c]) ch_high_mult_cnt[c]++;
                }
            }
        }

        // Print All Events Table
        std::cout << "All events (N = " << n_v8 << "):\n";
        std::cout << std::left << std::setw(6) << "n"
                  << std::right << std::setw(14) << "count_sig_X"
                  << std::setw(14) << "count_ctl_X"
                  << std::setw(14) << "count_sig_Y"
                  << std::setw(14) << "count_ctl_Y"
                  << "\n";
        std::cout << std::string(62, '-') << "\n";

        for (int n = 0; n < 10; ++n) {
            std::string n_str = (n == 9) ? ">=9" : std::to_string(n);
            std::cout << std::left << std::setw(6) << n_str
                      << std::right << std::setw(14) << all_sig_x[n]
                      << std::setw(14) << all_ctl_x[n]
                      << std::setw(14) << all_sig_y[n]
                      << std::setw(14) << all_ctl_y[n]
                      << "\n";
        }
        std::cout << std::string(62, '-') << "\n\n";

        // Print Restricted Events Table
        std::cout << "Restricted events (nX_sig >= 1 AND nY_sig >= 1; ctl: nX_ctl >= 1 AND nY_ctl >= 1):\n";
        std::cout << std::left << std::setw(6) << "n"
                  << std::right << std::setw(14) << "count_sig_X"
                  << std::setw(14) << "count_ctl_X"
                  << std::setw(14) << "count_sig_Y"
                  << std::setw(14) << "count_ctl_Y"
                  << "\n";
        std::cout << std::string(62, '-') << "\n";

        for (int n = 0; n < 10; ++n) {
            std::string n_str = (n == 9) ? ">=9" : std::to_string(n);
            std::cout << std::left << std::setw(6) << n_str
                      << std::right << std::setw(14) << rest_sig_x[n]
                      << std::setw(14) << rest_ctl_x[n]
                      << std::setw(14) << rest_sig_y[n]
                      << std::setw(14) << rest_ctl_y[n]
                      << "\n";
        }
        std::cout << std::string(62, '-') << "\n\n";

        // High multiplicity summary
        std::cout << "Number of events with nX_sig >= 4 OR nY_sig >= 4: " << n_high_mult_ev << "\n";
        std::cout << "Channel hit counts for these events:\n";
        for (int c = 0; c < 32; ++c) {
            std::cout << c << ": " << ch_high_mult_cnt[c];
            if (c < 31) std::cout << ", ";
        }
        std::cout << "\n\n";
    }

    // =========================================================================
    // PART C: ch29 Amplitude Spectrum, n = 10 p.e.
    // =========================================================================
    std::cout << "=======================================================================================================================================================\n";
    std::cout << "PART C: ch29 Amplitude Spectrum Diagnostics at n = 10 p.e. (Run v8 Good Time Only)\n";
    std::cout << "=======================================================================================================================================================\n";

    const std::vector<int> target_chs = {29, 28, 30, 21, 22, 23};

    // Histograms: binning 1 pe from 0 to 60 pe (60 bins: [0, 60])
    std::vector<TH1D*> h_sig(32, nullptr), h_ctl(32, nullptr), h_diff(32, nullptr);
    for (int c : target_chs) {
        h_sig[c]  = new TH1D(Form("h_sig_ch%d", c),  Form("ch%d Sig;p.e.;Entries", c), 60, 0.0, 60.0);
        h_ctl[c]  = new TH1D(Form("h_ctl_ch%d", c),  Form("ch%d Ctl;p.e.;Entries", c), 60, 0.0, 60.0);
        h_diff[c] = new TH1D(Form("h_diff_ch%d", c), Form("ch%d (Sig - Ctl);p.e.;Entries", c), 60, 0.0, 60.0);
    }

    long long n_ev_ch29_sig = 0;
    long long n_ev_ch29_ctl = 0;
    long long ch29_sig_coinc[32] = {};
    long long ch29_ctl_coinc[32] = {};

    for (Long64_t i = 0; i < n_v8; ++i) {
        bool is_ch29_sig = false;
        bool is_ch29_ctl = false;

        bool hit_s[32] = {};
        bool hit_c[32] = {};

        for (int c = 0; c < 32; ++c) {
            if (c == 3) continue;
            double thr10 = get_thr(c, 10.0);
            float ph_s = kept_events[i].ph_s[c];
            float ph_c = kept_events[i].ph_c[c];

            if (ph_s > thr10) {
                hit_s[c] = true;
                if (c == 29) is_ch29_sig = true;
            }
            if (ph_c > thr10) {
                hit_c[c] = true;
                if (c == 29) is_ch29_ctl = true;
            }
        }

        // Fill spectrum histograms for target channels
        for (int c : target_chs) {
            double thr10 = get_thr(c, 10.0);
            float ph_s = kept_events[i].ph_s[c];
            float ph_c = kept_events[i].ph_c[c];

            if (ph_s > thr10) {
                double pe_s = ((double)ph_s - calib[c].p0) / calib[c].G;
                h_sig[c]->Fill(pe_s);
            }
            if (ph_c > thr10) {
                double pe_c = ((double)ph_c - calib[c].p0) / calib[c].G;
                h_ctl[c]->Fill(pe_c);
            }
        }

        // Coincidence fractions for ch29
        if (is_ch29_sig) {
            n_ev_ch29_sig++;
            for (int d = 0; d < 32; ++d) {
                if (hit_s[d]) ch29_sig_coinc[d]++;
            }
        }
        if (is_ch29_ctl) {
            n_ev_ch29_ctl++;
            for (int d = 0; d < 32; ++d) {
                if (hit_c[d]) ch29_ctl_coinc[d]++;
            }
        }
    }

    // Compute diff and print bin contents
    for (int c : target_chs) {
        h_diff[c]->Add(h_sig[c], 1.0);
        h_diff[c]->Add(h_ctl[c], -1.0);

        long long entries_sig = (long long)h_sig[c]->GetEntries();
        long long entries_ctl = (long long)h_ctl[c]->GetEntries();

        std::cout << "Channel c = " << c << " (entries_sig = " << entries_sig << ", entries_ctl = " << entries_ctl << "):\n";
        std::cout << "bin_contents_sig_minus_ctl: ";
        for (int b = 1; b <= 60; ++b) {
            double val = h_diff[c]->GetBinContent(b);
            std::cout << (long long)val;
            if (b < 60) std::cout << " ";
        }
        std::cout << "\n\n";
    }

    // Print ch29 coincidence fractions
    std::cout << "Channel 29 Signal Window Coincidence Fractions (N_sig = " << n_ev_ch29_sig << "):\n";
    for (int d = 0; d < 32; ++d) {
        double frac = (n_ev_ch29_sig > 0) ? ((double)ch29_sig_coinc[d] / (double)n_ev_ch29_sig) : 0.0;
        std::cout << d << ": " << std::fixed << std::setprecision(4) << frac;
        if (d < 31) std::cout << ", ";
    }
    std::cout << "\n\n";

    std::cout << "Channel 29 Control Window Coincidence Fractions (N_ctl = " << n_ev_ch29_ctl << "):\n";
    for (int d = 0; d < 32; ++d) {
        double frac = (n_ev_ch29_ctl > 0) ? ((double)ch29_ctl_coinc[d] / (double)n_ev_ch29_ctl) : 0.0;
        std::cout << d << ": " << std::fixed << std::setprecision(4) << frac;
        if (d < 31) std::cout << ", ";
    }
    std::cout << "\n\n";

    // Overlay PNG normalized to unit area
    TCanvas c1("c1", "ch29 Spectrum Diagnostics", 900, 650);
    c1.SetRightMargin(0.08);
    c1.SetLeftMargin(0.12);
    c1.SetTopMargin(0.08);
    c1.SetBottomMargin(0.12);
    c1.SetGridx();
    c1.SetGridy();

    TLegend leg(0.65, 0.60, 0.90, 0.90);
    leg.SetBorderSize(1);
    leg.SetFillColor(kWhite);

    int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 2, kCyan + 2};

    double max_val = 0.0;
    for (size_t ic = 0; ic < target_chs.size(); ++ic) {
        int c = target_chs[ic];
        double integ = h_diff[c]->Integral();
        if (integ > 0) {
            h_diff[c]->Scale(1.0 / integ);
        }
        double m = h_diff[c]->GetMaximum();
        if (m > max_val) max_val = m;
    }

    for (size_t ic = 0; ic < target_chs.size(); ++ic) {
        int c = target_chs[ic];
        h_diff[c]->SetLineColor(colors[ic]);
        h_diff[c]->SetLineWidth(2);
        h_diff[c]->SetTitle("Normalized (Signal - Control) Spectra at 10 p.e.;Pulse Height [p.e.];Probability Density [1/p.e.]");
        h_diff[c]->GetYaxis()->SetRangeUser(-0.01, max_val * 1.15);
        h_diff[c]->GetXaxis()->SetRangeUser(0.0, 60.0);

        if (ic == 0) {
            h_diff[c]->Draw("HIST");
        } else {
            h_diff[c]->Draw("HIST SAME");
        }
        leg.AddEntry(h_diff[c], Form("Channel %d", c), "l");
    }
    leg.Draw();

    c1.SaveAs("outputs/fix17_spectra_ch29.png");
    std::cout << "Saved overlaid spectra figure: outputs/fix17_spectra_ch29.png\n\n";
}

