// analysis/macros/hod2018_amplitude.C
// Step 2 of ANALYSIS-FIX-8: hod2018 amplitude and saturation analysis on run v8
//
// 2a Lowest waveform value (mV) at imin over hod2018 hits in W2 (and control window),
//    and count of hits within 1 mV of lowest value.
// 2b Pulse height spectra for hod2018 (summed over 31 live channels): 25 mV bins from 50 mV to
//    ADC limit + 50 mV, for W2 and control window.
// 2c Q spectra (summed over live channels): 1000 mV·ns bins from 0 to 80000 plus overflow, W2 and control.
// 2d ph-cut scan: for each threshold T in {50, 100, 200, 300, 500, 700, 900} mV and for "saturated"
//    (within 1 mV of lowest value):
//    - >= 1 hit in 18X (ch 0-15) with ph > T in W2
//    - >= 1 hit in 18Y (ch 16-31) with ph > T in W2
//    - coincident 18X and 18Y (both > T in W2)
//    - accidental (control window, same formula)
//    - excess = coincident - accidental

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <filesystem>
#include <ctime>
#include <algorithm>

namespace fs = std::filesystem;

#include "../src/hodo_consts.h"
#include "../src/hodo_mapping.h"
#include "../src/hodo_paths.h"

constexpr double kPrePulse = 2500.0; // ns

struct HitData {
    int ch = -1;
    double vmin = 0.0;
    double ph = 0.0;
    double Q = 0.0;
};

struct EventData {
    std::vector<HitData> w2;
    std::vector<HitData> ctl;
};

void hod2018_amplitude(const char* input_file = "lectura_vx1740_ext_nim_1M_v8.root") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    std::cout << "=== hod2018_amplitude ===\n";
    std::cout << "Input file: " << input_file << "\n";
    std::cout << "Signal window W2:  [" << kSigTlo << ", " << kSigThi << "] ns\n";
    std::cout << "Control window:    [" << (kSigTlo + kCtlShift) << ", " << (kSigThi + kCtlShift) << "] ns\n";

    const std::string resolved_input = hodo::resolveRawFile(input_file);
    TFile* fin = TFile::Open(resolved_input.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input file " << resolved_input << "\n";
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(fin->Get("data_dig"));
    if (!tree) {
        std::cerr << "ERROR: data_dig tree not found in " << input_file << "\n";
        return;
    }

    std::array<std::vector<float>*, 32> cw{};
    std::vector<float>* tw = nullptr;
    tree->SetBranchStatus("*", 0);
    tree->SetCacheSize(256LL << 20);

    for (int c = 0; c < 32; ++c) {
        char bn[16]; std::snprintf(bn, sizeof(bn), "Ch%d", c);
        if (tree->GetBranch(bn)) {
            tree->SetBranchStatus(bn, 1);
            tree->SetBranchAddress(bn, &cw[c]);
        }
    }
    tree->SetBranchStatus("time", 1);
    tree->SetBranchAddress("time", &tw);

    const Long64_t nEnt = tree->GetEntries();
    std::cout << "Processing " << nEnt << " events...\n";

    std::vector<EventData> ev_data(nEnt);
    double lowest_vmin_w2 = 1e30;
    double lowest_vmin_ctl = 1e30;

    bool checked_step = false;

    // Pass 1: Extract hits in W2 and Control
    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (iev % 30000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
        if (!tw || tw->empty()) continue;
        const auto& tv = *tw;
        const int n_samp = (int)tv.size();
        if (!checked_step && n_samp > 1) {
            checked_step = true;
            if (std::fabs((tv[1] - tv[0]) - kDt_ns) > 1e-3) {
                std::fprintf(stderr, "Sampling step %.4f ns != kDt_ns %.4f ns; aborting\n", tv[1] - tv[0], kDt_ns);
                std::exit(1);
            }
        }

        for (int c = 0; c < 32; ++c) {
            if (!cw[c] || cw[c]->empty()) continue;
            const auto& wf = *cw[c];

            // Pedestal from pre-pulse
            double s_ped = 0.0; int n_ped = 0;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] < kPrePulse) { s_ped += wf[i]; n_ped++; }
            }
            if (n_ped < 2) continue;
            const double ped = s_ped / n_ped;

            // W2: [kSigTlo, kSigThi]
            double vmin_w2 = 1e30; int imin_w2 = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= kSigTlo && tv[i] <= kSigThi && wf[i] < vmin_w2) {
                    vmin_w2 = wf[i]; imin_w2 = i;
                }
            }
            if (imin_w2 >= 0) {
                double ph_w2 = ped - vmin_w2;
                if (ph_w2 > 50.0) {
                    // Local pedestal for Q
                    double ped_q = ped;
                    if (imin_w2 >= 12 && imin_w2 - 4 < n_samp) {
                        double s_loc = 0.0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin_w2 + j];
                        ped_q = s_loc / 9.0;
                    }
                    double Q = 0.0;
                    double tilo = tv[imin_w2] - kNpre * kDt_ns;
                    double tihi = tv[imin_w2] + kNpost18 * kDt_ns;
                    for (int i = 0; i < n_samp; ++i) {
                        if (tv[i] >= tilo && tv[i] <= tihi) Q += (ped_q - wf[i]) * kDt_ns;
                    }

                    ev_data[iev].w2.push_back({c, vmin_w2, ph_w2, Q});
                    if (vmin_w2 < lowest_vmin_w2) lowest_vmin_w2 = vmin_w2;
                }
            }

            // Control window: [kSigTlo + kCtlShift, kSigThi + kCtlShift]
            double vmin_ctl = 1e30; int imin_ctl = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= (kSigTlo + kCtlShift) && tv[i] <= (kSigThi + kCtlShift) && wf[i] < vmin_ctl) {
                    vmin_ctl = wf[i]; imin_ctl = i;
                }
            }
            if (imin_ctl >= 0) {
                double ph_ctl = ped - vmin_ctl;
                if (ph_ctl > 50.0) {
                    // Local pedestal for Q
                    double ped_q_c = ped;
                    if (imin_ctl >= 12 && imin_ctl - 4 < n_samp) {
                        double s_loc = 0.0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin_ctl + j];
                        ped_q_c = s_loc / 9.0;
                    }
                    double Q_ctl = 0.0;
                    double tilo_c = tv[imin_ctl] - kNpre * kDt_ns;
                    double tihi_c = tv[imin_ctl] + kNpost18 * kDt_ns;
                    for (int i = 0; i < n_samp; ++i) {
                        if (tv[i] >= tilo_c && tv[i] <= tihi_c) Q_ctl += (ped_q_c - wf[i]) * kDt_ns;
                    }

                    ev_data[iev].ctl.push_back({c, vmin_ctl, ph_ctl, Q_ctl});
                    if (vmin_ctl < lowest_vmin_ctl) lowest_vmin_ctl = vmin_ctl;
                }
            }
        }
    }
    std::cout << "  " << nEnt << "/" << nEnt << " done.\n\n";
    fin->Close();

    // -------------------------------------------------------------------------
    // 2a. Lowest waveform values and counts within 1 mV
    // -------------------------------------------------------------------------
    long long n_within_1mv_w2 = 0;
    long long n_within_1mv_ctl = 0;

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        for (const auto& h : ev_data[iev].w2) {
            if (h.vmin <= lowest_vmin_w2 + 1.0) n_within_1mv_w2++;
        }
        for (const auto& h : ev_data[iev].ctl) {
            if (h.vmin <= lowest_vmin_ctl + 1.0) n_within_1mv_ctl++;
        }
    }

    std::cout << "================================================================================\n";
    std::cout << "2a. Lowest Waveform Value at imin over hod2018 Hits (ph > 50 mV)\n";
    std::cout << "================================================================================\n";
    std::cout << "Window W2 [" << kSigTlo << ", " << kSigThi << "] ns:\n";
    std::cout << "  Lowest waveform value:       " << std::fixed << std::setprecision(2) << lowest_vmin_w2 << " mV\n";
    std::cout << "  Hits within 1.0 mV of lowest: " << n_within_1mv_w2 << " hits (vmin <= " << (lowest_vmin_w2 + 1.0) << " mV)\n";
    std::cout << "Control window [" << (kSigTlo + kCtlShift) << ", " << (kSigThi + kCtlShift) << "] ns:\n";
    std::cout << "  Lowest waveform value:       " << lowest_vmin_ctl << " mV\n";
    std::cout << "  Hits within 1.0 mV of lowest: " << n_within_1mv_ctl << " hits (vmin <= " << (lowest_vmin_ctl + 1.0) << " mV)\n\n";

    // -------------------------------------------------------------------------
    // 2b. Pulse Height Spectra (25 mV bins from 50 mV to ADC limit + 50 mV)
    // ADC limit: 1000 mV pulse height -> bins up to 1050 mV (40 bins) + overflow
    // -------------------------------------------------------------------------
    constexpr int kNPBins = 40; // 50 to 1050 mV
    std::vector<long long> ph_bins_w2(kNPBins, 0);
    std::vector<long long> ph_bins_ctl(kNPBins, 0);
    long long ph_overflow_w2 = 0;
    long long ph_overflow_ctl = 0;

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        for (const auto& h : ev_data[iev].w2) {
            int b = static_cast<int>((h.ph - 50.0) / 25.0);
            if (b >= 0 && b < kNPBins) ph_bins_w2[b]++;
            else if (b >= kNPBins) ph_overflow_w2++;
        }
        for (const auto& h : ev_data[iev].ctl) {
            int b = static_cast<int>((h.ph - 50.0) / 25.0);
            if (b >= 0 && b < kNPBins) ph_bins_ctl[b]++;
            else if (b >= kNPBins) ph_overflow_ctl++;
        }
    }

    std::cout << "================================================================================\n";
    std::cout << "2b. Pulse Height Spectra for hod2018 (Summed over Live Channels)\n";
    std::cout << "================================================================================\n";
    std::cout << std::setw(18) << std::left << "ph_bin [mV]"
              << std::setw(12) << std::right << "W2"
              << std::setw(12) << "Control"
              << std::setw(14) << "W2 - Control" << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    for (int b = 0; b < kNPBins; ++b) {
        double lo = 50.0 + 25.0 * b;
        double hi = lo + 25.0;
        char bname[32];
        std::snprintf(bname, sizeof(bname), "[%.0f, %.0f)", lo, hi);
        long long diff = ph_bins_w2[b] - ph_bins_ctl[b];
        std::cout << std::setw(18) << std::left << bname
                  << std::setw(12) << std::right << ph_bins_w2[b]
                  << std::setw(12) << ph_bins_ctl[b]
                  << std::setw(14) << diff << "\n";
    }
    std::cout << std::setw(18) << std::left << ">= 1050 (overflow)"
              << std::setw(12) << std::right << ph_overflow_w2
              << std::setw(12) << ph_overflow_ctl
              << std::setw(14) << (ph_overflow_w2 - ph_overflow_ctl) << "\n\n";

    // -------------------------------------------------------------------------
    // 2c. Q Spectra (1000 mV·ns bins from 0 to 80000 plus overflow)
    // -------------------------------------------------------------------------
    constexpr int kNQBins = 80; // 0 to 80000 mV*ns
    std::vector<long long> q_bins_w2(kNQBins, 0);
    std::vector<long long> q_bins_ctl(kNQBins, 0);
    long long q_underflow_w2 = 0, q_underflow_ctl = 0;
    long long q_overflow_w2 = 0, q_overflow_ctl = 0;

    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        for (const auto& h : ev_data[iev].w2) {
            if (h.Q < 0.0) q_underflow_w2++;
            else {
                int b = static_cast<int>(h.Q / 1000.0);
                if (b >= 0 && b < kNQBins) q_bins_w2[b]++;
                else if (b >= kNQBins) q_overflow_w2++;
            }
        }
        for (const auto& h : ev_data[iev].ctl) {
            if (h.Q < 0.0) q_underflow_ctl++;
            else {
                int b = static_cast<int>(h.Q / 1000.0);
                if (b >= 0 && b < kNQBins) q_bins_ctl[b]++;
                else if (b >= kNQBins) q_overflow_ctl++;
            }
        }
    }

    std::cout << "================================================================================\n";
    std::cout << "2c. Q Spectra for hod2018 (Summed over Live Channels) [mV·ns]\n";
    std::cout << "================================================================================\n";
    std::cout << std::setw(20) << std::left << "Q_bin [mV·ns]"
              << std::setw(12) << std::right << "W2"
              << std::setw(12) << "Control"
              << std::setw(14) << "W2 - Control" << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    if (q_underflow_w2 > 0 || q_underflow_ctl > 0) {
        std::cout << std::setw(20) << std::left << "< 0 (underflow)"
                  << std::setw(12) << std::right << q_underflow_w2
                  << std::setw(12) << q_underflow_ctl
                  << std::setw(14) << (q_underflow_w2 - q_underflow_ctl) << "\n";
    }
    for (int b = 0; b < kNQBins; ++b) {
        double lo = 1000.0 * b;
        double hi = lo + 1000.0;
        char bname[32];
        std::snprintf(bname, sizeof(bname), "[%.0f, %.0f)", lo, hi);
        long long diff = q_bins_w2[b] - q_bins_ctl[b];
        std::cout << std::setw(20) << std::left << bname
                  << std::setw(12) << std::right << q_bins_w2[b]
                  << std::setw(12) << q_bins_ctl[b]
                  << std::setw(14) << diff << "\n";
    }
    std::cout << std::setw(20) << std::left << ">= 80000 (overflow)"
              << std::setw(12) << std::right << q_overflow_w2
              << std::setw(12) << q_overflow_ctl
              << std::setw(14) << (q_overflow_w2 - q_overflow_ctl) << "\n\n";

    // -------------------------------------------------------------------------
    // 2d. ph-cut scan
    // Thresholds: {50, 100, 200, 300, 500, 700, 900} mV and "saturated" (vmin <= lowest + 1 mV)
    // -------------------------------------------------------------------------
    const std::vector<double> thresholds = {50.0, 100.0, 200.0, 300.0, 500.0, 700.0, 900.0};

    struct ScanRow {
        std::string label;
        long long n_18x_w2 = 0;
        long long n_18y_w2 = 0;
        long long coinc_w2 = 0;
        long long n_18x_ctl = 0;
        long long n_18y_ctl = 0;
        long long coinc_ctl = 0;
        long long excess = 0;
    };
    std::vector<ScanRow> scan_rows;

    for (double T : thresholds) {
        ScanRow row;
        row.label = Form("%.0f mV", T);
        for (Long64_t iev = 0; iev < nEnt; ++iev) {
            bool has_x_w2 = false, has_y_w2 = false;
            for (const auto& h : ev_data[iev].w2) {
                if (h.ph > T) {
                    if (h.ch >= 0 && h.ch <= 15) has_x_w2 = true;
                    else if (h.ch >= 16 && h.ch <= 31) has_y_w2 = true;
                }
            }
            if (has_x_w2) row.n_18x_w2++;
            if (has_y_w2) row.n_18y_w2++;
            if (has_x_w2 && has_y_w2) row.coinc_w2++;

            bool has_x_ctl = false, has_y_ctl = false;
            for (const auto& h : ev_data[iev].ctl) {
                if (h.ph > T) {
                    if (h.ch >= 0 && h.ch <= 15) has_x_ctl = true;
                    else if (h.ch >= 16 && h.ch <= 31) has_y_ctl = true;
                }
            }
            if (has_x_ctl) row.n_18x_ctl++;
            if (has_y_ctl) row.n_18y_ctl++;
            if (has_x_ctl && has_y_ctl) row.coinc_ctl++;
        }
        row.excess = row.coinc_w2 - row.coinc_ctl;
        scan_rows.push_back(row);
    }

    // Saturated cut (vmin <= lowest_vmin_w2 + 1.0 mV)
    {
        ScanRow row;
        row.label = "saturated";
        for (Long64_t iev = 0; iev < nEnt; ++iev) {
            bool has_x_w2 = false, has_y_w2 = false;
            for (const auto& h : ev_data[iev].w2) {
                if (h.vmin <= lowest_vmin_w2 + 1.0) {
                    if (h.ch >= 0 && h.ch <= 15) has_x_w2 = true;
                    else if (h.ch >= 16 && h.ch <= 31) has_y_w2 = true;
                }
            }
            if (has_x_w2) row.n_18x_w2++;
            if (has_y_w2) row.n_18y_w2++;
            if (has_x_w2 && has_y_w2) row.coinc_w2++;

            bool has_x_ctl = false, has_y_ctl = false;
            for (const auto& h : ev_data[iev].ctl) {
                if (h.vmin <= lowest_vmin_w2 + 1.0) {
                    if (h.ch >= 0 && h.ch <= 15) has_x_ctl = true;
                    else if (h.ch >= 16 && h.ch <= 31) has_y_ctl = true;
                }
            }
            if (has_x_ctl) row.n_18x_ctl++;
            if (has_y_ctl) row.n_18y_ctl++;
            if (has_x_ctl && has_y_ctl) row.coinc_ctl++;
        }
        row.excess = row.coinc_w2 - row.coinc_ctl;
        scan_rows.push_back(row);
    }

    std::cout << "================================================================================\n";
    std::cout << "2d. ph-Cut Scan in Window W2 [" << kSigTlo << ", " << kSigThi << "] ns\n";
    std::cout << "================================================================================\n";
    std::cout << std::setw(14) << std::left << "Cut (T)"
              << std::setw(14) << std::right << ">=1 18X (W2)"
              << std::setw(14) << ">=1 18Y (W2)"
              << std::setw(14) << "Coinc 18X&18Y"
              << std::setw(14) << "Accidental"
              << std::setw(12) << "Excess" << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    for (const auto& r : scan_rows) {
        std::cout << std::setw(14) << std::left << r.label
                  << std::setw(14) << std::right << r.n_18x_w2
                  << std::setw(14) << r.n_18y_w2
                  << std::setw(14) << r.coinc_w2
                  << std::setw(14) << r.coinc_ctl
                  << std::setw(12) << r.excess << "\n";
    }
    std::cout << "\nDetailed plane counts in control window [" << (kSigTlo + kCtlShift) << ", " << (kSigThi + kCtlShift) << "] ns:\n";
    std::cout << std::setw(14) << std::left << "Cut (T)"
              << std::setw(14) << std::right << ">=1 18X (Ctl)"
              << std::setw(14) << ">=1 18Y (Ctl)"
              << std::setw(14) << "Coinc 18X&18Y" << "\n";
    std::cout << "------------------------------------------------------------\n";
    for (const auto& r : scan_rows) {
        std::cout << std::setw(14) << std::left << r.label
                  << std::setw(14) << std::right << r.n_18x_ctl
                  << std::setw(14) << r.n_18y_ctl
                  << std::setw(14) << r.coinc_ctl << "\n";
    }
    std::cout << "================================================================================\n";
}

