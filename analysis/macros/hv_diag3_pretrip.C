// analysis/macros/hv_diag3_pretrip.C
// Task HV-DIAG-3 Part 2: pre-trip behaviour of all channels (v8, v10)

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TAxis.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>

#include "../src/hodo_paths.h"

struct TripDef {
    std::string name;
    std::string run;
    std::string filename;
    std::vector<int> group_channels;
    Long_t approx_start_ts; // microsecond timestamp to start search
    Long64_t start_entry;   // approximate entry to start scan
};

struct ChannelCalib {
    double p0 = 4.3;
    double G  = 4.0;
};

void hv_diag3_pretrip() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // 1. Load global calibration from outputs/pe_calibration_v8.txt
    std::map<int, ChannelCalib> calib_global;
    std::ifstream fcal("outputs/pe_calibration_v8.txt");
    if (fcal.is_open()) {
        std::string line;
        while (std::getline(fcal, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            int ch; std::string pl; double p0, G, s0, s1, chi2; int ndf, status;
            if (iss >> ch >> pl >> p0 >> G >> s0 >> s1 >> chi2 >> ndf >> status) {
                calib_global[ch] = {p0, G};
            }
        }
        fcal.close();
    }

    // 2. Load hourly calibration from outputs/hv_diag3_gain_vs_time.csv if available
    // key: run + "_" + date + "_" + hour + "_" + ch -> G
    std::map<std::string, double> hourly_G_map;
    std::map<std::string, double> hourly_p0_map;
    std::ifstream fhour("outputs/hv_diag3_gain_vs_time.csv");
    if (fhour.is_open()) {
        std::string line;
        std::getline(fhour, line); // header
        while (std::getline(fhour, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string run, date, hour;
            int ch, status;
            Long64_t n_ev;
            double T_live, G, p0, chi2ndf, N1pe, rate_05, Iproxy, trig_rate;
            char comma;
            if (std::getline(iss, run, ',') &&
                std::getline(iss, date, ',') &&
                std::getline(iss, hour, ',') &&
                (iss >> ch >> comma >> n_ev >> comma >> T_live >> comma >> G >> comma
                     >> p0 >> comma >> chi2ndf >> comma >> status >> comma >> N1pe >> comma
                     >> rate_05 >> comma >> Iproxy >> comma >> trig_rate)) {
                if (status == 0 && chi2ndf < 3.0 && G > 0.5 && G < 10.0) {
                    std::string key = run + "_" + date + "_" + hour + "_" + std::to_string(ch);
                    hourly_G_map[key] = G;
                    hourly_p0_map[key] = p0;
                }
            }
        }
        fhour.close();
    }

    // Define the three trips
    // Group 1: 0, 22, 31, 43, 47, 54, 58, 63
    // Group 2: 6, 20, 46, 50, 56, 59, 60, 62
    // Group 6: 1, 16, 17, 34, 53, 55, 57
    std::vector<TripDef> trips = {
        {"T1", "v8",  "lectura_vx1740_ext_nim_1M_v8.root",  {6, 20, 46, 50, 56, 59, 60, 62}, 0, 91000},
        {"T2", "v10", "lectura_vx1740_ext_nim_1M_v10.root", {0, 22, 31, 43, 47, 54, 58, 63, 6, 20, 46, 50, 56, 59, 60, 62}, 0, 30000},
        {"T3", "v10", "lectura_vx1740_ext_nim_1M_v10.root", {1, 16, 17, 34, 53, 55, 57}, 0, 128000}
    };

    // Group mapping for 2c
    auto get_group_str = [](int c) -> std::string {
        static const std::vector<int> g1 = {0, 22, 31, 43, 47, 54, 58, 63};
        static const std::vector<int> g2 = {6, 20, 46, 50, 56, 59, 60, 62};
        static const std::vector<int> g6 = {1, 16, 17, 34, 53, 55, 57};
        if (std::find(g1.begin(), g1.end(), c) != g1.end()) return "1";
        if (std::find(g2.begin(), g2.end(), c) != g2.end()) return "2";
        if (std::find(g6.begin(), g6.end(), c) != g6.end()) return "6";
        return "-";
    };

    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // -04:00 offset

    // Data structures for 2e PNG
    struct TripPlotData {
        std::vector<double> t_min;
        std::vector<double> rate_tripping;
        std::vector<double> rate_other_18;
    };
    std::map<std::string, TripPlotData> plot_data;

    std::cout << "========================================================================================\n";
    std::cout << "TASK HV-DIAG-3 PART 2: Behaviour Before Each HV Trip\n";
    std::cout << "========================================================================================\n\n";

    for (const auto& tr : trips) {
        std::string resolved_path = hodo::resolveRawFile(tr.filename);
        std::unique_ptr<TFile> fin(TFile::Open(resolved_path.c_str(), "READ"));
        if (!fin || fin->IsZombie()) {
            std::cerr << "Error opening " << resolved_path << "\n";
            continue;
        }
        TTree* tree = dynamic_cast<TTree*>(fin->Get("data_dig"));
        if (!tree) {
            std::cerr << "Error finding data_dig in " << resolved_path << "\n";
            continue;
        }

        Long64_t nEnt = tree->GetEntries();

        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("time", 1);
        tree->SetBranchStatus("timestamp", 1);

        std::vector<float>* tw = nullptr;
        Long_t ts = 0;
        tree->SetBranchAddress("time", &tw);
        tree->SetBranchAddress("timestamp", &ts);

        std::vector<float>* cw[64] = {};
        for (int c = 0; c < 64; ++c) {
            tree->SetBranchStatus(Form("Ch%d", c), 1);
            tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
        }
        tree->SetCacheSize(256LL << 20);

        // ---------------------------------------------------------------------
        // 2a. Find onset time
        // ---------------------------------------------------------------------
        int consec_zeros = 0;
        Long64_t last_pulse_entry = -1;
        Long_t last_pulse_ts = -1;
        Long64_t onset_entry = -1;
        Long_t onset_ts = -1;
        Long64_t potential_onset_entry = -1;
        Long_t potential_onset_ts = -1;

        for (Long64_t iev = tr.start_entry; iev < nEnt; ++iev) {
            tree->GetEntry(iev);
            if (!tw || tw->empty()) continue;
            const auto& tv = *tw;
            const int n_samp = (int)tv.size();

            bool group_has_pulse = false;
            for (int c : tr.group_channels) {
                if (!cw[c] || (int)cw[c]->size() != n_samp) continue;
                const auto& wf = *cw[c];

                double vmin = 1e30; int imin = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= 100.0 && tv[i] <= 3480.0 && wf[i] < vmin) {
                        vmin = wf[i]; imin = i;
                    }
                }
                if (imin >= 0) {
                    double ped_loc = 0.0;
                    if (imin >= 12 && imin - 4 < n_samp) {
                        double s_loc = 0.0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin + j];
                        ped_loc = s_loc / 9.0;
                    } else {
                        double s_pre = 0.0; int n_pre = 0;
                        for (int i = 0; i < n_samp; ++i) {
                            if (tv[i] < 2500.0) { s_pre += wf[i]; n_pre++; }
                        }
                        ped_loc = (n_pre > 0) ? (s_pre / n_pre) : 0.0;
                    }
                    double ph = ped_loc - vmin;

                    double thr = 3.0; // default for hod2019
                    if (c < 32) {
                        double p0 = calib_global.count(c) ? calib_global[c].p0 : 4.3;
                        double G  = calib_global.count(c) ? calib_global[c].G  : 4.0;
                        thr = p0 + 0.5 * G;
                    }
                    if (ph > thr) {
                        group_has_pulse = true;
                        break;
                    }
                }
            }

            if (group_has_pulse) {
                consec_zeros = 0;
                last_pulse_entry = iev;
                last_pulse_ts = ts;
            } else {
                if (consec_zeros == 0) {
                    potential_onset_entry = iev;
                    potential_onset_ts = ts;
                }
                consec_zeros++;
                if (consec_zeros == 100) {
                    onset_entry = potential_onset_entry;
                    onset_ts = potential_onset_ts;
                    break;
                }
            }
        }

        std::cout << "========================================================================================\n";
        std::cout << "PART 2 (a) Onset Timestamps for Trip " << tr.name << " (" << tr.run << ")\n";
        std::cout << "========================================================================================\n";
        std::cout << "Onset timestamp:              " << onset_ts << "  (event " << onset_entry << ")\n";
        std::cout << "Last pulse timestamp:         " << last_pulse_ts << "  (event " << last_pulse_entry << ")\n";

        double dt_sec_onset = boot_ref_epoch + ((double)onset_ts) * 1e-6 + tz_offset;
        time_t t_local_onset = (time_t)dt_sec_onset;
        struct tm tm_onset; gmtime_r(&t_local_onset, &tm_onset);
        char onset_buf[64];
        std::strftime(onset_buf, sizeof(onset_buf), "%Y-%m-%d %H:%M:%S", &tm_onset);

        double dt_sec_last = boot_ref_epoch + ((double)last_pulse_ts) * 1e-6 + tz_offset;
        time_t t_local_last = (time_t)dt_sec_last;
        struct tm tm_last; gmtime_r(&t_local_last, &tm_last);
        char last_buf[64];
        std::strftime(last_buf, sizeof(last_buf), "%Y-%m-%d %H:%M:%S", &tm_last);

        std::cout << "Onset local time:             " << onset_buf << "\n";
        std::cout << "Last pulse local time:        " << last_buf << "\n";
        std::cout << "========================================================================================\n\n";

        // ---------------------------------------------------------------------
        // 2b. 2-minute bins from onset - 120 min to onset + 10 min
        // 65 bins: [ -7200 s, +600 s ]
        // ---------------------------------------------------------------------
        constexpr int kNBins = 65; // 130 minutes / 2 minutes
        struct BinData {
            Long64_t n_events = 0;
            long long pulse_count[64] = {};
            double amp_sum[64] = {};
        };
        std::vector<BinData> bins(kNBins);

        // Find scan entry range covering [onset - 7200 s, onset + 600 s]
        Long64_t start_ev = onset_entry;
        while (start_ev > 0) {
            tree->GetEntry(start_ev);
            if (((double)(onset_ts - ts)) * 1e-6 > 7300.0) break;
            start_ev = std::max((Long64_t)0, start_ev - 500);
        }

        Long64_t end_ev = onset_entry;
        while (end_ev < nEnt - 1) {
            tree->GetEntry(end_ev);
            if (((double)(ts - onset_ts)) * 1e-6 > 700.0) break;
            end_ev = std::min(nEnt - 1, end_ev + 500);
        }

        for (Long64_t iev = start_ev; iev <= end_ev; ++iev) {
            tree->GetEntry(iev);
            if (!tw || tw->empty()) continue;

            double dt_s = ((double)(ts - onset_ts)) * 1e-6; // seconds relative to onset
            int b = (int)std::floor((dt_s + 7200.0) / 120.0);
            if (b < 0 || b >= kNBins) continue;

            bins[b].n_events++;

            // Date and hour key for hourly calibration lookup
            double t_sec = boot_ref_epoch + ((double)ts) * 1e-6 + tz_offset;
            time_t t_local = (time_t)t_sec;
            struct tm tm_val; gmtime_r(&t_local, &tm_val);
            char dh_buf[64];
            std::snprintf(dh_buf, sizeof(dh_buf), "%s_%04d-%02d-%02d_%02d:00",
                          tr.run.c_str(), tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday, tm_val.tm_hour);
            std::string hour_key(dh_buf);

            const auto& tv = *tw;
            const int n_samp = (int)tv.size();

            for (int c = 0; c < 64; ++c) {
                if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
                const auto& wf = *cw[c];

                double vmin = 1e30; int imin = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= 100.0 && tv[i] <= 3480.0 && wf[i] < vmin) {
                        vmin = wf[i]; imin = i;
                    }
                }
                if (imin >= 0) {
                    double ped_loc = 0.0;
                    if (imin >= 12 && imin - 4 < n_samp) {
                        double s_loc = 0.0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin + j];
                        ped_loc = s_loc / 9.0;
                    } else {
                        double s_pre = 0.0; int n_pre = 0;
                        for (int i = 0; i < n_samp; ++i) {
                            if (tv[i] < 2500.0) { s_pre += wf[i]; n_pre++; }
                        }
                        ped_loc = (n_pre > 0) ? (s_pre / n_pre) : 0.0;
                    }
                    double ph = ped_loc - vmin;

                    if (c < 32) {
                        // hod2018
                        std::string ch_key = hour_key + "_" + std::to_string(c);
                        double G  = hourly_G_map.count(ch_key) ? hourly_G_map[ch_key] :
                                    (calib_global.count(c) ? calib_global[c].G : 4.0);
                        double p0 = hourly_p0_map.count(ch_key) ? hourly_p0_map[ch_key] :
                                    (calib_global.count(c) ? calib_global[c].p0 : 4.3);
                        double thr = p0 + 0.5 * G;
                        if (ph > thr) {
                            bins[b].pulse_count[c]++;
                            bins[b].amp_sum[c] += (ph - p0);
                        }
                    } else {
                        // hod2019
                        if (ph > 3.0) {
                            bins[b].pulse_count[c]++;
                            bins[b].amp_sum[c] += ph;
                        }
                    }
                }
            }
        }

        // Save CSV for 2b
        std::string out_csv_path = Form("outputs/hv_diag3_pretrip_%s.csv", tr.name.c_str());
        std::ofstream out_csv(out_csv_path);
        out_csv << "bin,t_start_min,t_end_min,ch,pulse_count,amp_sum,live_time_s\n";
        for (int b = 0; b < kNBins; ++b) {
            double t_start_min = -120.0 + b * 2.0;
            double t_end_min   = -120.0 + (b + 1) * 2.0;
            double live_time   = (double)bins[b].n_events * 3.38e-6;
            for (int c = 0; c < 64; ++c) {
                out_csv << b << ","
                        << std::fixed << std::setprecision(1) << t_start_min << ","
                        << t_end_min << ","
                        << c << ","
                        << bins[b].pulse_count[c] << ","
                        << std::setprecision(4) << bins[b].amp_sum[c] << ","
                        << std::setprecision(6) << live_time << "\n";
            }
        }
        out_csv.close();
        std::cout << "Saved: " << out_csv_path << "\n";

        // ---------------------------------------------------------------------
        // 2c & 2d. Early vs Late comparison
        // Early: [onset-120, onset-60] min -> bins 0 to 29 (30 bins, 60 minutes)
        // Late:  [onset-10, onset] min    -> bins 55 to 59 (5 bins, 10 minutes)
        // ---------------------------------------------------------------------
        std::cout << "\n========================================================================================================================\n";
        std::cout << "PART 2 (c) Early vs Late Comparison for Trip " << tr.name << " (" << tr.run << ")\n";
        std::cout << "========================================================================================================================\n";
        std::cout << std::left << std::setw(6) << "ch"
                  << std::setw(8) << "group"
                  << std::right << std::setw(14) << "rate_early"
                  << std::setw(14) << "err_early"
                  << std::setw(14) << "rate_late"
                  << std::setw(14) << "err_late"
                  << std::setw(12) << "ratio"
                  << std::setw(16) << "Iproxy_early"
                  << std::setw(16) << "Iproxy_late" << "\n";
        std::cout << std::string(108, '-') << "\n";

        std::vector<int> significant_channels;

        for (int c = 0; c < 64; ++c) {
            if (c == 3) continue; // dead channel

            long long count_early = 0, count_late = 0;
            Long64_t ev_early = 0, ev_late = 0;
            double amp_early = 0.0, amp_late = 0.0;

            for (int b = 0; b <= 29; ++b) {
                count_early += bins[b].pulse_count[c];
                amp_early   += bins[b].amp_sum[c];
                ev_early    += bins[b].n_events;
            }
            for (int b = 55; b <= 59; ++b) {
                count_late += bins[b].pulse_count[c];
                amp_late   += bins[b].amp_sum[c];
                ev_late    += bins[b].n_events;
            }

            double T_live_early = (double)ev_early * 3.38e-6;
            double T_live_late  = (double)ev_late  * 3.38e-6;

            double rate_early = (T_live_early > 0.0) ? ((double)count_early / T_live_early) : 0.0;
            double err_early  = (T_live_early > 0.0) ? (std::sqrt((double)count_early) / T_live_early) : 0.0;
            double rate_late  = (T_live_late > 0.0) ? ((double)count_late / T_live_late) : 0.0;
            double err_late   = (T_live_late > 0.0) ? (std::sqrt((double)count_late) / T_live_late) : 0.0;

            double ratio = (rate_early > 0.0) ? (rate_late / rate_early) : 0.0;
            double Iproxy_early = (T_live_early > 0.0) ? (amp_early / T_live_early) : 0.0;
            double Iproxy_late  = (T_live_late > 0.0) ? (amp_late / T_live_late) : 0.0;

            std::cout << std::left << std::setw(6) << c
                      << std::setw(8) << get_group_str(c)
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(14) << rate_early
                      << std::setw(14) << err_early
                      << std::setw(14) << rate_late
                      << std::setw(14) << err_late
                      << std::setprecision(4)
                      << std::setw(12) << ratio
                      << std::setprecision(2)
                      << std::setw(16) << Iproxy_early
                      << std::setw(16) << Iproxy_late << "\n";

            // 2d condition: ratio > 1.5 && (rate_late - rate_early) > 3 * sqrt(err_early^2 + err_late^2)
            double combined_err = std::sqrt(err_early * err_early + err_late * err_late);
            if (ratio > 1.5 && (rate_late - rate_early) > 3.0 * combined_err) {
                significant_channels.push_back(c);
            }
        }
        std::cout << "========================================================================================================================\n\n";

        // ---------------------------------------------------------------------
        // 2d. List of significant channels
        // ---------------------------------------------------------------------
        std::cout << "PART 2 (d) Channels with ratio > 1.5 and rate_late - rate_early > 3*sigma for Trip " << tr.name << ":\n";
        if (significant_channels.empty()) {
            std::cout << "  None\n";
        } else {
            for (int c : significant_channels) {
                std::cout << "  Channel " << c << " (Group " << get_group_str(c) << ")\n";
            }
        }
        std::cout << "----------------------------------------------------------------------------------------\n\n";

        // Collect rate curves for 2e PNG
        TripPlotData tpd;
        for (int b = 0; b < kNBins; ++b) {
            double t_mid = -120.0 + (b + 0.5) * 2.0;
            double live_time = (double)bins[b].n_events * 3.38e-6;

            long long sum_trip = 0;
            for (int c : tr.group_channels) {
                sum_trip += bins[b].pulse_count[c];
            }

            long long sum_other_18 = 0;
            for (int c = 0; c < 32; ++c) {
                if (c == 3) continue;
                if (std::find(tr.group_channels.begin(), tr.group_channels.end(), c) == tr.group_channels.end()) {
                    sum_other_18 += bins[b].pulse_count[c];
                }
            }

            double r_trip = (live_time > 0.0) ? ((double)sum_trip / live_time) : 0.0;
            double r_other = (live_time > 0.0) ? ((double)sum_other_18 / live_time) : 0.0;

            tpd.t_min.push_back(t_mid);
            tpd.rate_tripping.push_back(r_trip);
            tpd.rate_other_18.push_back(r_other);
        }
        plot_data[tr.name] = tpd;
    }

    // -------------------------------------------------------------------------
    // 2e. Save PNG outputs/hv_diag3_pretrip.png
    // 6 pads: 3 rows x 2 columns
    // Left column: Tripping group(s) summed rate
    // Right column: Other hod2018 channels summed rate
    // -------------------------------------------------------------------------
    TCanvas c_pretrip("c_pretrip", "Pre-trip Rates", 1400, 1000);
    c_pretrip.Divide(2, 3);

    std::vector<std::string> trip_names = {"T1", "T2", "T3"};
    std::vector<std::string> trip_labels = {
        "T1 (v8, Group 2)",
        "T2 (v10, Groups 1+2)",
        "T3 (v10, Group 6)"
    };

    for (size_t i = 0; i < trip_names.size(); ++i) {
        std::string tn = trip_names[i];
        const auto& tpd = plot_data[tn];
        int np = (int)tpd.t_min.size();

        // Left pad: Tripping group rate
        c_pretrip.cd((int)i * 2 + 1);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.12);

        TGraph* gr_trip = new TGraph(np, tpd.t_min.data(), tpd.rate_tripping.data());
        gr_trip->SetTitle(Form("%s: Tripping Group(s) Summed Rate;Time relative to onset [min];Rate [Hz]", trip_labels[i].c_str()));
        gr_trip->SetMarkerStyle(20);
        gr_trip->SetMarkerSize(0.7);
        gr_trip->SetMarkerColor(kRed + 1);
        gr_trip->SetLineColor(kRed + 1);
        gr_trip->GetXaxis()->SetRangeUser(-120.0, 10.0);
        gr_trip->GetXaxis()->SetTitleSize(0.045);
        gr_trip->GetYaxis()->SetTitleSize(0.045);
        gr_trip->Draw("ALP");

        // Right pad: Other hod2018 channels rate
        c_pretrip.cd((int)i * 2 + 2);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.12);

        TGraph* gr_other = new TGraph(np, tpd.t_min.data(), tpd.rate_other_18.data());
        gr_other->SetTitle(Form("%s: Other hod2018 Channels Summed Rate;Time relative to onset [min];Rate [Hz]", trip_labels[i].c_str()));
        gr_other->SetMarkerStyle(21);
        gr_other->SetMarkerSize(0.7);
        gr_other->SetMarkerColor(kBlue + 1);
        gr_other->SetLineColor(kBlue + 1);
        gr_other->GetXaxis()->SetRangeUser(-120.0, 10.0);
        gr_other->GetXaxis()->SetTitleSize(0.045);
        gr_other->GetYaxis()->SetTitleSize(0.045);
        gr_other->Draw("ALP");
    }

    c_pretrip.SaveAs("outputs/hv_diag3_pretrip.png");
    std::cout << "Saved figure: outputs/hv_diag3_pretrip.png\n";
}

