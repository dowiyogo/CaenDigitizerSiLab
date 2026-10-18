// analysis/macros/hv_diag3_gain_vs_time.C
// Task HV-DIAG-3 Part 1: hourly gain of hod2018 for runs v8, v9, v10

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TDatime.h>
#include <TAxis.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <cmath>
#include <algorithm>

#include "../src/hodo_paths.h"

// 6-peak model from pe_scale_step2.C
// par[0] = p0, par[1] = G, par[2] = s0, par[3] = s1, par[4..9] = A1..A6
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

struct Interval {
    std::string run;
    std::string start_local;
    std::string end_local;
    std::string status;
    std::vector<int> channels;
    bool all_channels = false;
};

struct HourData {
    std::string run;
    std::string date;
    std::string hour;
    std::string key; // run + "_" + date + "_" + hour
    Long64_t n_events = 0;
    Long_t first_ts = 0;
    Long_t last_ts = 0;
    // For each channel 0..31: vector of pulse amplitudes (ph) in this hour
    std::vector<float> ph[32];
    TH1D* h_ph[32] = {};
};

struct FitRecord {
    std::string key;
    std::string run;
    std::string date;
    std::string hour;
    int ch = -1;
    Long64_t n_events = 0;
    double T_live = 0.0;
    double G = 0.0;
    double p0 = 0.0;
    double chi2ndf = 0.0;
    int status = -1;
    double N1pe = 0.0;
    double rate_05 = 0.0;
    double Iproxy = 0.0;
    double trig_rate = 0.0;
    bool is_ok = false;
};

static double compute_median(std::vector<double>& v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

void hv_diag3_gain_vs_time() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // 1. Read good_intervals_2026.csv
    std::vector<Interval> intervals;
    std::ifstream fcsv("runs/good_intervals_2026.csv");
    if (fcsv.is_open()) {
        std::string line;
        std::getline(fcsv, line); // header
        while (std::getline(fcsv, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            Interval iv;
            std::string chs_str;
            if (std::getline(iss, iv.run, ',') &&
                std::getline(iss, iv.start_local, ',') &&
                std::getline(iss, iv.end_local, ',') &&
                std::getline(iss, iv.status, ',') &&
                std::getline(iss, chs_str, ',')) {
                if (chs_str == "all") {
                    iv.all_channels = true;
                } else {
                    std::istringstream ch_iss(chs_str);
                    int c;
                    while (ch_iss >> c) iv.channels.push_back(c);
                }
                intervals.push_back(iv);
            }
        }
        fcsv.close();
    }

    // Read outputs/pe_calibration_v8.txt for G_v8_global
    std::map<int, double> G_v8_global_map;
    std::ifstream fcal("outputs/pe_calibration_v8.txt");
    if (fcal.is_open()) {
        std::string line;
        while (std::getline(fcal, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            int ch; std::string pl; double p0, G, s0, s1, chi2; int ndf, status;
            if (iss >> ch >> pl >> p0 >> G >> s0 >> s1 >> chi2 >> ndf >> status) {
                G_v8_global_map[ch] = G;
            }
        }
        fcal.close();
    }

    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // -04:00 offset

    struct RunFile {
        std::string run;
        std::string filename;
    };
    std::vector<RunFile> runs = {
        {"v8",  "lectura_vx1740_ext_nim_1M_v8.root"},
        {"v9",  "lectura_vx1740_ext_nim_1M_v9.root"},
        {"v10", "lectura_vx1740_ext_nim_1M_v10.root"}
    };

    std::map<std::string, HourData> hour_map;
    std::vector<std::string> hour_keys_order;

    for (const auto& rf : runs) {
        std::string resolved_path = hodo::resolveRawFile(rf.filename);
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
        std::cout << "Processing run " << rf.run << " (" << nEnt << " entries) from " << resolved_path << "...\n";

        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("time", 1);
        tree->SetBranchStatus("timestamp", 1);

        std::vector<float>* tw = nullptr;
        Long_t ts = 0;
        tree->SetBranchAddress("time", &tw);
        tree->SetBranchAddress("timestamp", &ts);

        std::vector<float>* cw[32] = {};
        for (int c = 0; c < 32; ++c) {
            if (c == 3) continue;
            tree->SetBranchStatus(Form("Ch%d", c), 1);
            tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
        }
        tree->SetCacheSize(256LL << 20);

        for (Long64_t iev = 0; iev < nEnt; ++iev) {
            tree->GetEntry(iev);
            if (iev % 50000 == 0) std::cout << "  " << iev << "/" << nEnt << "\r" << std::flush;
            if (!tw || tw->empty()) continue;

            double t_sec = boot_ref_epoch + ((double)ts) * 1e-6 + tz_offset;
            time_t t_local = (time_t)t_sec;
            struct tm tm_val;
            gmtime_r(&t_local, &tm_val);

            char date_buf[64], hour_buf[64], key_buf[128];
            std::snprintf(date_buf, sizeof(date_buf), "%04d-%02d-%02d",
                          tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday);
            std::snprintf(hour_buf, sizeof(hour_buf), "%02d:00", tm_val.tm_hour);
            std::snprintf(key_buf, sizeof(key_buf), "%s_%s_%s", rf.run.c_str(), date_buf, hour_buf);
            std::string key(key_buf);

            if (hour_map.find(key) == hour_map.end()) {
                HourData hd;
                hd.run = rf.run;
                hd.date = date_buf;
                hd.hour = hour_buf;
                hd.key = key;
                hd.first_ts = ts;
                for (int c = 0; c < 32; ++c) {
                    if (c == 3) continue;
                    hd.h_ph[c] = new TH1D(Form("h_ph_%s_ch%d", key.c_str(), c), "", 240, 0.0, 60.0);
                    hd.h_ph[c]->SetDirectory(nullptr);
                }
                hour_map[key] = std::move(hd);
                hour_keys_order.push_back(key);
            }

            HourData& hd = hour_map[key];
            hd.n_events++;
            hd.last_ts = ts;

            const auto& tv = *tw;
            const int n_samp = (int)tv.size();

            for (int c = 0; c < 32; ++c) {
                if (c == 3 || !cw[c] || (int)cw[c]->size() != n_samp) continue;
                const auto& wf = *cw[c];

                // Pulse finding in pretrigger range [100, 3480] ns
                double vmin = 1e30;
                int imin = -1;
                for (int i = 0; i < n_samp; ++i) {
                    if (tv[i] >= 100.0 && tv[i] <= 3480.0 && wf[i] < vmin) {
                        vmin = wf[i];
                        imin = i;
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
                    double ph_val = ped_loc - vmin;
                    hd.ph[c].push_back((float)ph_val);
                    hd.h_ph[c]->Fill(ph_val);
                }
            }
        }
        std::cout << "  " << nEnt << "/" << nEnt << " done.\n";
    }

    std::cout << "Total clock hours recorded across runs: " << hour_keys_order.size() << "\n";

    // Helper to check if channel c in hour (date, hour) is in HV_OFF
    auto is_channel_in_hv_off = [&](const std::string& run, const std::string& date, const std::string& hour, int c) -> bool {
        std::string dh = date + " " + hour; // e.g. "2026-05-30 11:00"
        for (const auto& iv : intervals) {
            if (iv.run == run && iv.status == "HV_OFF") {
                // Check if dh is between iv.start_local and iv.end_local
                // Since format is "YYYY-MM-DD HH:MM", lexical comparison works
                if (dh >= iv.start_local && dh <= iv.end_local) {
                    if (iv.all_channels) return true;
                    for (int ch_off : iv.channels) {
                        if (ch_off == c) return true;
                    }
                }
            }
        }
        return false;
    };

    // Fits and analysis
    std::vector<FitRecord> all_records;
    int total_attempted = 0;
    int total_ok = 0;

    for (const auto& k : hour_keys_order) {
        HourData& hd = hour_map[k];

        double dt_s = ((double)(hd.last_ts - hd.first_ts)) * 1e-6;
        double covered_time = 3600.0;
        if (dt_s < 3500.0 && dt_s > 0.0) covered_time = dt_s;
        double trig_rate = (covered_time > 0.0) ? ((double)hd.n_events / covered_time) : 0.0;
        double T_live = (double)hd.n_events * 3.38e-6;

        for (int c = 0; c < 32; ++c) {
            if (c == 3) continue;

            // Exclude hours in which channel is in an HV_OFF interval
            if (is_channel_in_hv_off(hd.run, hd.date, hd.hour, c)) {
                continue;
            }

            total_attempted++;

            FitRecord rec;
            rec.key = hd.key;
            rec.run = hd.run;
            rec.date = hd.date;
            rec.hour = hd.hour;
            rec.ch = c;
            rec.n_events = hd.n_events;
            rec.T_live = T_live;
            rec.trig_rate = trig_rate;

            TH1D* h = hd.h_ph[c];
            if (!h || h->GetEntries() < 20) {
                rec.status = -1;
                all_records.push_back(rec);
                continue;
            }

            double p0_start = 4.3;
            double G_start  = 4.0;
            double s0_start = 0.8;
            double s1_start = 0.3;

            double fit_lo = p0_start - G_start / 2.0;
            double fit_hi = p0_start + 5.5 * G_start;
            if (fit_lo < 0.25) fit_lo = 0.25;
            if (fit_hi > 59.0) fit_hi = 59.0;

            TF1 f_fit(Form("f_fit_%d_%s", c, k.c_str()), fit_6peaks, fit_lo, fit_hi, 10);
            f_fit.SetParameter(0, p0_start);
            f_fit.SetParameter(1, G_start);
            f_fit.SetParameter(2, s0_start);
            f_fit.SetParameter(3, s1_start);

            f_fit.SetParLimits(0, 0.5, 12.0);
            f_fit.SetParLimits(1, 0.5, 12.0);
            f_fit.SetParLimits(2, 0.1, 3.5);
            f_fit.SetParLimits(3, 0.0, 3.5);

            for (int n = 1; n <= 6; ++n) {
                double mu = p0_start + (n - 1) * G_start;
                double A_init = h->GetBinContent(h->FindBin(mu));
                if (A_init <= 0) A_init = 1.0;
                f_fit.SetParameter(3 + n, A_init);
                f_fit.SetParLimits(3 + n, 0.0, h->GetMaximum() * 2.0);
            }

            TFitResultPtr res = h->Fit(&f_fit, "QRN0S", "", fit_lo, fit_hi);
            int status = (int)res;

            double p0 = f_fit.GetParameter(0);
            double G  = f_fit.GetParameter(1);
            double s0 = f_fit.GetParameter(2);
            double s1 = f_fit.GetParameter(3);
            double chi2 = f_fit.GetChisquare();
            int ndf = f_fit.GetNDF();
            double chi2_ndf = (ndf > 0) ? (chi2 / ndf) : 999.0;

            rec.status = status;
            rec.p0 = p0;
            rec.G = G;
            rec.chi2ndf = chi2_ndf;

            if (status == 0 && chi2_ndf < 3.0 && G > 0.5 && G < 10.0) {
                rec.is_ok = true;
                total_ok++;
                double A2 = f_fit.GetParameter(5);
                double sigma2 = std::sqrt(s0 * s0 + s1 * s1);
                rec.N1pe = (A2 * std::sqrt(2.0 * M_PI) * sigma2) / 0.25;
            } else {
                rec.is_ok = false;
                rec.N1pe = 0.0;
            }

            // Pulse rate (> p0 + 0.5*G) and Iproxy
            double p0_use = rec.is_ok ? p0 : p0_start;
            double G_use  = rec.is_ok ? G  : G_start;

            long long n_pulses = 0;
            double sum_amp_p0 = 0.0;
            for (float val : hd.ph[c]) {
                if (val > p0_use + 0.5 * G_use) {
                    n_pulses++;
                    sum_amp_p0 += (val - p0_use);
                }
            }

            rec.rate_05 = (T_live > 0.0) ? ((double)n_pulses / T_live) : 0.0;
            rec.Iproxy  = (T_live > 0.0) ? (sum_amp_p0 / T_live) : 0.0;

            all_records.push_back(rec);
        }
    }

    // Save CSV
    std::ofstream out_csv("outputs/hv_diag3_gain_vs_time.csv");
    out_csv << "run,date,hour,ch,n_events,T_live,G,p0,chi2ndf,status,N1pe,rate_05,Iproxy,trig_rate\n";
    for (const auto& r : all_records) {
        out_csv << r.run << ","
                << r.date << ","
                << r.hour << ","
                << r.ch << ","
                << r.n_events << ","
                << std::fixed << std::setprecision(6) << r.T_live << ","
                << std::setprecision(4) << r.G << ","
                << r.p0 << ","
                << r.chi2ndf << ","
                << r.status << ","
                << std::setprecision(2) << r.N1pe << ","
                << std::setprecision(3) << r.rate_05 << ","
                << std::setprecision(3) << r.Iproxy << ","
                << std::setprecision(4) << r.trig_rate << "\n";
    }
    out_csv.close();
    std::cout << "Saved: outputs/hv_diag3_gain_vs_time.csv\n";

    // -------------------------------------------------------------------------
    // Report (a)
    // -------------------------------------------------------------------------
    std::cout << "\n========================================================================================\n";
    std::cout << "PART 1 (a) Fits Status OK and chi2/ndf < 3\n";
    std::cout << "========================================================================================\n";
    std::cout << "Total attempted fits: " << total_attempted << "\n";
    std::cout << "Fits with status OK (0) and chi2/ndf < 3: " << total_ok << "\n";
    std::cout << "Percentage: " << std::fixed << std::setprecision(2)
              << (total_attempted > 0 ? (100.0 * total_ok / total_attempted) : 0.0) << " %\n";
    std::cout << "========================================================================================\n\n";

    // -------------------------------------------------------------------------
    // Report (b)
    // Per channel: G from the hourly fits averaged over all OK hours of v8, next to G from outputs/pe_calibration_v8.txt.
    // Columns EXACTLY: ch  G_hourly_mean  G_hourly_rms  G_v8_global
    // -------------------------------------------------------------------------
    std::cout << "========================================================================================\n";
    std::cout << "PART 1 (b) Per-Channel Hourly G on Run v8 vs Global v8 G\n";
    std::cout << "========================================================================================\n";
    std::cout << std::left << std::setw(6) << "ch"
              << std::right << std::setw(16) << "G_hourly_mean"
              << std::setw(16) << "G_hourly_rms"
              << std::setw(16) << "G_v8_global" << "\n";
    std::cout << std::string(54, '-') << "\n";

    for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        std::vector<double> g_vals;
        for (const auto& r : all_records) {
            if (r.run == "v8" && r.ch == c && r.is_ok) {
                g_vals.push_back(r.G);
            }
        }
        double mean = 0.0, rms = 0.0;
        if (!g_vals.empty()) {
            double sum = 0.0;
            for (double g : g_vals) sum += g;
            mean = sum / g_vals.size();
            double s2 = 0.0;
            for (double g : g_vals) s2 += (g - mean) * (g - mean);
            rms = (g_vals.size() > 1) ? std::sqrt(s2 / (g_vals.size() - 1)) : 0.0;
        }
        double g_global = G_v8_global_map.count(c) ? G_v8_global_map[c] : 0.0;

        std::cout << std::left << std::setw(6) << c
                  << std::right << std::fixed << std::setprecision(4)
                  << std::setw(16) << mean
                  << std::setw(16) << rms
                  << std::setw(16) << g_global << "\n";
    }
    std::cout << "========================================================================================\n\n";

    // -------------------------------------------------------------------------
    // Report (c)
    // One table, one row per hour, columns EXACTLY:
    // run  date  hour  n_ch_ok  median_G_X  median_G_Y  median_rate_05_X  median_Iproxy_X  trig_rate
    // -------------------------------------------------------------------------
    std::cout << "===============================================================================================================\n";
    std::cout << "PART 1 (c) Hourly Medians Across Runs v8, v9, v10\n";
    std::cout << "===============================================================================================================\n";
    std::cout << std::left << std::setw(6) << "run"
              << std::setw(12) << "date"
              << std::setw(8) << "hour"
              << std::right << std::setw(9) << "n_ch_ok"
              << std::setw(13) << "median_G_X"
              << std::setw(13) << "median_G_Y"
              << std::setw(18) << "median_rate_05_X"
              << std::setw(18) << "median_Iproxy_X"
              << std::setw(12) << "trig_rate" << "\n";
    std::cout << std::string(109, '-') << "\n";

    // We will also collect points for PNG
    std::vector<double> plot_times;
    std::vector<double> plot_g_x, plot_g_y;
    std::vector<double> plot_trig_rate;
    std::vector<double> plot_rate05_x;

    for (const auto& k : hour_keys_order) {
        HourData& hd = hour_map[k];

        std::vector<double> gx, gy, rx, ix;
        int n_ok = 0;
        double cur_trig_rate = 0.0;

        for (const auto& r : all_records) {
            if (r.key == k) {
                cur_trig_rate = r.trig_rate;
                if (r.is_ok) {
                    n_ok++;
                    if (r.ch >= 0 && r.ch <= 15) {
                        gx.push_back(r.G);
                        rx.push_back(r.rate_05);
                        ix.push_back(r.Iproxy);
                    } else if (r.ch >= 16 && r.ch <= 31) {
                        gy.push_back(r.G);
                    }
                }
            }
        }

        double med_gx = compute_median(gx);
        double med_gy = compute_median(gy);
        double med_rx = compute_median(rx);
        double med_ix = compute_median(ix);

        std::cout << std::left << std::setw(6) << hd.run
                  << std::setw(12) << hd.date
                  << std::setw(8) << hd.hour
                  << std::right << std::setw(9) << n_ok
                  << std::fixed << std::setprecision(4)
                  << std::setw(13) << med_gx
                  << std::setw(13) << med_gy
                  << std::setprecision(2)
                  << std::setw(18) << med_rx
                  << std::setw(18) << med_ix
                  << std::setprecision(4)
                  << std::setw(12) << cur_trig_rate << "\n";

        // Convert date + hour to calendar timestamp
        struct tm tm_h{};
        int yr, mo, dy, hr;
        if (std::sscanf(hd.date.c_str(), "%d-%d-%d", &yr, &mo, &dy) == 3 &&
            std::sscanf(hd.hour.c_str(), "%d", &hr) == 1) {
            tm_h.tm_year = yr - 1900;
            tm_h.tm_mon  = mo - 1;
            tm_h.tm_mday = dy;
            tm_h.tm_hour = hr;
            tm_h.tm_min  = 0;
            tm_h.tm_sec  = 0;
            time_t t_plot = timegm(&tm_h); // UTC epoch of this local time string
            plot_times.push_back((double)t_plot);
            plot_g_x.push_back(med_gx);
            plot_g_y.push_back(med_gy);
            plot_trig_rate.push_back(cur_trig_rate);
            plot_rate05_x.push_back(med_rx);
        }
    }
    std::cout << "===============================================================================================================\n\n";

    // -------------------------------------------------------------------------
    // (d) Save PNG outputs/hv_diag3_gain_vs_time.png
    // Three pads sharing the time axis
    // -------------------------------------------------------------------------
    if (!plot_times.empty()) {
        TCanvas c_plot("c_gain_time", "Gain and Rates vs Time", 1200, 900);
        c_plot.Divide(1, 3, 0.001, 0.001);

        int np = (int)plot_times.size();
        auto gr_gx = std::make_unique<TGraph>(np, plot_times.data(), plot_g_x.data());
        auto gr_gy = std::make_unique<TGraph>(np, plot_times.data(), plot_g_y.data());
        auto gr_trig = std::make_unique<TGraph>(np, plot_times.data(), plot_trig_rate.data());
        auto gr_rx = std::make_unique<TGraph>(np, plot_times.data(), plot_rate05_x.data());

        // Pad 1: Median G
        c_plot.cd(1);
        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.10);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.04);
        gPad->SetGrid();

        gr_gx->SetTitle("hod2018 Hourly Diagnostics (Runs v8 - v10);;Median G [mV/p.e.]");
        gr_gx->SetMarkerStyle(20);
        gr_gx->SetMarkerSize(0.8);
        gr_gx->SetMarkerColor(kBlue + 1);
        gr_gx->SetLineColor(kBlue + 1);
        gr_gx->GetYaxis()->SetRangeUser(3.0, 6.0);
        gr_gx->GetXaxis()->SetTimeDisplay(1);
        gr_gx->GetXaxis()->SetTimeFormat("%m/%d %H:%M");
        gr_gx->GetXaxis()->SetTimeOffset(0);
        gr_gx->GetXaxis()->SetLabelSize(0.0); // hide x label to share axis
        gr_gx->GetYaxis()->SetTitleSize(0.06);
        gr_gx->GetYaxis()->SetTitleOffset(0.7);
        gr_gx->Draw("ALP");

        gr_gy->SetMarkerStyle(21);
        gr_gy->SetMarkerSize(0.8);
        gr_gy->SetMarkerColor(kRed + 1);
        gr_gy->SetLineColor(kRed + 1);
        gr_gy->Draw("LP SAME");

        TLegend leg1(0.75, 0.65, 0.94, 0.88);
        leg1.AddEntry(gr_gx.get(), "18X Median G", "lp");
        leg1.AddEntry(gr_gy.get(), "18Y Median G", "lp");
        leg1.Draw();

        // Pad 2: Trigger Rate
        c_plot.cd(2);
        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.02);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.04);
        gPad->SetGrid();

        gr_trig->SetTitle(";;Trigger rate [Hz]");
        gr_trig->SetMarkerStyle(20);
        gr_trig->SetMarkerSize(0.7);
        gr_trig->SetMarkerColor(kBlack);
        gr_trig->SetLineColor(kBlack);
        gr_trig->GetYaxis()->SetRangeUser(0.0, 0.8);
        gr_trig->GetXaxis()->SetTimeDisplay(1);
        gr_trig->GetXaxis()->SetTimeFormat("%m/%d %H:%M");
        gr_trig->GetXaxis()->SetTimeOffset(0);
        gr_trig->GetXaxis()->SetLabelSize(0.0);
        gr_trig->GetYaxis()->SetTitleSize(0.06);
        gr_trig->GetYaxis()->SetTitleOffset(0.7);
        gr_trig->Draw("ALP");

        // Pad 3: Median rate_05
        c_plot.cd(3);
        gPad->SetBottomMargin(0.18);
        gPad->SetTopMargin(0.02);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.04);
        gPad->SetGrid();

        gr_rx->SetTitle(";;Median Rate 18X (>0.5 p.e.) [Hz]");
        gr_rx->SetMarkerStyle(20);
        gr_rx->SetMarkerSize(0.8);
        gr_rx->SetMarkerColor(kGreen + 2);
        gr_rx->SetLineColor(kGreen + 2);
        gr_rx->GetXaxis()->SetTimeDisplay(1);
        gr_rx->GetXaxis()->SetTimeFormat("%m/%d %H:%M");
        gr_rx->GetXaxis()->SetTimeOffset(0);
        gr_rx->GetXaxis()->SetLabelSize(0.05);
        gr_rx->GetXaxis()->SetTitleSize(0.06);
        gr_rx->GetXaxis()->SetTitle("Local Time (UTC-4)");
        gr_rx->GetYaxis()->SetTitleSize(0.06);
        gr_rx->GetYaxis()->SetTitleOffset(0.7);
        gr_rx->Draw("ALP");

        c_plot.SaveAs("outputs/hv_diag3_gain_vs_time.png");
        std::cout << "Saved figure: outputs/hv_diag3_gain_vs_time.png\n";
    }
}
