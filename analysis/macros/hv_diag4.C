// analysis/macros/hv_diag4.C
// Task HV-DIAG-4: occupancy-corrected rates, morning light test, ch9 gain

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TFitResult.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <ctime>

#include "../src/hodo_paths.h"

// 6-peak model as in pe_scale_step2.C
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

// -----------------------------------------------------------------------------
// PART A: Occupancy-corrected rates from hv_diag3_pretrip CSV files
// -----------------------------------------------------------------------------
void run_part_a() {
    std::cout << "========================================================================================\n";
    std::cout << "TASK HV-DIAG-4 PART A: Occupancy-Corrected Pre-Trip Rates\n";
    std::cout << "========================================================================================\n";
    std::cout << "Notice: n_ev per bin is computed from live time as n_ev = round(live_time_s / 3.38e-6 s).\n\n";

    std::vector<std::pair<std::string, std::string>> trips = {
        {"T1", "outputs/hv_diag3_pretrip_T1.csv"},
        {"T2", "outputs/hv_diag3_pretrip_T2.csv"},
        {"T3", "outputs/hv_diag3_pretrip_T3.csv"}
    };

    struct BinDataA {
        double t_start = 0.0;
        double t_end = 0.0;
        int nev = 0;
        std::map<int, int> counts;
        std::map<int, double> rates_khz;
        std::map<int, bool> sat;
        double sum56_63 = 0.0;
        bool sum_sat = false;
    };

    for (const auto& item : trips) {
        std::string trip_name = item.first;
        std::string fn = item.second;

        std::map<int, BinDataA> bins;
        std::ifstream fin(fn);
        if (!fin.is_open()) {
            std::cerr << "Error opening " << fn << "\n";
            continue;
        }

        std::string line;
        std::getline(fin, line); // header
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            int b, ch, pc;
            double t_s, t_e, amp, lt;
            char comma;
            if (iss >> b >> comma >> t_s >> comma >> t_e >> comma >> ch >> comma >> pc >> comma >> amp >> comma >> lt) {
                if (bins.find(b) == bins.end()) {
                    bins[b].t_start = t_s;
                    bins[b].t_end   = t_e;
                    bins[b].nev     = (int)std::round(lt / 3.38e-6);
                }
                bins[b].counts[ch] = pc;
            }
        }
        fin.close();

        // Compute rate_true in kHz for each bin
        for (auto& pair : bins) {
            int b = pair.first;
            auto& bd = pair.second;
            int nev = bd.nev;

            std::vector<int> target_chs = {45, 56, 57, 58, 59, 60, 61, 62, 63};
            double s56_63 = 0.0;
            bool s_sat = false;

            for (int ch : target_chs) {
                int k = bd.counts.count(ch) ? bd.counts[ch] : 0;
                if (nev == 0) {
                    bd.rates_khz[ch] = 0.0;
                    bd.sat[ch] = false;
                } else if (k >= nev) {
                    bd.rates_khz[ch] = 0.0;
                    bd.sat[ch] = true;
                    if (ch >= 56 && ch <= 63) s_sat = true;
                } else {
                    double occ = (double)k / (double)nev;
                    double r = -std::log(1.0 - occ) / 3.38e-3; // in kHz
                    bd.rates_khz[ch] = r;
                    bd.sat[ch] = false;
                    if (ch >= 56 && ch <= 63) s56_63 += r;
                }
            }
            bd.sum56_63 = s56_63;
            bd.sum_sat = s_sat;
        }

        // Baseline mean in [onset-120, onset-60] min (bins 0 to 29)
        double base_sum = 0.0;
        int n_base = 0;
        for (int b = 0; b <= 29; ++b) {
            if (bins.count(b) && !bins[b].sum_sat) {
                base_sum += bins[b].sum56_63;
                n_base++;
            }
        }
        double mean_base = (n_base > 0) ? (base_sum / n_base) : 0.0;
        double thresh = 2.0 * mean_base;

        std::cout << "------------------------------------------------------------------------------------------------------------------------\n";
        std::cout << "Trip " << trip_name << " (Baseline mean sum56_63 in [-120, -60] min: " << std::fixed << std::setprecision(2)
                  << mean_base << " kHz, 2x mean: " << thresh << " kHz)\n";
        std::cout << "------------------------------------------------------------------------------------------------------------------------\n";
        std::cout << std::right << std::setw(6) << "t_min"
                  << std::setw(9) << "rate56"
                  << std::setw(9) << "rate57"
                  << std::setw(9) << "rate58"
                  << std::setw(9) << "rate59"
                  << std::setw(9) << "rate60"
                  << std::setw(9) << "rate61"
                  << std::setw(9) << "rate62"
                  << std::setw(9) << "rate63"
                  << std::setw(9) << "rate45"
                  << std::setw(11) << "sum56_63" << "\n";
        std::cout << std::string(98, '-') << "\n";

        auto fmt_rate = [](double r, bool sat) -> std::string {
            if (sat) return "      SAT";
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << std::setw(9) << r;
            return oss.str();
        };

        // Bins from onset-60 to onset+4 min:
        // bin 30 (t_start = -60.0) through bin 62 (t_start = 4.0)
        for (int b = 30; b <= 62; ++b) {
            if (!bins.count(b)) continue;
            const auto& bd = bins[b];
            std::cout << std::fixed << std::setprecision(1) << std::setw(6) << bd.t_start
                      << fmt_rate(bd.rates_khz.at(56), bd.sat.at(56))
                      << fmt_rate(bd.rates_khz.at(57), bd.sat.at(57))
                      << fmt_rate(bd.rates_khz.at(58), bd.sat.at(58))
                      << fmt_rate(bd.rates_khz.at(59), bd.sat.at(59))
                      << fmt_rate(bd.rates_khz.at(60), bd.sat.at(60))
                      << fmt_rate(bd.rates_khz.at(61), bd.sat.at(61))
                      << fmt_rate(bd.rates_khz.at(62), bd.sat.at(62))
                      << fmt_rate(bd.rates_khz.at(63), bd.sat.at(63))
                      << fmt_rate(bd.rates_khz.at(45), bd.sat.at(45));
            if (bd.sum_sat) {
                std::cout << "        SAT\n";
            } else {
                std::cout << std::fixed << std::setprecision(2) << std::setw(11) << bd.sum56_63 << "\n";
            }
        }
        std::cout << "------------------------------------------------------------------------------------------------------------------------\n";

        // First bin, counting back from onset (t=0, bin 59 downwards), after which sum56_63 stays above 2x mean
        // Check consecutive bins above threshold ending at or near onset
        int first_bin_above = -1;
        for (int b = 59; b >= 0; --b) {
            if (!bins.count(b)) break;
            if (bins[b].sum56_63 > thresh) {
                first_bin_above = b;
            } else {
                // dropped below
                break;
            }
        }

        std::cout << "Onset-rise evaluation for Trip " << trip_name << ":\n";
        if (first_bin_above != -1) {
            double min_before_onset = -bins[first_bin_above].t_start;
            std::cout << "  Counting back from onset (bin 59): sum56_63 is above 2x mean continuously back to bin "
                      << first_bin_above << " [" << bins[first_bin_above].t_start << ", " << bins[first_bin_above].t_end
                      << "] min (" << min_before_onset << " minutes before onset).\n";
        } else {
            // Check if bin 59 was just slightly below and bin 58 was above
            std::cout << "  Note: At bin 59 [-2.0, 0.0] min, sum56_63 = " << bins[59].sum56_63
                      << " kHz (threshold = " << thresh << " kHz).\n";
            int b_peak = -1;
            for (int b = 58; b >= 30; --b) {
                if (bins[b].sum56_63 > thresh) {
                    b_peak = b;
                    break;
                }
            }
            if (b_peak != -1) {
                std::cout << "  The last pre-trip bin exceeding 2x mean is bin " << b_peak << " ["
                          << bins[b_peak].t_start << ", " << bins[b_peak].t_end << "] min ("
                          << -bins[b_peak].t_start << " minutes before onset).\n";
            }
        }
        std::cout << "\n";
    }
}

// -----------------------------------------------------------------------------
// PART B: Morning light test (raw data)
// -----------------------------------------------------------------------------
struct ChannelCalib {
    double p0 = 4.3;
    double G  = 4.0;
};

void run_part_b() {
    std::cout << "========================================================================================\n";
    std::cout << "TASK HV-DIAG-4 PART B: Morning Light Test (07:00 - 12:00 Local, 10-Min Bins)\n";
    std::cout << "========================================================================================\n\n";

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

    // 2. Load hourly calibration from outputs/hv_diag3_gain_vs_time.csv
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

    struct MorningDay {
        std::string run;
        std::string date; // "YYYY-MM-DD"
        std::string raw_file;
        Long64_t est_start_ev;
        Long64_t est_end_ev;
    };

    std::vector<MorningDay> mornings = {
        {"v8",  "2026-05-29", "lectura_vx1740_ext_nim_1M_v8.root",  32000, 48000},
        {"v8",  "2026-05-31", "lectura_vx1740_ext_nim_1M_v8.root",  134000, 143814},
        {"v10", "2026-06-02", "lectura_vx1740_ext_nim_1M_v10.root", 77000, 91000},
        {"v10", "2026-06-04", "lectura_vx1740_ext_nim_1M_v10.root", 173000, 190000},
        {"v10", "2026-06-05", "lectura_vx1740_ext_nim_1M_v10.root", 223000, 238000}
    };

    // Helper to check if channel is inside an HV_OFF interval of good_intervals_2026.csv
    auto is_channel_hv_off = [](const std::string& run, const std::string& date, int hour, int min, int ch) -> bool {
        // v8, 2026-05-30 11:00 to 2026-05-31 00:59: 6 20 46 50 56 59 60 62
        if (run == "v8") {
            // On May 29: none
            // On May 31 07:00-12:00: none (HV_OFF ended at 00:59)
            return false;
        }
        if (run == "v10") {
            // v10, 2026-06-01 09:00 to 2026-06-02 12:59:
            // 0 6 20 22 31 43 46 47 50 54 56 58 59 60 62 63
            if (date == "2026-06-02") {
                // Entire morning 07:00-12:00 is inside this interval
                static const std::vector<int> off_june2 = {0, 6, 20, 22, 31, 43, 46, 47, 50, 54, 56, 58, 59, 60, 62, 63};
                return std::find(off_june2.begin(), off_june2.end(), ch) != off_june2.end();
            }
            // v10, 2026-06-03 09:00 to 2026-06-05 15:59:
            // 0 1 6 16 17 20 22 31 34 43 46 47 50 53 54 55 56 57 58 59 60 62 63
            if (date == "2026-06-04" || date == "2026-06-05") {
                // Entire morning 07:00-12:00 is inside this interval
                static const std::vector<int> off_june3_5 = {0, 1, 6, 16, 17, 20, 22, 31, 34, 43, 46, 47, 50, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63};
                return std::find(off_june3_5.begin(), off_june3_5.end(), ch) != off_june3_5.end();
            }
        }
        return false;
    };

    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // UTC-4

    // Channels needed for Part B:
    // 45, 52, 56, 57, 58, 59, 60, 61, 62, 63, and hod2018 ch 0-15
    std::vector<int> needed_chs;
    for (int c = 0; c < 16; ++c) needed_chs.push_back(c);
    std::vector<int> extra_chs = {45, 52, 56, 57, 58, 59, 60, 61, 62, 63};
    for (int c : extra_chs) needed_chs.push_back(c);

    for (const auto& md : mornings) {
        std::string resolved_path = hodo::resolveRawFile(md.raw_file);
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
        for (int c : needed_chs) {
            tree->SetBranchStatus(Form("Ch%d", c), 1);
            tree->SetBranchAddress(Form("Ch%d", c), &cw[c]);
        }
        tree->SetCacheSize(256LL << 20);

        // 30 bins of 10 minutes: 07:00 to 12:00
        constexpr int kNBinsB = 30;
        struct BinDataB {
            Long64_t n_ev = 0;
            int pulse_count[64] = {};
            bool ch_off[64] = {};
        };
        std::vector<BinDataB> bins(kNBinsB);

        // Pre-mark OFF channels for each bin
        for (int b = 0; b < kNBinsB; ++b) {
            int h = 7 + b / 6;
            int m = (b % 6) * 10;
            for (int c : needed_chs) {
                bins[b].ch_off[c] = is_channel_hv_off(md.run, md.date, h, m, c);
            }
        }

        // Scan entries covering 07:00 - 12:00
        Long64_t start_ev = std::max((Long64_t)0, md.est_start_ev);
        Long64_t end_ev   = std::min(nEnt - 1, md.est_end_ev);

        for (Long64_t iev = start_ev; iev <= end_ev; ++iev) {
            tree->GetEntry(iev);
            if (!tw || tw->empty()) continue;

            double t_sec = boot_ref_epoch + ((double)ts) * 1e-6 + tz_offset;
            time_t t_local = (time_t)t_sec;
            struct tm tm_val;
            gmtime_r(&t_local, &tm_val);

            char date_buf[32];
            std::snprintf(date_buf, sizeof(date_buf), "%04d-%02d-%02d",
                          tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday);
            if (std::string(date_buf) != md.date) continue;

            int hour = tm_val.tm_hour;
            int min  = tm_val.tm_min;
            if (hour < 7 || hour >= 12) continue;

            int b = (hour - 7) * 6 + (min / 10);
            if (b < 0 || b >= kNBinsB) continue;

            bins[b].n_ev++;

            // Hourly lookup key
            char dh_buf[64];
            std::snprintf(dh_buf, sizeof(dh_buf), "%s_%s_%02d:00", md.run.c_str(), md.date.c_str(), hour);
            std::string hour_key(dh_buf);

            const auto& tv = *tw;
            const int n_samp = (int)tv.size();

            for (int c : needed_chs) {
                if (c == 3 || bins[b].ch_off[c]) continue;
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
                        std::string ch_key = hour_key + "_" + std::to_string(c);
                        double G  = hourly_G_map.count(ch_key) ? hourly_G_map[ch_key] :
                                    (calib_global.count(c) ? calib_global[c].G : 4.0);
                        double p0 = hourly_p0_map.count(ch_key) ? hourly_p0_map[ch_key] :
                                    (calib_global.count(c) ? calib_global[c].p0 : 4.3);
                        thr = p0 + 0.5 * G;
                    }
                    if (ph > thr) {
                        bins[b].pulse_count[c]++;
                    }
                }
            }
        }

        // Print Day Table
        std::cout << "==================================================================================================================================================\n";
        std::cout << "PART B Table: Run " << md.run << ", Date " << md.date << " (07:00 - 12:00 Local)\n";
        std::cout << "==================================================================================================================================================\n";
        std::cout << std::left << std::setw(7) << "time"
                  << std::right << std::setw(7) << "n_ev"
                  << std::setw(9) << "rate45"
                  << std::setw(9) << "rate52"
                  << std::setw(9) << "rate56"
                  << std::setw(9) << "rate57"
                  << std::setw(9) << "rate58"
                  << std::setw(9) << "rate59"
                  << std::setw(9) << "rate60"
                  << std::setw(9) << "rate61"
                  << std::setw(9) << "rate62"
                  << std::setw(9) << "rate63"
                  << std::setw(14) << "sum18X_rate"
                  << std::setw(12) << "trig_rate" << "\n";
        std::cout << std::string(118, '-') << "\n";

        for (int b = 0; b < kNBinsB; ++b) {
            int h = 7 + b / 6;
            int m = (b % 6) * 10;
            char time_buf[16];
            std::snprintf(time_buf, sizeof(time_buf), "%02d:%02d", h, m);

            Long64_t nev = bins[b].n_ev;
            double trig_r = (double)nev / 600.0;

            auto get_rate_str = [&](int c) -> std::string {
                if (bins[b].ch_off[c]) return "      OFF";
                if (nev == 0) return "     0.00";
                int k = bins[b].pulse_count[c];
                if (k >= nev) return "      SAT";
                double occ = (double)k / (double)nev;
                double r_khz = -std::log(1.0 - occ) / 3.38e-3;
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(2) << std::setw(9) << r_khz;
                return oss.str();
            };

            // sum18X_rate
            double sum18X = 0.0;
            bool sum18X_sat = false;
            for (int c = 0; c < 16; ++c) {
                if (c == 3 || bins[b].ch_off[c]) continue;
                if (nev == 0) continue;
                int k = bins[b].pulse_count[c];
                if (k >= nev) {
                    sum18X_sat = true;
                } else {
                    double occ = (double)k / (double)nev;
                    sum18X += (-std::log(1.0 - occ) / 3.38e-3);
                }
            }

            std::cout << std::left << std::setw(7) << time_buf
                      << std::right << std::setw(7) << nev
                      << get_rate_str(45)
                      << get_rate_str(52)
                      << get_rate_str(56)
                      << get_rate_str(57)
                      << get_rate_str(58)
                      << get_rate_str(59)
                      << get_rate_str(60)
                      << get_rate_str(61)
                      << get_rate_str(62)
                      << get_rate_str(63);
            if (sum18X_sat) {
                std::cout << "           SAT";
            } else {
                std::cout << std::fixed << std::setprecision(2) << std::setw(14) << sum18X;
            }
            std::cout << std::fixed << std::setprecision(4) << std::setw(12) << trig_r << "\n";
        }
        std::cout << "==================================================================================================================================================\n\n";
    }
}

// -----------------------------------------------------------------------------
// PART C: Channel 9 gain refit
// -----------------------------------------------------------------------------
void run_part_c() {
    std::cout << "========================================================================================\n";
    std::cout << "TASK HV-DIAG-4 PART C: Channel 9 Control-Window Gain Refit (Run v8)\n";
    std::cout << "========================================================================================\n";

    std::string raw_path = hodo::resolveRawFile("lectura_vx1740_ext_nim_1M_v8.root");
    std::unique_ptr<TFile> fin(TFile::Open(raw_path.c_str(), "READ"));
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening " << raw_path << "\n";
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(fin->Get("data_dig"));
    if (!tree) {
        std::cerr << "Error finding data_dig in " << raw_path << "\n";
        return;
    }

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("time", 1);
    tree->SetBranchStatus("Ch9", 1);

    std::vector<float>* tw = nullptr;
    std::vector<float>* cw = nullptr;
    tree->SetBranchAddress("time", &tw);
    tree->SetBranchAddress("Ch9", &cw);

    TH1D h_ph("h_ph9", "Ch 9 Control-Window Spectrum (v8);ph [mV];Events / 0.25 mV", 240, 0.0, 60.0);
    h_ph.SetDirectory(nullptr);

    Long64_t nEnt = tree->GetEntries();
    std::cout << "Reading " << nEnt << " events for Channel 9...\n";
    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        tree->GetEntry(iev);
        if (!tw || !cw) continue;
        const auto& tv = *tw;
        const auto& wf = *cw;
        int n_samp = (int)tv.size();

        double vmin = 1e30; int imin = -1;
        for (int i = 0; i < n_samp; ++i) {
            if (tv[i] >= 6128.0 && tv[i] <= 6272.0 && wf[i] < vmin) {
                vmin = wf[i]; imin = i;
            }
        }
        if (imin >= 0) {
            double ped = 0.0;
            if (imin >= 12 && imin - 4 < n_samp) {
                double s = 0.0;
                for (int j = -12; j <= -4; ++j) s += wf[imin + j];
                ped = s / 9.0;
            } else {
                double s = 0.0; int np = 0;
                for (int i = 0; i < n_samp; ++i) if (tv[i] < 2500.0) { s += wf[i]; np++; }
                ped = (np > 0) ? (s / np) : 0.0;
            }
            h_ph.Fill(ped - vmin);
        }
    }
    std::cout << "Filled histogram: " << h_ph.GetEntries() << " entries.\n\n";

    struct FitResultC {
        double G_start;
        double p0;
        double G;
        double s0;
        double s1;
        double chi2;
        int ndf;
        double chi2_ndf;
        int status;
        std::vector<double> peak_pos;
    };

    auto do_fit = [&](double G_start, const char* fname, int color) -> TF1* {
        double p0_start = 4.3;
        double s0_start = 0.8;
        double s1_start = 0.3;

        double fit_lo = p0_start - G_start / 2.0;
        double fit_hi = p0_start + 5.5 * G_start;
        if (fit_lo < 0.25) fit_lo = 0.25;
        if (fit_hi > 59.0) fit_hi = 59.0;

        TF1* f_fit = new TF1(fname, fit_6peaks, fit_lo, fit_hi, 10);
        f_fit->SetParameter(0, p0_start);
        f_fit->SetParameter(1, G_start);
        f_fit->SetParameter(2, s0_start);
        f_fit->SetParameter(3, s1_start);

        f_fit->SetParLimits(0, 0.5, 12.0);
        f_fit->SetParLimits(1, 0.5, 12.0);
        f_fit->SetParLimits(2, 0.1, 3.5);
        f_fit->SetParLimits(3, 0.0, 3.5);

        for (int n = 1; n <= 6; ++n) {
            double mu = p0_start + (n - 1) * G_start;
            double A_init = h_ph.GetBinContent(h_ph.FindBin(mu));
            if (A_init <= 0) A_init = 1.0;
            f_fit->SetParameter(3 + n, A_init);
            f_fit->SetParLimits(3 + n, 0.0, h_ph.GetMaximum() * 2.0);
        }
        f_fit->SetLineColor(color);
        f_fit->SetLineWidth(2);

        TFitResultPtr res = h_ph.Fit(f_fit, "QRN0S", "", fit_lo, fit_hi);
        int status = (int)res;

        double p0 = f_fit->GetParameter(0);
        double G  = f_fit->GetParameter(1);
        double s0 = f_fit->GetParameter(2);
        double s1 = f_fit->GetParameter(3);
        double chi2 = f_fit->GetChisquare();
        int ndf = f_fit->GetNDF();
        double chi2_ndf = (ndf > 0) ? (chi2 / ndf) : 0.0;

        std::cout << "Fit with G_start = " << std::fixed << std::setprecision(1) << G_start << " mV:\n";
        std::cout << "  Status:    " << status << "\n";
        std::cout << "  p0:        " << std::setprecision(5) << p0 << " mV\n";
        std::cout << "  G:         " << std::setprecision(5) << G << " mV/pe\n";
        std::cout << "  s0:        " << std::setprecision(5) << s0 << " mV\n";
        std::cout << "  s1:        " << std::setprecision(5) << s1 << " mV\n";
        std::cout << "  chi2/ndf:  " << std::setprecision(4) << chi2_ndf << "  (" << std::setprecision(2) << chi2 << " / " << ndf << ")\n";
        std::cout << "  Peak positions (mV):\n";
        for (int n = 1; n <= 6; ++n) {
            double pos = p0 + (n - 1) * G;
            std::cout << "    Peak " << n << ": " << std::setprecision(5) << pos << " mV\n";
        }
        std::cout << "\n";

        return f_fit;
    };

    TF1* fit_47 = do_fit(4.7, "fit_ch9_G47", kRed + 1);
    TF1* fit_30 = do_fit(3.0, "fit_ch9_G30", kBlue + 1);

    // Save PNG
    TCanvas c("c_ch9", "Ch 9 Gain Fits", 900, 600);
    gPad->SetGrid();
    gPad->SetLogy();
    h_ph.SetLineColor(kBlack);
    h_ph.GetXaxis()->SetRangeUser(0.0, 45.0);
    h_ph.Draw("HIST");
    fit_47->Draw("SAME");
    fit_30->Draw("SAME");

    TLegend leg(0.55, 0.68, 0.88, 0.88);
    leg.SetBorderSize(1);
    leg.SetFillColor(kWhite);
    leg.AddEntry(&h_ph, "Ch 9 Control ph (v8)", "l");
    leg.AddEntry(fit_47, Form("Fit G_start=4.7 mV (G=%.2f mV/pe, #chi^{2}/ndf=%.1f)", fit_47->GetParameter(1), fit_47->GetChisquare()/fit_47->GetNDF()), "l");
    leg.AddEntry(fit_30, Form("Fit G_start=3.0 mV (G=%.2f mV/pe, #chi^{2}/ndf=%.1f)", fit_30->GetParameter(1), fit_30->GetChisquare()/fit_30->GetNDF()), "l");
    leg.Draw();

    c.SaveAs("outputs/hv_diag4_ch9_fits.png");
    std::cout << "Saved figure: outputs/hv_diag4_ch9_fits.png\n\n";
}

void hv_diag4() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kError;

    run_part_a();
    run_part_b();
    run_part_c();
}
