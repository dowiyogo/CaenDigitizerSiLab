// analysis/macros/hv_diag5_episodes.C
// Task HV-DIAG-5: complete inventory of dead-channel episodes, v7-v10

#include <TFile.h>
#include <TTree.h>
#include <TError.h>
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <ctime>

#include "../src/hodo_paths.h"

struct ChannelCalib {
    double p0 = 4.3;
    double G  = 4.0;
};

struct BinInfo {
    std::string run;
    std::string date;
    std::string time;
    int hour = 0;
    int min = 0;
    Long64_t n_ev = 0;
    int k[64] = {};
    double rate[64] = {};
    bool sat[64] = {};
    bool dead[64] = {};
};

struct GoodInterval {
    std::string run;
    std::string start_date;
    std::string start_time;
    std::string end_date;
    std::string end_time;
    time_t t_start = 0;
    time_t t_end = 0;
    std::string status;
    std::string channels;
    std::string reason;
};

struct Episode {
    int id = 0;
    std::string run;
    std::string start_date;
    std::string start_time;
    std::string end_date;
    std::string end_time;
    int start_bin = 0;
    int end_bin = 0;
    std::vector<int> channels;
    int n18 = 0;
    int n19 = 0;
    bool is_listed = false;
    Long64_t affected_events = 0;
};

time_t parse_dt(const std::string& date_str, const std::string& time_str) {
    struct tm tm_val = {};
    int y, m, d, hh, mm;
    if (std::sscanf((date_str + " " + time_str).c_str(), "%d-%d-%d %d:%d", &y, &m, &d, &hh, &mm) == 5) {
        tm_val.tm_year = y - 1900;
        tm_val.tm_mon  = m - 1;
        tm_val.tm_mday = d;
        tm_val.tm_hour = hh;
        tm_val.tm_min  = mm;
        tm_val.tm_sec  = 0;
        return timegm(&tm_val);
    }
    return 0;
}

void hv_diag5_episodes() {
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError; // Suppress TFile recovery warnings

    // 1. Read calibration parameters from outputs/pe_calibration_v8.txt
    std::map<int, ChannelCalib> calib;
    std::ifstream fcal("outputs/pe_calibration_v8.txt");
    if (!fcal.is_open()) {
        std::cerr << "Error: outputs/pe_calibration_v8.txt not found!\n";
        return;
    }
    std::string line;
    std::vector<double> g_0_15;
    while (std::getline(fcal, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int ch; std::string pl; double p0, G, s0, s1, chi2; int ndf, status;
        if (iss >> ch >> pl >> p0 >> G >> s0 >> s1 >> chi2 >> ndf >> status) {
            calib[ch] = {p0, G};
            if (ch >= 0 && ch <= 15 && ch != 9) g_0_15.push_back(G);
        }
    }
    fcal.close();

    std::sort(g_0_15.begin(), g_0_15.end());
    double med_G_0_15 = (g_0_15[g_0_15.size()/2 - 1] + g_0_15[g_0_15.size()/2]) / 2.0;

    double thr[64];
    for (int c = 0; c < 64; ++c) {
        if (c < 32) {
            if (c == 3) thr[c] = 34.3; // dead channel
            else if (c == 9) thr[c] = calib[9].p0 + 7.5 * med_G_0_15;
            else if (calib.count(c)) thr[c] = calib[c].p0 + 7.5 * calib[c].G;
            else thr[c] = 34.3;
        } else {
            thr[c] = 3.0; // hod2019 threshold: 3 mV
        }
    }

    // 2. Read runs/good_intervals_2026.csv
    std::vector<GoodInterval> listed_intervals;
    std::ifstream fgood("runs/good_intervals_2026.csv");
    if (fgood.is_open()) {
        std::string hline;
        std::getline(fgood, hline);
        while (std::getline(fgood, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            GoodInterval gi;
            std::string s_dt, e_dt;
            if (std::getline(iss, gi.run, ',') &&
                std::getline(iss, s_dt, ',') &&
                std::getline(iss, e_dt, ',') &&
                std::getline(iss, gi.status, ',') &&
                std::getline(iss, gi.channels, ',') &&
                std::getline(iss, gi.reason, ',')) {
                gi.start_date = s_dt.substr(0, 10);
                gi.start_time = s_dt.substr(11, 5);
                gi.end_date   = e_dt.substr(0, 10);
                gi.end_time   = e_dt.substr(11, 5);
                gi.t_start    = parse_dt(gi.start_date, gi.start_time);
                gi.t_end      = parse_dt(gi.end_date, gi.end_time);
                listed_intervals.push_back(gi);
            }
        }
        fgood.close();
    }

    // 3. Process runs v7, v8, v9, v10
    struct RunInfo {
        std::string run_name;
        std::string filename;
    };
    std::vector<RunInfo> runs = {
        {"v7",  "lectura_vx1740_ext_nim_1M_v7.root"},
        {"v8",  "lectura_vx1740_ext_nim_1M_v8.root"},
        {"v9",  "lectura_vx1740_ext_nim_1M_v9.root"},
        {"v10", "lectura_vx1740_ext_nim_1M_v10.root"}
    };

    const double boot_ref_epoch = 1778886208.8; // UTC epoch of 2026-05-15T19:03:28.8-04:00
    const double tz_offset = -4.0 * 3600.0;     // UTC-4 offset

    std::map<std::string, std::vector<BinInfo>> run_bins;
    std::map<std::string, std::vector<double>> run_Rc;

    for (const auto& rinfo : runs) {
        std::string resolved = hodo::resolveRawFile(rinfo.filename);
        TFile* f = TFile::Open(resolved.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening " << resolved << "\n";
            continue;
        }
        TTree* t = (TTree*)f->Get("data_dig");
        if (!t) {
            delete f;
            continue;
        }

        Long64_t nEnt = t->GetEntries();

        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("time", 1);
        t->SetBranchStatus("timestamp", 1);
        std::vector<float>* tw = nullptr;
        Long_t ts = 0;
        t->SetBranchAddress("time", &tw);
        t->SetBranchAddress("timestamp", &ts);

        std::vector<float>* cw[64] = {};
        for (int c = 0; c < 64; ++c) {
            t->SetBranchStatus(Form("Ch%d", c), 1);
            t->SetBranchAddress(Form("Ch%d", c), &cw[c]);
        }
        t->SetCacheSize(256LL << 20);

        std::map<std::string, BinInfo> bmap;
        std::vector<std::string> border;

        for (Long64_t iev = 0; iev < nEnt; ++iev) {
            t->GetEntry(iev);
            if (!tw || tw->empty()) continue;

            double t_sec = boot_ref_epoch + ((double)ts) * 1e-6 + tz_offset;
            time_t t_local = (time_t)t_sec;
            struct tm tm_val;
            gmtime_r(&t_local, &tm_val);

            int bmin = (tm_val.tm_min / 10) * 10;
            char date_buf[32], time_buf[16];
            std::snprintf(date_buf, sizeof(date_buf), "%04d-%02d-%02d",
                          tm_val.tm_year + 1900, tm_val.tm_mon + 1, tm_val.tm_mday);
            std::snprintf(time_buf, sizeof(time_buf), "%02d:%02d",
                          tm_val.tm_hour, bmin);

            std::string bkey = std::string(date_buf) + "_" + std::string(time_buf);
            if (bmap.find(bkey) == bmap.end()) {
                BinInfo bi;
                bi.run = rinfo.run_name;
                bi.date = date_buf;
                bi.time = time_buf;
                bi.hour = tm_val.tm_hour;
                bi.min = bmin;
                bmap[bkey] = bi;
                border.push_back(bkey);
            }

            BinInfo& binfo = bmap[bkey];
            binfo.n_ev++;

            const auto& tv = *tw;
            int nsamp = (int)tv.size();

            for (int c = 0; c < 64; ++c) {
                if (c == 3) continue; // permanently dead
                if (!cw[c] || (int)cw[c]->size() != nsamp) continue;
                const auto& wf = *cw[c];

                double vmin = 1e30;
                int imin = -1;
                for (int i = 0; i < nsamp; ++i) {
                    if (tv[i] >= 100.0 && tv[i] <= 3480.0 && wf[i] < vmin) {
                        vmin = wf[i];
                        imin = i;
                    }
                }
                if (imin >= 0) {
                    double ped_loc = 0.0;
                    if (imin >= 12 && imin - 4 < nsamp) {
                        double s_loc = 0.0;
                        for (int j = -12; j <= -4; ++j) s_loc += wf[imin + j];
                        ped_loc = s_loc / 9.0;
                    } else {
                        double s_pre = 0.0; int n_pre = 0;
                        for (int i = 0; i < nsamp; ++i) {
                            if (tv[i] < 2500.0) { s_pre += wf[i]; n_pre++; }
                        }
                        ped_loc = (n_pre > 0) ? (s_pre / n_pre) : 0.0;
                    }
                    double ph = ped_loc - vmin;
                    if (ph > thr[c]) {
                        binfo.k[c]++;
                    }
                }
            }
        }

        std::vector<BinInfo> bins;
        for (const auto& bk : border) {
            BinInfo& bi = bmap[bk];
            for (int c = 0; c < 64; ++c) {
                if (bi.n_ev == 0) {
                    bi.rate[c] = 0.0;
                    bi.sat[c] = false;
                } else if (bi.k[c] >= bi.n_ev) {
                    bi.sat[c] = true;
                    bi.rate[c] = 999999.0; // SAT
                } else {
                    double occ = (double)bi.k[c] / (double)bi.n_ev;
                    bi.rate[c] = -std::log(1.0 - occ) / (3.38e-6 * 1000.0); // kHz
                    bi.sat[c] = false;
                }
            }
            bins.push_back(bi);
        }
        delete f;

        // Compute R_c: median of 10-min rates between 00:00 and 05:59 local
        std::vector<double> Rc(64, 0.0);
        for (int c = 0; c < 64; ++c) {
            std::vector<double> r_night;
            for (const auto& b : bins) {
                if (b.hour >= 0 && b.hour <= 5 && !b.sat[c]) {
                    r_night.push_back(b.rate[c]);
                }
            }
            if (!r_night.empty()) {
                std::sort(r_night.begin(), r_night.end());
                double med = (r_night.size() % 2 == 1) ? r_night[r_night.size()/2]
                             : 0.5 * (r_night[r_night.size()/2 - 1] + r_night[r_night.size()/2]);
                if (med <= 0.001) {
                    // If zero for > 50% of night bins (due to multi-day trips), take median of positive night rates
                    std::vector<double> vpos;
                    for (double v : r_night) if (v > 0.001) vpos.push_back(v);
                    if (!vpos.empty()) {
                        std::sort(vpos.begin(), vpos.end());
                        med = (vpos.size() % 2 == 1) ? vpos[vpos.size()/2]
                              : 0.5 * (vpos[vpos.size()/2 - 1] + vpos[vpos.size()/2]);
                    }
                }
                Rc[c] = med;
            } else {
                // For v9 (daytime run, no bins in 00:00-05:59): median across all bins of v9
                std::vector<double> r_all;
                for (const auto& b : bins) if (!b.sat[c]) r_all.push_back(b.rate[c]);
                if (!r_all.empty()) {
                    std::sort(r_all.begin(), r_all.end());
                    Rc[c] = (r_all.size() % 2 == 1) ? r_all[r_all.size()/2]
                            : 0.5 * (r_all[r_all.size()/2 - 1] + r_all[r_all.size()/2]);
                }
            }
        }
        run_Rc[rinfo.run_name] = Rc;

        // Flag dead channels per bin
        int nBins = (int)bins.size();
        for (int i = 0; i < nBins; ++i) {
            for (int c = 0; c < 64; ++c) {
                if (Rc[c] > 0.1) {
                    bins[i].dead[c] = (bins[i].rate[c] < 0.1 * Rc[c]);
                } else {
                    bins[i].dead[c] = false;
                }
            }
        }
        run_bins[rinfo.run_name] = bins;
    }

    // 4. Save outputs/hv_diag5_rates_10min.csv
    std::ofstream frate("outputs/hv_diag5_rates_10min.csv");
    frate << "run,date,time,ch,n_ev,k,rate\n";
    for (const auto& rinfo : runs) {
        std::string rname = rinfo.run_name;
        for (const auto& b : run_bins[rname]) {
            for (int c = 0; c < 64; ++c) {
                frate << rname << "," << b.date << "," << b.time << ","
                      << c << "," << b.n_ev << "," << b.k[c] << ",";
                if (b.sat[c]) frate << "SAT\n";
                else frate << std::fixed << std::setprecision(2) << b.rate[c] << "\n";
            }
        }
    }
    frate.close();

    // 5. Detect episodes across runs
    std::vector<Episode> all_episodes;
    int ep_counter = 1;

    for (const auto& rinfo : runs) {
        std::string rname = rinfo.run_name;
        const auto& bins = run_bins[rname];
        int nBins = (int)bins.size();

        int i = 0;
        while (i < nBins) {
            std::vector<int> candidate_chs;
            for (int c = 0; c < 64; ++c) {
                if (bins[i].dead[c] && (i == 0 || !bins[i-1].dead[c])) {
                    candidate_chs.push_back(c);
                }
            }
            if (i + 1 < nBins) {
                for (int c = 0; c < 64; ++c) {
                    if (std::find(candidate_chs.begin(), candidate_chs.end(), c) != candidate_chs.end()) continue;
                    if (bins[i+1].dead[c] && !bins[i].dead[c]) {
                        candidate_chs.push_back(c);
                    }
                }
            }

            std::vector<int> stable_dead;
            for (int c : candidate_chs) {
                int start_c = (bins[i].dead[c] && (i == 0 || !bins[i-1].dead[c])) ? i : (i + 1);
                if (start_c + 1 < nBins) {
                    if (bins[start_c].dead[c] && bins[start_c + 1].dead[c]) {
                        stable_dead.push_back(c);
                    }
                }
            }

            if (stable_dead.size() >= 3) {
                int ep_start = i;
                bool any_at_i = false;
                for (int c : stable_dead) {
                    if (bins[i].dead[c] && (i == 0 || !bins[i-1].dead[c])) any_at_i = true;
                }
                if (!any_at_i) ep_start = i + 1;

                std::vector<int> all_ep_chs = stable_dead;

                int ep_end = nBins - 1;
                for (int j = ep_start + 2; j < nBins; ++j) {
                    int n_alive = 0;
                    for (int c : all_ep_chs) {
                        if (!bins[j].dead[c]) n_alive++;
                    }
                    if ((double)n_alive >= (double)all_ep_chs.size() / 2.0) {
                        ep_end = j - 1;
                        break;
                    }
                }

                // Must stay dead for >= 2 bins
                if (ep_end >= ep_start + 1) {
                    Episode ep;
                    ep.id = ep_counter++;
                    ep.run = rinfo.run_name;
                    ep.start_bin = ep_start;
                    ep.end_bin = ep_end;
                    ep.start_date = bins[ep_start].date;
                    ep.start_time = bins[ep_start].time;
                    ep.end_date = bins[ep_end].date;
                    ep.end_time = bins[ep_end].time;
                    ep.channels = all_ep_chs;
                    std::sort(ep.channels.begin(), ep.channels.end());
                    for (int c : ep.channels) {
                        if (c < 32) ep.n18++;
                        else ep.n19++;
                    }

                    // Check overlap with runs/good_intervals_2026.csv
                    time_t ep_t_start = parse_dt(ep.start_date, ep.start_time);
                    time_t ep_t_end   = parse_dt(ep.end_date, ep.end_time) + 600; // end of final 10-min bin

                    ep.is_listed = false;
                    for (const auto& gi : listed_intervals) {
                        if (gi.run == ep.run) {
                            if (ep_t_start <= gi.t_end + 60 && ep_t_end >= gi.t_start) {
                                ep.is_listed = true;
                                break;
                            }
                        }
                    }

                    // Calculate affected events for NOT listed episodes
                    ep.affected_events = 0;
                    if (!ep.is_listed) {
                        for (int b = ep.start_bin; b <= ep.end_bin; ++b) {
                            ep.affected_events += bins[b].n_ev;
                        }
                    }

                    all_episodes.push_back(ep);
                    i = ep_end + 1;
                } else {
                    i++;
                }
            } else {
                i++;
            }
        }
    }

    // 6. Save outputs/hv_diag5_episodes.csv
    std::ofstream fep("outputs/hv_diag5_episodes.csv");
    fep << "run,start,end,n18,n19,channels,status,affected_events\n";
    for (const auto& ep : all_episodes) {
        fep << ep.run << "," << ep.start_date << " " << ep.start_time << ","
            << ep.end_date << " " << ep.end_time << ","
            << ep.n18 << "," << ep.n19 << ",\"";
        for (size_t ci = 0; ci < ep.channels.size(); ++ci) {
            fep << ep.channels[ci] << (ci + 1 == ep.channels.size() ? "" : ",");
        }
        fep << "\"," << (ep.is_listed ? "listed" : "NOT listed") << "," << ep.affected_events << "\n";
    }
    fep.close();

    // -------------------------------------------------------------------------
    // REPORT SECTION 1: R_c per channel (one line per run: "c:R_c" for c = 0..63)
    // -------------------------------------------------------------------------
    std::cout << "========================================================================================\n";
    std::cout << "OUTPUT 1: Reference Rate R_c per channel (kHz, 00:00 - 05:59 Local Median)\n";
    std::cout << "========================================================================================\n";
    for (const auto& rinfo : runs) {
        std::string rname = rinfo.run_name;
        const auto& Rc = run_Rc[rname];
        std::cout << "run " << rname << ":\n";
        for (int c = 0; c < 64; ++c) {
            std::cout << c << ":" << std::fixed << std::setprecision(2) << Rc[c] << (c == 63 ? "" : " ");
        }
        std::cout << "\n\n";
    }

    // -------------------------------------------------------------------------
    // REPORT SECTION 2: Table of episodes
    // Columns EXACTLY: run  start  end  n18  n19  channels
    // -------------------------------------------------------------------------
    std::cout << "========================================================================================\n";
    std::cout << "OUTPUT 2: Table of Episodes\n";
    std::cout << "========================================================================================\n";
    std::cout << std::left << std::setw(5) << "run"
              << std::setw(18) << "start"
              << std::setw(18) << "end"
              << std::right << std::setw(5) << "n18"
              << std::setw(5) << "n19"
              << "  " << std::left << "channels\n";
    std::cout << std::string(80, '-') << "\n";
    for (const auto& ep : all_episodes) {
        std::string s_str = ep.start_date + " " + ep.start_time;
        std::string e_str = ep.end_date + " " + ep.end_time;
        std::string chs_str;
        for (size_t ci = 0; ci < ep.channels.size(); ++ci) {
            chs_str += std::to_string(ep.channels[ci]);
            if (ci + 1 < ep.channels.size()) chs_str += ",";
        }
        std::cout << std::left << std::setw(5) << ep.run
                  << std::setw(18) << s_str
                  << std::setw(18) << e_str
                  << std::right << std::setw(5) << ep.n18
                  << std::setw(5) << ep.n19
                  << "  " << std::left << chs_str << "\n";
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // REPORT SECTION 3 & 4: Listed / NOT listed status & Affected events
    // -------------------------------------------------------------------------
    std::cout << "========================================================================================\n";
    std::cout << "OUTPUT 3 & 4: Episode Overlap with good_intervals_2026.csv and Affected Events\n";
    std::cout << "========================================================================================\n";
    std::cout << std::left << std::setw(5) << "id"
              << std::setw(5) << "run"
              << std::setw(18) << "start"
              << std::setw(18) << "end"
              << std::setw(12) << "status"
              << std::right << std::setw(16) << "affected_events" << "\n";
    std::cout << std::string(74, '-') << "\n";
    for (const auto& ep : all_episodes) {
        std::string s_str = ep.start_date + " " + ep.start_time;
        std::string e_str = ep.end_date + " " + ep.end_time;
        std::cout << std::left << std::setw(5) << ep.id
                  << std::setw(5) << ep.run
                  << std::setw(18) << s_str
                  << std::setw(18) << e_str
                  << std::setw(12) << (ep.is_listed ? "listed" : "NOT listed")
                  << std::right << std::setw(16) << ep.affected_events << "\n";
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // REPORT SECTION 5: Mean rate ratio of every hod2018 channel that is NOT dead
    // One line "c:ratio" for c = 0..31
    // -------------------------------------------------------------------------
    std::cout << "========================================================================================\n";
    std::cout << "OUTPUT 5: Mean Rate Ratio (rate in episode / R_c) for hod2018 Channels NOT Dead\n";
    std::cout << "========================================================================================\n";
    for (const auto& ep : all_episodes) {
        std::string rname = ep.run;
        const auto& bins = run_bins[rname];
        const auto& Rc = run_Rc[rname];

        std::set<int> dead_set(ep.channels.begin(), ep.channels.end());
        int n_ep_bins = ep.end_bin - ep.start_bin + 1;

        std::cout << "Episode " << ep.id << " (" << ep.run << " "
                  << ep.start_date << " " << ep.start_time << " to "
                  << ep.end_date << " " << ep.end_time << "):\n";

        bool first = true;
        for (int c = 0; c < 32; ++c) {
            if (c == 3 || dead_set.count(c)) continue;
            if (Rc[c] <= 0.001) continue;

            double sum_rate = 0.0;
            for (int b = ep.start_bin; b <= ep.end_bin; ++b) {
                sum_rate += bins[b].rate[c];
            }
            double mean_rate = sum_rate / (double)n_ep_bins;
            double ratio = mean_rate / Rc[c];

            if (!first) std::cout << " ";
            std::cout << c << ":" << std::fixed << std::setprecision(2) << ratio;
            first = false;
        }
        std::cout << "\n\n";
    }
}

