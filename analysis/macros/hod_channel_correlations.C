// analysis/macros/hod_channel_correlations.C
// Measures channel-level and plane-level correlations in hodoscopes.
// Run on v8 and on v10, all entries, READ only.
//
// Usage: root -l -b -q 'analysis/macros/hod_channel_correlations.C("lectura_vx1740_ext_nim_1M_v8.root", "outputs/hod_channel_correlations_v8")'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
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

static inline int ch_to_plane(int c) {
    if (c >= 0 && c <= 15) return 0;   // 18X
    if (c >= 16 && c <= 31) return 1;  // 18Y
    if (c >= 32 && c <= 47) return 2;  // 19Y
    if (c >= 48 && c <= 63) return 3;  // 19X
    return -1;
}

static inline const char* plane_name(int p) {
    switch (p) {
        case 0: return "18X";
        case 1: return "18Y";
        case 2: return "19Y";
        case 3: return "19X";
        default: return "UNK";
    }
}

struct EvHitInfo {
    bool ch_hit_S[64] = {};
    float ch_t_S[64] = {};
    bool ch_hit_C[64] = {};
    bool pl_hit_S[4] = {};
    float pl_t_S[4] = {};
};

void hod_channel_correlations(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root",
                              const char* outdir_arg = "") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    fs::path outdir;
    if (outdir_arg && std::string(outdir_arg).length() > 0) {
        outdir = fs::path(outdir_arg);
    } else {
        std::time_t now = std::time(nullptr);
        char ds[32];
        std::strftime(ds, sizeof(ds), "%Y%m%d_%H%M%S", std::localtime(&now));
        outdir = fs::path("outputs") / ("hod_channel_correlations_" + std::string(ds));
    }
    fs::create_directories(outdir);

    std::cout << "=== hod_channel_correlations ===\n";
    std::cout << "Input: " << fname << "\n";
    std::cout << "Output dir: " << outdir << "\n";

    TFile* fin = TFile::Open(fname, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open file " << fname << "\n";
        return;
    }
    TTree* t = static_cast<TTree*>(fin->Get("data_dig"));
    if (!t) {
        std::cerr << "ERROR: data_dig tree not found in " << fname << "\n";
        fin->Close();
        return;
    }

    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("time", 1);
    std::vector<float>* tw = nullptr;
    t->SetBranchAddress("time", &tw);

    std::vector<float>* cw[64] = {nullptr};
    for (int c = 0; c < 64; ++c) {
        std::string bn = Form("Ch%d", c);
        t->SetBranchStatus(bn.c_str(), 1);
        t->SetBranchAddress(bn.c_str(), &cw[c]);
    }

    const Long64_t nEnt = t->GetEntries();
    std::cout << "Total entries: " << nEnt << "\n";

    // 1A accumulators
    long long N_S[64] = {};
    long long N_C[64] = {};
    long long N_plane_S[4] = {};
    long long N_plane_C[4] = {};

    // 1D average waveforms (-10..+30 samples -> 41 points)
    // Offset k in [-10, +30] -> index idx = k + 10
    double sum_19Y_tag[41] = {};
    long long cnt_19Y_tag[41] = {};
    double sum_19Xother_tag[41] = {};
    long long cnt_19Xother_tag[41] = {};
    double sum_19Y_ref[41] = {};
    long long cnt_19Y_ref[41] = {};

    std::vector<EvHitInfo> ev_records;
    ev_records.reserve(nEnt);

    // Main event loop
    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        if (iev % 30000 == 0) {
            std::cout << "  Event " << iev << "/" << nEnt << "\r" << std::flush;
        }
        t->GetEntry(iev);
        if (!tw || tw->empty()) {
            ev_records.emplace_back();
            continue;
        }

        const auto& tv = *tw;
        const int n_samp = static_cast<int>(tv.size());

        EvHitInfo ev_info;

        // Channel analysis
        double ph_S_all[64] = {};
        double t_S_all[64] = {};
        int imin_S_all[64] = {};
        double ped_all[64] = {};

        for (int c = 0; c < 64; ++c) {
            if (!cw[c] || static_cast<int>(cw[c]->size()) != n_samp) continue;
            const auto& wf = *cw[c];

            // ped = mean of samples with t < 2500 ns
            double s_ped = 0;
            int n_ped = 0;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] < 2500.0) {
                    s_ped += wf[i];
                    n_ped++;
                }
            }
            double ped = (n_ped > 0) ? (s_ped / n_ped) : 0.0;
            ped_all[c] = ped;

            // Signal window S = [4000, 4400] ns
            double min_S = 1e30;
            int imin_S = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= 4000.0 && tv[i] <= 4400.0 && wf[i] < min_S) {
                    min_S = wf[i];
                    imin_S = i;
                }
            }
            double ph_S = (imin_S >= 0) ? (ped - min_S) : 0.0;
            double t_S = (imin_S >= 0) ? tv[imin_S] : 0.0;
            ph_S_all[c] = ph_S;
            t_S_all[c] = t_S;
            imin_S_all[c] = imin_S;

            // Control window C = [6000, 6400] ns
            double min_C = 1e30;
            int imin_C = -1;
            for (int i = 0; i < n_samp; ++i) {
                if (tv[i] >= 6000.0 && tv[i] <= 6400.0 && wf[i] < min_C) {
                    min_C = wf[i];
                    imin_C = i;
                }
            }
            double ph_C = (imin_C >= 0) ? (ped - min_C) : 0.0;

            const double thr = (c < 32) ? 50.0 : 35.0;
            bool hit_S = (ph_S > thr);
            bool hit_C = (ph_C > thr);

            if (hit_S) {
                N_S[c]++;
                ev_info.ch_hit_S[c] = true;
                ev_info.ch_t_S[c] = static_cast<float>(t_S);
            }
            if (hit_C) {
                N_C[c]++;
                ev_info.ch_hit_C[c] = true;
            }
        }

        // Plane hit determinations
        // Window S
        for (int p = 0; p < 4; ++p) {
            int ch_start = p * 16;
            int ch_end = ch_start + 16;
            double max_ph = -1e30;
            int best_ch = -1;
            bool any_hit = false;

            for (int c = ch_start; c < ch_end; ++c) {
                if (ev_info.ch_hit_S[c]) any_hit = true;
                if (ph_S_all[c] > max_ph) {
                    max_ph = ph_S_all[c];
                    best_ch = c;
                }
            }
            if (any_hit) {
                ev_info.pl_hit_S[p] = true;
                ev_info.pl_t_S[p] = static_cast<float>(t_S_all[best_ch]);
                N_plane_S[p]++;
            }
        }

        // Window C plane hit
        for (int p = 0; p < 4; ++p) {
            int ch_start = p * 16;
            int ch_end = ch_start + 16;
            for (int c = ch_start; c < ch_end; ++c) {
                if (ev_info.ch_hit_C[c]) {
                    N_plane_C[p]++;
                    break;
                }
            }
        }

        // 1D Average waveforms
        // Tag: 19X channel with the largest ph_S > 35 mV (ch 48-63)
        double max_ph_19x = -1e30;
        int tag_ch_19x = -1;
        for (int c = 48; c < 64; ++c) {
            if (ph_S_all[c] > 35.0 && ph_S_all[c] > max_ph_19x) {
                max_ph_19x = ph_S_all[c];
                tag_ch_19x = c;
            }
        }

        if (tag_ch_19x >= 0) {
            int imin_tag = imin_S_all[tag_ch_19x];
            // (i) 16 19Y channels (32-47) with ph_S <= 35 mV
            for (int c = 32; c < 48; ++c) {
                if (ph_S_all[c] <= 35.0 && cw[c]) {
                    const auto& wf = *cw[c];
                    for (int k = -10; k <= 30; ++k) {
                        int s = imin_tag + k;
                        if (s >= 0 && s < n_samp) {
                            sum_19Y_tag[k + 10] += (wf[s] - ped_all[c]);
                            cnt_19Y_tag[k + 10]++;
                        }
                    }
                }
            }
            // (ii) Other 15 19X channels (48-63, c != tag_ch_19x) with ph_S <= 35 mV
            for (int c = 48; c < 64; ++c) {
                if (c != tag_ch_19x && ph_S_all[c] <= 35.0 && cw[c]) {
                    const auto& wf = *cw[c];
                    for (int k = -10; k <= 30; ++k) {
                        int s = imin_tag + k;
                        if (s >= 0 && s < n_samp) {
                            sum_19Xother_tag[k + 10] += (wf[s] - ped_all[c]);
                            cnt_19Xother_tag[k + 10]++;
                        }
                    }
                }
            }
        }

        // Reference: events with no channel 32-63 hit in S (ph_S <= 35 mV for all ch 32..63)
        bool any_hod19_hit = false;
        for (int c = 32; c < 64; ++c) {
            if (ph_S_all[c] > 35.0) {
                any_hod19_hit = true;
                break;
            }
        }
        if (!any_hod19_hit) {
            const int s_ref = 262;
            for (int c = 32; c < 48; ++c) {
                if (cw[c]) {
                    const auto& wf = *cw[c];
                    for (int k = -10; k <= 30; ++k) {
                        int s = s_ref + k;
                        if (s >= 0 && s < n_samp) {
                            sum_19Y_ref[k + 10] += (wf[s] - ped_all[c]);
                            cnt_19Y_ref[k + 10]++;
                        }
                    }
                }
            }
        }

        ev_records.push_back(ev_info);
    }
    std::cout << "  Event " << nEnt << "/" << nEnt << " done.\n";
    fin->Close();

    // 1B. Plane pairs (window S, |dt| < 100 ns)
    // 6 pairs: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    struct PlanePairResult {
        int pA = 0, pB = 0;
        long long N_A = 0, N_B = 0;
        double N_indep = 0;
        long long N_k[7] = {}; // k = -3..+3 -> index 0..6
    };
    std::vector<PlanePairResult> pp_results(6);
    const int pair_indices[6][2] = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };

    for (int pidx = 0; pidx < 6; ++pidx) {
        int pA = pair_indices[pidx][0];
        int pB = pair_indices[pidx][1];
        pp_results[pidx].pA = pA;
        pp_results[pidx].pB = pB;
        pp_results[pidx].N_A = N_plane_S[pA];
        pp_results[pidx].N_B = N_plane_S[pB];
        pp_results[pidx].N_indep = (nEnt > 0) ? (static_cast<double>(N_plane_S[pA]) * N_plane_S[pB] / nEnt) : 0.0;

        for (int k = -3; k <= 3; ++k) {
            long long cnt = 0;
            int i_start = (k < 0) ? -k : 0;
            int i_end = (k > 0) ? static_cast<int>(nEnt - 1 - k) : static_cast<int>(nEnt - 1);
            for (int i = i_start; i <= i_end; ++i) {
                int j = i + k;
                if (ev_records[i].pl_hit_S[pA] && ev_records[j].pl_hit_S[pB]) {
                    if (std::abs(ev_records[i].pl_t_S[pA] - ev_records[j].pl_t_S[pB]) < 100.0) {
                        cnt++;
                    }
                }
            }
            pp_results[pidx].N_k[k + 3] = cnt;
        }
    }

    // 1C. Channel pairs (window S, |dt| < 100 ns)
    std::vector<std::vector<long long>> Same(64, std::vector<long long>(64, 0));
    std::vector<std::vector<long long>> Mixed(64, std::vector<long long>(64, 0));

    // Same: both a and b hit at event i, a != b
    for (size_t i = 0; i < ev_records.size(); ++i) {
        std::vector<int> hits;
        for (int c = 0; c < 64; ++c) {
            if (ev_records[i].ch_hit_S[c]) hits.push_back(c);
        }
        for (size_t m = 0; m < hits.size(); ++m) {
            int a = hits[m];
            for (size_t n = 0; n < hits.size(); ++n) {
                if (m == n) continue;
                int b = hits[n];
                if (std::abs(ev_records[i].ch_t_S[a] - ev_records[i].ch_t_S[b]) < 100.0) {
                    Same[a][b]++;
                }
            }
        }
    }

    // Mixed: a hit at event i, b hit at event i+1
    for (size_t i = 0; i + 1 < ev_records.size(); ++i) {
        std::vector<int> hits_i, hits_next;
        for (int c = 0; c < 64; ++c) {
            if (ev_records[i].ch_hit_S[c]) hits_i.push_back(c);
            if (ev_records[i + 1].ch_hit_S[c]) hits_next.push_back(c);
        }
        for (int a : hits_i) {
            for (int b : hits_next) {
                if (std::abs(ev_records[i].ch_t_S[a] - ev_records[i + 1].ch_t_S[b]) < 100.0) {
                    Mixed[a][b]++;
                }
            }
        }
    }

    // Write CSVs for 64x64 matrices
    {
        std::ofstream fc_same(outdir / "same_matrix_64x64.csv");
        fc_same << "ch";
        for (int c = 0; c < 64; ++c) fc_same << "," << c;
        fc_same << "\n";
        for (int a = 0; a < 64; ++a) {
            fc_same << a;
            for (int b = 0; b < 64; ++b) fc_same << "," << Same[a][b];
            fc_same << "\n";
        }
    }
    {
        std::ofstream fc_mixed(outdir / "mixed_matrix_64x64.csv");
        fc_mixed << "ch";
        for (int c = 0; c < 64; ++c) fc_mixed << "," << c;
        fc_mixed << "\n";
        for (int a = 0; a < 64; ++a) {
            fc_mixed << a;
            for (int b = 0; b < 64; ++b) fc_mixed << "," << Mixed[a][b];
            fc_mixed << "\n";
        }
    }

    // 8x8 Group matrices: distinct channels only (a != b)
    long long GrpSame[8][8] = {};
    long long GrpMixed[8][8] = {};
    for (int g1 = 0; g1 < 8; ++g1) {
        for (int g2 = 0; g2 < 8; ++g2) {
            for (int a = g1 * 8; a < (g1 + 1) * 8; ++a) {
                for (int b = g2 * 8; b < (g2 + 1) * 8; ++b) {
                    if (a != b) {
                        GrpSame[g1][g2] += Same[a][b];
                        GrpMixed[g1][g2] += Mixed[a][b];
                    }
                }
            }
        }
    }

    // 1D Average waveforms computation
    std::vector<double> avg_19Y(41, 0.0), avg_19Xother(41, 0.0), avg_19Y_ref_val(41, 0.0);
    std::vector<double> v_offset(41);
    for (int k = -10; k <= 30; ++k) {
        int idx = k + 10;
        v_offset[idx] = k;
        avg_19Y[idx] = (cnt_19Y_tag[idx] > 0) ? (sum_19Y_tag[idx] / cnt_19Y_tag[idx]) : 0.0;
        avg_19Xother[idx] = (cnt_19Xother_tag[idx] > 0) ? (sum_19Xother_tag[idx] / cnt_19Xother_tag[idx]) : 0.0;
        avg_19Y_ref_val[idx] = (cnt_19Y_ref[idx] > 0) ? (sum_19Y_ref[idx] / cnt_19Y_ref[idx]) : 0.0;
    }

    // Plot 1D average waveforms
    {
        TCanvas c1("c_1d", "Average Waveforms Aligned on Hod2019X Tag", 900, 600);
        c1.SetGrid();
        c1.SetLeftMargin(0.12);
        c1.SetBottomMargin(0.12);

        TGraph g_19Y(41, v_offset.data(), avg_19Y.data());
        g_19Y.SetLineColor(kBlue + 1);
        g_19Y.SetLineWidth(2);
        g_19Y.SetTitle("Average Waveforms Aligned on Hod2019X Tag;Offset [samples];Average (wf - ped) [mV]");

        TGraph g_19Xother(41, v_offset.data(), avg_19Xother.data());
        g_19Xother.SetLineColor(kRed + 1);
        g_19Xother.SetLineWidth(2);

        TGraph g_19Y_ref(41, v_offset.data(), avg_19Y_ref_val.data());
        g_19Y_ref.SetLineColor(kBlack);
        g_19Y_ref.SetLineWidth(2);
        g_19Y_ref.SetLineStyle(2);

        g_19Y.Draw("AL");
        g_19Xother.Draw("L SAME");
        g_19Y_ref.Draw("L SAME");

        TLegend leg(0.55, 0.72, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(&g_19Y, "19Y channels (tag event)", "l");
        leg.AddEntry(&g_19Xother, "Other 19X channels (tag event)", "l");
        leg.AddEntry(&g_19Y_ref, "19Y channels (reference event)", "l");
        leg.Draw();

        c1.SaveAs((outdir / "average_waveforms_1D.png").string().c_str());
    }

    // Summary text printer
    auto print_summary = [&](std::ostream& os) {
        os << "=== Hodoscope Channel and Plane Correlations ===\n\n";
        os << "Input file: " << fname << "\n";
        os << "Total events: " << nEnt << "\n";
        os << "Output directory: " << outdir << "\n\n";

        os << "--- 1A. Per-Channel Hit Counts (Signal window S [4000, 4400] ns, Control window C [6000, 6400] ns) ---\n";
        os << std::setw(6) << std::left << "ch"
           << std::setw(8) << "plane"
           << std::setw(12) << "N_S"
           << std::setw(12) << "N_C"
           << std::setw(12) << "N_S - N_C" << "\n";
        os << "--------------------------------------------------\n";
        for (int c = 0; c < 64; ++c) {
            os << std::setw(6) << std::left << c
               << std::setw(8) << plane_name(ch_to_plane(c))
               << std::setw(12) << N_S[c]
               << std::setw(12) << N_C[c]
               << std::setw(12) << (N_S[c] - N_C[c]) << "\n";
        }
        os << "\n--- 1A (cont). Per-Plane Hit Counts ---\n";
        os << std::setw(8) << std::left << "Plane"
           << std::setw(12) << "N_S"
           << std::setw(12) << "N_C"
           << std::setw(12) << "N_S - N_C" << "\n";
        os << "--------------------------------------------\n";
        for (int p = 0; p < 4; ++p) {
            os << std::setw(8) << std::left << plane_name(p)
               << std::setw(12) << N_plane_S[p]
               << std::setw(12) << N_plane_C[p]
               << std::setw(12) << (N_plane_S[p] - N_plane_C[p]) << "\n";
        }

        os << "\n--- 1B. Plane Pairs (Window S, |dt| < 100 ns) ---\n";
        os << std::setw(12) << std::left << "Pair"
           << std::setw(10) << "N_A"
           << std::setw(10) << "N_B"
           << std::setw(12) << "N_indep"
           << std::setw(8) << "N(-3)"
           << std::setw(8) << "N(-2)"
           << std::setw(8) << "N(-1)"
           << std::setw(8) << "N(0)"
           << std::setw(8) << "N(+1)"
           << std::setw(8) << "N(+2)"
           << std::setw(8) << "N(+3)" << "\n";
        os << "--------------------------------------------------------------------------------------------------\n";
        for (int pidx = 0; pidx < 6; ++pidx) {
            const auto& r = pp_results[pidx];
            std::string pname = std::string(plane_name(r.pA)) + "-" + plane_name(r.pB);
            os << std::setw(12) << std::left << pname
               << std::setw(10) << r.N_A
               << std::setw(10) << r.N_B
               << std::setw(12) << std::fixed << std::setprecision(1) << r.N_indep;
            for (int k = 0; k < 7; ++k) {
                os << std::setw(8) << r.N_k[k];
            }
            os << "\n";
        }

        os << "\n--- 1C. 8x8 Group Matrices (Distinct channels only, a != b) ---\n";
        os << "\nGroup Matrix: Same\n";
        os << "grp ";
        for (int g = 0; g < 8; ++g) os << std::setw(8) << ("G" + std::to_string(g));
        os << "\n--------------------------------------------------------------------\n";
        for (int g1 = 0; g1 < 8; ++g1) {
            os << "G" << g1 << "  ";
            for (int g2 = 0; g2 < 8; ++g2) {
                os << std::setw(8) << GrpSame[g1][g2];
            }
            os << "\n";
        }

        os << "\nGroup Matrix: Mixed\n";
        os << "grp ";
        for (int g = 0; g < 8; ++g) os << std::setw(8) << ("G" + std::to_string(g));
        os << "\n--------------------------------------------------------------------\n";
        for (int g1 = 0; g1 < 8; ++g1) {
            os << "G" << g1 << "  ";
            for (int g2 = 0; g2 < 8; ++g2) {
                os << std::setw(8) << GrpMixed[g1][g2];
            }
            os << "\n";
        }

        os << "\nGroup Matrix: Same / Mixed\n";
        os << "grp ";
        for (int g = 0; g < 8; ++g) os << std::setw(8) << ("G" + std::to_string(g));
        os << "\n--------------------------------------------------------------------\n";
        for (int g1 = 0; g1 < 8; ++g1) {
            os << "G" << g1 << "  ";
            for (int g2 = 0; g2 < 8; ++g2) {
                if (GrpMixed[g1][g2] > 0) {
                    double r = static_cast<double>(GrpSame[g1][g2]) / GrpMixed[g1][g2];
                    os << std::setw(8) << std::fixed << std::setprecision(2) << r;
                } else {
                    os << std::setw(8) << "n/a";
                }
            }
            os << "\n";
        }

        os << "\n--- 1D. Average Waveforms Aligned on Hod2019X Tag (ph_S > 35 mV) ---\n";
        os << std::setw(8) << std::left << "offset"
           << std::setw(16) << "avg_19Y_tag"
           << std::setw(20) << "avg_19Xother_tag"
           << std::setw(16) << "avg_19Y_ref" << "\n";
        os << "------------------------------------------------------------\n";
        for (int k = -10; k <= 30; ++k) {
            int idx = k + 10;
            os << std::setw(8) << std::left << k
               << std::setw(16) << std::fixed << std::setprecision(3) << avg_19Y[idx]
               << std::setw(20) << avg_19Xother[idx]
               << std::setw(16) << avg_19Y_ref_val[idx] << "\n";
        }
    };

    print_summary(std::cout);
    {
        std::ofstream sf(outdir / "summary.txt");
        print_summary(sf);
    }
    std::cout << "\nResults written to: " << outdir << "\n";
}
