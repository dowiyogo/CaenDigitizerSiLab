// analysis/macros/hod_hit_structure.C
// Analysis of hodoscope hit timing, multiplicity, pulse heights, and channel correlations.
//
// Usage:
//   root -l -b -q 'analysis/macros/hod_hit_structure.C("lectura_vx1740_ext_nim_1M_v8.root", "", false)'
//   root -l -b -q 'analysis/macros/hod_hit_structure.C("lectura_vx1740_ext_nim_1M_v8.root", "", true)'

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static const std::array<int, 16> kH19X = {53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56};
static const std::array<int, 16> kH19Y = {37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40};

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

void hod_hit_structure(const char* fname = "lectura_vx1740_ext_nim_1M_v8.root",
                       const char* outdir_arg = "",
                       bool counts_only = false) {
    gROOT->SetBatch(kTRUE);

    fs::path outdir;
    if (outdir_arg && std::string(outdir_arg).length() > 0) {
        outdir = fs::path(outdir_arg);
    } else {
        std::time_t now = std::time(nullptr);
        char ds[32];
        std::strftime(ds, sizeof(ds), "%Y%m%d_%H%M%S", std::localtime(&now));
        outdir = fs::path("outputs") / ("hod_hit_structure_" + std::string(ds));
    }

    if (!counts_only) {
        fs::create_directories(outdir);
    }

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

    // Extract short run name (e.g. v1..v10)
    std::string sname = fname;
    std::string run_label = "v1";
    for (int v = 10; v >= 1; --v) {
        std::string tag = "v" + std::to_string(v);
        if (sname.find("_" + tag + ".root") != std::string::npos ||
            sname.find(tag + ".root") != std::string::npos) {
            run_label = tag;
            break;
        }
    }

    // Step 2 / Counts-only accumulators (events with >= 1 hit in S/C)
    long long ev_hit_18X_S = 0, ev_hit_18X_C = 0;
    long long ev_hit_18Y_S = 0, ev_hit_18Y_C = 0;
    long long ev_hit_G4_S = 0, ev_hit_G4_C = 0;
    long long ev_hit_G5_S = 0, ev_hit_G5_C = 0;
    long long ev_hit_G6_S = 0, ev_hit_G6_C = 0;
    long long ev_hit_G7_S = 0, ev_hit_G7_C = 0;

    // Step 1a accumulators: 25 bins of 16 ns across S [4000, 4400] and C [6000, 6400]
    // t_center_S = 4008 + 16 * b
    long long hits_18X_S[25] = {};
    long long hits_18X_C[25] = {};
    long long hits_18Y_S[25] = {};
    long long hits_18Y_C[25] = {};
    long long hits_G4_S[25] = {};
    long long hits_G5_S[25] = {};
    long long hits_G6_S[25] = {};
    long long hits_G7_S[25] = {};

    // Step 1b accumulators: multiplicity of hod2019 (ch 32..63 in S)
    // hit channels (1, 2, 3, 4, 5+) vs distinct groups hit (1, 2, 3, 4)
    // Row 0..4: n_ch = 1, 2, 3, 4, >=5
    // Col 0..3: n_grp = 1, 2, 3, 4
    long long mult_table[5][4] = {};

    // Step 1c accumulators: 25 mV bins from 35 to 1035 mV (40 bins) + overflow (bin 40)
    // Groups G4..G7 (index 0..3)
    long long ph_hist[4][41] = {};
    double lowest_wf_imin[4] = {1e30, 1e30, 1e30, 1e30};

    // Main event loop
    for (Long64_t iev = 0; iev < nEnt; ++iev) {
        if (!counts_only && iev % 30000 == 0) {
            std::cout << "  Event " << iev << "/" << nEnt << "\r" << std::flush;
        }
        t->GetEntry(iev);
        if (!tw || tw->empty()) continue;

        const auto& tv = *tw;
        const int n_samp = static_cast<int>(tv.size());

        bool pl_has_S[4] = {}; // 18X, 18Y, 19Y, 19X
        bool pl_has_C[4] = {};
        bool grp_has_S[8] = {}; // G0..G7
        bool grp_has_C[8] = {};

        int n_hod19_hit_ch = 0;

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
            double t_C = (imin_C >= 0) ? tv[imin_C] : 0.0;

            const double thr = (c < 32) ? 50.0 : 35.0;
            bool hit_S = (ph_S > thr);
            bool hit_C = (ph_C > thr);

            int pl = ch_to_plane(c);
            int grp = c / 8;

            if (hit_S) {
                pl_has_S[pl] = true;
                grp_has_S[grp] = true;
                if (c >= 32) n_hod19_hit_ch++;

                if (!counts_only) {
                    // 1a Hit times binning
                    int b_S = static_cast<int>((t_S - 4000.0) / 16.0);
                    if (b_S < 0) b_S = 0;
                    if (b_S > 24) b_S = 24;

                    if (pl == 0) hits_18X_S[b_S]++;
                    else if (pl == 1) hits_18Y_S[b_S]++;

                    if (grp == 4) hits_G4_S[b_S]++;
                    else if (grp == 5) hits_G5_S[b_S]++;
                    else if (grp == 6) hits_G6_S[b_S]++;
                    else if (grp == 7) hits_G7_S[b_S]++;

                    // 1c ph spectra and lowest wf
                    if (c >= 32) {
                        int g_idx = grp - 4; // 0..3 for G4..G7
                        int b_ph = static_cast<int>((ph_S - 35.0) / 25.0);
                        if (b_ph < 0) b_ph = 0;
                        if (b_ph >= 40) b_ph = 40; // overflow
                        ph_hist[g_idx][b_ph]++;

                        if (min_S < lowest_wf_imin[g_idx]) {
                            lowest_wf_imin[g_idx] = min_S;
                        }
                    }
                }
            }

            if (hit_C) {
                pl_has_C[pl] = true;
                grp_has_C[grp] = true;

                if (!counts_only) {
                    int b_C = static_cast<int>((t_C - 6000.0) / 16.0);
                    if (b_C < 0) b_C = 0;
                    if (b_C > 24) b_C = 24;

                    if (pl == 0) hits_18X_C[b_C]++;
                    else if (pl == 1) hits_18Y_C[b_C]++;
                }
            }
        }

        // Per-event accumulation for counts
        if (pl_has_S[0]) ev_hit_18X_S++;
        if (pl_has_C[0]) ev_hit_18X_C++;
        if (pl_has_S[1]) ev_hit_18Y_S++;
        if (pl_has_C[1]) ev_hit_18Y_C++;

        if (grp_has_S[4]) ev_hit_G4_S++;
        if (grp_has_C[4]) ev_hit_G4_C++;
        if (grp_has_S[5]) ev_hit_G5_S++;
        if (grp_has_C[5]) ev_hit_G5_C++;
        if (grp_has_S[6]) ev_hit_G6_S++;
        if (grp_has_C[6]) ev_hit_G6_C++;
        if (grp_has_S[7]) ev_hit_G7_S++;
        if (grp_has_C[7]) ev_hit_G7_C++;

        // 1b multiplicity
        if (!counts_only && n_hod19_hit_ch >= 1) {
            int n_grp_hit = 0;
            for (int g = 4; g <= 7; ++g) {
                if (grp_has_S[g]) n_grp_hit++;
            }
            int r_idx = (n_hod19_hit_ch >= 5) ? 4 : (n_hod19_hit_ch - 1);
            int c_idx = (n_grp_hit >= 1 && n_grp_hit <= 4) ? (n_grp_hit - 1) : 0;
            mult_table[r_idx][c_idx]++;
        }
    }
    fin->Close();

    if (counts_only) {
        // Print single formatted row for Step 2
        // Format: run, N_events, 18X_S, 18X_C, 18Y_S, 18Y_C, G4_S, G4_C, G5_S, G5_C, G6_S, G6_C, G7_S, G7_C
        std::cout << std::setw(6) << std::left << run_label
                  << std::setw(10) << nEnt
                  << std::setw(8) << ev_hit_18X_S
                  << std::setw(8) << ev_hit_18X_C
                  << std::setw(8) << ev_hit_18Y_S
                  << std::setw(8) << ev_hit_18Y_C
                  << std::setw(8) << ev_hit_G4_S
                  << std::setw(8) << ev_hit_G4_C
                  << std::setw(8) << ev_hit_G5_S
                  << std::setw(8) << ev_hit_G5_C
                  << std::setw(8) << ev_hit_G6_S
                  << std::setw(8) << ev_hit_G6_C
                  << std::setw(8) << ev_hit_G7_S
                  << std::setw(8) << ev_hit_G7_C << "\n";
        return;
    }

    std::cout << "  Event " << nEnt << "/" << nEnt << " done.\n";

    // 1d: Channel pairs from CSVs of hod_channel_correlations
    struct TopPair {
        int a = 0, b = 0;
        std::string grp_a, grp_b;
        std::string pl_a, pl_b;
        int k_a = -1, k_b = -1;
        long long Same = 0;
        long long Mixed = 0;
    };
    std::vector<TopPair> top_pairs;

    // Inverse maps for 19X / 19Y
    std::array<int, 64> inv19x{}, inv19y{};
    inv19x.fill(-1); inv19y.fill(-1);
    for (int k = 0; k < 16; ++k) {
        inv19x[kH19X[k]] = k;
        inv19y[kH19Y[k]] = k;
    }

    auto get_pl_k = [&](int c, std::string& pl, int& k) {
        if (inv19x[c] >= 0) { pl = "19X"; k = inv19x[c]; }
        else if (inv19y[c] >= 0) { pl = "19Y"; k = inv19y[c]; }
        else { pl = "UNK"; k = -1; }
    };

    // Locate CSV directory
    fs::path csv_dir;
    if (sname.find("v8") != std::string::npos) {
        csv_dir = "outputs/hod_channel_correlations_20261017_202452";
    } else if (sname.find("v10") != std::string::npos) {
        csv_dir = "outputs/hod_channel_correlations_20261017_202805";
    }

    std::vector<std::vector<long long>> mat_same(64, std::vector<long long>(64, 0));
    std::vector<std::vector<long long>> mat_mixed(64, std::vector<long long>(64, 0));

    if (!csv_dir.empty() && fs::exists(csv_dir / "same_matrix_64x64.csv")) {
        // Read Same
        std::ifstream fsame(csv_dir / "same_matrix_64x64.csv");
        std::string line;
        std::getline(fsame, line); // header
        while (std::getline(fsame, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string item;
            std::getline(ss, item, ',');
            int a = std::stoi(item);
            for (int b = 0; b < 64 && std::getline(ss, item, ','); ++b) {
                mat_same[a][b] = std::stoll(item);
            }
        }
        // Read Mixed
        std::ifstream fmixed(csv_dir / "mixed_matrix_64x64.csv");
        std::getline(fmixed, line); // header
        while (std::getline(fmixed, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string item;
            std::getline(ss, item, ',');
            int a = std::stoi(item);
            for (int b = 0; b < 64 && std::getline(ss, item, ','); ++b) {
                mat_mixed[a][b] = std::stoll(item);
            }
        }

        for (int a = 32; a < 64; ++a) {
            for (int b = a + 1; b < 64; ++b) {
                TopPair tp;
                tp.a = a;
                tp.b = b;
                tp.grp_a = "G" + std::to_string(a / 8);
                tp.grp_b = "G" + std::to_string(b / 8);
                get_pl_k(a, tp.pl_a, tp.k_a);
                get_pl_k(b, tp.pl_b, tp.k_b);
                tp.Same = mat_same[a][b];
                tp.Mixed = mat_mixed[a][b];
                top_pairs.push_back(tp);
            }
        }

        std::sort(top_pairs.begin(), top_pairs.end(), [](const TopPair& p1, const TopPair& p2) {
            if (p1.Same != p2.Same) return p1.Same > p2.Same;
            if (p1.a != p2.a) return p1.a < p2.a;
            return p1.b < p2.b;
        });

        if (top_pairs.size() > 30) {
            top_pairs.resize(30);
        }
    }

    // Summary printer
    auto print_summary = [&](std::ostream& os) {
        os << "=== Hodoscope Hit Structure and Accidental Scan ===\n\n";
        os << "Input file: " << fname << "\n";
        os << "Total events: " << nEnt << "\n";
        os << "Output directory: " << outdir << "\n\n";

        os << "--- 1a. Hit Times (16 ns bins over S [4000, 4400] ns and C [6000, 6400] ns) ---\n";
        os << std::setw(12) << std::left << "t_center_S"
           << std::setw(10) << "18X_S"
           << std::setw(10) << "18X_C"
           << std::setw(10) << "18Y_S"
           << std::setw(10) << "18Y_C"
           << std::setw(10) << "G4_S"
           << std::setw(10) << "G5_S"
           << std::setw(10) << "G6_S"
           << std::setw(10) << "G7_S" << "\n";
        os << "----------------------------------------------------------------------------------------\n";
        for (int b = 0; b < 25; ++b) {
            double t_c = 4008.0 + 16.0 * b;
            os << std::setw(12) << std::left << std::fixed << std::setprecision(0) << t_c
               << std::setw(10) << hits_18X_S[b]
               << std::setw(10) << hits_18X_C[b]
               << std::setw(10) << hits_18Y_S[b]
               << std::setw(10) << hits_18Y_C[b]
               << std::setw(10) << hits_G4_S[b]
               << std::setw(10) << hits_G5_S[b]
               << std::setw(10) << hits_G6_S[b]
               << std::setw(10) << hits_G7_S[b] << "\n";
        }

        os << "\n--- 1b. hod2019 Multiplicity in Window S (events with >= 1 hit in ch 32-63) ---\n";
        os << "Hit channels \\ Groups hit:  G_1       G_2       G_3       G_4       Total\n";
        os << "--------------------------------------------------------------------------\n";
        const char* r_labels[5] = {"1", "2", "3", "4", "5+"};
        long long col_totals[4] = {};
        long long grand_total = 0;
        for (int r = 0; r < 5; ++r) {
            long long r_tot = 0;
            os << std::setw(28) << std::left << r_labels[r];
            for (int c = 0; c < 4; ++c) {
                os << std::setw(10) << mult_table[r][c];
                r_tot += mult_table[r][c];
                col_totals[c] += mult_table[r][c];
            }
            os << std::setw(10) << r_tot << "\n";
            grand_total += r_tot;
        }
        os << "--------------------------------------------------------------------------\n";
        os << std::setw(28) << std::left << "Total";
        for (int c = 0; c < 4; ++c) {
            os << std::setw(10) << col_totals[c];
        }
        os << std::setw(10) << grand_total << "\n";

        os << "\n--- 1c. hod2019 Pulse-Height Spectra per Group (Window S) ---\n";
        os << "Lowest wf value at imin over all hits [mV]:\n";
        os << "  G4 (ch 32-39): " << std::fixed << std::setprecision(2) << lowest_wf_imin[0] << " mV\n";
        os << "  G5 (ch 40-47): " << std::fixed << std::setprecision(2) << lowest_wf_imin[1] << " mV\n";
        os << "  G6 (ch 48-55): " << std::fixed << std::setprecision(2) << lowest_wf_imin[2] << " mV\n";
        os << "  G7 (ch 56-63): " << std::fixed << std::setprecision(2) << lowest_wf_imin[3] << " mV\n\n";

        os << std::setw(18) << std::left << "ph_bin [mV]"
           << std::setw(10) << "G4"
           << std::setw(10) << "G5"
           << std::setw(10) << "G6"
           << std::setw(10) << "G7" << "\n";
        os << "------------------------------------------------------------\n";
        for (int b = 0; b < 40; ++b) {
            double lo = 35.0 + 25.0 * b;
            double hi = lo + 25.0;
            char bname[32];
            std::snprintf(bname, sizeof(bname), "[%.0f, %.0f)", lo, hi);
            os << std::setw(18) << std::left << bname
               << std::setw(10) << ph_hist[0][b]
               << std::setw(10) << ph_hist[1][b]
               << std::setw(10) << ph_hist[2][b]
               << std::setw(10) << ph_hist[3][b] << "\n";
        }
        os << std::setw(18) << std::left << ">= 1035 (overflow)"
           << std::setw(10) << ph_hist[0][40]
           << std::setw(10) << ph_hist[1][40]
           << std::setw(10) << ph_hist[2][40]
           << std::setw(10) << ph_hist[3][40] << "\n";

        os << "\n--- 1d. Top 30 hod2019 Channel Pairs (a < b within ch 32..63) with Largest Same ---\n";
        os << std::setw(6) << std::left << "a"
           << std::setw(6) << "b"
           << std::setw(8) << "grp_a"
           << std::setw(8) << "grp_b"
           << std::setw(8) << "pl_a"
           << std::setw(6) << "k_a"
           << std::setw(8) << "pl_b"
           << std::setw(6) << "k_b"
           << std::setw(10) << "Same"
           << std::setw(10) << "Mixed" << "\n";
        os << "----------------------------------------------------------------------\n";
        for (const auto& tp : top_pairs) {
            os << std::setw(6) << std::left << tp.a
               << std::setw(6) << tp.b
               << std::setw(8) << tp.grp_a
               << std::setw(8) << tp.grp_b
               << std::setw(8) << tp.pl_a
               << std::setw(6) << tp.k_a
               << std::setw(8) << tp.pl_b
               << std::setw(6) << tp.k_b
               << std::setw(10) << tp.Same
               << std::setw(10) << tp.Mixed << "\n";
        }
    };

    print_summary(std::cout);
    {
        std::ofstream sf(outdir / "summary.txt");
        print_summary(sf);
    }
    std::cout << "\nResults written to: " << outdir << "\n";
}

