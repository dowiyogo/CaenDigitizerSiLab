#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

constexpr int kNChannels = 64;
constexpr int kNHod2018 = 32;
constexpr int kNHod2019 = 32;
constexpr int kNPlane = 16;
constexpr int kNBase = 50;
constexpr int kNCuts = 4;
constexpr int kLoose = 0;
constexpr int kMedium = 1;
constexpr int kStrict = 2;
constexpr int kVeryStrict = 3;

struct CutDef {
    std::string name;
    double nsigma = 0.0;
    double min_adc = 0.0;
};

const std::array<CutDef, kNCuts> kCuts = {{
    {"loose", 4.0, 5.0},
    {"medium", 5.0, 10.0},
    {"strict", 6.0, 15.0},
    {"very_strict", 5.0, 20.0},
}};

struct WaveformMetrics {
    bool valid = false;
    double baseline = 0.0;
    double sigma = 0.0;
    double vmin = 0.0;
    double vmax = 0.0;
    double amp_neg = 0.0;
    double amp_pos = 0.0;
    int sample_min = -1;
};

struct ChannelStats {
    bool branch_present = false;
    long long entries = 0;
    long long hits = 0;
    double sum_baseline = 0.0;
    double sum_baseline2 = 0.0;
    double sum_sigma = 0.0;
    double max_amp_neg = 0.0;
    double sum_amp_neg_hits = 0.0;
    double min_vmin = std::numeric_limits<double>::infinity();
    double max_vmax = -std::numeric_limits<double>::infinity();
    double sum_sample_min_hits = 0.0;
};

struct CutAccumulator {
    long long empty = 0;
    long long x_only = 0;
    long long y_only = 0;
    long long xy_candidate = 0;
    long long clean_xy = 0;
    long long near_clean_xy = 0;
    long long multi_small_xy = 0;
    long long multi_large_xy = 0;
    long double sum_nx = 0.0;
    long double sum_ny = 0.0;
    long double sum_nhits = 0.0;
    long long hod2019_total_hits = 0;
    int hod2019_channels_with_hits = 0;
    double hod2019_max_amp_neg = 0.0;
    std::array<long long, kNPlane> x_occ{};
    std::array<long long, kNPlane> y_occ{};
    std::array<std::array<long long, kNPlane>, kNPlane> xy_clean{};
    std::array<std::array<long long, kNPlane>, kNPlane> xy_near{};
    std::array<std::array<long long, kNPlane>, kNPlane> xy_all{};
};

struct CutSummary {
    int cut_index = 0;
    std::string cut_name;
    long long empty = 0;
    long long x_only = 0;
    long long y_only = 0;
    long long xy_candidate = 0;
    long long clean_xy = 0;
    long long near_clean_xy = 0;
    long long multi_small_xy = 0;
    long long multi_large_xy = 0;
    double frac_empty = 0.0;
    double frac_xy_candidate = 0.0;
    double frac_clean_xy = 0.0;
    double frac_near_clean_xy = 0.0;
    double frac_multi_large_xy = 0.0;
    double mean_nx = 0.0;
    double mean_ny = 0.0;
    double mean_nhits = 0.0;
    long long total_hits_x = 0;
    long long total_hits_y = 0;
    int dominant_x_channel = -1;
    int dominant_y_channel = -1;
    double dominant_x_fraction = 0.0;
    double dominant_y_fraction = 0.0;
    long long hod2019_total_hits = 0;
    double hod2019_max_amp_neg = 0.0;
    std::string hod2019_status;
    double quality_score = 0.0;
};

struct RunResult {
    int run_index = 0;
    std::string run_name;
    std::string filename;
    long long total_events = 0;
    bool ok = false;
    std::array<CutAccumulator, kNCuts> acc{};
    std::array<std::array<ChannelStats, kNChannels>, kNCuts> channel_stats{};
    std::array<CutSummary, kNCuts> summaries{};
    int best_cut_index = kMedium;
    double best_score = -1.0e99;
    std::string reason;
    bool recommended_for_analysis = false;
};

static std::string detector_name(int ch) {
    return (ch < 32) ? "hod2018" : "hod2019";
}

static std::string plane_name(int ch) {
    if (ch < 16) return "X";
    if (ch < 32) return "Y";
    return "HOD2019";
}

static double safe_fraction(long long num, long long den) {
    return den > 0 ? static_cast<double>(num) / static_cast<double>(den) : 0.0;
}

static std::string sanitize_name(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (unsigned char c : text) {
        if (std::isalnum(c) || c == '_' || c == '-') out.push_back(static_cast<char>(c));
        else out.push_back('_');
    }
    while (!out.empty() && out.back() == '_') out.pop_back();
    if (out.empty()) out = "run";
    return out;
}

static WaveformMetrics analyze_waveform(const std::vector<float>& v) {
    WaveformMetrics m;
    if (v.empty()) return m;

    m.valid = true;
    const int nbase = std::min<int>(kNBase, v.size());

    double sum = 0.0;
    for (int i = 0; i < nbase; ++i) sum += v[i];
    m.baseline = sum / static_cast<double>(nbase);

    if (nbase > 1) {
        double sum2 = 0.0;
        for (int i = 0; i < nbase; ++i) {
            const double d = v[i] - m.baseline;
            sum2 += d * d;
        }
        m.sigma = std::sqrt(sum2 / static_cast<double>(nbase - 1));
    }

    auto minmax = std::minmax_element(v.begin(), v.end());
    m.vmin = *minmax.first;
    m.vmax = *minmax.second;
    m.amp_neg = m.baseline - m.vmin;
    m.amp_pos = m.vmax - m.baseline;
    m.sample_min = static_cast<int>(std::distance(v.begin(), minmax.first));
    return m;
}

static std::array<bool, kNCuts> evaluate_hits(const std::vector<float>& v, const WaveformMetrics& m) {
    std::array<bool, kNCuts> hits{};
    if (!m.valid || v.size() < 2) return hits;

    std::array<double, kNCuts> delta{};
    std::array<bool, kNCuts> amp_ok{};
    std::array<bool, kNCuts> consec_ok{};

    for (int c = 0; c < kNCuts; ++c) {
        delta[c] = std::max(kCuts[c].nsigma * m.sigma, kCuts[c].min_adc);
        amp_ok[c] = (m.amp_neg > delta[c]);
    }

    for (size_t i = 1; i < v.size(); ++i) {
        for (int c = 0; c < kNCuts; ++c) {
            if (consec_ok[c] || !amp_ok[c]) continue;
            const double threshold = m.baseline - delta[c];
            if (v[i - 1] < threshold && v[i] < threshold) consec_ok[c] = true;
        }
    }

    for (int c = 0; c < kNCuts; ++c) hits[c] = amp_ok[c] && consec_ok[c];
    return hits;
}

static void update_channel_stats(ChannelStats& s, const WaveformMetrics& m, bool hit) {
    if (!m.valid) return;
    ++s.entries;
    s.sum_baseline += m.baseline;
    s.sum_baseline2 += m.baseline * m.baseline;
    s.sum_sigma += m.sigma;
    s.max_amp_neg = std::max(s.max_amp_neg, m.amp_neg);
    s.min_vmin = std::min(s.min_vmin, m.vmin);
    s.max_vmax = std::max(s.max_vmax, m.vmax);
    if (hit) {
        ++s.hits;
        s.sum_amp_neg_hits += m.amp_neg;
        s.sum_sample_min_hits += m.sample_min;
    }
}

static double mean_baseline(const ChannelStats& s) {
    return s.entries > 0 ? s.sum_baseline / static_cast<double>(s.entries) : 0.0;
}

static double mean_sigma(const ChannelStats& s) {
    return s.entries > 0 ? s.sum_sigma / static_cast<double>(s.entries) : 0.0;
}

static double rms_baseline(const ChannelStats& s) {
    if (s.entries <= 0) return 0.0;
    const double mean = s.sum_baseline / static_cast<double>(s.entries);
    const double mean2 = s.sum_baseline2 / static_cast<double>(s.entries);
    return std::sqrt(std::max(0.0, mean2 - mean * mean));
}

static double mean_amp_hits(const ChannelStats& s) {
    return s.hits > 0 ? s.sum_amp_neg_hits / static_cast<double>(s.hits) : 0.0;
}

static double mean_sample_min_hits(const ChannelStats& s) {
    return s.hits > 0 ? s.sum_sample_min_hits / static_cast<double>(s.hits) : 0.0;
}

static std::string channel_status(const RunResult& run, int cut_index, int ch) {
    const ChannelStats& s = run.channel_stats[cut_index][ch];
    if (!s.branch_present || s.entries == 0) return "missing_branch";

    const double frac = safe_fraction(s.hits, s.entries);
    if (ch < 32) {
        if (frac > 0.5) return "suspicious_noise";
        if (s.hits == 0) return "dead";
        if (frac <= 0.001) return "weak";
        return "active";
    }

    const long long strict_hits = run.channel_stats[kStrict][ch].hits;
    const double strict_max_amp = run.channel_stats[kStrict][ch].max_amp_neg;
    if (strict_hits == 0 && strict_max_amp < 5.0) return "hod2019_no_signal";
    if (s.hits == 0) return "no_hits_this_cut";
    if (strict_hits == 0) return "noise_or_outliers";
    return "possible_signal_or_outlier";
}

static std::string hod2019_status_for_cut(const RunResult& run, int cut_index) {
    const CutAccumulator& current = run.acc[cut_index];
    const CutAccumulator& strict = run.acc[kStrict];
    const double strict_hit_frac = safe_fraction(strict.hod2019_total_hits, run.total_events * kNHod2019);

    if (strict.hod2019_total_hits <= 2 && strict.hod2019_max_amp_neg < 5.0) {
        return "no_signal";
    }
    if (current.hod2019_total_hits > 0 && strict.hod2019_total_hits == 0) {
        return "noise_or_outliers";
    }
    if (strict.hod2019_total_hits < 10 || strict_hit_frac < 1.0e-4) {
        return "noise_or_outliers";
    }
    if (strict.hod2019_channels_with_hits >= 4 && strict.hod2019_max_amp_neg >= 15.0) {
        return "possible_signal";
    }
    return "noise_or_outliers";
}

static void configure_hist_1d(TH1* h, const std::string& x_title, const std::string& y_title) {
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(x_title.c_str());
    h->GetYaxis()->SetTitle(y_title.c_str());
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->SetLineWidth(2);
}

static void configure_hist_2d(TH2* h, const std::string& x_title, const std::string& y_title,
                              const std::string& z_title) {
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(x_title.c_str());
    h->GetYaxis()->SetTitle(y_title.c_str());
    h->GetZaxis()->SetTitle(z_title.c_str());
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetZaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
}

static void save_hist1_png(TH1* h, const fs::path& path, const std::string& option = "HIST") {
    TCanvas c(("c_" + std::string(h->GetName())).c_str(), "", 1100, 760);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.13);
    h->Draw(option.c_str());
    c.SaveAs(path.string().c_str());
}

static void save_hist2_png(TH2* h, const fs::path& path, const std::string& option = "COLZ") {
    TCanvas c(("c_" + std::string(h->GetName())).c_str(), "", 1000, 840);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.16);
    c.SetBottomMargin(0.12);
    h->Draw(option.c_str());
    c.SaveAs(path.string().c_str());
}

static int max_channel_index(const std::array<long long, kNPlane>& counts) {
    return static_cast<int>(std::distance(counts.begin(), std::max_element(counts.begin(), counts.end())));
}

static long long sum_counts(const std::array<long long, kNPlane>& counts) {
    long long total = 0;
    for (long long v : counts) total += v;
    return total;
}

static double quality_score(const CutSummary& s) {
    double penalty_hod2019 = 0.0;
    if (s.hod2019_status == "possible_signal") penalty_hod2019 = 0.15;
    if (s.hod2019_status == "noise_or_outliers") penalty_hod2019 = 0.05;
    if (s.hod2019_total_hits > 0 && s.frac_multi_large_xy > 0.2) penalty_hod2019 += 0.05;

    return 3.0 * s.frac_near_clean_xy
         + 1.0 * s.frac_clean_xy
         - 2.0 * s.frac_multi_large_xy
         - std::abs(s.dominant_x_fraction - 0.25)
         - std::abs(s.dominant_y_fraction - 0.25)
         - penalty_hod2019;
}

static std::string reason_for_summary(const CutSummary& s) {
    if (s.dominant_x_fraction > 0.6 || s.dominant_y_fraction > 0.6) return "dominated_by_single_channel";
    if (s.frac_multi_large_xy > 0.3) return "too_noisy";
    if (s.frac_xy_candidate < 0.001 || s.frac_near_clean_xy < 0.0005) return "too_few_hits";
    if (s.frac_near_clean_xy > 0.0 && s.frac_multi_large_xy < 0.2) return "reasonable";
    return "mixed_quality";
}

static void finalize_run_summaries(RunResult& run) {
    for (int c = 0; c < kNCuts; ++c) {
        CutAccumulator& a = run.acc[c];

        for (int ch = 32; ch < 64; ++ch) {
            if (run.channel_stats[c][ch].hits > 0) ++a.hod2019_channels_with_hits;
            a.hod2019_max_amp_neg = std::max(a.hod2019_max_amp_neg, run.channel_stats[c][ch].max_amp_neg);
        }

        CutSummary s;
        s.cut_index = c;
        s.cut_name = kCuts[c].name;
        s.empty = a.empty;
        s.x_only = a.x_only;
        s.y_only = a.y_only;
        s.xy_candidate = a.xy_candidate;
        s.clean_xy = a.clean_xy;
        s.near_clean_xy = a.near_clean_xy;
        s.multi_small_xy = a.multi_small_xy;
        s.multi_large_xy = a.multi_large_xy;
        s.frac_empty = safe_fraction(a.empty, run.total_events);
        s.frac_xy_candidate = safe_fraction(a.xy_candidate, run.total_events);
        s.frac_clean_xy = safe_fraction(a.clean_xy, run.total_events);
        s.frac_near_clean_xy = safe_fraction(a.near_clean_xy, run.total_events);
        s.frac_multi_large_xy = safe_fraction(a.multi_large_xy, run.total_events);
        s.mean_nx = run.total_events > 0 ? static_cast<double>(a.sum_nx / run.total_events) : 0.0;
        s.mean_ny = run.total_events > 0 ? static_cast<double>(a.sum_ny / run.total_events) : 0.0;
        s.mean_nhits = run.total_events > 0 ? static_cast<double>(a.sum_nhits / run.total_events) : 0.0;
        s.total_hits_x = sum_counts(a.x_occ);
        s.total_hits_y = sum_counts(a.y_occ);
        s.dominant_x_channel = max_channel_index(a.x_occ);
        s.dominant_y_channel = max_channel_index(a.y_occ);
        s.dominant_x_fraction = safe_fraction(a.x_occ[s.dominant_x_channel], s.total_hits_x);
        s.dominant_y_fraction = safe_fraction(a.y_occ[s.dominant_y_channel], s.total_hits_y);
        s.hod2019_total_hits = a.hod2019_total_hits;
        s.hod2019_max_amp_neg = a.hod2019_max_amp_neg;
        s.hod2019_status = hod2019_status_for_cut(run, c);
        s.quality_score = quality_score(s);
        run.summaries[c] = s;

        if (s.quality_score > run.best_score) {
            run.best_score = s.quality_score;
            run.best_cut_index = c;
        }
    }
    run.reason = reason_for_summary(run.summaries[run.best_cut_index]);
}

static RunResult process_run(int run_index, const fs::path& filename) {
    RunResult run;
    run.run_index = run_index;
    run.filename = filename.string();
    run.run_name = sanitize_name(filename.stem().string());

    std::cout << "\n=== Run " << run.run_index << ": " << run.filename << " ===\n";

    std::unique_ptr<TFile> file(TFile::Open(filename.string().c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: no pude abrir " << filename << "\n";
        return run;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("data_dig"));
    if (!tree) {
        std::cerr << "ERROR: no encontre TTree data_dig en " << filename << "\n";
        return run;
    }

    run.total_events = tree->GetEntries();
    std::cout << "Entradas: " << run.total_events << "\n";

    tree->SetBranchStatus("*", 0);
    tree->SetCacheSize(64LL * 1024LL * 1024LL);

    std::array<std::vector<float>*, kNChannels> channels{};
    std::array<bool, kNChannels> branch_present{};
    for (int ch = 0; ch < kNChannels; ++ch) {
        const std::string branch_name = "Ch" + std::to_string(ch);
        if (!tree->GetBranch(branch_name.c_str())) {
            std::cerr << "ADVERTENCIA: falta rama " << branch_name << " en " << filename << "\n";
            continue;
        }
        branch_present[ch] = true;
        tree->SetBranchStatus(branch_name.c_str(), 1);
        if (tree->SetBranchAddress(branch_name.c_str(), &channels[ch]) < 0) {
            std::cerr << "ADVERTENCIA: no pude conectar " << branch_name << " en " << filename << "\n";
            branch_present[ch] = false;
        }
        for (int c = 0; c < kNCuts; ++c) run.channel_stats[c][ch].branch_present = branch_present[ch];
    }

    if (tree->GetBranch("event")) tree->SetBranchStatus("event", 1);
    if (tree->GetBranch("timestamp")) tree->SetBranchStatus("timestamp", 1);
    tree->AddBranchToCache("*", kTRUE);

    for (long long entry = 0; entry < run.total_events; ++entry) {
        tree->GetEntry(entry);
        if (entry % 1000 == 0) {
            std::cout << "Run " << run.run_index << " procesando evento "
                      << entry << " / " << run.total_events << "\r" << std::flush;
        }

        std::array<std::vector<int>, kNCuts> x_hits;
        std::array<std::vector<int>, kNCuts> y_hits;
        for (int c = 0; c < kNCuts; ++c) {
            x_hits[c].reserve(kNPlane);
            y_hits[c].reserve(kNPlane);
        }

        for (int ch = 0; ch < kNChannels; ++ch) {
            if (!branch_present[ch] || channels[ch] == nullptr) continue;
            const std::vector<float>& waveform = *channels[ch];
            const WaveformMetrics m = analyze_waveform(waveform);
            if (!m.valid) continue;
            const std::array<bool, kNCuts> hits = evaluate_hits(waveform, m);

            for (int c = 0; c < kNCuts; ++c) {
                update_channel_stats(run.channel_stats[c][ch], m, hits[c]);
                if (!hits[c]) continue;

                if (ch < 16) {
                    x_hits[c].push_back(ch);
                    ++run.acc[c].x_occ[ch];
                } else if (ch < 32) {
                    y_hits[c].push_back(ch - 16);
                    ++run.acc[c].y_occ[ch - 16];
                } else {
                    ++run.acc[c].hod2019_total_hits;
                }
            }
        }

        for (int c = 0; c < kNCuts; ++c) {
            CutAccumulator& a = run.acc[c];
            const int nX = static_cast<int>(x_hits[c].size());
            const int nY = static_cast<int>(y_hits[c].size());
            const int nHits = nX + nY;
            a.sum_nx += nX;
            a.sum_ny += nY;
            a.sum_nhits += nHits;

            if (nX == 0 && nY == 0) ++a.empty;
            if (nX > 0 && nY == 0) ++a.x_only;
            if (nX == 0 && nY > 0) ++a.y_only;
            if (nX > 0 && nY > 0) {
                ++a.xy_candidate;
                for (int x : x_hits[c]) {
                    for (int y : y_hits[c]) {
                        ++a.xy_all[x][y];
                    }
                }
            }
            if (nX == 1 && nY == 1) {
                ++a.clean_xy;
                ++a.xy_clean[x_hits[c][0]][y_hits[c][0]];
            }
            if (nX >= 1 && nX <= 2 && nY >= 1 && nY <= 2) {
                ++a.near_clean_xy;
                for (int x : x_hits[c]) {
                    for (int y : y_hits[c]) {
                        ++a.xy_near[x][y];
                    }
                }
            }
            if (nX >= 1 && nX <= 3 && nY >= 1 && nY <= 3) ++a.multi_small_xy;
            if (nX > 3 || nY > 3) ++a.multi_large_xy;
        }
    }

    std::cout << "\nRun " << run.run_index << " terminado.\n";
    finalize_run_summaries(run);
    run.ok = true;
    return run;
}

static void write_csvs(const fs::path& csvdir, const std::vector<RunResult>& runs) {
    std::ofstream run_csv(csvdir / "run_summary.csv");
    run_csv << "run_index,run_name,filename,total_events,cut_name,"
            << "empty,X_only,Y_only,XY_candidate,clean_XY,near_clean_XY,multi_small_XY,multi_large_XY,"
            << "frac_empty,frac_XY_candidate,frac_clean_XY,frac_near_clean_XY,frac_multi_large_XY,"
            << "mean_nX,mean_nY,mean_nHits,hod2018_total_hits_X,hod2018_total_hits_Y,"
            << "dominant_X_channel,dominant_X_fraction,dominant_Y_channel,dominant_Y_fraction,"
            << "hod2019_total_hits,hod2019_max_amp_neg,hod2019_status,quality_score\n";
    run_csv << std::setprecision(10);

    std::ofstream channel_csv(csvdir / "channel_summary_all_runs.csv");
    channel_csv << "run_index,run_name,cut_name,channel,detector,plane,entries,hits,frac_hits,"
                << "mean_baseline,mean_sigma,rms_baseline,max_amp_neg,mean_amp_neg_hits,"
                << "min_vmin,max_vmax,mean_sample_min_hits,status\n";
    channel_csv << std::setprecision(10);

    std::ofstream class_csv(csvdir / "event_class_counts_all_runs.csv");
    class_csv << "run_index,run_name,cut_name,class,count,fraction\n";
    class_csv << std::setprecision(10);

    std::ofstream xy_csv(csvdir / "xy_map_summary_all_runs.csv");
    xy_csv << "run_index,run_name,cut_name,map_type,x_channel,y_channel,count\n";

    for (const RunResult& run : runs) {
        if (!run.ok) continue;
        for (int c = 0; c < kNCuts; ++c) {
            const CutSummary& s = run.summaries[c];
            run_csv << run.run_index << ',' << run.run_name << ',' << run.filename << ','
                    << run.total_events << ',' << s.cut_name << ','
                    << s.empty << ',' << s.x_only << ',' << s.y_only << ','
                    << s.xy_candidate << ',' << s.clean_xy << ',' << s.near_clean_xy << ','
                    << s.multi_small_xy << ',' << s.multi_large_xy << ','
                    << s.frac_empty << ',' << s.frac_xy_candidate << ',' << s.frac_clean_xy << ','
                    << s.frac_near_clean_xy << ',' << s.frac_multi_large_xy << ','
                    << s.mean_nx << ',' << s.mean_ny << ',' << s.mean_nhits << ','
                    << s.total_hits_x << ',' << s.total_hits_y << ','
                    << s.dominant_x_channel << ',' << s.dominant_x_fraction << ','
                    << s.dominant_y_channel << ',' << s.dominant_y_fraction << ','
                    << s.hod2019_total_hits << ',' << s.hod2019_max_amp_neg << ','
                    << s.hod2019_status << ',' << s.quality_score << '\n';

            const std::array<std::pair<std::string, long long>, 8> class_rows = {{
                {"empty", s.empty},
                {"X_only", s.x_only},
                {"Y_only", s.y_only},
                {"XY_candidate", s.xy_candidate},
                {"clean_XY", s.clean_xy},
                {"near_clean_XY", s.near_clean_xy},
                {"multi_small_XY", s.multi_small_xy},
                {"multi_large_XY", s.multi_large_xy},
            }};
            for (const auto& row : class_rows) {
                class_csv << run.run_index << ',' << run.run_name << ',' << s.cut_name << ','
                          << row.first << ',' << row.second << ','
                          << safe_fraction(row.second, run.total_events) << '\n';
            }

            for (int ch = 0; ch < kNChannels; ++ch) {
                const ChannelStats& cs = run.channel_stats[c][ch];
                const double min_vmin = std::isfinite(cs.min_vmin) ? cs.min_vmin : 0.0;
                const double max_vmax = std::isfinite(cs.max_vmax) ? cs.max_vmax : 0.0;
                channel_csv << run.run_index << ',' << run.run_name << ',' << s.cut_name << ','
                            << ch << ',' << detector_name(ch) << ',' << plane_name(ch) << ','
                            << cs.entries << ',' << cs.hits << ',' << safe_fraction(cs.hits, cs.entries) << ','
                            << mean_baseline(cs) << ',' << mean_sigma(cs) << ',' << rms_baseline(cs) << ','
                            << cs.max_amp_neg << ',' << mean_amp_hits(cs) << ','
                            << min_vmin << ',' << max_vmax << ',' << mean_sample_min_hits(cs) << ','
                            << channel_status(run, c, ch) << '\n';
            }

            for (int x = 0; x < kNPlane; ++x) {
                for (int y = 0; y < kNPlane; ++y) {
                    xy_csv << run.run_index << ',' << run.run_name << ',' << s.cut_name
                           << ",clean_XY," << x << ',' << y << ',' << run.acc[c].xy_clean[x][y] << '\n';
                    xy_csv << run.run_index << ',' << run.run_name << ',' << s.cut_name
                           << ",near_clean_XY," << x << ',' << y << ',' << run.acc[c].xy_near[x][y] << '\n';
                    xy_csv << run.run_index << ',' << run.run_name << ',' << s.cut_name
                           << ",all_combinations," << x << ',' << y << ',' << run.acc[c].xy_all[x][y] << '\n';
                }
            }
        }
    }
}

static void write_recommendation_csv(const fs::path& csvdir, std::vector<RunResult>& runs) {
    int best_run_pos = -1;
    double best_score = -1.0e99;
    for (size_t i = 0; i < runs.size(); ++i) {
        if (!runs[i].ok) continue;
        if (runs[i].best_score > best_score) {
            best_score = runs[i].best_score;
            best_run_pos = static_cast<int>(i);
        }
    }
    if (best_run_pos >= 0) runs[best_run_pos].recommended_for_analysis = true;

    std::ofstream rec(csvdir / "recommendation_summary.csv");
    rec << "run_index,run_name,best_cut,quality_score,reason,recommended_for_analysis\n";
    rec << std::setprecision(10);
    for (RunResult& run : runs) {
        if (!run.ok) continue;
        std::string reason = run.reason;
        if (run.recommended_for_analysis) reason = "best_available";
        rec << run.run_index << ',' << run.run_name << ','
            << kCuts[run.best_cut_index].name << ','
            << run.best_score << ',' << reason << ','
            << (run.recommended_for_analysis ? "true" : "false") << '\n';
    }
}

static const CutSummary* find_summary(const RunResult& run, const std::string& cut_name) {
    for (const CutSummary& s : run.summaries) {
        if (s.cut_name == cut_name) return &s;
    }
    return nullptr;
}

static void make_comparison_graph(const std::vector<RunResult>& runs, const fs::path& png_path,
                                  TFile& root_out, const std::string& name,
                                  const std::string& title, const std::string& y_title,
                                  double (*value_fn)(const CutSummary&)) {
    TCanvas canvas(("c_" + name).c_str(), "", 1100, 760);
    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.04);
    canvas.SetBottomMargin(0.13);
    TMultiGraph mg;
    mg.SetTitle((title + ";Run index;" + y_title).c_str());
    const int colors[kNCuts] = {kAzure + 2, kOrange + 7, kGreen + 3, kRed + 1};
    TLegend legend(0.70, 0.70, 0.92, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    std::vector<std::unique_ptr<TGraph>> graphs;
    graphs.reserve(kNCuts);

    for (int c = 0; c < kNCuts; ++c) {
        auto graph = std::make_unique<TGraph>();
        graph->SetName((name + "_" + kCuts[c].name).c_str());
        graph->SetLineColor(colors[c]);
        graph->SetMarkerColor(colors[c]);
        graph->SetMarkerStyle(20 + c);
        graph->SetLineWidth(2);
        int p = 0;
        for (const RunResult& run : runs) {
            if (!run.ok) continue;
            graph->SetPoint(p++, run.run_index, value_fn(run.summaries[c]));
        }
        mg.Add(graph.get(), "LP");
        legend.AddEntry(graph.get(), kCuts[c].name.c_str(), "lp");
        graphs.push_back(std::move(graph));
    }

    mg.Draw("A");
    mg.GetXaxis()->SetNdivisions(static_cast<int>(runs.size()), false);
    mg.GetYaxis()->SetTitleOffset(1.2);
    legend.Draw();
    canvas.SaveAs(png_path.string().c_str());
    root_out.cd();
    for (auto& g : graphs) g->Write();
}

static double get_frac_near(const CutSummary& s) { return s.frac_near_clean_xy; }
static double get_frac_multi_large(const CutSummary& s) { return s.frac_multi_large_xy; }
static double get_quality(const CutSummary& s) { return s.quality_score; }
static double get_hod2019_hits(const CutSummary& s) { return static_cast<double>(s.hod2019_total_hits); }

static void write_run_cut_hists(const RunResult& run, int cut_index, const fs::path& pngdir, TFile& root_out) {
    const std::string suffix = (cut_index == kStrict) ? (run.run_name + "_strict") : run.run_name;
    const std::string hprefix = "h_run" + std::to_string(run.run_index) + "_" + kCuts[cut_index].name + "_";
    const CutAccumulator& a = run.acc[cut_index];

    auto hx = std::make_unique<TH1D>((hprefix + "x_occupancy").c_str(),
                                    ("hod2018 X occupancy " + run.run_name + " " + kCuts[cut_index].name + ";X channel;Hits").c_str(),
                                    16, -0.5, 15.5);
    auto hy = std::make_unique<TH1D>((hprefix + "y_occupancy").c_str(),
                                    ("hod2018 Y occupancy " + run.run_name + " " + kCuts[cut_index].name + ";Y channel;Hits").c_str(),
                                    16, -0.5, 15.5);
    auto hnear = std::make_unique<TH2D>((hprefix + "xy_near_clean").c_str(),
                                       ("hod2018 near_clean map " + run.run_name + " " + kCuts[cut_index].name + ";X channel;Y channel;Counts").c_str(),
                                       16, -0.5, 15.5, 16, -0.5, 15.5);
    configure_hist_1d(hx.get(), "X channel", "Hits");
    configure_hist_1d(hy.get(), "Y channel", "Hits");
    configure_hist_2d(hnear.get(), "X channel", "Y channel", "Counts");
    hx->SetFillColor(kAzure - 9);
    hx->SetLineColor(kAzure + 2);
    hy->SetFillColor(kOrange - 3);
    hy->SetLineColor(kOrange + 7);

    for (int i = 0; i < 16; ++i) {
        hx->SetBinContent(i + 1, a.x_occ[i]);
        hy->SetBinContent(i + 1, a.y_occ[i]);
    }
    for (int x = 0; x < 16; ++x) {
        for (int y = 0; y < 16; ++y) {
            hnear->SetBinContent(x + 1, y + 1, a.xy_near[x][y]);
        }
    }

    save_hist1_png(hx.get(), pngdir / ("hod2018_x_occupancy_" + suffix + ".png"));
    save_hist1_png(hy.get(), pngdir / ("hod2018_y_occupancy_" + suffix + ".png"));
    save_hist2_png(hnear.get(), pngdir / ("hod2018_xy_map_near_clean_" + suffix + ".png"));
    root_out.cd();
    hx->Write();
    hy->Write();
    hnear->Write();

    if (cut_index == kMedium) {
        auto hall = std::make_unique<TH2D>((hprefix + "xy_all_combinations").c_str(),
                                          ("hod2018 all combinations " + run.run_name + " medium;X channel;Y channel;Counts").c_str(),
                                          16, -0.5, 15.5, 16, -0.5, 15.5);
        configure_hist_2d(hall.get(), "X channel", "Y channel", "Counts");
        for (int x = 0; x < 16; ++x) {
            for (int y = 0; y < 16; ++y) {
                hall->SetBinContent(x + 1, y + 1, a.xy_all[x][y]);
            }
        }
        save_hist2_png(hall.get(), pngdir / ("hod2018_xy_map_all_combinations_" + run.run_name + ".png"));
        root_out.cd();
        hall->Write();
    }
}

static void make_event_classes_medium(const std::vector<RunResult>& runs, const fs::path& pngdir, TFile& root_out) {
    auto h_empty = std::make_unique<TH1D>("h_compare_event_empty_medium", "Event class fractions, medium cut;Run index;Fraction", runs.size(), 0.5, runs.size() + 0.5);
    auto h_clean = std::make_unique<TH1D>("h_compare_event_clean_medium", "Event class fractions, medium cut;Run index;Fraction", runs.size(), 0.5, runs.size() + 0.5);
    auto h_near = std::make_unique<TH1D>("h_compare_event_near_medium", "Event class fractions, medium cut;Run index;Fraction", runs.size(), 0.5, runs.size() + 0.5);
    auto h_large = std::make_unique<TH1D>("h_compare_event_large_medium", "Event class fractions, medium cut;Run index;Fraction", runs.size(), 0.5, runs.size() + 0.5);
    std::array<TH1D*, 4> hists = {h_empty.get(), h_clean.get(), h_near.get(), h_large.get()};
    const std::array<std::string, 4> labels = {"empty", "clean_XY", "near_clean_XY", "multi_large_XY"};
    const int colors[4] = {kGray + 1, kAzure + 2, kGreen + 3, kRed + 1};

    double ymax = 0.0;
    for (size_t i = 0; i < runs.size(); ++i) {
        const RunResult& run = runs[i];
        const CutSummary& s = run.summaries[kMedium];
        const std::array<double, 4> vals = {s.frac_empty, s.frac_clean_xy, s.frac_near_clean_xy, s.frac_multi_large_xy};
        for (int j = 0; j < 4; ++j) {
            hists[j]->SetBinContent(i + 1, vals[j]);
            hists[j]->GetXaxis()->SetBinLabel(i + 1, std::to_string(run.run_index).c_str());
            ymax = std::max(ymax, vals[j]);
        }
    }

    for (int j = 0; j < 4; ++j) {
        configure_hist_1d(hists[j], "Run index", "Fraction");
        hists[j]->SetFillColor(colors[j]);
        hists[j]->SetLineColor(colors[j]);
        hists[j]->SetBarWidth(0.18);
        hists[j]->SetBarOffset(0.12 + 0.19 * j);
        hists[j]->SetMaximum(std::max(0.05, ymax * 1.25));
    }

    TCanvas canvas("c_compare_event_classes_by_run_medium", "", 1100, 760);
    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.04);
    canvas.SetBottomMargin(0.13);
    h_empty->Draw("BAR");
    h_clean->Draw("BAR SAME");
    h_near->Draw("BAR SAME");
    h_large->Draw("BAR SAME");
    TLegend legend(0.68, 0.70, 0.92, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    for (int j = 0; j < 4; ++j) legend.AddEntry(hists[j], labels[j].c_str(), "f");
    legend.Draw();
    canvas.SaveAs((pngdir / "compare_event_classes_by_run_medium.png").string().c_str());

    root_out.cd();
    for (TH1D* h : hists) h->Write();
}

static void make_mean_multiplicity_plot(const std::vector<RunResult>& runs, const fs::path& pngdir, TFile& root_out) {
    const int colors[kNCuts] = {kAzure + 2, kOrange + 7, kGreen + 3, kRed + 1};
    const std::array<std::string, 3> panels = {"mean_nX", "mean_nY", "mean_nHits"};
    TCanvas canvas("c_compare_mean_multiplicity_by_run", "", 1100, 1200);
    canvas.Divide(1, 3);
    std::vector<std::unique_ptr<TGraph>> graphs;

    for (int panel = 0; panel < 3; ++panel) {
        canvas.cd(panel + 1);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.13);
        TMultiGraph* mg = new TMultiGraph();
        mg->SetName(("mg_" + panels[panel]).c_str());
        mg->SetTitle(("Mean multiplicity;Run index;" + panels[panel]).c_str());
        TLegend* legend = new TLegend(0.72, 0.68, 0.92, 0.90);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (int c = 0; c < kNCuts; ++c) {
            auto graph = std::make_unique<TGraph>();
            graph->SetName(("g_" + panels[panel] + "_" + kCuts[c].name).c_str());
            graph->SetLineColor(colors[c]);
            graph->SetMarkerColor(colors[c]);
            graph->SetMarkerStyle(20 + c);
            graph->SetLineWidth(2);
            int p = 0;
            for (const RunResult& run : runs) {
                const CutSummary& s = run.summaries[c];
                const double value = (panel == 0) ? s.mean_nx : (panel == 1 ? s.mean_ny : s.mean_nhits);
                graph->SetPoint(p++, run.run_index, value);
            }
            mg->Add(graph.get(), "LP");
            legend->AddEntry(graph.get(), kCuts[c].name.c_str(), "lp");
            graphs.push_back(std::move(graph));
        }
        mg->Draw("A");
        mg->GetXaxis()->SetNdivisions(static_cast<int>(runs.size()), false);
        legend->Draw();
        root_out.cd();
        mg->Write();
    }
    canvas.SaveAs((pngdir / "compare_mean_multiplicity_by_run.png").string().c_str());
    root_out.cd();
    for (auto& g : graphs) g->Write();
}

static void make_dominant_fraction_plot(const std::vector<RunResult>& runs, const fs::path& pngdir, TFile& root_out) {
    TCanvas canvas("c_compare_dominant_channel_fraction_by_run", "", 1100, 760);
    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.04);
    canvas.SetBottomMargin(0.13);
    TMultiGraph mg;
    mg.SetTitle("Dominant channel fraction, medium cut;Run index;Dominant fraction");
    auto gx = std::make_unique<TGraph>();
    auto gy = std::make_unique<TGraph>();
    gx->SetName("g_dominant_x_fraction_medium");
    gy->SetName("g_dominant_y_fraction_medium");
    gx->SetLineColor(kAzure + 2);
    gx->SetMarkerColor(kAzure + 2);
    gx->SetMarkerStyle(20);
    gx->SetLineWidth(2);
    gy->SetLineColor(kOrange + 7);
    gy->SetMarkerColor(kOrange + 7);
    gy->SetMarkerStyle(21);
    gy->SetLineWidth(2);
    int p = 0;
    for (const RunResult& run : runs) {
        const CutSummary& s = run.summaries[kMedium];
        gx->SetPoint(p, run.run_index, s.dominant_x_fraction);
        gy->SetPoint(p, run.run_index, s.dominant_y_fraction);
        ++p;
    }
    mg.Add(gx.get(), "LP");
    mg.Add(gy.get(), "LP");
    mg.Draw("A");
    mg.GetXaxis()->SetNdivisions(static_cast<int>(runs.size()), false);
    TLegend legend(0.70, 0.78, 0.92, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(gx.get(), "X medium", "lp");
    legend.AddEntry(gy.get(), "Y medium", "lp");
    legend.Draw();
    canvas.SaveAs((pngdir / "compare_dominant_channel_fraction_by_run.png").string().c_str());
    root_out.cd();
    gx->Write();
    gy->Write();
}

static void make_hod2019_max_amp_plot(const std::vector<RunResult>& runs, const fs::path& pngdir, TFile& root_out) {
    auto h = std::make_unique<TH1D>("h_compare_hod2019_max_amp_neg_by_run",
                                   "hod2019 max amp_{neg} by run;Run index;max amp_{neg} [ADC]",
                                   runs.size(), 0.5, runs.size() + 0.5);
    configure_hist_1d(h.get(), "Run index", "max amp_{neg} [ADC]");
    h->SetFillColor(kRed - 9);
    h->SetLineColor(kRed + 1);
    for (size_t i = 0; i < runs.size(); ++i) {
        h->SetBinContent(i + 1, runs[i].summaries[kStrict].hod2019_max_amp_neg);
        h->GetXaxis()->SetBinLabel(i + 1, std::to_string(runs[i].run_index).c_str());
    }
    save_hist1_png(h.get(), pngdir / "compare_hod2019_max_amp_neg_by_run.png");
    root_out.cd();
    h->Write();
}

static void make_all_plots(const fs::path& outdir, const std::vector<RunResult>& runs) {
    const fs::path pngdir = outdir / "png";
    TFile root_out((outdir / "compare_hv_runs_histograms.root").string().c_str(), "RECREATE");
    if (root_out.IsZombie()) {
        std::cerr << "ADVERTENCIA: no pude crear ROOT de histogramas\n";
        return;
    }

    make_event_classes_medium(runs, pngdir, root_out);
    make_comparison_graph(runs, pngdir / "compare_near_clean_fraction_by_run.png", root_out,
                          "g_compare_near_clean_fraction", "near_clean_XY fraction by run", "Fraction", get_frac_near);
    make_comparison_graph(runs, pngdir / "compare_multi_large_fraction_by_run.png", root_out,
                          "g_compare_multi_large_fraction", "multi_large_XY fraction by run", "Fraction", get_frac_multi_large);
    make_mean_multiplicity_plot(runs, pngdir, root_out);
    make_dominant_fraction_plot(runs, pngdir, root_out);
    make_comparison_graph(runs, pngdir / "compare_quality_score_by_run.png", root_out,
                          "g_compare_quality_score", "Quality score by run", "quality_score", get_quality);
    make_comparison_graph(runs, pngdir / "compare_hod2019_total_hits_by_run.png", root_out,
                          "g_compare_hod2019_total_hits", "hod2019 total hits by run", "Hits", get_hod2019_hits);
    make_hod2019_max_amp_plot(runs, pngdir, root_out);

    for (const RunResult& run : runs) {
        if (!run.ok) continue;
        write_run_cut_hists(run, kMedium, pngdir, root_out);
        write_run_cut_hists(run, kStrict, pngdir, root_out);
    }

    root_out.Close();
}

static void write_readme(const fs::path& outdir) {
    std::ofstream readme(outdir / "README.md");
    readme << "# Comparacion de runs HV\n\n";
    readme << "Analisis comparativo de runs CAEN VX1740 con hod2018 en Ch0-Ch31 y diagnostico de hod2019 en Ch32-Ch63.\n\n";
    readme << "## Comandos\n\n";
    readme << "Compilacion:\n\n";
    readme << "```bash\n";
    readme << "g++ -O2 -std=c++17 compare_hv_runs.cpp $(root-config --cflags --libs) -o compare_hv_runs\n";
    readme << "```\n\n";
    readme << "Ejecucion:\n\n";
    readme << "```bash\n";
    readme << "./compare_hv_runs outputs/compare_hv_runs \\\n";
    readme << "  lectura_vx1740_ext_nim_1M.root \\\n";
    readme << "  lectura_vx1740_ext_nim_1M_v2.root \\\n";
    readme << "  lectura_vx1740_ext_nim_1M_v3.root \\\n";
    readme << "  lectura_vx1740_ext_nim_1M_v4.root\n";
    readme << "```\n\n";
    readme << "Tambien se puede correr `./run_compare_hv_runs.sh`.\n\n";
    readme << "## Archivos\n\n";
    readme << "- `csv/run_summary.csv`: resumen por run y corte.\n";
    readme << "- `csv/channel_summary_all_runs.csv`: resumen por canal, run y corte.\n";
    readme << "- `csv/event_class_counts_all_runs.csv`: conteos de clases. Las clases `clean_XY`, `near_clean_XY` y `multi_small_XY` son inclusivas/subconjuntos.\n";
    readme << "- `csv/xy_map_summary_all_runs.csv`: mapas clean, near_clean y all_combinations.\n";
    readme << "- `csv/recommendation_summary.csv`: mejor corte por run y recomendacion global.\n";
    readme << "- `png/`: graficos comparativos y mapas por run.\n";
    readme << "- `compare_hv_runs_histograms.root`: histogramas y graficos principales en ROOT.\n";
    readme << "- `summary.txt`: lectura resumida y recomendacion final.\n\n";
    readme << "## Como interpretar\n\n";
    readme << "`near_clean_XY` es la muestra fisicamente mas util: permite 1 o 2 hits por plano y evita gran parte de los eventos muy poblados. `multi_large_XY` indica multiplicidad alta, probablemente ruido, crosstalk, pickup o umbral demasiado bajo.\n\n";
    readme << "`quality_score` recompensa `near_clean_XY` y `clean_XY`, penaliza `multi_large_XY`, dominancia extrema de un canal y actividad ruidosa en hod2019. Es un ranking automatico simple; debe leerse junto con ocupancias y mapas.\n\n";
    readme << "Para analisis posterior se prefiere `medium` o `strict`. `loose` sirve como diagnostico de pulsos marginales y ruido.\n\n";
    readme << "Los ROOT originales se abren solo en modo lectura y no se modifican.\n";
}

static void write_summary(const fs::path& outdir, const std::vector<RunResult>& runs) {
    std::ofstream summary(outdir / "summary.txt");
    summary << std::setprecision(8);
    summary << "Comparacion de runs HV\n";
    summary << "======================\n\n";
    summary << "Archivos analizados:\n";
    for (const RunResult& run : runs) {
        summary << "  Run " << run.run_index << " (" << run.run_name << "): "
                << run.filename << ", eventos=" << run.total_events << "\n";
    }
    summary << "\n";

    summary << "Tabla por run y corte:\n";
    summary << "run cut events frac_clean frac_near frac_multi_large mean_nHits domX domY hod2019 quality\n";
    for (const RunResult& run : runs) {
        if (!run.ok) continue;
        for (const CutSummary& s : run.summaries) {
            summary << run.run_index << ' ' << s.cut_name << ' ' << run.total_events << ' '
                    << s.frac_clean_xy << ' ' << s.frac_near_clean_xy << ' '
                    << s.frac_multi_large_xy << ' ' << s.mean_nhits << ' '
                    << "X" << s.dominant_x_channel << "(" << s.dominant_x_fraction << ") "
                    << "Y" << s.dominant_y_channel << "(" << s.dominant_y_fraction << ") "
                    << s.hod2019_status << ' ' << s.quality_score << "\n";
        }
    }
    summary << "\n";

    std::vector<const RunResult*> ranking;
    for (const RunResult& run : runs) if (run.ok) ranking.push_back(&run);
    std::sort(ranking.begin(), ranking.end(), [](const RunResult* a, const RunResult* b) {
        return a->best_score > b->best_score;
    });

    summary << "Ranking de runs por mejor corte:\n";
    for (size_t i = 0; i < ranking.size(); ++i) {
        const RunResult* run = ranking[i];
        summary << "  " << (i + 1) << ". Run " << run->run_index << " " << run->run_name
                << " best_cut=" << kCuts[run->best_cut_index].name
                << " quality_score=" << run->best_score
                << " reason=" << (run->recommended_for_analysis ? "best_available" : run->reason)
                << "\n";
    }
    summary << "\n";

    if (!ranking.empty()) {
        const RunResult* best = ranking.front();
        summary << "Recomendacion final:\n";
        summary << "El run recomendado para analisis es Run " << best->run_index
                << " (" << best->run_name << ") usando corte "
                << kCuts[best->best_cut_index].name << ".\n";
        summary << "Para analisis posterior se recomienda usar medium o strict; si el ranking eligio loose, revisar manualmente mapas y multiplicidades antes de aceptarlo.\n\n";
    }

    if (!runs.empty()) {
        const RunResult& highest_hv = runs.back();
        const CutSummary& s = highest_hv.summaries[kMedium];
        if (s.frac_multi_large_xy > 0.3 || s.mean_nhits > 8.0) {
            summary << "ADVERTENCIA: el run de mayor HV muestra multiplicidad elevada en medium; probablemente esta dominado por ruido/crosstalk/pickup.\n";
        }
    }

    for (const RunResult& run : runs) {
        if (!run.ok) continue;
        const CutSummary& loose = run.summaries[kLoose];
        const CutSummary& strict = run.summaries[kStrict];
        if (loose.frac_near_clean_xy > 0.01 && strict.frac_near_clean_xy < 0.2 * loose.frac_near_clean_xy) {
            summary << "ADVERTENCIA: Run " << run.run_index
                    << " tiene muchas senales en loose que no sobreviven a strict; probablemente son ruido o pulsos marginales.\n";
        }
    }

    bool hod2019_useful = false;
    for (const RunResult& run : runs) {
        if (run.summaries[kStrict].hod2019_status == "possible_signal") hod2019_useful = true;
    }
    if (!hod2019_useful) {
        summary << "ADVERTENCIA: hod2019 no entrega senal util consistente; no usar estos runs para tracks 3D.\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Uso:\n  " << argv[0] << " outputs/compare_hv_runs archivo1.root [archivo2.root ...]\n";
        return 1;
    }

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    const fs::path outdir = argv[1];
    const fs::path csvdir = outdir / "csv";
    const fs::path pngdir = outdir / "png";
    fs::create_directories(csvdir);
    fs::create_directories(pngdir);

    std::vector<RunResult> runs;
    runs.reserve(argc - 2);
    for (int i = 2; i < argc; ++i) {
        runs.push_back(process_run(i - 1, fs::path(argv[i])));
    }

    if (std::none_of(runs.begin(), runs.end(), [](const RunResult& r) { return r.ok; })) {
        std::cerr << "ERROR: ningun run pudo procesarse.\n";
        return 2;
    }

    write_recommendation_csv(csvdir, runs);
    write_csvs(csvdir, runs);
    make_all_plots(outdir, runs);
    write_readme(outdir);
    write_summary(outdir, runs);

    std::cout << "\nSalidas escritas en " << outdir << "\n";
    std::cout << "Resumen: " << (outdir / "summary.txt") << "\n";
    std::cout << "PNGs: " << pngdir << "\n";
    return 0;
}
