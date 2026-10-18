#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "hodo_paths.h"

constexpr int kNChannels = 64;

struct ChannelConfig {
    int baseline_start = 0;
    int baseline_end = 2500;
    int signal_start = 0;
    int signal_end = 0;
    double nsigma = 0.0;
    double abs_threshold_adc = 0.0;
    int min_consecutive = 2;
};

struct EffectiveWindow {
    int start = 0;
    int end = -1;
    bool valid = false;
    bool from_time_axis = false;
};

struct WaveformResult {
    bool valid = false;
    double baseline = 0.0;
    double sigma = 0.0;
    double vmin = 0.0;
    int tmin_index = -1;
    double tmin = -1.0;
    double amp_neg = 0.0;
    double threshold_delta = 0.0;
    bool hit_neg = false;
};

struct ChannelStats {
    int channel = -1;
    bool branch_exists = false;
    ChannelConfig cfg;

    long long entries = 0;
    long long hits_neg = 0;

    double sum_baseline = 0.0;
    double sum_sigma = 0.0;
    double sum_amp_neg_all = 0.0;
    double sum_amp_neg_hits = 0.0;
    double max_amp_neg = 0.0;

    double sum_tmin_hits = 0.0;
    double sum_tmin2_hits = 0.0;
    double min_tmin_hits = std::numeric_limits<double>::infinity();
    double max_tmin_hits = -std::numeric_limits<double>::infinity();
};

static ChannelConfig config_for_channel(int ch) {
    ChannelConfig cfg;
    cfg.baseline_start = 0;
    cfg.baseline_end = 2500;
    cfg.min_consecutive = 2;

    if (ch < 32) {
        cfg.signal_start = 3500;
        cfg.signal_end = 7000;
        cfg.nsigma = 5.0;
        cfg.abs_threshold_adc = 20.0;
    } else {
        cfg.signal_start = 3000;
        cfg.signal_end = 5000;
        cfg.nsigma = 4.0;
        cfg.abs_threshold_adc = 3.0;
    }

    return cfg;
}

static const char* detector_name(int ch) {
    return (ch < 32) ? "hod2018" : "hod2019";
}

static EffectiveWindow clamp_window(int requested_start, int requested_end, int n_samples) {
    EffectiveWindow w;
    if (n_samples <= 0) return w;

    int start = std::max(0, requested_start);
    int end = std::max(0, requested_end);
    if (start > end) std::swap(start, end);

    if (start >= n_samples) start = n_samples - 1;
    if (end >= n_samples) end = n_samples - 1;
    if (start > end) return w;

    w.start = start;
    w.end = end;
    w.valid = true;
    return w;
}

static EffectiveWindow window_from_time_axis(int requested_start, int requested_end,
                                             int n_samples,
                                             const std::vector<float>* time_axis) {
    if (!time_axis || static_cast<int>(time_axis->size()) != n_samples) {
        return clamp_window(requested_start, requested_end, n_samples);
    }

    int start = std::max(0, requested_start);
    int end = std::max(0, requested_end);
    if (start > end) std::swap(start, end);

    EffectiveWindow w;
    for (int i = 0; i < n_samples; ++i) {
        const double t = time_axis->at(i);
        if (t >= start && t <= end) {
            if (!w.valid) {
                w.start = i;
                w.valid = true;
            }
            w.end = i;
        }
    }
    w.from_time_axis = w.valid;
    return w;
}

static double mean_in_window(const std::vector<float>& v, const EffectiveWindow& w) {
    double sum = 0.0;
    int n = 0;
    for (int i = w.start; i <= w.end; ++i) {
        sum += v[i];
        ++n;
    }
    return n > 0 ? sum / static_cast<double>(n) : 0.0;
}

static double sigma_in_window(const std::vector<float>& v, const EffectiveWindow& w, double mean) {
    const int n = w.end - w.start + 1;
    if (n <= 1) return 0.0;

    double sum2 = 0.0;
    for (int i = w.start; i <= w.end; ++i) {
        const double d = v[i] - mean;
        sum2 += d * d;
    }
    return std::sqrt(sum2 / static_cast<double>(n - 1));
}

static bool has_consecutive_below(const std::vector<float>& v, const EffectiveWindow& w,
                                  double threshold, int min_consecutive) {
    if (!w.valid || min_consecutive <= 0) return false;

    int consecutive = 0;
    for (int i = w.start; i <= w.end; ++i) {
        if (v[i] < threshold) {
            ++consecutive;
            if (consecutive >= min_consecutive) return true;
        } else {
            consecutive = 0;
        }
    }
    return false;
}

static WaveformResult analyze_waveform_windowed(const std::vector<float>& v,
                                                const std::vector<float>* time_axis,
                                                const ChannelConfig& cfg) {
    WaveformResult r;
    if (v.empty()) return r;

    const int n_samples = static_cast<int>(v.size());
    const EffectiveWindow baseline_window = window_from_time_axis(cfg.baseline_start, cfg.baseline_end,
                                                                  n_samples, time_axis);
    const EffectiveWindow signal_window = window_from_time_axis(cfg.signal_start, cfg.signal_end,
                                                                n_samples, time_axis);
    if (!baseline_window.valid || !signal_window.valid) return r;

    r.valid = true;
    r.baseline = mean_in_window(v, baseline_window);
    r.sigma = sigma_in_window(v, baseline_window, r.baseline);

    r.vmin = std::numeric_limits<double>::infinity();
    r.tmin_index = signal_window.start;
    for (int i = signal_window.start; i <= signal_window.end; ++i) {
        if (v[i] < r.vmin) {
            r.vmin = v[i];
            r.tmin_index = i;
        }
    }
    if (time_axis && static_cast<int>(time_axis->size()) == n_samples) {
        r.tmin = time_axis->at(r.tmin_index);
    } else {
        r.tmin = static_cast<double>(r.tmin_index);
    }

    r.amp_neg = r.baseline - r.vmin;
    r.threshold_delta = std::max(cfg.nsigma * r.sigma, cfg.abs_threshold_adc);
    const double threshold = r.baseline - r.threshold_delta;

    r.hit_neg = (r.amp_neg > r.threshold_delta) &&
                has_consecutive_below(v, signal_window, threshold, cfg.min_consecutive);
    return r;
}

static void update_stats(ChannelStats& stats, const WaveformResult& r) {
    if (!r.valid) return;

    ++stats.entries;
    stats.sum_baseline += r.baseline;
    stats.sum_sigma += r.sigma;
    stats.sum_amp_neg_all += r.amp_neg;
    stats.max_amp_neg = std::max(stats.max_amp_neg, r.amp_neg);

    if (r.hit_neg) {
        ++stats.hits_neg;
        stats.sum_amp_neg_hits += r.amp_neg;
        stats.sum_tmin_hits += r.tmin;
        stats.sum_tmin2_hits += r.tmin * r.tmin;
        stats.min_tmin_hits = std::min(stats.min_tmin_hits, r.tmin);
        stats.max_tmin_hits = std::max(stats.max_tmin_hits, r.tmin);
    }
}

static std::string join_channels(const std::vector<int>& channels) {
    std::ostringstream os;
    for (size_t i = 0; i < channels.size(); ++i) {
        if (i) os << ';';
        os << channels[i];
    }
    return os.str();
}

static std::string leaf_value_or_blank(TLeaf* leaf) {
    if (!leaf) return "";
    return std::to_string(leaf->GetValueLong64(0));
}

static double fraction(long long num, long long den) {
    return den > 0 ? static_cast<double>(num) / static_cast<double>(den) : 0.0;
}

static double mean_or_zero(double sum, long long n) {
    return n > 0 ? sum / static_cast<double>(n) : 0.0;
}

static double std_tmin_hits(const ChannelStats& s) {
    if (s.hits_neg <= 1) return 0.0;
    const double mean = s.sum_tmin_hits / static_cast<double>(s.hits_neg);
    const double mean2 = s.sum_tmin2_hits / static_cast<double>(s.hits_neg);
    return std::sqrt(std::max(0.0, mean2 - mean * mean));
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Uso:\n"
                  << "  " << argv[0] << " lectura_vx1740_ext_nim_1M_v4.root [max_events]\n";
        return 1;
    }

    const std::string filename = hodo::resolveRawFile(argv[1]);
    long long max_events = -1;
    if (argc >= 3) max_events = std::atoll(argv[2]);

    gROOT->SetBatch(kTRUE);

    std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: no pude abrir " << filename << "\n";
        return 2;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("data_dig"));
    if (!tree) {
        std::cerr << "ERROR: no encontre el TTree data_dig\n";
        return 3;
    }

    long long nentries = tree->GetEntries();
    if (max_events >= 0 && max_events < nentries) nentries = max_events;

    std::cout << "Archivo: " << filename << "\n";
    std::cout << "Eventos en data_dig: " << tree->GetEntries() << "\n";
    std::cout << "Eventos a procesar: " << nentries << "\n";

    tree->SetBranchStatus("*", 0);
    tree->SetCacheSize(64LL * 1024LL * 1024LL);

    std::array<std::vector<float>*, kNChannels> channels{};
    std::vector<float>* time_axis = nullptr;
    std::array<ChannelStats, kNChannels> stats{};

    for (int ch = 0; ch < kNChannels; ++ch) {
        stats[ch].channel = ch;
        stats[ch].cfg = config_for_channel(ch);

        const std::string branch_name = "Ch" + std::to_string(ch);
        TBranch* branch = tree->GetBranch(branch_name.c_str());
        if (!branch) {
            std::cerr << "ADVERTENCIA: falta rama " << branch_name << "\n";
            continue;
        }

        stats[ch].branch_exists = true;
        tree->SetBranchStatus(branch_name.c_str(), 1);
        if (tree->SetBranchAddress(branch_name.c_str(), &channels[ch]) < 0) {
            std::cerr << "ADVERTENCIA: no pude conectar rama " << branch_name << "\n";
            stats[ch].branch_exists = false;
        }
    }

    TLeaf* event_leaf = nullptr;
    TLeaf* timestamp_leaf = nullptr;

    if (tree->GetBranch("event")) {
        tree->SetBranchStatus("event", 1);
        event_leaf = tree->GetLeaf("event");
    } else {
        std::cerr << "ADVERTENCIA: falta rama event\n";
    }

    if (tree->GetBranch("timestamp")) {
        tree->SetBranchStatus("timestamp", 1);
        timestamp_leaf = tree->GetLeaf("timestamp");
    } else {
        std::cerr << "ADVERTENCIA: falta rama timestamp\n";
    }

    if (tree->GetBranch("time")) {
        tree->SetBranchStatus("time", 1);
        if (tree->SetBranchAddress("time", &time_axis) < 0) {
            std::cerr << "ADVERTENCIA: no pude conectar rama time; usare indices de muestra\n";
            time_axis = nullptr;
        }
    } else {
        std::cerr << "ADVERTENCIA: falta rama time; usare indices de muestra\n";
    }

    tree->AddBranchToCache("*", kTRUE);

    std::ofstream event_csv("event_multiplicity_windowed.csv");
    event_csv << "entry,event,timestamp,"
              << "n_hit_neg_total,"
              << "n_hit_neg_hod2018,"
              << "n_hit_neg_hod2019,"
              << "n_hit_neg_hod2018_x_or_y_rough,"
              << "n_hit_neg_hod2019_x_or_y_rough\n";

    std::ofstream candidate_csv("candidate_events_windowed.csv");
    candidate_csv << "entry,event,timestamp,"
                  << "n_hit_hod2018,n_hit_hod2019,"
                  << "channels_hod2018,channels_hod2019,"
                  << "max_amp_hod2018,max_amp_hod2019\n";

    long long total_hits_hod2018 = 0;
    long long total_hits_hod2019 = 0;
    long long events_with_hod2018 = 0;
    long long events_with_hod2019 = 0;
    long long events_with_both = 0;

    for (long long entry = 0; entry < nentries; ++entry) {
        tree->GetEntry(entry);

        if (entry % 1000 == 0) {
            std::cout << "Procesando evento " << entry << " / " << nentries << "\r" << std::flush;
        }

        std::vector<int> hit_channels_hod2018;
        std::vector<int> hit_channels_hod2019;
        double max_amp_hod2018 = 0.0;
        double max_amp_hod2019 = 0.0;

        for (int ch = 0; ch < kNChannels; ++ch) {
            if (!stats[ch].branch_exists) continue;
            if (channels[ch] == nullptr || channels[ch]->empty()) continue;

            const std::vector<float>* axis = nullptr;
            if (time_axis && time_axis->size() == channels[ch]->size()) axis = time_axis;
            const WaveformResult r = analyze_waveform_windowed(*channels[ch], axis, stats[ch].cfg);
            if (!r.valid) continue;

            update_stats(stats[ch], r);

            if (!r.hit_neg) continue;

            if (ch < 32) {
                hit_channels_hod2018.push_back(ch);
                max_amp_hod2018 = std::max(max_amp_hod2018, r.amp_neg);
            } else {
                hit_channels_hod2019.push_back(ch);
                max_amp_hod2019 = std::max(max_amp_hod2019, r.amp_neg);
            }
        }

        const int n_hod2018 = static_cast<int>(hit_channels_hod2018.size());
        const int n_hod2019 = static_cast<int>(hit_channels_hod2019.size());
        const int n_total = n_hod2018 + n_hod2019;

        total_hits_hod2018 += n_hod2018;
        total_hits_hod2019 += n_hod2019;
        if (n_hod2018 > 0) ++events_with_hod2018;
        if (n_hod2019 > 0) ++events_with_hod2019;
        if (n_hod2018 > 0 && n_hod2019 > 0) ++events_with_both;

        const std::string event_value = leaf_value_or_blank(event_leaf);
        const std::string timestamp_value = leaf_value_or_blank(timestamp_leaf);

        event_csv << entry << ','
                  << event_value << ','
                  << timestamp_value << ','
                  << n_total << ','
                  << n_hod2018 << ','
                  << n_hod2019 << ','
                  << n_hod2018 << ','
                  << n_hod2019 << '\n';

        if (n_hod2018 > 0 && n_hod2019 > 0) {
            candidate_csv << entry << ','
                          << event_value << ','
                          << timestamp_value << ','
                          << n_hod2018 << ','
                          << n_hod2019 << ','
                          << join_channels(hit_channels_hod2018) << ','
                          << join_channels(hit_channels_hod2019) << ','
                          << max_amp_hod2018 << ','
                          << max_amp_hod2019 << '\n';
        }
    }

    std::cout << "\nProcesamiento terminado.\n";
    event_csv.close();
    candidate_csv.close();

    std::ofstream summary_csv("active_channels_windowed_summary.csv");
    summary_csv << "channel,detector,entries,hits_neg,frac_hits_neg,"
                << "mean_baseline,mean_sigma,mean_amp_neg_all,mean_amp_neg_hits,max_amp_neg,"
                << "mean_tmin_hits,std_tmin_hits,min_tmin_hits,max_tmin_hits,"
                << "baseline_start,baseline_end,signal_start,signal_end,"
                << "nsigma,abs_threshold_adc,min_consecutive\n";
    summary_csv << std::setprecision(10);

    std::vector<ChannelStats> sorted(stats.begin(), stats.end());
    std::sort(sorted.begin(), sorted.end(), [](const ChannelStats& a, const ChannelStats& b) {
        if (a.hits_neg != b.hits_neg) return a.hits_neg > b.hits_neg;
        return a.channel < b.channel;
    });

    for (const ChannelStats& s : stats) {
        const double min_t = (s.hits_neg > 0) ? s.min_tmin_hits : -1.0;
        const double max_t = (s.hits_neg > 0) ? s.max_tmin_hits : -1.0;
        summary_csv << s.channel << ','
                    << detector_name(s.channel) << ','
                    << s.entries << ','
                    << s.hits_neg << ','
                    << fraction(s.hits_neg, s.entries) << ','
                    << mean_or_zero(s.sum_baseline, s.entries) << ','
                    << mean_or_zero(s.sum_sigma, s.entries) << ','
                    << mean_or_zero(s.sum_amp_neg_all, s.entries) << ','
                    << mean_or_zero(s.sum_amp_neg_hits, s.hits_neg) << ','
                    << s.max_amp_neg << ','
                    << mean_or_zero(s.sum_tmin_hits, s.hits_neg) << ','
                    << std_tmin_hits(s) << ','
                    << min_t << ','
                    << max_t << ','
                    << s.cfg.baseline_start << ','
                    << s.cfg.baseline_end << ','
                    << s.cfg.signal_start << ','
                    << s.cfg.signal_end << ','
                    << s.cfg.nsigma << ','
                    << s.cfg.abs_threshold_adc << ','
                    << s.cfg.min_consecutive << '\n';
    }
    summary_csv.close();

    std::cout << "\nResumen global:\n";
    std::cout << "  total eventos procesados: " << nentries << "\n";
    std::cout << "  hits totales hod2018: " << total_hits_hod2018 << "\n";
    std::cout << "  hits totales hod2019: " << total_hits_hod2019 << "\n";
    std::cout << "  eventos con al menos un hit hod2018: " << events_with_hod2018 << "\n";
    std::cout << "  eventos con al menos un hit hod2019: " << events_with_hod2019 << "\n";
    std::cout << "  eventos con hits en ambos hodoscopios: " << events_with_both << "\n";

    std::cout << "\nTabla ordenada por hits_neg:\n";
    std::cout << std::setw(8) << "channel"
              << std::setw(10) << "detector"
              << std::setw(12) << "entries"
              << std::setw(12) << "hits_neg"
              << std::setw(14) << "frac_hits"
              << std::setw(14) << "mean_base"
              << std::setw(12) << "sigma"
              << std::setw(14) << "mean_amp"
              << std::setw(14) << "max_amp"
              << std::setw(12) << "mean_tmin"
              << "\n";

    std::cout << std::fixed << std::setprecision(5);
    for (const ChannelStats& s : sorted) {
        std::cout << std::setw(8) << s.channel
                  << std::setw(10) << detector_name(s.channel)
                  << std::setw(12) << s.entries
                  << std::setw(12) << s.hits_neg
                  << std::setw(14) << fraction(s.hits_neg, s.entries)
                  << std::setw(14) << mean_or_zero(s.sum_baseline, s.entries)
                  << std::setw(12) << mean_or_zero(s.sum_sigma, s.entries)
                  << std::setw(14) << mean_or_zero(s.sum_amp_neg_hits, s.hits_neg)
                  << std::setw(14) << s.max_amp_neg
                  << std::setw(12) << mean_or_zero(s.sum_tmin_hits, s.hits_neg)
                  << "\n";
    }

    std::cout << "\nArchivos escritos:\n";
    std::cout << "  active_channels_windowed_summary.csv\n";
    std::cout << "  event_multiplicity_windowed.csv\n";
    std::cout << "  candidate_events_windowed.csv\n";

    return 0;
}
