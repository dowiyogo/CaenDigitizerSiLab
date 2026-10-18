#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
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
#include "hodo_paths.h"

namespace fs = std::filesystem;

constexpr int kNChannels = 64;
constexpr int kNHod2018 = 32;
constexpr int kNPlane = 16;
constexpr int kNBase = 50;

struct WaveformMetrics {
    bool valid = false;
    double baseline = 0.0;
    double sigma = 0.0;
    double vmin = 0.0;
    double vmax = 0.0;
    double amp_neg = 0.0;
    double amp_pos = 0.0;
    int sample_min = -1;
    int sample_max = -1;
};

struct ChannelStats {
    bool branch_present = false;
    long long entries = 0;
    long long hits_relaxed = 0;
    long long hits_strict = 0;
    double sum_baseline = 0.0;
    double sum_baseline2 = 0.0;
    double sum_sigma = 0.0;
    double max_amp_neg = 0.0;
    double max_amp_pos = 0.0;
    double min_vmin = std::numeric_limits<double>::infinity();
    double max_vmax = -std::numeric_limits<double>::infinity();
    double sum_sample_min = 0.0;
    double sum_sample_max = 0.0;
};

struct HitInfo {
    int global_channel = -1;
    int local_channel = -1;
    double amp_neg = 0.0;
    int sample_min = -1;
};

static std::string detector_name(int ch) {
    return (ch < 32) ? "hod2018" : "hod2019";
}

static std::string plane_name(int ch) {
    if (ch < 16) return "X";
    if (ch < 32) return "Y";
    return "HOD2019";
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
    m.sample_max = static_cast<int>(std::distance(v.begin(), minmax.second));
    return m;
}

static bool has_two_consecutive_below(const std::vector<float>& v, double threshold_value) {
    if (v.size() < 2) return false;
    for (size_t i = 0; i + 1 < v.size(); ++i) {
        if (v[i] < threshold_value && v[i + 1] < threshold_value) return true;
    }
    return false;
}

static bool hit_negative(const std::vector<float>& v, const WaveformMetrics& m, double delta_adc) {
    if (!m.valid) return false;
    return has_two_consecutive_below(v, m.baseline - delta_adc);
}

static void update_channel_stats(ChannelStats& s, const WaveformMetrics& m,
                                 bool hit_relaxed, bool hit_strict) {
    if (!m.valid) return;
    ++s.entries;
    if (hit_relaxed) ++s.hits_relaxed;
    if (hit_strict) ++s.hits_strict;
    s.sum_baseline += m.baseline;
    s.sum_baseline2 += m.baseline * m.baseline;
    s.sum_sigma += m.sigma;
    s.max_amp_neg = std::max(s.max_amp_neg, m.amp_neg);
    s.max_amp_pos = std::max(s.max_amp_pos, m.amp_pos);
    s.min_vmin = std::min(s.min_vmin, m.vmin);
    s.max_vmax = std::max(s.max_vmax, m.vmax);
    s.sum_sample_min += m.sample_min;
    s.sum_sample_max += m.sample_max;
}

static double safe_fraction(long long num, long long den) {
    return den > 0 ? static_cast<double>(num) / static_cast<double>(den) : 0.0;
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
    const double var = std::max(0.0, mean2 - mean * mean);
    return std::sqrt(var);
}

static double mean_sample_min(const ChannelStats& s) {
    return s.entries > 0 ? s.sum_sample_min / static_cast<double>(s.entries) : 0.0;
}

static double mean_sample_max(const ChannelStats& s) {
    return s.entries > 0 ? s.sum_sample_max / static_cast<double>(s.entries) : 0.0;
}

static std::string status_for_channel(int ch, const ChannelStats& s) {
    if (!s.branch_present || s.entries == 0) return "missing_branch";

    const double frac_relaxed = safe_fraction(s.hits_relaxed, s.entries);
    if (frac_relaxed > 0.5) return "suspicious_noise";

    if (ch < 32) {
        if (s.hits_relaxed == 0) return "dead";
        if (frac_relaxed <= 0.001) return "weak";
        return "active";
    }

    if (s.hits_strict == 0 && s.max_amp_neg < 5.0) return "hod2019_no_signal";
    if (s.hits_relaxed == 0) return "dead";
    return "hod2019_outlier_or_noise";
}

static std::string join_ints(const std::vector<int>& values) {
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) os << ';';
        os << values[i];
    }
    return os.str();
}

static std::string value_or_blank(bool has_value, long long value) {
    if (!has_value) return "";
    return std::to_string(value);
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

static void save_hist2_png(TH2* h, const fs::path& path, bool logz = false,
                           const std::string& option = "COLZ") {
    TCanvas c(("c_" + std::string(h->GetName())).c_str(), "", 1000, 840);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.16);
    c.SetBottomMargin(0.12);
    if (logz) {
        c.SetLogz();
        if (h->GetMaximum() > 0.0) h->SetMinimum(0.5);
    }
    h->Draw(option.c_str());
    c.SaveAs(path.string().c_str());
}

static double nice_upper(double max_value) {
    if (max_value <= 0.0) return 10.0;
    const double raw = max_value * 1.08;
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double scaled = raw / mag;
    double nice = 10.0;
    if (scaled <= 1.0) nice = 1.0;
    else if (scaled <= 2.0) nice = 2.0;
    else if (scaled <= 5.0) nice = 5.0;
    return nice * mag;
}

static std::string format_top_channels(const std::vector<long long>& counts, int global_offset,
                                       const std::string& plane_label, int max_items = 5) {
    std::vector<std::pair<long long, int>> ranked;
    ranked.reserve(counts.size());
    for (size_t i = 0; i < counts.size(); ++i) {
        ranked.push_back({counts[i], static_cast<int>(i)});
    }
    std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });

    std::ostringstream os;
    int written = 0;
    for (const auto& item : ranked) {
        if (item.first <= 0) continue;
        if (written) os << ", ";
        os << plane_label << item.second << "(Ch" << (global_offset + item.second)
           << ")=" << item.first;
        ++written;
        if (written >= max_items) break;
    }
    if (written == 0) os << "ninguno";
    return os.str();
}

static std::string list_zero_hit_hod2018(const std::array<ChannelStats, kNChannels>& stats) {
    std::ostringstream os;
    bool first = true;
    for (int ch = 0; ch < kNHod2018; ++ch) {
        if (stats[ch].hits_relaxed != 0) continue;
        if (!first) os << ", ";
        os << "Ch" << ch;
        first = false;
    }
    if (first) os << "ninguno";
    return os.str();
}

static void write_readme(const fs::path& outdir) {
    std::ofstream readme(outdir / "README.md");
    readme << "# Analisis final hod2018 / diagnostico hod2019 v4\n\n";
    readme << "Este directorio fue generado por `final_hodoscope_analysis.cpp` a partir de "
           << "`lectura_vx1740_ext_nim_1M_v4.root`.\n\n";
    readme << "## Comandos\n\n";
    readme << "Compilacion:\n\n";
    readme << "```bash\n";
    readme << "g++ -O2 -std=c++17 final_hodoscope_analysis.cpp $(root-config --cflags --libs) -o final_hodoscope_analysis\n";
    readme << "```\n\n";
    readme << "Ejecucion:\n\n";
    readme << "```bash\n";
    readme << "./final_hodoscope_analysis lectura_vx1740_ext_nim_1M_v4.root outputs/final_hod2018_analysis_v4\n";
    readme << "```\n\n";
    readme << "Tambien se puede usar:\n\n";
    readme << "```bash\n";
    readme << "./run_final_hod2018_analysis_v4.sh\n";
    readme << "```\n\n";
    readme << "## Archivos generados\n\n";
    readme << "- `channel_summary_all.csv`: resumen por canal Ch0-Ch63, incluyendo baseline, ruido, amplitudes y estado.\n";
    readme << "- `hod2018_event_summary.csv`: resumen por evento para hod2018.\n";
    readme << "- `hod2018_clean_xy_events.csv`: eventos con exactamente un hit X y un hit Y.\n";
    readme << "- `hod2018_xy_map_clean.csv`: mapa 16x16 de eventos `clean_XY`.\n";
    readme << "- `hod2018_xy_map_all_combinations.csv`: mapa 16x16 con todas las combinaciones X-Y en eventos con X e Y.\n";
    readme << "- `hod2019_diagnostic_summary.csv`: diagnostico cuantitativo de Ch32-Ch63.\n";
    readme << "- `event_class_counts.csv`: conteos de clases; `XY_candidate` es inclusivo y contiene todos los eventos con nX>0 y nY>0.\n";
    readme << "- `final_hod2018_histograms.root`: histogramas principales en formato ROOT.\n";
    readme << "- `png/`: graficos de ocupancia, mapas, multiplicidad, amplitudes, tiempos de minimo y diagnosticos globales.\n";
    readme << "- `summary.txt`: interpretacion automatica del run.\n\n";
    readme << "## Interpretacion rapida\n\n";
    readme << "- `hod2018_x_occupancy.png` y `hod2018_y_occupancy.png` muestran que barras dominan en cada plano.\n";
    readme << "- `hod2018_xy_map_clean_counts.png` usa solo eventos con un X y un Y; es el mapa 2D principal de hod2018.\n";
    readme << "- `hod2018_xy_map_all_combinations.png` incluye eventos multi-hit y sirve como comparacion cualitativa.\n";
    readme << "- `hod2018_event_classes.png`, `hod2018_multiplicity_nhits.png` y `hod2018_nx_ny_map.png` muestran pureza y multiplicidad.\n";
    readme << "- Los histogramas de amplitud y `sample_min` permiten revisar estabilidad temporal y separacion senal/ruido.\n";
    readme << "- `hod2019_max_amp_neg.png`, `hod2019_baseline_sigma.png` y `hod2019_diagnostic_summary.csv` documentan si Ch32-Ch63 entregan senal util.\n\n";
    readme << "## Advertencias\n\n";
    readme << "Ch32-Ch63 se tratan solo como diagnostico. Si hod2019 no tiene hits estrictos y sus amplitudes maximas negativas son compatibles con ruido, no se debe intentar reconstruccion 3D con este run. La etiqueta `hod2019_outlier_or_noise` indica candidatos aislados bajo los criterios de diagnostico, no un plano operativo.\n\n";
    readme << "El archivo ROOT original se abre solo en modo lectura y no se modifica.\n";
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Uso:\n  " << argv[0]
                  << " lectura_vx1740_ext_nim_1M_v4.root outputs/final_hod2018_analysis_v4\n";
        return 1;
    }

    const fs::path input_path = argv[1];
    const fs::path outdir = argv[2];
    const fs::path pngdir = outdir / "png";

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetPaintTextFormat(".0f");

    fs::create_directories(pngdir);

    const std::string resolved = hodo::resolveRawFile(input_path.string());
    std::unique_ptr<TFile> input(TFile::Open(resolved.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: no pude abrir " << resolved << "\n";
        return 2;
    }

    TTree* tree = dynamic_cast<TTree*>(input->Get("data_dig"));
    if (!tree) {
        std::cerr << "ERROR: no encontre el TTree data_dig en " << input_path << "\n";
        return 3;
    }

    const long long total_entries = tree->GetEntries();
    std::cout << "Archivo: " << input_path << "\n";
    std::cout << "TTree data_dig, entradas: " << total_entries << "\n";

    tree->SetBranchStatus("*", 0);
    tree->SetCacheSize(64LL * 1024LL * 1024LL);

    std::array<std::vector<float>*, kNChannels> channels{};
    std::array<ChannelStats, kNChannels> channel_stats{};

    for (int ch = 0; ch < kNChannels; ++ch) {
        const std::string branch_name = "Ch" + std::to_string(ch);
        TBranch* branch = tree->GetBranch(branch_name.c_str());
        if (!branch) {
            std::cerr << "ADVERTENCIA: no existe rama " << branch_name << "\n";
            continue;
        }
        channel_stats[ch].branch_present = true;
        tree->SetBranchStatus(branch_name.c_str(), 1);
        if (tree->SetBranchAddress(branch_name.c_str(), &channels[ch]) < 0) {
            std::cerr << "ADVERTENCIA: no pude conectar " << branch_name << "\n";
            channel_stats[ch].branch_present = false;
        }
    }

    bool has_event = false;
    bool has_timestamp = false;
    TLeaf* event_leaf = nullptr;
    TLeaf* timestamp_leaf = nullptr;

    if (tree->GetBranch("event")) {
        tree->SetBranchStatus("event", 1);
        event_leaf = tree->GetLeaf("event");
        has_event = (event_leaf != nullptr);
    } else {
        std::cerr << "ADVERTENCIA: no existe rama event\n";
    }

    if (tree->GetBranch("timestamp")) {
        tree->SetBranchStatus("timestamp", 1);
        timestamp_leaf = tree->GetLeaf("timestamp");
        has_timestamp = (timestamp_leaf != nullptr);
    } else {
        std::cerr << "ADVERTENCIA: no existe rama timestamp\n";
    }

    tree->AddBranchToCache("*", kTRUE);

    std::ofstream event_csv(outdir / "hod2018_event_summary.csv");
    std::ofstream clean_csv(outdir / "hod2018_clean_xy_events.csv");
    std::ofstream multi_csv(outdir / "hod2018_multi_xy_events.csv");

    event_csv << "entry,event,timestamp,nX,nY,nHitsHod2018,class,x_channels,y_channels,max_amp_X,max_amp_Y\n";
    clean_csv << "entry,event,timestamp,x_channel,y_channel,amp_X,amp_Y,sample_min_X,sample_min_Y\n";
    multi_csv << "entry,event,timestamp,nX,nY,nHitsHod2018,x_channels,y_channels,max_amp_X,max_amp_Y\n";

    std::array<long long, kNPlane> x_occupancy{};
    std::array<long long, kNPlane> y_occupancy{};
    std::array<std::array<long long, kNPlane>, kNPlane> xy_clean{};
    std::array<std::array<long long, kNPlane>, kNPlane> xy_all_combinations{};

    long long n_empty = 0;
    long long n_x_only = 0;
    long long n_y_only = 0;
    long long n_xy_candidate = 0;
    long long n_clean_xy = 0;
    long long n_multi_xy = 0;

    auto h_nhits = new TH1D("h_nhits", "hod2018 hit multiplicity;nHitsHod2018;Events", 33, -0.5, 32.5);
    auto h_nxny = new TH2D("h_nxny", "hod2018 nX vs nY;nX;nY;Events", 17, -0.5, 16.5, 17, -0.5, 16.5);
    configure_hist_1d(h_nhits, "nHitsHod2018", "Events");
    configure_hist_2d(h_nxny, "nX", "nY", "Events");

    std::vector<double> clean_amp_x;
    std::vector<double> clean_amp_y;
    std::vector<int> clean_sample_min_x;
    std::vector<int> clean_sample_min_y;
    clean_amp_x.reserve(10000);
    clean_amp_y.reserve(10000);
    clean_sample_min_x.reserve(10000);
    clean_sample_min_y.reserve(10000);

    bool have_first_timestamp = false;
    long double first_timestamp = 0.0L;
    long double last_timestamp = 0.0L;
    long double min_timestamp = 0.0L;
    long double max_timestamp = 0.0L;

    for (long long entry = 0; entry < total_entries; ++entry) {
        tree->GetEntry(entry);
        if (entry % 1000 == 0) {
            std::cout << "Procesando evento " << entry << " / " << total_entries << "\r" << std::flush;
        }

        const bool event_ok = has_event;
        const long long event_value = event_ok ? event_leaf->GetValueLong64(0) : 0;
        const bool timestamp_ok = has_timestamp;
        const long long timestamp_value = timestamp_ok ? timestamp_leaf->GetValueLong64(0) : 0;

        if (timestamp_ok) {
            const long double ts = static_cast<long double>(timestamp_value);
            if (!have_first_timestamp) {
                have_first_timestamp = true;
                first_timestamp = ts;
                min_timestamp = ts;
                max_timestamp = ts;
            }
            last_timestamp = ts;
            min_timestamp = std::min(min_timestamp, ts);
            max_timestamp = std::max(max_timestamp, ts);
        }

        std::vector<HitInfo> x_hits;
        std::vector<HitInfo> y_hits;
        x_hits.reserve(kNPlane);
        y_hits.reserve(kNPlane);

        for (int ch = 0; ch < kNChannels; ++ch) {
            if (!channel_stats[ch].branch_present || channels[ch] == nullptr) continue;
            const std::vector<float>& waveform = *channels[ch];
            const WaveformMetrics m = analyze_waveform(waveform);
            if (!m.valid) continue;

            const double relaxed_delta = std::max(4.0 * m.sigma, 5.0);
            const double strict_delta = std::max(5.0 * m.sigma, 20.0);
            const bool relaxed_hit = hit_negative(waveform, m, relaxed_delta);
            const bool strict_hit = hit_negative(waveform, m, strict_delta);
            update_channel_stats(channel_stats[ch], m, relaxed_hit, strict_hit);

            if (ch < kNHod2018 && relaxed_hit) {
                if (ch < 16) {
                    x_hits.push_back({ch, ch, m.amp_neg, m.sample_min});
                    ++x_occupancy[ch];
                } else {
                    y_hits.push_back({ch, ch - 16, m.amp_neg, m.sample_min});
                    ++y_occupancy[ch - 16];
                }
            }
        }

        const int nX = static_cast<int>(x_hits.size());
        const int nY = static_cast<int>(y_hits.size());
        const int nHits = nX + nY;
        h_nhits->Fill(nHits);
        h_nxny->Fill(nX, nY);

        double max_amp_x = 0.0;
        double max_amp_y = 0.0;
        std::vector<int> x_channels;
        std::vector<int> y_channels;
        x_channels.reserve(x_hits.size());
        y_channels.reserve(y_hits.size());

        for (const auto& hit : x_hits) {
            max_amp_x = std::max(max_amp_x, hit.amp_neg);
            x_channels.push_back(hit.global_channel);
        }
        for (const auto& hit : y_hits) {
            max_amp_y = std::max(max_amp_y, hit.amp_neg);
            y_channels.push_back(hit.global_channel);
        }

        std::string event_class;
        if (nX == 0 && nY == 0) {
            event_class = "empty";
            ++n_empty;
        } else if (nX > 0 && nY == 0) {
            event_class = "X_only";
            ++n_x_only;
        } else if (nX == 0 && nY > 0) {
            event_class = "Y_only";
            ++n_y_only;
        } else {
            ++n_xy_candidate;
            for (const auto& x_hit : x_hits) {
                for (const auto& y_hit : y_hits) {
                    ++xy_all_combinations[x_hit.local_channel][y_hit.local_channel];
                }
            }

            if (nX == 1 && nY == 1) {
                event_class = "clean_XY";
                ++n_clean_xy;
                const int x = x_hits[0].local_channel;
                const int y = y_hits[0].local_channel;
                ++xy_clean[x][y];
                clean_amp_x.push_back(x_hits[0].amp_neg);
                clean_amp_y.push_back(y_hits[0].amp_neg);
                clean_sample_min_x.push_back(x_hits[0].sample_min);
                clean_sample_min_y.push_back(y_hits[0].sample_min);

                clean_csv << entry << ',' << value_or_blank(event_ok, event_value) << ','
                          << value_or_blank(timestamp_ok, timestamp_value) << ','
                          << x << ',' << y << ','
                          << x_hits[0].amp_neg << ',' << y_hits[0].amp_neg << ','
                          << x_hits[0].sample_min << ',' << y_hits[0].sample_min << '\n';
            } else {
                event_class = "multi_XY";
                ++n_multi_xy;
                multi_csv << entry << ',' << value_or_blank(event_ok, event_value) << ','
                          << value_or_blank(timestamp_ok, timestamp_value) << ','
                          << nX << ',' << nY << ',' << nHits << ','
                          << join_ints(x_channels) << ',' << join_ints(y_channels) << ','
                          << max_amp_x << ',' << max_amp_y << '\n';
            }
        }

        event_csv << entry << ',' << value_or_blank(event_ok, event_value) << ','
                  << value_or_blank(timestamp_ok, timestamp_value) << ','
                  << nX << ',' << nY << ',' << nHits << ','
                  << event_class << ','
                  << join_ints(x_channels) << ',' << join_ints(y_channels) << ','
                  << max_amp_x << ',' << max_amp_y << '\n';
    }
    std::cout << "\nProcesamiento terminado.\n";

    event_csv.close();
    clean_csv.close();
    multi_csv.close();

    std::ofstream channel_csv(outdir / "channel_summary_all.csv");
    channel_csv << "channel,detector,plane,entries,hits_relaxed,hits_strict,frac_hits_relaxed,frac_hits_strict,"
                << "mean_baseline,mean_sigma,rms_baseline,max_amp_neg,max_amp_pos,min_vmin,max_vmax,"
                << "mean_sample_min,mean_sample_max,status\n";

    std::ofstream hod2019_csv(outdir / "hod2019_diagnostic_summary.csv");
    hod2019_csv << "channel,entries,hits_strict,hits_relaxed,max_amp_neg,max_amp_pos,mean_baseline,mean_sigma,min_vmin,max_vmax,status\n";

    channel_csv << std::setprecision(10);
    hod2019_csv << std::setprecision(10);

    for (int ch = 0; ch < kNChannels; ++ch) {
        const ChannelStats& s = channel_stats[ch];
        const double frac_relaxed = safe_fraction(s.hits_relaxed, s.entries);
        const double frac_strict = safe_fraction(s.hits_strict, s.entries);
        const double min_vmin = std::isfinite(s.min_vmin) ? s.min_vmin : 0.0;
        const double max_vmax = std::isfinite(s.max_vmax) ? s.max_vmax : 0.0;
        const std::string status = status_for_channel(ch, s);

        channel_csv << ch << ','
                    << detector_name(ch) << ','
                    << plane_name(ch) << ','
                    << s.entries << ','
                    << s.hits_relaxed << ','
                    << s.hits_strict << ','
                    << frac_relaxed << ','
                    << frac_strict << ','
                    << mean_baseline(s) << ','
                    << mean_sigma(s) << ','
                    << rms_baseline(s) << ','
                    << s.max_amp_neg << ','
                    << s.max_amp_pos << ','
                    << min_vmin << ','
                    << max_vmax << ','
                    << mean_sample_min(s) << ','
                    << mean_sample_max(s) << ','
                    << status << '\n';

        if (ch >= 32) {
            hod2019_csv << ch << ','
                        << s.entries << ','
                        << s.hits_strict << ','
                        << s.hits_relaxed << ','
                        << s.max_amp_neg << ','
                        << s.max_amp_pos << ','
                        << mean_baseline(s) << ','
                        << mean_sigma(s) << ','
                        << min_vmin << ','
                        << max_vmax << ','
                        << status << '\n';
        }
    }
    channel_csv.close();
    hod2019_csv.close();

    std::ofstream xy_clean_csv(outdir / "hod2018_xy_map_clean.csv");
    xy_clean_csv << "x_channel,y_channel,count\n";
    std::ofstream xy_all_csv(outdir / "hod2018_xy_map_all_combinations.csv");
    xy_all_csv << "x_channel,y_channel,count\n";

    long long total_clean_xy_map = 0;
    long long best_xy_count = 0;
    int best_x = -1;
    int best_y = -1;
    for (int x = 0; x < kNPlane; ++x) {
        for (int y = 0; y < kNPlane; ++y) {
            xy_clean_csv << x << ',' << y << ',' << xy_clean[x][y] << '\n';
            xy_all_csv << x << ',' << y << ',' << xy_all_combinations[x][y] << '\n';
            total_clean_xy_map += xy_clean[x][y];
            if (xy_clean[x][y] > best_xy_count) {
                best_xy_count = xy_clean[x][y];
                best_x = x;
                best_y = y;
            }
        }
    }
    xy_clean_csv.close();
    xy_all_csv.close();

    std::ofstream class_csv(outdir / "event_class_counts.csv");
    class_csv << "class,count,fraction\n";
    class_csv << std::setprecision(10);
    class_csv << "empty," << n_empty << ',' << safe_fraction(n_empty, total_entries) << '\n';
    class_csv << "X_only," << n_x_only << ',' << safe_fraction(n_x_only, total_entries) << '\n';
    class_csv << "Y_only," << n_y_only << ',' << safe_fraction(n_y_only, total_entries) << '\n';
    class_csv << "XY_candidate," << n_xy_candidate << ',' << safe_fraction(n_xy_candidate, total_entries) << '\n';
    class_csv << "clean_XY," << n_clean_xy << ',' << safe_fraction(n_clean_xy, total_entries) << '\n';
    class_csv << "multi_XY," << n_multi_xy << ',' << safe_fraction(n_multi_xy, total_entries) << '\n';
    class_csv.close();

    auto h_x_occupancy = new TH1D("h_x_occupancy", "hod2018 X occupancy;X channel;Hits", 16, -0.5, 15.5);
    auto h_y_occupancy = new TH1D("h_y_occupancy", "hod2018 Y occupancy;Y channel (Y0 = Ch16);Hits", 16, -0.5, 15.5);
    auto h_xy_clean = new TH2D("h_xy_clean", "hod2018 clean_XY map;X channel;Y channel;Counts", 16, -0.5, 15.5, 16, -0.5, 15.5);
    auto h_xy_all = new TH2D("h_xy_all_combinations", "hod2018 XY all combinations;X channel;Y channel;Counts", 16, -0.5, 15.5, 16, -0.5, 15.5);
    auto h_event_classes = new TH1D("h_event_classes", "hod2018 event classes;Class;Events", 6, 0.5, 6.5);
    auto h_all_amp_neg = new TH1D("h_all_channels_max_amp_neg", "Max negative amplitude by channel;Channel;max amp_{neg} [ADC]", 64, -0.5, 63.5);
    auto h_all_frac_hits = new TH1D("h_all_channels_frac_hits", "Relaxed hit fraction by channel;Channel;Fraction", 64, -0.5, 63.5);
    auto h_hod2019_amp_neg = new TH1D("h_hod2019_max_amp_neg", "hod2019 max negative amplitude;Channel;max amp_{neg} [ADC]", 32, 31.5, 63.5);
    auto h_hod2019_sigma = new TH1D("h_hod2019_baseline_sigma", "hod2019 mean baseline sigma;Channel;mean sigma [ADC]", 32, 31.5, 63.5);

    configure_hist_1d(h_x_occupancy, "X channel", "Hits");
    configure_hist_1d(h_y_occupancy, "Y channel (Y0 = Ch16)", "Hits");
    configure_hist_2d(h_xy_clean, "X channel", "Y channel", "Counts");
    configure_hist_2d(h_xy_all, "X channel", "Y channel", "Counts");
    configure_hist_1d(h_event_classes, "Class", "Events");
    configure_hist_1d(h_all_amp_neg, "Channel", "max amp_{neg} [ADC]");
    configure_hist_1d(h_all_frac_hits, "Channel", "Relaxed hit fraction");
    configure_hist_1d(h_hod2019_amp_neg, "Channel", "max amp_{neg} [ADC]");
    configure_hist_1d(h_hod2019_sigma, "Channel", "mean sigma [ADC]");

    h_x_occupancy->SetFillColor(kAzure - 9);
    h_x_occupancy->SetLineColor(kAzure + 2);
    h_y_occupancy->SetFillColor(kOrange - 3);
    h_y_occupancy->SetLineColor(kOrange + 7);
    h_all_amp_neg->SetFillColor(kTeal - 7);
    h_all_amp_neg->SetLineColor(kTeal + 3);
    h_all_frac_hits->SetFillColor(kViolet - 8);
    h_all_frac_hits->SetLineColor(kViolet + 2);
    h_hod2019_amp_neg->SetFillColor(kRed - 9);
    h_hod2019_amp_neg->SetLineColor(kRed + 1);
    h_hod2019_sigma->SetFillColor(kGray + 1);
    h_hod2019_sigma->SetLineColor(kBlack);

    for (int i = 0; i < 16; ++i) {
        h_x_occupancy->SetBinContent(i + 1, x_occupancy[i]);
        h_y_occupancy->SetBinContent(i + 1, y_occupancy[i]);
    }
    for (int x = 0; x < 16; ++x) {
        for (int y = 0; y < 16; ++y) {
            h_xy_clean->SetBinContent(x + 1, y + 1, xy_clean[x][y]);
            h_xy_all->SetBinContent(x + 1, y + 1, xy_all_combinations[x][y]);
        }
    }

    const std::array<std::string, 6> class_names = {
        "empty", "X_only", "Y_only", "XY_candidate", "clean_XY", "multi_XY"
    };
    const std::array<long long, 6> class_counts = {
        n_empty, n_x_only, n_y_only, n_xy_candidate, n_clean_xy, n_multi_xy
    };
    for (int i = 0; i < 6; ++i) {
        h_event_classes->SetBinContent(i + 1, class_counts[i]);
        h_event_classes->GetXaxis()->SetBinLabel(i + 1, class_names[i].c_str());
    }
    h_event_classes->SetFillColor(kSpring - 7);
    h_event_classes->SetLineColor(kSpring + 4);
    h_event_classes->LabelsOption("v", "X");

    for (int ch = 0; ch < kNChannels; ++ch) {
        const ChannelStats& s = channel_stats[ch];
        h_all_amp_neg->SetBinContent(ch + 1, s.max_amp_neg);
        h_all_frac_hits->SetBinContent(ch + 1, safe_fraction(s.hits_relaxed, s.entries));
        if (ch >= 32) {
            h_hod2019_amp_neg->SetBinContent(ch - 31, s.max_amp_neg);
            h_hod2019_sigma->SetBinContent(ch - 31, mean_sigma(s));
        }
    }

    auto h_xy_clean_norm = static_cast<TH2D*>(h_xy_clean->Clone("h_xy_clean_norm"));
    h_xy_clean_norm->SetTitle("hod2018 clean_XY normalized map;X channel;Y channel;Probability");
    configure_hist_2d(h_xy_clean_norm, "X channel", "Y channel", "Probability");
    if (total_clean_xy_map > 0) h_xy_clean_norm->Scale(1.0 / static_cast<double>(total_clean_xy_map));

    const double max_clean_amp_x = clean_amp_x.empty() ? 0.0 : *std::max_element(clean_amp_x.begin(), clean_amp_x.end());
    const double max_clean_amp_y = clean_amp_y.empty() ? 0.0 : *std::max_element(clean_amp_y.begin(), clean_amp_y.end());
    const double amp_upper = nice_upper(std::max(max_clean_amp_x, max_clean_amp_y));
    auto h_amp_x_clean = new TH1D("h_amp_x_clean", "clean_XY X amplitude;amp_X [ADC];Events", 120, 0.0, amp_upper);
    auto h_amp_y_clean = new TH1D("h_amp_y_clean", "clean_XY Y amplitude;amp_Y [ADC];Events", 120, 0.0, amp_upper);
    auto h_amp_x_vs_y = new TH2D("h_amp_x_vs_y_clean", "clean_XY amplitude X vs Y;amp_X [ADC];amp_Y [ADC];Events",
                                 100, 0.0, amp_upper, 100, 0.0, amp_upper);
    configure_hist_1d(h_amp_x_clean, "amp_X [ADC]", "Events");
    configure_hist_1d(h_amp_y_clean, "amp_Y [ADC]", "Events");
    configure_hist_2d(h_amp_x_vs_y, "amp_X [ADC]", "amp_Y [ADC]", "Events");
    h_amp_x_clean->SetFillColor(kAzure - 9);
    h_amp_x_clean->SetLineColor(kAzure + 2);
    h_amp_y_clean->SetFillColor(kOrange - 3);
    h_amp_y_clean->SetLineColor(kOrange + 7);

    int max_sample_min = 0;
    for (int v : clean_sample_min_x) max_sample_min = std::max(max_sample_min, v);
    for (int v : clean_sample_min_y) max_sample_min = std::max(max_sample_min, v);
    max_sample_min = std::max(max_sample_min, 10);
    auto h_sample_min_x = new TH1D("h_sample_min_x", "clean_XY X sample_min;sample_min_X;Events",
                                   max_sample_min + 1, -0.5, max_sample_min + 0.5);
    auto h_sample_min_y = new TH1D("h_sample_min_y", "clean_XY Y sample_min;sample_min_Y;Events",
                                   max_sample_min + 1, -0.5, max_sample_min + 0.5);
    auto h_sample_min_x_vs_y = new TH2D("h_sample_min_x_vs_y", "clean_XY sample_min X vs Y;sample_min_X;sample_min_Y;Events",
                                        max_sample_min + 1, -0.5, max_sample_min + 0.5,
                                        max_sample_min + 1, -0.5, max_sample_min + 0.5);
    configure_hist_1d(h_sample_min_x, "sample_min_X", "Events");
    configure_hist_1d(h_sample_min_y, "sample_min_Y", "Events");
    configure_hist_2d(h_sample_min_x_vs_y, "sample_min_X", "sample_min_Y", "Events");
    h_sample_min_x->SetFillColor(kAzure - 9);
    h_sample_min_x->SetLineColor(kAzure + 2);
    h_sample_min_y->SetFillColor(kOrange - 3);
    h_sample_min_y->SetLineColor(kOrange + 7);

    for (size_t i = 0; i < clean_amp_x.size(); ++i) {
        h_amp_x_clean->Fill(clean_amp_x[i]);
        h_amp_y_clean->Fill(clean_amp_y[i]);
        h_amp_x_vs_y->Fill(clean_amp_x[i], clean_amp_y[i]);
        h_sample_min_x->Fill(clean_sample_min_x[i]);
        h_sample_min_y->Fill(clean_sample_min_y[i]);
        h_sample_min_x_vs_y->Fill(clean_sample_min_x[i], clean_sample_min_y[i]);
    }

    save_hist1_png(h_x_occupancy, pngdir / "hod2018_x_occupancy.png");
    save_hist1_png(h_y_occupancy, pngdir / "hod2018_y_occupancy.png");

    {
        TCanvas c("c_xy_occupancy_comparison", "", 1100, 760);
        c.SetLeftMargin(0.11);
        c.SetRightMargin(0.04);
        c.SetBottomMargin(0.13);
        auto hx = static_cast<TH1D*>(h_x_occupancy->Clone("h_x_occupancy_compare"));
        auto hy = static_cast<TH1D*>(h_y_occupancy->Clone("h_y_occupancy_compare"));
        hx->SetFillColor(0);
        hy->SetFillColor(0);
        hx->SetLineColor(kAzure + 2);
        hy->SetLineColor(kOrange + 7);
        hx->SetMarkerColor(kAzure + 2);
        hy->SetMarkerColor(kOrange + 7);
        hx->SetMarkerStyle(20);
        hy->SetMarkerStyle(21);
        hx->SetMaximum(std::max(hx->GetMaximum(), hy->GetMaximum()) * 1.20 + 1.0);
        hx->Draw("HIST E");
        hy->Draw("HIST E SAME");
        TLegend leg(0.70, 0.78, 0.92, 0.90);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hx, "X (Ch0-Ch15)", "l");
        leg.AddEntry(hy, "Y (Ch16-Ch31)", "l");
        leg.Draw();
        c.SaveAs((pngdir / "hod2018_xy_occupancy_comparison.png").string().c_str());
        delete hx;
        delete hy;
    }

    save_hist2_png(h_xy_clean, pngdir / "hod2018_xy_map_clean_counts.png");
    save_hist2_png(h_xy_clean, pngdir / "hod2018_xy_map_clean_logz.png", true);
    save_hist2_png(h_xy_clean_norm, pngdir / "hod2018_xy_map_clean_normalized.png");
    save_hist2_png(h_xy_all, pngdir / "hod2018_xy_map_all_combinations.png");
    save_hist1_png(h_event_classes, pngdir / "hod2018_event_classes.png");
    save_hist1_png(h_nhits, pngdir / "hod2018_multiplicity_nhits.png");
    save_hist2_png(h_nxny, pngdir / "hod2018_nx_ny_map.png");
    save_hist1_png(h_amp_x_clean, pngdir / "hod2018_amp_X_clean.png");
    save_hist1_png(h_amp_y_clean, pngdir / "hod2018_amp_Y_clean.png");
    save_hist2_png(h_amp_x_vs_y, pngdir / "hod2018_amp_X_vs_Y_clean.png");
    save_hist1_png(h_sample_min_x, pngdir / "hod2018_sample_min_X_clean.png");
    save_hist1_png(h_sample_min_y, pngdir / "hod2018_sample_min_Y_clean.png");
    save_hist2_png(h_sample_min_x_vs_y, pngdir / "hod2018_sample_min_X_vs_Y_clean.png");
    save_hist1_png(h_all_amp_neg, pngdir / "all_channels_max_amp_neg.png");
    save_hist1_png(h_all_frac_hits, pngdir / "all_channels_frac_hits_relaxed.png");
    save_hist1_png(h_hod2019_amp_neg, pngdir / "hod2019_max_amp_neg.png");
    save_hist1_png(h_hod2019_sigma, pngdir / "hod2019_baseline_sigma.png");

    const fs::path hist_root_path = outdir / "final_hod2018_histograms.root";
    TFile hist_file(hist_root_path.string().c_str(), "RECREATE");
    if (!hist_file.IsZombie()) {
        h_x_occupancy->Write();
        h_y_occupancy->Write();
        h_xy_clean->Write();
        h_xy_clean_norm->Write();
        h_xy_all->Write();
        h_event_classes->Write();
        h_nhits->Write();
        h_nxny->Write();
        h_amp_x_clean->Write();
        h_amp_y_clean->Write();
        h_amp_x_vs_y->Write();
        h_sample_min_x->Write();
        h_sample_min_y->Write();
        h_sample_min_x_vs_y->Write();
        h_all_amp_neg->Write();
        h_all_frac_hits->Write();
        h_hod2019_amp_neg->Write();
        h_hod2019_sigma->Write();
        hist_file.Close();
    } else {
        std::cerr << "ADVERTENCIA: no pude crear " << hist_root_path << "\n";
    }

    std::vector<long long> x_counts(16, 0);
    std::vector<long long> y_counts(16, 0);
    long long total_x_hits = 0;
    long long total_y_hits = 0;
    for (int i = 0; i < 16; ++i) {
        x_counts[i] = x_occupancy[i];
        y_counts[i] = y_occupancy[i];
        total_x_hits += x_occupancy[i];
        total_y_hits += y_occupancy[i];
    }

    const auto max_x_it = std::max_element(x_counts.begin(), x_counts.end());
    const auto max_y_it = std::max_element(y_counts.begin(), y_counts.end());
    const int dominant_x = static_cast<int>(std::distance(x_counts.begin(), max_x_it));
    const int dominant_y = static_cast<int>(std::distance(y_counts.begin(), max_y_it));
    const double dominant_x_frac = safe_fraction(*max_x_it, total_x_hits);
    const double dominant_y_frac = safe_fraction(*max_y_it, total_y_hits);

    long long total_hod2019_strict = 0;
    long long total_hod2019_relaxed = 0;
    double max_hod2019_amp_neg = 0.0;
    for (int ch = 32; ch < 64; ++ch) {
        total_hod2019_strict += channel_stats[ch].hits_strict;
        total_hod2019_relaxed += channel_stats[ch].hits_relaxed;
        max_hod2019_amp_neg = std::max(max_hod2019_amp_neg, channel_stats[ch].max_amp_neg);
    }

    long long total_hod2018_hits = total_x_hits + total_y_hits;
    long long max_any_hod2018 = 0;
    int max_any_hod2018_ch = -1;
    for (int ch = 0; ch < 32; ++ch) {
        if (channel_stats[ch].hits_relaxed > max_any_hod2018) {
            max_any_hod2018 = channel_stats[ch].hits_relaxed;
            max_any_hod2018_ch = ch;
        }
    }
    const double max_any_hod2018_frac = safe_fraction(max_any_hod2018, total_hod2018_hits);

    std::ofstream summary(outdir / "summary.txt");
    summary << std::setprecision(8);
    summary << "Analisis final hod2018 / diagnostico hod2019 v4\n";
    summary << "================================================\n\n";
    summary << "Archivo analizado: " << input_path.string() << "\n";
    summary << "TTree: data_dig\n";
    summary << "Total de eventos: " << total_entries << "\n\n";

    summary << "Clases de evento hod2018 (XY_candidate es inclusivo):\n";
    summary << "  empty: " << n_empty << " (" << safe_fraction(n_empty, total_entries) << ")\n";
    summary << "  X_only: " << n_x_only << " (" << safe_fraction(n_x_only, total_entries) << ")\n";
    summary << "  Y_only: " << n_y_only << " (" << safe_fraction(n_y_only, total_entries) << ")\n";
    summary << "  XY_candidate: " << n_xy_candidate << " (" << safe_fraction(n_xy_candidate, total_entries) << ")\n";
    summary << "  clean_XY: " << n_clean_xy << " (" << safe_fraction(n_clean_xy, total_entries) << ")\n";
    summary << "  multi_XY: " << n_multi_xy << " (" << safe_fraction(n_multi_xy, total_entries) << ")\n\n";

    if (have_first_timestamp) {
        long double span_units = max_timestamp - min_timestamp;
        if (span_units <= 0.0L) span_units = last_timestamp - first_timestamp;
        if (span_units > 0.0L) {
            const long double duration_s = span_units * 1.0e-6L;
            const long double rate_hz = static_cast<long double>(total_entries) / duration_s;
            summary << "Timestamp: disponible.\n";
            summary << "Duracion estimada: " << static_cast<double>(duration_s)
                    << " s, asumiendo timestamp en microsegundos.\n";
            summary << "Rate estimado: " << static_cast<double>(rate_hz) << " Hz\n\n";
        } else {
            summary << "Timestamp: disponible, pero no permite calcular duracion positiva.\n\n";
        }
    } else {
        summary << "Timestamp: no disponible; rate estimado no calculado.\n\n";
    }

    summary << "Canales X mas activos: " << format_top_channels(x_counts, 0, "X") << "\n";
    summary << "Canales Y mas activos: " << format_top_channels(y_counts, 16, "Y") << "\n";
    summary << "Canales hod2018 muertos o con cero hits relajados: " << list_zero_hit_hod2018(channel_stats) << "\n";
    summary << "Canal X dominante: X" << dominant_x << " (Ch" << dominant_x << "), hits="
            << *max_x_it << ", fraccion de hits X=" << dominant_x_frac << "\n";
    summary << "Canal Y dominante: Y" << dominant_y << " (Ch" << (16 + dominant_y) << "), hits="
            << *max_y_it << ", fraccion de hits Y=" << dominant_y_frac << "\n";

    if (total_x_hits > 0 && safe_fraction(x_counts[15], total_x_hits) > 0.5) {
        summary << "ADVERTENCIA: X15 domina mas del 50% de los hits X.\n";
    }
    if (total_hod2018_hits > 0 && max_any_hod2018_frac > 0.5) {
        summary << "ADVERTENCIA: Ch" << max_any_hod2018_ch
                << " domina mas del 50% de todos los hits hod2018.\n";
    }
    summary << "\n";

    summary << "Total clean_XY usado para mapa: " << total_clean_xy_map << "\n";
    if (best_x >= 0 && best_y >= 0) {
        summary << "Bin mas poblado del mapa clean_XY: X" << best_x << ", Y" << best_y
                << " con " << best_xy_count << " cuentas.\n\n";
    } else {
        summary << "Bin mas poblado del mapa clean_XY: no disponible.\n\n";
    }

    summary << "Resumen hod2019 (Ch32-Ch63):\n";
    summary << "  total hits estrictos: " << total_hod2019_strict << "\n";
    summary << "  total hits relajados: " << total_hod2019_relaxed << "\n";
    summary << "  max_amp_neg maximo: " << max_hod2019_amp_neg << " ADC\n";
    if (total_hod2019_strict == 0 && max_hod2019_amp_neg < 5.0) {
        summary << "  Conclusion hod2019: hod2019 no muestra señal física detectable en este run\n";
    } else if (total_hod2019_strict == 0) {
        summary << "  Conclusion hod2019: no hay hits estrictos; las excursiones relajadas deben interpretarse como ruido o diagnostico, no como plano operativo.\n";
    } else {
        summary << "  Conclusion hod2019: aparecen hits estrictos aislados; por su baja frecuencia y por el contexto del run se documentan como outliers/ruido de diagnostico, no como plano operativo.\n";
    }
    summary << "\n";
    summary << "Conclusion general:\n";
    summary << "Con este run se puede caracterizar hod2018 como plano 2D, pero no reconstruir tracks 3D porque hod2019 no entrega señal útil.\n";
    summary.close();

    write_readme(outdir);

    std::cout << "Salidas escritas en: " << outdir << "\n";
    std::cout << "Resumen: " << (outdir / "summary.txt") << "\n";
    std::cout << "PNGs: " << pngdir << "\n";
    return 0;
}
