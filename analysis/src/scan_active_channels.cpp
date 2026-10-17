#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

struct ChannelStats {
    int ch = -1;
    bool exists = false;

    long long entries = 0;
    long long hits_neg = 0;
    long long hits_pos = 0;

    double sum_baseline = 0.0;
    double sum_sigma = 0.0;

    double max_amp_neg = 0.0;
    double max_amp_pos = 0.0;

    double min_vmin = 1e99;
    double max_vmax = -1e99;
};

static double mean_first_n(const std::vector<float>& v, int n) {
    if (v.empty()) return 0.0;
    n = std::min<int>(n, v.size());

    double s = 0.0;
    for (int i = 0; i < n; ++i) s += v[i];

    return s / n;
}

static double sigma_first_n(const std::vector<float>& v, int n, double mu) {
    if (v.size() < 2) return 0.0;
    n = std::min<int>(n, v.size());

    if (n < 2) return 0.0;

    double s2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = v[i] - mu;
        s2 += d * d;
    }

    return std::sqrt(s2 / (n - 1));
}

int main(int argc, char** argv) {
    const int NCH = 64;
    const int NBASE = 50;

    if (argc < 2) {
        std::cerr << "Uso:\n"
                  << "  " << argv[0]
                  << " archivo.root [max_events] [nsigma] [abs_threshold_adc]\n\n"
                  << "Ejemplo:\n"
                  << "  " << argv[0]
                  << " lectura_vx1740_ext_nim_1M.root 20000 5 20\n";
        return 1;
    }

    std::string filename = argv[1];

    long long max_events = -1;
    if (argc >= 3) max_events = std::atoll(argv[2]);

    double nsigma = 5.0;
    if (argc >= 4) nsigma = std::atof(argv[3]);

    double abs_threshold_adc = 20.0;
    if (argc >= 5) abs_threshold_adc = std::atof(argv[4]);

    std::cout << "Archivo: " << filename << "\n";
    std::cout << "max_events: " << max_events << "\n";
    std::cout << "nsigma: " << nsigma << "\n";
    std::cout << "abs_threshold_adc: " << abs_threshold_adc << "\n";

    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: no pude abrir " << filename << "\n";
        return 1;
    }

    TTree* t = dynamic_cast<TTree*>(f->Get("data_dig"));
    if (!t) {
        std::cerr << "ERROR: no encontré el TTree data_dig\n";
        f->Close();
        return 1;
    }

    std::cout << "Entradas en data_dig: " << t->GetEntries() << "\n";

    std::vector<float>* ch[NCH] = {nullptr};

    int event = 0;
    enum class TimestampKind {
        None,
        Long,
        Long64,
        ULong64
    };
    TimestampKind timestamp_kind = TimestampKind::None;
    Long_t timestamp_long = 0;
    Long64_t timestamp_long64 = 0;
    ULong64_t timestamp_ulong64 = 0;
    std::vector<float>* time = nullptr;

    std::vector<ChannelStats> stats(NCH);
    for (int i = 0; i < NCH; ++i) {
        stats[i].ch = i;

        std::string brname = "Ch" + std::to_string(i);
        TBranch* br = t->GetBranch(brname.c_str());

        if (!br) {
            std::cerr << "ADVERTENCIA: no existe rama " << brname << "\n";
            stats[i].exists = false;
            continue;
        }

        stats[i].exists = true;
        t->SetBranchAddress(brname.c_str(), &ch[i]);
    }

    if (t->GetBranch("event")) {
        t->SetBranchAddress("event", &event);
    } else {
        std::cerr << "ADVERTENCIA: no existe rama event\n";
    }

    if (t->GetBranch("timestamp")) {
        TBranch* ts_branch = t->GetBranch("timestamp");
        TLeaf* ts_leaf = ts_branch ? ts_branch->GetLeaf("timestamp") : nullptr;
        std::string ts_type = ts_leaf ? ts_leaf->GetTypeName() : "";

        if (ts_type == "ULong64_t" || ts_type == "unsigned long long") {
            t->SetBranchAddress("timestamp", &timestamp_ulong64);
            timestamp_kind = TimestampKind::ULong64;
        } else if (ts_type == "Long64_t" || ts_type == "long long") {
            t->SetBranchAddress("timestamp", &timestamp_long64);
            timestamp_kind = TimestampKind::Long64;
        } else {
            t->SetBranchAddress("timestamp", &timestamp_long);
            timestamp_kind = TimestampKind::Long;
            if (!ts_type.empty()) {
                std::cerr << "ADVERTENCIA: rama timestamp tipo " << ts_type
                          << ", usando Long_t compatible\n";
            }
        }
    } else {
        std::cerr << "ADVERTENCIA: no existe rama timestamp\n";
    }

    if (t->GetBranch("time")) {
        t->SetBranchAddress("time", &time);
    } else {
        std::cerr << "ADVERTENCIA: no existe rama time\n";
    }

    long long nentries = t->GetEntries();
    if (max_events > 0 && max_events < nentries) nentries = max_events;

    std::ofstream mult("event_multiplicity.csv");
    mult << "entry,event,timestamp,n_hit_neg,n_hit_pos\n";

    for (long long ie = 0; ie < nentries; ++ie) {
        t->GetEntry(ie);

        unsigned long long timestamp_out = 0;
        if (timestamp_kind == TimestampKind::ULong64) {
            timestamp_out = static_cast<unsigned long long>(timestamp_ulong64);
        } else if (timestamp_kind == TimestampKind::Long64) {
            timestamp_out = static_cast<unsigned long long>(timestamp_long64);
        } else if (timestamp_kind == TimestampKind::Long) {
            timestamp_out = static_cast<unsigned long long>(timestamp_long);
        }

        if (ie % 1000 == 0) {
            std::cout << "Procesando evento " << ie << " / " << nentries << "\r" << std::flush;
        }

        int n_hit_neg_event = 0;
        int n_hit_pos_event = 0;

        for (int ic = 0; ic < NCH; ++ic) {
            if (!stats[ic].exists) continue;
            if (!ch[ic]) continue;
            if (ch[ic]->empty()) continue;

            const auto& v = *ch[ic];

            double baseline = mean_first_n(v, NBASE);
            double sigma = sigma_first_n(v, NBASE, baseline);

            auto mm = std::minmax_element(v.begin(), v.end());
            double vmin = *mm.first;
            double vmax = *mm.second;

            double amp_neg = baseline - vmin;
            double amp_pos = vmax - baseline;

            bool hit_neg = (amp_neg > nsigma * sigma) && (amp_neg > abs_threshold_adc);
            bool hit_pos = (amp_pos > nsigma * sigma) && (amp_pos > abs_threshold_adc);

            stats[ic].entries++;
            stats[ic].sum_baseline += baseline;
            stats[ic].sum_sigma += sigma;

            stats[ic].max_amp_neg = std::max(stats[ic].max_amp_neg, amp_neg);
            stats[ic].max_amp_pos = std::max(stats[ic].max_amp_pos, amp_pos);

            stats[ic].min_vmin = std::min(stats[ic].min_vmin, vmin);
            stats[ic].max_vmax = std::max(stats[ic].max_vmax, vmax);

            if (hit_neg) {
                stats[ic].hits_neg++;
                n_hit_neg_event++;
            }

            if (hit_pos) {
                stats[ic].hits_pos++;
                n_hit_pos_event++;
            }
        }

        mult << ie << ","
             << event << ","
               << timestamp_out << ","
             << n_hit_neg_event << ","
             << n_hit_pos_event << "\n";
    }

    std::cout << "\nProcesamiento terminado.\n";
    mult.close();

    std::vector<ChannelStats> sorted = stats;
    std::sort(sorted.begin(), sorted.end(),
              [](const ChannelStats& a, const ChannelStats& b) {
                  return a.hits_neg > b.hits_neg;
              });

    std::ofstream csv("active_channels_summary.csv");
    csv << "channel,entries,hits_neg,frac_hits_neg,hits_pos,frac_hits_pos,"
        << "mean_baseline,mean_sigma,max_amp_neg,max_amp_pos,min_vmin,max_vmax\n";

    std::cout << "\nResumen por canal, ordenado por hits negativos:\n";
    std::cout << std::setw(8) << "channel"
              << std::setw(10) << "entries"
              << std::setw(12) << "hits_neg"
              << std::setw(14) << "frac_hits_neg"
              << std::setw(12) << "hits_pos"
              << std::setw(14) << "frac_hits_pos"
              << std::setw(15) << "mean_baseline"
              << std::setw(12) << "mean_sigma"
              << std::setw(12) << "max_amp_neg"
              << std::setw(12) << "max_amp_pos"
              << std::setw(12) << "min_vmin"
              << std::setw(12) << "max_vmax"
              << "\n";

    std::vector<int> channels_probablemente_activos;
    std::vector<int> channels_sin_senal;
    std::vector<int> channels_ruidosos_o_saturados;

    for (const auto& s : sorted) {
        double frac_neg = (s.entries > 0) ? double(s.hits_neg) / s.entries : 0.0;
        double frac_pos = (s.entries > 0) ? double(s.hits_pos) / s.entries : 0.0;
        double mean_baseline = (s.entries > 0) ? s.sum_baseline / s.entries : 0.0;
        double mean_sigma = (s.entries > 0) ? s.sum_sigma / s.entries : 0.0;

        std::cout << std::setw(8) << s.ch
                  << std::setw(10) << s.entries
                  << std::setw(12) << s.hits_neg
                  << std::setw(14) << std::fixed << std::setprecision(5) << frac_neg
                  << std::setw(12) << s.hits_pos
                  << std::setw(14) << std::fixed << std::setprecision(5) << frac_pos
                  << std::setw(15) << std::fixed << std::setprecision(2) << mean_baseline
                  << std::setw(12) << std::fixed << std::setprecision(2) << mean_sigma
                  << std::setw(12) << std::fixed << std::setprecision(2) << s.max_amp_neg
                  << std::setw(12) << std::fixed << std::setprecision(2) << s.max_amp_pos
                  << std::setw(12) << std::fixed << std::setprecision(2)
                  << ((s.entries > 0) ? s.min_vmin : 0.0)
                  << std::setw(12) << std::fixed << std::setprecision(2)
                  << ((s.entries > 0) ? s.max_vmax : 0.0)
                  << "\n";

        csv << s.ch << ","
            << s.entries << ","
            << s.hits_neg << ","
            << frac_neg << ","
            << s.hits_pos << ","
            << frac_pos << ","
            << mean_baseline << ","
            << mean_sigma << ","
            << s.max_amp_neg << ","
            << s.max_amp_pos << ","
            << ((s.entries > 0) ? s.min_vmin : 0.0) << ","
            << ((s.entries > 0) ? s.max_vmax : 0.0) << "\n";

        if (s.exists) {
            if (s.hits_neg == 0) {
                channels_sin_senal.push_back(s.ch);
            }
            if (frac_neg > 0.001) {
                channels_probablemente_activos.push_back(s.ch);
            }
            if (frac_neg > 0.5) {
                channels_ruidosos_o_saturados.push_back(s.ch);
            }
        }
    }

    csv.close();

    std::cout << "\nArchivos generados:\n";
    std::cout << "  active_channels_summary.csv\n";
    std::cout << "  event_multiplicity.csv\n";

    std::cout << "\nInterpretación rápida:\n";
    std::cout << "  - 'probablemente_activo': frac_hits_neg > 0.001\n";
    std::cout << "  - 'sin_senal_detectada': ningún cruce significativo\n";
    std::cout << "  - 'posiblemente_ruidoso_o_saturado': frac_hits_neg > 0.5\n";
    std::cout << "\nCanales probablemente activos: ";
    if (channels_probablemente_activos.empty()) {
        std::cout << "(ninguno)";
    } else {
        for (size_t i = 0; i < channels_probablemente_activos.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << channels_probablemente_activos[i];
        }
    }
    std::cout << "\n";

    std::cout << "Canales sin señal detectada: ";
    if (channels_sin_senal.empty()) {
        std::cout << "(ninguno)";
    } else {
        for (size_t i = 0; i < channels_sin_senal.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << channels_sin_senal[i];
        }
    }
    std::cout << "\n";

    std::cout << "Canales posiblemente ruidosos o saturados: ";
    if (channels_ruidosos_o_saturados.empty()) {
        std::cout << "(ninguno)";
    } else {
        for (size_t i = 0; i < channels_ruidosos_o_saturados.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << channels_ruidosos_o_saturados[i];
        }
    }
    std::cout << "\n";
    std::cout << "\nSi tus pulsos fueran positivos, mira hits_pos y frac_hits_pos.\n";

    f->Close();
    return 0;
}