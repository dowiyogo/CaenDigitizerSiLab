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
#include <sstream>
#include <string>
#include <vector>
#include <map>

struct ChannelStats {
    int ch = -1;
    char plane = 'X'; // 'X' or 'Y'
    long long entries = 0;
    long long hits = 0;
    double sum_baseline = 0.0;
    double sum_sigma = 0.0;
    double max_amp_neg = 0.0;
    double min_vmin = 1e99;
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

static std::string join_channels(const std::vector<int>& v) {
    if (v.empty()) return "";
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ";";
        oss << v[i];
    }
    return oss.str();
}

int main(int argc, char** argv) {
    const int NBASE = 50;
    const int NCH_TOTAL = 64;
    const int HOD2018_START = 0;
    const int HOD2018_END = 31; // inclusive

    if (argc < 2) {
        std::cerr << "Uso:\n  " << argv[0] << " archivo.root [max_events]\n";
        return 1;
    }
    std::string filename = argv[1];
    long long max_events = -1;
    if (argc >= 3) max_events = std::atoll(argv[2]);

    std::cout << "Archivo: " << filename << "\n";
    std::cout << "max_events: " << max_events << "\n";

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

    // Branch pointers
    std::vector<float>* ch[NCH_TOTAL] = {nullptr};
    int event = 0;

    // Robust timestamp handling (types may vary)
    enum class TimestampKind { None, Long, Long64, ULong64 };
    TimestampKind timestamp_kind = TimestampKind::None;
    Long_t timestamp_long = 0;
    Long64_t timestamp_long64 = 0;
    ULong64_t timestamp_ulong64 = 0;

    // Set branch addresses for Ch0..Ch31
    for (int i = 0; i <= HOD2018_END; ++i) {
        std::string br = "Ch" + std::to_string(i);
        TBranch* b = t->GetBranch(br.c_str());
        if (!b) {
            std::cerr << "ADVERTENCIA: no existe rama " << br << "\n";
            continue;
        }
        t->SetBranchAddress(br.c_str(), &ch[i]);
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

    long long nentries = t->GetEntries();
    if (max_events > 0 && max_events < nentries) nentries = max_events;

    // Prepare outputs
    std::ofstream evcsv("hod2018_event_summary.csv");
    evcsv << "entry,event,timestamp,nX,nY,nHits,class,x_channels,y_channels,max_amp_X,max_amp_Y\n";

    std::ofstream chcsv("hod2018_channel_summary.csv");
    chcsv << "channel,plane,entries,hits,frac_hits,mean_baseline,mean_sigma,max_amp_neg,min_vmin\n";

    std::ofstream xycsv("hod2018_xy_map.csv");
    xycsv << "x_channel,y_channel,count\n";

    // stats for channels 0..31
    std::vector<ChannelStats> stats(32);
    for (int i = 0; i <= 31; ++i) {
        stats[i].ch = i;
        stats[i].plane = (i <= 15) ? 'X' : 'Y';
    }

    // xy map counts 16x16
    std::vector<std::vector<long long>> xy_map(16, std::vector<long long>(16, 0));

    // event class counters
    long long n_empty = 0, n_X_only = 0, n_Y_only = 0, n_XY_candidate = 0, n_clean_XY = 0, n_multi_XY = 0;

    for (long long ie = 0; ie < nentries; ++ie) {
        t->GetEntry(ie);
        if (ie % 1000 == 0) {
            std::cout << "Procesando evento " << ie << " / " << nentries << "\r" << std::flush;
        }

        unsigned long long timestamp_out = 0;
        if (timestamp_kind == TimestampKind::ULong64) timestamp_out = static_cast<unsigned long long>(timestamp_ulong64);
        else if (timestamp_kind == TimestampKind::Long64) timestamp_out = static_cast<unsigned long long>(timestamp_long64);
        else if (timestamp_kind == TimestampKind::Long) timestamp_out = static_cast<unsigned long long>(timestamp_long);

        std::vector<int> X_chs; X_chs.reserve(16);
        std::vector<int> Y_chs; Y_chs.reserve(16);
        double max_amp_X = 0.0;
        double max_amp_Y = 0.0;

        for (int ic = HOD2018_START; ic <= HOD2018_END; ++ic) {
            if (!ch[ic]) continue;
            if (ch[ic]->empty()) continue;
            const auto& v = *ch[ic];
            double baseline = mean_first_n(v, NBASE);
            double sigma = sigma_first_n(v, NBASE, baseline);
            auto mm = std::minmax_element(v.begin(), v.end());
            double vmin = *mm.first;
            double amp_neg = baseline - vmin;

            double thr = std::max(4.0 * sigma, 5.0);
            // require at least 2 consecutive samples below baseline - thr
            bool consec_ok = false;
            for (size_t k = 0; k + 1 < v.size(); ++k) {
                if (v[k] < baseline - thr && v[k+1] < baseline - thr) { consec_ok = true; break; }
            }

            bool is_hit = (amp_neg > thr) && consec_ok;

            // update per-channel stats (only consider events where waveform existed)
            if (ic <= 31) {
                ChannelStats& s = stats[ic];
                s.entries++;
                s.sum_baseline += baseline;
                s.sum_sigma += sigma;
                s.max_amp_neg = std::max(s.max_amp_neg, amp_neg);
                s.min_vmin = std::min(s.min_vmin, vmin);
                if (is_hit) {
                    s.hits++;
                    if (ic <= 15) {
                        X_chs.push_back(ic);
                        max_amp_X = std::max(max_amp_X, amp_neg);
                    } else {
                        Y_chs.push_back(ic - 16);
                        max_amp_Y = std::max(max_amp_Y, amp_neg);
                    }
                }
            }
        }

        int nX = X_chs.size();
        int nY = Y_chs.size();
        int nHits = nX + nY;
        std::string cls;
        if (nX == 0 && nY == 0) { cls = "empty"; n_empty++; }
        else if (nX > 0 && nY == 0) { cls = "X_only"; n_X_only++; }
        else if (nX == 0 && nY > 0) { cls = "Y_only"; n_Y_only++; }
        else { // both >0
            n_XY_candidate++;
            if (nX == 1 && nY == 1) { cls = "clean_XY"; n_clean_XY++; }
            else if (nX + nY > 2) { cls = "multi_XY"; n_multi_XY++; }
            else { cls = "XY_candidate"; }
        }

        // fill xy map only for clean_XY (nX==1 && nY==1)
        if (cls == "clean_XY") {
            int x = X_chs[0];
            int y = Y_chs[0];
            int x_idx = x; // 0..15
            int y_idx = y; // already 0..15 (we stored ic-16)
            if (x_idx >=0 && x_idx < 16 && y_idx >= 0 && y_idx < 16) xy_map[x_idx][y_idx]++;
        }

        // write event summary
        evcsv << ie << "," << event << "," << timestamp_out << ","
              << nX << "," << nY << "," << nHits << "," << cls << ","
              << join_channels(X_chs) << "," << join_channels(Y_chs) << ","
              << max_amp_X << "," << max_amp_Y << "\n";
    }

    std::cout << "\nProcesamiento terminado.\n";

    // write channel summary
    for (int i = 0; i <= 31; ++i) {
        ChannelStats& s = stats[i];
        double frac = (s.entries > 0) ? double(s.hits) / s.entries : 0.0;
        double mean_bas = (s.entries > 0) ? s.sum_baseline / s.entries : 0.0;
        double mean_sig = (s.entries > 0) ? s.sum_sigma / s.entries : 0.0;
        chcsv << s.ch << "," << s.plane << "," << s.entries << "," << s.hits << ","
              << frac << "," << mean_bas << "," << mean_sig << ","
              << s.max_amp_neg << "," << ((s.entries>0)?s.min_vmin:0.0) << "\n";
    }

    // write xy map CSV
    for (int x = 0; x < 16; ++x) {
        for (int y = 0; y < 16; ++y) {
            xycsv << x << "," << y << "," << xy_map[x][y] << "\n";
        }
    }

    evcsv.close(); chcsv.close(); xycsv.close();

    long long total = nentries;
    std::cout << "Total eventos procesados: " << total << "\n";
    std::cout << "empty: " << n_empty << "\n";
    std::cout << "X_only: " << n_X_only << "\n";
    std::cout << "Y_only: " << n_Y_only << "\n";
    std::cout << "XY_candidate: " << n_XY_candidate << "\n";
    std::cout << "clean_XY: " << n_clean_XY << "\n";
    std::cout << "multi_XY: " << n_multi_XY << "\n";

    // lists of channels X and Y ordered by hits
    std::vector<ChannelStats> sortedX, sortedY;
    for (int i = 0; i <= 15; ++i) sortedX.push_back(stats[i]);
    for (int i = 16; i <= 31; ++i) sortedY.push_back(stats[i]);
    auto cmp = [](const ChannelStats& a, const ChannelStats& b){ return a.hits > b.hits; };
    std::sort(sortedX.begin(), sortedX.end(), cmp);
    std::sort(sortedY.begin(), sortedY.end(), cmp);

    std::cout << "\nCanales X ordenados por hits:\n";
    for (auto &s: sortedX) std::cout << "Ch" << s.ch << ": " << s.hits << "\n";
    std::cout << "\nCanales Y ordenados por hits:\n";
    for (auto &s: sortedY) std::cout << "Ch" << s.ch << ": " << s.hits << "\n";

    std::cout << "\nFracción clean_XY / total: ";
    std::cout << (double(n_clean_XY) / double(std::max<long long>(1,total))) << "\n";
    std::cout << "Fracción XY_candidate / total: ";
    std::cout << (double(n_XY_candidate) / double(std::max<long long>(1,total))) << "\n";

    std::cout << "\nArchivos generados:\n"
              << "  hod2018_event_summary.csv\n"
              << "  hod2018_channel_summary.csv\n"
              << "  hod2018_xy_map.csv\n";

    f->Close();
    return 0;
}
