#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct ChannelRow {
    int channel = -1;
    std::string plane;
    long long entries = 0;
    long long hits = 0;
    double frac_hits = 0.0;
    double mean_baseline = 0.0;
    double mean_sigma = 0.0;
    double max_amp_neg = 0.0;
    double min_vmin = 0.0;
};

struct EventRow {
    long long entry = 0;
    long long event = 0;
    unsigned long long timestamp = 0;
    int nX = 0;
    int nY = 0;
    int nHits = 0;
    std::string cls;
    std::string x_channels;
    std::string y_channels;
    double max_amp_X = 0.0;
    double max_amp_Y = 0.0;
};

static std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::stringstream ss(line);
    while (std::getline(ss, field, ',')) {
        fields.push_back(trim(field));
    }
    if (!line.empty() && line.back() == ',') {
        fields.emplace_back("");
    }
    return fields;
}

static bool is_empty_or_header(const std::string& line) {
    std::string s = trim(line);
    if (s.empty()) return true;
    if (!s.empty() && (s[0] == '#')) return true;
    return false;
}

static bool parse_long_long(const std::string& s, long long& out) {
    try {
        size_t idx = 0;
        out = std::stoll(trim(s), &idx);
        return idx > 0;
    } catch (...) {
        return false;
    }
}

static bool parse_ulong_long(const std::string& s, unsigned long long& out) {
    try {
        size_t idx = 0;
        out = std::stoull(trim(s), &idx);
        return idx > 0;
    } catch (...) {
        return false;
    }
}

static bool parse_int(const std::string& s, int& out) {
    long long tmp = 0;
    if (!parse_long_long(s, tmp)) return false;
    out = static_cast<int>(tmp);
    return true;
}

static bool parse_double(const std::string& s, double& out) {
    try {
        size_t idx = 0;
        out = std::stod(trim(s), &idx);
        return idx > 0;
    } catch (...) {
        return false;
    }
}

static bool load_channel_summary(const fs::path& path, std::vector<ChannelRow>& rows, std::string& err) {
    std::ifstream in(path);
    if (!in) {
        err = "no pude abrir " + path.string();
        return false;
    }

    std::string line;
    bool header_seen = false;
    while (std::getline(in, line)) {
        if (is_empty_or_header(line)) continue;
        auto fields = split_csv_line(line);
        if (fields.size() < 9) continue;

        ChannelRow row;
        if (!parse_int(fields[0], row.channel)) continue;
        row.plane = fields[1];
        parse_long_long(fields[2], row.entries);
        parse_long_long(fields[3], row.hits);
        parse_double(fields[4], row.frac_hits);
        parse_double(fields[5], row.mean_baseline);
        parse_double(fields[6], row.mean_sigma);
        parse_double(fields[7], row.max_amp_neg);
        parse_double(fields[8], row.min_vmin);
        rows.push_back(row);
        header_seen = true;
    }

    if (!header_seen || rows.empty()) {
        err = "no encontré filas válidas en " + path.string();
        return false;
    }
    return true;
}

static bool load_event_summary(const fs::path& path, std::vector<EventRow>& rows, std::map<std::string, long long>& class_counts, std::string& err) {
    std::ifstream in(path);
    if (!in) {
        err = "no pude abrir " + path.string();
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (is_empty_or_header(line)) continue;
        auto fields = split_csv_line(line);
        if (fields.size() < 11) continue;

        EventRow row;
        if (!parse_long_long(fields[0], row.entry)) continue;
        parse_long_long(fields[1], row.event);

        unsigned long long ts = 0;
        if (parse_ulong_long(fields[2], ts)) row.timestamp = ts;

        parse_int(fields[3], row.nX);
        parse_int(fields[4], row.nY);
        parse_int(fields[5], row.nHits);
        row.cls = fields[6];
        row.x_channels = fields[7];
        row.y_channels = fields[8];
        parse_double(fields[9], row.max_amp_X);
        parse_double(fields[10], row.max_amp_Y);

        rows.push_back(row);
        class_counts[row.cls]++;
    }

    if (rows.empty()) {
        err = "no encontré filas válidas en " + path.string();
        return false;
    }
    return true;
}

static bool load_xy_map(const fs::path& path, std::vector<std::vector<long long>>& counts, long long& total_clean_xy, std::string& err) {
    std::ifstream in(path);
    if (!in) {
        err = "no pude abrir " + path.string();
        return false;
    }

    counts.assign(16, std::vector<long long>(16, 0));
    total_clean_xy = 0;

    std::string line;
    while (std::getline(in, line)) {
        if (is_empty_or_header(line)) continue;
        auto fields = split_csv_line(line);
        if (fields.size() < 3) continue;

        int x = -1;
        int y = -1;
        long long c = 0;
        if (!parse_int(fields[0], x)) continue;
        if (!parse_int(fields[1], y)) continue;
        if (!parse_long_long(fields[2], c)) continue;
        if (x < 0 || x >= 16 || y < 0 || y >= 16) continue;
        counts[x][y] = c;
        total_clean_xy += c;
    }
    return true;
}

static void draw_text(TCanvas* c, const std::string& text, double x, double y, double size = 0.035) {
    TLatex latex;
    latex.SetNDC(kTRUE);
    latex.SetTextSize(size);
    latex.DrawLatex(x, y, text.c_str());
    c->Modified();
    c->Update();
}

static void save_canvas(TCanvas* c, const fs::path& path) {
    c->SaveAs(path.string().c_str());
}

static TH1D* make_bar_hist(const std::string& name, const std::string& title, int nbins, double xmin, double xmax, const std::string& xlab, const std::string& ylab) {
    TH1D* h = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(xlab.c_str());
    h->GetYaxis()->SetTitle(ylab.c_str());
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    return h;
}

int main(int argc, char** argv) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetFrameLineWidth(1);

    fs::path event_csv = "hod2018_event_summary.csv";
    fs::path channel_csv = "hod2018_channel_summary.csv";
    fs::path xy_csv = "hod2018_xy_map.csv";

    if (argc >= 2) event_csv = argv[1];
    if (argc >= 3) channel_csv = argv[2];
    if (argc >= 4) xy_csv = argv[3];

    const fs::path out_dir = fs::path("outputs") / "hod2018_diagnostics";
    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "ERROR: no pude crear " << out_dir.string() << ": " << ec.message() << "\n";
        return 1;
    }

    if (!fs::exists(event_csv)) {
        std::cerr << "ERROR: no existe " << event_csv.string() << "\n";
        return 1;
    }
    if (!fs::exists(channel_csv)) {
        std::cerr << "ERROR: no existe " << channel_csv.string() << "\n";
        return 1;
    }
    if (!fs::exists(xy_csv)) {
        std::cerr << "ERROR: no existe " << xy_csv.string() << "\n";
        return 1;
    }

    std::vector<ChannelRow> channels;
    std::vector<EventRow> events;
    std::map<std::string, long long> class_counts;
    std::vector<std::vector<long long>> xy_counts;
    long long total_clean_xy = 0;

    std::string err;
    if (!load_channel_summary(channel_csv, channels, err)) {
        std::cerr << "ERROR: " << err << "\n";
        return 1;
    }
    if (!load_event_summary(event_csv, events, class_counts, err)) {
        std::cerr << "ERROR: " << err << "\n";
        return 1;
    }
    if (!load_xy_map(xy_csv, xy_counts, total_clean_xy, err)) {
        std::cerr << "ERROR: " << err << "\n";
        return 1;
    }

    std::map<int, ChannelRow> by_channel;
    for (const auto& row : channels) by_channel[row.channel] = row;

    long long total_events = static_cast<long long>(events.size());
    auto class_get = [&](const std::string& c) -> long long {
        auto it = class_counts.find(c);
        return (it == class_counts.end()) ? 0LL : it->second;
    };

    const long long n_empty = class_get("empty");
    const long long n_X_only = class_get("X_only");
    const long long n_Y_only = class_get("Y_only");
    const long long n_XY_candidate = class_get("XY_candidate") + class_get("clean_XY") + class_get("multi_XY");
    const long long n_clean_XY = class_get("clean_XY");
    const long long n_multi_XY = class_get("multi_XY");

    std::vector<ChannelRow> x_rows;
    std::vector<ChannelRow> y_rows;
    for (const auto& row : channels) {
        if (row.plane == "X") x_rows.push_back(row);
        if (row.plane == "Y") y_rows.push_back(row);
    }
    auto by_hits_desc = [](const ChannelRow& a, const ChannelRow& b) {
        if (a.hits != b.hits) return a.hits > b.hits;
        return a.channel < b.channel;
    };
    std::sort(x_rows.begin(), x_rows.end(), by_hits_desc);
    std::sort(y_rows.begin(), y_rows.end(), by_hits_desc);

    long long total_hits_x = 0;
    long long total_hits_y = 0;
    int top_x_channel = -1;
    int top_y_channel = -1;
    long long top_x_hits = -1;
    long long top_y_hits = -1;
    std::vector<int> dead_channels;

    for (const auto& row : channels) {
        if (row.plane == "X") total_hits_x += row.hits;
        if (row.plane == "Y") total_hits_y += row.hits;
        if (row.hits == 0) dead_channels.push_back(row.channel);
    }
    for (const auto& row : x_rows) {
        if (row.hits > top_x_hits) {
            top_x_hits = row.hits;
            top_x_channel = row.channel;
        }
    }
    for (const auto& row : y_rows) {
        if (row.hits > top_y_hits) {
            top_y_hits = row.hits;
            top_y_channel = row.channel;
        }
    }

    auto fraction = [](long long a, long long b) -> double {
        return (b > 0) ? static_cast<double>(a) / static_cast<double>(b) : 0.0;
    };

    // Canvas helpers
    auto make_canvas = [](const std::string& name, int w = 1100, int h = 800) {
        TCanvas* c = new TCanvas(name.c_str(), name.c_str(), w, h);
        c->SetMargin(0.10, 0.06, 0.11, 0.08);
        return c;
    };

    // 1) Channel occupancy plots
    TH1D* hx = make_bar_hist("hx", "HOD2018 X occupancy;X channel;Hits", 16, -0.5, 15.5, "X channel", "Hits");
    TH1D* hy = make_bar_hist("hy", "HOD2018 Y occupancy;Y channel (Ch16=0);Hits", 16, -0.5, 15.5, "Y channel (Ch16=0)", "Hits");
    TH1D* hxy_x = make_bar_hist("hxy_x", "X occupancy;X channel;Hits", 16, -0.5, 15.5, "X channel", "Hits");
    TH1D* hxy_y = make_bar_hist("hxy_y", "Y occupancy;Y channel (Ch16=0);Hits", 16, -0.5, 15.5, "Y channel (Ch16=0)", "Hits");
    TH1D* hfrac = make_bar_hist("hfrac", "Frac. hits by channel;Channel;frac_hits", 32, -0.5, 31.5, "Channel", "frac_hits");
    TH1D* hamp = make_bar_hist("hamp", "Max negative amplitude by channel;Channel;max_amp_neg", 32, -0.5, 31.5, "Channel", "max_amp_neg");

    for (const auto& row : channels) {
        if (row.channel >= 0 && row.channel < 32) {
            hfrac->SetBinContent(row.channel + 1, row.frac_hits);
            hamp->SetBinContent(row.channel + 1, row.max_amp_neg);
        }
        if (row.plane == "X" && row.channel >= 0 && row.channel < 16) {
            hx->SetBinContent(row.channel + 1, row.hits);
            hxy_x->SetBinContent(row.channel + 1, row.hits);
        }
        if (row.plane == "Y" && row.channel >= 16 && row.channel < 32) {
            int yidx = row.channel - 16;
            hy->SetBinContent(yidx + 1, row.hits);
            hxy_y->SetBinContent(yidx + 1, row.hits);
        }
    }

    {
        TCanvas* c = make_canvas("c_x_occ");
        hx->SetFillColor(kAzure + 1);
        hx->SetLineColor(kAzure + 4);
        hx->Draw("hist");
        save_canvas(c, out_dir / "hod2018_x_occupancy.png");
    }
    {
        TCanvas* c = make_canvas("c_y_occ");
        hy->SetFillColor(kOrange + 1);
        hy->SetLineColor(kOrange + 7);
        hy->Draw("hist");
        save_canvas(c, out_dir / "hod2018_y_occupancy.png");
    }
    {
        TCanvas* c = make_canvas("c_xy_occ");
        hx->SetFillColor(kAzure + 1);
        hx->SetLineColor(kAzure + 4);
        hy->SetFillColorAlpha(kOrange + 1, 0.45);
        hy->SetLineColor(kOrange + 7);
        hx->Draw("hist");
        hy->Draw("hist same");
        TLegend* leg = new TLegend(0.72, 0.75, 0.92, 0.90);
        leg->AddEntry(hx, "X (Ch0-15)", "f");
        leg->AddEntry(hy, "Y (Ch16-31 -> 0-15)", "f");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();
        save_canvas(c, out_dir / "hod2018_xy_occupancy_comparison.png");
    }
    {
        TCanvas* c = make_canvas("c_frac");
        hfrac->SetFillColor(kGreen + 1);
        hfrac->SetLineColor(kGreen + 3);
        hfrac->Draw("hist");
        save_canvas(c, out_dir / "hod2018_frac_hits_by_channel.png");
    }
    {
        TCanvas* c = make_canvas("c_amp");
        hamp->SetFillColor(kMagenta + 1);
        hamp->SetLineColor(kMagenta + 3);
        hamp->Draw("hist");
        save_canvas(c, out_dir / "hod2018_max_amp_neg_by_channel.png");
    }

    // 2) XY map plots
    TH2D* hxy = new TH2D("hxy", "HOD2018 XY map;X channel;Y channel (0-15)", 16, -0.5, 15.5, 16, -0.5, 15.5);
    TH2D* hxy_norm = new TH2D("hxy_norm", "HOD2018 XY normalized map;X channel;Y channel (0-15)", 16, -0.5, 15.5, 16, -0.5, 15.5);
    for (int x = 0; x < 16; ++x) {
        for (int y = 0; y < 16; ++y) {
            hxy->SetBinContent(x + 1, y + 1, xy_counts[x][y]);
            hxy_norm->SetBinContent(x + 1, y + 1, fraction(xy_counts[x][y], total_clean_xy));
        }
    }

    {
        TCanvas* c = make_canvas("c_xy_counts");
        hxy->SetContour(50);
        hxy->Draw("COLZ TEXT");
        save_canvas(c, out_dir / "hod2018_xy_map_counts.png");
    }
    {
        TCanvas* c = make_canvas("c_xy_logz");
        c->SetLogz(true);
        hxy->SetContour(50);
        hxy->Draw("COLZ TEXT");
        save_canvas(c, out_dir / "hod2018_xy_map_logz.png");
    }
    {
        TCanvas* c = make_canvas("c_xy_norm");
        hxy_norm->SetContour(50);
        hxy_norm->Draw("COLZ TEXT");
        save_canvas(c, out_dir / "hod2018_xy_map_normalized.png");
    }

    // 3) Event summary plots
    TH1D* h_nhits = new TH1D("h_nhits", "Event multiplicity;nHits;Events", 33, -0.5, 32.5);
    TH2D* h_nxny = new TH2D("h_nxny", "nX vs nY;nX;nY", 17, -0.5, 16.5, 17, -0.5, 16.5);
    TH1D* h_maxx = new TH1D("h_maxx", "max_amp_X for events with max_amp_X>0;max_amp_X;Events", 100, 0.0, 250.0);
    TH1D* h_maxy = new TH1D("h_maxy", "max_amp_Y for events with max_amp_Y>0;max_amp_Y;Events", 100, 0.0, 250.0);
    TH1D* h_classes = new TH1D("h_classes", "Event classes;Class;Events", 6, -0.5, 5.5);
    h_classes->GetXaxis()->SetBinLabel(1, "empty");
    h_classes->GetXaxis()->SetBinLabel(2, "X_only");
    h_classes->GetXaxis()->SetBinLabel(3, "Y_only");
    h_classes->GetXaxis()->SetBinLabel(4, "XY_candidate");
    h_classes->GetXaxis()->SetBinLabel(5, "clean_XY");
    h_classes->GetXaxis()->SetBinLabel(6, "multi_XY");
    h_classes->SetStats(kFALSE);

    for (const auto& ev : events) {
        h_nhits->Fill(ev.nHits);
        h_nxny->Fill(ev.nX, ev.nY);
        if (ev.max_amp_X > 0.0) h_maxx->Fill(ev.max_amp_X);
        if (ev.max_amp_Y > 0.0) h_maxy->Fill(ev.max_amp_Y);

        if (ev.nX == 0 && ev.nY == 0) {
            h_classes->Fill(0.0);
        } else if (ev.nX > 0 && ev.nY == 0) {
            h_classes->Fill(1.0);
        } else if (ev.nX == 0 && ev.nY > 0) {
            h_classes->Fill(2.0);
        } else {
            h_classes->Fill(3.0);
            if (ev.nX == 1 && ev.nY == 1) {
                h_classes->Fill(4.0);
            }
            if (ev.nX > 0 && ev.nY > 0 && (ev.nX + ev.nY) > 2) {
                h_classes->Fill(5.0);
            }
        }
    }

    {
        TCanvas* c = make_canvas("c_classes");
        h_classes->SetFillColor(kTeal + 1);
        h_classes->SetLineColor(kTeal + 4);
        h_classes->Draw("hist");
        save_canvas(c, out_dir / "hod2018_event_classes.png");
    }
    {
        TCanvas* c = make_canvas("c_mult");
        h_nhits->SetFillColor(kBlue - 7);
        h_nhits->SetLineColor(kBlue + 2);
        h_nhits->Draw("hist");
        save_canvas(c, out_dir / "hod2018_multiplicity.png");
    }
    {
        TCanvas* c = make_canvas("c_nxny");
        h_nxny->SetContour(50);
        h_nxny->Draw("COLZ TEXT");
        save_canvas(c, out_dir / "hod2018_nx_ny_map.png");
    }
    {
        TCanvas* c = make_canvas("c_maxx");
        h_maxx->SetFillColor(kRed - 7);
        h_maxx->SetLineColor(kRed + 1);
        h_maxx->Draw("hist");
        save_canvas(c, out_dir / "hod2018_max_amp_X.png");
    }
    {
        TCanvas* c = make_canvas("c_maxy");
        h_maxy->SetFillColor(kViolet - 7);
        h_maxy->SetLineColor(kViolet + 1);
        h_maxy->Draw("hist");
        save_canvas(c, out_dir / "hod2018_max_amp_Y.png");
    }

    // Summary text
    std::ostringstream summary;
    summary << "Total events: " << total_events << "\n";
    summary << "empty: " << n_empty << "\n";
    summary << "X_only: " << n_X_only << "\n";
    summary << "Y_only: " << n_Y_only << "\n";
    summary << "XY_candidate: " << n_XY_candidate << "\n";
    summary << "clean_XY: " << n_clean_XY << "\n";
    summary << "multi_XY: " << n_multi_XY << "\n";
    summary << "clean_XY / total: " << std::fixed << std::setprecision(6) << fraction(n_clean_XY, total_events) << "\n";
    summary << "XY_candidate / total: " << std::fixed << std::setprecision(6) << fraction(n_XY_candidate, total_events) << "\n";
    summary << "clean_XY used in map: " << total_clean_xy << "\n";
    summary << "Most active X channel: " << top_x_channel << " (hits=" << top_x_hits << ")\n";
    summary << "Most active Y channel: " << top_y_channel << " (hits=" << top_y_hits << ")\n";

    summary << "Channels with zero hits: ";
    if (dead_channels.empty()) {
        summary << "none\n";
    } else {
        for (size_t i = 0; i < dead_channels.size(); ++i) {
            if (i) summary << ",";
            summary << dead_channels[i];
        }
        summary << "\n";
    }

    if (total_hits_x > 0 && by_channel.count(15) && by_channel[15].hits > 0 && fraction(by_channel[15].hits, total_hits_x) > 0.5) {
        summary << "Comment: X15 dominates more than 50% of X hits.\n";
    } else {
        summary << "Comment: X15 does not dominate more than 50% of X hits.\n";
    }

    bool any_y_dominates = false;
    int dom_y_channel = -1;
    double dom_y_frac = 0.0;
    for (const auto& row : y_rows) {
        double frac = fraction(row.hits, total_hits_y);
        if (frac > dom_y_frac) {
            dom_y_frac = frac;
            dom_y_channel = row.channel;
        }
        if (frac > 0.5) any_y_dominates = true;
    }
    if (any_y_dominates && dom_y_channel >= 0) {
        summary << "Comment: Y channel " << dom_y_channel << " dominates more than 50% of Y hits.\n";
    } else {
        summary << "Comment: no Y channel dominates more than 50% of Y hits.\n";
    }

    if (!dead_channels.empty()) {
        summary << "Comment: dead channels detected.\n";
    } else {
        summary << "Comment: no dead channels detected.\n";
    }

    std::ofstream summary_out(out_dir / "summary.txt");
    if (!summary_out) {
        std::cerr << "ERROR: no pude escribir " << (out_dir / "summary.txt").string() << "\n";
        return 1;
    }
    summary_out << summary.str();
    summary_out.close();

    std::cout << "Gráficos guardados en outputs/hod2018_diagnostics/\n";
    return 0;
}
