#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "hodo_mapping.h"

namespace fs = std::filesystem;

constexpr double kRadToDeg = 180.0 / M_PI;

static const auto& kHod2019XChannelsByIndex = kH19X;
static const auto& kHod2019YChannelsByIndex = kH19Y;


static double coord_mm_from_index(int idx, double pitch) {
    return (static_cast<double>(idx) + 0.5 - 8.0) * pitch;
}

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
    while (std::getline(ss, field, ',')) fields.push_back(trim(field));
    if (!line.empty() && line.back() == ',') fields.emplace_back("");
    return fields;
}

static std::vector<int> parse_channel_list(const std::string& text) {
    std::vector<int> channels;
    const std::string clean = trim(text);
    if (clean.empty()) return channels;

    std::stringstream ss(clean);
    std::string item;
    while (std::getline(ss, item, ';')) {
        item = trim(item);
        if (item.empty()) continue;
        try {
            channels.push_back(std::stoi(item));
        } catch (...) {
            std::cerr << "ADVERTENCIA: no pude parsear canal '" << item << "'\n";
        }
    }
    return channels;
}

static double parse_double_or_zero(const std::string& text) {
    try {
        return std::stod(trim(text));
    } catch (...) {
        return 0.0;
    }
}

static std::string join_ints(const std::vector<int>& values) {
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) os << ';';
        os << values[i];
    }
    return os.str();
}

static std::string format_double_or_blank(bool valid, double value) {
    if (!valid) return "";
    std::ostringstream os;
    os << std::setprecision(10) << value;
    return os.str();
}

static double mean_coord_from_indices(const std::vector<int>& indices, double pitch) {
    if (indices.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for (int idx : indices) sum += coord_mm_from_index(idx, pitch);
    return sum / static_cast<double>(indices.size());
}

struct CandidateRow {
    std::string entry;
    std::string event;
    std::string timestamp;
    std::vector<int> channels_hod2018;
    std::vector<int> channels_hod2019;
    double max_amp_hod2018 = 0.0;
    double max_amp_hod2019 = 0.0;
};

struct TrackReco {
    CandidateRow row;

    std::vector<int> ch_x2018;
    std::vector<int> ch_y2018;
    std::vector<int> ch_x2019;
    std::vector<int> ch_y2019;

    std::vector<int> idx_x2018;
    std::vector<int> idx_y2018;
    std::vector<int> idx_x2019;
    std::vector<int> idx_y2019;

    bool has_x2018 = false;
    bool has_y2018 = false;
    bool has_x2019 = false;
    bool has_y2019 = false;
    bool has_xz = false;
    bool has_yz = false;
    bool has_full3d = false;

    double x2018_mm = std::numeric_limits<double>::quiet_NaN();
    double y2018_mm = std::numeric_limits<double>::quiet_NaN();
    double x2019_mm = std::numeric_limits<double>::quiet_NaN();
    double y2019_mm = std::numeric_limits<double>::quiet_NaN();

    double mx = std::numeric_limits<double>::quiet_NaN();
    double my = std::numeric_limits<double>::quiet_NaN();
    double theta_x_deg = std::numeric_limits<double>::quiet_NaN();
    double theta_y_deg = std::numeric_limits<double>::quiet_NaN();
    double theta_3d_deg = std::numeric_limits<double>::quiet_NaN();
};

static bool load_candidates(const fs::path& input, std::vector<CandidateRow>& rows) {
    std::ifstream in(input);
    if (!in) {
        std::cerr << "ERROR: no existe o no pude abrir " << input << "\n";
        return false;
    }

    std::string line;
    bool first_line = true;
    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        if (first_line) {
            first_line = false;
            continue;
        }

        const auto fields = split_csv_line(line);
        if (fields.size() < 9) {
            std::cerr << "ADVERTENCIA: linea con pocas columnas, saltando: " << line << "\n";
            continue;
        }

        CandidateRow row;
        row.entry = fields[0];
        row.event = fields[1];
        row.timestamp = fields[2];
        row.channels_hod2018 = parse_channel_list(fields[5]);
        row.channels_hod2019 = parse_channel_list(fields[6]);
        row.max_amp_hod2018 = parse_double_or_zero(fields[7]);
        row.max_amp_hod2019 = parse_double_or_zero(fields[8]);
        rows.push_back(std::move(row));
    }
    return true;
}

static TrackReco reconstruct_candidate(const CandidateRow& row, double pitch, double z_sep_mm) {
    TrackReco r;
    r.row = row;

    for (int ch : row.channels_hod2018) {
        const int ix = hod2018_x_index(ch);
        const int iy = hod2018_y_index(ch);
        if (ix >= 0) {
            r.ch_x2018.push_back(ch);
            r.idx_x2018.push_back(ix);
        }
        if (iy >= 0) {
            r.ch_y2018.push_back(ch);
            r.idx_y2018.push_back(iy);
        }
    }

    for (int ch : row.channels_hod2019) {
        const int ix = hod2019_x_index(ch);
        const int iy = hod2019_y_index(ch);
        if (ix >= 0) {
            r.ch_x2019.push_back(ch);
            r.idx_x2019.push_back(ix);
        }
        if (iy >= 0) {
            r.ch_y2019.push_back(ch);
            r.idx_y2019.push_back(iy);
        }
    }

    r.has_x2018 = !r.idx_x2018.empty();
    r.has_y2018 = !r.idx_y2018.empty();
    r.has_x2019 = !r.idx_x2019.empty();
    r.has_y2019 = !r.idx_y2019.empty();
    r.has_xz = r.has_x2018 && r.has_x2019;
    r.has_yz = r.has_y2018 && r.has_y2019;
    r.has_full3d = r.has_xz && r.has_yz;

    if (r.has_x2018) r.x2018_mm = mean_coord_from_indices(r.idx_x2018, pitch);
    if (r.has_y2018) r.y2018_mm = mean_coord_from_indices(r.idx_y2018, pitch);
    if (r.has_x2019) r.x2019_mm = mean_coord_from_indices(r.idx_x2019, pitch);
    if (r.has_y2019) r.y2019_mm = mean_coord_from_indices(r.idx_y2019, pitch);

    if (r.has_xz) {
        r.mx = (r.x2019_mm - r.x2018_mm) / z_sep_mm;
        r.theta_x_deg = std::atan(r.mx) * kRadToDeg;
    }
    if (r.has_yz) {
        r.my = (r.y2019_mm - r.y2018_mm) / z_sep_mm;
        r.theta_y_deg = std::atan(r.my) * kRadToDeg;
    }
    if (r.has_full3d) {
        r.theta_3d_deg = std::atan(std::sqrt(r.mx * r.mx + r.my * r.my)) * kRadToDeg;
    }

    return r;
}

static void configure_hist_1d(TH1D* h, const std::string& x_title, const std::string& y_title) {
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(x_title.c_str());
    h->GetYaxis()->SetTitle(y_title.c_str());
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->SetLineWidth(2);
}

static void configure_hist_2d(TH2D* h, const std::string& x_title, const std::string& y_title) {
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(x_title.c_str());
    h->GetYaxis()->SetTitle(y_title.c_str());
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
}

static void save_hist1_png(TH1D* h, const fs::path& path) {
    TCanvas c(("c_" + std::string(h->GetName())).c_str(), "", 1100, 760);
    c.SetLeftMargin(0.11);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.18);
    h->Draw("HIST");
    c.SaveAs(path.string().c_str());
}

static void save_hist2_png(TH2D* h, const fs::path& path) {
    TCanvas c(("c_" + std::string(h->GetName())).c_str(), "", 900, 820);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.16);
    c.SetBottomMargin(0.12);
    h->Draw("COLZ");
    c.SaveAs(path.string().c_str());
}

static void draw_projection_lines(const std::vector<TrackReco>& tracks,
                                  bool use_x,
                                  double active_size_mm,
                                  double z_sep_mm,
                                  const fs::path& path) {
    const double half = active_size_mm * 0.5;
    const std::string axis = use_x ? "X" : "Y";
    TH2D frame(("frame_" + axis).c_str(),
               (axis + "Z projected candidate tracks;" + axis + " [mm];Z [mm]").c_str(),
               10, -half, half, 10, -2.0, z_sep_mm + 2.0);
    frame.SetStats(kFALSE);

    TCanvas c(("c_" + axis + "z_projected_tracks").c_str(), "", 900, 900);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.12);
    frame.Draw();

    std::vector<std::unique_ptr<TLine>> lines;
    lines.reserve(tracks.size() + 2);

    auto plane0 = std::make_unique<TLine>(-half, 0.0, half, 0.0);
    plane0->SetLineColor(kGray + 2);
    plane0->SetLineWidth(2);
    plane0->Draw();
    lines.push_back(std::move(plane0));

    auto plane1 = std::make_unique<TLine>(-half, z_sep_mm, half, z_sep_mm);
    plane1->SetLineColor(kGray + 2);
    plane1->SetLineWidth(2);
    plane1->SetLineStyle(2);
    plane1->Draw();
    lines.push_back(std::move(plane1));

    for (const TrackReco& t : tracks) {
        if (use_x && !t.has_xz) continue;
        if (!use_x && !t.has_yz) continue;

        const double a = use_x ? t.x2018_mm : t.y2018_mm;
        const double b = use_x ? t.x2019_mm : t.y2019_mm;
        auto line = std::make_unique<TLine>(a, 0.0, b, z_sep_mm);
        line->SetLineColor(use_x ? kAzure + 2 : kOrange + 7);
        line->SetLineWidth(1);
        line->Draw();
        lines.push_back(std::move(line));
    }

    c.SaveAs(path.string().c_str());
}

static void write_reconstructed_csv(const fs::path& outdir, const std::vector<TrackReco>& tracks) {
    std::ofstream out(outdir / "reconstructed_tracks_interleaved.csv");
    out << "entry,event,timestamp,"
        << "ch_x2018,ch_y2018,ch_x2019,ch_y2019,"
        << "idx_x2018,idx_y2018,idx_x2019,idx_y2019,"
        << "has_xz,has_yz,has_full3d,"
        << "x2018_mm,y2018_mm,x2019_mm,y2019_mm,"
        << "mx,my,theta_x_deg,theta_y_deg,theta_3d_deg,"
        << "max_amp_hod2018,max_amp_hod2019\n";

    for (const TrackReco& t : tracks) {
        out << t.row.entry << ','
            << t.row.event << ','
            << t.row.timestamp << ','
            << join_ints(t.ch_x2018) << ','
            << join_ints(t.ch_y2018) << ','
            << join_ints(t.ch_x2019) << ','
            << join_ints(t.ch_y2019) << ','
            << join_ints(t.idx_x2018) << ','
            << join_ints(t.idx_y2018) << ','
            << join_ints(t.idx_x2019) << ','
            << join_ints(t.idx_y2019) << ','
            << (t.has_xz ? 1 : 0) << ','
            << (t.has_yz ? 1 : 0) << ','
            << (t.has_full3d ? 1 : 0) << ','
            << format_double_or_blank(t.has_x2018, t.x2018_mm) << ','
            << format_double_or_blank(t.has_y2018, t.y2018_mm) << ','
            << format_double_or_blank(t.has_x2019, t.x2019_mm) << ','
            << format_double_or_blank(t.has_y2019, t.y2019_mm) << ','
            << format_double_or_blank(t.has_xz, t.mx) << ','
            << format_double_or_blank(t.has_yz, t.my) << ','
            << format_double_or_blank(t.has_xz, t.theta_x_deg) << ','
            << format_double_or_blank(t.has_yz, t.theta_y_deg) << ','
            << format_double_or_blank(t.has_full3d, t.theta_3d_deg) << ','
            << t.row.max_amp_hod2018 << ','
            << t.row.max_amp_hod2019 << '\n';
    }
}

static void write_hod2019_physical_mapping_csv(const fs::path& outdir) {
    std::ofstream out(outdir / "hod2019_physical_mapping.csv");
    out << "plane,index,channel\n";
    for (size_t i = 0; i < kHod2019XChannelsByIndex.size(); ++i) {
        out << "X," << i << ',' << kHod2019XChannelsByIndex[i] << '\n';
    }
    for (size_t i = 0; i < kHod2019YChannelsByIndex.size(); ++i) {
        out << "Y," << i << ',' << kHod2019YChannelsByIndex[i] << '\n';
    }
}

static void write_summary(const fs::path& outdir,
                          const std::vector<TrackReco>& tracks,
                          double active_size_mm,
                          double z_sep_mm,
                          double pitch,
                          long long n_xz,
                          long long n_yz,
                          long long n_full3d,
                          long long n_xz_only,
                          long long n_yz_only,
                          long long n_no_projection) {
    std::ofstream summary(outdir / "summary.txt");
    summary << "Track reconstruction with interleaved hod2019 mapping\n";
    summary << "====================================================\n\n";
    summary << "Total candidates: " << tracks.size() << "\n";
    summary << "With XZ projection: " << n_xz << "\n";
    summary << "With YZ projection: " << n_yz << "\n";
    summary << "Full 3D: " << n_full3d << "\n";
    summary << "Only XZ: " << n_xz_only << "\n";
    summary << "Only YZ: " << n_yz_only << "\n";
    summary << "No reconstructible projection: " << n_no_projection << "\n\n";

    summary << "Geometry:\n";
    summary << "  active_size_mm = " << active_size_mm << "\n";
    summary << "  z_sep_mm = " << z_sep_mm << "\n";
    summary << "  pitch = " << pitch << "\n\n";

    summary << "Mapping used:\n";
    summary << "  hod2018 X: X0=Ch0, X1=Ch1, X2=Ch2, X3=Ch3, X4=Ch4, X5=Ch5, X6=Ch6, X7=Ch7, X8=Ch8, X9=Ch9, X10=Ch10, X11=Ch11, X12=Ch12, X13=Ch13, X14=Ch14, X15=Ch15\n";
    summary << "  hod2018 Y: Y0=Ch16, Y1=Ch17, Y2=Ch18, Y3=Ch19, Y4=Ch20, Y5=Ch21, Y6=Ch22, Y7=Ch23, Y8=Ch24, Y9=Ch25, Y10=Ch26, Y11=Ch27, Y12=Ch28, Y13=Ch29, Y14=Ch30, Y15=Ch31\n";
    summary << "  hod2019 X literal order: X0=Ch53, X1=Ch61, X2=Ch54, X3=Ch62, X4=Ch55, X5=Ch63, X6=Ch52, X7=Ch60, X8=Ch51, X9=Ch59, X10=Ch50, X11=Ch58, X12=Ch49, X13=Ch57, X14=Ch48, X15=Ch56\n";
    summary << "  hod2019 Y literal order: Y0=Ch37, Y1=Ch45, Y2=Ch38, Y3=Ch46, Y4=Ch39, Y5=Ch47, Y6=Ch36, Y7=Ch44, Y8=Ch35, Y9=Ch43, Y10=Ch34, Y11=Ch42, Y12=Ch33, Y13=Ch41, Y14=Ch32, Y15=Ch40\n";
    summary << "  Physical coordinate: coord_mm = (index + 0.5 - 8.0) * pitch\n\n";

    if (n_full3d == 0) {
        summary << "Comment: full3D = 0, so no event has both XZ and YZ projections simultaneously.\n";
    }
    summary << "Comment: these are candidate tracks because hod2019 has few hits and should be treated as a sparse/diagnostic plane.\n";
}

int main(int argc, char** argv) {
    const fs::path input = (argc >= 2) ? fs::path(argv[1]) : fs::path("candidate_events_windowed.csv");
    const fs::path outdir = (argc >= 3) ? fs::path(argv[2]) : fs::path("outputs/track_reconstruction_interleaved");
    const double active_size_mm = (argc >= 4) ? std::stod(argv[3]) : 33.0;
    const double z_sep_mm = (argc >= 5) ? std::stod(argv[4]) : 98.4;
    const double pitch = active_size_mm / 16.0;

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    std::vector<CandidateRow> candidates;
    if (!load_candidates(input, candidates)) return 2;

    fs::create_directories(outdir);

    std::vector<TrackReco> tracks;
    tracks.reserve(candidates.size());
    for (const CandidateRow& row : candidates) {
        tracks.push_back(reconstruct_candidate(row, pitch, z_sep_mm));
    }

    long long n_xz = 0;
    long long n_yz = 0;
    long long n_full3d = 0;
    long long n_xz_only = 0;
    long long n_yz_only = 0;
    long long n_no_projection = 0;

    const double half = active_size_mm * 0.5;
    auto h_theta_x = std::make_unique<TH1D>("h_theta_x", "XZ angle;#theta_{x} [deg];Candidate tracks", 100, -30.0, 30.0);
    auto h_theta_y = std::make_unique<TH1D>("h_theta_y", "YZ angle;#theta_{y} [deg];Candidate tracks", 100, -30.0, 30.0);
    auto h_theta_3d = std::make_unique<TH1D>("h_theta_3d", "3D angle;#theta_{3D} [deg];Candidate tracks", 100, 0.0, 35.0);
    auto h_hod2018_xy = std::make_unique<TH2D>("h_hod2018_xy_points", "hod2018 XY points;X2018 [mm];Y2018 [mm];Events", 16, -half, half, 16, -half, half);
    auto h_hod2019_xy = std::make_unique<TH2D>("h_hod2019_xy_points", "hod2019 XY points;X2019 [mm];Y2019 [mm];Events", 16, -half, half, 16, -half, half);
    auto h_x_occ = std::make_unique<TH1D>("h_hod2019_x_interleaved_occupancy", "hod2019 X interleaved occupancy;X index;Hits", 16, -0.5, 15.5);
    auto h_y_occ = std::make_unique<TH1D>("h_hod2019_y_interleaved_occupancy", "hod2019 Y interleaved occupancy;Y index;Hits", 16, -0.5, 15.5);
    auto h_categories = std::make_unique<TH1D>("h_projection_category_counts", "Projection categories;Category;Events", 4, 0.5, 4.5);

    configure_hist_1d(h_theta_x.get(), "#theta_{x} [deg]", "Candidate tracks");
    configure_hist_1d(h_theta_y.get(), "#theta_{y} [deg]", "Candidate tracks");
    configure_hist_1d(h_theta_3d.get(), "#theta_{3D} [deg]", "Candidate tracks");
    configure_hist_1d(h_x_occ.get(), "X index", "Hits");
    configure_hist_1d(h_y_occ.get(), "Y index", "Hits");
    configure_hist_1d(h_categories.get(), "Category", "Events");
    configure_hist_2d(h_hod2018_xy.get(), "X2018 [mm]", "Y2018 [mm]");
    configure_hist_2d(h_hod2019_xy.get(), "X2019 [mm]", "Y2019 [mm]");

    for (int i = 0; i < 16; ++i) {
        h_x_occ->GetXaxis()->SetBinLabel(i + 1, (std::to_string(i) + "/Ch" + std::to_string(kHod2019XChannelsByIndex[i])).c_str());
        h_y_occ->GetXaxis()->SetBinLabel(i + 1, (std::to_string(i) + "/Ch" + std::to_string(kHod2019YChannelsByIndex[i])).c_str());
    }
    h_x_occ->LabelsOption("v", "X");
    h_y_occ->LabelsOption("v", "X");

    h_categories->GetXaxis()->SetBinLabel(1, "no_projection");
    h_categories->GetXaxis()->SetBinLabel(2, "xz_only");
    h_categories->GetXaxis()->SetBinLabel(3, "yz_only");
    h_categories->GetXaxis()->SetBinLabel(4, "full3d");

    for (const TrackReco& t : tracks) {
        if (t.has_xz) {
            ++n_xz;
            h_theta_x->Fill(t.theta_x_deg);
        }
        if (t.has_yz) {
            ++n_yz;
            h_theta_y->Fill(t.theta_y_deg);
        }
        if (t.has_full3d) {
            ++n_full3d;
            h_theta_3d->Fill(t.theta_3d_deg);
            h_categories->Fill(4);
        } else if (t.has_xz) {
            ++n_xz_only;
            h_categories->Fill(2);
        } else if (t.has_yz) {
            ++n_yz_only;
            h_categories->Fill(3);
        } else {
            ++n_no_projection;
            h_categories->Fill(1);
        }

        if (t.has_x2018 && t.has_y2018) h_hod2018_xy->Fill(t.x2018_mm, t.y2018_mm);
        if (t.has_x2019 && t.has_y2019) h_hod2019_xy->Fill(t.x2019_mm, t.y2019_mm);
        for (int idx : t.idx_x2019) h_x_occ->Fill(idx);
        for (int idx : t.idx_y2019) h_y_occ->Fill(idx);
    }

    write_reconstructed_csv(outdir, tracks);
    write_hod2019_physical_mapping_csv(outdir);
    write_summary(outdir, tracks, active_size_mm, z_sep_mm, pitch,
                  n_xz, n_yz, n_full3d, n_xz_only, n_yz_only, n_no_projection);

    draw_projection_lines(tracks, true, active_size_mm, z_sep_mm, outdir / "xz_projected_tracks.png");
    draw_projection_lines(tracks, false, active_size_mm, z_sep_mm, outdir / "yz_projected_tracks.png");
    save_hist1_png(h_theta_x.get(), outdir / "theta_x.png");
    save_hist1_png(h_theta_y.get(), outdir / "theta_y.png");
    save_hist1_png(h_theta_3d.get(), outdir / "theta_3d.png");
    save_hist2_png(h_hod2018_xy.get(), outdir / "hod2018_xy_points.png");
    save_hist2_png(h_hod2019_xy.get(), outdir / "hod2019_xy_points.png");
    save_hist1_png(h_x_occ.get(), outdir / "hod2019_x_interleaved_occupancy.png");
    save_hist1_png(h_y_occ.get(), outdir / "hod2019_y_interleaved_occupancy.png");
    save_hist1_png(h_categories.get(), outdir / "projection_category_counts.png");

    TFile fout((outdir / "track_reconstruction_interleaved.root").string().c_str(), "RECREATE");
    if (!fout.IsZombie()) {
        h_theta_x->Write();
        h_theta_y->Write();
        h_theta_3d->Write();
        h_hod2018_xy->Write();
        h_hod2019_xy->Write();
        h_x_occ->Write();
        h_y_occ->Write();
        h_categories->Write();
        fout.Close();
    }

    std::cout << "Total candidates: " << tracks.size() << "\n";
    std::cout << "With XZ: " << n_xz << "\n";
    std::cout << "With YZ: " << n_yz << "\n";
    std::cout << "Full 3D: " << n_full3d << "\n";
    std::cout << "Output directory: " << outdir << "\n";
    return 0;
}
