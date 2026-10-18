#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <iterator>
#include <vector>
#include "hodo_paths.h"

namespace fs = std::filesystem;
constexpr int NCH = 64;
constexpr int NRUN = 9;
constexpr double PRE_SECONDS = 3.2e-6;
constexpr std::array<int, 16> H19X = {53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
constexpr std::array<int, 16> H19Y = {37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};
const std::array<std::string, NRUN> RUN_NAMES = {
    "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9"
};
const std::array<std::string, NRUN> RUN_FILES = {
    "lectura_vx1740_ext_nim_1M.root",
    "lectura_vx1740_ext_nim_1M_v2.root",
    "lectura_vx1740_ext_nim_1M_v3.root",
    "lectura_vx1740_ext_nim_1M_v4.root",
    "lectura_vx1740_ext_nim_1M_v5.root",
    "lectura_vx1740_ext_nim_1M_v6.root",
    "lectura_vx1740_ext_nim_1M_v7.root",
    "lectura_vx1740_ext_nim_1M_v8.root",
    "lectura_vx1740_ext_nim_1M_v9.root"
};

struct Window { int first = 0; int last = -1; };
struct RunData {
    long long entries = 0;
    std::array<long long, NCH> down{};
    std::array<long long, NCH> up{};
    std::array<std::unique_ptr<TH1D>, NCH> amp_hist;
    std::array<long long, NCH> muon_hits{};
    std::array<long long, NCH> accidental_hits{};
    std::array<long long, NCH> zero_sigma_events{};
    std::array<std::vector<double>, NCH> muon_amp;
};
struct FitResult {
    double mpv = std::numeric_limits<double>::quiet_NaN();
    double mpv_err = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double sigma_err = std::numeric_limits<double>::quiet_NaN();
    long long n = 0;
    int status = -1;
};
struct ChannelResult {
    std::string verdict;
    std::string reason;
    double false_rate = 0;
    double noise_z = 0;
    double zero_sigma_fraction = 0;
    double slope = 0;
    double slope_z = 0;
    bool one_pe = false;
    std::array<double, NRUN> dcr{};
    std::array<double, NRUN> dcr_error{};
    std::array<double, NRUN> accidental{};
    std::array<double, NRUN> occupancy{};
    bool active_v7 = false;
    bool active_v8 = false;
    FitResult fit_v7;
    FitResult fit_v8;
};

static std::string detector(int ch) { return ch < 32 ? "hod2018" : "hod2019"; }
static std::string axis(int ch) {
    if (ch < 16) return "X";
    if (ch < 32) return "Y";
    for (int k = 0; k < 16; ++k) {
        if (H19X[k] == ch) return "X";
        if (H19Y[k] == ch) return "Y";
    }
    return "?";
}
static int physical_k(int ch) {
    if (ch < 16) return ch;
    if (ch < 32) return ch - 16;
    for (int k = 0; k < 16; ++k) {
        if (H19X[k] == ch || H19Y[k] == ch) return k;
    }
    return -1;
}

static Window find_window(const std::vector<float>& time, double low, double high) {
    auto lo = std::lower_bound(time.begin(), time.end(), static_cast<float>(low));
    auto hi = std::upper_bound(time.begin(), time.end(), static_cast<float>(high));
    if (lo == time.end() || hi == time.begin() || lo >= hi) return {};
    return {static_cast<int>(lo - time.begin()), static_cast<int>(hi - time.begin()) - 1};
}

static double median(std::vector<double>& v) {
    const size_t n = v.size();
    if (!n) return 0;
    auto mid = v.begin() + n / 2;
    std::nth_element(v.begin(), mid, v.end());
    double m = *mid;
    if ((n & 1U) == 0) {
        auto lower = std::max_element(v.begin(), mid);
        m = 0.5 * (m + *lower);
    }
    return m;
}

static std::pair<long long, std::vector<double>> count_excursions(
    const std::vector<float>& wf, Window w, double pedestal, double cut,
    bool negative, TH1D* hist) {
    long long count = 0;
    std::vector<double> amplitudes;
    int i = w.first;
    while (i <= w.last) {
        const double excursion = negative ? pedestal - wf[i] : wf[i] - pedestal;
        if (excursion <= cut) { ++i; continue; }
        const int start = i;
        double peak = excursion;
        while (i <= w.last) {
            const double a = negative ? pedestal - wf[i] : wf[i] - pedestal;
            if (a <= cut) break;
            peak = std::max(peak, a);
            ++i;
        }
        if (i - start >= 2) {
            ++count;
            amplitudes.push_back(peak);
            if (hist) hist->Fill(peak);
        }
        if (i == start) ++i;
    }
    return {count, std::move(amplitudes)};
}

static double max_ph(const std::vector<float>& wf, Window w, double pedestal) {
    double peak = 0;
    for (int i = w.first; i <= w.last; ++i) peak = std::max(peak, pedestal - wf[i]);
    return peak;
}

static uint32_t png_crc32(const unsigned char* data, size_t size) {
    uint32_t crc = 0xffffffffU;
    for (size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) crc = (crc >> 1) ^ (0xedb88320U & (-(crc & 1U)));
    }
    return crc ^ 0xffffffffU;
}

static void set_png_300_dpi(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return;
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)), {});
    if (bytes.size() < 33 || std::string(bytes.begin() + 12, bytes.begin() + 16) != "IHDR") return;
    const uint32_t ppm = 11811;
    std::vector<unsigned char> chunk = {0,0,0,9,'p','H','Y','s',
        static_cast<unsigned char>(ppm >> 24), static_cast<unsigned char>(ppm >> 16),
        static_cast<unsigned char>(ppm >> 8), static_cast<unsigned char>(ppm),
        static_cast<unsigned char>(ppm >> 24), static_cast<unsigned char>(ppm >> 16),
        static_cast<unsigned char>(ppm >> 8), static_cast<unsigned char>(ppm), 1};
    const uint32_t crc = png_crc32(chunk.data() + 4, 13);
    chunk.push_back(static_cast<unsigned char>(crc >> 24));
    chunk.push_back(static_cast<unsigned char>(crc >> 16));
    chunk.push_back(static_cast<unsigned char>(crc >> 8));
    chunk.push_back(static_cast<unsigned char>(crc));
    bytes.insert(bytes.begin() + 33, chunk.begin(), chunk.end());
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
}

static bool analyze_run(const fs::path& filename, int run_index, RunData& out) {
    const std::string resolved = hodo::resolveRawFile(filename.string());
    std::unique_ptr<TFile> file(TFile::Open(resolved.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open " << resolved << "\n";
        return false;
    }
    auto* tree = dynamic_cast<TTree*>(file->Get("data_dig"));
    if (!tree) { std::cerr << "ERROR: missing data_dig in " << filename << "\n"; return false; }
    out.entries = tree->GetEntries();
    tree->SetBranchStatus("*", 0);
    tree->SetCacheSize(96LL * 1024 * 1024);
    std::array<std::vector<float>*, NCH> channels{};
    std::vector<float>* time = nullptr;
    for (int ch = 0; ch < NCH; ++ch) {
        const std::string name = "Ch" + std::to_string(ch);
        if (!tree->GetBranch(name.c_str())) { std::cerr << "ERROR: missing " << name << "\n"; return false; }
        tree->SetBranchStatus(name.c_str(), 1);
        if (tree->SetBranchAddress(name.c_str(), &channels[ch]) < 0) return false;
    }
    if (!tree->GetBranch("time")) { std::cerr << "ERROR: missing time\n"; return false; }
    tree->SetBranchStatus("time", 1);
    if (tree->SetBranchAddress("time", &time) < 0) return false;

    std::array<std::string, NCH> hist_names;
    for (int ch = 0; ch < NCH; ++ch) {
        hist_names[ch] = "preamp_" + RUN_NAMES[run_index] + "_ch" + std::to_string(ch);
        out.amp_hist[ch] = std::make_unique<TH1D>(hist_names[ch].c_str(),
            ("Pretrigger pulse amplitudes " + RUN_NAMES[run_index] + " Ch" + std::to_string(ch) + ";Amplitude [ADC];Pulses").c_str(),
            200, 0, 400);
        out.amp_hist[ch]->SetDirectory(nullptr);
    }

    Window pre, muon, accidental;
    std::vector<double> base_samples, deviations;
    for (long long entry = 0; entry < out.entries; ++entry) {
        tree->GetEntry(entry);
        if (!time || time->empty()) continue;
        if (entry == 0) {
            pre = find_window(*time, 200, 3400);
            muon = find_window(*time, 4000, 4400);
            accidental = find_window(*time, 6000, 6400);
            if (pre.last < pre.first || muon.last < muon.first || accidental.last < accidental.first) {
                std::cerr << "ERROR: a requested time window is absent in " << filename << "\n";
                return false;
            }
        }
        for (int ch = 0; ch < NCH; ++ch) {
            if (!channels[ch] || channels[ch]->empty()) continue;
            const auto& wf = *channels[ch];
            const int last = std::min<int>(pre.last, wf.size() - 1);
            const int first = std::min<int>(pre.first, last);
            base_samples.clear();
            base_samples.reserve(last - first + 1);
            for (int s = first; s <= last; ++s) base_samples.push_back(wf[s]);
            const double ped = median(base_samples);
            deviations.clear();
            deviations.reserve(base_samples.size());
            for (double sample : base_samples) deviations.push_back(std::abs(sample - ped));
            const double sigma = 1.4826 * median(deviations);
            if (sigma == 0.0) ++out.zero_sigma_events[ch];
            const double cut = 4.0 * sigma;
            Window pre_ch{first, std::min<int>(pre.last, wf.size() - 1)};
            Window mu_ch{std::min<int>(muon.first, static_cast<int>(wf.size()) - 1), std::min<int>(muon.last, static_cast<int>(wf.size()) - 1)};
            Window acc_ch{std::min<int>(accidental.first, static_cast<int>(wf.size()) - 1), std::min<int>(accidental.last, static_cast<int>(wf.size()) - 1)};
            auto down = count_excursions(wf, pre_ch, ped, cut, true, out.amp_hist[ch].get());
            auto up = count_excursions(wf, pre_ch, ped, cut, false, nullptr);
            out.down[ch] += down.first;
            out.up[ch] += up.first;

            if (run_index == 6 || run_index == 7) {
                const double threshold = ch < 32 ? 50.0 : 35.0;
                const double mu_ph = max_ph(wf, mu_ch, ped);
                const double acc_ph = max_ph(wf, acc_ch, ped);
                if (mu_ph > threshold) {
                    ++out.muon_hits[ch];
                    out.muon_amp[ch].push_back(mu_ph);
                }
                if (acc_ph > threshold) ++out.accidental_hits[ch];
            }
        }
        if (entry && entry % 20000 == 0) {
            std::cout << RUN_NAMES[run_index] << " " << entry << "/" << out.entries << "\r" << std::flush;
        }
    }
    std::cout << RUN_NAMES[run_index] << " processed " << out.entries << " entries\n";
    return true;
}

static FitResult fit_landau(const std::vector<double>& values, double threshold, const std::string& name) {
    FitResult result;
    result.n = values.size();
    if (values.size() < 25) return result;
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    const double q995 = sorted[std::min(sorted.size() - 1, static_cast<size_t>(0.995 * sorted.size()))];
    const double xmax = std::max(threshold * 1.25, q995 * 1.15);
    TH1D hist((name + "_hist").c_str(), "", 160, 0, xmax);
    hist.SetDirectory(nullptr);
    for (double v : values) if (v < xmax) hist.Fill(v);
    TF1 fit((name + "_landau").c_str(), "landau", threshold, xmax);
    fit.SetParameters(hist.GetMaximum(), std::max(threshold, hist.GetXaxis()->GetBinCenter(hist.GetMaximumBin())), 8.0);
    auto fit_ptr = hist.Fit(&fit, "LQ0SNR");
    result.status = static_cast<int>(fit_ptr);
    if (fit_ptr.Get() && fit_ptr->IsValid()) {
        result.mpv = fit.GetParameter(1);
        result.mpv_err = fit.GetParError(1);
        result.sigma = fit.GetParameter(2);
        result.sigma_err = fit.GetParError(2);
    }
    return result;
}

static bool has_one_pe_peak(const TH1D& h) {
    if (h.GetEntries() < 20) return false;
    std::vector<double> smoothed(h.GetNbinsX() + 1, 0);
    for (int b = 2; b < h.GetNbinsX(); ++b)
        smoothed[b] = (h.GetBinContent(b - 1) + 2 * h.GetBinContent(b) + h.GetBinContent(b + 1)) / 4.0;
    const double needed = std::max(10.0, 0.02 * h.GetEntries());
    for (int b = 2; b < h.GetNbinsX(); ++b) {
        const double amplitude = h.GetXaxis()->GetBinCenter(b);
        const double neighbor = std::max(smoothed[b - 1], smoothed[b + 1]);
        if (amplitude >= 8.0 && smoothed[b] >= needed && smoothed[b] >= 1.5 * neighbor &&
            smoothed[b] > smoothed[b - 1] && smoothed[b] >= smoothed[b + 1]) return true;
    }
    return false;
}

static std::pair<double, double> weighted_slope(const std::array<double, NRUN>& y,
                                                const std::array<double, NRUN>& err) {
    double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < NRUN; ++i) {
        const double x = i + 1.0;
        const double e = std::max(err[i], 1e-12);
        const double w = 1.0 / (e * e);
        sw += w; sx += w * x; sy += w * y[i]; sxx += w * x * x; sxy += w * x * y[i];
    }
    const double delta = sw * sxx - sx * sx;
    if (delta <= 0) return {0, 0};
    const double slope = (sw * sxy - sx * sy) / delta;
    const double slope_err = std::sqrt(sw / delta);
    return {slope, slope_err > 0 ? slope / slope_err : 0};
}

static std::string trend_text(double slope, double z) {
    const char* sign = slope > 0 ? "increasing" : (slope < 0 ? "decreasing" : "flat");
    std::ostringstream os;
    os << sign << " (z=" << std::fixed << std::setprecision(2) << z << ")";
    return os.str();
}

static std::string verdict_color(const std::string& verdict) {
    return verdict;
}

static void draw_figures(const fs::path& outdir, const std::array<RunData, NRUN>& runs,
                         const std::array<ChannelResult, NCH>& result) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    constexpr int colors[] = {kGreen + 2, kOrange + 7, kRed + 1, kGray + 2};
    constexpr const char* verdicts[] = {"active", "live SiPM, no muon light", "no SiPM signal (electrical suspect)", "ambiguous"};
    TCanvas c1("c_dcr_channel", "", 1204, 928);
    c1.SetLeftMargin(0.12); c1.SetBottomMargin(0.13); c1.SetRightMargin(0.04);
    TH1D frame("dcr_channel_frame", "DCR per channel, v8;Channel;DCR [Hz]", 64, -0.5, 63.5);
    frame.SetMinimum(1e5); frame.SetMaximum(1e7); frame.SetStats(0);
    frame.GetXaxis()->SetLabelSize(0.027); frame.GetYaxis()->SetTitleOffset(1.25);
    frame.Draw(); c1.SetLogy();
    std::array<std::unique_ptr<TGraph>, 4> groups;
    for (auto& g : groups) g = std::make_unique<TGraph>();
    for (int ch = 0; ch < NCH; ++ch) {
        int group = result[ch].verdict == verdicts[0] ? 0 : result[ch].verdict == verdicts[1] ? 1 : result[ch].verdict == verdicts[2] ? 2 : 3;
        double y = std::max(result[ch].dcr[7], 0.01);
        groups[group]->SetPoint(groups[group]->GetN(), ch, y);
    }
    TLegend leg(0.60, 0.68, 0.96, 0.94); leg.SetBorderSize(0); leg.SetFillStyle(0);
    for (int i = 0; i < 4; ++i) {
        groups[i]->SetMarkerStyle(20); groups[i]->SetMarkerSize(1.15); groups[i]->SetMarkerColor(colors[i]);
        groups[i]->SetLineColor(colors[i]);
        auto* entry = leg.AddEntry(groups[i].get(), verdicts[i], "p");
        entry->SetMarkerColor(colors[i]); entry->SetMarkerStyle(20);
        if (groups[i]->GetN() > 0) groups[i]->Draw("P SAME");
    }
    for (double x : {15.5, 31.5}) { TLine line(x, 1e5, x, 1e7); line.SetLineStyle(2); line.SetLineColor(kBlack); line.Draw(); }
    leg.Draw(); c1.SaveAs((outdir / "dcr_vs_channel_v8.png").string().c_str());

    TCanvas c2("c_dcr_run", "", 1204, 928);
    c2.SetLeftMargin(0.12); c2.SetBottomMargin(0.12); c2.SetRightMargin(0.04); c2.SetLogy();
    TH1D fr2("dcr_run_frame", "hod2019 DCR vs run;Run;DCR [Hz]", 9, 0.5, 9.5);
    fr2.SetMinimum(1e4); fr2.SetMaximum(1e7); fr2.SetStats(0);
    for (int r = 0; r < NRUN; ++r) fr2.GetXaxis()->SetBinLabel(r + 1, RUN_NAMES[r].c_str());
    fr2.GetXaxis()->SetLabelSize(0.045); fr2.Draw();
    std::array<std::unique_ptr<TGraph>, 32> curves;
    for (int i = 0; i < 32; ++i) {
        const int ch = 32 + i;
        curves[i] = std::make_unique<TGraph>(NRUN);
        for (int r = 0; r < NRUN; ++r) curves[i]->SetPoint(r, r + 1, std::max(result[ch].dcr[r], 0.01));
        curves[i]->SetLineColor(kAzure + (i % 9) + 1); curves[i]->SetMarkerColor(kAzure + (i % 9) + 1);
        curves[i]->SetMarkerStyle(result[ch].active_v8 ? 20 : 24);
        curves[i]->SetLineStyle(result[ch].active_v8 ? 1 : 2); curves[i]->SetLineWidth(1);
        curves[i]->Draw("LP SAME");
    }
    TLegend l2(0.72, 0.82, 0.94, 0.94); l2.SetBorderSize(0); l2.SetFillStyle(0);
    auto active_key = new TGraph(); active_key->SetMarkerStyle(20); active_key->SetLineStyle(1); active_key->SetMarkerColor(kBlack); active_key->SetLineColor(kBlack);
    auto silent_key = new TGraph(); silent_key->SetMarkerStyle(24); silent_key->SetLineStyle(2); silent_key->SetMarkerColor(kBlack); silent_key->SetLineColor(kBlack);
    l2.AddEntry(active_key, "Muon-active", "lp"); l2.AddEntry(silent_key, "Muon-silent", "lp"); l2.Draw();
    c2.SaveAs((outdir / "dcr_vs_run.png").string().c_str());
    delete active_key; delete silent_key;

    TCanvas c3("c_mpv", "", 1204, 928);
    c3.SetLeftMargin(0.12); c3.SetBottomMargin(0.12); c3.SetRightMargin(0.04);
    TH1D fr3("mpv_frame", "Muon-active strip pulse-height MPV, v8;Physical strip k;Landau MPV [ADC]", 16, -0.5, 15.5);
    fr3.SetMinimum(0); fr3.SetStats(0); fr3.GetXaxis()->SetTitleOffset(1.15); fr3.Draw();
    TGraphErrors h18, h19;
    int n18 = 0, n19 = 0;
    for (int ch = 0; ch < NCH; ++ch) {
        if (!result[ch].active_v8 || !std::isfinite(result[ch].fit_v8.mpv)) continue;
        const int k = physical_k(ch);
        if (ch < 32) h18.SetPoint(n18, k - 0.12, result[ch].fit_v8.mpv), h18.SetPointError(n18++, 0, result[ch].fit_v8.mpv_err);
        else h19.SetPoint(n19, k + 0.12, result[ch].fit_v8.mpv), h19.SetPointError(n19++, 0, result[ch].fit_v8.mpv_err);
    }
    h18.SetMarkerStyle(21); h18.SetMarkerColor(kBlue + 1); h18.SetLineColor(kBlue + 1); h18.SetMarkerSize(1.2);
    h19.SetMarkerStyle(22); h19.SetMarkerColor(kRed + 1); h19.SetLineColor(kRed + 1); h19.SetMarkerSize(1.2);
    h18.Draw("P SAME"); h19.Draw("P SAME");
    TLegend l3(0.72, 0.82, 0.94, 0.94); l3.SetBorderSize(0); l3.SetFillStyle(0);
    l3.AddEntry(&h18, "hod2018", "p"); l3.AddEntry(&h19, "hod2019", "p"); l3.Draw();
    c3.SaveAs((outdir / "mpv_per_channel_v8.png").string().c_str());

    std::vector<int> silent, active;
    for (int ch : H19X) if (!result[ch].active_v8) silent.push_back(ch);
    for (int ch : H19Y) if (!result[ch].active_v8 && silent.size() < 2) silent.push_back(ch);
    for (int ch : H19X) if (result[ch].active_v8) active.push_back(ch);
    for (int ch : H19Y) if (result[ch].active_v8 && active.size() < 2) active.push_back(ch);
    TCanvas c4("c_preamp_examples", "", 1204, 928); c4.Divide(2, 2, 0.001, 0.001);
    for (int pad = 0; pad < 4; ++pad) {
        c4.cd(pad + 1); gPad->SetLogy();
        const bool is_silent = pad < 2;
        const int ch = is_silent ? (pad < static_cast<int>(silent.size()) ? silent[pad] : -1)
                                 : (pad - 2 < static_cast<int>(active.size()) ? active[pad - 2] : -1);
        if (ch < 0) continue;
        TH1D* h = runs[7].amp_hist[ch].get();
        h->SetTitle((std::string(is_silent ? "Muon-silent" : "Muon-active") + " hod2019 Ch" + std::to_string(ch) + ";Pulse amplitude [ADC];Pulses").c_str());
        double probability = 0.995, q = 0;
        h->GetQuantiles(1, &q, &probability);
        h->GetXaxis()->SetRangeUser(0, std::max(20.0, 1.2 * q));
        h->SetLineColor(is_silent ? kRed + 1 : kBlue + 1); h->Draw("HIST");
    }
    c4.SaveAs((outdir / "pretrigger_amplitude_examples.png").string().c_str());
    for (const char* name : {"dcr_vs_channel_v8.png", "dcr_vs_run.png", "mpv_per_channel_v8.png", "pretrigger_amplitude_examples.png"})
        set_png_300_dpi(outdir / name);
}

static void write_cable_comparison(const fs::path& path, const std::array<ChannelResult, NCH>& result) {
    std::ofstream out(path);
    out << "Cable report comparison\nObserved set means verdict = no SiPM signal (electrical suspect).\n";
    out << "Cable pins: Lateral OP AMP 4 and 8; Lateral Opuesta OP AMP 1.\n";
    out << "Consistency means exact equality between predicted and observed k sets.\n\n";
    constexpr std::array<int, 3> amps = {4, 8, 1};
    constexpr std::array<int, 3> lateral = {1, 1, 0};
    for (const std::string plane : {"hod2018_X", "hod2018_Y", "hod2019_X", "hod2019_Y"}) {
        const bool h19 = plane.find("hod2019") == 0;
        const bool is_x = plane.back() == 'X';
        std::vector<int> observed;
        for (int k = 0; k < 16; ++k) {
            int ch = h19 ? (is_x ? H19X[k] : H19Y[k]) : (is_x ? k : 16 + k);
            if (result[ch].verdict == "no SiPM signal (electrical suspect)") observed.push_back(k);
        }
        auto list = [](const std::vector<int>& v) {
            std::ostringstream s; s << "{";
            for (size_t i = 0; i < v.size(); ++i) { if (i) s << ","; s << v[i]; }
            s << "}"; return s.str();
        };
        out << plane << " observed electrical-suspect k=" << list(observed) << "\n";
        std::vector<std::string> rules = {"mirrored_offsets"};
        if (!h19) rules.push_back("identity_offsets");
        for (const auto& rule : rules) {
            for (int lateral_face_a = 1; lateral_face_a >= 0; --lateral_face_a) {
                std::vector<int> predicted;
                for (int p = 0; p < 3; ++p) {
                    int offset = rule == "mirrored_offsets" ? std::array<int, 8>{5,6,7,4,3,2,1,0}[amps[p] - 1] : amps[p] - 1;
                    int face_a = lateral[p] ? lateral_face_a : !lateral_face_a;
                    predicted.push_back(2 * offset + (face_a ? 0 : 1));
                }
                std::sort(predicted.begin(), predicted.end());
                const bool consistent = predicted == observed;
                std::vector<int> predicted_missing, observed_unpredicted;
                for (int k : predicted) if (std::find(observed.begin(), observed.end(), k) == observed.end()) predicted_missing.push_back(k);
                for (int k : observed) if (std::find(predicted.begin(), predicted.end(), k) == predicted.end()) observed_unpredicted.push_back(k);
                out << "  face assignment Lateral=" << (lateral_face_a ? "A" : "B") << ", rule=" << rule
                    << ", predicted=" << list(predicted) << ": " << (consistent ? "CONSISTENT" : "EXCLUDED");
                if (!consistent) out << "; predicted without observed support=" << list(predicted_missing)
                                     << "; observed strips not predicted=" << list(observed_unpredicted);
                out << "\n";
            }
        }
        if (h19) out << "  identity_offsets was not specified as a hod2019 candidate rule.\n";
        out << "\n";
    }
}

int main(int argc, char** argv) {
    const fs::path outdir = argc > 1 ? argv[1] : "outputs/channel_health_20261016";
    fs::create_directories(outdir);
    std::array<RunData, NRUN> runs;
    for (int r = 0; r < NRUN; ++r) {
        if (!analyze_run(RUN_FILES[r], r, runs[r])) return 2;
    }

    std::array<ChannelResult, NCH> result;
    const int v7 = 6, v8 = 7;
    for (int ch = 0; ch < NCH; ++ch) {
        auto& cr = result[ch];
        long long down_total = 0, up_total = 0;
        for (int r = 0; r < NRUN; ++r) {
            const double exposure = runs[r].entries * PRE_SECONDS;
            cr.dcr[r] = exposure > 0 ? runs[r].down[ch] / exposure : 0;
            cr.dcr_error[r] = exposure > 0 ? std::sqrt(std::max<long long>(runs[r].down[ch], 1)) / exposure : 0;
            cr.accidental[r] = runs[r].entries > 0 ? static_cast<double>(runs[r].accidental_hits[ch]) / runs[r].entries : 0;
            cr.occupancy[r] = runs[r].entries > 0 ? static_cast<double>(runs[r].muon_hits[ch]) / runs[r].entries : 0;
            down_total += runs[r].down[ch]; up_total += runs[r].up[ch];
        }
        const long long n_events = std::accumulate(runs.begin(), runs.end(), 0LL,
            [](long long s, const RunData& d) { return s + d.entries; });
        const double exposure = n_events * PRE_SECONDS;
        long long zero_sigma = 0;
        for (const auto& run : runs) zero_sigma += run.zero_sigma_events[ch];
        cr.zero_sigma_fraction = n_events > 0 ? static_cast<double>(zero_sigma) / n_events : 0;
        const double expected = exposure > 0 ? up_total / exposure : 0;
        cr.false_rate = expected;
        cr.noise_z = down_total + up_total > 0 ? (down_total - up_total) / std::sqrt(static_cast<double>(down_total + up_total)) : 0;
        auto [slope, z] = weighted_slope(cr.dcr, cr.dcr_error);
        cr.slope = slope; cr.slope_z = z;
        cr.active_v7 = runs[v7].entries && cr.occupancy[v7] > 3.0 * cr.accidental[v7];
        cr.active_v8 = runs[v8].entries && cr.occupancy[v8] > 3.0 * cr.accidental[v8];
        const double threshold = ch < 32 ? 50 : 35;
        cr.fit_v7 = fit_landau(runs[v7].muon_amp[ch], threshold, "fit_v7_ch" + std::to_string(ch));
        cr.fit_v8 = fit_landau(runs[v8].muon_amp[ch], threshold, "fit_v8_ch" + std::to_string(ch));
        cr.one_pe = has_one_pe_peak(*runs[v8].amp_hist[ch]);
        if (cr.active_v8) {
            cr.verdict = "active";
            cr.reason = "v8 muon occupancy exceeds 3x the equal-width accidental occupancy";
        } else if (cr.noise_z <= 2.0 && std::abs(cr.slope_z) < 2.0) {
            cr.verdict = "no SiPM signal (electrical suspect)";
            cr.reason = "downward pretrigger excursions consistent with upward control; no significant DCR trend";
        } else if (cr.noise_z > 3.0 && cr.slope_z > 2.0) {
            cr.verdict = "live SiPM, no muon light";
            cr.reason = "downward pretrigger excess over upward control and significant increasing DCR";
        } else {
            cr.verdict = "ambiguous";
            cr.reason = "criteria not jointly met: noise_z=" + std::to_string(cr.noise_z) +
                        ", slope_z=" + std::to_string(cr.slope_z) +
                        ", sigma_zero_fraction=" + std::to_string(100.0 * cr.zero_sigma_fraction) + "%";
        }
    }

    std::ofstream csv(outdir / "channel_verdicts.csv");
    csv << "ch,hodoscope,axis,k,muon_active,occupancy_v8,occupancy_v7,accidental_v8,MPV_v8,MPV_v8_error,sigma_v8,sigma_v8_error,N_entries_v8,MPV_v7,MPV_v7_error,sigma_v7,N_entries_v7,DCR_v8_Hz,DCR_v1_Hz,DCR_v2_Hz,DCR_v3_Hz,DCR_v4_Hz,DCR_v5_Hz,DCR_v6_Hz,DCR_v7_Hz,DCR_v9_Hz,DCR_trend_slope_Hz_per_rank,DCR_trend_z,false_pulse_rate_Hz,noise_excess_z,sigma_zero_fraction_v1_v9,1pe_peak_seen,verdict,reason\n";
    for (int ch = 0; ch < NCH; ++ch) {
        const auto& r = result[ch];
        csv << ch << ',' << detector(ch) << ',' << axis(ch) << ',' << physical_k(ch) << ',' << (r.active_v8 ? "y" : "n") << ','
            << r.occupancy[7] << ',' << r.occupancy[6] << ',' << r.accidental[7] << ','
            << r.fit_v8.mpv << ',' << r.fit_v8.mpv_err << ',' << r.fit_v8.sigma << ',' << r.fit_v8.sigma_err << ',' << r.fit_v8.n << ','
            << r.fit_v7.mpv << ',' << r.fit_v7.mpv_err << ',' << r.fit_v7.sigma << ',' << r.fit_v7.n << ','
            << r.dcr[7];
        for (int run : {0,1,2,3,4,5,6,8}) csv << ',' << r.dcr[run];
        csv << ',' << r.slope << ',' << r.slope_z << ',' << r.false_rate << ',' << r.noise_z << ',' << r.zero_sigma_fraction << ','
            << (r.one_pe ? "y" : "n") << ',' << '"' << r.verdict << "\",\"" << r.reason << "\"\n";
    }

    std::ofstream summary(outdir / "summary_channel_health.txt");
    summary << std::fixed << std::setprecision(6);
    summary << "Channel health analysis\n=======================\n\nInput files and entries (v10 excluded):\n";
    for (int r = 0; r < NRUN; ++r) summary << "  " << fs::absolute(RUN_FILES[r]) << "  " << runs[r].entries << "\n";
    summary << "\nWindows: pretrigger 200-3400 ns; muon 4000-4400 ns; accidental 6000-6400 ns.\n"
            << "DCR exposure per event uses 3.2e-6 s. Pedestal=median and sigma=MAD*1.4826 in the pretrigger window.\n"
            << "Pulse = contiguous excursion beyond 4 sigma with at least 2 samples.\n"
            << "Muon hit = maximum pedestal-subtracted negative amplitude above 50 ADC (hod2018) or 35 ADC (hod2019).\n"
            << "Muon active when occupancy > 3x accidental occupancy.\n"
            << "False-pulse control uses upward excursions in the same pretrigger samples.\n"
            << "Events with sigma_noise=0 are counted using the literal zero threshold and tracked explicitly below; such events do not provide a resolved noise scale.\n"
            << "Trend is weighted linear fit of raw DCR against voltage rank v1..v9; z= slope/error; v9 is approximately v8.\n"
            << "Noise excess z=(N_down-N_up)/sqrt(N_down+N_up), pooled over v1..v9.\n"
            << "1pe_peak_seen is an operational candidate flag: >=20 pulses and a smoothed local maximum at amplitude >=8 ADC, at least 1.5x either neighbor, and >=max(10 counts, 2% of entries); ADC-to-PE calibration is unavailable.\n"
            << "Landau MPV and sigma are ROOT landau fits to event peak amplitude above threshold; N_entries counts amplitudes above threshold.\n\n";
    for (const std::string hodo : {"hod2018", "hod2019"}) {
        for (const std::string ax : {"X", "Y"}) {
            std::array<int, 4> counts{};
            for (int ch = 0; ch < NCH; ++ch) if (detector(ch) == hodo && axis(ch) == ax) {
                const auto& v = result[ch].verdict;
                if (v == "active") ++counts[0]; else if (v == "live SiPM, no muon light") ++counts[1];
                else if (v == "no SiPM signal (electrical suspect)") ++counts[2]; else ++counts[3];
            }
            summary << hodo << ' ' << ax << ": active=" << counts[0] << ", live_SiPM_no_muon=" << counts[1]
                    << ", electrical_suspect=" << counts[2] << ", ambiguous=" << counts[3] << "\n";
        }
    }
    int h19_live_silent = 0, h19_dead_silent = 0, h19_ambig_silent = 0;
    for (int ch = 32; ch < 64; ++ch) if (!result[ch].active_v8) {
        if (result[ch].verdict == "live SiPM, no muon light") ++h19_live_silent;
        else if (result[ch].verdict == "no SiPM signal (electrical suspect)") ++h19_dead_silent;
        else ++h19_ambig_silent;
    }
    summary << "\nhod2019 muon-silent strips: live-SiPM=" << h19_live_silent << ", electrical-suspect=" << h19_dead_silent
            << ", ambiguous=" << h19_ambig_silent << "\n";
    long long zero_sigma_total = 0, channel_event_total = 0;
    for (const auto& run : runs) for (int ch = 0; ch < NCH; ++ch) {
        zero_sigma_total += run.zero_sigma_events[ch];
        channel_event_total += run.entries;
    }
    summary << "Pretrigger channel-event windows with sigma_noise=0: " << zero_sigma_total << "/" << channel_event_total
            << " (" << (channel_event_total ? 100.0 * zero_sigma_total / channel_event_total : 0.0) << "%)\n";
    std::vector<double> h18_pool, h19_pool;
    for (int ch = 0; ch < NCH; ++ch) if (result[ch].active_v8) {
        const auto& values = runs[7].muon_amp[ch];
        auto& pool = ch < 32 ? h18_pool : h19_pool;
        pool.insert(pool.end(), values.begin(), values.end());
    }
    auto pooled18 = fit_landau(h18_pool, 50, "pooled_hod2018_v8");
    auto pooled19 = fit_landau(h19_pool, 35, "pooled_hod2019_v8");
    summary << "\nPooled v8 Landau fits over muon-active strips:\n"
            << "  hod2018: MPV=" << pooled18.mpv << " +/- " << pooled18.mpv_err << " ADC, sigma=" << pooled18.sigma << " +/- " << pooled18.sigma_err << ", N=" << pooled18.n << ", fit_status=" << pooled18.status << "\n"
            << "  hod2019: MPV=" << pooled19.mpv << " +/- " << pooled19.mpv_err << " ADC, sigma=" << pooled19.sigma << " +/- " << pooled19.sigma_err << ", N=" << pooled19.n << ", fit_status=" << pooled19.status << "\n";
    long long all_down = 0, all_up = 0;
    for (const auto& run : runs) for (int ch = 0; ch < NCH; ++ch) { all_down += run.down[ch]; all_up += run.up[ch]; }
    const double total_events = std::accumulate(runs.begin(), runs.end(), 0.0, [](double s, const RunData& d){ return s + d.entries; });
    summary << "\nFalse-pulse control from upward polarity: N_up=" << all_up << ", N_down=" << all_down
            << ", pooled upward rate per channel=" << all_up / (NCH * total_events * PRE_SECONDS) << " Hz"
            << ", all-channel sum=" << all_up / (total_events * PRE_SECONDS) << " Hz.\n";
    summary << "No physics cause is inferred from these channel classifications.\n";

    write_cable_comparison(outdir / "cable_report_comparison.txt", result);
    draw_figures(outdir, runs, result);
    TFile hist_file((outdir / "pretrigger_amplitude_hists.root").string().c_str(), "RECREATE");
    for (auto& run : runs) for (auto& hist : run.amp_hist) hist->Write();
    hist_file.Close();
    std::cout << "Output directory: " << fs::absolute(outdir) << "\n";
    return 0;
}
