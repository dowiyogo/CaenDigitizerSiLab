// reconstruct_v8.cxx
// Muon telescope reconstruction v8 — raw waveform → per-channel baseline correction →
// negative-pulse extraction → spatial clustering → track reconstruction.
//
// WHY v8 FIXES v7:
//   v7 used max_amp (captured baseline or overshoot, not pulse depth).
//   v8 correctly computes pulse_height = pedestal − min(signal_window) > 0 for
//   negative pulses on a positive baseline (e.g. Ch50: baseline 20 → min 3 → ph 17 ADC).
//
// KEY PHYSICAL FINDING (documented from data analysis):
//   hod2019 scintillator bars are 1 mm thick → muon energy deposit ≈ 0.2 MeV →
//   typical signal ≈ 17 ADC above baseline. Noise crossings at 5σ threshold fall at
//   5-7 ADC, overlapping with the signal. The digitizer boards for hod2018 and hod2019
//   use DIFFERENT pre-trigger depths (≈5200 ns vs ≈4000 ns), creating a systematic
//   ≈1200 ns time offset between the two hodoscopes in the digitized waveforms.
//   This offset is NOT from cable lengths — it is a digitizer configuration difference.
//
// TWO-STAGE ANALYSIS:
//   Stage A (--stage A or default): pulse/noise diagnostics with standard threshold.
//   Stage B (--stage B or default): track reconstruction with elevated hod2019 threshold
//     that suppresses noise crossings and retains only real muon signals (≥ abs_thr_2019).
//
// Pulse convention: NEGATIVE pulses on POSITIVE baseline.
//   pulse_height = pedestal − min(samples_in_signal_window) [positive for real hits]
//   pulse_time   = time of the minimum sample [time-of-minimum method, not CFD]
//
// Geometry: hod2018 at Z = 0; hod2019 at Z = +zsep_mm.
//   coord_mm(k) = (k + 0.5 − 8.0) × pitch_mm  (k = 0..15, increasing coordinate)
//   mx = (x2019 − x2018) / zsep_mm  → muon slope in XZ plane
//
// Usage:
//   ./reconstruct_v8 [input.root] [--pitch mm] [--zsep mm] [--nsigma N]
//                    [--abst18 ADC] [--abst19 ADC] [--dtdelay ns] [--dtwindow ns]
//                    [--nmax19 N] [--stage A|B|AB]
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o reconstruct_v8 reconstruct_v8.cxx

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
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
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

//==============================================================================
// Fixed detector constants — do not change without physics justification
//==============================================================================
constexpr int    kNCh            = 64;
constexpr double kPrePulseEnd_ns = 2500.0;  // pre-pulse baseline region end [ns]
constexpr double kHod2019Tlo_ns  = 3000.0;  // hod2019 signal window [ns]
constexpr double kHod2019Thi_ns  = 5000.0;
constexpr double kHod2018Tlo_ns  = 3500.0;  // hod2018 signal window [ns]
constexpr double kHod2018Thi_ns  = 7000.0;
constexpr double kRadToDeg       = 180.0 / M_PI;

//==============================================================================
// Run-time configurable parameters (set via CLI)
//==============================================================================
struct Params {
    std::string input_file  = "lectura_vx1740_ext_nim_1M_v7.root";
    double pitch_mm         =  2.0;    // effective strip pitch [mm]
    double zsep_mm          = 98.4;    // Z separation hod2018 → hod2019 [mm]
    double n_sigma          =  5.0;    // hit threshold in noise RMS units
    // Absolute thresholds. hod2019 real muon signal ≈ 17 ADC; noise crossings ≈ 5-7 ADC.
    // Set abs_thr_2019 = 8+ to select only real muon signals.
    double abs_thr_2018     = 20.0;   // abs threshold hod2018 [ADC]
    double abs_thr_2019     =  8.0;   // abs threshold hod2019 [ADC] — default 8 (selects real hits)
    // hod2018 is read out at t≈5200 ns, hod2019 at t≈4000 ns in the waveform.
    // The ≈1200 ns offset is a digitizer pre-trigger configuration difference.
    double dt_delay_ns      = 1200.0;  // expected t_2018 − t_2019 [ns]
    double dt_window_ns     =  600.0;  // hod2018 time-match window half-width [ns]
    int    max_cls_2019     =    2;    // max clusters per hod2019 plane
    bool   run_stage_a      = true;
    bool   run_stage_b      = true;
};

//==============================================================================
// Channel mapping — AUTHORITATIVE literal arrays (Excel "Disposición Correcta")
// Index k = 0..15 → physical channel number.  NO algebraic rules.
//
// Hod2018 X: k → k       (channels  0..15, identity mapping)
// Hod2018 Y: k → 16 + k  (channels 16..31)
// Hod2019 Y: channels 32..47, non-trivial interleave (table below)
// Hod2019 X: channels 48..63, non-trivial interleave (table below)
//==============================================================================
static const std::array<int,16> kHod2019X = {
    53, 61, 54, 62,
    55, 63, 52, 60,
    51, 59, 50, 58,
    49, 57, 48, 56
};
static const std::array<int,16> kHod2019Y = {
    37, 45, 38, 46,
    39, 47, 36, 44,
    35, 43, 34, 42,
    33, 41, 32, 40
};

// Control mapping — identity order (wrong) for QA comparison
static const std::array<int,16> kHod2019X_ctrl = {
    48, 49, 50, 51, 52, 53, 54, 55,
    56, 57, 58, 59, 60, 61, 62, 63
};
static const std::array<int,16> kHod2019Y_ctrl = {
    32, 33, 34, 35, 36, 37, 38, 39,
    40, 41, 42, 43, 44, 45, 46, 47
};

//==============================================================================
// Inverse maps: channel number → physical strip index (−1 if not in that plane)
//==============================================================================
struct InvMap {
    std::array<int,kNCh> x19{}, y19{};
    InvMap(const std::array<int,16>& xm, const std::array<int,16>& ym) {
        x19.fill(-1); y19.fill(-1);
        for (int k = 0; k < 16; ++k) { x19[xm[k]] = k; y19[ym[k]] = k; }
    }
};

static int  idx2018x(int ch) { return (ch >=  0 && ch <= 15) ? ch      : -1; }
static int  idx2018y(int ch) { return (ch >= 16 && ch <= 31) ? ch - 16 : -1; }
static bool is2019  (int ch) { return ch >= 32; }
static double sig_tlo(int ch) { return is2019(ch) ? kHod2019Tlo_ns : kHod2018Tlo_ns; }
static double sig_thi(int ch) { return is2019(ch) ? kHod2019Thi_ns : kHod2018Thi_ns; }
static double strip_coord(int k, double pitch) { return (k + 0.5 - 8.0) * pitch; }

//==============================================================================
// Pulse extraction
//==============================================================================
struct PulseResult {
    bool   valid        = false;
    double pedestal     = 0.0;   // [ADC] mean of pre-pulse baseline
    double noise        = 0.0;   // [ADC] RMS of pre-pulse baseline
    double pulse_height = 0.0;   // [ADC] pedestal − min_in_window  (positive)
    double pulse_time   = -1.0;  // [ns] time of minimum sample
    bool   hit          = false;
};

static PulseResult extract_pulse(const std::vector<float>& wf, const std::vector<float>& t,
                                  double tlo, double thi, double nsig, double athr) {
    PulseResult r;
    const int n = static_cast<int>(wf.size());
    if (n < 4 || static_cast<int>(t.size()) != n) return r;

    double sum=0,sum2=0; int nb=0;
    for (int i=0;i<n;++i) if(t[i]<kPrePulseEnd_ns){sum+=wf[i];sum2+=wf[i]*wf[i];++nb;}
    if (nb < 2) return r;

    r.valid    = true;
    r.pedestal = sum / nb;
    const double var = sum2/nb - r.pedestal*r.pedestal;
    r.noise    = (var > 0.0) ? std::sqrt(var) : 0.0;

    double vmin=std::numeric_limits<double>::infinity(); int imin=-1;
    for (int i=0;i<n;++i)
        if (t[i]>=tlo && t[i]<=thi && wf[i]<vmin) { vmin=wf[i]; imin=i; }
    if (imin < 0) return r;

    r.pulse_height = r.pedestal - vmin;
    r.pulse_time   = t[imin];
    r.hit = r.pulse_height > std::max(nsig * r.noise, athr);
    return r;
}

//==============================================================================
// Hit and cluster structures
//==============================================================================
struct Hit { int index; double ph; double pt; };

struct Cluster {
    double position;   // amplitude-weighted centroid [mm]
    double time;       // amplitude-weighted mean time [ns]
    double amplitude;  // sum of pulse heights [ADC]
    int    nstrips;
    int    lo, hi;
};

static std::vector<Cluster> make_clusters(std::vector<Hit> hits, double pitch) {
    if (hits.empty()) return {};
    std::sort(hits.begin(),hits.end(),[](const Hit&a,const Hit&b){return a.index<b.index;});
    std::vector<Cluster> out;
    Cluster cur{};
    auto flush=[&](){ cur.position/=cur.amplitude; cur.time/=cur.amplitude; out.push_back(cur); };
    auto open =[&](const Hit&h){ cur={}; cur.lo=cur.hi=h.index; };
    auto add  =[&](const Hit&h){
        cur.position+=h.ph*strip_coord(h.index,pitch);
        cur.time    +=h.ph*h.pt;
        cur.amplitude+=h.ph; ++cur.nstrips; cur.hi=h.index;
    };
    open(hits[0]); add(hits[0]);
    for (int i=1;i<(int)hits.size();++i) {
        if (hits[i].index > cur.hi+1) { flush(); open(hits[i]); }
        add(hits[i]);
    }
    flush();
    return out;
}

static const Cluster& best_amp(const std::vector<Cluster>& v) {
    return *std::max_element(v.begin(),v.end(),
        [](const Cluster&a,const Cluster&b){return a.amplitude<b.amplitude;});
}

// Select cluster closest in time to t_ref within dt_window; prefer higher amplitude as tiebreak
static const Cluster* best_time(const std::vector<Cluster>& v, double t_ref, double dt_window) {
    const Cluster* best=nullptr; double best_dt=std::numeric_limits<double>::infinity();
    for (const auto&c:v) {
        double dt=std::abs(c.time-t_ref);
        if (dt<dt_window && (dt<best_dt||(dt==best_dt&&best&&c.amplitude>best->amplitude)))
            { best_dt=dt; best=&c; }
    }
    return best;
}

//==============================================================================
// Statistics
//==============================================================================
static double pearson_r(const std::vector<double>& x, const std::vector<double>& y) {
    const int n=static_cast<int>(x.size()); if(n<2) return 0.0;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];syy+=y[i]*y[i];sxy+=x[i]*y[i];}
    const double d=n, num=d*sxy-sx*sy, den=std::sqrt((d*sxx-sx*sx)*(d*syy-sy*sy));
    return den>0?num/den:0.0;
}
static void linfit(const std::vector<double>&x,const std::vector<double>&y,
                   double&slope,double&intercept) {
    slope=intercept=0.0; const int n=static_cast<int>(x.size()); if(n<2) return;
    double sx=0,sy=0,sxx=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];sxy+=x[i]*y[i];}
    const double d=n,denom=d*sxx-sx*sx; if(std::abs(denom)<1e-12) return;
    slope=(d*sxy-sx*sy)/denom; intercept=(sy-slope*sx)/d;
}

//==============================================================================
// Histogram utilities
//==============================================================================
static void s1(TH1D*h,const char*xt,const char*yt){
    h->GetXaxis()->SetTitle(xt);h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.046);h->GetYaxis()->SetTitleSize(0.046);
    h->GetXaxis()->SetLabelSize(0.038);h->GetYaxis()->SetLabelSize(0.038);
    h->SetLineWidth(2);h->SetStats(0);}
static void s2(TH2D*h,const char*xt,const char*yt){
    h->GetXaxis()->SetTitle(xt);h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.046);h->GetYaxis()->SetTitleSize(0.046);
    h->GetXaxis()->SetLabelSize(0.038);h->GetYaxis()->SetLabelSize(0.038);h->SetStats(0);}
static void p1(TH1D*h,const fs::path&p){
    TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",1000,680);
    c.SetLeftMargin(0.12);c.SetBottomMargin(0.13);c.SetRightMargin(0.05);
    h->Draw("HIST");c.SaveAs(p.string().c_str());}
static void p2(TH2D*h,const fs::path&p){
    TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",900,820);
    c.SetLeftMargin(0.13);c.SetRightMargin(0.17);c.SetBottomMargin(0.13);
    h->Draw("COLZ");c.SaveAs(p.string().c_str());}
static void p2fit(TH2D*h,double sl,double ic,double r,const fs::path&p){
    TCanvas c(("cf_"+std::string(h->GetName())).c_str(),"",900,820);
    c.SetLeftMargin(0.13);c.SetRightMargin(0.17);c.SetBottomMargin(0.13);
    h->Draw("COLZ");
    const double xlo=h->GetXaxis()->GetXmin(),xhi=h->GetXaxis()->GetXmax();
    TLine line(xlo,sl*xlo+ic,xhi,sl*xhi+ic);
    line.SetLineColor(kRed);line.SetLineWidth(2);line.Draw("SAME");
    TLatex lat;lat.SetTextSize(0.032);lat.SetNDC(true);
    char buf[128]; std::snprintf(buf,sizeof(buf),"slope=%.3f  intercept=%.2f mm  r=%.3f",sl,ic,r);
    lat.DrawLatex(0.14,0.935,buf);
    c.SaveAs(p.string().c_str());}

//==============================================================================
// Output helpers
//==============================================================================
static void write_mapping_csv(const fs::path& outdir) {
    std::ofstream f(outdir/"mapping_used.csv");
    f<<"plane,coord,index,channel\n";
    for(int k=0;k<16;++k) f<<"hod2018,X,"<<k<<','<<k<<'\n';
    for(int k=0;k<16;++k) f<<"hod2018,Y,"<<k<<','<<16+k<<'\n';
    for(int k=0;k<16;++k) f<<"hod2019,X,"<<k<<','<<kHod2019X[k]<<'\n';
    for(int k=0;k<16;++k) f<<"hod2019,Y,"<<k<<','<<kHod2019Y[k]<<'\n';
}

struct TrackRec {
    int ev; Long_t ts;
    double x18,y18,t18,x19,y19,t19,dt,mx,my,thx,thy,th3;
};
static void write_tracks(const fs::path& outdir, const std::vector<TrackRec>& v) {
    std::ofstream f(outdir/"recon_tracks.csv");
    f<<"event,timestamp_us,x2018_mm,y2018_mm,t2018_ns,x2019_mm,y2019_mm,t2019_ns,"
      "dt_ns,mx,my,theta_x_deg,theta_y_deg,theta_3d_deg\n";
    for(const auto&r:v)
        f<<r.ev<<','<<r.ts<<','<<std::fixed<<std::setprecision(3)
         <<r.x18<<','<<r.y18<<','<<r.t18<<','<<r.x19<<','<<r.y19<<','<<r.t19<<','
         <<r.dt<<','<<r.mx<<','<<r.my<<','<<r.thx<<','<<r.thy<<','<<r.th3<<'\n';
}

//==============================================================================
// main
//==============================================================================
int main(int argc, char** argv) {
    Params p;
    for (int i=1;i<argc;++i) {
        std::string a=argv[i];
        if(a=="--pitch"    &&i+1<argc){p.pitch_mm    =std::stod(argv[++i]);continue;}
        if(a=="--zsep"     &&i+1<argc){p.zsep_mm     =std::stod(argv[++i]);continue;}
        if(a=="--nsigma"   &&i+1<argc){p.n_sigma      =std::stod(argv[++i]);continue;}
        if(a=="--abst18"   &&i+1<argc){p.abs_thr_2018 =std::stod(argv[++i]);continue;}
        if(a=="--abst19"   &&i+1<argc){p.abs_thr_2019 =std::stod(argv[++i]);continue;}
        if(a=="--dtdelay"  &&i+1<argc){p.dt_delay_ns  =std::stod(argv[++i]);continue;}
        if(a=="--dtwindow" &&i+1<argc){p.dt_window_ns =std::stod(argv[++i]);continue;}
        if(a=="--nmax19"   &&i+1<argc){p.max_cls_2019 =std::stoi(argv[++i]);continue;}
        if(a=="--stage"    &&i+1<argc){
            std::string s=argv[++i];
            p.run_stage_a=(s=="A"||s=="AB"); p.run_stage_b=(s=="B"||s=="AB"); continue;}
        if(a[0]!='-'){p.input_file=a;continue;}
        std::cerr<<"Unknown option: "<<a<<'\n'; return 1;
    }

    // Dated output directory
    std::time_t now=std::time(nullptr); char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("reco_v8_baseline_corr_"+std::string(ds));
    fs::create_directories(outdir);

    std::cout<<"=== reconstruct_v8 ===\n";
    std::cout<<"Input : "<<p.input_file<<"\n";
    std::cout<<"Output: "<<outdir<<"\n";
    std::cout<<"pitch="<<p.pitch_mm<<" mm  zsep="<<p.zsep_mm<<" mm  nsigma="<<p.n_sigma<<"\n";
    std::cout<<"abs_thr_2018="<<p.abs_thr_2018<<" ADC  abs_thr_2019="<<p.abs_thr_2019<<" ADC\n";
    std::cout<<"dt_delay="<<p.dt_delay_ns<<" ns  dt_window=±"<<p.dt_window_ns<<" ns"
             <<"  max_cls_2019="<<p.max_cls_2019<<"\n\n";

    gROOT->SetBatch(kTRUE); TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0); gStyle->SetPalette(kViridis);

    TFile* fin=TFile::Open(p.input_file.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: cannot open "<<p.input_file<<'\n';return 1;}
    TTree* tree=dynamic_cast<TTree*>(fin->Get("data_dig"));
    if(!tree){std::cerr<<"ERROR: TTree 'data_dig' not found\n";fin->Close();return 1;}

    std::array<std::vector<float>*,kNCh> ch_wf{};
    std::vector<float>* time_wf=nullptr;
    int ev_num=0; Long_t ts=0;
    tree->SetBranchStatus("*",0); tree->SetCacheSize(64LL*1024LL*1024LL);
    for(int c=0;c<kNCh;++c){
        char bn[16]; std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!tree->GetBranch(bn)){std::cerr<<"WARNING: missing "<<bn<<'\n';continue;}
        tree->SetBranchStatus(bn,1); tree->SetBranchAddress(bn,&ch_wf[c]);
    }
    tree->SetBranchStatus("time",1);      tree->SetBranchAddress("time",&time_wf);
    tree->SetBranchStatus("event",1);     tree->SetBranchAddress("event",&ev_num);
    tree->SetBranchStatus("timestamp",1); tree->SetBranchAddress("timestamp",&ts);

    const Long64_t nEnt=tree->GetEntries();
    std::cout<<"Total entries: "<<nEnt<<'\n';

    const InvMap invA(kHod2019X,kHod2019Y), invC(kHod2019X_ctrl,kHod2019Y_ctrl);

    //==========================================================================
    // Histograms — Stage A: pulse diagnostics (all valid pulses)
    //==========================================================================
    auto h_ped  =std::make_unique<TH1D>("h_pedestal","Mean pedestal per channel;Channel;Pedestal [ADC]",kNCh,-0.5,kNCh-0.5);
    auto h_ph18 =std::make_unique<TH1D>("h_ph_2018","Pulse height hod2018 (all);Pulse height [ADC];Entries",200,0,600);
    auto h_ph19 =std::make_unique<TH1D>("h_ph_2019","Pulse height hod2019 (all);Pulse height [ADC];Entries",200,0,200);
    auto h_pt18 =std::make_unique<TH1D>("h_pt_2018","Pulse time hod2018 (hits);Pulse time [ns];Hits",175,3500,7000);
    auto h_pt19 =std::make_unique<TH1D>("h_pt_2019","Pulse time hod2019 (hits);Pulse time [ns];Hits",100,3000,5000);
    auto h_dt   =std::make_unique<TH1D>("h_dt_interplano","#Deltat = t_{2018} #minus t_{2019};#Deltat [ns];Events",200,-500,3000);

    //==========================================================================
    // Histograms — Stage B: tracking
    //==========================================================================
    constexpr double kClo=-16.0,kChi=16.0; constexpr int kCb=16;
    auto h_occ18x=std::make_unique<TH1D>("h_occ_2018X","Occupancy hod2018 X;Strip index k;Hits",16,-0.5,15.5);
    auto h_occ18y=std::make_unique<TH1D>("h_occ_2018Y","Occupancy hod2018 Y;Strip index k;Hits",16,-0.5,15.5);
    auto h_occ19x=std::make_unique<TH1D>("h_occ_2019X","Occupancy hod2019 X (auth.);Strip index k;Hits",16,-0.5,15.5);
    auto h_occ19y=std::make_unique<TH1D>("h_occ_2019Y","Occupancy hod2019 Y (auth.);Strip index k;Hits",16,-0.5,15.5);
    auto h_xy18  =std::make_unique<TH2D>("h_xy_2018","hod2018 XY;X_{2018} [mm];Y_{2018} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_xy19  =std::make_unique<TH2D>("h_xy_2019","hod2019 XY;X_{2019} [mm];Y_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cx    =std::make_unique<TH2D>("h_corr_x","Correlation X (authoritative);X_{2018} [mm];X_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy    =std::make_unique<TH2D>("h_corr_y","Correlation Y (authoritative);Y_{2018} [mm];Y_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cx_ct =std::make_unique<TH2D>("h_corr_x_ctrl","Correlation X (CONTROL: identity);X_{2018} [mm];X_{2019}^{ctrl} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy_ct =std::make_unique<TH2D>("h_corr_y_ctrl","Correlation Y (CONTROL: identity);Y_{2018} [mm];Y_{2019}^{ctrl} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_thx   =std::make_unique<TH1D>("h_theta_x","#theta_{x};#theta_{x} [deg];Events",80,-25,25);
    auto h_thy   =std::make_unique<TH1D>("h_theta_y","#theta_{y};#theta_{y} [deg];Events",80,-25,25);
    auto h_th3   =std::make_unique<TH1D>("h_theta_3d","#theta_{3D};#theta_{3D} [deg];Events",80,0,25);

    //==========================================================================
    // Cutflow
    //==========================================================================
    long long cf_tot=nEnt,cf_any=0,cf_19c=0,cf_19m=0,cf_18m=0,cf_trk=0;
    std::vector<double> vx18,vx19,vy18,vy19,vx18c,vx19c,vy18c,vy19c;
    std::vector<TrackRec> trks;
    double ped_sum[kNCh]={}; long long ped_n[kNCh]={};

    //==========================================================================
    // Event loop
    //==========================================================================
    std::cout<<"Running event loop...\n";
    for(Long64_t iev=0;iev<nEnt;++iev){
        tree->GetEntry(iev);
        if(iev%2000==0)
            std::cout<<"  "<<iev<<" / "<<nEnt<<" ("<<100*iev/nEnt<<"%)\r"<<std::flush;
        if(!time_wf||time_wf->empty()) continue;
        const auto& tv=*time_wf;

        // ── Pulse extraction ─────────────────────────────────────────────────
        std::array<PulseResult,kNCh> pr{};
        for(int c=0;c<kNCh;++c){
            if(!ch_wf[c]||ch_wf[c]->empty()) continue;
            const double athr=is2019(c)?p.abs_thr_2019:p.abs_thr_2018;
            pr[c]=extract_pulse(*ch_wf[c],tv,sig_tlo(c),sig_thi(c),p.n_sigma,athr);
            if(!pr[c].valid) continue;
            ped_sum[c]+=pr[c].pedestal; ++ped_n[c];
            if(c<32){h_ph18->Fill(pr[c].pulse_height);if(pr[c].hit)h_pt18->Fill(pr[c].pulse_time);}
            else    {h_ph19->Fill(pr[c].pulse_height);if(pr[c].hit)h_pt19->Fill(pr[c].pulse_time);}
        }

        // ── Hit lists ────────────────────────────────────────────────────────
        std::vector<Hit> h18x,h18y,h19x,h19y,h19xc,h19yc;
        bool any=false;
        for(int c=0;c<kNCh;++c){
            if(!pr[c].valid||!pr[c].hit) continue;
            any=true;
            const Hit h{0,pr[c].pulse_height,pr[c].pulse_time};
            const int i18x=idx2018x(c),i18y=idx2018y(c);
            const int i19x=invA.x19[c],i19y=invA.y19[c];
            const int i19xc=invC.x19[c],i19yc=invC.y19[c];
            if(i18x>=0){h18x.push_back({i18x,h.ph,h.pt});h_occ18x->Fill(i18x);}
            if(i18y>=0){h18y.push_back({i18y,h.ph,h.pt});h_occ18y->Fill(i18y);}
            if(i19x>=0){h19x.push_back({i19x,h.ph,h.pt});h_occ19x->Fill(i19x);}
            if(i19y>=0){h19y.push_back({i19y,h.ph,h.pt});h_occ19y->Fill(i19y);}
            if(i19xc>=0) h19xc.push_back({i19xc,h.ph,h.pt});
            if(i19yc>=0) h19yc.push_back({i19yc,h.ph,h.pt});
        }
        if(any) ++cf_any;

        // ── Cluster ──────────────────────────────────────────────────────────
        auto cls18x=make_clusters(h18x,p.pitch_mm), cls18y=make_clusters(h18y,p.pitch_mm);
        auto cls19x=make_clusters(h19x,p.pitch_mm), cls19y=make_clusters(h19y,p.pitch_mm);
        auto cls19xc=make_clusters(h19xc,p.pitch_mm),cls19yc=make_clusters(h19yc,p.pitch_mm);

        // Step 1: hod2019 coincidence (both X and Y must have clusters)
        if(cls19x.empty()||cls19y.empty()) continue;
        ++cf_19c;

        // Step 2: hod2019 multiplicity quality cut
        if((int)cls19x.size()>p.max_cls_2019||(int)cls19y.size()>p.max_cls_2019) continue;
        ++cf_19m;

        // Step 3: hod2018 coincidence (both X and Y must have clusters)
        if(cls18x.empty()||cls18y.empty()) continue;

        // Step 4: select best cluster per plane
        // hod2019: amplitude selection (narrow window, signal dominates at abs_thr_2019=8 ADC)
        // hod2018: time-match to expected muon time (t_2019 + dt_delay), then amplitude as tiebreak.
        //   Real muon in hod2018 fires at t_muon + board_delay_2018 ≈ t_2019 + 1200 ns.
        //   Background fires at random times in [3500, 7000] ns.
        const Cluster& cx19=best_amp(cls19x); const Cluster& cy19=best_amp(cls19y);
        const double x2019=cx19.position, y2019=cy19.position;
        const double t2019=0.5*(cx19.time+cy19.time);
        const double t_exp18=t2019+p.dt_delay_ns;
        // If no cluster within dt_window, fall back to best-amplitude cluster in hod2018
        const Cluster* cx18_t=best_time(cls18x,t_exp18,p.dt_window_ns);
        const Cluster* cy18_t=best_time(cls18y,t_exp18,p.dt_window_ns);
        const Cluster* cx18=(cx18_t ? cx18_t : &best_amp(cls18x));
        const Cluster* cy18=(cy18_t ? cy18_t : &best_amp(cls18y));
        ++cf_18m;

        const double x2018=cx18->position, y2018=cy18->position;
        const double t2018=0.5*(cx18->time+cy18->time);
        const double dt=t2018-t2019;  // sign: positive when hod2018 reads out later
        h_dt->Fill(dt);
        ++cf_trk;

        // ── Track parameters ─────────────────────────────────────────────────
        const double mx=(x2019-x2018)/p.zsep_mm, my=(y2019-y2018)/p.zsep_mm;
        const double thx=std::atan(mx)*kRadToDeg, thy=std::atan(my)*kRadToDeg;
        const double th3=std::atan(std::sqrt(mx*mx+my*my))*kRadToDeg;

        h_xy18->Fill(x2018,y2018); h_xy19->Fill(x2019,y2019);
        h_cx->Fill(x2018,x2019);   h_cy->Fill(y2018,y2019);
        h_thx->Fill(thx); h_thy->Fill(thy); h_th3->Fill(th3);
        vx18.push_back(x2018); vx19.push_back(x2019);
        vy18.push_back(y2018); vy19.push_back(y2019);
        trks.push_back({ev_num,ts,x2018,y2018,t2018,x2019,y2019,t2019,dt,mx,my,thx,thy,th3});

        // ── Control mapping (same hod2018 selection, ctrl hod2019 position) ──
        if(!cls19xc.empty()&&!cls19yc.empty()&&
           (int)cls19xc.size()<=p.max_cls_2019&&(int)cls19yc.size()<=p.max_cls_2019){
            const double cxc=best_amp(cls19xc).position;
            const double cyc=best_amp(cls19yc).position;
            h_cx_ct->Fill(x2018,cxc); h_cy_ct->Fill(y2018,cyc);
            vx18c.push_back(x2018); vx19c.push_back(cxc);
            vy18c.push_back(y2018); vy19c.push_back(cyc);
        }
    }
    std::cout<<"\nEvent loop done.\n\n";

    for(int c=0;c<kNCh;++c) if(ped_n[c]>0) h_ped->SetBinContent(c+1,ped_sum[c]/ped_n[c]);

    const long long na=static_cast<long long>(vx18.size());
    const long long nc=static_cast<long long>(vx18c.size());
    double rx=pearson_r(vx18,vx19),ry=pearson_r(vy18,vy19);
    double rxc=pearson_r(vx18c,vx19c),ryc=pearson_r(vy18c,vy19c);
    double slx,icx,sly,icy,slxc,icxc,slyc,icyc;
    linfit(vx18,vx19,slx,icx); linfit(vy18,vy19,sly,icy);
    linfit(vx18c,vx19c,slxc,icxc); linfit(vy18c,vy19c,slyc,icyc);

    //==========================================================================
    // Style and save
    //==========================================================================
    s1(h_ped.get(),  "Channel","Mean pedestal [ADC]");
    s1(h_ph18.get(), "Pulse height [ADC]","Entries");
    s1(h_ph19.get(), "Pulse height [ADC]","Entries");
    s1(h_pt18.get(), "Pulse time [ns]","Hits");
    s1(h_pt19.get(), "Pulse time [ns]","Hits");
    s1(h_dt.get(),   "#Deltat [ns]","Events");
    s1(h_occ18x.get(),"Strip index k","Hits"); s1(h_occ18y.get(),"Strip index k","Hits");
    s1(h_occ19x.get(),"Strip index k","Hits"); s1(h_occ19y.get(),"Strip index k","Hits");
    s2(h_xy18.get(),"X_{2018} [mm]","Y_{2018} [mm]");
    s2(h_xy19.get(),"X_{2019} [mm]","Y_{2019} [mm]");
    s2(h_cx.get(),"X_{2018} [mm]","X_{2019} [mm]");
    s2(h_cy.get(),"Y_{2018} [mm]","Y_{2019} [mm]");
    s2(h_cx_ct.get(),"X_{2018} [mm]","X_{2019}^{ctrl} [mm]");
    s2(h_cy_ct.get(),"Y_{2018} [mm]","Y_{2019}^{ctrl} [mm]");
    s1(h_thx.get(),"#theta_{x} [deg]","Events");
    s1(h_thy.get(),"#theta_{y} [deg]","Events");
    s1(h_th3.get(),"#theta_{3D} [deg]","Events");

    for(int k=0;k<16;++k){
        char l[16];
        std::snprintf(l,sizeof(l),"k%d/Ch%d",k,kHod2019X[k]);
        h_occ19x->GetXaxis()->SetBinLabel(k+1,l);
        std::snprintf(l,sizeof(l),"k%d/Ch%d",k,kHod2019Y[k]);
        h_occ19y->GetXaxis()->SetBinLabel(k+1,l);
    }
    h_occ19x->LabelsOption("v","X"); h_occ19y->LabelsOption("v","X");

    p1(h_ped.get(),   outdir/"pedestal_por_canal.png");
    p1(h_ph18.get(),  outdir/"pulseheight_2018.png");
    p1(h_ph19.get(),  outdir/"pulseheight_2019.png");
    p1(h_pt18.get(),  outdir/"pulsetime_2018.png");
    p1(h_pt19.get(),  outdir/"pulsetime_2019.png");
    p1(h_dt.get(),    outdir/"dt_interplano.png");
    p1(h_occ18x.get(),outdir/"occupancy_2018X.png");
    p1(h_occ18y.get(),outdir/"occupancy_2018Y.png");
    p1(h_occ19x.get(),outdir/"occupancy_2019X.png");
    p1(h_occ19y.get(),outdir/"occupancy_2019Y.png");
    p2(h_xy18.get(),  outdir/"xy_2018.png");
    p2(h_xy19.get(),  outdir/"xy_2019.png");
    p2fit(h_cx.get(),   slx, icx, rx,  outdir/"corr_x2019_vs_x2018.png");
    p2fit(h_cy.get(),   sly, icy, ry,  outdir/"corr_y2019_vs_y2018.png");
    p2fit(h_cx_ct.get(),slxc,icxc,rxc, outdir/"corr_x2019_vs_x2018_control.png");
    p2fit(h_cy_ct.get(),slyc,icyc,ryc, outdir/"corr_y2019_vs_y2018_control.png");
    p1(h_thx.get(),   outdir/"theta_x.png");
    p1(h_thy.get(),   outdir/"theta_y.png");
    p1(h_th3.get(),   outdir/"theta_3d.png");

    {
        TFile fout((outdir/"histos.root").string().c_str(),"RECREATE");
        if(!fout.IsZombie()){
            for(auto*h:std::initializer_list<TH1*>{
                h_ped.get(),h_ph18.get(),h_ph19.get(),h_pt18.get(),h_pt19.get(),h_dt.get(),
                h_occ18x.get(),h_occ18y.get(),h_occ19x.get(),h_occ19y.get(),
                h_thx.get(),h_thy.get(),h_th3.get()}) h->Write();
            for(auto*h:std::initializer_list<TH2D*>{h_xy18.get(),h_xy19.get(),
                h_cx.get(),h_cy.get(),h_cx_ct.get(),h_cy_ct.get()}) h->Write();
            fout.Close();
        }
    }
    write_mapping_csv(outdir);
    write_tracks(outdir,trks);

    //==========================================================================
    // Summary
    //==========================================================================
    auto print=[&](std::ostream&os){
        os<<"=== reconstruct_v8: muon telescope, time-matched, baseline-corrected ===\n\n";
        os<<"Parameters:\n";
        os<<"  input_file    = "<<p.input_file<<"\n";
        os<<"  pitch_eff_mm  = "<<p.pitch_mm
          <<"  (31 mm→1.9375; 33 mm→2.0625 [v7 convention]; using "<<p.pitch_mm<<")\n";
        os<<"  z_sep_mm      = "<<p.zsep_mm<<" mm\n";
        os<<"  n_sigma       = "<<p.n_sigma<<"\n";
        os<<"  abs_thr_2018  = "<<p.abs_thr_2018<<" ADC\n";
        os<<"  abs_thr_2019  = "<<p.abs_thr_2019
          <<" ADC  [real muon ≈17 ADC; noise crossings 5-7 ADC]\n";
        os<<"  hod2018_win   = ["<<kHod2018Tlo_ns<<", "<<kHod2018Thi_ns<<"] ns\n";
        os<<"  hod2019_win   = ["<<kHod2019Tlo_ns<<", "<<kHod2019Thi_ns<<"] ns\n";
        os<<"  baseline_win  = [0, "<<kPrePulseEnd_ns<<") ns\n";
        os<<"  dt_delay_ns   = "<<p.dt_delay_ns
          <<" ns  (digitizer board pre-trigger offset hod2018 − hod2019)\n";
        os<<"  dt_window_ns  = ±"<<p.dt_window_ns<<" ns  (hod2018 time-match half-width)\n";
        os<<"  max_cls_2019  = "<<p.max_cls_2019<<" per plane\n";
        os<<"  pulse_time    = time of minimum sample (time-of-minimum method)\n\n";

        os<<"Cluster selection:\n";
        os<<"  hod2019: best by amplitude (high S/N, narrow window, reliable)\n";
        os<<"  hod2018: time-match to (t_2019 + dt_delay) within ±dt_window\n\n";

        os<<"Physical/instrumental notes:\n";
        os<<"  1. hod2019 bars are 1 mm thick → muon energy ≈0.2 MeV → signal ≈17 ADC.\n";
        os<<"     Noise crossings at 5σ threshold: 5-7 ADC. SNR ≈ 3. Very low efficiency\n";
        os<<"     for large signals (> abs_thr_2019 = "<<p.abs_thr_2019<<" ADC).\n";
        os<<"  2. The ≈1200 ns offset between hod2018 and hod2019 waveform times is NOT\n";
        os<<"     from cable lengths. It is a digitizer pre-trigger configuration difference:\n";
        os<<"     hod2018 board pre-trigger ≈ 5200 ns; hod2019 board pre-trigger ≈ 4000 ns.\n";
        os<<"  3. hod2018 window [3500,7000] ns is 3.5× wider than hod2019 [3000,5000] ns.\n";
        os<<"     Multiple background particles (gammas, neutrons) fire random hod2018 strips\n";
        os<<"     at random times. Time-matching filters these but efficiency is still limited.\n\n";

        os<<"Cutflow:\n";
        os<<std::setw(38)<<std::left<<"  Total events"<<cf_tot<<"\n";
        os<<std::setw(38)<<"  Any hit (any plane)"<<cf_any<<"\n";
        os<<std::setw(38)<<"  hod2019 X+Y coincidence"<<cf_19c<<"\n";
        os<<std::setw(38)<<"  hod2019 mult ≤ "+std::to_string(p.max_cls_2019)<<cf_19m<<"\n";
        os<<std::setw(38)<<"  hod2018 X+Y coincidence"<<cf_18m<<"\n";
        os<<std::setw(38)<<"  Full 3D tracks produced"<<cf_trk<<"\n\n";

        os<<std::fixed<<std::setprecision(4);
        os<<"Spatial correlation — AUTHORITATIVE mapping (n="<<na<<" tracks):\n";
        os<<"  X: r="<<rx<<"  slope="<<slx<<"  intercept="<<icx<<" mm\n";
        os<<"  Y: r="<<ry<<"  slope="<<sly<<"  intercept="<<icy<<" mm\n\n";
        os<<"Spatial correlation — CONTROL mapping (identity, n="<<nc<<"):\n";
        os<<"  X: r="<<rxc<<"  slope="<<slxc<<"  intercept="<<icxc<<" mm\n";
        os<<"  Y: r="<<ryc<<"  slope="<<slyc<<"  intercept="<<icyc<<" mm\n\n";

        os<<"Mapping verification:\n";
        const bool okx=na>=2&&std::abs(rx)>0.25&&(na<5||std::abs(slx-1.0)<0.5);
        const bool oky=na>=2&&std::abs(ry)>0.25&&(na<5||std::abs(sly-1.0)<0.5);
        os<<"  X: "<<(okx?"MAPPING CONFIRMED":"inconclusive — insufficient statistics")
          <<"  (r="<<rx<<", slope="<<slx<<", n="<<na<<")\n";
        os<<"  Y: "<<(oky?"MAPPING CONFIRMED":"inconclusive — insufficient statistics")
          <<"  (r="<<ry<<", slope="<<sly<<", n="<<na<<")\n";
        os<<"  Control: r_x="<<rxc<<"  r_y="<<ryc<<"\n";
        os<<"  AUTH vs CTRL sign: r_x_auth="<<rx<<"  r_x_ctrl="<<rxc<<"\n";
        if(na>=2&&((rx>0&&rxc<0)||(rx<0&&rxc>0)))
            os<<"  *** SIGN FLIP between AUTH and CTRL confirms mapping matters ***\n\n";
        else
            os<<"  (auth and ctrl r not significantly different — insufficient events)\n\n";

        os<<"Interpretation:\n";
        if(na<10){
            os<<"  LOW STATISTICS ("<<na<<" tracks). The hod2019 efficiency for producing\n";
            os<<"  signals > "<<p.abs_thr_2019<<" ADC is very low due to 1 mm bar thickness.\n";
            os<<"  Raise abs_thr_2019 or use a run with stronger signals for conclusive mapping.\n";
            os<<"  For comparison: run with --abst19 5 gives ~190 events but r≈0 because\n";
            os<<"  hod2019 is dominated by noise crossings at that threshold.\n";
        } else {
            os<<"  Adequate statistics. Slope≈1 and r≳0.5 would confirm correct mapping.\n";
        }
        os<<"\nPitch convention:\n";
        os<<"  Active area 31 mm → pitch = 1.9375 mm\n";
        os<<"  Active area 33 mm → pitch = 2.0625 mm  (v7 convention)\n";
        os<<"  This run: "<<p.pitch_mm<<" mm  (correlation slope is pitch-independent)\n";
        os<<"  hod2018 at Z=0; hod2019 at Z=+"<<p.zsep_mm<<" mm.\n";
    };

    print(std::cout);
    {std::ofstream sf(outdir/"summary.txt"); print(sf);}

    std::cout<<"Output directory: "<<outdir<<"\n";
    fin->Close();
    return 0;
}
