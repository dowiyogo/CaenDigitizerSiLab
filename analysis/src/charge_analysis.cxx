// charge_analysis.cxx
// Análisis de carga integrada Q = Σ(pedestal − V(t)) × Δt_sample [mV·ns]
// para hod2018 y hod2019 en ambos archivos v7.root y v8.root.
//
// Ventaja vs pulse_height: un spike de 1 muestra tiene ph alto pero Q pequeño
// → la integración rechaza mejor los crossings de ruido estrecho (DCR, EMI).
//
// Integración adaptativa: para cada evento/canal, integra en
//   [t_min − N_pre·Δt, t_min + N_post·Δt] donde t_min es el tiempo del mínimo.
// Δt_sample = 16.0 ns/muestra (VX1740, 62.5 MS/s). N_pre=3, N_post=8 → 12 muestras.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o charge_analysis charge_analysis.cxx

#include <TBranch.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <THStack.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "hodo_consts.h"

namespace fs = std::filesystem;

//==============================================================================
// Constantes del detector
//==============================================================================
constexpr double kPre_ns  = 2500.0;   // fin región baseline [ns]
constexpr int    kNCh     = 64;

static const double kTlo[4] = {3500, 3500, 3000, 3000};  // ventana señal [ns]
static const double kThi[4] = {7000, 7000, 5000, 5000};

// Planos: 0=hod2018X, 1=hod2018Y, 2=hod2019X, 3=hod2019Y
constexpr int P18X=0,P18Y=1,P19X=2,P19Y=3;
static const char* kPname[4]={"hod2018X","hod2018Y","hod2019X","hod2019Y"};

// Mapeo canales → plano e índice k
static int ch_plane(int c){
    if(c>=0&&c<=15)return P18X; if(c>=16&&c<=31)return P18Y;
    if(c>=32&&c<=47)return P19Y; if(c>=48&&c<=63)return P19X; return -1;}
static int ch_k(int c){
    if(c>=0&&c<=15)return c; if(c>=16&&c<=31)return c-16;
    static const std::array<int,16> h19y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};
    static const std::array<int,16> h19x={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
    for(int k=0;k<16;k++)if(h19y[k]==c)return k;
    for(int k=0;k<16;k++)if(h19x[k]==c)return k;
    return -1;}

// Umbral alto para Step 0 (selección de señales reales para timing)
constexpr double kThrHigh18 = 15.0;   // mV — hod2018
constexpr double kThrHigh19_v7 = 17.0; // mV — hod2019 v7
constexpr double kThrHigh19_v8 = 35.0; // mV, pulse height from pedestal; = 72 ADC counts

//==============================================================================
// Resultado de extracción de pulso + carga integrada
//==============================================================================
struct CQ {
    bool   valid  = false;
    double ped    = 0;       // pedestal [mV]
    double noise  = 0;       // RMS baseline [mV]
    double ph     = 0;       // pulse height [mV]
    double pt     = -1;      // tiempo del mínimo [ns]
    double Q      = 0;       // carga integrada [mV·ns]
    double sigma_Q = 0;      // error en Q por ruido: noise × Δt × √N_samples
    int    nw     = 0;       // muestras en la ventana de integración
    // r = Q / (ph × N × Δt) — indicador de forma del pulso
    double r      = 0;
};

static CQ extract_cq(const std::vector<float>& wf, const std::vector<float>& t,
                      double tlo, double thi, int kNpost, bool use_loc) {
    CQ r;
    const int n = (int)wf.size();
    if (n < 4 || (int)t.size() != n) return r;

    // Baseline
    double s=0, s2=0; int nb=0;
    for(int i=0;i<n;i++) if(t[i]<kPre_ns){s+=wf[i];s2+=wf[i]*wf[i];nb++;}
    if(nb<2) return r;
    r.valid=true; r.ped=s/nb;
    const double var=s2/nb-r.ped*r.ped;
    r.noise=(var>0)?std::sqrt(var):0.0;

    // Mínimo en ventana de señal
    double vmin=1e30; int imin=-1;
    for(int i=0;i<n;i++) if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0) return r;
    r.ph = r.ped - vmin;
    r.pt = t[imin];

    // Baseline for Q integration
    double ped_Q = r.ped;
    if (use_loc && imin >= 12 && imin - 4 < n) {
        double s_loc = 0;
        for (int j = -12; j <= -4; ++j) s_loc += wf[imin + j];
        ped_Q = s_loc / 9.0;
    }

    // Ventana de integración adaptativa: [t_min - kNpre×Δt, t_min + kNpost×Δt]
    const double t_int_lo = r.pt - kNpre * kDt_ns;
    const double t_int_hi = r.pt + kNpost * kDt_ns;
    double Q = 0; int nw = 0;
    for(int i=0;i<n;i++) if(t[i]>=t_int_lo&&t[i]<=t_int_hi){
        Q += (ped_Q - wf[i]) * kDt_ns; nw++;}
    r.Q = Q;
    r.nw = nw;
    r.sigma_Q = r.noise * kDt_ns * std::sqrt(std::max(nw,1));
    // Ratio de forma (0=ruido estrecho, 1=rectangular, 0.5=triangular)
    if(r.ph > 0.5 && nw > 0)
        r.r = r.Q / (r.ph * nw * kDt_ns);
    return r;
}

//==============================================================================
// Ajuste gaussiano y Landau
//==============================================================================
struct GausFit{double mean=0,sigma=0,chi2ndf=0;bool ok=false;};
static GausFit fit_gaus(TH1D*h,double xlo,double xhi){
    GausFit g;
    if(h->Integral(h->FindBin(xlo),h->FindBin(xhi))<5)return g;
    TF1 f("_fg","gaus",xlo,xhi);
    f.SetParameters(h->GetMaximum(),h->GetMean(),h->GetRMS());
    if(h->Fit(&f,"RQNS")==0){
        g.ok=true;g.mean=f.GetParameter(1);g.sigma=std::abs(f.GetParameter(2));
        g.chi2ndf=f.GetChisquare()/(f.GetNDF()>0?f.GetNDF():1);}
    return g;}

struct LandauFit{
    double mpv=0,sigma=0;
    int status=-1;
    double chi2ndf=0;
    bool at_limit=false;
    bool ok=false;
};

static LandauFit fit_landau_Q(TH1D*h,double xlo,double xhi){
    LandauFit l;
    if(h->Integral(h->FindBin(xlo),h->FindBin(xhi))<10)return l;
    const double mpv0=h->GetBinCenter(h->GetMaximumBin());
    TF1 f("_fl","TMath::Landau(x,[0],[1],0)",xlo,xhi);
    f.SetParameters(std::max(mpv0,xlo+1),std::max((xhi-xlo)*0.1,10.0));
    f.SetParLimits(0,xlo,xhi);f.SetParLimits(1,1.0,xhi-xlo);
    l.status = h->Fit(&f,"RQNS");
    if(l.status==0){
        l.ok=true;
        l.mpv=f.GetParameter(0);
        l.sigma=f.GetParameter(1);
        l.chi2ndf=f.GetNDF()>0?f.GetChisquare()/f.GetNDF():0.0;
        double lim0_lo,lim0_hi,lim1_lo,lim1_hi;
        f.GetParLimits(0,lim0_lo,lim0_hi);
        f.GetParLimits(1,lim1_lo,lim1_hi);
        auto near_lim = [](double val, double lo, double hi) {
            double span = hi - lo;
            return span > 0 && ((std::abs(val - lo) <= 0.001 * span) || (std::abs(val - hi) <= 0.001 * span));
        };
        if(near_lim(l.mpv,lim0_lo,lim0_hi) || near_lim(l.sigma,lim1_lo,lim1_hi))
            l.at_limit = true;
    }
    return l;
}

//==============================================================================
// Helpers de gráficas
//==============================================================================
static void sl(TH1D*h,const char*x,const char*y){
    h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    h->GetXaxis()->SetTitleSize(0.05);h->GetYaxis()->SetTitleSize(0.05);h->SetStats(0);}
static void sl2(TH2D*h,const char*x,const char*y){
    h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    h->GetXaxis()->SetTitleSize(0.05);h->GetYaxis()->SetTitleSize(0.05);h->SetStats(0);}

//==============================================================================
// Datos acumulados por archivo y por plano
//==============================================================================
struct PlaneQ {
    // Histogramas de carga integrada
    TH1D* h_Q    = nullptr;   // max(Q_c) por evento
    TH1D* h_Qall = nullptr;   // todos los Q_c > -1000 (incluye pedestal completo)
    // Histogramas de pulse_time al umbral alto (Step 0)
    TH1D* h_pt_hi = nullptr;
    GausFit   gaus_ped;
    LandauFit lan_Q;
    double sigma_Q_ped = 0;   // del fit gaussiano al pedestal
    double umbral_Q = 0;      // valle en espacio Q
    bool has_valley = false;  // valle interior encontrado
    double eff_nat_Q = 0;     // eficiencia al umbral Q natural [%]
    // Eficiencia vs umbral
    std::vector<double> Q_per_event;   // max(Q) por evento (para eficiencia)
    // Estimador de carga robusto (pulsos sobre umbral en ventana de muones)
    std::vector<double> signal_Q;
    long long n_sig_Q = 0;
    std::vector<double> nopulse_Q;    // Q de canales sin pulso (ph < 5 mV en [2500,8000] ns)
    double sigma_Q_noise = 0;         // RMS de nopulse_Q (F1)
    int blo = 0, bhi = 0, min_b = 0;  // bins de búsqueda de valle (F2)
    double med_Q = 0;
    double q10_Q = 0;
    double q90_Q = 0;
};

struct FileQ {
    std::string label;
    Long64_t n_events = 0;
    PlaneQ planes[4];
    // Per-channel (solo hod2019 X e Y, y hod2018 X e Y para v8)
    TH1D* h_Q_ch[kNCh] = {};
    // Step 3: Q vs ph para canal más activo (Ch50 hod2019X)
    TH2D* h_QvsPH = nullptr;
    TH1D* h_ratio_r = nullptr;  // ratio r = Q/(ph×N×Δt)
    // Ruido por canal
    double noise[kNCh] = {}, sigma_Q_ch[kNCh] = {};
    long long n_ped[kNCh] = {};
};

//==============================================================================
// Procesa un archivo
//==============================================================================
static void process_file_Q(FileQ& fq, const std::string& fname,
                             double thr19_hi, bool make_per_channel) {
    TFile* fin = TFile::Open(fname.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: "<<fname<<'\n';return;}
    TTree* t = static_cast<TTree*>(fin->Get("data_dig"));
    if(!t){std::cerr<<"ERROR: no data_dig\n";fin->Close();return;}

    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw = nullptr;
    t->SetBranchStatus("*",0); t->SetCacheSize(256LL<<20);
    for(int c=0;c<kNCh;c++){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!t->GetBranch(bn))continue;
        t->SetBranchStatus(bn,1);t->SetBranchAddress(bn,&cw[c]);}
    t->SetBranchStatus("time",1);t->SetBranchAddress("time",&tw);

    const Long64_t nEnt = t->GetEntries();
    fq.n_events = nEnt;
    std::cout<<"  ["<<fq.label<<"] "<<nEnt<<" eventos...\n";

    // Rango de Q por plano
    const double Qlo18=-2000, Qhi18=25000, Qlo19=-500, Qhi19=5000;

    for(int pl=0;pl<4;pl++){
        auto& p = fq.planes[pl];
        const bool is19=(pl>=P19X);
        const double Qlo=is19?Qlo19:Qlo18, Qhi=is19?Qhi19:Qhi18;
        const std::string nm=fq.label+"_"+kPname[pl];
        p.h_Q    =new TH1D((nm+"_Q").c_str(),   (std::string(kPname[pl])+" "+fq.label+";Q_{max} [mV#upointns];Entries").c_str(),200,Qlo,Qhi);
        p.h_Qall =new TH1D((nm+"_Qall").c_str(), "",200,Qlo,Qhi);
        p.h_pt_hi=new TH1D((nm+"_pth").c_str(),  (std::string(kPname[pl])+" t_{pulso} (umbral alto);t [ns];Hits").c_str(),100,2500,8000);
    }
    if(make_per_channel)
        for(int c=0;c<kNCh;c++){
            const int pl=ch_plane(c);if(pl<0)continue;
            const bool is19=(pl>=P19X);
            const double Qlo=is19?Qlo19:Qlo18,Qhi=is19?Qhi19:Qhi18;
            char bn[32];std::snprintf(bn,sizeof(bn),"%s_ch%d",fq.label.c_str(),c);
            fq.h_Q_ch[c]=new TH1D(bn,bn,200,Qlo,Qhi);}
    // Q vs ph para Ch50 (hod2019X k=10)
    fq.h_QvsPH =new TH2D((fq.label+"_Qvph").c_str(),
        ("Q vs PH Ch50 hod2019X "+fq.label+";Pulse height [mV];Q [mV#upointns];Eventos").c_str(),
        100,0,100,100,Qlo19,Qhi19);
    fq.h_ratio_r=new TH1D((fq.label+"_r").c_str(),
        ("Ratio r = Q/(ph#timesN#timesdt) "+fq.label+";r;Hits").c_str(),
        100,-0.1,1.5);

    bool checked_step = false;
    for(Long64_t iev=0;iev<nEnt;iev++){
        t->GetEntry(iev);
        if(iev%30000==0)std::cout<<"    "<<iev<<"/"<<nEnt<<"\r"<<std::flush;
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;
        if (!checked_step && tv.size() > 1) {
            checked_step = true;
            if (std::fabs((tv[1] - tv[0]) - kDt_ns) > 1e-3) {
                std::fprintf(stderr, "Sampling step %.4f ns != kDt_ns %.4f ns; aborting\n", tv[1]-tv[0], kDt_ns);
                std::exit(1);
            }
        }

        // Máximo Q por plano para este evento
        double maxQ[4]={-1e9,-1e9,-1e9,-1e9};

        for(int c=0;c<kNCh;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const int pl=ch_plane(c);if(pl<0)continue;
            const bool is19=(c>=32);
            const int npost=is19?kNpost19:kNpost18;
            const bool use_loc=is19?kQLocalPed19:kQLocalPed18;
            const CQ cq=extract_cq(*cw[c],tv,kTlo[pl],kThi[pl],npost,use_loc);
            if(!cq.valid)continue;

            // Estimador de carga robusto: pulsos con ph > umbral de estación y pt en ventana de muones
            const double ph_thr = is19 ? kPh19_mV : kPh18_mV;
            if (cq.ph > ph_thr && cq.pt >= kMuonTlo && cq.pt <= kMuonThi) {
                fq.planes[pl].signal_Q.push_back(cq.Q);
            }

            // F1: Q noise over channel-events with ph < 5 mV in [2500, 8000] ns
            double vmin_np = 1e30;
            for (int i = 0; i < (int)tv.size(); ++i) {
                if (tv[i] >= 2500.0 && tv[i] <= 8000.0 && (*cw[c])[i] < vmin_np)
                    vmin_np = (*cw[c])[i];
            }
            if ((cq.ped - vmin_np) < 5.0) {
                const int i0 = 262; // sample at 4200 ns
                double ped_Q_np = cq.ped;
                if (use_loc && i0 >= 12 && i0 - 4 < (int)tv.size()) {
                    double s_loc = 0;
                    for (int k = -12; k <= -4; ++k) s_loc += (*cw[c])[i0 + k];
                    ped_Q_np = s_loc / 9.0;
                }
                double Q_noise = 0;
                for (int k = i0 - 3; k <= i0 + npost; ++k) {
                    if (k >= 0 && k < (int)tv.size())
                        Q_noise += (ped_Q_np - (*cw[c])[k]) * kDt_ns;
                }
                fq.planes[pl].nopulse_Q.push_back(Q_noise);
            }

            // Acumula ruido
            fq.noise[c]+=cq.noise; fq.n_ped[c]++;

            // Histogramas de Q por plano
            fq.planes[pl].h_Qall->Fill(cq.Q);
            if(cq.Q>maxQ[pl]) maxQ[pl]=cq.Q;

            // Timing a umbral alto (Step 0)
            const double thr_hi=(pl>=P19X)?thr19_hi:kThrHigh18;
            if(cq.ph>thr_hi && cq.pt>0) fq.planes[pl].h_pt_hi->Fill(cq.pt);

            // Per canal
            if(make_per_channel && fq.h_Q_ch[c]) fq.h_Q_ch[c]->Fill(cq.Q);

            // Q vs ph y ratio r para Ch50
            if(c==50){
                fq.h_QvsPH->Fill(cq.ph, cq.Q);
                if(cq.ph>5 && cq.nw>0)
                    fq.h_ratio_r->Fill(cq.r);
            }
        }
        // Llena histograma max(Q) por plano
        for(int pl=0;pl<4;pl++)
            if(maxQ[pl]>-1e8) fq.planes[pl].h_Q->Fill(maxQ[pl]);
        // Guarda max(Q) para eficiencia
        for(int pl=0;pl<4;pl++)
            fq.planes[pl].Q_per_event.push_back(maxQ[pl]);
    }
    std::cout<<"    "<<nEnt<<"/"<<nEnt<<" ✓\n";
    fin->Close();

    // Normaliza ruido
    for(int c=0;c<kNCh;c++) if(fq.n_ped[c]>0){
        fq.noise[c]/=fq.n_ped[c];
        const int nw_c = kNpre + (c >= 32 ? kNpost19 : kNpost18) + 1;
        fq.sigma_Q_ch[c]=fq.noise[c]*kDt_ns*std::sqrt(nw_c);}
}

//==============================================================================
// Calcula umbral natural Q y fits
//==============================================================================
static void compute_fits(FileQ& fq){
    for(int pl=0;pl<4;pl++){
        auto& p=fq.planes[pl];
        TH1D* h=p.h_Q;
        if(!h||h->Integral()<5)continue;

        // F1: Pedestal noise of Q: RMS over channel-events with ph < 5 mV in [2500,8000] ns
        if(!p.nopulse_Q.empty()){
            double sum=0, sum2=0;
            for(double q:p.nopulse_Q) sum += q;
            double mean = sum / p.nopulse_Q.size();
            for(double q:p.nopulse_Q) sum2 += (q - mean) * (q - mean);
            p.sigma_Q_noise = std::sqrt(sum2 / p.nopulse_Q.size());
        }
        p.sigma_Q_ped = p.sigma_Q_noise;

        // Valle: mínimo en [umbral_ped, Q_max/2]
        const double Q_max=h->GetBinCenter(h->FindLastBinAbove(1.0));
        int blo=h->FindBin(p.sigma_Q_ped*3), bhi=h->FindBin(Q_max*0.5);
        int min_b = -1;
        double min_v = 1e30;
        for(int b=blo;b<=bhi;b++){
            double c=h->GetBinContent(b);
            if(c>0 && c<min_v){
                min_v=c;
                min_b=b;
            }
        }
        p.blo = blo; p.bhi = bhi; p.min_b = min_b;
        // F2: treat a minimum within 2 bins of blo or bhi as "no valley"
        if(min_b > blo + 2 && min_b < bhi - 2){
            p.has_valley = true;
            p.umbral_Q = h->GetBinCenter(min_b);
        } else {
            p.has_valley = false;
            p.umbral_Q = 0;
        }

        // Landau en la región de señal (solo si hay valle real)
        if(p.has_valley && p.umbral_Q>0 && Q_max>p.umbral_Q)
            p.lan_Q=fit_landau_Q(h,p.umbral_Q,std::min(Q_max,p.h_Q->GetXaxis()->GetXmax()-1));

        // Eficiencia al umbral Q natural
        const long long N=fq.n_events;
        if(N>0 && p.has_valley && p.umbral_Q>0){
            long long cnt=0;
            for(double Qmx:p.Q_per_event) if(Qmx>p.umbral_Q)cnt++;
            p.eff_nat_Q=100.0*cnt/N;
        } else {
            p.eff_nat_Q=0;
        }

        // Estimador de carga robusto (percentiles)
        p.n_sig_Q = (long long)p.signal_Q.size();
        if(!p.signal_Q.empty()){
            std::sort(p.signal_Q.begin(), p.signal_Q.end());
            auto get_pct = [&](double pct) {
                double idx = pct * (p.signal_Q.size() - 1);
                size_t i0 = static_cast<size_t>(std::floor(idx));
                size_t i1 = static_cast<size_t>(std::ceil(idx));
                double frac = idx - i0;
                return p.signal_Q[i0] * (1.0 - frac) + p.signal_Q[i1] * frac;
            };
            p.med_Q = get_pct(0.50);
            p.q10_Q = get_pct(0.10);
            p.q90_Q = get_pct(0.90);
        }
    }
}

//==============================================================================
// Figuras
//==============================================================================
static void make_figures(FileQ& fq, const FileQ& fq_other,
                          bool is_primary, const fs::path& out){
    const std::string L=fq.label;
    const int col=(L=="v7")?kBlue+1:kOrange+7;

    //── Step 0: Pulse time al umbral alto ─────────────────────────────────────
    if(is_primary){  // solo para v8 (más estadística)
        TCanvas c0("c_pth","",1200,900);c0.Divide(2,2);
        for(int pl=0;pl<4;pl++){
            c0.cd(pl+1);gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);
            TH1D* h=fq.planes[pl].h_pt_hi;
            h->SetLineColor(col);h->SetLineWidth(2);h->SetStats(0);
            h->Draw("HIST");
            // Fit gaussiana para determinar t_peak
            if(h->Integral()>5){
                GausFit gf=fit_gaus(h,h->GetMean()-h->GetRMS()*2,h->GetMean()+h->GetRMS()*2);
                if(gf.ok){
                    TF1 fg("_fg","gaus",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
                    fg.SetParameters(h->GetMaximum(),gf.mean,gf.sigma);
                    fg.SetLineColor(kRed);fg.SetLineWidth(2);fg.Draw("SAME");
                    TLatex lt;lt.SetNDC();lt.SetTextSize(0.038);
                    char buf[60];std::snprintf(buf,sizeof(buf),"t_{{peak}}=%.0f  #sigma=%.0f ns",gf.mean,gf.sigma);
                    lt.DrawLatex(0.15,0.92,buf);}
            }
        }
        c0.SaveAs((out/"pulsetime_highthreshold_all_planes_v8.png").string().c_str());
    }

    //── Step 2: Histogramas de carga por plano (v7 y v8 superpuestos) ─────────
    // Se genera en make_comparison_Q

    //── Step 5: Per-canal Q hod2019 X e Y ────────────────────────────────────
    if(is_primary){
        static const std::array<int,16> h19x={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
        static const std::array<int,16> h19y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};
        auto make_4x4=[&](const std::array<int,16>& chs,const std::string&nm,const fs::path&path){
            TCanvas cc(("cc"+nm).c_str(),"",1800,1200);cc.Divide(4,4);
            for(int k=0;k<16;k++){
                int c=chs[k];cc.cd(k+1);
                gPad->SetLogy();gPad->SetLeftMargin(0.15);gPad->SetBottomMargin(0.15);
                if(!fq.h_Q_ch[c])continue;
                fq.h_Q_ch[c]->SetLineColor(col);fq.h_Q_ch[c]->SetLineWidth(1);
                fq.h_Q_ch[c]->SetStats(0);
                char tt[24];std::snprintf(tt,sizeof(tt),"k=%d Ch%d",k,c);
                fq.h_Q_ch[c]->SetTitle(tt);fq.h_Q_ch[c]->Draw("HIST");
                if(fq.sigma_Q_ch[c]>0){
                    double vl=fq.sigma_Q_ch[c]*3;
                    gPad->Update();
                    TLine ln(vl,gPad->GetUymin()<0?std::pow(10,gPad->GetUymin()):0.5,
                              vl,std::pow(10,gPad->GetUymax()*0.9));
                    ln.SetLineColor(kRed);ln.SetLineStyle(2);ln.Draw();}}
            cc.SaveAs(path.string().c_str());};
        make_4x4(h19x,"19x_"+L,out/("charge_per_channel_hod2019X_"+L+".png"));
        make_4x4(h19y,"19y_"+L,out/("charge_per_channel_hod2019Y_"+L+".png"));
        // hod2018 X e Y
        std::array<int,16> h18x,h18y;
        for(int k=0;k<16;k++){h18x[k]=k;h18y[k]=16+k;}
        make_4x4(h18x,"18x_"+L,out/("charge_per_channel_hod2018X_"+L+".png"));
        make_4x4(h18y,"18y_"+L,out/("charge_per_channel_hod2018Y_"+L+".png"));
    }

    //── Step 3: Q vs ph + ratio r ────────────────────────────────────────────
    if(is_primary){
        {sl2(fq.h_QvsPH,"Pulse height [mV]","Q [mV#upointns]");
         TCanvas c3("c3","",900,820);c3.SetLeftMargin(0.13);c3.SetRightMargin(0.17);c3.SetBottomMargin(0.13);
         fq.h_QvsPH->Draw("COLZ");
         // Línea de ratio r=0.5 (señal tipo Landau)
         const double Nw = kNpre + kNpost19 + 1; const double dt = kDt_ns;
         // Q = r × ph × Nw × dt → Q = 0.5×ph×Nw×dt → línea en espacio QvsPH
         double ph_max=fq.h_QvsPH->GetXaxis()->GetXmax();
         TLine lr(0,0,ph_max,0.5*ph_max*Nw*dt);
         lr.SetLineColor(kRed);lr.SetLineWidth(2);lr.Draw("SAME");
         TLatex lt;lt.SetNDC();lt.SetTextSize(0.032);lt.SetTextColor(kRed);
         lt.DrawLatex(0.50,0.93,"r=0.5 (pulso SiPM)");
         c3.SaveAs((out/"charge_vs_ph_hod2019X_v8.png").string().c_str());}
        {sl(fq.h_ratio_r,"r = Q/(ph#timesN#timesdt)","Hits (ph>5 mV)");
         fq.h_ratio_r->SetLineColor(col);fq.h_ratio_r->SetLineWidth(2);
         TCanvas cr("cr","",900,680);cr.SetLeftMargin(0.12);cr.SetBottomMargin(0.13);
         fq.h_ratio_r->Draw("HIST");
         // Línea en r=0.5
         gPad->Update();
         TLine lr(0.5,0,0.5,fq.h_ratio_r->GetMaximum()*0.9);
         lr.SetLineColor(kRed);lr.SetLineWidth(2);lr.Draw();
         cr.SaveAs((out/"ratio_r_distribution_v8.png").string().c_str());}
    }
}

//==============================================================================
// Figura comparación charge_histograms_all_planes.png (8 paneles)
//==============================================================================
static void make_charge_histos(const FileQ& v7, const FileQ& v8, const fs::path& out){
    TCanvas cc("c_cq","",1600,1800);cc.Divide(2,4);
    const FileQ* fqs[2]={&v7,&v8};
    const int cols[2]={kBlue+1,kOrange+7};
    for(int pl=0;pl<4;pl++){
        for(int v=0;v<2;v++){
            cc.cd(pl*2+v+1);
            gPad->SetLogy();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
            const FileQ& fq=*fqs[v];
            TH1D* h=fq.planes[pl].h_Q;
            if(!h){gPad->Clear();continue;}
            h->SetLineColor(cols[v]);h->SetLineWidth(2);h->SetStats(0);
            sl(h,"Q_{max} [mV#upointns]","Entradas");
            h->Draw("HIST");
            // Línea umbral Q natural
            const double uQ=fq.planes[pl].umbral_Q;
            if(uQ>0){gPad->Update();
                TLine lu(uQ,std::pow(10,gPad->GetUymin()),uQ,std::pow(10,gPad->GetUymax()*0.9));
                lu.SetLineColor(kRed);lu.SetLineStyle(2);lu.SetLineWidth(2);lu.Draw();}
            // Fit pedestal (gaussiana)
            if(fq.planes[pl].gaus_ped.ok){
                TF1 fg("_fp","gaus",h->GetXaxis()->GetXmin(),uQ);
                fg.SetParameters(h->GetBinContent(h->FindBin(0)),
                                  fq.planes[pl].gaus_ped.mean,
                                  fq.planes[pl].gaus_ped.sigma);
                fg.SetLineColor(kGreen+2);fg.SetLineWidth(2);fg.Draw("SAME");}
            // Fit Landau en señal
            if(fq.planes[pl].lan_Q.ok){
                TF1 fl("_fl","TMath::Landau(x,[0],[1],0)",uQ,h->GetXaxis()->GetXmax());
                fl.SetParameters(fq.planes[pl].lan_Q.mpv,fq.planes[pl].lan_Q.sigma);
                fl.SetLineColor(kMagenta);fl.SetLineWidth(2);fl.Draw("SAME");}
            TLatex lt;lt.SetNDC();lt.SetTextSize(0.038);
            char buf[80];
            std::snprintf(buf,sizeof(buf),"%s %s",kPname[pl],fqs[v]->label.c_str());
            lt.DrawLatex(0.15,0.93,buf);
            if(fq.planes[pl].gaus_ped.ok){
                std::snprintf(buf,sizeof(buf),"#sigma_Q=%.0f",fq.planes[pl].gaus_ped.sigma);
                lt.DrawLatex(0.55,0.88,buf);}
            if(fq.planes[pl].lan_Q.ok){
                std::snprintf(buf,sizeof(buf),"MPV=%.0f",fq.planes[pl].lan_Q.mpv);
                lt.DrawLatex(0.55,0.84,buf);}
        }
    }
    cc.SaveAs((out/"charge_histograms_all_planes.png").string().c_str());
}

//==============================================================================
// Figura de eficiencia en Q
//==============================================================================
static void make_eff_Q(const FileQ& v7, const FileQ& v8, const fs::path& out){
    std::vector<double> Q_thr_v8={0,100,200,500,1000,1500,2000,3000};
    TCanvas cc("c_effQ","",900,680);cc.SetLeftMargin(0.12);cc.SetBottomMargin(0.13);
    const FileQ* fqs[2]={&v7,&v8};
    const int cols[2]={kBlue+1,kOrange+7};
    TH1D frame("fr","ε(Q_{thr}) hod2019X;Q_{umbral} [mV#upointns];#varepsilon",
               100,0,3500);
    frame.SetMinimum(0);frame.SetMaximum(1.05);frame.SetStats(0);frame.Draw();
    for(int v=0;v<2;v++){
        const FileQ&fq=*fqs[v];
        const auto& Qvec=fq.planes[P19X].Q_per_event;
        const long long N=fq.n_events;
        if(N==0||Qvec.empty())continue;
        std::vector<double> xs,ys,eys;
        for(double Qt:Q_thr_v8){
            long long cnt=0;
            for(double Qmx:Qvec)if(Qmx>Qt)cnt++;
            double eps=(double)cnt/N;
            xs.push_back(Qt);ys.push_back(eps);
            eys.push_back(std::sqrt(eps*(1-eps)/N));}
        TGraphErrors* g=new TGraphErrors((int)xs.size(),xs.data(),ys.data(),
            std::vector<double>(xs.size(),0).data(),eys.data());
        g->SetLineColor(cols[v]);g->SetMarkerColor(cols[v]);
        g->SetMarkerStyle(20+v);g->SetLineWidth(2);g->Draw("PL SAME");
        // Marca umbral natural en v8
        if(v==1){double uQ=fq.planes[P19X].umbral_Q;
            TLine lu(uQ,0,uQ,0.9);lu.SetLineColor(kRed);lu.SetLineStyle(2);lu.Draw("SAME");}
    }
    TLegend lg(0.55,0.65,0.88,0.78);lg.SetBorderSize(0);lg.SetFillStyle(0);
    lg.AddEntry((TObject*)nullptr,"v7","");lg.AddEntry((TObject*)nullptr,"v8","");lg.Draw();
    cc.SaveAs((out/"efficiency_Q_vs_threshold_hod2019X.png").string().c_str());
}

//==============================================================================
// Figura resumen de 6 paneles (charge_figures_of_merit.png)
//==============================================================================
static void make_summary_fig(const FileQ& v7, const FileQ& v8, const fs::path& out){
    TCanvas cc("c_fom","",2400,1600);cc.Divide(3,2);

    //── [1,1] hod2018X v8 con fits ────────────────────────────────────────────
    cc.cd(1);gPad->SetLogy();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
    {TH1D*h=v8.planes[P18X].h_Q;if(h){h->SetLineColor(kOrange+7);h->SetLineWidth(2);h->SetStats(0);sl(h,"Q_{max} [mV#upointns]","Entradas");h->Draw("HIST");
     if(v8.planes[P18X].gaus_ped.ok){TF1 fg("_f1","gaus",h->GetXaxis()->GetXmin(),v8.planes[P18X].umbral_Q);fg.SetParameters(h->GetBinContent(h->FindBin(0)),v8.planes[P18X].gaus_ped.mean,v8.planes[P18X].gaus_ped.sigma);fg.SetLineColor(kGreen+2);fg.SetLineWidth(2);fg.Draw("SAME");}
     if(v8.planes[P18X].lan_Q.ok){TF1 fl("_fl1","TMath::Landau(x,[0],[1],0)",v8.planes[P18X].umbral_Q,h->GetXaxis()->GetXmax());fl.SetParameters(v8.planes[P18X].lan_Q.mpv,v8.planes[P18X].lan_Q.sigma);fl.SetLineColor(kMagenta);fl.SetLineWidth(2);fl.Draw("SAME");}
     TLatex lt;lt.SetNDC();lt.SetTextSize(0.042);lt.DrawLatex(0.15,0.93,"hod2018X v8");
     if(v8.planes[P18X].lan_Q.ok){char buf[40];std::snprintf(buf,sizeof(buf),"MPV=%.0f mV·ns",v8.planes[P18X].lan_Q.mpv);lt.DrawLatex(0.40,0.85,buf);}}}

    //── [1,2] hod2019X v8 + flecha señal ─────────────────────────────────────
    cc.cd(2);gPad->SetLogy();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
    {TH1D*h=v8.planes[P19X].h_Q;if(h){h->SetLineColor(kOrange+7);h->SetLineWidth(2);h->SetStats(0);sl(h,"Q_{max} [mV#upointns]","Entradas");h->Draw("HIST");
     const double uQ=v8.planes[P19X].umbral_Q;
     if(uQ>0){gPad->Update();TLine lu(uQ,std::pow(10,gPad->GetUymin()),uQ,std::pow(10,gPad->GetUymax()*0.9));lu.SetLineColor(kRed);lu.SetLineStyle(2);lu.SetLineWidth(2);lu.Draw();}
     TLatex lt;lt.SetNDC();lt.SetTextSize(0.042);lt.DrawLatex(0.15,0.93,"hod2019X v8");}}

    //── [2,1] hod2019X v7 vs v8 normalizados ─────────────────────────────────
    cc.cd(3);gPad->SetLogy();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
    {TH1D* h7=(TH1D*)v7.planes[P19X].h_Q->Clone("_h7n");
     TH1D* h8=(TH1D*)v8.planes[P19X].h_Q->Clone("_h8n");
     if(h7->Integral()>0)h7->Scale(1.0/h7->Integral());
     if(h8->Integral()>0)h8->Scale(1.0/h8->Integral());
     h7->SetLineColor(kBlue+1);h7->SetLineWidth(2);h7->SetStats(0);
     h8->SetLineColor(kOrange+7);h8->SetLineWidth(2);
     sl(h7,"Q_{max} [mV#upointns]","Densidad");h7->Draw("HIST");h8->Draw("HIST SAME");
     TLegend lg(0.55,0.75,0.90,0.90);lg.SetBorderSize(0);lg.SetFillStyle(0);
     lg.AddEntry(h7,"v7","l");lg.AddEntry(h8,"v8","l");lg.Draw();
     TLatex lt;lt.SetNDC();lt.SetTextSize(0.042);lt.DrawLatex(0.15,0.93,"hod2019X v7 vs v8");}

    //── [2,2] Q vs ph TH2D ────────────────────────────────────────────────────
    cc.cd(4);gPad->SetLeftMargin(0.14);gPad->SetRightMargin(0.18);gPad->SetBottomMargin(0.14);
    {sl2(v8.h_QvsPH,"Pulse height [mV]","Q [mV#upointns]");
     v8.h_QvsPH->Draw("COLZ");
     TLatex lt;lt.SetNDC();lt.SetTextSize(0.042);lt.DrawLatex(0.15,0.93,"Ch50 hod2019X v8");}

    //── [3,1] Eficiencia vs Q para hod2019X ──────────────────────────────────
    cc.cd(5);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
    {TH1D fr("fr5","ε vs Q_{thr} hod2019X;Q [mV#upointns];#varepsilon",100,0,3000);
     fr.SetMinimum(0);fr.SetMaximum(1.05);fr.SetStats(0);fr.Draw();
     const FileQ* fqs[2]={&v7,&v8};const int cols[2]={kBlue+1,kOrange+7};
     const std::vector<double> Qt={0,100,200,500,1000,1500,2000};
     for(int v=0;v<2;v++){
         const auto& Qvec=fqs[v]->planes[P19X].Q_per_event;
         const long long N=fqs[v]->n_events;
         if(!N||Qvec.empty())continue;
         std::vector<double>xs,ys;
         for(double qt:Qt){long long cnt=0;for(double q:Qvec)if(q>qt)cnt++;xs.push_back(qt);ys.push_back((double)cnt/N);}
         TGraph*g=new TGraph((int)xs.size(),xs.data(),ys.data());
         g->SetLineColor(cols[v]);g->SetLineWidth(2);g->SetMarkerStyle(20+v);g->SetMarkerColor(cols[v]);g->Draw("PL SAME");}
     TLegend lg(0.55,0.65,0.88,0.78);lg.SetBorderSize(0);lg.SetFillStyle(0);
     lg.AddEntry((TObject*)nullptr,"v7","");lg.AddEntry((TObject*)nullptr,"v8","");lg.Draw();}

    //── [3,2] SNR bar plot horizontal ────────────────────────────────────────
    cc.cd(6);gPad->SetLeftMargin(0.25);gPad->SetBottomMargin(0.13);
    {// SNR_Q = MPV_Q / sigma_Q_pedestal; SNR_ph from previous analysis
     static const double snr_ph[4]={19.7,19.7,14.0,14.0};  // approx (ph/noise)
     TH1D hsnr("hsnr","SNR hod2018X;SNR;Plano",8,0.5,8.5);
     const char* labels[8]={"18X-ph","18X-Q","18Y-ph","18Y-Q","19X-ph","19X-Q","19Y-ph","19Y-Q"};
     // Fill SNR values (use v8 as primary)
     for(int pl=0;pl<4;pl++){
         hsnr.SetBinContent(pl*2+1, snr_ph[pl]);  // SNR_ph
         double snrQ=0;
         if(v8.planes[pl].lan_Q.ok&&v8.planes[pl].sigma_Q_noise>0)
             snrQ=v8.planes[pl].lan_Q.mpv/v8.planes[pl].sigma_Q_noise;
         hsnr.SetBinContent(pl*2+2, snrQ);
     }
     for(int i=0;i<8;i++)hsnr.GetXaxis()->SetBinLabel(i+1,labels[i]);
     hsnr.SetStats(0);hsnr.SetFillColor(kAzure+2);hsnr.SetTitle("SNR: pulse_height vs carga;SNR;");
     hsnr.Draw("HBAR");
     TLatex lt;lt.SetNDC();lt.SetTextSize(0.038);lt.DrawLatex(0.15,0.93,"SNR v8 (MPV/#sigma_{Q noise})");}

    cc.SaveAs((out/"charge_figures_of_merit.png").string().c_str());
}

//==============================================================================
// main
//==============================================================================
int main(int argc, char** argv){
    std::string fv7="lectura_vx1740_ext_nim_1M_v7.root";
    std::string fv8="lectura_vx1740_ext_nim_1M_v8.root";
    for(int i=1;i<argc;i++){
        std::string a=argv[i];
        if(a=="--v7"&&i+1<argc){fv7=argv[++i];continue;}
        if(a=="--v8"&&i+1<argc){fv8=argv[++i];continue;}
    }

    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("charge_analysis_"+std::string(ds));
    fs::create_directories(outdir);

    std::cout<<"=== charge_analysis ===\n";
    std::cout<<"v7: "<<fv7<<"\nv8: "<<fv8<<"\nOut: "<<outdir<<"\n";
    std::cout<<"Δt_sample="<<kDt_ns<<" ns  N_pre="<<kNpre
             <<"  N_post18="<<kNpost18<<"  N_post19="<<kNpost19<<"\n\n";

    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    FileQ fv7_data, fv8_data;
    fv7_data.label="v7"; fv8_data.label="v8";

    std::cout<<"[1/2] Procesando v7...\n";
    process_file_Q(fv7_data, fv7, kThrHigh19_v7, false);
    std::cout<<"[2/2] Procesando v8...\n";
    process_file_Q(fv8_data, fv8, kThrHigh19_v8, true);   // per-channel solo v8
    std::cout<<"Pasadas completas.\n\n";

    // Fits y umbrales
    compute_fits(fv7_data);
    compute_fits(fv8_data);

    // Genera figuras
    make_figures(fv7_data, fv8_data, false, outdir);
    make_figures(fv8_data, fv7_data, true,  outdir);
    make_charge_histos(fv7_data, fv8_data, outdir);
    make_eff_Q(fv7_data, fv8_data, outdir);
    make_summary_fig(fv7_data, fv8_data, outdir);

    // ROOT file
    {TFile fout((outdir/"histos_charge.root").string().c_str(),"RECREATE");
     if(!fout.IsZombie()){
         for(auto* fq : {&fv7_data,&fv8_data}){
             for(int pl=0;pl<4;pl++){
                 if(fq->planes[pl].h_Q)    fq->planes[pl].h_Q->Write();
                 if(fq->planes[pl].h_pt_hi)fq->planes[pl].h_pt_hi->Write();}
             if(fq->h_QvsPH)  fq->h_QvsPH->Write();
             if(fq->h_ratio_r)fq->h_ratio_r->Write();}
         fout.Close();}}

    // Tabla summary
    auto print_table=[&](std::ostream& os){
        os<<"=== Análisis de carga integrada Q ===\n\n";
        os<<"Parámetros de integración:\n";
        os<<"  Δt_sample = "<<kDt_ns<<" ns  N_pre="<<kNpre
          <<"  N_post18="<<kNpost18<<" (hod2018)  N_post19="<<kNpost19<<" (hod2019)\n\n";
        const int kW = 40;
        os<<std::setw(kW)<<std::left<<"Métrica"
          <<std::setw(14)<<"18X-v7"<<std::setw(14)<<"18X-v8"
          <<std::setw(14)<<"19X-v7"<<std::setw(14)<<"19X-v8\n";
        os<<"---------------------------------------------------------------------------------------------\n";
        auto fmt=[](double v,int p=1)->std::string{if(v==0||!std::isfinite(v))return "N/A";
            char b[20];std::snprintf(b,sizeof(b),"%.*f",p,v);return b;};

        os<<std::setw(kW)<<"σ_Q noise [mV·ns] (no-pulse channels)";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++)
            os<<std::setw(14)<<fmt((v==0?fv7_data:fv8_data).planes[pp].sigma_Q_noise);
        os<<"\n";

        os<<std::setw(kW)<<"N pulsos señal";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++)
            os<<std::setw(14)<<(v==0?fv7_data.planes[pp].n_sig_Q:fv8_data.planes[pp].n_sig_Q);
        os<<"\n";

        os<<std::setw(kW)<<"Mediana Q [mV·ns]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++)
            os<<std::setw(14)<<fmt(v==0?fv7_data.planes[pp].med_Q:fv8_data.planes[pp].med_Q);
        os<<"  (Q of pulses above the ph threshold; not a gain estimator)\n";

        os<<std::setw(kW)<<"Percentil 10 Q [mV·ns]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++)
            os<<std::setw(14)<<fmt(v==0?fv7_data.planes[pp].q10_Q:fv8_data.planes[pp].q10_Q);
        os<<"\n";

        os<<std::setw(kW)<<"Percentil 90 Q [mV·ns]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++)
            os<<std::setw(14)<<fmt(v==0?fv7_data.planes[pp].q90_Q:fv8_data.planes[pp].q90_Q);
        os<<"\n";

        os<<std::setw(kW)<<"MPV_Q Landau [mV·ns]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.lan_Q.ok || p.lan_Q.status != 0 || p.lan_Q.at_limit)
                os<<std::setw(14)<<"n/a";
            else
                os<<std::setw(14)<<fmt(p.lan_Q.mpv);
        }
        os<<"\n";

        os<<std::setw(kW)<<"Landau fit status";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.has_valley) os<<std::setw(14)<<"no fit";
            else os<<std::setw(14)<<p.lan_Q.status;
        }
        os<<"\n";

        os<<std::setw(kW)<<"Landau chi2/ndf";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.lan_Q.ok || p.lan_Q.status != 0) os<<std::setw(14)<<"n/a";
            else os<<std::setw(14)<<fmt(p.lan_Q.chi2ndf, 2);
        }
        os<<"\n";

        os<<std::setw(kW)<<"Landau limit flag";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.lan_Q.ok || p.lan_Q.status != 0) os<<std::setw(14)<<"n/a";
            else if(p.lan_Q.at_limit) os<<std::setw(14)<<"AT_LIMIT";
            else os<<std::setw(14)<<"interior";
        }
        os<<"\n";

        // SNR_Q
        os<<std::setw(kW)<<"SNR_Q = MPV/σ_Q";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.lan_Q.ok || p.lan_Q.status != 0 || p.lan_Q.at_limit || p.sigma_Q_noise <= 0)
                os<<std::setw(14)<<"n/a";
            else
                os<<std::setw(14)<<fmt(p.lan_Q.mpv / p.sigma_Q_noise);
        }
        os<<"\n";

        os<<std::setw(kW)<<"umbral_Q natural [mV·ns]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.has_valley)
                os<<std::setw(14)<<"no valley";
            else
                os<<std::setw(14)<<fmt(p.umbral_Q);
        }
        os<<"\n";

        os<<std::setw(kW)<<"  valley bins [blo,bhi,min]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            char bbuf[32];
            std::snprintf(bbuf, sizeof(bbuf), "[%d,%d,%d]", p.blo, p.bhi, p.min_b);
            os<<std::setw(14)<<bbuf;
        }
        os<<"\n";

        os<<std::setw(kW)<<"ε(umbral_Q_nat) [%]";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            if(!p.has_valley)
                os<<std::setw(14)<<"n/a";
            else
                os<<std::setw(14)<<fmt(p.eff_nat_Q, 2);
        }
        os<<"\n\n";

        os<<"Búsqueda de valle (umbral natural):\n";
        for(int pp:{P18X,P19X})for(int v=0;v<2;v++){
            const auto& p = (v==0?fv7_data:fv8_data).planes[pp];
            const char* name = (pp==P18X) ? (v==0?"18X-v7":"18X-v8") : (v==0?"19X-v7":"19X-v8");
            os<<"  "<<std::setw(8)<<std::left<<name<<": ";
            if(p.has_valley) os<<fmt(p.umbral_Q)<<" mV·ns";
            else os<<"no valley";
            os<<" (blo="<<p.blo<<", bhi="<<p.bhi<<", min="<<p.min_b;
            if(p.min_b <= p.blo + 2 || p.min_b >= p.bhi - 2) os<<", edge guard -> no valley";
            os<<")\n";
        }
        os<<"\n";
    };
    print_table(std::cout);
    {std::ofstream sf(outdir/"summary_charge.txt");print_table(sf);}
    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
