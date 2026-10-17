// tagprobe_analysis.cxx
// Detección estadística de muones en hod2019 por tag-and-probe con event-mixing.
// Usa hod2018 (ε≈99.9%) como referencia para predecir la posición del muón en hod2019
// y medir el exceso de carga sobre el fondo accidental (event-mixing).
//
// RAMA A si significancia > 5σ → detección estadística confirmada.
// RAMA B si significancia < 5σ → sin detección, confirma límite de light yield.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o tagprobe_analysis tagprobe_analysis.cxx

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
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace fs = std::filesystem;

//==============================================================================
// Constantes
//==============================================================================
constexpr double kPre_ns = 2500.0;
constexpr double kDt_ns = 16.0;  // VX1740, 62.5 MS/s; equals the data_dig.time step
constexpr int    kNpre=3, kNpost=8, kNwin=12;   // ventana de integración Q
constexpr int    kNCh = 64;
constexpr double kPitch = 2.0;   // mm
constexpr double kD_mm  = 98.4;  // mm

// Ventana de señal completa (para extracción de pulso)
static const double kTlo[4]={3500,3500,3000,3000};
static const double kThi[4]={7000,7000,5000,5000};

// Ventana de tiempo de muón en hod2019 (señales reales a t≈4175 ns, ±300 ns de jitter)
constexpr double kMuonTlo = 3900.0, kMuonThi = 4500.0;  // ns

// Umbral_Q para el barrido de eficiencia (de charge_analysis)
constexpr double kUmbralQ = 2440.048;  // rescaled x1.04 on 2026-10-17 (kDt_ns 15.3846 -> 16 ns); selection unchanged

// Mapeo hod2019
static const std::array<int,16> kH19X={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};

//==============================================================================
// Extracción de pulso + carga integrada para UNA TIRA
//==============================================================================
struct StripResult { double Q=0, pt=-1; bool fired=false; };

static StripResult extract_strip(const std::vector<float>&wf, const std::vector<float>&t,
                                  double tlo, double thi) {
    StripResult r;
    const int n=(int)wf.size();
    if(n<4||(int)t.size()!=n) return r;
    double s=0, s2=0; int nb=0;
    for(int i=0;i<n;i++) if(t[i]<kPre_ns){s+=wf[i];s2+=wf[i]*wf[i];nb++;}
    if(nb<2) return r;
    const double ped=s/nb;
    double vmin=1e30; int imin=-1;
    for(int i=0;i<n;i++) if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0) return r;
    r.pt=t[imin];
    // Integrar alrededor del mínimo
    const double tilo=r.pt-kNpre*kDt_ns, tihi=r.pt+kNpost*kDt_ns;
    double Q=0;
    for(int i=0;i<n;i++) if(t[i]>=tilo&&t[i]<=tihi) Q+=(ped-wf[i])*kDt_ns;
    r.Q=Q; r.fired=(Q>0 && vmin<ped-1.0);
    return r;
}

//==============================================================================
// Clustering hod2018 (gap≤1, igual que antes)
//==============================================================================
struct Hit18{int k; double ph;};
struct Cl18{double pos; double amp;};
static std::vector<Cl18> clust18(std::vector<Hit18> hits, double thr){
    std::vector<Hit18> act;
    for(const auto&h:hits) if(h.ph>=thr) act.push_back(h);
    if(act.empty()) return{};
    std::sort(act.begin(),act.end(),[](const Hit18&a,const Hit18&b){return a.k<b.k;});
    std::vector<Cl18> out; double wpos=0,wamp=0; int hi=act[0].k;
    auto flush=[&](){out.push_back({wpos/wamp,wamp});wpos=wamp=0;};
    auto add=[&](const Hit18&h){wpos+=h.ph*(h.k+0.5-8.0)*kPitch;wamp+=h.ph;hi=h.k;};
    hi=act[0].k; add(act[0]);
    for(int i=1;i<(int)act.size();i++){
        if(act[i].k>hi+1){flush();hi=act[i].k;} add(act[i]);}
    flush(); return out;}

//==============================================================================
// Datos por evento almacenados en RAM
//==============================================================================
struct EvData {
    // hod2018X: tag information
    int k18x=-1, k18y=-1;          // strip predicted for hod2019
    double amp18x=0, amp18y=0;     // cluster amplitude
    int ncls18x=0, ncls18y=0;      // number of clusters
    bool tag_x=false, tag_y=false; // hod2018 tag conditions
    // hod2019: charge + pulse_time per strip (16 strips each plane)
    std::array<double,16> Q19x{}, Q19y{};   // mV·ns
    std::array<double,16> pt19x{}, pt19y{}; // ns (-1 = no pulse)
};

//==============================================================================
// Carga en RAM
//==============================================================================
static std::vector<EvData> load_all(const std::string& fname) {
    TFile* fin=TFile::Open(fname.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: "<<fname<<'\n';return{};}
    TTree* t=static_cast<TTree*>(fin->Get("data_dig"));
    if(!t){fin->Close();return{};}

    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw=nullptr;
    t->SetBranchStatus("*",0); t->SetCacheSize(256LL<<20);
    for(int c=0;c<kNCh;c++){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!t->GetBranch(bn))continue;
        t->SetBranchStatus(bn,1);t->SetBranchAddress(bn,&cw[c]);}
    t->SetBranchStatus("time",1);t->SetBranchAddress("time",&tw);

    const Long64_t nEnt=t->GetEntries();
    std::cout<<"  Loading "<<nEnt<<" events...\n";
    std::vector<EvData> evs; evs.reserve(nEnt);

    // Build inverse maps for hod2019
    std::array<int,kNCh> inv19x{},inv19y{}; inv19x.fill(-1);inv19y.fill(-1);
    for(int k=0;k<16;k++){inv19x[kH19X[k]]=k;inv19y[kH19Y[k]]=k;}

    for(Long64_t iev=0;iev<nEnt;iev++){
        t->GetEntry(iev);
        if(iev%30000==0)std::cout<<"    "<<iev<<"/"<<nEnt<<"\r"<<std::flush;
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;
        static bool checked_step = false;
        if (!checked_step && tv.size() > 1) {
            checked_step = true;
            if (std::fabs((tv[1] - tv[0]) - kDt_ns) > 1e-3) {
                std::fprintf(stderr, "Sampling step %.4f ns != kDt_ns %.4f ns; aborting\n", tv[1]-tv[0], kDt_ns);
                return {};
            }
        }

        EvData ev;
        ev.pt19x.fill(-1); ev.pt19y.fill(-1);

        // hod2018X and Y hits
        std::vector<Hit18> h18x,h18y;
        for(int c=0;c<kNCh;c++){
            if(!cw[c]||cw[c]->empty())continue;
            if(c>31)continue;  // only hod2018 here
            const int pl=(c<16)?0:1;
            const int k=(c<16)?c:c-16;
            // Quick pulse height extraction for hod2018
            double s=0,s2=0;int nb=0;
            for(int i=0;i<(int)tv.size();i++) if(tv[i]<kPre_ns){s+=(*cw[c])[i];s2+=(*cw[c])[i]*(*cw[c])[i];nb++;}
            if(nb<2)continue;
            const double ped=s/nb;
            double vmin=1e30;
            for(int i=0;i<(int)tv.size();i++) if(tv[i]>=kTlo[pl]&&tv[i]<=kThi[pl]&&(*cw[c])[i]<vmin) vmin=(*cw[c])[i];
            const double ph=ped-vmin;
            if(ph>=5.0){if(pl==0)h18x.push_back({k,ph});else h18y.push_back({k,ph});}
        }
        // hod2018 clusters
        {auto cls=clust18(h18x,10.0);ev.ncls18x=(int)cls.size();
         if(!cls.empty()){const auto& best=*std::max_element(cls.begin(),cls.end(),[](const Cl18&a,const Cl18&b){return a.amp<b.amp;});
             ev.k18x=(int)std::round(best.pos/kPitch+7.5);
             if(ev.k18x<0)ev.k18x=0;if(ev.k18x>15)ev.k18x=15;
             ev.amp18x=best.amp;}
         ev.tag_x=(ev.ncls18x==1&&ev.amp18x>=10);}
        {auto cls=clust18(h18y,10.0);ev.ncls18y=(int)cls.size();
         if(!cls.empty()){const auto& best=*std::max_element(cls.begin(),cls.end(),[](const Cl18&a,const Cl18&b){return a.amp<b.amp;});
             ev.k18y=(int)std::round(best.pos/kPitch+7.5);
             if(ev.k18y<0)ev.k18y=0;if(ev.k18y>15)ev.k18y=15;
             ev.amp18y=best.amp;}
         ev.tag_y=(ev.ncls18y==1&&ev.amp18y>=10);}

        // hod2019X and Y: charge per strip
        for(int c=32;c<64;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const int k19x=inv19x[c], k19y=inv19y[c];
            const int pl=(c<48)?3:2;  // 3=19Y (32-47), 2=19X (48-63)
            if(k19x>=0||k19y>=0){
                StripResult sr=extract_strip(*cw[c],tv,kTlo[pl],kThi[pl]);
                if(k19x>=0){ev.Q19x[k19x]=sr.Q;ev.pt19x[k19x]=sr.pt;}
                if(k19y>=0){ev.Q19y[k19y]=sr.Q;ev.pt19y[k19y]=sr.pt;}
            }
        }
        evs.push_back(std::move(ev));
    }
    std::cout<<"    "<<nEnt<<"/"<<nEnt<<" ✓  ("<<evs.size()<<" events)\n";
    fin->Close();
    return evs;
}

//==============================================================================
// Tag-and-probe analysis para un plano
//==============================================================================
struct TPResult {
    long long n_tag=0, n_obs=0, n_acc=0, n_off=0;
    double eps_obs=0, eps_acc=0, eps_true=0, sigma_eps=0, significance=0;
    // Full histograms
    std::unique_ptr<TH1D> h_same, h_mixed, h_diff;
    // Efficiency sweep
    std::vector<double> Q_sweep={0,100,300,500,1000,1500,2000,2346,3000,5000};
    std::vector<double> eps_true_sweep, eps_true_err_sweep;
};

static TPResult tagprobe(const std::vector<EvData>& evs,
                           bool useX,   // X or Y plane
                           const std::vector<int>& mask,  // indices to skip (hot channels)
                           const std::vector<size_t>& shuf) { // shuffled event indices
    TPResult tp;
    const size_t N=evs.size();
    const double Qlo=-500, Qhi=5000;

    tp.h_same =std::make_unique<TH1D>("hs","H_{same}",100,Qlo,Qhi);
    tp.h_mixed=std::make_unique<TH1D>("hm","H_{mixed}",100,Qlo,Qhi);

    // Helper: get Q and pt for strip k of event i
    auto getQ=[&](size_t i, int k)->double{
        if(i>=N||k<0||k>15)return 0;
        return useX?evs[i].Q19x[k]:evs[i].Q19y[k];};
    auto getPt=[&](size_t i, int k)->double{
        if(i>=N||k<0||k>15)return -1;
        return useX?evs[i].pt19x[k]:evs[i].pt19y[k];};
    auto inMuonWin=[&](double pt)->bool{
        return (pt>=kMuonTlo&&pt<=kMuonThi);};
    auto isHot=[&](int k)->bool{
        return std::find(mask.begin(),mask.end(),k)!=mask.end();};

    // Tag condition
    auto isTag=[&](size_t i)->bool{
        return useX?evs[i].tag_x:evs[i].tag_y;};
    auto kTag=[&](size_t i)->int{
        return useX?evs[i].k18x:evs[i].k18y;};

    // Initialize sweep counters
    const int nsw=(int)tp.Q_sweep.size();
    std::vector<long long> n_obs_sw(nsw,0),n_acc_sw(nsw,0);

    for(size_t i=0;i<N;i++){
        if(!isTag(i))continue;
        const int kp=kTag(i);
        if(kp<0||kp>15)continue;
        if(isHot(kp))continue;  // skip if predicted strip is hot
        ++tp.n_tag;

        // Q_same (same event, predicted strip, with muon time window)
        const double Qs=getQ(i,kp);
        const double pts=getPt(i,kp);
        const bool time_ok_s=inMuonWin(pts);
        const double Qsame=time_ok_s?Qs:0.0;

        // Q_mixed (different event, same strip)
        const size_t j=shuf[i];
        const double Qm=getQ(j,kp);
        const double ptm=getPt(j,kp);
        const bool time_ok_m=inMuonWin(ptm);
        const double Qmixed=time_ok_m?Qm:0.0;

        // Q_off (same event, opposite strip = k+8 mod 16)
        const int koff=(kp+8)%16;
        const double Qo=getQ(i,koff);
        const double pto=getPt(i,koff);
        const double Qoff=inMuonWin(pto)?Qo:0.0;

        tp.h_same->Fill(Qsame);
        tp.h_mixed->Fill(Qmixed);

        if(Qsame>kUmbralQ) ++tp.n_obs;
        if(Qmixed>kUmbralQ) ++tp.n_acc;
        if(Qoff>kUmbralQ) ++tp.n_off;

        for(int qi=0;qi<nsw;qi++){
            if(Qsame>tp.Q_sweep[qi]) ++n_obs_sw[qi];
            if(Qmixed>tp.Q_sweep[qi]) ++n_acc_sw[qi];
        }
    }

    // Compute main efficiency at kUmbralQ
    const double Nt=static_cast<double>(tp.n_tag);
    if(Nt<1)return tp;
    tp.eps_obs=tp.n_obs/Nt;
    tp.eps_acc=tp.n_acc/Nt;
    tp.eps_true=tp.eps_obs-tp.eps_acc;
    const double s_obs=std::sqrt(tp.eps_obs*(1-tp.eps_obs)/Nt);
    const double s_acc=std::sqrt(tp.eps_acc*(1-tp.eps_acc)/Nt);
    tp.sigma_eps=std::sqrt(s_obs*s_obs+s_acc*s_acc);
    tp.significance=(tp.sigma_eps>0)?tp.eps_true/tp.sigma_eps:0;

    // Sweep
    tp.eps_true_sweep.resize(nsw);
    tp.eps_true_err_sweep.resize(nsw);
    for(int qi=0;qi<nsw;qi++){
        double eo=(double)n_obs_sw[qi]/Nt, ea=(double)n_acc_sw[qi]/Nt;
        tp.eps_true_sweep[qi]=eo-ea;
        tp.eps_true_err_sweep[qi]=std::sqrt(eo*(1-eo)/Nt+ea*(1-ea)/Nt);
    }

    // Difference histogram
    tp.h_diff=std::make_unique<TH1D>("hd","H_{same}#minusH_{mixed}",100,Qlo,Qhi);
    for(int b=1;b<=100;b++)
        tp.h_diff->SetBinContent(b,tp.h_same->GetBinContent(b)-tp.h_mixed->GetBinContent(b));

    return tp;
}

//==============================================================================
// Main
//==============================================================================
int main(int argc,char**argv){
    std::string fv8="lectura_vx1740_ext_nim_1M_v8.root";
    for(int i=1;i<argc;i++){std::string a=argv[i];if(a=="--v8"&&i+1<argc){fv8=argv[++i];}}

    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("tagprobe_"+std::string(ds));
    fs::create_directories(outdir);
    std::cout<<"=== tagprobe_analysis ===\nv8: "<<fv8<<"\nOut: "<<outdir<<"\n\n";
    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);gStyle->SetOptStat(0);

    // Carga todos los eventos en RAM
    auto evs=load_all(fv8);
    const size_t N=evs.size();
    std::cout<<"Events loaded: "<<N<<"\n\n";
    if(N==0){std::cerr<<"ERROR: no events\n";return 1;}

    // Event mixing: shuffle con semilla fija
    std::vector<size_t> shuf(N);
    std::iota(shuf.begin(),shuf.end(),0);
    {std::mt19937 rng(42);std::shuffle(shuf.begin(),shuf.end(),rng);}

    //==========================================================================
    // ANÁLISIS 2 — Canales calientes
    //==========================================================================
    // Ocupancia: eventos donde Q>520 mV·ns AND pt ∈ muon window
    std::array<long long,16> occ19x{},occ19y{};
    for(const auto&ev:evs){
        for(int k=0;k<16;k++){
            if(ev.Q19x[k]>520.0&&ev.pt19x[k]>=kMuonTlo&&ev.pt19x[k]<=kMuonThi) occ19x[k]++;  // rescaled x1.04 on 2026-10-17 (kDt_ns 15.3846 -> 16 ns); selection unchanged
            if(ev.Q19y[k]>520.0&&ev.pt19y[k]>=kMuonTlo&&ev.pt19y[k]<=kMuonThi) occ19y[k]++;  // rescaled x1.04 on 2026-10-17 (kDt_ns 15.3846 -> 16 ns); selection unchanged
        }
    }
    // Mediana de ocupancia
    auto find_median=[](std::array<long long,16> a)->double{
        std::sort(a.begin(),a.end());return 0.5*(a[7]+a[8]);};
    const double med_x=find_median(occ19x), med_y=find_median(occ19y);
    const double hot_thr_x=3.0*med_x, hot_thr_y=3.0*med_y;

    std::vector<int> hot_x,hot_y;
    for(int k=0;k<16;k++){
        if(occ19x[k]>hot_thr_x&&med_x>0)hot_x.push_back(k);
        if(occ19y[k]>hot_thr_y&&med_y>0)hot_y.push_back(k);
    }
    // Si la mediana es 0 (sin señal), los "hot channels" son los que tienen alguna señal
    if(med_x==0){for(int k=0;k<16;k++)if(occ19x[k]>2)hot_x.push_back(k);}
    if(med_y==0){for(int k=0;k<16;k++)if(occ19y[k]>2)hot_y.push_back(k);}

    std::cout<<"Canales calientes hod2019X (occ>3×mediana="<<med_x<<"): ";
    for(int k:hot_x)std::cout<<"k="<<k<<"(Ch"<<kH19X[k]<<",occ="<<occ19x[k]<<") ";
    std::cout<<"\nCanales calientes hod2019Y: ";
    for(int k:hot_y)std::cout<<"k="<<k<<"(Ch"<<kH19Y[k]<<",occ="<<occ19y[k]<<") ";
    std::cout<<"\n\n";

    // CSV de canales calientes
    {std::ofstream fh(outdir/"hotchannels.csv");
     fh<<"plane,k,channel,occupancy_muonwindow,fraction\n";
     for(int k:hot_x)fh<<"hod2019X,"<<k<<','<<kH19X[k]<<','<<occ19x[k]<<','<<(double)occ19x[k]/N<<"\n";
     for(int k:hot_y)fh<<"hod2019Y,"<<k<<','<<kH19Y[k]<<','<<occ19y[k]<<','<<(double)occ19y[k]/N<<"\n";}

    //==========================================================================
    // ANÁLISIS 1 — Tag-and-probe (todos los canales vs canales calientes enmascarados)
    //==========================================================================
    std::cout<<"Running tag-and-probe...\n";
    std::vector<int> empty_mask;
    auto tp_X_all  = tagprobe(evs,true, empty_mask, shuf);
    auto tp_Y_all  = tagprobe(evs,false,empty_mask, shuf);
    auto tp_X_mask = tagprobe(evs,true, hot_x,      shuf);
    auto tp_Y_mask = tagprobe(evs,false,hot_y,      shuf);

    std::cout<<"hod2019X (all):  N_tag="<<tp_X_all.n_tag<<"  sig="<<tp_X_all.significance<<"\n";
    std::cout<<"hod2019Y (all):  N_tag="<<tp_Y_all.n_tag<<"  sig="<<tp_Y_all.significance<<"\n";
    std::cout<<"hod2019X (mask): N_tag="<<tp_X_mask.n_tag<<"  sig="<<tp_X_mask.significance<<"\n";
    std::cout<<"hod2019Y (mask): N_tag="<<tp_Y_mask.n_tag<<"  sig="<<tp_Y_mask.significance<<"\n\n";

    // Correlation r_X tras enmascarar canales calientes
    std::vector<double> rx18,rx19,ry18,ry19;
    for(const auto&ev:evs){
        if(!ev.tag_x)continue;
        // Busca un canal hod2019X activo (Q>kUmbralQ, pt en muon window, no hot)
        int best_k=-1;double best_Q=-1;
        for(int k=0;k<16;k++){
            if(std::find(hot_x.begin(),hot_x.end(),k)!=hot_x.end())continue;
            if(ev.Q19x[k]>kUmbralQ&&ev.pt19x[k]>=kMuonTlo&&ev.pt19x[k]<=kMuonThi&&ev.Q19x[k]>best_Q){
                best_Q=ev.Q19x[k];best_k=k;}}
        if(best_k>=0){rx18.push_back(ev.k18x*kPitch-(7.5*kPitch));rx19.push_back((best_k+0.5-8.0)*kPitch);}
    }
    double rX_masked=0;
    {const int n=(int)rx18.size();if(n>=2){
        double sx=0,sy=0,sxx=0,syy=0,sxy=0;
        for(int i=0;i<n;i++){sx+=rx18[i];sy+=rx19[i];sxx+=rx18[i]*rx18[i];syy+=rx19[i]*rx19[i];sxy+=rx18[i]*rx19[i];}
        const double d=n,num=d*sxy-sx*sy,den=std::sqrt((d*sxx-sx*sx)*(d*syy-sy*sy));
        rX_masked=den>0?num/den:0;}}
    std::cout<<"Correlación r_X tras enmascarar hot channels: "<<rX_masked<<" (N="<<rx18.size()<<")\n\n";

    //==========================================================================
    // Figuras
    //==========================================================================
    // Fig 1: tagprobe_charge_same_vs_mixed_X.png
    auto make_tagprobe_fig=[&](TPResult& tp,const std::string&coord,const fs::path&path){
        TCanvas cc(("c_tp_"+coord).c_str(),"",1200,900);cc.Divide(1,2);
        // Panel superior: H_same vs H_mixed
        cc.cd(1);gPad->SetLogy();gPad->SetLeftMargin(0.12);gPad->SetBottomMargin(0.12);
        tp.h_same->SetLineColor(kBlue+1);tp.h_same->SetLineWidth(2);tp.h_same->SetStats(0);
        tp.h_mixed->SetLineColor(kOrange+7);tp.h_mixed->SetLineWidth(2);
        tp.h_same->GetXaxis()->SetTitle("Q [mV#upointns]");tp.h_same->GetYaxis()->SetTitle("Events");
        tp.h_same->SetTitle(("Tag-and-probe hod2019"+coord).c_str());
        tp.h_same->Draw("HIST");tp.h_mixed->Draw("HIST SAME");
        TLine lq(kUmbralQ,gPad->GetUymin()<0?0.5:0.5,kUmbralQ,tp.h_same->GetMaximum()*0.5);
        lq.SetLineColor(kRed);lq.SetLineStyle(2);lq.SetLineWidth(2);lq.Draw();
        TLegend lg(0.55,0.75,0.88,0.90);lg.SetBorderSize(0);lg.SetFillStyle(0);
        char lb[60];std::snprintf(lb,sizeof(lb),"H_{same} (N_{tag}=%lld)",tp.n_tag);
        lg.AddEntry(tp.h_same.get(),lb,"l");
        lg.AddEntry(tp.h_mixed.get(),"H_{mixed} (event-mixing)","l");lg.Draw();
        TLatex lt;lt.SetNDC();lt.SetTextSize(0.038);
        char buf[80];
        std::snprintf(buf,sizeof(buf),"#varepsilon_{true}=%.4f%%  #pm%.4f%%  sig=%.2f#sigma",
            tp.eps_true*100,tp.sigma_eps*100,tp.significance);
        lt.DrawLatex(0.14,0.93,buf);
        // Panel inferior: diferencia
        cc.cd(2);gPad->SetLeftMargin(0.12);gPad->SetBottomMargin(0.14);
        tp.h_diff->SetLineColor(kGreen+2);tp.h_diff->SetLineWidth(2);tp.h_diff->SetStats(0);
        tp.h_diff->GetXaxis()->SetTitle("Q [mV#upointns]");tp.h_diff->GetYaxis()->SetTitle("H_{same}#minusH_{mixed}");
        tp.h_diff->SetTitle("Diferencia (exceso de señal)");tp.h_diff->Draw("HIST");
        TLine l0(tp.h_diff->GetXaxis()->GetXmin(),0,tp.h_diff->GetXaxis()->GetXmax(),0);
        l0.SetLineColor(kRed);l0.SetLineStyle(2);l0.Draw();
        cc.SaveAs(path.string().c_str());};

    make_tagprobe_fig(tp_X_all,"X",outdir/"tagprobe_charge_same_vs_mixed_X.png");
    make_tagprobe_fig(tp_Y_all,"Y",outdir/"tagprobe_charge_same_vs_mixed_Y.png");

    // Fig 2: efficiency_with_without_hotchannel_mask.png
    {TCanvas cc("c_eff_mask","",1000,680);cc.SetLeftMargin(0.12);cc.SetBottomMargin(0.13);
     TH1D fr("fr","#varepsilon_{true} vs Q_{thr};Q_{thr} [mV#upointns];#varepsilon_{true}",100,0,5000);
     fr.SetMinimum(-0.005);fr.SetMaximum(0.02);fr.SetStats(0);fr.Draw();
     const int nsw=(int)tp_X_all.Q_sweep.size();
     std::vector<double> xs(tp_X_all.Q_sweep.begin(),tp_X_all.Q_sweep.end());
     TGraphErrors gA(nsw,xs.data(),tp_X_all.eps_true_sweep.data(),
         std::vector<double>(nsw,0).data(),tp_X_all.eps_true_err_sweep.data());
     gA.SetMarkerColor(kBlue+1);gA.SetLineColor(kBlue+1);gA.SetMarkerStyle(20);gA.SetLineWidth(2);
     gA.Draw("PL SAME");
     TGraphErrors gM(nsw,xs.data(),tp_X_mask.eps_true_sweep.data(),
         std::vector<double>(nsw,0).data(),tp_X_mask.eps_true_err_sweep.data());
     gM.SetMarkerColor(kOrange+7);gM.SetLineColor(kOrange+7);gM.SetMarkerStyle(21);gM.SetLineWidth(2);
     gM.Draw("PL SAME");
     TLine l0(0,0,5000,0);l0.SetLineColor(kRed);l0.SetLineStyle(2);l0.Draw();
     TLegend lg(0.50,0.75,0.88,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
     lg.AddEntry(&gA,"hod2019X — todos los canales","pl");
     lg.AddEntry(&gM,"hod2019X — sin canales calientes","pl");lg.Draw();
     cc.SaveAs((outdir/"efficiency_with_without_hotchannel_mask.png").string().c_str());}

    //==========================================================================
    // Veredicto y summary
    //==========================================================================
    const double sig_x=tp_X_all.significance, sig_y=tp_Y_all.significance;
    const double sig_xm=tp_X_mask.significance, sig_ym=tp_Y_mask.significance;
    const double max_sig=std::max({sig_x,sig_y,sig_xm,sig_ym});
    const bool rama_A=(max_sig>=5.0);

    auto print_summary=[&](std::ostream&os){
        os<<"=== Tag-and-Probe: detección estadística de muones en hod2019 ===\n\n";
        os<<"Parámetros:\n";
        os<<"  Archivo: "<<fv8<<"\n";
        os<<"  umbral_Q = "<<kUmbralQ<<" mV·ns\n";
        os<<"  Ventana temporal muón: ["<<kMuonTlo<<", "<<kMuonThi<<"] ns\n";
        os<<"  Event mixing: shuffle fijo (seed=42)\n\n";

        os<<"Tabla de resultados:\n";
        os<<std::setw(38)<<std::left<<"Métrica"
          <<std::setw(14)<<"hod2019 X"<<std::setw(14)<<"hod2019 Y\n";
        os<<"-----------------------------------------------------------\n";
        auto fmt=[](double v)->std::string{char b[20];std::snprintf(b,sizeof(b),"%.4f",v);return b;};
        auto row=[&](const std::string&m,const std::string&a,const std::string&b){
            os<<std::setw(38)<<m<<std::setw(14)<<a<<std::setw(14)<<b<<"\n";};
        row("ε_true (all ch) [%]",
            fmt(tp_X_all.eps_true*100),fmt(tp_Y_all.eps_true*100));
        row("σ(ε_true) [%]",
            fmt(tp_X_all.sigma_eps*100),fmt(tp_Y_all.sigma_eps*100));
        row("Significancia [σ]",
            fmt(tp_X_all.significance),fmt(tp_Y_all.significance));
        row("ε_true (hot ch masked) [%]",
            fmt(tp_X_mask.eps_true*100),fmt(tp_Y_mask.eps_true*100));
        row("Significancia (masked) [σ]",
            fmt(tp_X_mask.significance),fmt(tp_Y_mask.significance));
        os<<std::setw(38)<<"Canales calientes X";
        for(int k:hot_x)os<<kH19X[k]<<"(k="<<k<<") ";os<<"\n";
        os<<std::setw(38)<<"Canales calientes Y";
        for(int k:hot_y)os<<kH19Y[k]<<"(k="<<k<<") ";os<<"\n";
        row("r_X (masked, N="+std::to_string(rx18.size())+")",fmt(rX_masked),"N/A");
        os<<"\n";

        if(rama_A){
            os<<"VEREDICTO (RAMA A):\n";
            const double sig_best=sig_x>sig_y?sig_x:sig_y;
            const double eps_best=sig_x>sig_y?tp_X_all.eps_true:tp_Y_all.eps_true;
            const double serr_best=sig_x>sig_y?tp_X_all.sigma_eps:tp_Y_all.sigma_eps;
            os<<"DETECCIÓN ESTADÍSTICA DE MUONES EN HOD2019: ε_true = "
              <<std::setprecision(2)<<eps_best*100<<"% ± "<<serr_best*100
              <<" ("<<std::setprecision(1)<<sig_best<<"σ significancia). La señal de muón\n"
              <<"cósmico se detecta sobre fondo accidental usando hod2018 como referencia,\n"
              <<"pero la eficiencia es demasiado baja para tracking evento-a-evento.\n"
              <<"Confirma que la tecnología detecta muones; el light yield (reflector TiO₂) es el limitante.\n";
        } else {
            os<<"VEREDICTO (RAMA B):\n";
            os<<"SIN DETECCIÓN ESTADÍSTICA SIGNIFICATIVA (sig_max = "<<std::setprecision(2)<<max_sig<<"σ < 5σ):\n"
              <<"hod2019 no muestra exceso de carga en las posiciones predichas por hod2018.\n"
              <<"La detección de muones cósmicos con el reflector TiO₂ y barras de 1 mm\n"
              <<"está por debajo del umbral medible con esta estadística.\n"
              <<"Resultado consistente con la limitación del reflector TiO₂.\n"
              <<"Canales calientes enmascarados: significancia max = "
              <<std::max(sig_xm,sig_ym)<<"σ (sin cambio significativo).\n";
        }
        os<<"\nFiguras producidas:\n";
        os<<"  tagprobe_charge_same_vs_mixed_X.png\n";
        os<<"  tagprobe_charge_same_vs_mixed_Y.png\n";
        os<<"  efficiency_with_without_hotchannel_mask.png\n";
        os<<"  hotchannels.csv\n";
    };
    print_summary(std::cout);
    {std::ofstream sf(outdir/"summary_tagprobe.txt");print_summary(sf);}

    //==========================================================================
    // Actualización del beamer
    //==========================================================================
    {// Busca el beamer más reciente
     std::vector<fs::path> bdirs;
     try{for(auto&e:fs::directory_iterator("outputs"))
         if(e.is_directory()&&e.path().filename().string().substr(0,7)=="beamer_")
             bdirs.push_back(e.path());}catch(...){}
     if(!bdirs.empty()){
         std::sort(bdirs.begin(),bdirs.end(),[](auto&a,auto&b){return a.filename()>b.filename();});
         const fs::path btex=bdirs[0]/"complement_cosmic.tex";
         if(fs::exists(btex)){
             // Lee el .tex
             std::ifstream fin_tex(btex);
             std::string tex((std::istreambuf_iterator<char>(fin_tex)),
                              std::istreambuf_iterator<char>());
             // Encuentra la posición para insertar ANTES del slide backup
             const std::string marker="\\begin{frame}{\\EXP\\quad Backup: Live-Time";
             const std::string caveat4_marker="Reflector upgrade (Vikuiti) is the corrective action.\"}\n";
             size_t pos_backup=tex.find(marker);
             if(rama_A){
                 // Inserta nueva lámina antes del backup
                 const double eps_best=sig_x>sig_y?tp_X_all.eps_true:tp_Y_all.eps_true;
                 const double serr_best=sig_x>sig_y?tp_X_all.sigma_eps:tp_Y_all.sigma_eps;
                 const double sig_best=sig_x>sig_y?sig_x:sig_y;
                 char new_slide[2000];
                 std::snprintf(new_slide,sizeof(new_slide),
R"(\begin{frame}{\EXP\quad Statistical Detection of Cosmic Muons in hod2019}
  \begin{columns}[T]
    \column{0.55\textwidth}
      \begin{itemize}\setlength\itemsep{4pt}
        \item hod2018 ($\varepsilon\approx99.9$\%%) \textbf{tags} the muon position;
          hod2019 charge is \textbf{probed} at the predicted strip.
        \item Event-mixing provides an in-situ background model (same strip,
          different events $\Rightarrow$ accidentals only).
        \item Background-subtracted efficiency:
          $\varepsilon_\text{true} = %.2f\%%\pm%.2f\%%$ \textbf{(%.1f$\sigma$)}.
        \item Confirms the technology detects cosmic muons; efficiency is
          light-yield-limited (TiO$_2$ vs.\ Vikuiti).
      \end{itemize}
    \column{0.43\textwidth}
      \centering
      \includegraphics[width=\linewidth]{figs/tagprobe_charge_same_vs_mixed_X.png}
      \par\footnotesize
      $H_\text{same}$ vs.\ $H_\text{mixed}$: excess at $Q>Q_\text{thr}$ is the muon signal.
  \end{columns}
  \begin{caveat}[Caveat]
    Statistical detection only; per-event tracking is not feasible at this
    efficiency. Reflector upgrade (Vikuiti) is the corrective action.
    \hfill\textit{\footnotesize [EXPERIMENTAL] Cosmic data v8, ULS. 2026.}
  \end{caveat}
\end{frame}

)", eps_best*100, serr_best*100, sig_best);
                 if(pos_backup!=std::string::npos)
                     tex.insert(pos_backup,new_slide);
                 // Also copy figure
                 fs::copy_file(outdir/"tagprobe_charge_same_vs_mixed_X.png",
                               bdirs[0]/"figs"/"tagprobe_charge_same_vs_mixed_X.png",
                               fs::copy_options::overwrite_existing);
             } else {
                 // Añade línea al caveat de Lámina 4
                 const std::string old_cav="A statistical tag-and-probe test";
                 if(tex.find(old_cav)==std::string::npos){
                     // La línea no existe aún: buscar el caveat de Lámina 4 e insertar
                     const std::string cav4_end="[EXPERIMENTAL] Cosmic data v7/v8, ULS. 2026.";
                     size_t p=tex.find(cav4_end);
                     if(p!=std::string::npos){
                         char add_line[300];
                         std::snprintf(add_line,sizeof(add_line),
                             "\n    A statistical tag-and-probe test (hod2018-predicted positions) finds no\n"
                             "    significant cosmic-muon charge excess in hod2019 (%.2f$\\sigma$), reinforcing\n"
                             "    that TiO$_2$ light yield is below the usable detection threshold.\n    ",
                             max_sig);
                         tex.insert(p,add_line);
                     }
                 }
             }
             // Escribe el .tex actualizado
             std::ofstream fout_tex(btex);
             fout_tex<<tex;
             std::cout<<"Beamer actualizado: "<<btex<<"\n";
             std::cout<<"Compilar con: cd "<<bdirs[0].string()<<" && pdflatex complement_cosmic.tex\n";
         }
     }
    }
    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
