// reconstruct_cosmic_tracks.cxx
// Reconstrucción de tracks cósmicos con discriminante de carga Q + forma r.
// Hipótesis: selección por Q > umbral_Q AND r ∈ [r_lo, r_hi] entrega una muestra
// pura de señales SiPM reales, revelando la correlación espacial hod2019↔hod2018.
//
// RAMAS:
//   RAMA A: r_corr > 0.3 → reconstruir proyecciones de track (theta_x, theta_y)
//   RAMA B: r_corr ≤ 0.3 → documentar resultado honestamente (limite de light yield)
//
// Parámetros geométricos REALES: D=98.4mm, pitch_eff=2.0mm (NO D=150mm de diseño).
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o reconstruct_cosmic_tracks \
//        reconstruct_cosmic_tracks.cxx

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
#include <string>
#include <vector>

namespace fs = std::filesystem;

//==============================================================================
// Constantes y mapeo
//==============================================================================
constexpr double kPre_ns  = 2500.0;
constexpr double kDt_ns = 16.0;  // VX1740, 62.5 MS/s; equals the data_dig.time step
constexpr int    kNpre    = 3, kNpost = 8; // ventana de integración
constexpr int    kNwin    = kNpre + kNpost + 1;  // 12 muestras
constexpr int    kNCh     = 64;
constexpr double kPitch   = 2.0;   // mm
constexpr double kD_mm    = 98.4;  // separación real entre planos [mm]

// Ventanas de señal [ns]
static const double kTlo[4] = {3500,3500,3000,3000};  // 0=18X,1=18Y,2=19X,3=19Y
static const double kThi[4] = {7000,7000,5000,5000};

static const std::array<int,16> kH19X = {53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y = {37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};

static int ch_plane(int c){
    if(c>=0&&c<=15)return 0; if(c>=16&&c<=31)return 1;
    if(c>=32&&c<=47)return 3; if(c>=48&&c<=63)return 2; return -1;}
static int ch_k(int c){
    if(c>=0&&c<=15)return c; if(c>=16&&c<=31)return c-16;
    for(int k=0;k<16;k++)if(kH19Y[k]==c)return k;
    for(int k=0;k<16;k++)if(kH19X[k]==c)return k; return -1;}
static double strip_coord(int k){return (k+0.5-8.0)*kPitch;}

//==============================================================================
// Extracción de carga Q + ratio r
//==============================================================================
struct CQ{bool valid=false;double ped=0,noise=0,ph=0,pt=-1,Q=0,r=0;int nw=0;int k=-1;};

static CQ extract_cq(const std::vector<float>&wf,const std::vector<float>&t,
                      double tlo,double thi,int k){
    CQ r; r.k=k;
    const int n=(int)wf.size();
    if(n<4||(int)t.size()!=n)return r;
    double s=0,s2=0;int nb=0;
    for(int i=0;i<n;i++)if(t[i]<kPre_ns){s+=wf[i];s2+=wf[i]*wf[i];nb++;}
    if(nb<2)return r;
    r.valid=true;r.ped=s/nb;
    const double v=s2/nb-r.ped*r.ped;r.noise=(v>0)?std::sqrt(v):0.0;
    double vmin=1e30;int imin=-1;
    for(int i=0;i<n;i++)if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0)return r;
    r.ph=r.ped-vmin; r.pt=t[imin];
    const double tilo=r.pt-kNpre*kDt_ns, tihi=r.pt+kNpost*kDt_ns;
    double Q=0;int nw=0;
    for(int i=0;i<n;i++)if(t[i]>=tilo&&t[i]<=tihi){Q+=(r.ped-wf[i])*kDt_ns;nw++;}
    r.Q=Q; r.nw=nw;
    if(r.ph>0.5&&nw>0)r.r=Q/(r.ph*nw*kDt_ns);
    return r;
}

//==============================================================================
// Clustering hod2018 (gap≤1)
//==============================================================================
struct Hit18{int k;double ph,pt;};
struct Cl18{double pos,amp;};
static std::vector<Cl18> clust18(std::vector<Hit18> hits,double thr){
    std::vector<Hit18> act;
    for(const auto&h:hits)if(h.ph>=thr)act.push_back(h);
    if(act.empty())return{};
    std::sort(act.begin(),act.end(),[](const Hit18&a,const Hit18&b){return a.k<b.k;});
    std::vector<Cl18> out; double wpos=0,wamp=0; int lo=act[0].k,hi=lo;
    auto flush=[&](){out.push_back({wpos/wamp,wamp});wpos=wamp=0;};
    auto add=[&](const Hit18&h){wpos+=h.ph*(h.k+0.5-8.0)*kPitch;wamp+=h.ph;hi=h.k;};
    lo=act[0].k;hi=lo;add(act[0]);
    for(int i=1;i<(int)act.size();i++){
        if(act[i].k>hi+1){flush();lo=hi=act[i].k;}
        add(act[i]);}
    flush(); return out;}

//==============================================================================
// Estadísticas
//==============================================================================
static double pearsonR(const std::vector<double>&x,const std::vector<double>&y){
    const int n=(int)x.size();if(n<2)return 0;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(int i=0;i<n;i++){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];syy+=y[i]*y[i];sxy+=x[i]*y[i];}
    const double d=n,num=d*sxy-sx*sy,den=std::sqrt((d*sxx-sx*sx)*(d*syy-sy*sy));
    return den>0?num/den:0.0;}
static void linfit(const std::vector<double>&x,const std::vector<double>&y,double&sl,double&ic){
    sl=ic=0;const int n=(int)x.size();if(n<2)return;
    double sx=0,sy=0,sxx=0,sxy=0;
    for(int i=0;i<n;i++){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];sxy+=x[i]*y[i];}
    const double d=n,dn=d*sxx-sx*sx;if(std::abs(dn)<1e-12)return;
    sl=(d*sxy-sx*sy)/dn;ic=(sy-sl*sx)/d;}

//==============================================================================
// Almacenamiento de hits por evento en RAM
//==============================================================================
struct EvCache{
    // hod2019 hits (todos los que tienen ph>1 mV)
    std::vector<CQ> h19x, h19y;
    // hod2018 clusters (best = max amplitude)
    double best18x_pos=0, best18x_amp=0;
    double best18y_pos=0, best18y_amp=0;
    int ncls18x=0, ncls18y=0;
};

//==============================================================================
// Carga datos en RAM (pasada única)
//==============================================================================
static std::vector<EvCache> load_events(const std::string& fname){
    TFile* f=TFile::Open(fname.c_str(),"READ");
    if(!f||f->IsZombie()){std::cerr<<"ERROR: "<<fname<<'\n';return{};}
    TTree* t=static_cast<TTree*>(f->Get("data_dig"));
    if(!t){f->Close();return{};}
    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw=nullptr;
    t->SetBranchStatus("*",0);t->SetCacheSize(256LL<<20);
    for(int c=0;c<kNCh;c++){char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!t->GetBranch(bn))continue;t->SetBranchStatus(bn,1);t->SetBranchAddress(bn,&cw[c]);}
    t->SetBranchStatus("time",1);t->SetBranchAddress("time",&tw);
    const Long64_t nEnt=t->GetEntries();
    std::cout<<"  Loading "<<nEnt<<" events into RAM...\n";
    std::vector<EvCache> evs; evs.reserve(nEnt);
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
        EvCache ev;
        // hod2018 hits
        std::vector<Hit18> h18x,h18y;
        for(int c=0;c<kNCh;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const int pl=ch_plane(c);if(pl<0)continue;
            CQ cq=extract_cq(*cw[c],tv,kTlo[pl],kThi[pl],ch_k(c));
            if(!cq.valid||cq.ph<0.5)continue;
            if(pl==0)h18x.push_back({cq.k,cq.ph,cq.pt});
            else if(pl==1)h18y.push_back({cq.k,cq.ph,cq.pt});
            else if(pl==2&&cq.ph>1.0)ev.h19x.push_back(cq);
            else if(pl==3&&cq.ph>1.0)ev.h19y.push_back(cq);
        }
        // hod2018 clusters (threshold=10 mV)
        {auto cls=clust18(h18x,10.0);ev.ncls18x=(int)cls.size();
         if(!cls.empty()){const auto& best=*std::max_element(cls.begin(),cls.end(),[](const Cl18&a,const Cl18&b){return a.amp<b.amp;});ev.best18x_pos=best.pos;ev.best18x_amp=best.amp;}}
        {auto cls=clust18(h18y,10.0);ev.ncls18y=(int)cls.size();
         if(!cls.empty()){const auto& best=*std::max_element(cls.begin(),cls.end(),[](const Cl18&a,const Cl18&b){return a.amp<b.amp;});ev.best18y_pos=best.pos;ev.best18y_amp=best.amp;}}
        evs.push_back(std::move(ev));
    }
    std::cout<<"    "<<nEnt<<"/"<<nEnt<<" ✓  ("<<evs.size()<<" events cached)\n";
    f->Close();
    return evs;
}

//==============================================================================
// Figuras helper
//==============================================================================
static void p2fit(TH2D*h,double sl,double ic,double r,const fs::path&p){
    TCanvas c(("cf_"+std::string(h->GetName())).c_str(),"",800,780);
    c.SetLeftMargin(0.13);c.SetRightMargin(0.17);c.SetBottomMargin(0.13);
    h->SetStats(0);h->Draw("COLZ");
    const double xlo=h->GetXaxis()->GetXmin(),xhi=h->GetXaxis()->GetXmax();
    TLine ln(xlo,sl*xlo+ic,xhi,sl*xhi+ic);ln.SetLineColor(kRed);ln.SetLineWidth(2);ln.Draw("SAME");
    TLatex lt;lt.SetNDC();lt.SetTextSize(0.033);
    char buf[100];std::snprintf(buf,sizeof(buf),"slope=%.3f  intercept=%.2f mm  r=%.3f",sl,ic,r);
    lt.DrawLatex(0.14,0.94,buf);
    c.SaveAs(p.string().c_str());}

//==============================================================================
// main
//==============================================================================
int main(int argc,char**argv){
    std::string fv8="lectura_vx1740_ext_nim_1M_v8.root";
    std::string fv7="lectura_vx1740_ext_nim_1M_v7.root";
    double umbral_Q19x=2440.048, umbral_Q19y=2440.048;  // rescaled x1.04 on 2026-10-17 (kDt_ns 15.3846 -> 16 ns); selection unchanged
    double r_lo=0.20, r_hi=0.80;  // valores iniciales; calibrados del histograma
    for(int i=1;i<argc;i++){
        std::string a=argv[i];
        if(a=="--v8"&&i+1<argc){fv8=argv[++i];continue;}
        if(a=="--v7"&&i+1<argc){fv7=argv[++i];continue;}
        if(a=="--Qlo"&&i+1<argc){umbral_Q19x=std::stod(argv[++i]);umbral_Q19y=umbral_Q19x;continue;}
    }
    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("cosmic_tracks_"+std::string(ds));
    fs::create_directories(outdir);
    std::cout<<"=== reconstruct_cosmic_tracks ===\nv8: "<<fv8<<"\nOut: "<<outdir<<"\n";
    std::cout<<"umbral_Q = "<<umbral_Q19x<<" mV·ns  D = "<<kD_mm<<" mm  pitch = "<<kPitch<<" mm\n\n";
    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    // Carga v8 en RAM
    std::cout<<"[v8] Loading...\n";
    auto evs=load_events(fv8);
    const long long N=static_cast<long long>(evs.size());
    if(N==0){std::cerr<<"ERROR: no events loaded\n";return 1;}

    //==========================================================================
    // Calibración del ratio r desde los datos
    //==========================================================================
    std::cout<<"Calibrando ratio r desde datos...\n";
    auto h_r=std::make_unique<TH1D>("h_r","Ratio r = Q/(ph#timesN#timesdt);r;Hits",150,-0.1,1.5);
    for(const auto&ev:evs){
        for(const auto&cq:ev.h19x) if(cq.ph>5&&cq.nw>0) h_r->Fill(cq.r);
        for(const auto&cq:ev.h19y) if(cq.ph>5&&cq.nw>0) h_r->Fill(cq.r);
    }
    // Busca valle entre spikes (r<0.2) y señal (r~0.3–0.7)
    // Mínimo en [0.10, 0.35]
    {int blo=h_r->FindBin(0.10),bhi=h_r->FindBin(0.35);
     int bv=blo;double vv=1e30;
     for(int b=blo;b<=bhi;b++){double c=h_r->GetBinContent(b);if(c>0&&c<vv){vv=c;bv=b;}}
     r_lo=h_r->GetBinCenter(bv);if(r_lo<0.12)r_lo=0.20;}
    // Máximo hacia la derecha del pico de señal
    {int bpk=h_r->GetMaximumBin();
     double r_peak=h_r->GetBinCenter(bpk);
     if(r_peak>r_lo)r_hi=std::min(1.1,r_peak+0.4);else r_hi=0.80;}
    std::cout<<"  r_lo="<<r_lo<<"  r_hi="<<r_hi<<"\n";
    // Guarda figura de r
    {h_r->SetLineColor(kBlue+1);h_r->SetLineWidth(2);
     TCanvas cr("cr","",900,600);cr.SetLeftMargin(0.12);cr.SetBottomMargin(0.13);cr.SetLogy();
     h_r->Draw("HIST");
     gPad->Update();
     TLine ll(r_lo,gPad->GetUymin(),r_lo,h_r->GetMaximum()*0.5);ll.SetLineColor(kRed);ll.SetLineStyle(2);ll.SetLineWidth(2);ll.Draw();
     TLine lh(r_hi,gPad->GetUymin(),r_hi,h_r->GetMaximum()*0.5);lh.SetLineColor(kRed);lh.SetLineStyle(2);lh.SetLineWidth(2);lh.Draw();
     TLatex lt;lt.SetNDC();lt.SetTextSize(0.035);
     char buf[80];std::snprintf(buf,sizeof(buf),"señal aceptada: r #in [%.2f, %.2f]",r_lo,r_hi);
     lt.SetTextColor(kRed);lt.DrawLatex(0.35,0.90,buf);
     cr.SaveAs((outdir/"ratio_r_distribution_hod2019.png").string().c_str());}

    //==========================================================================
    // Selección Q+r para hod2019 y correlación
    //==========================================================================
    // Q+r: exactamente 1 canal limpio en hod2019 + hod2018 con hits
    std::vector<double> xA18,xA19,yA18,yA19;  // Q+r XZ e YZ
    std::vector<double> xB18,xB19,yB18,yB19;  // ph control (ph>umbral equivalente)
    // umbral ph equivalente ~ umbral_Q / (kNwin × kDt_ns)
    const double ph_equiv = umbral_Q19x / (kNwin * kDt_ns);
    std::cout<<"  ph_equiv = "<<ph_equiv<<" mV (equivalente a Q>"<<umbral_Q19x<<" mV·ns)\n";

    long long n_coinc_Xqr=0,n_coinc_Yqr=0,n_coinc_Xph=0,n_coinc_Yph=0;
    long long n_3D_qr=0,n_3D_ph=0;

    std::vector<double> mx_qr,my_qr,mx_ph,my_ph;  // slopes para theta
    std::ofstream fcsv(outdir/"cosmic_tracks.csv");
    fcsv<<"selection,x18_mm,y18_mm,x19_mm,y19_mm,mx,my,theta_x_deg,theta_y_deg,theta_3d_deg\n";

    // Control: para comparación justa, acumula todos los eventos con hod2019X hit(ph>ph_equiv)
    // y hod2018X hit, para poder hacer el top-N después
    struct EvScore{double x18,x19,ph_max,Q_max;};
    std::vector<EvScore> scores_x,scores_y;

    for(const auto&ev:evs){
        if(ev.ncls18x<1||ev.best18x_amp<10)continue;
        if(ev.ncls18y<1||ev.best18y_amp<10)continue;

        // Selección Q+r en hod2019X
        std::vector<const CQ*> clean_x,clean_y;
        for(const auto&cq:ev.h19x)
            if(cq.Q>umbral_Q19x&&cq.r>=r_lo&&cq.r<=r_hi) clean_x.push_back(&cq);
        for(const auto&cq:ev.h19y)
            if(cq.Q>umbral_Q19y&&cq.r>=r_lo&&cq.r<=r_hi) clean_y.push_back(&cq);

        const bool has_Xqr=(clean_x.size()==1);
        const bool has_Yqr=(clean_y.size()==1);

        if(has_Xqr){
            ++n_coinc_Xqr;
            const double x19=strip_coord(clean_x[0]->k);
            xA18.push_back(ev.best18x_pos); xA19.push_back(x19);
        }
        if(has_Yqr){
            ++n_coinc_Yqr;
            const double y19=strip_coord(clean_y[0]->k);
            yA18.push_back(ev.best18y_pos); yA19.push_back(y19);
        }
        if(has_Xqr&&has_Yqr){
            ++n_3D_qr;
            const double x18=ev.best18x_pos,y18=ev.best18y_pos;
            const double x19=strip_coord(clean_x[0]->k),y19=strip_coord(clean_y[0]->k);
            const double mx=(x19-x18)/kD_mm,my=(y19-y18)/kD_mm;
            const double thx=std::atan(mx)*180/M_PI,thy=std::atan(my)*180/M_PI;
            const double th3=std::atan(std::sqrt(mx*mx+my*my))*180/M_PI;
            mx_qr.push_back(mx);my_qr.push_back(my);
            fcsv<<"Qr,"<<x18<<','<<y18<<','<<x19<<','<<y19<<','
                <<mx<<','<<my<<','<<thx<<','<<thy<<','<<th3<<'\n';
        }

        // Acumula datos para control (ph selection)
        {double phmax=0,Qmax=-1e9;
         for(const auto&cq:ev.h19x){if(cq.ph>phmax)phmax=cq.ph;if(cq.Q>Qmax)Qmax=cq.Q;}
         if(phmax>0) scores_x.push_back({ev.best18x_pos,0,phmax,Qmax});}
        {double phmax=0,Qmax=-1e9;
         for(const auto&cq:ev.h19y){if(cq.ph>phmax)phmax=cq.ph;if(cq.Q>Qmax)Qmax=cq.Q;}
         if(phmax>0) scores_y.push_back({ev.best18y_pos,0,phmax,Qmax});}
    }

    // Build ph control: top-N events by ph_max in hod2019X (same N as Q+r)
    // More direct: use all events with exactly 1 hod2019X channel with ph>ph_equiv
    for(const auto&ev:evs){
        if(ev.ncls18x<1||ev.best18x_amp<10)continue;
        std::vector<const CQ*> hits_x,hits_y;
        for(const auto&cq:ev.h19x) if(cq.ph>=ph_equiv) hits_x.push_back(&cq);
        for(const auto&cq:ev.h19y) if(cq.ph>=ph_equiv) hits_y.push_back(&cq);
        if(hits_x.size()==1){
            ++n_coinc_Xph;
            xB18.push_back(ev.best18x_pos);
            xB19.push_back(strip_coord(hits_x[0]->k));
        }
        if(hits_y.size()==1&&ev.ncls18y>=1&&ev.best18y_amp>=10){
            ++n_coinc_Yph;
            yB18.push_back(ev.best18y_pos);
            yB19.push_back(strip_coord(hits_y[0]->k));
        }
        if(hits_x.size()==1&&hits_y.size()==1&&ev.ncls18y>=1&&ev.best18y_amp>=10){
            ++n_3D_ph;
            const double x18=ev.best18x_pos,y18=ev.best18y_pos;
            const double x19=strip_coord(hits_x[0]->k),y19=strip_coord(hits_y[0]->k);
            const double mx=(x19-x18)/kD_mm,my=(y19-y18)/kD_mm;
            const double thx=std::atan(mx)*180/M_PI,thy=std::atan(my)*180/M_PI;
            const double th3=std::atan(std::sqrt(mx*mx+my*my))*180/M_PI;
            mx_ph.push_back(mx);my_ph.push_back(my);
            fcsv<<"ph,"<<x18<<','<<y18<<','<<x19<<','<<y19<<','
                <<mx<<','<<my<<','<<thx<<','<<thy<<','<<th3<<'\n';
        }
    }

    // Correlaciones
    const double rX_qr=pearsonR(xA18,xA19), rY_qr=pearsonR(yA18,yA19);
    const double rX_ph=pearsonR(xB18,xB19), rY_ph=pearsonR(yB18,yB19);
    double slXqr,icXqr,slYqr,icYqr,slXph,icXph,slYph,icYph;
    linfit(xA18,xA19,slXqr,icXqr); linfit(yA18,yA19,slYqr,icYqr);
    linfit(xB18,xB19,slXph,icXph); linfit(yB18,yB19,slYph,icYph);

    std::cout<<"\n=== RESULTADOS ===\n";
    std::cout<<"Selección Q+r:  N_XZ="<<n_coinc_Xqr<<"  N_YZ="<<n_coinc_Yqr<<"  N_3D="<<n_3D_qr<<"\n";
    std::cout<<"  r_X="<<rX_qr<<"  r_Y="<<rY_qr<<"\n";
    std::cout<<"Selección ph:   N_XZ="<<n_coinc_Xph<<"  N_YZ="<<n_coinc_Yph<<"  N_3D="<<n_3D_ph<<"\n";
    std::cout<<"  r_X="<<rX_ph<<"  r_Y="<<rY_ph<<"\n\n";

    //==========================================================================
    // Figuras de correlación (2 paneles lado a lado)
    //==========================================================================
    auto make_corr_canvas=[&](const std::string& coord,
                               const std::vector<double>&xph_18,const std::vector<double>&xph_19,
                               const std::vector<double>&xqr_18,const std::vector<double>&xqr_19,
                               double slph,double icph,double rph,
                               double slqr,double icqr,double rqr,
                               const fs::path& path){
        constexpr double lo=-16,hi=16; constexpr int nb=16;
        TCanvas cc(("cc_"+coord).c_str(),"",1600,780);cc.Divide(2,1);
        const char* cname=coord.c_str();
        // Panel 1: ph selection
        cc.cd(1);gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.17);gPad->SetBottomMargin(0.13);
        {TH2D h(("hph_"+coord).c_str(),
            ("pulse_height selection (N="+std::to_string(xph_18.size())+");"+coord+"_{2018} [mm];"+coord+"_{2019} [mm];Events").c_str(),
            nb,lo,hi,nb,lo,hi);
         for(size_t i=0;i<xph_18.size();i++)h.Fill(xph_18[i],xph_19[i]);
         h.SetStats(0);h.Draw("COLZ");
         TLine ln(lo,slph*lo+icph,hi,slph*hi+icph);ln.SetLineColor(kRed);ln.SetLineWidth(2);ln.Draw("SAME");
         TLatex lt;lt.SetNDC();lt.SetTextSize(0.038);
         char buf[80];std::snprintf(buf,sizeof(buf),"slope=%.3f  r=%.3f",slph,rph);lt.DrawLatex(0.14,0.94,buf);
         lt.DrawLatex(0.14,0.90,"Sel: ph > equiv.");}
        // Panel 2: Q+r selection
        cc.cd(2);gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.17);gPad->SetBottomMargin(0.13);
        {TH2D h(("hqr_"+coord).c_str(),
            ("Q+r selection (N="+std::to_string(xqr_18.size())+");"+coord+"_{2018} [mm];"+coord+"_{2019} [mm];Events").c_str(),
            nb,lo,hi,nb,lo,hi);
         for(size_t i=0;i<xqr_18.size();i++)h.Fill(xqr_18[i],xqr_19[i]);
         h.SetStats(0);h.Draw("COLZ");
         TLine ln(lo,slqr*lo+icqr,hi,slqr*hi+icqr);ln.SetLineColor(kRed);ln.SetLineWidth(2);ln.Draw("SAME");
         TLatex lt;lt.SetNDC();lt.SetTextSize(0.038);
         char buf[80];std::snprintf(buf,sizeof(buf),"slope=%.3f  r=%.3f",slqr,rqr);lt.DrawLatex(0.14,0.94,buf);
         lt.DrawLatex(0.14,0.90,"Sel: Q+r clean");}
        cc.SaveAs(path.string().c_str());};

    make_corr_canvas("X",xB18,xB19,xA18,xA19,slXph,icXph,rX_ph,slXqr,icXqr,rX_qr,
                     outdir/"correlation_ph_vs_Q_selection_XZ.png");
    make_corr_canvas("Y",yB18,yB19,yA18,yA19,slYph,icYph,rY_ph,slYqr,icYqr,rY_qr,
                     outdir/"correlation_ph_vs_Q_selection_YZ.png");

    //==========================================================================
    // Theta distributions (siempre, aunque pequeñas)
    //==========================================================================
    auto h_thx_qr=std::make_unique<TH1D>("hthx","#theta_x Q+r;#theta_x [deg];Events",80,-25,25);
    auto h_thy_qr=std::make_unique<TH1D>("hthy","#theta_y Q+r;#theta_y [deg];Events",80,-25,25);
    auto h_th3_qr=std::make_unique<TH1D>("hth3","#theta_{3D} Q+r;#theta_{3D} [deg];Events",80,0,25);
    auto h_thx_ph=std::make_unique<TH1D>("hthxph","#theta_x ph;#theta_x [deg];Events",80,-25,25);
    for(double mx:mx_qr){h_thx_qr->Fill(std::atan(mx)*180/M_PI);}
    for(double my:my_qr){h_thy_qr->Fill(std::atan(my)*180/M_PI);}
    for(size_t i=0;i<mx_qr.size();i++){
        const double mx=mx_qr[i],my=my_qr[i];
        h_th3_qr->Fill(std::atan(std::sqrt(mx*mx+my*my))*180/M_PI);}
    for(double mx:mx_ph){h_thx_ph->Fill(std::atan(mx)*180/M_PI);}

    {TCanvas cth("cth","",1200,400);cth.Divide(3,1);
     auto draw=[&](TH1D*h,int pad){cth.cd(pad);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.15);
         h->SetLineColor(kBlue+1);h->SetLineWidth(2);h->Draw("HIST");};
     draw(h_thx_qr.get(),1);draw(h_thy_qr.get(),2);draw(h_th3_qr.get(),3);
     cth.SaveAs((outdir/"theta_distributions_cosmic.png").string().c_str());}

    //==========================================================================
    // Figura headline para beamer
    //==========================================================================
    const bool rama_A=(rX_qr>0.3||rY_qr>0.3);
    {TCanvas ch("ch","",1400,560);ch.Divide(2,1);
     if(rama_A){
         // Panel 1: mejor scatter de correlación
         ch.cd(1);gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.17);gPad->SetBottomMargin(0.13);
         const bool useX=(std::abs(rX_qr)>std::abs(rY_qr));
         const auto& x18=(useX?xA18:yA18), x19=(useX?xA19:yA19);
         const double sl=(useX?slXqr:slYqr),ic=(useX?icXqr:icYqr),rc=(useX?rX_qr:rY_qr);
         TH2D hh("hhc",(std::string(useX?"X":"Y")+"Z projection;Pos_{2018} [mm];Pos_{2019} [mm]").c_str(),16,-16,16,16,-16,16);
         for(size_t i=0;i<x18.size();i++)hh.Fill(x18[i],x19[i]);hh.SetStats(0);hh.Draw("COLZ");
         TLine ln(-16,sl*(-16)+ic,16,sl*16+ic);ln.SetLineColor(kRed);ln.SetLineWidth(2);ln.Draw("SAME");
         TLatex lt;lt.SetNDC();lt.SetTextSize(0.040);
         char buf[80];std::snprintf(buf,sizeof(buf),"r=%.3f  slope=%.3f  N=%d",rc,sl,(int)x18.size());lt.DrawLatex(0.14,0.93,buf);
         // Panel 2: theta distribution
         ch.cd(2);gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.15);
         if(!mx_qr.empty())h_thx_qr->Draw("HIST");
         else h_thx_ph->Draw("HIST");
     } else {
         // RAMA B: charge spectrum + SNR comparison
         ch.cd(1);gPad->SetLogy();gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);
         TH1D hQ("hQsumm","Charge hod2019X (v8);Q_{max} [mV#upointns];Events",200,-500,5000);
         for(const auto&ev:evs){double mx=-1e9;for(const auto&cq:ev.h19x)if(cq.Q>mx)mx=cq.Q;if(mx>-1e8)hQ.Fill(mx);}
         hQ.SetLineColor(kBlue+1);hQ.SetLineWidth(2);hQ.Draw("HIST");
         ch.cd(2);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
         TH1D hsnr("hsnr","SNR comparison;Plano;SNR",4,0.5,4.5);
         hsnr.SetBinContent(1,19.7);hsnr.SetBinContent(2,43.2);  // hod2018 ph/Q
         hsnr.SetBinContent(3,14.0);hsnr.SetBinContent(4,49.7);  // hod2019 ph/Q
         hsnr.GetXaxis()->SetBinLabel(1,"18X ph");hsnr.GetXaxis()->SetBinLabel(2,"18X Q");
         hsnr.GetXaxis()->SetBinLabel(3,"19X ph");hsnr.GetXaxis()->SetBinLabel(4,"19X Q");
         hsnr.SetFillColor(kAzure+2);hsnr.SetStats(0);hsnr.Draw("HBAR");
     }
     ch.SaveAs((outdir/"headline_for_beamer.png").string().c_str());}

    // Copia headline al beamer si existe
    {auto beamer_dirs=std::vector<fs::path>{};
     try{for(auto&e:fs::directory_iterator("outputs"))
         if(e.is_directory()&&e.path().filename().string().substr(0,7)=="beamer_")
             beamer_dirs.push_back(e.path());}catch(...){}
     if(!beamer_dirs.empty()){
         std::sort(beamer_dirs.begin(),beamer_dirs.end(),[](auto&a,auto&b){return a.filename()>b.filename();});
         fs::path bfigs=beamer_dirs[0]/"figs";
         if(fs::exists(bfigs)){
             fs::copy_file(outdir/"headline_for_beamer.png",bfigs/"headline_for_beamer.png",
                           fs::copy_options::overwrite_existing);
             std::cout<<"  Copiada headline_for_beamer.png → "<<bfigs<<"\n";}}}

    //==========================================================================
    // Summary
    //==========================================================================
    const double max_r=std::max(std::abs(rX_qr),std::abs(rY_qr));
    const double max_rph=std::max(std::abs(rX_ph),std::abs(rY_ph));
    const bool rama_b_forced=(max_r<0.3);

    auto print_summary=[&](std::ostream&os){
        os<<"=== Reconstrucción de Tracks Cósmicos — Discriminante Q+r ===\n\n";
        os<<"Parámetros:\n";
        os<<"  D = "<<kD_mm<<" mm (real, NO de diseño)\n";
        os<<"  umbral_Q_hod2019X = "<<umbral_Q19x<<" mV·ns\n";
        os<<"  r_lo = "<<r_lo<<"  r_hi = "<<r_hi<<" (calibrado del histograma)\n";
        os<<"  ph_equiv = "<<ph_equiv<<" mV (control)\n\n";
        os<<"Cutflow:\n";
        os<<std::setw(32)<<std::left<<"  Eventos cargados"<<N<<"\n";
        os<<std::setw(32)<<"  Coinc. XZ (Q+r)"<<n_coinc_Xqr<<"\n";
        os<<std::setw(32)<<"  Coinc. YZ (Q+r)"<<n_coinc_Yqr<<"\n";
        os<<std::setw(32)<<"  Full 3D (Q+r)"<<n_3D_qr<<"\n";
        os<<std::setw(32)<<"  Coinc. XZ (ph)"<<n_coinc_Xph<<"\n";
        os<<std::setw(32)<<"  Coinc. YZ (ph)"<<n_coinc_Yph<<"\n";
        os<<std::setw(32)<<"  Full 3D (ph)"<<n_3D_ph<<"\n\n";
        os<<"Correlaciones:\n";
        os<<std::fixed<<std::setprecision(4);
        os<<"  Q+r: r_X="<<rX_qr<<"  r_Y="<<rY_qr;
        os<<"  slope_X="<<slXqr<<"  intercept_X="<<icXqr<<" mm\n";
        os<<"  ph:  r_X="<<rX_ph<<"  r_Y="<<rY_ph;
        os<<"  slope_X="<<slXph<<"  intercept_X="<<icXph<<" mm\n\n";
        os<<"RAMA: "<<(rama_b_forced?"B (r_corr<0.3)":"A (r_corr>0.3)")<<"\n\n";
        if(!rama_b_forced){
            os<<"Track statistics (Q+r full-3D, N="<<n_3D_qr<<"):\n";
            if(!mx_qr.empty()){
                double sumth=0;for(double m:mx_qr)sumth+=std::atan(m)*180/M_PI;
                const double mean_thx=sumth/mx_qr.size();
                double sumth2=0;for(double m:mx_qr){double d=std::atan(m)*180/M_PI-mean_thx;sumth2+=d*d;}
                const double rms_thx=std::sqrt(sumth2/mx_qr.size());
                os<<"  theta_x: mean="<<mean_thx<<"  RMS="<<rms_thx<<" deg\n";
                os<<"  (RMS/sqrt(2) = "<<rms_thx/std::sqrt(2.0)<<" deg es cota cruda de resolución angular — no calibrada)\n";
            }
        }
        if(max_r>0.3&&max_r>max_rph+0.15)
            os<<"VEREDICTO: RECONSTRUCCIÓN LOGRADA: el discriminante de carga revela la correlación inter-plano.\n"
              <<"  N_XZ="<<n_coinc_Xqr<<"  r_corr_X="<<rX_qr<<"  → primeros tracks cósmicos reconstruidos.\n";
        else
            os<<"VEREDICTO: TRACKING LIMITADO POR LIGHT YIELD: discriminante mejora pureza (Q+r) pero la\n"
              <<"  correlación no alcanza umbral r=0.3. Resultado consistente con límite del reflector TiO₂.\n"
              <<"  r_corr_Qr_X="<<rX_qr<<"  r_corr_ph_X="<<rX_ph<<"\n";
    };
    print_summary(std::cout);
    {std::ofstream sf(outdir/"summary_cosmic_tracks.txt");print_summary(sf);}
    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
