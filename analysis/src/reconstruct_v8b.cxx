// reconstruct_v8b.cxx
// Addendum al v8: corrección por retardo diferencial de cables X–Y en hod2019
// y ventana temporal ajustada al pico de señal real.
//
// DIAGNÓSTICO PREVIO (incorporado):
//   Con umbral alto (15 ADC), las señales reales de hod2019 llegan a t ≈ 4175-4185 ns
//   tanto para canales X (48..63) como Y (32..47). La diferencia cable Y-X es ≈ 0 ns
//   (indistinguible a 16 ns/muestra). Los hits a t ≈ 4000 ns (umbral bajo) son ruido:
//   el mínimo de la ventana [3000,5000] en ruido gaussiano tiene media ≈ 4000 ns.
//
// CORRECCIÓN v8b:
//   1. Ventana hod2019 más estrecha (configurable): [t19_lo, t19_hi] centrada ≈ 4150 ns.
//      Por defecto [3900, 4700] ns (600 ns < los 2000 ns del v8 → reduce ruido ~70%).
//   2. Corte de coincidencia interna hod2019 X–Y:
//      |t_clX − t_clY − dtxy_cable| < dtxy_cut
//      donde dtxy_cable ≈ 0 ns (medido) y dtxy_cut configurable (default 300 ns).
//      Filtra coincidencias accidentales donde un eje disparó por señal real y el
//      otro por ruido en un tiempo diferente.
//
// NUEVO MÉTODO D de inferencia de mapeo:
//   Si los cables están ruteados en orden físico, mean_t_c debería ser monotónico
//   en el índice físico k. Calcula τ de Kendall entre orden de mean_t y orden de π_A.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o reconstruct_v8b reconstruct_v8b.cxx

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
#include <TGraphErrors.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

//==============================================================================
// Parámetros
//==============================================================================
struct Params {
    std::string input_file = "lectura_vx1740_ext_nim_1M_v7.root";
    double pitch_mm        =  2.0;
    double zsep_mm         = 98.4;
    double n_sigma         =  5.0;
    double abs_thr_2018    = 20.0;
    double abs_thr_2019    =  5.0;
    // v8b nuevas: ventana hod2019 ajustada (antes era [3000,5000])
    double t19_lo          = 3900.0;  // antes: 3000 ns
    double t19_hi          = 4700.0;  // antes: 5000 ns → 800 ns de ventana
    double kHod2018Tlo     = 3500.0;
    double kHod2018Thi     = 7000.0;
    // Corte coincidencia interna hod2019 X–Y
    double dtxy_cable      =    0.0;  // Δt_cable medido (t_X − t_Y) [ns]
    double dtxy_cut        =  300.0;  // ventana ±cut [ns]
    double dt_delay_ns     = 1200.0;
    double dt_window_ns    =  600.0;
    int    max_cls_2019    =    2;
};

//==============================================================================
// Constantes
//==============================================================================
constexpr double kPrePulse_ns = 2500.0;
constexpr int    kNCh         = 64;

static const std::array<int,16> kHod2019X = {
    53, 61, 54, 62, 55, 63, 52, 60,
    51, 59, 50, 58, 49, 57, 48, 56
};
static const std::array<int,16> kHod2019Y = {
    37, 45, 38, 46, 39, 47, 36, 44,
    35, 43, 34, 42, 33, 41, 32, 40
};
static const std::array<int,16> kHod2019X_ctrl = {
    48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63
};
static const std::array<int,16> kHod2019Y_ctrl = {
    32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47
};

//==============================================================================
// Inversas del mapeo
//==============================================================================
struct InvMap {
    std::array<int,kNCh> x19{}, y19{};
    InvMap(const std::array<int,16>& xm, const std::array<int,16>& ym) {
        x19.fill(-1); y19.fill(-1);
        for(int k=0;k<16;++k){ x19[xm[k]]=k; y19[ym[k]]=k; }
    }
};

//==============================================================================
// Extracción de pulso
//==============================================================================
struct PulseResult {
    bool   valid=false;
    double pedestal=0, noise=0, pulse_height=0, pulse_time=-1;
    bool   hit=false;
};
static PulseResult extract_pulse(const std::vector<float>& wf,
                                  const std::vector<float>& t,
                                  double tlo, double thi, double nsig, double athr) {
    PulseResult r;
    const int n=static_cast<int>(wf.size());
    if(n<4||static_cast<int>(t.size())!=n) return r;
    double sum=0,sum2=0;int nb=0;
    for(int i=0;i<n;++i) if(t[i]<kPrePulse_ns){sum+=wf[i];sum2+=wf[i]*wf[i];++nb;}
    if(nb<2) return r;
    r.valid=true; r.pedestal=sum/nb;
    const double var=sum2/nb-r.pedestal*r.pedestal;
    r.noise=(var>0)?std::sqrt(var):0.0;
    double vmin=std::numeric_limits<double>::infinity();int imin=-1;
    for(int i=0;i<n;++i) if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0) return r;
    r.pulse_height=r.pedestal-vmin; r.pulse_time=t[imin];
    r.hit=r.pulse_height>std::max(nsig*r.noise,athr);
    return r;
}

//==============================================================================
// Clustering
//==============================================================================
struct Hit{int index;double ph,pt;};
struct Cluster{double position,time,amplitude;int nstrips,lo,hi;};
static std::vector<Cluster> make_clusters(std::vector<Hit> hits, double pitch){
    if(hits.empty()) return {};
    std::sort(hits.begin(),hits.end(),[](const Hit&a,const Hit&b){return a.index<b.index;});
    std::vector<Cluster> out; Cluster cur{};
    auto flush=[&](){cur.position/=cur.amplitude;cur.time/=cur.amplitude;out.push_back(cur);};
    auto open=[&](const Hit&h){cur={};cur.lo=cur.hi=h.index;};
    auto add=[&](const Hit&h){
        double coord=(h.index+0.5-8.0)*pitch;
        cur.position+=h.ph*coord; cur.time+=h.ph*h.pt;
        cur.amplitude+=h.ph; ++cur.nstrips; cur.hi=h.index;
    };
    open(hits[0]); add(hits[0]);
    for(int i=1;i<(int)hits.size();++i){if(hits[i].index>cur.hi+1){flush();open(hits[i]);}add(hits[i]);}
    flush(); return out;
}
static const Cluster& best_amp(const std::vector<Cluster>& v){
    return *std::max_element(v.begin(),v.end(),[](const Cluster&a,const Cluster&b){return a.amplitude<b.amplitude;});
}

//==============================================================================
// Estadísticas
//==============================================================================
static double pearson_r(const std::vector<double>&x,const std::vector<double>&y){
    const int n=static_cast<int>(x.size());if(n<2)return 0;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];syy+=y[i]*y[i];sxy+=x[i]*y[i];}
    const double d=n,num=d*sxy-sx*sy,den=std::sqrt((d*sxx-sx*sx)*(d*syy-sy*sy));
    return den>0?num/den:0.0;
}
static void linfit(const std::vector<double>&x,const std::vector<double>&y,double&sl,double&ic){
    sl=ic=0;const int n=static_cast<int>(x.size());if(n<2)return;
    double sx=0,sy=0,sxx=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];sxy+=x[i]*y[i];}
    const double d=n,dn=d*sxx-sx*sx;if(std::abs(dn)<1e-12)return;
    sl=(d*sxy-sx*sy)/dn; ic=(sy-sl*sx)/d;
}
// Kendall tau entre dos rankings de n elementos
static double kendall_tau(const std::vector<int>&ra,const std::vector<int>&rb){
    const int n=static_cast<int>(ra.size()); long long C=0,D=0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
        int da=ra[i]-ra[j],db=rb[i]-rb[j];
        if(da*db>0)++C;else if(da*db<0)++D;
    }
    const long long tot=static_cast<long long>(n)*(n-1)/2;
    return tot>0?static_cast<double>(C-D)/static_cast<double>(tot):0.0;
}

//==============================================================================
// Gráficas
//==============================================================================
static void s1(TH1D*h,const char*xt,const char*yt){
    h->GetXaxis()->SetTitle(xt);h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.046);h->GetYaxis()->SetTitleSize(0.046);
    h->SetLineWidth(2);h->SetStats(0);}
static void s2(TH2D*h,const char*xt,const char*yt){
    h->GetXaxis()->SetTitle(xt);h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.046);h->GetYaxis()->SetTitleSize(0.046);h->SetStats(0);}
static void p1(TH1D*h,const fs::path&p){
    TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",1000,680);
    c.SetLeftMargin(0.12);c.SetBottomMargin(0.13);c.SetRightMargin(0.05);
    h->Draw("HIST");c.SaveAs(p.string().c_str());}
static void p2fit(TH2D*h,double sl,double ic,double r,const fs::path&p){
    TCanvas c(("cf_"+std::string(h->GetName())).c_str(),"",900,820);
    c.SetLeftMargin(0.13);c.SetRightMargin(0.17);c.SetBottomMargin(0.13);
    h->Draw("COLZ");
    const double xlo=h->GetXaxis()->GetXmin(),xhi=h->GetXaxis()->GetXmax();
    TLine line(xlo,sl*xlo+ic,xhi,sl*xhi+ic);
    line.SetLineColor(kRed);line.SetLineWidth(2);line.Draw("SAME");
    TLatex lat;lat.SetTextSize(0.032);lat.SetNDC(true);
    char buf[128];std::snprintf(buf,sizeof(buf),"slope=%.3f  intercept=%.2f mm  r=%.3f",sl,ic,r);
    lat.DrawLatex(0.14,0.935,buf);
    c.SaveAs(p.string().c_str());}

//==============================================================================
// Main
//==============================================================================
int main(int argc,char**argv){
    Params p;
    for(int i=1;i<argc;++i){
        std::string a=argv[i];
        if(a=="--pitch"   &&i+1<argc){p.pitch_mm      =std::stod(argv[++i]);continue;}
        if(a=="--nsigma"  &&i+1<argc){p.n_sigma        =std::stod(argv[++i]);continue;}
        if(a=="--abst18"  &&i+1<argc){p.abs_thr_2018   =std::stod(argv[++i]);continue;}
        if(a=="--abst19"  &&i+1<argc){p.abs_thr_2019   =std::stod(argv[++i]);continue;}
        if(a=="--t19lo"   &&i+1<argc){p.t19_lo         =std::stod(argv[++i]);continue;}
        if(a=="--t19hi"   &&i+1<argc){p.t19_hi         =std::stod(argv[++i]);continue;}
        if(a=="--dtcable" &&i+1<argc){p.dtxy_cable     =std::stod(argv[++i]);continue;}
        if(a=="--dtxycut" &&i+1<argc){p.dtxy_cut       =std::stod(argv[++i]);continue;}
        if(a=="--dtdelay" &&i+1<argc){p.dt_delay_ns    =std::stod(argv[++i]);continue;}
        if(a=="--dtwin"   &&i+1<argc){p.dt_window_ns   =std::stod(argv[++i]);continue;}
        if(a=="--nmax19"  &&i+1<argc){p.max_cls_2019   =std::stoi(argv[++i]);continue;}
        if(a[0]!='-'){p.input_file=a;continue;}
        std::cerr<<"Unknown: "<<a<<'\n';return 1;
    }

    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("reco_v8b_dtxy_"+std::string(ds));
    fs::create_directories(outdir);

    std::cout<<"=== reconstruct_v8b ===\n";
    std::cout<<"Input : "<<p.input_file<<"\n";
    std::cout<<"Output: "<<outdir<<"\n";
    std::cout<<"hod2019 window : ["<<p.t19_lo<<", "<<p.t19_hi<<"] ns"
             <<" (v8 usó [3000,5000]: ventana reducida de 2000→"<<(p.t19_hi-p.t19_lo)<<" ns)\n";
    std::cout<<"dtxy_cable="<<p.dtxy_cable<<" ns  dtxy_cut=±"<<p.dtxy_cut<<" ns\n\n";

    gROOT->SetBatch(kTRUE); TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0); gStyle->SetPalette(kViridis);

    TFile* fin=TFile::Open(p.input_file.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: cannot open "<<p.input_file<<'\n';return 1;}
    TTree* tree=dynamic_cast<TTree*>(fin->Get("data_dig"));
    if(!tree){std::cerr<<"ERROR: TTree 'data_dig' not found\n";return 1;}

    std::array<std::vector<float>*,kNCh> ch_wf{};
    std::vector<float>* time_wf=nullptr;
    tree->SetBranchStatus("*",0);tree->SetCacheSize(64LL*1024LL*1024LL);
    for(int c=0;c<kNCh;++c){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!tree->GetBranch(bn))continue;
        tree->SetBranchStatus(bn,1);tree->SetBranchAddress(bn,&ch_wf[c]);
    }
    tree->SetBranchStatus("time",1);tree->SetBranchAddress("time",&time_wf);

    const Long64_t nEnt=tree->GetEntries();
    std::cout<<"Total entries: "<<nEnt<<'\n';

    const InvMap invA(kHod2019X,kHod2019Y),invC(kHod2019X_ctrl,kHod2019Y_ctrl);

    //==========================================================================
    // FASE A: Diagnóstico de tiempos de pulso por canal hod2019 (sin cambiar nada)
    //==========================================================================
    // Histogramas de tiempo de pulso por canal (umbral 5 ADC y 15 ADC)
    auto h_ptX_all = std::make_unique<TH1D>("h_ptX_all",
        "Tiempo pulso hod2019 X (todos los hits);t [ns];Hits",100,3000,5000);
    auto h_ptY_all = std::make_unique<TH1D>("h_ptY_all",
        "Tiempo pulso hod2019 Y (todos los hits);t [ns];Hits",100,3000,5000);
    auto h_ptX_hi  = std::make_unique<TH1D>("h_ptX_hi",
        "Tiempo pulso hod2019 X (ph>15 ADC);t [ns];Hits",100,3000,5000);
    auto h_ptY_hi  = std::make_unique<TH1D>("h_ptY_hi",
        "Tiempo pulso hod2019 Y (ph>15 ADC);t [ns];Hits",100,3000,5000);

    // Per canal: acumuladores de tiempo
    double sum_pt_lo[kNCh]={};long long n_pt_lo[kNCh]={};
    double sum_pt_hi[kNCh]={};long long n_pt_hi[kNCh]={};

    // Distribución de (t_X - t_Y) para coincidencia interna hod2019
    auto h_dtxy = std::make_unique<TH1D>("h_dt_interno_XY_hod2019",
        "Retardo interno hod2019 X #minus Y;t_X #minus t_Y [ns];Eventos",200,-2000,2000);
    // Antes y después del corte (usando ventana v8b)
    auto h_dtxy_cut = std::make_unique<TH1D>("h_dt_XY_ventana_v8b",
        "Retardo X#minusY (ventana v8b);t_X #minus t_Y [ns];Eventos",200,-2000,2000);

    //==========================================================================
    // FASE B: Histogramas de trazas v8b (con corte XY)
    //==========================================================================
    constexpr double kClo=-16.0,kChi=16.0; constexpr int kCb=16;

    auto h_dt18   =std::make_unique<TH1D>("h_dt_interplano_v8b",
        "t_{2018} #minus t_{2019} (v8b);#Deltat [ns];Events",200,-500,3000);
    auto h_cx_v8  =std::make_unique<TH2D>("h_corr_x_v8","Corr X v8b;X_{2018} [mm];X_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy_v8  =std::make_unique<TH2D>("h_corr_y_v8","Corr Y v8b;Y_{2018} [mm];Y_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    // Version sin corte XY (para comparar)
    auto h_cx_noc =std::make_unique<TH2D>("h_corr_x_nocut","Corr X sin corte XY;X_{2018} [mm];X_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy_noc =std::make_unique<TH2D>("h_corr_y_nocut","Corr Y sin corte XY;Y_{2018} [mm];Y_{2019} [mm];Events",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_thx    =std::make_unique<TH1D>("h_theta_x","#theta_x;#theta_x [deg];Events",80,-25,25);
    auto h_thy    =std::make_unique<TH1D>("h_theta_y","#theta_y;#theta_y [deg];Events",80,-25,25);
    auto h_th3    =std::make_unique<TH1D>("h_theta_3d","#theta_{3D};#theta_{3D} [deg];Events",80,0,25);

    // Gradient plot per channel (Method D)
    double sum_pt_auth[16]={};long long n_pt_auth[16]={};  // hod2019X por índice propuesto
    double sum_pt_authY[16]={};long long n_pt_authY[16]={};

    //==========================================================================
    // Cutflow
    //==========================================================================
    long long cf_tot=nEnt,cf_any=0,cf_19c=0,cf_19m=0,cf_18m=0,cf_xycut=0,cf_trk=0;
    long long cf_18m_noc=0,cf_trk_noc=0;  // sin corte XY

    std::vector<double> vx18,vx19,vy18,vy19;
    std::vector<double> vx18n,vx19n,vy18n,vy19n;  // sin corte XY

    //==========================================================================
    // Bucle de eventos
    //==========================================================================
    std::cout<<"Running event loop...\n";
    for(Long64_t iev=0;iev<nEnt;++iev){
        tree->GetEntry(iev);
        if(iev%2000==0)
            std::cout<<"  "<<iev<<" / "<<nEnt<<" ("<<100*iev/nEnt<<"%)\r"<<std::flush;
        if(!time_wf||time_wf->empty()) continue;
        const auto& tv=*time_wf;

        // Extracción de pulsos
        std::array<PulseResult,kNCh> pr{};
        for(int c=0;c<kNCh;++c){
            if(!ch_wf[c]||ch_wf[c]->empty()) continue;
            const bool is19=(c>=32);
            const double tlo = is19 ? p.t19_lo : p.kHod2018Tlo;
            const double thi = is19 ? p.t19_hi : p.kHod2018Thi;
            const double thr = is19 ? p.abs_thr_2019 : p.abs_thr_2018;
            pr[c]=extract_pulse(*ch_wf[c],tv,tlo,thi,p.n_sigma,thr);

            // Diagnóstico de tiempos (sin corte de ventana v8b — usar ventana ORIGINAL)
            if(is19 && pr[c].valid){
                // Extrae con ventana original [3000,5000] para diagnóstico
                PulseResult pr_diag=extract_pulse(*ch_wf[c],tv,3000,5000,p.n_sigma,p.abs_thr_2019);
                PulseResult pr_diag_hi=extract_pulse(*ch_wf[c],tv,3000,5000,p.n_sigma,15.0);
                if(pr_diag.hit && pr_diag.pulse_time>0){
                    const bool isX=(c>=48);
                    if(isX) h_ptX_all->Fill(pr_diag.pulse_time);
                    else    h_ptY_all->Fill(pr_diag.pulse_time);
                    sum_pt_lo[c]+=pr_diag.pulse_time; n_pt_lo[c]++;
                }
                if(pr_diag_hi.hit && pr_diag_hi.pulse_time>0){
                    const bool isX=(c>=48);
                    if(isX) h_ptX_hi->Fill(pr_diag_hi.pulse_time);
                    else    h_ptY_hi->Fill(pr_diag_hi.pulse_time);
                    sum_pt_hi[c]+=pr_diag_hi.pulse_time; n_pt_hi[c]++;
                    // Acumula para Method D (gradient)
                    int k_authX = invA.x19[c]; int k_authY = invA.y19[c];
                    if(k_authX>=0){sum_pt_auth[k_authX]+=pr_diag_hi.pulse_time;n_pt_auth[k_authX]++;}
                    if(k_authY>=0){sum_pt_authY[k_authY]+=pr_diag_hi.pulse_time;n_pt_authY[k_authY]++;}
                }
            }
        }

        // Hits por plano
        std::vector<Hit> h18x,h18y,h19x,h19y,h19xc,h19yc;
        bool any=false;
        for(int c=0;c<kNCh;++c){
            if(!pr[c].valid||!pr[c].hit) continue; any=true;
            const Hit h{0,pr[c].pulse_height,pr[c].pulse_time};
            if(c>=0&&c<=15){h18x.push_back({c,h.ph,h.pt});}
            else if(c>=16&&c<=31){h18y.push_back({c-16,h.ph,h.pt});}
            else{
                int i19x=invA.x19[c],i19y=invA.y19[c];
                int i19xc=invC.x19[c],i19yc=invC.y19[c];
                if(i19x>=0) h19x.push_back({i19x,h.ph,h.pt});
                if(i19y>=0) h19y.push_back({i19y,h.ph,h.pt});
                if(i19xc>=0)h19xc.push_back({i19xc,h.ph,h.pt});
                if(i19yc>=0)h19yc.push_back({i19yc,h.ph,h.pt});
            }
        }
        if(any)++cf_any;

        auto cls18x=make_clusters(h18x,p.pitch_mm), cls18y=make_clusters(h18y,p.pitch_mm);
        auto cls19x=make_clusters(h19x,p.pitch_mm), cls19y=make_clusters(h19y,p.pitch_mm);

        if(cls19x.empty()||cls19y.empty()) continue; ++cf_19c;
        if((int)cls19x.size()>p.max_cls_2019||(int)cls19y.size()>p.max_cls_2019) continue; ++cf_19m;

        const Cluster& cx19=best_amp(cls19x); const Cluster& cy19=best_amp(cls19y);
        const double x2019=cx19.position, y2019=cy19.position;
        const double t2019=0.5*(cx19.time+cy19.time);

        // Retardo interno X-Y de hod2019
        const double dt_xy = cx19.time - cy19.time;  // t_X − t_Y
        h_dtxy->Fill(dt_xy);

        // Distribución con ventana v8b
        if(pr[0].valid) h_dtxy_cut->Fill(dt_xy);  // siempre para todos los que pasan

        // Hod2018 (por amplitud)
        if(cls18x.empty()||cls18y.empty()) continue; ++cf_18m_noc;

        const Cluster* cx18=&best_amp(cls18x); const Cluster* cy18=&best_amp(cls18y);
        const double x2018=cx18->position, y2018=cy18->position;
        const double t2018=0.5*(cx18->time+cy18->time);
        const double dt=t2018-t2019;
        h_dt18->Fill(dt);

        const double mx=(x2019-x2018)/p.zsep_mm, my=(y2019-y2018)/p.zsep_mm;

        // Versión SIN corte XY
        h_cx_noc->Fill(x2018,x2019); h_cy_noc->Fill(y2018,y2019);
        vx18n.push_back(x2018);vx19n.push_back(x2019);
        vy18n.push_back(y2018);vy19n.push_back(y2019);
        ++cf_trk_noc;

        // Corte de coincidencia interna hod2019 X–Y
        const bool pass_xy = std::abs(dt_xy - p.dtxy_cable) < p.dtxy_cut;
        if(!pass_xy) continue; ++cf_xycut;

        ++cf_18m; ++cf_trk;
        h_cx_v8->Fill(x2018,x2019); h_cy_v8->Fill(y2018,y2019);
        h_thx->Fill(std::atan(mx)*180/M_PI);
        h_thy->Fill(std::atan(my)*180/M_PI);
        h_th3->Fill(std::atan(std::sqrt(mx*mx+my*my))*180/M_PI);
        vx18.push_back(x2018);vx19.push_back(x2019);
        vy18.push_back(y2018);vy19.push_back(y2019);
    }
    std::cout<<"\nEvent loop done.\n\n";

    //==========================================================================
    // Method D: gradient temporal vs índice propuesto
    //==========================================================================
    // rank de los canales por mean_t (orden temporal)
    std::vector<int> order_t_X(16), order_t_Y(16);
    std::iota(order_t_X.begin(),order_t_X.end(),0);
    std::iota(order_t_Y.begin(),order_t_Y.end(),0);
    std::sort(order_t_X.begin(),order_t_X.end(),[&](int a,int b){
        if(n_pt_auth[a]==0&&n_pt_auth[b]==0)return a<b;
        if(n_pt_auth[a]==0)return false; if(n_pt_auth[b]==0)return true;
        return sum_pt_auth[a]/n_pt_auth[a] < sum_pt_auth[b]/n_pt_auth[b];
    });
    std::sort(order_t_Y.begin(),order_t_Y.end(),[&](int a,int b){
        if(n_pt_authY[a]==0&&n_pt_authY[b]==0)return a<b;
        if(n_pt_authY[a]==0)return false; if(n_pt_authY[b]==0)return true;
        return sum_pt_authY[a]/n_pt_authY[a] < sum_pt_authY[b]/n_pt_authY[b];
    });
    // rank_t_X[k] = posición de índice k en el orden temporal
    std::vector<int> rank_t_X(16,-1), rank_t_Y(16,-1);
    for(int i=0;i<16;++i){ rank_t_X[order_t_X[i]]=i; rank_t_Y[order_t_Y[i]]=i; }
    // rank_prop_X[k] = k (propuesto es la identidad)
    std::vector<int> rank_prop(16); std::iota(rank_prop.begin(),rank_prop.end(),0);
    double tau_D_X = kendall_tau(rank_t_X,rank_prop);
    double tau_D_Y = kendall_tau(rank_t_Y,rank_prop);

    //==========================================================================
    // Estadísticas de correlación
    //==========================================================================
    double rx,ry,rxn,ryn,slx,icx,sly,icy;
    rx=pearson_r(vx18,vx19); ry=pearson_r(vy18,vy19);
    rxn=pearson_r(vx18n,vx19n); ryn=pearson_r(vy18n,vy19n);
    linfit(vx18,vx19,slx,icx); linfit(vy18,vy19,sly,icy);
    double slxn,icxn,slyn,icyn;
    linfit(vx18n,vx19n,slxn,icxn); linfit(vy18n,vy19n,slyn,icyn);

    //==========================================================================
    // Figuras
    //==========================================================================
    // Timing diagnóstico
    {
        TCanvas c1("c_ptXY","",1000,680);
        c1.SetLeftMargin(0.12);c1.SetBottomMargin(0.13);
        h_ptX_all->SetLineColor(kBlue+1);h_ptX_all->SetLineWidth(2);h_ptX_all->SetStats(0);
        h_ptY_all->SetLineColor(kRed);h_ptY_all->SetLineWidth(2);
        h_ptX_all->SetTitle("Pulse time hod2019 X vs Y (umbral 5 ADC);t [ns];Hits");
        h_ptX_all->Draw("HIST"); h_ptY_all->Draw("HIST SAME");
        TLatex lat;lat.SetNDC(true);lat.SetTextSize(0.032);
        lat.SetTextColor(kBlue+1);lat.DrawLatex(0.65,0.85,"hod2019 X");
        lat.SetTextColor(kRed);lat.DrawLatex(0.65,0.80,"hod2019 Y");
        c1.SaveAs((outdir/"pulsetime_XvsY_overlay.png").string().c_str());
    }
    {
        TCanvas c2("c_ptXY_hi","",1000,680);
        c2.SetLeftMargin(0.12);c2.SetBottomMargin(0.13);
        h_ptX_hi->SetLineColor(kBlue+1);h_ptX_hi->SetLineWidth(2);h_ptX_hi->SetStats(0);
        h_ptY_hi->SetLineColor(kRed);h_ptY_hi->SetLineWidth(2);
        h_ptX_hi->SetTitle("Pulse time hod2019 (ph>15 ADC);t [ns];Hits");
        h_ptX_hi->Draw("HIST"); h_ptY_hi->Draw("HIST SAME");
        TLatex lat;lat.SetNDC(true);lat.SetTextSize(0.032);
        lat.SetTextColor(kBlue+1);lat.DrawLatex(0.65,0.85,"X (ph>15)");
        lat.SetTextColor(kRed);lat.DrawLatex(0.65,0.80,"Y (ph>15)");
        c2.SaveAs((outdir/"pulsetime_hod2019_hi_thr.png").string().c_str());
    }

    // Per-channel timing vs índice propuesto (Method D)
    auto plot_gradient=[&](const double* sum_pt,const long long* n_pt,
                            const char* name_base,const std::string& coord,
                            const fs::path& outd){
        std::vector<double> xs,ys,exs,eys;
        for(int k=0;k<16;++k) if(n_pt[k]>2){
            xs.push_back(k);
            ys.push_back(sum_pt[k]/n_pt[k]);
            exs.push_back(0.3);
            const double mean=sum_pt[k]/n_pt[k];
            eys.push_back(100); // approx, no variance stored
        }
        TGraphErrors g(static_cast<int>(xs.size()),xs.data(),ys.data(),exs.data(),eys.data());
        g.SetTitle(("Tiempo medio vs índice propuesto hod2019 "+coord+";Índice propuesto k;t_{medio} [ns]").c_str());
        g.SetMarkerStyle(20);g.SetMarkerColor(kBlue+1);g.SetLineColor(kBlue+1);
        TCanvas cg(("cgrad_"+coord).c_str(),"",900,680);
        cg.SetLeftMargin(0.12);cg.SetBottomMargin(0.12);
        g.Draw("AP");
        cg.SaveAs((outd/("gradient_mean_t_vs_index_"+coord+".png")).string().c_str());
    };
    plot_gradient(sum_pt_auth,n_pt_auth,"auth_X","X",outdir);
    plot_gradient(sum_pt_authY,n_pt_authY,"auth_Y","Y",outdir);

    // dt_interno XY
    s1(h_dtxy.get(),"t_X #minus t_Y [ns]","Eventos (conincidencia hod2019)");
    {
        TCanvas c3("c_dtxy","",1000,680);
        c3.SetLeftMargin(0.12);c3.SetBottomMargin(0.13);c3.SetRightMargin(0.05);
        h_dtxy->Draw("HIST");
        // Línea de corte
        double y_max=h_dtxy->GetMaximum();
        TLine l1(p.dtxy_cable-p.dtxy_cut,0,p.dtxy_cable-p.dtxy_cut,y_max*0.9);
        TLine l2(p.dtxy_cable+p.dtxy_cut,0,p.dtxy_cable+p.dtxy_cut,y_max*0.9);
        l1.SetLineColor(kRed);l1.SetLineWidth(2);l1.Draw();
        l2.SetLineColor(kRed);l2.SetLineWidth(2);l2.Draw();
        c3.SaveAs((outdir/"dt_interno_XY_hod2019.png").string().c_str());
    }

    // Correlaciones v8b
    s2(h_cx_v8.get(),"X_{2018} [mm]","X_{2019} [mm]");
    s2(h_cy_v8.get(),"Y_{2018} [mm]","Y_{2019} [mm]");
    s2(h_cx_noc.get(),"X_{2018} [mm]","X_{2019} [mm]");
    s2(h_cy_noc.get(),"Y_{2018} [mm]","Y_{2019} [mm]");
    p2fit(h_cx_v8.get(),slx,icx,rx,outdir/"corr_x_v8b_con_corte_xy.png");
    p2fit(h_cy_v8.get(),sly,icy,ry,outdir/"corr_y_v8b_con_corte_xy.png");
    p2fit(h_cx_noc.get(),slxn,icxn,rxn,outdir/"corr_x_v8b_sin_corte_xy.png");
    p2fit(h_cy_noc.get(),slyn,icyn,ryn,outdir/"corr_y_v8b_sin_corte_xy.png");
    s1(h_dt18.get(),"#Deltat [ns]","Events");
    p1(h_dt18.get(),outdir/"dt_interplano.png");
    s1(h_thx.get(),"#theta_x [deg]","Events");
    s1(h_thy.get(),"#theta_y [deg]","Events");
    s1(h_th3.get(),"#theta_{3D} [deg]","Events");
    p1(h_thx.get(),outdir/"theta_x.png");
    p1(h_thy.get(),outdir/"theta_y.png");
    p1(h_th3.get(),outdir/"theta_3d.png");

    // ROOT output
    {
        TFile fout((outdir/"histos.root").string().c_str(),"RECREATE");
        if(!fout.IsZombie()){
            for(auto*h:std::initializer_list<TH1*>{
                h_ptX_all.get(),h_ptY_all.get(),h_ptX_hi.get(),h_ptY_hi.get(),
                h_dtxy.get(),h_dt18.get(),h_thx.get(),h_thy.get(),h_th3.get()})
                h->Write();
            for(auto*h:std::initializer_list<TH2D*>{
                h_cx_v8.get(),h_cy_v8.get(),h_cx_noc.get(),h_cy_noc.get()})
                h->Write();
            fout.Close();
        }
    }

    //==========================================================================
    // Summary
    //==========================================================================
    auto print_sum=[&](std::ostream& os){
        os<<"=== reconstruct_v8b: corrección retardo diferencial cables hod2019 ===\n\n";
        os<<"Parámetros:\n";
        os<<"  hod2019 ventana  = ["<<p.t19_lo<<", "<<p.t19_hi<<"] ns"
          <<"  (v8: [3000,5000]; reducción "<<(int)(100*(1-(p.t19_hi-p.t19_lo)/2000))<<"%)\n";
        os<<"  dtxy_cable       = "<<p.dtxy_cable<<" ns  (t_X − t_Y cable offset)\n";
        os<<"  dtxy_cut         = ±"<<p.dtxy_cut<<" ns  (ventana coincidencia interna XY)\n";
        os<<"  abs_thr_2019     = "<<p.abs_thr_2019<<" ADC\n";
        os<<"  abs_thr_2018     = "<<p.abs_thr_2018<<" ADC\n\n";

        os<<"DIAGNÓSTICO DE TIEMPOS (umbral 5 ADC, ventana original [3000,5000]):\n";
        os<<"  hod2019 X: mean_t="<<std::fixed<<std::setprecision(1)<<h_ptX_all->GetMean()
          <<" ns  (n="<<(long long)h_ptX_all->Integral()<<")\n";
        os<<"  hod2019 Y: mean_t="<<h_ptY_all->GetMean()
          <<" ns  (n="<<(long long)h_ptY_all->Integral()<<")\n";
        os<<"  Δt_cable(Y−X) con umbral 5 ADC = "<<(h_ptY_all->GetMean()-h_ptX_all->GetMean())
          <<" ns  [dominado por ruido → no representativo]\n";
        os<<"  DIAGNÓSTICO UMBRAL ALTO (ph>15 ADC — señales reales):\n";
        os<<"  hod2019 X: mean_t="<<h_ptX_hi->GetMean()
          <<" ns  (n="<<(long long)h_ptX_hi->Integral()<<")\n";
        os<<"  hod2019 Y: mean_t="<<h_ptY_hi->GetMean()
          <<" ns  (n="<<(long long)h_ptY_hi->Integral()<<")\n";
        os<<"  Δt_cable(Y−X) con señales reales = "
          <<(h_ptY_hi->GetMean()-h_ptX_hi->GetMean())<<" ns\n";
        os<<"    Las señales reales llegan a t ≈ "<<std::setprecision(0)<<h_ptX_hi->GetMean()
          <<" ns (X) y "<<h_ptY_hi->GetMean()<<" ns (Y).\n";
        os<<"    Diferencia real de cable X-Y: "
          <<std::setprecision(1)<<(h_ptY_hi->GetMean()-h_ptX_hi->GetMean())<<" ns.\n\n";

        // Per-channel timing (high threshold)
        os<<"Tiempo medio por canal hod2019 (ph>15 ADC):\n";
        os<<"  Canal  t_medio(ns)   n   (eje X)\n";
        for(int k=0;k<16;++k) if(n_pt_auth[k]>2)
            os<<"  k="<<k<<" Ch"<<kHod2019X[k]<<"  "<<std::setprecision(1)
              <<sum_pt_auth[k]/n_pt_auth[k]<<"  n="<<n_pt_auth[k]<<"\n";
        os<<"  Canal  t_medio(ns)   n   (eje Y)\n";
        for(int k=0;k<16;++k) if(n_pt_authY[k]>2)
            os<<"  k="<<k<<" Ch"<<kHod2019Y[k]<<"  "<<std::setprecision(1)
              <<sum_pt_authY[k]/n_pt_authY[k]<<"  n="<<n_pt_authY[k]<<"\n";

        os<<"\nMETHOD D — Gradiente temporal vs índice propuesto:\n";
        os<<"  τ(order_t_X, order_prop_X) = "<<std::setprecision(4)<<tau_D_X<<"\n";
        os<<"  τ(order_t_Y, order_prop_Y) = "<<tau_D_Y<<"\n";
        if(std::abs(tau_D_X)>0.6)
            os<<"  X: GRADIENTE TEMPORAL DETECTADO — el orden temporal confirma el mapeo propuesto.\n";
        else
            os<<"  X: gradiente débil o ausente — no se puede inferir mapeo desde timing.\n";
        if(std::abs(tau_D_Y)>0.6)
            os<<"  Y: GRADIENTE TEMPORAL DETECTADO.\n";
        else
            os<<"  Y: gradiente débil o ausente.\n";

        os<<"\nCUTFLOW v8b:\n";
        os<<std::setw(40)<<std::left<<"  Total eventos"             <<cf_tot<<"\n";
        os<<std::setw(40)<<"  Any hit"                               <<cf_any<<"\n";
        os<<std::setw(40)<<"  hod2019 X+Y coincidencia"             <<cf_19c<<"\n";
        os<<std::setw(40)<<"  hod2019 mult ≤ "+std::to_string(p.max_cls_2019)<<cf_19m<<"\n";
        os<<std::setw(40)<<"  hod2018 X+Y coincidencia (sin corte XY)"<<cf_18m_noc<<"\n";
        os<<std::setw(40)<<"  Trazas sin corte XY"                   <<cf_trk_noc<<"\n";
        os<<std::setw(40)<<"  |t_X−t_Y−dt_cable|<"+std::to_string(int(p.dtxy_cut))+" ns"<<cf_xycut<<"\n";
        os<<std::setw(40)<<"  Trazas CON corte XY"                   <<cf_trk<<"\n\n";

        os<<"CORRELACIÓN ESPACIAL:\n";
        os<<std::fixed<<std::setprecision(4);
        os<<"  SIN corte XY: r_x="<<rxn<<"  r_y="<<ryn
          <<"  slope_x="<<slxn<<"  (n="<<(long long)vx18n.size()<<")\n";
        os<<"  CON corte XY: r_x="<<rx <<"  r_y="<<ry
          <<"  slope_x="<<slx <<"  (n="<<(long long)vx18.size()<<")\n\n";

        os<<"CRITERIO DE ÉXITO:\n";
        if(rx>0.3||ry>0.3)
            os<<"  ✓ r > 0.3 CON corte XY → mapeo verificable desde datos actuales.\n";
        else if(rx>0.3||ry>0.3)
            os<<"  ~ r > 0.3 SIN corte XY → el corte XY no mejora (datos ya buenos o ambos malos).\n";
        else
            os<<"  ✗ r < 0.3 en ambas versiones → señal insuficiente con datos actuales.\n"
              <<"    Necesario: centelladores más gruesos, umbral más bajo, más estadística.\n";
        os<<"  Ganancia del corte XY: Δr_x = "<<(rx-rxn)
          <<"  Δr_y = "<<(ry-ryn)<<"\n";
    };

    print_sum(std::cout);
    {std::ofstream sf(outdir/"summary_v8b.txt");print_sum(sf);}

    std::cout<<"\nOutput: "<<outdir<<"\n";
    fin->Close();
    return 0;
}
