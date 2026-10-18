// analyze_highgain.cxx
// Análisis del archivo v8.root (mayor ganancia MPPC, 143k eventos).
// Diseño de pasada única: lee el archivo UNA VEZ, almacena datos relevantes en RAM,
// luego hace el barrido de umbrales y toda la inferencia en memoria.
//
// Fases:
//   0 — Caracterización: espectros ph, timing, ocupancia (en la pasada única).
//   1 — Barrido de umbral en memoria (11 valores): maximiza r(x2019,x2018).
//   2 — Inferencia de mapeo con umbral óptimo (Métodos A, B, C).
//   3 — Reconstrucción de trayectorias + corte interno XY hod2019.
//   4 — Tabla comparativa v7 vs v8.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o analyze_highgain analyze_highgain.cxx

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVectorD.h>

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
#include <random>
#include <string>
#include <vector>
#include "hodo_paths.h"

namespace fs = std::filesystem;

//==============================================================================
// Configuración
//==============================================================================
struct Params {
    std::string input_file = "lectura_vx1740_ext_nim_1M_v8.root";
    double pitch_mm     =  2.0;
    double zsep_mm      = 98.4;
    double n_sigma      =  5.0;
    double abs_thr_2018 = 20.0;
    double t19_lo=3000,t19_hi=5000,t18_lo=3500,t18_hi=7000;
    double dtxy_cut     = 300.0;
    int    max_cls      =    2;    // max clusters por coordenada
    std::vector<double> sweep = {3,5,8,10,12,15,20,25,30,40,50};
};

//==============================================================================
// Detector
//==============================================================================
constexpr double kPre=2500.0; constexpr int kNCh=64;
static const std::array<int,16> kH19X={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};
static const std::array<int,16> kH19Xc={48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};
static const std::array<int,16> kH19Yc={32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};

struct InvMap{
    std::array<int,kNCh> x{},y{};
    InvMap(const std::array<int,16>&xm,const std::array<int,16>&ym){
        x.fill(-1);y.fill(-1);
        for(int k=0;k<16;++k){x[xm[k]]=k;y[ym[k]]=k;}
    }
};

//==============================================================================
// Datos por evento (almacenados en RAM tras la pasada única)
//==============================================================================
struct ChHit { int k; double ph, pt; };  // índice físico (0..15), altura, tiempo

struct EvData {
    std::vector<ChHit> h18x,h18y,h19x,h19y;
    std::vector<ChHit> h19xc,h19yc;  // control mapping
    int ev_num; Long_t ts;
};

//==============================================================================
// Extracción de pulso
//==============================================================================
struct PR{bool valid=false;double ped=0,noise=0,ph=0,pt=-1;};
static PR extr(const std::vector<float>&wf,const std::vector<float>&t,
               double tlo,double thi,double nsig,double athr){
    PR r; const int n=static_cast<int>(wf.size());
    if(n<4||static_cast<int>(t.size())!=n)return r;
    double s=0,s2=0;int nb=0;
    for(int i=0;i<n;++i)if(t[i]<kPre){s+=wf[i];s2+=wf[i]*wf[i];++nb;}
    if(nb<2)return r;
    r.valid=true;r.ped=s/nb;
    const double var=s2/nb-r.ped*r.ped;
    r.noise=(var>0)?std::sqrt(var):0.0;
    double vmin=1e30;int imin=-1;
    for(int i=0;i<n;++i)if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0)return r;
    r.ph=r.ped-vmin; r.pt=t[imin]; return r;
}

//==============================================================================
// Clustering
//==============================================================================
struct Cluster{double pos,time,amp;int ns,lo,hi;};
static std::vector<Cluster> clust(std::vector<ChHit> hits,double pitch,double thr){
    std::vector<ChHit> active;
    for(const auto&h:hits)if(h.ph>=thr)active.push_back(h);
    if(active.empty())return{};
    std::sort(active.begin(),active.end(),[](const ChHit&a,const ChHit&b){return a.k<b.k;});
    std::vector<Cluster> out;Cluster cur{};
    auto flush=[&](){cur.pos/=cur.amp;cur.time/=cur.amp;out.push_back(cur);};
    auto open=[&](const ChHit&h){cur={};cur.lo=cur.hi=h.k;};
    auto add=[&](const ChHit&h){
        cur.pos+=h.ph*(h.k+0.5-8.0);cur.time+=h.ph*h.pt;
        cur.amp+=h.ph;++cur.ns;cur.hi=h.k;
    };
    open(active[0]);add(active[0]);
    for(int i=1;i<(int)active.size();++i){
        if(active[i].k>cur.hi+1){flush();open(active[i]);}add(active[i]);}
    flush();
    for(auto&c:out)c.pos*=pitch;
    return out;
}
static const Cluster& bamp(const std::vector<Cluster>&v){
    return *std::max_element(v.begin(),v.end(),[](const Cluster&a,const Cluster&b){return a.amp<b.amp;});
}

//==============================================================================
// Estadísticas
//==============================================================================
static double pR(const std::vector<double>&x,const std::vector<double>&y){
    const int n=static_cast<int>(x.size());if(n<2)return 0;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];syy+=y[i]*y[i];sxy+=x[i]*y[i];}
    const double d=n,num=d*sxy-sx*sy,den=std::sqrt((d*sxx-sx*sx)*(d*syy-sy*sy));
    return den>0?num/den:0.0;}
static void lf(const std::vector<double>&x,const std::vector<double>&y,double&sl,double&ic){
    sl=ic=0;const int n=static_cast<int>(x.size());if(n<2)return;
    double sx=0,sy=0,sxx=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];sxy+=x[i]*y[i];}
    const double d=n,dn=d*sxx-sx*sx;if(std::abs(dn)<1e-12)return;
    sl=(d*sxy-sx*sy)/dn;ic=(sy-sl*sx)/d;}
static double kT(const std::vector<int>&a,const std::vector<int>&b){
    const int n=static_cast<int>(a.size());long long C=0,D=0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
        int da=a[i]-a[j],db=b[i]-b[j];
        if(da*db>0)++C;else if(da*db<0)++D;}
    const long long tot=(long long)n*(n-1)/2;
    return tot>0?(double)(C-D)/(double)tot:0.0;}

//==============================================================================
// Fiedler
//==============================================================================
static std::vector<double> fiedler(const long long co[16][16],const long long hn[16]){
    const int n=16;TMatrixDSym L(n);
    double A[16][16]={};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j){
        if(i==j)continue;double d=std::sqrt((double)hn[i]*(double)hn[j]);
        A[i][j]=(d>0.5)?(double)co[i][j]/d:0.0;}
    for(int i=0;i<n;++i){double di=0;for(int j=0;j<n;++j)di+=A[i][j];
        for(int j=0;j<n;++j)L[i][j]=(i==j)?di:-A[i][j];}
    TMatrixDSymEigen e(L);
    const TVectorD&ev=e.GetEigenValues();const TMatrixD&evec=e.GetEigenVectors();
    int fc=-1;double sm=1e30,s2=1e30;
    for(int k=0;k<n;++k)if(ev(k)<sm)sm=ev(k);
    for(int k=0;k<n;++k)if(ev(k)>sm+1e-9&&ev(k)<s2){s2=ev(k);fc=k;}
    std::vector<double> fv(n,0);
    if(fc>=0)for(int i=0;i<n;++i)fv[i]=evec[i][fc];
    return fv;}

//==============================================================================
// Plot helpers
//==============================================================================
static void s1(TH1D*h,const char*x,const char*y){
    h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    h->GetXaxis()->SetTitleSize(0.046);h->GetYaxis()->SetTitleSize(0.046);
    h->SetLineWidth(2);h->SetStats(0);}
static void s2(TH2D*h,const char*x,const char*y){
    h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    h->GetXaxis()->SetTitleSize(0.046);h->GetYaxis()->SetTitleSize(0.046);h->SetStats(0);}
static void p1(TH1D*h,const fs::path&p){
    TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",1000,680);
    c.SetLeftMargin(0.12);c.SetBottomMargin(0.13);c.SetRightMargin(0.05);
    h->Draw("HIST");c.SaveAs(p.string().c_str());}
static void p2f(TH2D*h,double sl,double ic,double r,const fs::path&p){
    TCanvas c(("cf_"+std::string(h->GetName())).c_str(),"",900,820);
    c.SetLeftMargin(0.13);c.SetRightMargin(0.17);c.SetBottomMargin(0.13);
    h->Draw("COLZ");
    const double xlo=h->GetXaxis()->GetXmin(),xhi=h->GetXaxis()->GetXmax();
    TLine ln(xlo,sl*xlo+ic,xhi,sl*xhi+ic);ln.SetLineColor(kRed);ln.SetLineWidth(2);ln.Draw("SAME");
    TLatex lat;lat.SetTextSize(0.029);lat.SetNDC(true);
    char buf[128];std::snprintf(buf,sizeof(buf),"slope=%.3f  intercept=%.2f mm  r=%.3f",sl,ic,r);
    lat.DrawLatex(0.14,0.94,buf);c.SaveAs(p.string().c_str());}

//==============================================================================
// Sweep en memoria: dada lista de eventos en RAM, prueba un umbral
//==============================================================================
struct SweepRes{double thr;long long n=0;double rx=0,ry=0;};

static SweepRes sweep_in_mem(const std::vector<EvData>&evs,double thr,
                              double pitch,int max_cls){
    SweepRes res; res.thr=thr;
    std::vector<double> vx18,vx19,vy18,vy19;
    for(const auto&ev:evs){
        auto c18x=clust(ev.h18x,pitch,20.0),c18y=clust(ev.h18y,pitch,20.0);
        auto c19x=clust(ev.h19x,pitch,thr),  c19y=clust(ev.h19y,pitch,thr);
        if(c19x.empty()||c19y.empty()||c18x.empty()||c18y.empty())continue;
        if((int)c19x.size()>max_cls||(int)c19y.size()>max_cls)continue;
        if((int)c18x.size()>max_cls||(int)c18y.size()>max_cls)continue;
        ++res.n;
        vx18.push_back(bamp(c18x).pos);vx19.push_back(bamp(c19x).pos);
        vy18.push_back(bamp(c18y).pos);vy19.push_back(bamp(c19y).pos);
    }
    res.rx=pR(vx18,vx19); res.ry=pR(vy18,vy19);
    return res;
}

//==============================================================================
// main
//==============================================================================
int main(int argc,char**argv){
    Params p;
    for(int i=1;i<argc;++i){
        std::string a=argv[i];
        if(a=="--input"  &&i+1<argc){p.input_file=argv[++i];continue;}
        if(a=="--pitch"  &&i+1<argc){p.pitch_mm=std::stod(argv[++i]);continue;}
        if(a=="--abst18" &&i+1<argc){p.abs_thr_2018=std::stod(argv[++i]);continue;}
        if(a=="--dtcut"  &&i+1<argc){p.dtxy_cut=std::stod(argv[++i]);continue;}
        if(a[0]!='-'){p.input_file=a;continue;}
    }

    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("reco_v8_highgain_"+std::string(ds));
    for(auto s:{"char","mapping","tracks"})fs::create_directories(outdir/s);

    std::cout<<"=== analyze_highgain (single-pass) ===\n";
    std::cout<<"Input : "<<p.input_file<<"\n";
    std::cout<<"Output: "<<outdir<<"\n\n";

    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    const std::string resolved = hodo::resolveRawFile(p.input_file);
    TFile* fin=TFile::Open(resolved.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR\n";return 1;}
    TTree* tree=dynamic_cast<TTree*>(fin->Get("data_dig"));
    if(!tree){std::cerr<<"ERROR: no data_dig\n";return 1;}

    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw=nullptr;
    Long_t ts=0;int evn=0;
    tree->SetBranchStatus("*",0);tree->SetCacheSize(256LL<<20);
    for(int c=0;c<kNCh;++c){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!tree->GetBranch(bn))continue;
        tree->SetBranchStatus(bn,1);tree->SetBranchAddress(bn,&cw[c]);
    }
    tree->SetBranchStatus("time",1);     tree->SetBranchAddress("time",&tw);
    tree->SetBranchStatus("event",1);    tree->SetBranchAddress("event",&evn);
    tree->SetBranchStatus("timestamp",1);tree->SetBranchAddress("timestamp",&ts);

    const Long64_t nEnt=tree->GetEntries();
    std::cout<<"Total entries: "<<nEnt<<"\nReading into RAM...\n";

    const InvMap invA(kH19X,kH19Y),invC(kH19Xc,kH19Yc);

    //==========================================================================
    // Fase 0 histogramas
    //==========================================================================
    auto h_ph19x=std::make_unique<TH1D>("h19x","Pulse height hod2019 X;[ADC];Entries",100,0,100);
    auto h_ph19y=std::make_unique<TH1D>("h19y","Pulse height hod2019 Y;[ADC];Entries",100,0,100);
    auto h_ph18x=std::make_unique<TH1D>("h18x","Pulse height hod2018 X;[ADC];Entries",100,0,400);
    auto h_ph18y=std::make_unique<TH1D>("h18y","Pulse height hod2018 Y;[ADC];Entries",100,0,400);
    auto h_pt19x=std::make_unique<TH1D>("pt19x","Pulse time hod2019 X (ph>15);[ns];Hits",100,3000,5500);
    auto h_pt19y=std::make_unique<TH1D>("pt19y","Pulse time hod2019 Y (ph>15);[ns];Hits",100,3000,5500);
    auto h_occ19x=std::make_unique<TH1D>("occ19x","Ocupancia hod2019 X;k;Hits",16,-0.5,15.5);
    auto h_occ19y=std::make_unique<TH1D>("occ19y","Ocupancia hod2019 Y;k;Hits",16,-0.5,15.5);

    // Acumuladores para Fase 0 y Método A/B
    double sum_ph19x=0,sum_ph19y=0;
    long long n_ph19x=0,n_ph19y=0;
    double max_ph[kNCh]={};

    // Método A datos de sesión (umbral temporal mínimo para almacenar)
    constexpr double kMinPh = 1.0;  // umbral mínimo para almacenar en RAM

    // Almacenamiento en RAM
    std::vector<EvData> all_evs;
    all_evs.reserve(nEnt);

    //==========================================================================
    // PASADA ÚNICA
    //==========================================================================
    for(Long64_t iev=0;iev<nEnt;++iev){
        tree->GetEntry(iev);
        if(iev%20000==0)std::cout<<"  "<<iev<<"/"<<nEnt<<"\r"<<std::flush;
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;

        EvData ev; ev.ev_num=evn; ev.ts=ts;
        bool any=false;

        for(int c=0;c<kNCh;++c){
            if(!cw[c]||cw[c]->empty())continue;
            const bool is19=(c>=32);
            const double tlo=is19?p.t19_lo:p.t18_lo,thi=is19?p.t19_hi:p.t18_hi;
            // Extrae con umbral mínimo (para histogramas y almacenamiento)
            PR pr=extr(*cw[c],tv,tlo,thi,p.n_sigma,0.5);
            if(!pr.valid)continue;

            const double ph=pr.ph;
            if(ph>max_ph[c])max_ph[c]=ph;

            // Fase 0 histogramas
            if(c>=48){h_ph19x->Fill(ph);sum_ph19x+=ph;n_ph19x++;
                int k=invA.x[c]; if(k>=0&&ph>5.0)h_occ19x->Fill(k);
                if(ph>15&&pr.pt>0)h_pt19x->Fill(pr.pt);}
            else if(c>=32){h_ph19y->Fill(ph);sum_ph19y+=ph;n_ph19y++;
                int k=invA.y[c]; if(k>=0&&ph>5.0)h_occ19y->Fill(k);
                if(ph>15&&pr.pt>0)h_pt19y->Fill(pr.pt);}
            else if(c<=15)h_ph18x->Fill(ph);
            else h_ph18y->Fill(ph);

            // Almacena en RAM si supera umbral mínimo
            if(ph<kMinPh)continue;
            any=true;

            const ChHit h{0,ph,pr.pt};  // índice k se asigna debajo
            if(c>=48){
                int k=invA.x[c];if(k>=0){ev.h19x.push_back({k,ph,pr.pt});}
                int kc=invC.x[c];if(kc>=0)ev.h19xc.push_back({kc,ph,pr.pt});
            }else if(c>=32){
                int k=invA.y[c];if(k>=0){ev.h19y.push_back({k,ph,pr.pt});}
                int kc=invC.y[c];if(kc>=0)ev.h19yc.push_back({kc,ph,pr.pt});
            }else if(c<=15){ev.h18x.push_back({c,ph,pr.pt});}
            else{ev.h18y.push_back({c-16,ph,pr.pt});}
        }
        if(any)all_evs.push_back(std::move(ev));
    }
    fin->Close();
    std::cout<<"\nPassada completa. Eventos en RAM: "<<all_evs.size()<<"\n\n";

    const double ph19x_mean=n_ph19x>0?sum_ph19x/n_ph19x:0;
    const double ph19y_mean=n_ph19y>0?sum_ph19y/n_ph19y:0;

    //==========================================================================
    // FASE 0 — Resultados
    //==========================================================================
    std::cout<<"=== FASE 0: Caracterización ===\n";
    std::cout<<std::setw(25)<<"Subconjunto"<<std::setw(14)<<"v7[ADC]"<<std::setw(14)<<"v8[ADC]"<<std::setw(10)<<"ratio\n";
    std::cout<<"-------------------------------------------------------\n";
    std::cout<<std::setw(25)<<"hod2019 X (ph_medio)"<<std::setw(14)<<"4.4"
              <<std::setw(14)<<std::fixed<<std::setprecision(2)<<ph19x_mean
              <<std::setw(10)<<ph19x_mean/4.4<<"\n";
    std::cout<<std::setw(25)<<"hod2019 Y (ph_medio)"<<std::setw(14)<<"4.4"
              <<std::setw(14)<<ph19y_mean<<std::setw(10)<<ph19y_mean/4.4<<"\n";
    double mx19x=0,mx19y=0;
    for(int c=48;c<64;++c)mx19x=std::max(mx19x,max_ph[c]);
    for(int c=32;c<48;++c)mx19y=std::max(mx19y,max_ph[c]);
    std::cout<<"  max_ph hod2019 X: "<<mx19x<<" ADC   Y: "<<mx19y<<" ADC\n";
    std::cout<<"  Timing (ph>15): X="<<h_pt19x->GetMean()<<" ns  Y="<<h_pt19y->GetMean()<<" ns\n";
    const bool ganancia_ok=(ph19x_mean>6.0||mx19x>50.0);
    std::cout<<">>> GANANCIA "<<(ganancia_ok?"MEJORÓ":"NO MEJORÓ")
             <<" (ph_medio="<<ph19x_mean<<" ADC, max="<<mx19x<<" ADC)\n";
    std::cout<<"    Distribución sesgada: max >> media → mayoría de eventos sin señal real en hod2019.\n\n";

    // Umbral natural (valle del espectro)
    double umbral_nat=10.0;
    {int vb=2;double vv=1e30;
     for(int b=2;b<=40;++b){double v=h_ph19x->GetBinContent(b);if(v>0&&v<vv){vv=v;vb=b;}}
     umbral_nat=h_ph19x->GetBinCenter(vb);}
    std::cout<<"  Umbral natural estimado: "<<umbral_nat<<" ADC\n\n";

    // Guarda PNGs Fase 0
    s1(h_ph19x.get(),"Pulse height [ADC]","Entries");
    s1(h_ph19y.get(),"Pulse height [ADC]","Entries");
    s1(h_ph18x.get(),"Pulse height [ADC]","Entries");
    s1(h_ph18y.get(),"Pulse height [ADC]","Entries");
    {TCanvas cc("c_ph","",1200,900);cc.Divide(2,2);
     for(int i=0;i<4;++i){cc.cd(i+1);cc.GetPad(i+1)->SetLeftMargin(0.13);cc.GetPad(i+1)->SetBottomMargin(0.13);
         TH1D* h=(i==0)?h_ph19x.get():(i==1)?h_ph19y.get():(i==2)?h_ph18x.get():h_ph18y.get();
         h->Draw("HIST");}
     cc.SaveAs((outdir/"char"/"char_pulseheight_v8_highgain.png").string().c_str());}
    {TCanvas cc("c_occ","",1200,500);cc.Divide(2,1);
     cc.cd(1);s1(h_occ19x.get(),"Strip k","Hits");h_occ19x->Draw("HIST");
     cc.cd(2);s1(h_occ19y.get(),"Strip k","Hits");h_occ19y->Draw("HIST");
     cc.SaveAs((outdir/"char"/"char_occupancy_v8_highgain.png").string().c_str());}
    {TCanvas cc("c_pt","",1200,500);cc.Divide(2,1);
     cc.cd(1);s1(h_pt19x.get(),"t [ns]","Hits");h_pt19x->Draw("HIST");
     cc.cd(2);s1(h_pt19y.get(),"t [ns]","Hits");h_pt19y->Draw("HIST");
     cc.SaveAs((outdir/"char"/"char_pulsetime_v8_highgain.png").string().c_str());}

    //==========================================================================
    // FASE 1 — Barrido en memoria
    //==========================================================================
    std::cout<<"=== FASE 1: Barrido de umbral (en RAM) ===\n";
    std::cout<<std::setw(10)<<"abst19"<<std::setw(10)<<"N_clean"<<std::setw(10)<<"r_x"<<std::setw(10)<<"r_y\n";
    std::cout<<"----------------------------------------------\n";

    std::vector<SweepRes> sweep_results;
    double abst19_opt=5.0;double best_r=0;

    for(double thr:p.sweep){
        SweepRes sr=sweep_in_mem(all_evs,thr,p.pitch_mm,p.max_cls);
        sweep_results.push_back(sr);
        std::cout<<std::setw(10)<<thr<<std::setw(10)<<sr.n
                  <<std::setw(10)<<std::fixed<<std::setprecision(4)<<sr.rx
                  <<std::setw(10)<<sr.ry<<"\n";
        if(sr.n>=50&&std::max(std::abs(sr.rx),std::abs(sr.ry))>best_r){
            best_r=std::max(std::abs(sr.rx),std::abs(sr.ry));abst19_opt=thr;}
    }
    std::cout<<">>> abst19_opt = "<<abst19_opt<<" ADC  (max|r| = "<<best_r<<")\n\n";

    //==========================================================================
    // FASE 2 — Inferencia de mapeo (con abst19_opt)
    //==========================================================================
    std::cout<<"=== FASE 2: Inferencia de mapeo (abst19_opt="<<abst19_opt<<") ===\n";

    long long nRX[16]={},nRY[16]={},nAX[16]={},nAY[16]={};
    double sXR[16]={},sX2R[16]={},sYR[16]={},sY2R[16]={};
    long long coX[16][16]={},hnX[16]={};
    long long coY[16][16]={},hnY[16]={};
    long long occ19X[16]={},occ19Y[16]={};
    std::vector<std::pair<int,double>> scX,scY;  // (k_idx_19, ref_18 coord)

    for(const auto&ev:all_evs){
        auto c18x=clust(ev.h18x,p.pitch_mm,p.abs_thr_2018);
        auto c18y=clust(ev.h18y,p.pitch_mm,p.abs_thr_2018);
        auto c19x=clust(ev.h19x,p.pitch_mm,abst19_opt);
        auto c19y=clust(ev.h19y,p.pitch_mm,abst19_opt);

        // Ocupancia y co-activación
        for(const auto&cl:c19x){occ19X[cl.lo]++;hnX[cl.lo]++;}
        for(const auto&cl:c19y){occ19Y[cl.lo]++;hnY[cl.lo]++;}
        for(size_t i=0;i<c19x.size();++i)for(size_t j=i+1;j<c19x.size();++j){
            coX[c19x[i].lo][c19x[j].lo]++;coX[c19x[j].lo][c19x[i].lo]++;}
        for(size_t i=0;i<c19y.size();++i)for(size_t j=i+1;j<c19y.size();++j){
            coY[c19y[i].lo][c19y[j].lo]++;coY[c19y[j].lo][c19y[i].lo]++;}

        // Método A relajado: 1 activo en hod2019X + lead hod2018X (ph>30)
        if(c19x.size()==1){
            int ki=c19x[0].lo;
            int best18=-1;double bph=0;
            for(const auto&h:ev.h18x)if(h.ph>bph){bph=h.ph;best18=h.k;}
            if(best18>=0&&bph>=30.0){
                double x18=(best18+0.5-8.0)*p.pitch_mm;
                nRX[ki]++;sXR[ki]+=x18;sX2R[ki]+=x18*x18;
                scX.push_back({ki,x18});}
            if(c18x.size()==1){double x18=(ev.h18x[0].k+0.5-8.0)*p.pitch_mm;
                nAX[ki]++;}}
        if(c19y.size()==1){
            int ki=c19y[0].lo;
            int best18=-1;double bph=0;
            for(const auto&h:ev.h18y)if(h.ph>bph){bph=h.ph;best18=h.k;}
            if(best18>=0&&bph>=30.0){
                double y18=(best18+0.5-8.0)*p.pitch_mm;
                nRY[ki]++;sYR[ki]+=y18;sY2R[ki]+=y18*y18;
                scY.push_back({ki,y18});}
            if(c18y.size()==1)nAY[c19y[0].lo]++;}
    }

    // Construye ordenamiento π_A
    double xhat[16]={},errX[16]={};
    double yhat[16]={},errY[16]={};
    std::vector<int> ordAX(16),ordAY(16);
    std::iota(ordAX.begin(),ordAX.end(),0);
    std::iota(ordAY.begin(),ordAY.end(),0);
    for(int c=0;c<16;++c){
        if(nRX[c]>0){xhat[c]=sXR[c]/nRX[c];double v=sX2R[c]/nRX[c]-xhat[c]*xhat[c];
            errX[c]=(v>0)?std::sqrt(v)/std::sqrt((double)nRX[c]):1e9;}else{xhat[c]=1e9;errX[c]=1e9;}
        if(nRY[c]>0){yhat[c]=sYR[c]/nRY[c];double v=sY2R[c]/nRY[c]-yhat[c]*yhat[c];
            errY[c]=(v>0)?std::sqrt(v)/std::sqrt((double)nRY[c]):1e9;}else{yhat[c]=1e9;errY[c]=1e9;}
    }
    std::stable_sort(ordAX.begin(),ordAX.end(),[&](int a,int b){
        if(nRX[a]==0&&nRX[b]==0)return a<b;if(nRX[a]==0)return false;if(nRX[b]==0)return true;
        return xhat[a]<xhat[b];});
    std::stable_sort(ordAY.begin(),ordAY.end(),[&](int a,int b){
        if(nRY[a]==0&&nRY[b]==0)return a<b;if(nRY[a]==0)return false;if(nRY[b]==0)return true;
        return yhat[a]<yhat[b];});

    // Rangos para Kendall τ
    std::vector<int> rAX(16,-1),rAY(16,-1),r0X(16,-1),r0Y(16,-1),rBX(16,-1),rBY(16,-1);
    for(int k=0;k<16;++k){rAX[ordAX[k]]=k;rAY[ordAY[k]]=k;}
    for(int k=0;k<16;++k){r0X[kH19X[k]-48]=k;r0Y[kH19Y[k]-32]=k;}

    // Método B: Fiedler
    auto fvX=fiedler(coX,hnX),fvY=fiedler(coY,hnY);
    std::vector<int> ordBX(16),ordBY(16);
    std::iota(ordBX.begin(),ordBX.end(),0);std::iota(ordBY.begin(),ordBY.end(),0);
    std::sort(ordBX.begin(),ordBX.end(),[&](int a,int b){return fvX[a]<fvX[b];});
    std::sort(ordBY.begin(),ordBY.end(),[&](int a,int b){return fvY[a]<fvY[b];});
    for(int k=0;k<16;++k){rBX[ordBX[k]]=k;rBY[ordBY[k]]=k;}

    const double tau_AX=kT(rAX,r0X),tau_AY=kT(rAY,r0Y);
    const double tau_BX=kT(rBX,r0X),tau_BY=kT(rBY,r0Y);
    const double tau_ABX=kT(rAX,rBX),tau_ABY=kT(rAY,rBY);

    // Pearson r con π₀ vs π_A en eventos Método A
    std::vector<double> vx18p0,vx19p0,vx18pA,vx19pA;
    std::vector<double> vy18p0,vy19p0,vy18pA,vy19pA;
    for(const auto&[ki,x18]:scX){
        int k0=r0X[ki];if(k0>=0){vx18p0.push_back(x18);vx19p0.push_back((k0+0.5-8.0)*p.pitch_mm);}
        int kA=rAX[ki];if(kA>=0){vx18pA.push_back(x18);vx19pA.push_back((kA+0.5-8.0)*p.pitch_mm);}}
    for(const auto&[ki,y18]:scY){
        int k0=r0Y[ki];if(k0>=0){vy18p0.push_back(y18);vy19p0.push_back((k0+0.5-8.0)*p.pitch_mm);}
        int kA=rAY[ki];if(kA>=0){vy18pA.push_back(y18);vy19pA.push_back((kA+0.5-8.0)*p.pitch_mm);}}
    const double r_pi0X=pR(vx18p0,vx19p0),r_piAX=pR(vx18pA,vx19pA);
    const double r_pi0Y=pR(vy18p0,vy19p0),r_piAY=pR(vy18pA,vy19pA);

    // Método C: rugosidad
    auto rg=[&](const long long*occ,const std::vector<int>&ord)->double{
        double R=0;for(int k=0;k<15;++k){double d=occ[ord[k+1]]-occ[ord[k]];R+=d*d;}return R;};
    std::vector<int> ord0X(16),ord0Y(16);
    for(int k=0;k<16;++k){ord0X[k]=kH19X[k]-48;ord0Y[k]=kH19Y[k]-32;}
    std::mt19937 rng(42);std::vector<int> perm(16);std::iota(perm.begin(),perm.end(),0);
    double randRX=0,randRY=0;
    for(int t=0;t<500;++t){std::shuffle(perm.begin(),perm.end(),rng);
        randRX+=rg(occ19X,perm);randRY+=rg(occ19Y,perm);}
    randRX/=500;randRY/=500;

    // Imprime tabla Método A — X
    std::cout<<"Tabla Método A — hod2019 X:\n";
    std::cout<<std::setw(6)<<"Canal"<<std::setw(9)<<"x̂_c"<<std::setw(9)<<"±err"
              <<std::setw(7)<<"n_rel"<<std::setw(7)<<"n_cl"<<std::setw(6)<<"π_A"<<std::setw(6)<<"π₀"<<"\n";
    for(int k=0;k<16;++k){
        int c=ordAX[k];int ch=c+48;int k0=r0X[c];
        std::cout<<std::setw(6)<<("Ch"+std::to_string(ch));
        if(nRX[c]>0)std::cout<<std::fixed<<std::setprecision(1)
            <<std::setw(9)<<xhat[c]<<std::setw(9)<<errX[c];
        else std::cout<<std::setw(9)<<"N/A"<<std::setw(9)<<"---";
        // Ambigüedad
        bool amb=false;if(nRX[c]>0)for(int c2=0;c2<16;++c2){
            if(c2==c||nRX[c2]==0)continue;
            if(std::abs(xhat[c]-xhat[c2])<2*(errX[c]+errX[c2])){amb=true;break;}}
        std::cout<<std::setw(7)<<nRX[c]<<std::setw(7)<<nAX[c]<<std::setw(6)<<k<<std::setw(6)<<k0
                  <<(amb?" AMB":"")<<"\n";
    }
    std::cout<<"\nKendall τ (X): τ(A,0)="<<std::fixed<<std::setprecision(4)<<tau_AX
             <<"  τ(B,0)="<<tau_BX<<"  τ(A,B)="<<tau_ABX<<"\n";
    std::cout<<"Kendall τ (Y): τ(A,0)="<<tau_AY<<"  τ(B,0)="<<tau_BY<<"  τ(A,B)="<<tau_ABY<<"\n";
    std::cout<<"Pearson r (Mét.A): r(π₀)_X="<<r_pi0X<<"  r(π_A)_X="<<r_piAX
             <<"  r(π₀)_Y="<<r_pi0Y<<"  r(π_A)_Y="<<r_piAY<<"\n";
    std::cout<<"Rugosidad: R(π₀)_X="<<std::setprecision(0)<<rg(occ19X,ord0X)
             <<"  R(π_A)_X="<<rg(occ19X,ordAX)<<"  R(rand)_X="<<randRX<<"\n\n";

    auto verdict=[](double tA,double tB,const std::string& ax)->std::string{
        const double best=std::max(std::abs(tA),std::abs(tB));
        if(best>0.75&&tA*tB>0)return "MAPEO "+ax+" CONFIRMADO (τ_A="+std::to_string(tA).substr(0,6)+")";
        if(std::abs(tA)<0.3)return "MAPEO "+ax+" INCIERTO — estadística/SNR insuficiente";
        return "MAPEO "+ax+" PARCIALMENTE CORREGIDO";
    };


    //==========================================================================
    // FASE 3 — Reconstrucción con abst19_opt + corte XY
    //==========================================================================
    std::cout<<"=== FASE 3: Reconstrucción ===\n";
    constexpr double kClo=-16,kChi=16;constexpr int kCb=16;
    long long cf_19c=0,cf_19m=0,cf_18c=0,cf_xy=0,cf_trk=0;
    const long long n_total=nEnt;

    auto h_dt=std::make_unique<TH1D>("h_dt","t_{2018}#minust_{2019};[ns];Events",200,-500,3000);
    auto h_dtxy=std::make_unique<TH1D>("h_dtxy","t_X#minust_Y hod2019;[ns];Events",200,-2000,2000);
    auto h_cx=std::make_unique<TH2D>("h_cx","Corr X auth;X_{2018}[mm];X_{2019}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy=std::make_unique<TH2D>("h_cy","Corr Y auth;Y_{2018}[mm];Y_{2019}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cx_c=std::make_unique<TH2D>("h_cxc","Corr X ctrl;X_{2018}[mm];X_{2019}^{ctrl}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy_c=std::make_unique<TH2D>("h_cyc","Corr Y ctrl;Y_{2018}[mm];Y_{2019}^{ctrl}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_thx=std::make_unique<TH1D>("h_thx","#theta_x;[deg];Events",80,-25,25);
    auto h_thy=std::make_unique<TH1D>("h_thy","#theta_y;[deg];Events",80,-25,25);
    auto h_th3=std::make_unique<TH1D>("h_th3","#theta_{3D};[deg];Events",80,0,25);
    auto h_xy18=std::make_unique<TH2D>("h_xy18","XY 2018;X[mm];Y[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_xy19=std::make_unique<TH2D>("h_xy19","XY 2019;X[mm];Y[mm]",kCb,kClo,kChi,kCb,kClo,kChi);

    std::vector<double> vx18,vx19,vy18,vy19,vxc,vyc,vxc19,vyc19;
    std::ofstream fcsv(outdir/"recon_tracks.csv");
    fcsv<<"event,ts,x18,y18,t18,x19,y19,t19,dt,mx,my,thx,thy,th3\n";

    for(const auto&ev:all_evs){
        auto c18x=clust(ev.h18x,p.pitch_mm,p.abs_thr_2018);
        auto c18y=clust(ev.h18y,p.pitch_mm,p.abs_thr_2018);
        auto c19x=clust(ev.h19x,p.pitch_mm,abst19_opt);
        auto c19y=clust(ev.h19y,p.pitch_mm,abst19_opt);
        auto c19xc=clust(ev.h19xc,p.pitch_mm,abst19_opt);
        auto c19yc=clust(ev.h19yc,p.pitch_mm,abst19_opt);

        if(c19x.empty()||c19y.empty())continue;++cf_19c;
        if((int)c19x.size()>p.max_cls||(int)c19y.size()>p.max_cls)continue;++cf_19m;
        if(c18x.empty()||c18y.empty())continue;++cf_18c;

        const Cluster& cx19=bamp(c19x),&cy19=bamp(c19y);
        const Cluster& cx18=bamp(c18x),&cy18=bamp(c18y);
        const double x19=cx19.pos,y19=cy19.pos,x18=cx18.pos,y18=cy18.pos;
        const double t19=0.5*(cx19.time+cy19.time),t18=0.5*(cx18.time+cy18.time);
        const double dt=t18-t19,dt_xy=cx19.time-cy19.time;
        h_dt->Fill(dt);h_dtxy->Fill(dt_xy);
        h_cx->Fill(x18,x19);h_cy->Fill(y18,y19);
        vx18.push_back(x18);vx19.push_back(x19);
        vy18.push_back(y18);vy19.push_back(y19);
        if(!c19xc.empty()&&!c19yc.empty()&&(int)c19xc.size()<=p.max_cls&&(int)c19yc.size()<=p.max_cls){
            h_cx_c->Fill(x18,bamp(c19xc).pos);h_cy_c->Fill(y18,bamp(c19yc).pos);
            vxc.push_back(x18);vxc19.push_back(bamp(c19xc).pos);
            vyc.push_back(y18);vyc19.push_back(bamp(c19yc).pos);}

        const bool pxy=std::abs(dt_xy)<p.dtxy_cut;
        if(!pxy)continue;++cf_xy;++cf_trk;
        const double mx=(x19-x18)/p.zsep_mm,my=(y19-y18)/p.zsep_mm;
        const double thx=std::atan(mx)*180/M_PI,thy=std::atan(my)*180/M_PI;
        const double th3=std::atan(std::sqrt(mx*mx+my*my))*180/M_PI;
        h_thx->Fill(thx);h_thy->Fill(thy);h_th3->Fill(th3);
        h_xy18->Fill(x18,y18);h_xy19->Fill(x19,y19);
        fcsv<<ev.ev_num<<','<<ev.ts<<','<<std::fixed<<std::setprecision(3)
            <<x18<<','<<y18<<','<<t18<<','<<x19<<','<<y19<<','<<t19<<','
            <<dt<<','<<mx<<','<<my<<','<<thx<<','<<thy<<','<<th3<<"\n";
    }

    const double rx=pR(vx18,vx19),ry=pR(vy18,vy19);
    const double rxc=pR(vxc,vxc19),ryc=pR(vyc,vyc19);
    double slx,icx,sly,icy;lf(vx18,vx19,slx,icx);lf(vy18,vy19,sly,icy);
    double slxc,icxc,slyc,icyc;lf(vxc,vxc19,slxc,icxc);lf(vyc,vyc19,slyc,icyc);

    std::cout<<"Cutflow: tot="<<n_total<<" → 19coinc="<<cf_19c<<" → 19mult="<<cf_19m
             <<" → 18coinc="<<cf_18c<<" → XYcut="<<cf_xy<<" → full3D="<<cf_trk<<"\n";
    std::cout<<"Auth:  r_x="<<std::fixed<<std::setprecision(4)<<rx<<"  r_y="<<ry
             <<"  slope_x="<<slx<<"  n="<<(long long)vx18.size()<<"\n";
    std::cout<<"Ctrl:  r_x="<<rxc<<"  r_y="<<ryc<<"\n";

    //==========================================================================
    // FASE 4 — Comparación
    //==========================================================================
    const SweepRes& best_sw=*std::max_element(sweep_results.begin(),sweep_results.end(),
        [](const SweepRes&a,const SweepRes&b){return std::max(std::abs(a.rx),std::abs(a.ry))<std::max(std::abs(b.rx),std::abs(b.ry));});

    std::cout<<"\n=== FASE 4: Comparativa v7 vs v8 ===\n";
    std::cout<<std::setw(30)<<std::left<<"Métrica"<<std::setw(20)<<"v7.root (17k)"<<std::setw(20)<<"v8.root (140k)\n";
    std::cout<<"----------------------------------------------------------------------\n";
    auto row=[&](const std::string&m,const std::string&v7,const std::string&v8){
        std::cout<<std::setw(30)<<m<<std::setw(20)<<v7<<std::setw(20)<<v8<<"\n";};
    row("ph_medio hod2019X [ADC]","4.4",std::to_string(ph19x_mean).substr(0,5));
    row("max_ph hod2019X [ADC]","~198",std::to_string((int)mx19x));
    row("abst19_opt [ADC]","N/A",std::to_string((int)abst19_opt));
    row("N_clean (abst19_opt)","~52",std::to_string(best_sw.n));
    row("r_x (tracks + XY cut)","0.069",std::to_string(rx).substr(0,6));
    row("r_y (tracks + XY cut)","−0.059",std::to_string(ry).substr(0,6));
    row("τ_A_X","0.12",std::to_string(tau_AX).substr(0,6));
    row("τ_A_Y","−0.17",std::to_string(tau_AY).substr(0,6));


    //==========================================================================
    // Guarda figuras y archivos
    //==========================================================================
    // Mapping CSVs
    {std::ofstream f(outdir/"mapping_used.csv");f<<"plane,coord,index,channel\n";
     for(int k=0;k<16;++k)f<<"hod2019,X,"<<k<<','<<kH19X[k]<<'\n';
     for(int k=0;k<16;++k)f<<"hod2019,Y,"<<k<<','<<kH19Y[k]<<'\n';}
    {std::ofstream f(outdir/"mapping_inferred.csv");f<<"plane,coord,index,channel\n";
     for(int k=0;k<16;++k)f<<"hod2019,X,"<<k<<','<<(ordAX[k]+48)<<'\n';
     for(int k=0;k<16;++k)f<<"hod2019,Y,"<<k<<','<<(ordAY[k]+32)<<'\n';}

    // Scatter plots
    {auto hs0=std::make_unique<TH2D>("hs0","Scatter X π₀;x18[mm];x19[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
     auto hsA=std::make_unique<TH2D>("hsA","Scatter X π_A;x18[mm];x19[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
     for(size_t i=0;i<vx18p0.size();++i)hs0->Fill(vx18p0[i],vx19p0[i]);
     for(size_t i=0;i<vx18pA.size();++i)hsA->Fill(vx18pA[i],vx19pA[i]);
     double sl,ic;lf(vx18p0,vx19p0,sl,ic);p2f(hs0.get(),sl,ic,r_pi0X,outdir/"mapping"/"scatter_x2019_vs_x2018_pi0.png");
     lf(vx18pA,vx19pA,sl,ic);p2f(hsA.get(),sl,ic,r_piAX,outdir/"mapping"/"scatter_x2019_vs_x2018_piA.png");}

    // Posición estimada
    {std::vector<double> xs,ys,exs,eys;
     for(int k=0;k<16;++k){int c=ordAX[k];if(nRX[c]<2)continue;
         xs.push_back(k);ys.push_back(xhat[c]);exs.push_back(0.3);eys.push_back(errX[c]);}
     TGraphErrors g(static_cast<int>(xs.size()),xs.data(),ys.data(),exs.data(),eys.data());
     g.SetTitle("Posición estimada hod2019 X;Índice π_A;x̂_c [mm]");
     g.SetMarkerStyle(20);g.SetMarkerColor(kBlue+1);g.SetLineColor(kBlue+1);
     TCanvas cg("cg","",900,680);cg.SetLeftMargin(0.12);cg.SetBottomMargin(0.12);
     g.Draw("AP");
     TLine ref(-0.5,(0+0.5-8.0)*p.pitch_mm,15.5,(15+0.5-8.0)*p.pitch_mm);
     ref.SetLineColor(kRed);ref.SetLineStyle(2);ref.SetLineWidth(2);ref.Draw("SAME");
     cg.SaveAs((outdir/"mapping"/"posicion_estimada_X.png").string().c_str());}

    // Correlaciones
    s2(h_cx.get(),"X_{2018}[mm]","X_{2019}[mm]");s2(h_cy.get(),"Y_{2018}[mm]","Y_{2019}[mm]");
    s2(h_cx_c.get(),"X_{2018}[mm]","X_{2019}^{ctrl}[mm]");s2(h_cy_c.get(),"Y_{2018}[mm]","Y_{2019}^{ctrl}[mm]");
    p2f(h_cx.get(),slx,icx,rx,outdir/"tracks"/"corr_x2019_vs_x2018.png");
    p2f(h_cy.get(),sly,icy,ry,outdir/"tracks"/"corr_y2019_vs_y2018.png");
    p2f(h_cx_c.get(),slxc,icxc,rxc,outdir/"tracks"/"corr_x2019_vs_x2018_ctrl.png");
    p2f(h_cy_c.get(),slyc,icyc,ryc,outdir/"tracks"/"corr_y2019_vs_y2018_ctrl.png");
    s1(h_thx.get(),"#theta_x[deg]","Events");s1(h_thy.get(),"#theta_y[deg]","Events");
    s1(h_th3.get(),"#theta_{3D}[deg]","Events");
    p1(h_thx.get(),outdir/"tracks"/"theta_x.png");p1(h_thy.get(),outdir/"tracks"/"theta_y.png");
    p1(h_th3.get(),outdir/"tracks"/"theta_3d.png");
    s1(h_dt.get(),"#Deltat[ns]","Events");s1(h_dtxy.get(),"t_X#minust_Y[ns]","Events");
    p1(h_dt.get(),outdir/"tracks"/"dt_interplano.png");p1(h_dtxy.get(),outdir/"tracks"/"dt_interno_XY.png");
    {TCanvas c2("c2","",900,820);c2.SetLeftMargin(0.13);c2.SetRightMargin(0.17);c2.SetBottomMargin(0.13);
     s2(h_xy18.get(),"X[mm]","Y[mm]");h_xy18->Draw("COLZ");c2.SaveAs((outdir/"tracks"/"xy_2018.png").string().c_str());}
    {TCanvas c2("c2b","",900,820);c2.SetLeftMargin(0.13);c2.SetRightMargin(0.17);c2.SetBottomMargin(0.13);
     s2(h_xy19.get(),"X[mm]","Y[mm]");h_xy19->Draw("COLZ");c2.SaveAs((outdir/"tracks"/"xy_2019.png").string().c_str());}

    // ROOT
    {TFile fout((outdir/"histos.root").string().c_str(),"RECREATE");
     if(!fout.IsZombie()){
         for(auto*h:std::initializer_list<TH1*>{
             h_ph19x.get(),h_ph19y.get(),h_ph18x.get(),h_ph18y.get(),
             h_pt19x.get(),h_pt19y.get(),h_occ19x.get(),h_occ19y.get(),
             h_dt.get(),h_dtxy.get(),h_thx.get(),h_thy.get(),h_th3.get()})h->Write();
         for(auto*h:std::initializer_list<TH2D*>{h_cx.get(),h_cy.get(),
             h_cx_c.get(),h_cy_c.get(),h_xy18.get(),h_xy19.get()})h->Write();
         fout.Close();}}

    // summary_highgain.txt
    {std::ofstream sf(outdir/"summary_highgain.txt");
     sf<<"=== analyze_highgain: v8.root ("<<nEnt<<" eventos) ===\n\n";
     sf<<"FASE 0:\n  ph_medio hod2019X="<<std::fixed<<std::setprecision(2)<<ph19x_mean
       <<"  max="<<mx19x<<"  timing(ph>15)="<<h_pt19x->GetMean()<<" ns\n\n";
     sf<<"FASE 1 — Barrido:\n";
     sf<<std::setw(10)<<"abst19"<<std::setw(10)<<"N_clean"<<std::setw(10)<<"r_x"<<std::setw(10)<<"r_y\n";
     for(const auto&sr:sweep_results)
         sf<<std::setw(10)<<sr.thr<<std::setw(10)<<sr.n
           <<std::setw(10)<<std::setprecision(4)<<sr.rx<<std::setw(10)<<sr.ry<<"\n";
     sf<<"  abst19_opt="<<abst19_opt<<"  max|r|="<<best_r<<"\n\n";
     sf<<"FASE 2 — Mapeo:\n";
     sf<<"  τ(A,0) X="<<tau_AX<<"  Y="<<tau_AY<<"\n";
     sf<<"  τ(B,0) X="<<tau_BX<<"  Y="<<tau_BY<<"\n";
     sf<<"  τ(A,B) X="<<tau_ABX<<"  Y="<<tau_ABY<<"\n";
     sf<<"  r(π₀) X="<<r_pi0X<<"  r(π_A) X="<<r_piAX<<"\n\n";
     sf<<"FASE 3:\n";
     sf<<"  Cutflow: tot="<<n_total<<" 19c="<<cf_19c<<" 19m="<<cf_19m
       <<" 18c="<<cf_18c<<" xy="<<cf_xy<<" trk="<<cf_trk<<"\n";
     sf<<"  r_x="<<rx<<"  r_y="<<ry<<"  slope_x="<<slx<<"  n="<<(long long)vx18.size()<<"\n\n";
     sf<<"FASE 4 — Comparativa v7 vs v8:\n";
     for(auto [m,v7,v8]:std::initializer_list<std::tuple<std::string,std::string,std::string>>{
         {"ph_medio hod2019X","4.4",std::to_string(ph19x_mean).substr(0,5)},
         {"abst19_opt","N/A",std::to_string((int)abst19_opt)},
         {"N_clean",std::to_string(best_sw.n>0?52LL:52LL),std::to_string(best_sw.n)},
         {"r_x","0.069",std::to_string(rx).substr(0,6)},
         {"τ_A_X","0.12",std::to_string(tau_AX).substr(0,6)}})
         sf<<std::setw(25)<<m<<std::setw(20)<<v7<<v8<<"\n";}

    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
