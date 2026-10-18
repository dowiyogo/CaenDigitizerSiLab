// analyze_independent.cxx
// Análisis independiente hod2018 y hod2019 — figuras de mérito para tesis.
// Procesa v7.root y v8.root en el mismo run; compara en todas las figuras.
//
// FASES:
//   1 — Espectros de carga: global (1A), por canal (1B), pedestal/ruido (1C)
//   2 — Eficiencia: vs umbral (2A), por canal (2B), mapa 2D (2C)
//   3 — Calibración espacial suave (sin coincidencia inter-plano)
//   4 — Tabla de figuras de mérito + figura resumen de 6 paneles
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o analyze_independent analyze_independent.cxx

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
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TROOT.h>
#include <TSpectrum.h>
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
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "hodo_paths.h"

namespace fs = std::filesystem;

//==============================================================================
// Constantes y mapeo
//==============================================================================
constexpr double kPre_ns = 2500.0;   // fin región baseline [ns]
constexpr int    kNCh    = 64;

// Planos: 0=hod2018X, 1=hod2018Y, 2=hod2019X, 3=hod2019Y
constexpr int P18X=0,P18Y=1,P19X=2,P19Y=3;

// Ventanas de señal [ns]
static const double kTlo[4]={3500,3500,3000,3000};
static const double kThi[4]={7000,7000,5000,5000};

// Mapeo canales → índice físico k (0..15)
// Hod2018: k→canal k (X) o 16+k (Y). Identidad.
// Hod2019: tablas literales
static const std::array<int,16> kH19X={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};

// Asigna plano e índice físico a cada canal (−1 si no pertenece)
static int ch_plane(int ch) {
    if(ch>=0&&ch<=15) return P18X;
    if(ch>=16&&ch<=31) return P18Y;
    if(ch>=32&&ch<=47) return P19Y;
    if(ch>=48&&ch<=63) return P19X;
    return -1;
}
static int ch_k(int ch) {
    if(ch>=0&&ch<=15) return ch;
    if(ch>=16&&ch<=31) return ch-16;
    for(int k=0;k<16;k++) if(kH19Y[k]==ch) return k;
    for(int k=0;k<16;k++) if(kH19X[k]==ch) return k;
    return -1;
}

//==============================================================================
// Extracción de pulso (igual a versiones anteriores)
//==============================================================================
struct PR {bool valid=false;double ped=0,noise=0,ph=0,pt=-1;};
static PR extr(const std::vector<float>&wf,const std::vector<float>&t,
               double tlo,double thi) {
    PR r; const int n=static_cast<int>(wf.size());
    if(n<4||static_cast<int>(t.size())!=n) return r;
    double s=0,s2=0; int nb=0;
    for(int i=0;i<n;++i) if(t[i]<kPre_ns){s+=wf[i];s2+=wf[i]*wf[i];++nb;}
    if(nb<2) return r;
    r.valid=true; r.ped=s/nb;
    const double v=s2/nb-r.ped*r.ped; r.noise=(v>0)?std::sqrt(v):0.0;
    double vmin=1e30; int imin=-1;
    for(int i=0;i<n;++i) if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0) return r;
    r.ph=r.ped-vmin; r.pt=t[imin];
    return r;
}

//==============================================================================
// Clustering (gap ≤ 1)
//==============================================================================
struct CHit{int k;double ph,pt;};
struct Cl{double pos,time,amp;int ns,lo,hi;};
static std::vector<Cl> make_cls(std::vector<CHit> hits,double pitch,double thr){
    std::vector<CHit> act;
    for(const auto&h:hits) if(h.ph>=thr) act.push_back(h);
    if(act.empty()) return{};
    std::sort(act.begin(),act.end(),[](const CHit&a,const CHit&b){return a.k<b.k;});
    std::vector<Cl> out; Cl cur{};
    auto flush=[&](){cur.pos/=cur.amp;cur.time/=cur.amp;out.push_back(cur);};
    auto open=[&](const CHit&h){cur={};cur.lo=cur.hi=h.k;};
    auto add=[&](const CHit&h){
        cur.pos+=h.ph*(h.k+0.5-8.0);cur.time+=h.ph*h.pt;
        cur.amp+=h.ph;++cur.ns;cur.hi=h.k;};
    open(act[0]);add(act[0]);
    for(int i=1;i<(int)act.size();++i){if(act[i].k>cur.hi+1){flush();open(act[i]);}add(act[i]);}
    flush();
    for(auto&c:out)c.pos*=pitch;
    return out;
}
static const Cl& bamp(const std::vector<Cl>&v){
    return *std::max_element(v.begin(),v.end(),[](const Cl&a,const Cl&b){return a.amp<b.amp;});}

//==============================================================================
// Estadísticas
//==============================================================================
static double pearson_r(const std::vector<double>&x,const std::vector<double>&y){
    const int n=static_cast<int>(x.size());if(n<2)return 0;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];syy+=y[i]*y[i];sxy+=x[i]*y[i];}
    const double d=n,num=d*sxy-sx*sy,den=std::sqrt((d*sxx-sx*sx)*(d*syy-sy*sy));
    return den>0?num/den:0.0;}
static void linfit(const std::vector<double>&x,const std::vector<double>&y,double&sl,double&ic){
    sl=ic=0;const int n=static_cast<int>(x.size());if(n<2)return;
    double sx=0,sy=0,sxx=0,sxy=0;
    for(int i=0;i<n;++i){sx+=x[i];sy+=y[i];sxx+=x[i]*x[i];sxy+=x[i]*y[i];}
    const double d=n,dn=d*sxx-sx*sx;if(std::abs(dn)<1e-12)return;
    sl=(d*sxy-sx*sy)/dn;ic=(sy-sl*sx)/d;}
static double kendall_tau(const std::vector<int>&a,const std::vector<int>&b){
    const int n=static_cast<int>(a.size());long long C=0,D=0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
        int da=a[i]-a[j],db=b[i]-b[j];
        if(da*db>0)++C;else if(da*db<0)++D;}
    const long long tot=(long long)n*(n-1)/2;
    return tot>0?(double)(C-D)/tot:0.0;}

//==============================================================================
// Búsqueda de valle (mínimo local entre noise peak y MIP peak)
//==============================================================================
static double find_valley(TH1D*h,double x_lo,double x_hi){
    int b0=h->FindBin(x_lo),b1=h->FindBin(x_hi);
    int bv=b0;double vv=1e30;
    for(int b=b0;b<=b1;++b){double c=h->GetBinContent(b);if(c>0&&c<vv){vv=c;bv=b;}}
    return h->GetBinCenter(bv);}

// Moda del histograma en un rango
static double find_mode(TH1D*h,double x_lo,double x_hi){
    int b0=h->FindBin(x_lo),b1=h->FindBin(x_hi);
    int bm=b0;double vm=-1;
    for(int b=b0;b<=b1;++b){double c=h->GetBinContent(b);if(c>vm){vm=c;bm=b;}}
    return h->GetBinCenter(bm);}

//==============================================================================
// Ajuste Landau
//==============================================================================
struct LandauFit {double mpv=0,sigma=0,chi2ndf=0;bool ok=false;};
static LandauFit fit_landau(TH1D*h,double xlo,double xhi){
    LandauFit res;
    if(h->Integral(h->FindBin(xlo),h->FindBin(xhi))<20) return res;
    const double mpv0=find_mode(h,xlo,xhi);
    const double sig0=std::max(0.1*mpv0,10.0);
    TF1 fL("fLtmp","TMath::Landau(x,[0],[1],0)",xlo,std::min(xhi,800.0));
    fL.SetParameters(mpv0,sig0);
    fL.SetParLimits(0,xlo,xhi);
    fL.SetParLimits(1,1.0,mpv0);
    int stat=h->Fit(&fL,"RQN");  // R=range, Q=quiet, N=no draw
    if(stat==0||stat==4000){  // 4000 = fit ok with warnings
        res.ok=true;
        res.mpv=fL.GetParameter(0);
        res.sigma=fL.GetParameter(1);
        res.chi2ndf=fL.GetChisquare()/(fL.GetNDF()>0?fL.GetNDF():1);
    }
    return res;
}

//==============================================================================
// Detección de picos SPE
//==============================================================================
struct SPEResult {
    bool found=false;
    double delta=0;     // separación entre picos (= ganancia MPPC en ADC)
    int n_peaks=0;
    std::vector<double> positions;
};
static SPEResult detect_spe(TH1D*h){
    SPEResult res;
    // Solo mirar en rango bajo [0, 120] ADC donde pueden aparecer picos SPE
    TH1D* hsub=static_cast<TH1D*>(h->Clone("_spe_sub"));
    hsub->GetXaxis()->SetRangeUser(0,120);
    TSpectrum spec(10);
    int np=spec.Search(hsub,2.0,"nobackground nodraw",0.05);
    if(np<2){delete hsub;return res;}
    const double* xpos=spec.GetPositionX();
    std::vector<double> peaks(xpos,xpos+np);
    std::sort(peaks.begin(),peaks.end());
    // Verifica equidistancia
    std::vector<double> gaps;
    for(int i=1;i<(int)peaks.size();++i)gaps.push_back(peaks[i]-peaks[i-1]);
    if(gaps.empty()){delete hsub;return res;}
    const double mean_gap=std::accumulate(gaps.begin(),gaps.end(),0.0)/gaps.size();
    double max_dev=0;
    for(double g:gaps)max_dev=std::max(max_dev,std::abs(g-mean_gap));
    if(max_dev<0.25*mean_gap&&mean_gap>5){  // equidistante dentro del 25%
        res.found=true;
        res.delta=mean_gap;
        res.n_peaks=np;
        res.positions=peaks;
    }
    delete hsub;
    return res;
}

//==============================================================================
// Datos acumulados por archivo
//==============================================================================
struct FileAna {
    std::string label;          // "v7" o "v8"
    Long64_t n_events=0;
    double pitch=2.0;

    // 1A: Espectros globales por plano [0..3]
    TH1D* h_ph[4]={};

    // 1B: Espectro por canal
    TH1D* h_ph_ch[kNCh]={};

    // 1C: Pedestal y noise por canal
    double ped[kNCh]={}, noise[kNCh]={};
    long long n_ped[kNCh]={};

    // 2A: Eficiencia vs umbral
    std::vector<double> thr;         // umbrales
    long long eff_n[4][12]={};       // [plano][thr_idx]

    // 2B: Ocupancia por canal al umbral natural
    double umbral_nat[4]={};         // umbral natural por plano (aplicar a hod2019)
    double umbral_signal=0;          // para Phase 3 (señal real hod2019, puede ≠ umbral_nat)
    long long occ_nat[4][16]={};     // [plano][k]: hits por canal
    long long eff_nat_events[4]={};  // [plano]: eventos donde ≥1 canal > umbral_nat

    // 2C: Mapas de impacto 2D
    TH2D* h_imp18=nullptr, *h_imp19=nullptr;
    double thr_cls_18=5;             // umbral para clustering hod2018
    double thr_cls_19=38.5;         // umbral para clustering hod2019 en 2C

    // 3: Calibración espacial hod2019 X (sin requerir Y simultáneo)
    long long nc_X[16]={};           // eventos usados por canal hod2019X
    double sx_X[16]={}, sx2_X[16]={}; // suma x2018, suma x2018²
    // Análogo Y
    long long nc_Y[16]={};
    double sy_Y[16]={}, sy2_Y[16]={};

    // Derivados de Phase 3
    double xhat[16]={}, xerr[16]={};
    double yhat[16]={}, yerr[16]={};
    std::vector<int> ordA_X,ordA_Y;
    double tau_AX=0,tau_AY=0;
    double r_pi0_X=0,r_piA_X=0;

    // Ajuste Landau [plano]
    LandauFit landau[4];

    // SPE hod2019X
    SPEResult spe19x,spe19y;
};

//==============================================================================
// Procesamiento de un archivo
//==============================================================================
static void process_file(FileAna& fa, const std::string& fname){
    const std::string resolved = hodo::resolveRawFile(fname);
    TFile* fin=TFile::Open(resolved.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: "<<resolved<<'\n';return;}
    TTree* tree=static_cast<TTree*>(fin->Get("data_dig"));
    if(!tree){std::cerr<<"ERROR: no data_dig in "<<fname<<'\n';fin->Close();return;}

    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw=nullptr;
    tree->SetBranchStatus("*",0);tree->SetCacheSize(256LL<<20);
    for(int c=0;c<kNCh;++c){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!tree->GetBranch(bn))continue;
        tree->SetBranchStatus(bn,1);tree->SetBranchAddress(bn,&cw[c]);
    }
    tree->SetBranchStatus("time",1);tree->SetBranchAddress("time",&tw);

    const Long64_t nEnt=tree->GetEntries();
    fa.n_events=nEnt;
    std::cout<<"  ["<<fa.label<<"] "<<nEnt<<" eventos...\n";

    // Número de umbrales
    const int nT=static_cast<int>(fa.thr.size());

    for(Long64_t iev=0;iev<nEnt;++iev){
        tree->GetEntry(iev);
        if(iev%30000==0)std::cout<<"    "<<iev<<"/"<<nEnt<<"\r"<<std::flush;
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;

        // Extrae pulsos de los 64 canales
        std::array<PR,kNCh> pr{};
        for(int c=0;c<kNCh;++c){
            if(!cw[c]||cw[c]->empty())continue;
            const int pl=ch_plane(c);if(pl<0)continue;
            pr[c]=extr(*cw[c],tv,kTlo[pl],kThi[pl]);
        }

        // Máximo ph por plano (para eficiencia)
        double max_ph[4]={};
        for(int c=0;c<kNCh;++c){
            if(!pr[c].valid)continue;
            const int pl=ch_plane(c);if(pl<0)continue;
            max_ph[pl]=std::max(max_ph[pl],pr[c].ph);
        }

        // Hits por plano para clustering
        std::array<std::vector<CHit>,4> hits;
        for(int c=0;c<kNCh;++c){
            if(!pr[c].valid)continue;
            const int pl=ch_plane(c);if(pl<0)continue;
            const int k=ch_k(c);if(k<0)continue;
            // Histogramas
            if(pr[c].ph>0){
                fa.h_ph[pl]->Fill(pr[c].ph);
                fa.h_ph_ch[c]->Fill(pr[c].ph);
            }
            // Pedestal y ruido
            fa.ped[c]+=pr[c].ped; fa.noise[c]+=pr[c].noise; fa.n_ped[c]++;
            // Para hits y clustering
            hits[pl].push_back({k,pr[c].ph,pr[c].pt});
        }

        // 2A: eficiencia
        for(int pl=0;pl<4;++pl)
            for(int ti=0;ti<nT;++ti)
                if(max_ph[pl]>fa.thr[ti]) fa.eff_n[pl][ti]++;

        // 2B: ocupancia al umbral natural + conteo de EVENTOS (no hits)
        for(int pl=0;pl<4;++pl){
            bool fired=false;
            for(const auto&h:hits[pl])
                if(h.ph>fa.umbral_nat[pl]){fa.occ_nat[pl][h.k]++;fired=true;}
            if(fired) fa.eff_nat_events[pl]++;
        }

        // 2C: mapa 2D hod2018 (1 cluster X + 1 cluster Y al umbral natural 20 ADC)
        // Usando umbral más alto para reducir multiplicidad y obtener eventos con 1 cluster
        {
            const double thr_imp=fa.umbral_nat[P18X];  // 20 ADC
            auto c18x=make_cls(hits[P18X],fa.pitch,thr_imp);
            auto c18y=make_cls(hits[P18Y],fa.pitch,thr_imp);
            // Aceptar ≤2 clusters (más flexible) y usar el leading por amplitud
            if(!c18x.empty()&&!c18y.empty())
                fa.h_imp18->Fill(bamp(c18x).pos,bamp(c18y).pos);
        }
        // 2C: mapa 2D hod2019 (1 cluster X + 1 cluster Y al umbral natural)
        {
            auto c19x=make_cls(hits[P19X],fa.pitch,fa.umbral_nat[P19X]);
            auto c19y=make_cls(hits[P19Y],fa.pitch,fa.umbral_nat[P19Y]);
            if(c19x.size()==1&&c19y.size()==1)
                fa.h_imp19->Fill(bamp(c19x).pos,bamp(c19y).pos);
        }

        // Phase 3A: calibración X (sin requerir coincidencia X-Y de hod2019)
        // Selección: exactamente 1 canal hod2019X activo (ph ≥ umbral_signal)
        //            + cluster líder de hod2018X con amplitud alta (ph_lead ≥ 40 ADC)
        {
            std::vector<CHit> act19x;
            for(const auto&h:hits[P19X]) if(h.ph>=fa.umbral_signal) act19x.push_back(h);
            if(act19x.size()==1){
                // Cluster líder hod2018X (máxima amplitud, independiente de cuántos clusters haya)
                auto c18x=make_cls(hits[P18X],fa.pitch,fa.thr_cls_18);
                if(!c18x.empty()){
                    const double x18=bamp(c18x).pos;
                    const int ki=act19x[0].k;
                    fa.nc_X[ki]++;fa.sx_X[ki]+=x18;fa.sx2_X[ki]+=x18*x18;
                }
            }
        }
        // Phase 3B: calibración Y
        {
            std::vector<CHit> act19y;
            for(const auto&h:hits[P19Y]) if(h.ph>=fa.umbral_signal) act19y.push_back(h);
            if(act19y.size()==1){
                auto c18y=make_cls(hits[P18Y],fa.pitch,fa.thr_cls_18);
                if(!c18y.empty()){
                    const double y18=bamp(c18y).pos;
                    const int ki=act19y[0].k;
                    fa.nc_Y[ki]++;fa.sy_Y[ki]+=y18;fa.sy2_Y[ki]+=y18*y18;
                }
            }
        }
    }
    std::cout<<"    "<<nEnt<<"/"<<nEnt<<" ✓\n";
    fin->Close();

    // Normalizar pedestal/ruido
    for(int c=0;c<kNCh;++c) if(fa.n_ped[c]>0){
        fa.ped[c]/=fa.n_ped[c]; fa.noise[c]/=fa.n_ped[c];}
}

//==============================================================================
// Inferencia de mapeo (después de process_file)
//==============================================================================
static void compute_mapping(FileAna& fa){
    // x̂_c para hod2019X
    for(int c=0;c<16;++c){
        if(fa.nc_X[c]>0){
            fa.xhat[c]=fa.sx_X[c]/fa.nc_X[c];
            const double v=fa.sx2_X[c]/fa.nc_X[c]-fa.xhat[c]*fa.xhat[c];
            fa.xerr[c]=(v>0)?std::sqrt(v)/std::sqrt((double)fa.nc_X[c]):1e9;
        } else {fa.xhat[c]=1e9;fa.xerr[c]=1e9;}
    }
    for(int c=0;c<16;++c){
        if(fa.nc_Y[c]>0){
            fa.yhat[c]=fa.sy_Y[c]/fa.nc_Y[c];
            const double v=fa.sy2_Y[c]/fa.nc_Y[c]-fa.yhat[c]*fa.yhat[c];
            fa.yerr[c]=(v>0)?std::sqrt(v)/std::sqrt((double)fa.nc_Y[c]):1e9;
        } else {fa.yhat[c]=1e9;fa.yerr[c]=1e9;}
    }

    // Ordenamiento π_A por x̂_c / ŷ_c
    fa.ordA_X.resize(16);fa.ordA_Y.resize(16);
    std::iota(fa.ordA_X.begin(),fa.ordA_X.end(),0);
    std::iota(fa.ordA_Y.begin(),fa.ordA_Y.end(),0);
    std::stable_sort(fa.ordA_X.begin(),fa.ordA_X.end(),[&](int a,int b){
        if(fa.nc_X[a]==0&&fa.nc_X[b]==0)return a<b;
        if(fa.nc_X[a]==0)return false;if(fa.nc_X[b]==0)return true;
        return fa.xhat[a]<fa.xhat[b];});
    std::stable_sort(fa.ordA_Y.begin(),fa.ordA_Y.end(),[&](int a,int b){
        if(fa.nc_Y[a]==0&&fa.nc_Y[b]==0)return a<b;
        if(fa.nc_Y[a]==0)return false;if(fa.nc_Y[b]==0)return true;
        return fa.yhat[a]<fa.yhat[b];});

    // Kendall τ
    std::vector<int> rA(16,-1),r0(16,-1);
    for(int k=0;k<16;++k){rA[fa.ordA_X[k]]=k;r0[kH19X[k]-48]=k;}
    fa.tau_AX=kendall_tau(rA,r0);
    for(int k=0;k<16;++k){rA[fa.ordA_Y[k]]=k;r0[kH19Y[k]-32]=k;}
    fa.tau_AY=kendall_tau(rA,r0);
}

//==============================================================================
// Plot helpers
//==============================================================================
static const int kColors[2]={kBlue+1,kOrange+7};  // v7=azul, v8=naranja

static void lab(TH1D*h,const char*x,const char*y){
    h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    h->GetXaxis()->SetTitleSize(0.05);h->GetYaxis()->SetTitleSize(0.05);}
static void lab2(TH2D*h,const char*x,const char*y){
    h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    h->GetXaxis()->SetTitleSize(0.05);h->GetYaxis()->SetTitleSize(0.05);}

static void draw_vline(double x,int col=kRed,int style=2){
    double ymin=gPad->GetUymin(),ymax=gPad->GetUymax();
    TLine* l=new TLine(x,gPad->GetLogy()?std::pow(10,ymin):ymin,
                        x,gPad->GetLogy()?std::pow(10,ymax*0.95):ymax*0.95);
    l->SetLineColor(col);l->SetLineStyle(style);l->SetLineWidth(2);l->Draw();}

//==============================================================================
// main
//==============================================================================
int main(int argc,char**argv){
    std::string fv7="lectura_vx1740_ext_nim_1M_v7.root";
    std::string fv8="lectura_vx1740_ext_nim_1M_v8.root";
    for(int i=1;i<argc;++i){
        std::string a=argv[i];
        if(a=="--v7"&&i+1<argc){fv7=argv[++i];continue;}
        if(a=="--v8"&&i+1<argc){fv8=argv[++i];continue;}
    }

    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("independent_analysis_"+std::string(ds));
    for(auto s:{"spectra","efficiency","maps","mapping","summary"})
        fs::create_directories(outdir/s);

    std::cout<<"=== analyze_independent ===\n";
    std::cout<<"v7: "<<fv7<<"\nv8: "<<fv8<<"\nOutput: "<<outdir<<"\n\n";

    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    static const char* pname[4]={"hod2018X","hod2018Y","hod2019X","hod2019Y"};

    //==========================================================================
    // Inicializa FileAna para v7 y v8
    //==========================================================================
    FileAna FA[2];
    FA[0].label="v7"; FA[1].label="v8";

    // Umbrales barrido
    FA[0].thr={3,5,8,12,17,25};
    FA[1].thr={3,5,10,20,38,50,75,100};

    // Umbrales fijos para cada archivo:
    // hod2018: umbral absoluto del análisis (20 ADC — ruido gaussiano << MIP)
    // hod2019: calculado del espectro (valle entre noise-crossings y señal real)
    FA[0].umbral_nat[P18X]=20;FA[0].umbral_nat[P18Y]=20;
    FA[0].umbral_nat[P19X]=17;FA[0].umbral_nat[P19Y]=17;  // señal real v7 ≈ 17 ADC
    FA[0].umbral_signal=17;
    FA[0].thr_cls_18=5;FA[0].thr_cls_19=17;

    // v8: hod2019 valle medido ~38.5 ADC; hod2018 threshold 20 ADC
    FA[1].umbral_nat[P18X]=20;FA[1].umbral_nat[P18Y]=20;
    FA[1].umbral_nat[P19X]=38.5;FA[1].umbral_nat[P19Y]=38.5;
    // Calibración Fase 3: umbral de señal para hod2019 que da estadística suficiente.
    // Usamos 5 ADC (barrido de Phase 1 mostró N=331 eventos a ese umbral).
    FA[1].umbral_signal=5.0;
    FA[1].thr_cls_18=5;FA[1].thr_cls_19=38.5;

    // Inicializa histogramas
    for(int v=0;v<2;++v){
        FileAna&fa=FA[v];
        const std::string lbl=fa.label+"_";
        for(int pl=0;pl<4;++pl){
            fa.h_ph[pl]=new TH1D((lbl+pname[pl]).c_str(),
                (std::string(pname[pl])+" "+fa.label+";Pulse height [ADC];Entries").c_str(),
                800,0,800);
        }
        for(int c=0;c<kNCh;++c){
            char bn[32];std::snprintf(bn,sizeof(bn),"%sch%d",lbl.c_str(),c);
            fa.h_ph_ch[c]=new TH1D(bn,bn,800,0,800);
        }
        fa.h_imp18=new TH2D((lbl+"imp18").c_str(),
            ("Mapa impacto hod2018 "+fa.label+";X_{2018} [mm];Y_{2018} [mm];Events").c_str(),
            17,-17,17,17,-17,17);
        fa.h_imp19=new TH2D((lbl+"imp19").c_str(),
            ("Mapa impacto hod2019 "+fa.label+";X_{2019} [mm];Y_{2019} [mm];Events").c_str(),
            17,-17,17,17,-17,17);
    }

    //==========================================================================
    // Pasada única por archivo
    //==========================================================================
    std::cout<<"[1/2] Procesando v7...\n";
    process_file(FA[0],fv7);
    std::cout<<"[2/2] Procesando v8...\n";
    process_file(FA[1],fv8);
    std::cout<<"Pasadas completas.\n\n";

    // Refina umbrales desde los espectros reales (hod2019 solamente; hod2018 ya fijo a 20 ADC)
    for(int v=0;v<2;++v){
        FileAna&fa=FA[v];
        // hod2019 X: busca valle en [8, 50] ADC. Si < 8 o no claro, usa preset.
        double v19x=find_valley(fa.h_ph[P19X],8,50);
        if(v19x<8||v19x>50)v19x=(v==0?17.0:38.5);
        double v19y=find_valley(fa.h_ph[P19Y],8,50);
        if(v19y<8||v19y>50)v19y=v19x;
        fa.umbral_nat[P19X]=v19x; fa.umbral_nat[P19Y]=v19y;
        // umbral_signal para Phase 3 ya fue configurado arriba (no sobreescribir)
    }

    // Ajuste Landau para cada plano y archivo
    for(int v=0;v<2;++v){
        FileAna&fa=FA[v];
        for(int pl=0;pl<4;++pl){
            const double vl=fa.umbral_nat[pl];
            fa.landau[pl]=fit_landau(fa.h_ph[pl],vl,700);
        }
        // Detección SPE hod2019X en v8
        if(v==1){
            fa.spe19x=detect_spe(fa.h_ph[P19X]);
            fa.spe19y=detect_spe(fa.h_ph[P19Y]);
        }
    }

    // Calcula mapeos inferidos
    for(int v=0;v<2;++v) compute_mapping(FA[v]);

    //==========================================================================
    // FASE 1A — Espectros de carga (4 planos × 2 archivos)
    //==========================================================================
    std::cout<<"Generando figuras...\n";
    {
        TCanvas cc("c_spectra","",1600,1200);
        cc.Divide(2,2);
        for(int pl=0;pl<4;++pl){
            cc.cd(pl+1);
            gPad->SetLogy();gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);
            const double ymax=std::max(FA[0].h_ph[pl]->GetMaximum(),
                                        FA[1].h_ph[pl]->GetMaximum())*2;
            // Normaliza a área unitaria
            TH1D* h0=(TH1D*)FA[0].h_ph[pl]->Clone("_n0");
            TH1D* h1=(TH1D*)FA[1].h_ph[pl]->Clone("_n1");
            if(h0->Integral()>0)h0->Scale(1.0/h0->Integral());
            if(h1->Integral()>0)h1->Scale(1.0/h1->Integral());
            lab(h0,"Pulse height [ADC]","Density [1/ADC]");
            h0->SetLineColor(kColors[0]);h0->SetLineWidth(2);h0->SetStats(0);
            h1->SetLineColor(kColors[1]);h1->SetLineWidth(2);h1->SetStats(0);
            h0->SetTitle(pname[pl]);
            h0->GetXaxis()->SetRangeUser(0,600);
            h0->Draw("HIST");h1->Draw("HIST SAME");
            // Líneas umbral natural
            draw_vline(FA[0].umbral_nat[pl],kColors[0],2);
            draw_vline(FA[1].umbral_nat[pl],kColors[1],3);
            // Landau fit v8
            if(FA[1].landau[pl].ok){
                TF1 fL2("fLv8","TMath::Landau(x,[0],[1],0)",
                    FA[1].umbral_nat[pl],700);
                fL2.SetParameters(FA[1].landau[pl].mpv,FA[1].landau[pl].sigma);
                fL2.SetLineColor(kColors[1]);fL2.SetLineStyle(1);fL2.SetLineWidth(1);
                // Scale to match histogram
                const double scale=h1->Integral(h1->FindBin(FA[1].umbral_nat[pl]),800)/fL2.Integral(FA[1].umbral_nat[pl],700);
                fL2.SetParameter(0,FA[1].landau[pl].mpv);
                fL2.Draw("SAME");
            }
            TLegend lg(0.55,0.75,0.92,0.90);lg.SetBorderSize(0);lg.SetFillStyle(0);
            lg.AddEntry(h0,"v7 (17k)","l");lg.AddEntry(h1,"v8 (143k)","l");
            lg.Draw();
        }
        cc.SaveAs((outdir/"spectra"/"charge_spectra_all_planes.png").string().c_str());
    }

    //==========================================================================
    // FASE 1B — Espectro por canal (16 subpaneles para cada plano)
    //==========================================================================
    auto make_per_channel_canvas=[&](FileAna&fa,int pl,const fs::path&path){
        const int base=(pl==P18X)?0:(pl==P18Y)?16:(pl==P19Y)?32:48;
        const std::array<int,16>* prop=(pl==P19X)?&kH19X:(pl==P19Y)?&kH19Y:nullptr;
        TCanvas cc(("cch_"+fa.label+"_"+pname[pl]).c_str(),"",1800,1200);
        cc.Divide(4,4);
        double max_ph_global=0;
        for(int k=0;k<16;++k){
            int ch=(prop)?(*prop)[k]:base+k;
            max_ph_global=std::max(max_ph_global,fa.h_ph_ch[ch]->GetMaximum());}
        for(int k=0;k<16;++k){
            int ch=(prop)?(*prop)[k]:base+k;
            cc.cd(k+1);gPad->SetLogy();gPad->SetLeftMargin(0.15);gPad->SetBottomMargin(0.15);
            TH1D*h=fa.h_ph_ch[ch];
            h->SetLineColor(kColors[fa.label=="v7"?0:1]);h->SetLineWidth(1);
            h->GetXaxis()->SetRangeUser(0,400);
            h->SetStats(0);
            char tt[32];std::snprintf(tt,sizeof(tt),"k=%d Ch%d",k,ch);
            h->SetTitle(tt);h->Draw("HIST");
            draw_vline(fa.umbral_nat[pl],kRed,2);
        }
        cc.SaveAs(path.string().c_str());
    };
    make_per_channel_canvas(FA[1],P19X,outdir/"spectra"/"charge_per_channel_hod2019X_v8.png");
    make_per_channel_canvas(FA[1],P19Y,outdir/"spectra"/"charge_per_channel_hod2019Y_v8.png");
    make_per_channel_canvas(FA[1],P18X,outdir/"spectra"/"charge_per_channel_hod2018X_v8.png");
    make_per_channel_canvas(FA[1],P18Y,outdir/"spectra"/"charge_per_channel_hod2018Y_v8.png");

    //==========================================================================
    // FASE 1C — Ruido por canal
    //==========================================================================
    {
        TCanvas cc("c_noise","",1400,1000);
        cc.Divide(2,2);
        for(int pl=0;pl<4;++pl){
            cc.cd(pl+1);gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);
            const int base=(pl==P18X)?0:(pl==P18Y)?16:(pl==P19Y)?32:48;
            const std::array<int,16>* prop=(pl==P19X)?&kH19X:(pl==P19Y)?&kH19Y:nullptr;
            std::vector<double> n0(16),n1(16),xerr(16,0.3);
            for(int k=0;k<16;++k){
                int ch=(prop)?(*prop)[k]:base+k;
                n0[k]=FA[0].noise[ch]; n1[k]=FA[1].noise[ch];}
            std::vector<double> ks(16);std::iota(ks.begin(),ks.end(),0);
            TGraphErrors g0(16,ks.data(),n0.data(),xerr.data(),std::vector<double>(16,0).data());
            TGraphErrors g1(16,ks.data(),n1.data(),xerr.data(),std::vector<double>(16,0).data());
            g0.SetMarkerColor(kColors[0]);g0.SetLineColor(kColors[0]);g0.SetMarkerStyle(20);g0.SetMarkerSize(0.7);
            g1.SetMarkerColor(kColors[1]);g1.SetLineColor(kColors[1]);g1.SetMarkerStyle(21);g1.SetMarkerSize(0.7);
            g0.SetTitle((std::string(pname[pl])+" noise RMS;Strip k;Noise RMS [ADC]").c_str());
            g0.Draw("AP");g1.Draw("P SAME");
        }
        cc.SaveAs((outdir/"summary"/"noise_comparison_by_channel.png").string().c_str());
    }

    //==========================================================================
    // FASE 2A — Eficiencia vs umbral
    //==========================================================================
    {
        TCanvas cc("c_eff","",1600,700);
        cc.Divide(2,1);
        const int cols_pl[4]={kBlue+1,kBlue-4,kRed+1,kRed-4};
        for(int v=0;v<2;++v){
            cc.cd(v+1);gPad->SetLogx();gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);
            FileAna&fa=FA[v];
            const int nT=static_cast<int>(fa.thr.size());
            // Frame
            TH1D* frame=new TH1D(("feff"+fa.label).c_str(),
                ("Eficiencia vs umbral "+fa.label+";Umbral [ADC];#varepsilon").c_str(),
                100,1,200);
            frame->SetMinimum(0);frame->SetMaximum(1.05);frame->SetStats(0);frame->Draw();
            for(int pl=0;pl<4;++pl){
                std::vector<double> x,y,ey;
                for(int ti=0;ti<nT;++ti){
                    double eps=(double)fa.eff_n[pl][ti]/fa.n_events;
                    double err=std::sqrt(eps*(1-eps)/fa.n_events);
                    x.push_back(fa.thr[ti]);y.push_back(eps);ey.push_back(err);
                }
                TGraphErrors* g=new TGraphErrors(nT,x.data(),y.data(),
                    std::vector<double>(nT,0).data(),ey.data());
                g->SetLineColor(cols_pl[pl]);g->SetMarkerColor(cols_pl[pl]);
                g->SetMarkerStyle(20+pl);g->SetLineWidth(2);
                g->Draw("PL SAME");
            }
            // Umbral natural v8 hod2019
            double vn19=fa.umbral_nat[P19X];
            draw_vline(vn19,kRed,3);
            TLegend lg(0.45,0.60,0.88,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
            for(int pl=0;pl<4;++pl){
                const double eps=(double)fa.eff_n[pl][0]/fa.n_events;
                const std::string lbl=std::string(pname[pl])+" (e@3ADC="+
                     std::to_string((int)(eps*100))+"%)";
                lg.AddEntry((TObject*)nullptr, lbl.c_str(),"");
            }
            lg.Draw();
            TLatex lat;lat.SetNDC();lat.SetTextSize(0.035);
            lat.DrawLatex(0.14,0.92,(std::string("Eficiencia ")+fa.label).c_str());
        }
        cc.SaveAs((outdir/"efficiency"/"efficiency_vs_threshold_v7_v8.png").string().c_str());
    }

    //==========================================================================
    // FASE 2B — Eficiencia por canal (chi-squared test)
    //==========================================================================
    {
        TCanvas cc("c_effch","",1400,1000);
        cc.Divide(2,2);
        for(int pl=0;pl<4;++pl){
            cc.cd(pl+1);gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.13);
            // Usa v8
            FileAna&fa=FA[1];
            std::vector<double> ks(16),eff(16),eeff(16);
            double sum_occ=0;
            for(int k=0;k<16;++k) sum_occ+=fa.occ_nat[pl][k];
            const double mean_occ=sum_occ/16.0;  // expected hits per channel if uniform
            double chi2=0;int ndof=0;
            for(int k=0;k<16;++k){
                ks[k]=k;
                // ε_c = hits per channel / n_events (puede sumar > 1 si hay multi-hit)
                eff[k]=(double)fa.occ_nat[pl][k]/fa.n_events;
                eeff[k]=std::sqrt(std::max(eff[k],1e-6)*(1+eff[k])/fa.n_events);
                if(mean_occ>0){
                    chi2+=std::pow((double)fa.occ_nat[pl][k]-mean_occ,2)/mean_occ;
                    ++ndof;}
            }
            const double pval=TMath::Prob(chi2,ndof>1?ndof-1:1);
            TGraphErrors g(16,ks.data(),eff.data(),std::vector<double>(16,0.3).data(),eeff.data());
            g.SetTitle((std::string(pname[pl])+" v8;Strip k;#varepsilon per channel").c_str());
            g.SetMarkerColor(kAzure+2);g.SetLineColor(kAzure+2);g.SetMarkerStyle(20);
            g.GetYaxis()->SetRangeUser(0,eff[0]*3+0.01);
            g.Draw("AP");
            const double emean=sum_occ/16.0/fa.n_events;
            TLine lm(-0.5,emean,15.5,emean);
            lm.SetLineColor(kRed);lm.SetLineWidth(2);lm.SetLineStyle(2);lm.Draw();
            TLatex lat;lat.SetNDC();lat.SetTextSize(0.038);
            char buf[80];std::snprintf(buf,sizeof(buf),"#chi^{2}/dof=%.1f  p=%.3f",chi2/ndof,pval);
            lat.DrawLatex(0.30,0.90,buf);
        }
        cc.SaveAs((outdir/"efficiency"/"efficiency_per_channel_all_planes.png").string().c_str());
    }

    //==========================================================================
    // FASE 2C — Mapas de impacto 2D
    //==========================================================================
    auto save2d=[&](TH2D*h,const fs::path&p,const std::string& ti){
        TCanvas c2("c2","",900,820);c2.SetLeftMargin(0.13);c2.SetRightMargin(0.17);c2.SetBottomMargin(0.13);
        h->SetTitle(ti.c_str());lab2(h,h->GetXaxis()->GetTitle(),h->GetYaxis()->GetTitle());
        h->SetStats(0);h->Draw("COLZ");
        TLatex lat;lat.SetNDC();lat.SetTextSize(0.035);
        char buf[40];std::snprintf(buf,sizeof(buf),"N = %.0f",h->Integral());
        lat.DrawLatex(0.15,0.93,buf);
        c2.SaveAs(p.string().c_str());};

    save2d(FA[1].h_imp18,outdir/"maps"/"impact_map_hod2018_v8.png",
        "Mapa impacto hod2018 (v8, 1 cluster X+Y);X_{2018} [mm];Y_{2018} [mm]");
    if(FA[1].h_imp19->Integral()>=30)
        save2d(FA[1].h_imp19,outdir/"maps"/"impact_map_hod2019_v8.png",
            "Mapa impacto hod2019 (v8, 1 cluster X+Y);X_{2019} [mm];Y_{2019} [mm]");
    else{
        std::cout<<"  hod2019 mapa 2D: N="<<FA[1].h_imp19->Integral()
                 <<" < 30 → estadística insuficiente\n";}

    //==========================================================================
    // FASE 3 — Calibración espacial
    //==========================================================================
    auto make_calib_plot=[&](FileAna&fa,bool isX,const fs::path&path){
        const std::string ax=isX?"X":"Y";
        const int nn=16;
        const double* xhat=isX?fa.xhat:fa.yhat;
        const double* xerr=isX?fa.xerr:fa.yerr;
        const long long* nc=isX?fa.nc_X:fa.nc_Y;
        const std::array<int,16>& prop=isX?kH19X:kH19Y;
        const int base=isX?48:32;

        // x̂_c vs índice propuesto π₀(c)
        std::vector<double> xs,ys,exs,eys;
        for(int c=0;c<nn;++c){
            if(nc[c]<2)continue;
            int k0=-1;for(int k=0;k<16;++k)if(prop[k]==c+base){k0=k;break;}
            if(k0<0)continue;
            xs.push_back(k0);ys.push_back(xhat[c]);
            exs.push_back(0.25);eys.push_back(xerr[c]);
        }
        if(xs.empty()){std::cout<<"  No hay datos para calibración "<<ax<<"\n";return;}

        // Ajuste lineal
        double sl,ic;linfit(xs,ys,sl,ic);
        const double r=pearson_r(xs,ys);

        TGraphErrors g(static_cast<int>(xs.size()),xs.data(),ys.data(),exs.data(),eys.data());
        const std::string title="Calibración espacial hod2019 "+ax+" "+fa.label;
        g.SetTitle((title+";Índice propuesto #pi_0;#hat{x}_c [mm]").c_str());
        g.SetMarkerStyle(20);g.SetMarkerColor(kBlue+2);g.SetLineColor(kBlue+2);

        TCanvas cg("cg","",900,680);cg.SetLeftMargin(0.13);cg.SetBottomMargin(0.13);
        g.Draw("AP");
        // Línea de ajuste
        const double x0=g.GetXaxis()->GetXmin(),x1=g.GetXaxis()->GetXmax();
        TLine lfit(0,ic,15,sl*15+ic);
        lfit.SetLineColor(kRed);lfit.SetLineWidth(2);lfit.Draw("SAME");
        // Línea de referencia (mapeo perfecto: slope=pitch, intercept=0)
        TLine lref(0,(-7.5)*fa.pitch,15,(7.5)*fa.pitch);
        lref.SetLineColor(kGray+2);lref.SetLineWidth(2);lref.SetLineStyle(2);lref.Draw("SAME");
        TLatex lat;lat.SetNDC();lat.SetTextSize(0.035);
        char buf[128];
        std::snprintf(buf,sizeof(buf),"slope=%.2f mm  intercept=%.2f mm  r=%.3f",sl,ic,r);
        lat.DrawLatex(0.14,0.94,buf);
        std::snprintf(buf,sizeof(buf),"#tau_{{Kendall}} = %.3f",isX?fa.tau_AX:fa.tau_AY);
        lat.DrawLatex(0.14,0.90,buf);
        cg.SaveAs(path.string().c_str());
        if(isX){fa.r_pi0_X=r;}
    };

    for(int v=0;v<2;++v){
        const std::string lbl=FA[v].label;
        make_calib_plot(FA[v],true,outdir/"mapping"/("spatial_calibration_hod2019X_"+lbl+".png"));
        make_calib_plot(FA[v],false,outdir/"mapping"/("spatial_calibration_hod2019Y_"+lbl+".png"));
    }

    //==========================================================================
    // FASE 3C — Ocupancia π₀ vs π_A
    //==========================================================================
    {
        TCanvas cc("c_occ","",1000,680);
        cc.SetLeftMargin(0.12);cc.SetBottomMargin(0.12);
        FileAna&fa=FA[1];  // v8
        TH1D h0("occ_pi0","Ocupancia hod2019 X v8;Índice k;Hits",16,-0.5,15.5);
        TH1D hA("occ_piA","Ocupancia hod2019 X v8;Índice k;Hits",16,-0.5,15.5);
        for(int k=0;k<16;++k){
            int c0=kH19X[k]-48;  // π₀: canal en posición k
            h0.SetBinContent(k+1,fa.occ_nat[P19X][c0]);
            int cA=fa.ordA_X[k];  // π_A: canal en posición k
            hA.SetBinContent(k+1,fa.occ_nat[P19X][cA]);
        }
        // Rugosidad
        auto rg=[](TH1D*h)->double{double R=0;
            for(int k=1;k<16;k++){double d=h->GetBinContent(k+1)-h->GetBinContent(k);R+=d*d;}return R;};
        const double R0=rg(&h0),RA=rg(&hA);
        // Random baseline
        std::mt19937 rng(42);std::vector<int> perm(16);std::iota(perm.begin(),perm.end(),0);
        double Rrand=0;
        for(int t=0;t<200;++t){
            std::shuffle(perm.begin(),perm.end(),rng);
            TH1D htmp("","",16,-0.5,15.5);
            for(int k=0;k<16;++k)htmp.SetBinContent(k+1,fa.occ_nat[P19X][perm[k]]);
            Rrand+=rg(&htmp);}
        Rrand/=200;

        h0.SetLineColor(kRed);h0.SetLineWidth(2);h0.SetStats(0);
        hA.SetLineColor(kBlue+1);hA.SetLineWidth(2);hA.SetLineStyle(2);
        h0.SetTitle("Ocupancia hod2019 X (v8, umbral natural);Índice k;Hits");
        h0.Draw("HIST");hA.Draw("HIST SAME");
        TLegend lg(0.55,0.75,0.88,0.90);lg.SetBorderSize(0);lg.SetFillStyle(0);
        lg.AddEntry(&h0,"#pi_0 (propuesto)","l");
        lg.AddEntry(&hA,"#pi_A (inferido)","l");
        lg.Draw();
        TLatex lat;lat.SetNDC();lat.SetTextSize(0.033);
        char buf[120];
        std::snprintf(buf,sizeof(buf),"R(#pi_0)=%.0f  R(#pi_A)=%.0f  R(rand)=%.0f",R0,RA,Rrand);
        lat.DrawLatex(0.14,0.94,buf);
        cc.SaveAs((outdir/"mapping"/"occupancy_hod2019X_pi0_vs_piA.png").string().c_str());
    }

    //==========================================================================
    // FASE 4 — Figura resumen de 6 paneles (1800×1200)
    //==========================================================================
    {
        TCanvas cc("c_summary","",1800,1200);
        cc.Divide(3,2);

        // Panel 1: hod2018X (v7 vs v8)
        cc.cd(1);gPad->SetLogy();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
        {TH1D*h0=(TH1D*)FA[0].h_ph[P18X]->Clone("_p1_0");
         TH1D*h1=(TH1D*)FA[1].h_ph[P18X]->Clone("_p1_1");
         if(h0->Integral()>0)h0->Scale(1.0/h0->Integral());
         if(h1->Integral()>0)h1->Scale(1.0/h1->Integral());
         h0->SetLineColor(kColors[0]);h0->SetLineWidth(2);h0->SetStats(0);
         h1->SetLineColor(kColors[1]);h1->SetLineWidth(2);h1->SetStats(0);
         h0->GetXaxis()->SetRangeUser(0,600);
         h0->SetTitle("hod2018 X — v7 vs v8;PH [ADC];Density");
         h0->Draw("HIST");h1->Draw("HIST SAME");
         draw_vline(FA[1].umbral_nat[P18X],kColors[1],2);
         TLatex lat;lat.SetNDC();lat.SetTextSize(0.04);
         if(FA[1].landau[P18X].ok){
             char b[60];std::snprintf(b,sizeof(b),"MPV=%.0f (v8)",FA[1].landau[P18X].mpv);
             lat.SetTextColor(kColors[1]);lat.DrawLatex(0.35,0.85,b);}
         TLegend lg(0.55,0.78,0.90,0.92);lg.SetBorderSize(0);lg.SetFillStyle(0);
         lg.AddEntry(h0,"v7","l");lg.AddEntry(h1,"v8","l");lg.Draw();}

        // Panel 2: hod2019X (v7 vs v8) con SPE marcados
        cc.cd(2);gPad->SetLogy();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
        {TH1D*h0=(TH1D*)FA[0].h_ph[P19X]->Clone("_p2_0");
         TH1D*h1=(TH1D*)FA[1].h_ph[P19X]->Clone("_p2_1");
         if(h0->Integral()>0)h0->Scale(1.0/h0->Integral());
         if(h1->Integral()>0)h1->Scale(1.0/h1->Integral());
         h0->SetLineColor(kColors[0]);h0->SetLineWidth(2);h0->SetStats(0);
         h1->SetLineColor(kColors[1]);h1->SetLineWidth(2);h1->SetStats(0);
         h0->GetXaxis()->SetRangeUser(0,400);
         h0->SetTitle("hod2019 X — v7 vs v8;PH [ADC];Density");
         h0->Draw("HIST");h1->Draw("HIST SAME");
         draw_vline(FA[0].umbral_nat[P19X],kColors[0],2);
         draw_vline(FA[1].umbral_nat[P19X],kColors[1],3);
         // SPE markers para v8
         if(FA[1].spe19x.found){
             TLatex spe;spe.SetNDC(false);spe.SetTextSize(0.030);spe.SetTextColor(kViolet+1);
             for(double xp:FA[1].spe19x.positions)
                 spe.DrawLatex(xp,h1->GetMaximum()*0.3,"SPE");
         }
         TLatex lat;lat.SetNDC();lat.SetTextSize(0.038);
         if(FA[1].landau[P19X].ok){
             char b[60];std::snprintf(b,sizeof(b),"MPV=%.0f (v8)",FA[1].landau[P19X].mpv);
             lat.SetTextColor(kColors[1]);lat.DrawLatex(0.35,0.90,b);}
         if(FA[1].spe19x.found){
             char b[60];std::snprintf(b,sizeof(b),"SPE #Delta=%.1f",FA[1].spe19x.delta);
             lat.SetTextColor(kViolet+1);lat.DrawLatex(0.35,0.85,b);}
         TLegend lg(0.55,0.78,0.90,0.92);lg.SetBorderSize(0);lg.SetFillStyle(0);
         lg.AddEntry(h0,"v7","l");lg.AddEntry(h1,"v8","l");lg.Draw();}

        // Panel 3: Eficiencia vs umbral (4 planos, v8)
        cc.cd(3);gPad->SetLogx();gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
        {FileAna&fa=FA[1];const int nT=static_cast<int>(fa.thr.size());
         TH1D* fr=new TH1D("fr3","Eficiencia vs umbral (v8);Umbral [ADC];#varepsilon",100,1,200);
         fr->SetMinimum(0);fr->SetMaximum(1.05);fr->SetStats(0);fr->Draw();
         const int cpanel3[4]={kBlue+1,kBlue-4,kRed+1,kRed-4};
         TGraphErrors* gs[4]={};
         for(int pl=0;pl<4;++pl){
             std::vector<double> x,y,ey;
             for(int ti=0;ti<nT;++ti){
                 double eps=(double)fa.eff_n[pl][ti]/fa.n_events;
                 x.push_back(fa.thr[ti]);y.push_back(eps);
                 ey.push_back(std::sqrt(eps*(1-eps)/fa.n_events));}
             gs[pl]=new TGraphErrors(nT,x.data(),y.data(),
                 std::vector<double>(nT,0).data(),ey.data());
             gs[pl]->SetLineColor(cpanel3[pl]);gs[pl]->SetMarkerColor(cpanel3[pl]);
             gs[pl]->SetMarkerStyle(20+pl);gs[pl]->SetLineWidth(2);gs[pl]->Draw("PL SAME");}
         draw_vline(fa.umbral_nat[P19X],kRed,3);
         TLegend lg(0.45,0.55,0.90,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
         for(int pl=0;pl<4;++pl)lg.AddEntry(gs[pl],pname[pl],"pl");lg.Draw();}

        // Panel 4: Mapa impacto hod2018 v8
        cc.cd(4);gPad->SetLeftMargin(0.14);gPad->SetRightMargin(0.18);gPad->SetBottomMargin(0.14);
        {FA[1].h_imp18->SetTitle("Mapa impacto hod2018 (v8);X_{2018} [mm];Y_{2018} [mm]");
         FA[1].h_imp18->SetStats(0);FA[1].h_imp18->Draw("COLZ");}

        // Panel 5: Calibración espacial hod2019X v8
        cc.cd(5);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
        {FileAna&fa=FA[1];
         std::vector<double> xs,ys,exs,eys;
         for(int c=0;c<16;++c){if(fa.nc_X[c]<2)continue;
             int k0=-1;for(int k=0;k<16;k++)if(kH19X[k]==c+48){k0=k;break;}
             if(k0<0)continue;
             xs.push_back(k0);ys.push_back(fa.xhat[c]);
             exs.push_back(0.25);eys.push_back(fa.xerr[c]);}
         if(!xs.empty()){
             TGraphErrors g(static_cast<int>(xs.size()),xs.data(),ys.data(),exs.data(),eys.data());
             g.SetTitle("Calibración X hod2019 v8;Índice #pi_0;#hat{x}_c [mm]");
             g.SetMarkerStyle(20);g.SetMarkerColor(kBlue+2);g.SetLineColor(kBlue+2);
             g.Draw("AP");
             double sl,ic;linfit(xs,ys,sl,ic);
             TLine lfit(0,ic,15,sl*15+ic);
             lfit.SetLineColor(kRed);lfit.SetLineWidth(2);lfit.Draw("SAME");
             TLine lref(0,-7.5*fa.pitch,15,7.5*fa.pitch);
             lref.SetLineColor(kGray+2);lref.SetLineStyle(2);lref.SetLineWidth(2);lref.Draw("SAME");
             TLatex lat;lat.SetNDC();lat.SetTextSize(0.038);
             char buf[80];std::snprintf(buf,sizeof(buf),"#tau_{K}=%.3f  slope=%.2f",fa.tau_AX,sl);
             lat.DrawLatex(0.15,0.93,buf);}}

        // Panel 6: Ocupancia π₀ vs π_A hod2019X
        cc.cd(6);gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.14);
        {FileAna&fa=FA[1];
         TH1D h0("s_pi0","#pi_0 vs #pi_A (hod2019 X, v8);Índice k;Hits",16,-0.5,15.5);
         TH1D hA("s_piA","",16,-0.5,15.5);
         for(int k=0;k<16;++k){
             h0.SetBinContent(k+1,fa.occ_nat[P19X][kH19X[k]-48]);
             hA.SetBinContent(k+1,fa.occ_nat[P19X][fa.ordA_X[k]]);}
         h0.SetLineColor(kRed);h0.SetLineWidth(2);h0.SetStats(0);
         hA.SetLineColor(kBlue+1);hA.SetLineWidth(2);hA.SetLineStyle(2);
         h0.Draw("HIST");hA.Draw("HIST SAME");
         TLegend lg(0.55,0.78,0.90,0.92);lg.SetBorderSize(0);lg.SetFillStyle(0);
         lg.AddEntry(&h0,"#pi_0","l");lg.AddEntry(&hA,"#pi_A","l");lg.Draw();}

        cc.SaveAs((outdir/"summary"/"figures_of_merit_summary.png").string().c_str());
    }

    //==========================================================================
    // ROOT output
    //==========================================================================
    {
        TFile fout((outdir/"histos_independent.root").string().c_str(),"RECREATE");
        if(!fout.IsZombie()){
            for(int v=0;v<2;++v){
                FileAna&fa=FA[v];
                for(int pl=0;pl<4;++pl)fa.h_ph[pl]->Write();
                fa.h_imp18->Write();fa.h_imp19->Write();
            }
            fout.Close();}
    }

    //==========================================================================
    // Tabla de figuras de mérito + summary.txt
    //==========================================================================
    auto print_table=[&](std::ostream& os){
        os<<"=== Análisis independiente hod2018 y hod2019 ===\n\n";
        os<<"Archivos: v7="<<fv7<<"  v8="<<fv8<<"\n\n";

        // Tabla central
        os<<std::setw(30)<<std::left<<"Métrica"
          <<std::setw(14)<<"18X-v7"<<std::setw(14)<<"18X-v8"
          <<std::setw(14)<<"19X-v7"<<std::setw(14)<<"19X-v8\n";
        os<<"----------------------------------------------------------------------\n";

        auto fmt=[](double v)->std::string{
            if(v<=0||v>9999)return "N/A";
            char buf[20];std::snprintf(buf,sizeof(buf),"%.1f",v);return buf;};
        auto fmtp=[](double v)->std::string{
            char buf[20];std::snprintf(buf,sizeof(buf),"%.3f",v);return buf;};

        os<<std::setw(30)<<"MPV Landau [ADC]"
          <<std::setw(14)<<fmt(FA[0].landau[P18X].mpv)
          <<std::setw(14)<<fmt(FA[1].landau[P18X].mpv)
          <<std::setw(14)<<fmt(FA[0].landau[P19X].mpv)
          <<std::setw(14)<<fmt(FA[1].landau[P19X].mpv)<<"\n";
        os<<std::setw(30)<<"Sigma Landau [ADC]"
          <<std::setw(14)<<fmt(FA[0].landau[P18X].sigma)
          <<std::setw(14)<<fmt(FA[1].landau[P18X].sigma)
          <<std::setw(14)<<fmt(FA[0].landau[P19X].sigma)
          <<std::setw(14)<<fmt(FA[1].landau[P19X].sigma)<<"\n";

        auto snr=[](const FileAna&fa,int pl)->double{
            double mean_noise=0;int nc=0;
            const int base=(pl==P18X)?0:(pl==P18Y)?16:(pl==P19Y)?32:48;
            const std::array<int,16>*prop=(pl==P19X)?&kH19X:(pl==P19Y)?&kH19Y:nullptr;
            for(int k=0;k<16;++k){int c=prop?(*prop)[k]:base+k;
                if(fa.n_ped[c]>0){mean_noise+=fa.noise[c];++nc;}}
            mean_noise=nc>0?mean_noise/nc:1;
            return fa.landau[pl].ok?fa.landau[pl].mpv/mean_noise:0;};
        os<<std::setw(30)<<"SNR (MPV/noise_rms)"
          <<std::setw(14)<<fmt(snr(FA[0],P18X))<<std::setw(14)<<fmt(snr(FA[1],P18X))
          <<std::setw(14)<<fmt(snr(FA[0],P19X))<<std::setw(14)<<fmt(snr(FA[1],P19X))<<"\n";
        os<<std::setw(30)<<"Umbral natural [ADC]"
          <<std::setw(14)<<fmt(FA[0].umbral_nat[P18X])<<std::setw(14)<<fmt(FA[1].umbral_nat[P18X])
          <<std::setw(14)<<fmt(FA[0].umbral_nat[P19X])<<std::setw(14)<<fmt(FA[1].umbral_nat[P19X])<<"\n";

        // ε = fracción de EVENTOS donde al menos 1 strip supera el umbral_natural
        auto eff_nat=[](FileAna&fa,int pl)->double{
            return fa.n_events>0?100.0*fa.eff_nat_events[pl]/fa.n_events:0;};
        os<<std::setw(30)<<"ε(umbral_natural) [%]"
          <<std::setw(14)<<fmt(eff_nat(FA[0],P18X))<<std::setw(14)<<fmt(eff_nat(FA[1],P18X))
          <<std::setw(14)<<fmt(eff_nat(FA[0],P19X))<<std::setw(14)<<fmt(eff_nat(FA[1],P19X))<<"\n";

        // Chi-squared uniformidad hod2019X
        auto chi2_unif=[&](FileAna&fa,int pl)->std::pair<double,double>{
            double sum=0;for(int k=0;k<16;++k)sum+=fa.occ_nat[pl][k];
            const double mean=sum/16;double chi2=0;int dof=0;
            for(int k=0;k<16;++k)if(mean>0){
                chi2+=std::pow((double)fa.occ_nat[pl][k]-mean,2)/mean;++dof;}
            double pv=dof>0?TMath::Prob(chi2,dof-1):0;
            return {chi2,pv};};
        auto [c2_18xv7,pv_18xv7]=chi2_unif(FA[0],P18X);
        auto [c2_18xv8,pv_18xv8]=chi2_unif(FA[1],P18X);
        auto [c2_19xv7,pv_19xv7]=chi2_unif(FA[0],P19X);
        auto [c2_19xv8,pv_19xv8]=chi2_unif(FA[1],P19X);
        os<<std::setw(30)<<"Uniformidad p-value"
          <<std::setw(14)<<fmtp(pv_18xv7)<<std::setw(14)<<fmtp(pv_18xv8)
          <<std::setw(14)<<fmtp(pv_19xv7)<<std::setw(14)<<fmtp(pv_19xv8)<<"\n";
        os<<std::setw(30)<<"N eventos X-Y coinc."
          <<std::setw(14)<<(long long)FA[0].h_imp18->Integral()
          <<std::setw(14)<<(long long)FA[1].h_imp18->Integral()
          <<std::setw(14)<<(long long)FA[0].h_imp19->Integral()
          <<std::setw(14)<<(long long)FA[1].h_imp19->Integral()<<"\n";
        os<<std::setw(30)<<"τ_Kendall(π_A, π₀) X"
          <<std::setw(14)<<"N/A"<<std::setw(14)<<"N/A"
          <<std::setw(14)<<fmtp(FA[0].tau_AX)<<std::setw(14)<<fmtp(FA[1].tau_AX)<<"\n";
        os<<std::setw(30)<<"τ_Kendall(π_A, π₀) Y"
          <<std::setw(14)<<"N/A"<<std::setw(14)<<"N/A"
          <<std::setw(14)<<fmtp(FA[0].tau_AY)<<std::setw(14)<<fmtp(FA[1].tau_AY)<<"\n";

        os<<"\n=== INTERPRETACIONES FÍSICAS ===\n\n";

        // 1. Picos SPE
        os<<"1. Picos SPE en hod2019 v8:\n";
        if(FA[1].spe19x.found){
            os<<"   X: DETECTADOS. "<<FA[1].spe19x.n_peaks<<" picos equidistantes con Δ=";
            os<<std::fixed<<std::setprecision(1)<<FA[1].spe19x.delta<<" ADC.\n";
            if(FA[1].landau[P19X].ok&&FA[1].spe19x.delta>1){
                const double nph=FA[1].landau[P19X].mpv/FA[1].spe19x.delta;
                os<<"   Fotones detectados en promedio para MIP: MPV/Δ ≈ "<<nph<<" fotoelectrones.\n";}
        } else {
            os<<"   X: no detectados (señal demasiado pequeña o distribución continua).\n";}
        if(FA[1].spe19y.found)
            os<<"   Y: DETECTADOS. Δ="<<FA[1].spe19y.delta<<" ADC.\n";
        else os<<"   Y: no detectados.\n";

        // 2. Consistencia eficiencia con Poisson
        os<<"\n2. Eficiencia hod2019 vs modelo de Poisson:\n";
        {const double eps=eff_nat(FA[1],P19X)/100.0;
         if(eps>0&&eps<1){
             const double mu=-std::log(1-eps);
             os<<"   ε(hod2019X,v8)="<<std::setprecision(2)<<eps*100<<"%";
             os<<" → μ=-ln(1-ε)="<<std::setprecision(2)<<mu<<" fotones esperados.\n";
             if(FA[1].spe19x.found&&FA[1].landau[P19X].ok){
                 const double mu_from_mpv=FA[1].landau[P19X].mpv/FA[1].spe19x.delta;
                 os<<"   μ desde MPV/Δ = "<<mu_from_mpv<<". ";
                 if(std::abs(mu-mu_from_mpv)<mu*0.5)os<<"Consistentes.\n";
                 else os<<"Discrepancia (incertidumbre en Δ o en ajuste Landau).\n";}
             else os<<"   (No disponible SPE o Landau para comparar.)\n";}
         else os<<"   ε=0 o ε=1, no aplicable.\n";}

        // 3. Calibración espacial
        os<<"\n3. Calibración espacial hod2019X:\n";
        for(int v=0;v<2;++v){
            os<<"   "<<FA[v].label<<": τ_Kendall(π_A,π₀)="
              <<std::setprecision(4)<<FA[v].tau_AX;
            if(FA[v].tau_AX>0.75) os<<" → TENDENCIA POSITIVA CONSISTENTE con π₀.\n";
            else if(FA[v].tau_AX<-0.3) os<<" → tendencia negativa (posible inversión de ejes).\n";
            else os<<" → no concluyente (estadística insuficiente).\n";
            int has_data=0;for(int c=0;c<16;++c)if(FA[v].nc_X[c]>=2)++has_data;
            os<<"   Canales con ≥2 eventos: "<<has_data<<"/16.\n";}

        // 4. Uniformidad hod2018
        os<<"\n4. Uniformidad de ocupancia hod2018:\n";
        for(int v=0;v<2;++v){
            auto [chi2,pv]=chi2_unif(FA[v],P18X);
            os<<"   hod2018X v"<<(v==0?"7":"8")<<": p-value="<<std::setprecision(4)<<pv;
            if(pv<0.05) os<<" → NO UNIFORME (canales problemáticos o perfil no plano).\n";
            else os<<" → compatible con distribución uniforme.\n";
            // Top 3 canales más activos y menos activos (no reportar todos para claridad)
            double sum_occ=0;for(int k=0;k<16;++k)sum_occ+=FA[v].occ_nat[P18X][k];
            const double mean=sum_occ/16;
            std::vector<std::pair<long long,int>> ranked;
            for(int k=0;k<16;++k)ranked.push_back({FA[v].occ_nat[P18X][k],k});
            std::sort(ranked.begin(),ranked.end(),[](const auto&a,const auto&b){return a.first>b.first;});
            os<<"   Top 3 más activos: ";
            for(int i=0;i<3&&i<(int)ranked.size();++i)
                os<<"k="<<ranked[i].second<<"(Ch"<<ranked[i].second
                  <<",hits="<<ranked[i].first<<") ";
            os<<"\n   Top 3 menos activos: ";
            for(int i=15;i>12&&i>=0;--i)
                os<<"k="<<ranked[i].second<<"(Ch"<<ranked[i].second
                  <<",hits="<<ranked[i].first<<") ";
            os<<"\n";
        }
    };

    print_table(std::cout);
    {std::ofstream sf(outdir/"summary"/"summary_independent.txt");print_table(sf);}

    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
