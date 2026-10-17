// reconstruct_corrected.cxx
// Primer tracking real con ventana temporal CORRECTA para ambos hodoscopios.
// Condición nueva verificada: muones en [4000,4400] ns para hod2018 Y hod2019 (Δt=0).
// El análisis anterior usaba t=3510 ns (fondo) → fallaba por construcción.
//
// Paso 1: tag-and-probe corregido (confirma señal si sig>5σ).
// Paso 2: tracking completo con los 4 planos.
// Paso 3: figuras de correlación espacial, ángulos, mapa 2D.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o reconstruct_corrected reconstruct_corrected.cxx

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
#include <cstdlib>
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
// Parámetros congelados (NO cambiar sin evidencia nueva)
//==============================================================================
// Ventana del muón verificada: ambos hodoscopios a t≈4170 ns
constexpr double kMuonTlo =  4000.0;  // ns
constexpr double kMuonThi =  4400.0;  // ns
constexpr double kPrePulse =  2500.0;  // ns (región de baseline)
constexpr double kDt_ns = 16.0;  // VX1740, 62.5 MS/s; equals the data_dig.time step
constexpr int    kNpre=3, kNpost=8, kNwin=12; // integración Q

// Umbrales verificados
constexpr double kPh18 = 50.0;  // mV, pulse height from pedestal; = 102 ADC counts
constexpr double kPh19 = 35.0;  // mV, pulse height from pedestal; = 72 ADC counts

// Geometría real de la medición
constexpr double kPitch = 2.0;   // mm
constexpr double kD_mm  = 98.4;  // mm (separación entre planos)
constexpr int    kNCh   = 64;

// Mapeo autoritativo (literal, sin reglas algebraicas)
static const std::array<int,16> kH19X={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
static const std::array<int,16> kH19Y={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};

//==============================================================================
// Hit en la ventana del muón: extrae el pulso de máxima amplitud en [4000,4400] ns
//==============================================================================
struct MuonHit {
    bool   valid  = false;
    int    k      = -1;     // índice físico (0..15)
    double ph     = 0;      // pulse height [mV]
    double pt     = -1;     // tiempo del mínimo [ns]
    double Q      = 0;      // carga integrada [mV·ns]
    double pos    = 0;      // posición física [mm]
};

static MuonHit extract_muon_hit(const std::vector<float>& wf,
                                  const std::vector<float>& t,
                                  int k) {
    MuonHit h; h.k = k;
    const int n = (int)wf.size();
    if (n < 4 || (int)t.size() != n) return h;
    double s = 0; int nb = 0;
    for (int i=0;i<n;i++) if(t[i]<kPrePulse){s+=wf[i];nb++;}
    if (nb < 2) return h;
    const double ped = s/nb;
    double vmin = 1e30; int imin = -1;
    for (int i=0;i<n;i++)
        if (t[i]>=kMuonTlo && t[i]<=kMuonThi && wf[i]<vmin) {vmin=wf[i];imin=i;}
    if (imin < 0) return h;
    h.valid = true;
    h.ph  = ped - vmin;
    h.pt  = t[imin];
    // Carga integrada alrededor del mínimo
    const double tilo = h.pt - kNpre*kDt_ns, tihi = h.pt + kNpost*kDt_ns;
    for (int i=0;i<n;i++) if(t[i]>=tilo && t[i]<=tihi) h.Q += (ped-wf[i])*kDt_ns;
    h.pos = (k + 0.5 - 8.0) * kPitch;
    return h;
}

//==============================================================================
// Datos por evento (RAM)
//==============================================================================
struct EvData {
    // Hits válidos por plano: máxima amplitud que pasa el umbral en [4000,4400]
    MuonHit best18x, best18y, best19x, best19y;
    // Todos los hits válidos (para clustering)
    std::vector<MuonHit> hits18x, hits18y, hits19x, hits19y;
};

//==============================================================================
// Clustering (gap≤1 en índice k)
//==============================================================================
struct Cluster {double pos=0,amp=0,time=0;int lo=0,hi=0,n=0;};
static std::vector<Cluster> make_cls(const std::vector<MuonHit>& hits, double ph_thr){
    std::vector<MuonHit> act;
    for (const auto& h:hits) if(h.valid && h.ph>=ph_thr) act.push_back(h);
    if (act.empty()) return {};
    std::sort(act.begin(),act.end(),[](const MuonHit&a,const MuonHit&b){return a.k<b.k;});
    std::vector<Cluster> out; Cluster cur{};
    auto flush=[&](){cur.pos/=cur.amp;cur.time/=cur.amp;out.push_back(cur);};
    auto open=[&](const MuonHit&h){cur={};cur.lo=cur.hi=h.k;};
    auto add=[&](const MuonHit&h){cur.pos+=h.ph*h.pos;cur.time+=h.ph*h.pt;cur.amp+=h.ph;++cur.n;cur.hi=h.k;};
    open(act[0]);add(act[0]);
    for(int i=1;i<(int)act.size();i++){if(act[i].k>cur.hi+1){flush();open(act[i]);}add(act[i]);}
    flush(); return out;
}

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
// Main
//==============================================================================
int main(int argc,char**argv){
    std::string fname="lectura_vx1740_ext_nim_1M_v8.root";
    if(argc>1)fname=argv[1];

    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("corrected_tracking_"+std::string(ds));
    fs::create_directories(outdir);

    std::cout<<"=== reconstruct_corrected ===\n";
    std::cout<<"File: "<<fname<<"\nOut: "<<outdir<<"\n";
    std::cout<<"Muon window: ["<<kMuonTlo<<", "<<kMuonThi<<"] ns (BOTH hodoscopes)\n";
    std::cout<<"Thresholds: hod2018 ph>"<<kPh18<<" mV,  hod2019 ph>"<<kPh19<<" mV\n\n";

    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    TFile* fin=TFile::Open(fname.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: "<<fname<<'\n';return 1;}
    TTree* tree=dynamic_cast<TTree*>(fin->Get("data_dig"));
    if(!tree){std::cerr<<"ERROR: no data_dig\n";return 1;}

    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw=nullptr;
    tree->SetBranchStatus("*",0);tree->SetCacheSize(256LL<<20);
    for(int c=0;c<kNCh;c++){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!tree->GetBranch(bn))continue;
        tree->SetBranchStatus(bn,1);tree->SetBranchAddress(bn,&cw[c]);}
    tree->SetBranchStatus("time",1);tree->SetBranchAddress("time",&tw);
    Long_t ts_val = 0;
    if(tree->GetBranch("timestamp")){
        tree->SetBranchStatus("timestamp",1);
        tree->SetBranchAddress("timestamp",&ts_val);
    }

    // Inverse maps for hod2019
    std::array<int,kNCh> inv19x{},inv19y{};inv19x.fill(-1);inv19y.fill(-1);
    for(int k=0;k<16;k++){inv19x[kH19X[k]]=k;inv19y[kH19Y[k]]=k;}

    const Long64_t nEnt=tree->GetEntries();
    std::cout<<"Loading "<<nEnt<<" events into RAM...\n";
    std::vector<EvData> evs; evs.reserve(nEnt);

    Long_t ts_first = 0, ts_last = 0;
    bool checked_step = false;
    for(Long64_t iev=0;iev<nEnt;iev++){
        tree->GetEntry(iev);
        if(iev%30000==0)std::cout<<"  "<<iev<<"/"<<nEnt<<"\r"<<std::flush;
        if(iev==0) ts_first = ts_val;
        if(iev==nEnt-1) ts_last = ts_val;
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;
        if (!checked_step && tv.size() > 1) {
            checked_step = true;
            if (std::fabs((tv[1] - tv[0]) - kDt_ns) > 1e-3) {
                std::fprintf(stderr, "Sampling step %.4f ns != kDt_ns %.4f ns; aborting\n", tv[1]-tv[0], kDt_ns);
                std::exit(1);
            }
        }
        EvData ev;
        for(int c=0;c<kNCh;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const bool is19=(c>=32);
            int k=-1;
            if(c>=0&&c<=15)  k=c;
            else if(c>=16&&c<=31) k=c-16;
            else if(inv19x[c]>=0) k=inv19x[c];
            else if(inv19y[c]>=0) k=inv19y[c];
            if(k<0)continue;
            MuonHit h=extract_muon_hit(*cw[c],tv,k);
            if(!h.valid)continue;
            // Store in appropriate plane
            if(c>=0&&c<=15)  ev.hits18x.push_back(h);
            else if(c>=16&&c<=31) ev.hits18y.push_back(h);
            else if(inv19x[c]>=0) ev.hits19x.push_back(h);
            else if(inv19y[c]>=0) ev.hits19y.push_back(h);
        }
        evs.push_back(std::move(ev));
    }
    std::cout<<"  "<<nEnt<<"/"<<nEnt<<" ✓  ("<<evs.size()<<" events cached)\n\n";
    fin->Close();

    //==========================================================================
    // PASO 1: Tag-and-probe corregido
    //==========================================================================
    std::cout<<"=== PASO 1: Tag-and-probe con ventana correcta ===\n";
    std::vector<size_t> shuf(evs.size());
    std::iota(shuf.begin(),shuf.end(),0);
    {std::mt19937 rng(42);std::shuffle(shuf.begin(),shuf.end(),rng);}

    long long n_tag=0, n_same=0, n_mixed=0;
    for(size_t i=0;i<evs.size();i++){
        // Tag: hod2018X has exactly 1 cluster with ph>kPh18 in muon window
        auto cls18x=make_cls(evs[i].hits18x,kPh18);
        if(cls18x.size()!=1)continue;
        const int k_pred=(int)std::round(cls18x[0].pos/kPitch+7.5);
        if(k_pred<0||k_pred>15)continue;
        ++n_tag;
        // Probe same-event: hod2019X at k_pred with ph>kPh19
        bool same_hit=false;
        for(const auto& h:evs[i].hits19x)
            if(h.k==k_pred && h.ph>=kPh19){same_hit=true;break;}
        if(same_hit)++n_same;
        // Probe mixed-event
        const size_t j=shuf[i];
        bool mixed_hit=false;
        for(const auto& h:evs[j].hits19x)
            if(h.k==k_pred && h.ph>=kPh19){mixed_hit=true;break;}
        if(mixed_hit)++n_mixed;
    }
    const double Nt=static_cast<double>(n_tag);
    const double eps_obs=n_tag>0?(double)n_same/Nt:0;
    const double eps_acc=n_tag>0?(double)n_mixed/Nt:0;
    const double eps_true=eps_obs-eps_acc;
    const double so=std::sqrt(eps_obs*(1-eps_obs)/Nt);
    const double sa=std::sqrt(eps_acc*(1-eps_acc)/Nt);
    const double sig_tp=std::sqrt(so*so+sa*sa)>0?eps_true/std::sqrt(so*so+sa*sa):0;

    std::cout<<"N_tag="<<n_tag<<"  N_same="<<n_same<<"  N_mixed="<<n_mixed<<"\n";
    std::cout<<"ε_obs="<<eps_obs*100<<"% ε_acc="<<eps_acc*100<<"% ε_true="<<eps_true*100<<"%\n";
    std::cout<<"Significancia tag-and-probe: "<<sig_tp<<"σ\n";
    if(sig_tp<5.0)std::cout<<"[!] sig<5σ — tracking puede ser limitado por estadística o geometría\n";
    else std::cout<<"[✓] sig>5σ — señal confirmada, proceder al tracking completo\n";
    std::cout<<"\n";

    //==========================================================================
    // PASO 2: Tracking completo
    //==========================================================================
    std::cout<<"=== PASO 2: Tracking completo ===\n";

    long long cf_any=0,cf_18c=0,cf_19c=0,cf_XZ=0,cf_YZ=0,cf_3D=0;
    std::vector<double> vx18,vx19,vy18,vy19,vmx,vmy;

    // Histogramas
    constexpr double kClo=-16,kChi=16;constexpr int kCb=16;
    auto h_cx=std::make_unique<TH2D>("h_cx","Correlacion X;X_{2018}[mm];X_{2019}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_cy=std::make_unique<TH2D>("h_cy","Correlacion Y;Y_{2018}[mm];Y_{2019}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);
    auto h_thx=std::make_unique<TH1D>("h_thx","#theta_x;#theta_x [deg];Events",80,-25,25);
    auto h_thy=std::make_unique<TH1D>("h_thy","#theta_y;#theta_y [deg];Events",80,-25,25);
    auto h_th3=std::make_unique<TH1D>("h_th3","#theta_{3D};#theta_{3D} [deg];Events",80,0,25);
    auto h_xy=std::make_unique<TH2D>("h_xy","Track 2D map;X_{2018}[mm];Y_{2018}[mm]",kCb,kClo,kChi,kCb,kClo,kChi);

    std::ofstream fcsv(outdir/"cosmic_tracks_corrected.csv");
    fcsv<<"x2018_mm,y2018_mm,t18_ns,x2019_mm,y2019_mm,t19_ns,mx,my,theta_x_deg,theta_y_deg,theta_3d_deg\n";

    for(auto& ev:evs){
        // Cluster cada plano con los umbrales verificados
        auto cls18x=make_cls(ev.hits18x,kPh18);
        auto cls18y=make_cls(ev.hits18y,kPh18);
        auto cls19x=make_cls(ev.hits19x,kPh19);
        auto cls19y=make_cls(ev.hits19y,kPh19);

        bool has18=(cls18x.size()>=1||cls18y.size()>=1);
        bool has19=(cls19x.size()>=1||cls19y.size()>=1);
        if(has18||has19)++cf_any;
        if(cls18x.size()>=1&&cls18y.size()>=1)++cf_18c;
        if(cls19x.size()>=1&&cls19y.size()>=1)++cf_19c;

        // Selecciona el cluster de mayor amplitud por plano
        auto best=[](const std::vector<Cluster>&v)->const Cluster*{
            if(v.empty())return nullptr;
            return &*std::max_element(v.begin(),v.end(),[](const Cluster&a,const Cluster&b){return a.amp<b.amp;});};

        const Cluster* cx18=best(cls18x); const Cluster* cy18=best(cls18y);
        const Cluster* cx19=best(cls19x); const Cluster* cy19=best(cls19y);

        const bool has_XZ=(cx18&&cx19); if(has_XZ)++cf_XZ;
        const bool has_YZ=(cy18&&cy19); if(has_YZ)++cf_YZ;
        const bool has_3D=(has_XZ&&has_YZ); if(has_3D)++cf_3D;

        if(has_XZ){vx18.push_back(cx18->pos);vx19.push_back(cx19->pos);h_cx->Fill(cx18->pos,cx19->pos);}
        if(has_YZ){vy18.push_back(cy18->pos);vy19.push_back(cy19->pos);h_cy->Fill(cy18->pos,cy19->pos);}
        if(has_3D){
            const double mx=(cx19->pos-cx18->pos)/kD_mm;
            const double my=(cy19->pos-cy18->pos)/kD_mm;
            const double thx=std::atan(mx)*180/M_PI;
            const double thy=std::atan(my)*180/M_PI;
            const double th3=std::atan(std::sqrt(mx*mx+my*my))*180/M_PI;
            vmx.push_back(mx);vmy.push_back(my);
            h_thx->Fill(thx);h_thy->Fill(thy);h_th3->Fill(th3);
            h_xy->Fill(cx18->pos,cy18->pos);
            fcsv<<std::fixed<<std::setprecision(3)
                <<cx18->pos<<','<<cy18->pos<<','<<cx18->time<<','
                <<cx19->pos<<','<<cy19->pos<<','<<cx19->time<<','
                <<mx<<','<<my<<','<<thx<<','<<thy<<','<<th3<<'\n';
        }
    }

    // Correlaciones
    const double rx=pearsonR(vx18,vx19),ry=pearsonR(vy18,vy19);
    double slx,icx,sly,icy;
    linfit(vx18,vx19,slx,icx);linfit(vy18,vy19,sly,icy);

    std::cout<<"Cutflow:\n";
    std::cout<<"  Eventos con algún hit en ventana muón: "<<cf_any<<"\n";
    std::cout<<"  hod2018 X+Y coincidencia:              "<<cf_18c<<"\n";
    std::cout<<"  hod2019 X+Y coincidencia:              "<<cf_19c<<"\n";
    std::cout<<"  Trazas XZ (x2018+x2019):               "<<cf_XZ<<"\n";
    std::cout<<"  Trazas YZ (y2018+y2019):               "<<cf_YZ<<"\n";
    std::cout<<"  Trazas 3D completas:                   "<<cf_3D<<"\n\n";
    std::cout<<"Correlación espacial:\n";
    std::cout<<"  X: r="<<rx<<"  slope="<<slx<<"  intercept="<<icx<<" mm  (n="<<vx18.size()<<")\n";
    std::cout<<"  Y: r="<<ry<<"  slope="<<sly<<"  intercept="<<icy<<" mm  (n="<<vy18.size()<<")\n\n";

    // Ángulos
    if(!vmx.empty()){
        double sumth=0;for(double m:vmx)sumth+=std::atan(m)*180/M_PI;
        const double mean_thx=sumth/vmx.size();
        double sumth2=0;for(double m:vmx){double d=std::atan(m)*180/M_PI-mean_thx;sumth2+=d*d;}
        const double rms_thx=std::sqrt(sumth2/vmx.size());
        std::cout<<"theta_x: mean="<<mean_thx<<" deg  RMS="<<rms_thx<<" deg\n";
    }

    //==========================================================================
    // PASO 3: Figuras
    //==========================================================================
    auto p2fit=[&](TH2D*h,double sl,double ic,double r,const fs::path&p){
        TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",900,820);
        c.SetLeftMargin(0.13);c.SetRightMargin(0.17);c.SetBottomMargin(0.13);
        h->SetStats(0);h->Draw("COLZ");
        const double xlo=h->GetXaxis()->GetXmin(),xhi=h->GetXaxis()->GetXmax();
        TLine ln(xlo,sl*xlo+ic,xhi,sl*xhi+ic);ln.SetLineColor(kRed);ln.SetLineWidth(2);ln.Draw("SAME");
        TLatex lt;lt.SetNDC();lt.SetTextSize(0.033);
        char buf[100];std::snprintf(buf,sizeof(buf),"slope=%.3f  intercept=%.2f mm  r=%.3f  N=%d",sl,ic,r,(int)std::max(vx18.size(),vy18.size()));
        lt.DrawLatex(0.14,0.935,buf);
        c.SaveAs(p.string().c_str());};

    // FIG A: Correlación espacial
    h_cx->SetTitle("Correlacion espacial X;X_{2018}[mm];X_{2019}[mm];Events");
    h_cy->SetTitle("Correlacion espacial Y;Y_{2018}[mm];Y_{2019}[mm];Events");
    {TCanvas cc("cc","",1600,780);cc.Divide(2,1);
     cc.cd(1);gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.17);gPad->SetBottomMargin(0.13);
     h_cx->Draw("COLZ");
     TLine lnx(kClo,slx*kClo+icx,kChi,slx*kChi+icx);lnx.SetLineColor(kRed);lnx.SetLineWidth(2);lnx.Draw("SAME");
     TLatex ltx;ltx.SetNDC();ltx.SetTextSize(0.036);
     char bx[80];std::snprintf(bx,sizeof(bx),"slope=%.3f  r=%.3f  N=%d",slx,rx,(int)vx18.size());ltx.DrawLatex(0.14,0.935,bx);
     cc.cd(2);gPad->SetLeftMargin(0.13);gPad->SetRightMargin(0.17);gPad->SetBottomMargin(0.13);
     h_cy->Draw("COLZ");
     TLine lny(kClo,sly*kClo+icy,kChi,sly*kChi+icy);lny.SetLineColor(kRed);lny.SetLineWidth(2);lny.Draw("SAME");
     TLatex lty;lty.SetNDC();lty.SetTextSize(0.036);
     char by[80];std::snprintf(by,sizeof(by),"slope=%.3f  r=%.3f  N=%d",sly,ry,(int)vy18.size());lty.DrawLatex(0.14,0.935,by);
     cc.SaveAs((outdir/"fig_spatial_correlation_XY.png").string().c_str());}

    // FIG B: Distribuciones angulares
    {TCanvas ca("ca","",1200,450);ca.Divide(3,1);
     for(int i=0;i<3;i++){ca.cd(i+1);gPad->SetLeftMargin(0.13);gPad->SetBottomMargin(0.14);}
     ca.cd(1);h_thx->SetLineColor(kBlue+1);h_thx->SetLineWidth(2);h_thx->SetStats(0);h_thx->Draw("HIST");
     ca.cd(2);h_thy->SetLineColor(kBlue+1);h_thy->SetLineWidth(2);h_thy->SetStats(0);h_thy->Draw("HIST");
     ca.cd(3);h_th3->SetLineColor(kBlue+1);h_th3->SetLineWidth(2);h_th3->SetStats(0);h_th3->Draw("HIST");
     ca.SaveAs((outdir/"fig_angular_distributions.png").string().c_str());}

    // FIG C: Mapa 2D de tracks 3D
    {TCanvas cm("cm","",900,820);cm.SetLeftMargin(0.13);cm.SetRightMargin(0.17);cm.SetBottomMargin(0.13);
     h_xy->SetTitle("2D impact map (3D tracks);X_{2018}[mm];Y_{2018}[mm];Events");
     h_xy->SetStats(0);h_xy->Draw("COLZ");
     cm.SaveAs((outdir/"fig_track_2Dmap.png").string().c_str());}

    // FIG D: N_tracks vs umbral hod2019
    {std::vector<double> ph_thrs={10,15,20,25,30,35,40,50,75,100};
     std::vector<double> n_xz_v,n_3d_v;
     for(double thr:ph_thrs){
         long long nxz=0,n3d=0;
         for(auto&ev:evs){
             auto cx18=make_cls(ev.hits18x,kPh18);
             auto cx19=make_cls(ev.hits19x,thr);
             auto cy18=make_cls(ev.hits18y,kPh18);
             auto cy19=make_cls(ev.hits19y,thr);
             if(!cx18.empty()&&!cx19.empty())++nxz;
             if(!cx18.empty()&&!cx19.empty()&&!cy18.empty()&&!cy19.empty())++n3d;
         }
         n_xz_v.push_back(nxz);n_3d_v.push_back(n3d);}
     TCanvas ct("ct","",900,620);ct.SetLeftMargin(0.12);ct.SetBottomMargin(0.13);
     int np=(int)ph_thrs.size();
     TGraph gxz(np,ph_thrs.data(),n_xz_v.data()),g3d(np,ph_thrs.data(),n_3d_v.data());
     gxz.SetTitle("N tracks vs umbral hod2019;ph_{thr,hod2019} [mV];N tracks");
     gxz.SetMarkerStyle(20);gxz.SetLineColor(kBlue+1);gxz.SetMarkerColor(kBlue+1);gxz.SetLineWidth(2);
     g3d.SetMarkerStyle(21);g3d.SetLineColor(kRed+1);g3d.SetMarkerColor(kRed+1);g3d.SetLineWidth(2);
     gxz.Draw("APL");g3d.Draw("PL SAME");
     TLegend lg(0.55,0.75,0.88,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
     lg.AddEntry(&gxz,"N_XZ","pl");lg.AddEntry(&g3d,"N_3D","pl");lg.Draw();
     ct.SaveAs((outdir/"fig_tracks_vs_threshold.png").string().c_str());}

    // ROOT file
    {TFile fout((outdir/"histos.root").string().c_str(),"RECREATE");
     if(!fout.IsZombie()){
         for(auto*h:std::initializer_list<TH1*>{h_thx.get(),h_thy.get(),h_th3.get()})h->Write();
         for(auto*h:std::initializer_list<TH2D*>{h_cx.get(),h_cy.get(),h_xy.get()})h->Write();
         fout.Close();}}

    //==========================================================================
    // PASO 4: Summary
    //==========================================================================
    const double dur_h = (nEnt > 0 && ts_last > ts_first) ? (double)(ts_last - ts_first) / 3.6e9 : 0.0;
    auto print=[&](std::ostream&os){
        os<<"=== Primer Tracking Real con Ventana Temporal Correcta ===\n\n";
        os<<"Ventana muón: ["<<kMuonTlo<<", "<<kMuonThi<<"] ns (AMBOS hodoscopios, Δt=0 verificado)\n";
        os<<"Umbral hod2018: ph>"<<kPh18<<" mV  |  hod2019: ph>"<<kPh19<<" mV\n\n";
        os<<"PASO 1 — Tag-and-probe corregido:\n";
        os<<"  N_tag="<<n_tag<<"  N_same="<<n_same<<"  N_mixed="<<n_mixed<<"\n";
        os<<"  ε_true="<<eps_true*100<<"% ± "<<std::sqrt(so*so+sa*sa)*100<<"%\n";
        os<<"  Significancia: "<<sig_tp<<"σ\n\n";
        os<<"PASO 2 — Cutflow tracking:\n";
        os<<"  Eventos con hit en ventana muón: "<<cf_any<<"\n";
        os<<"  hod2018 X+Y:                     "<<cf_18c<<"\n";
        os<<"  hod2019 X+Y:                     "<<cf_19c<<"\n";
        os<<"  Trazas XZ:                       "<<cf_XZ<<"\n";
        os<<"  Trazas YZ:                       "<<cf_YZ<<"\n";
        os<<"  Trazas 3D completas:             "<<cf_3D<<"\n\n";
        os<<"Correlación espacial (muon window [4000,4400] ns):\n";
        os<<std::fixed<<std::setprecision(4);
        os<<"  X: r="<<rx<<"  slope="<<slx<<"  intercept="<<icx<<" mm\n";
        os<<"  Y: r="<<ry<<"  slope="<<sly<<"  intercept="<<icy<<" mm\n\n";
        // Comparison with expected
        constexpr double kFBoth = 0.0167; // fraction of bench triggers crossing both stations; geometric MC with dN/dOmega ~ cos^3(theta), D = 98.4 mm (Muographers 2026 proceeding, section "Fraction of triggers that cross the stations")
        const double expected = kFBoth * nEnt;
        os<<"N_3D observado: "<<cf_3D<<"  vs esperado: "<<std::fixed<<std::setprecision(0)<<expected<<" (del MC geométrico)\n";
        os<<"Ratio obs/exp: "<<std::setprecision(4)<<(expected>0?(double)cf_3D/expected:0)<<"\n\n";
        // Citable text
        os<<"TEXTO CITABLE:\n";
        os<<'"';
        os<<"Using the corrected muon time window [4000,4400] ns for hod2018 ";
        os<<"(eliminating the t≈3510 ns background), ";
        os<<"we reconstruct N_3D="<<cf_3D<<" three-dimensional cosmic-muon tracks ";
        os<<"in run v8 (143k events, "<<std::setprecision(2)<<dur_h<<" h), ";
        os<<"with spatial correlation r_X="<<std::setprecision(3)<<rx;
        os<<" (slope="<<slx<<").";
        os<<'"';
        os<<"\n\n";
        // Mapping assessment
        os<<"Evaluación del mapeo (pendiente≈1 → mapeo correcto):\n";
        if(cf_3D>=10){
            os<<"  Pendiente X="<<std::setprecision(3)<<slx<<"\n";
            os<<"  Pendiente Y="<<sly<<"\n";
        } else {
            os<<"  Insuficientes tracks para evaluación de mapeo (N_3D="<<cf_3D<<"<10)\n";}
    };
    print(std::cout);
    {std::ofstream sf(outdir/"summary_tracking_corrected.txt");print(sf);}

    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
