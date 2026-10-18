// multirun_analysis.cxx
// Análisis multi-run: barrido de sobretensión MPPC (v1..v10).
// Por cada archivo: ventana per-run, espectros de carga, tag-and-probe.
// Combinación global: H_same + H_mixed en unidades de sigma del pedestal → significancia máxima.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o multirun_analysis multirun_analysis.cxx

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

#include "hodo_consts.h"
#include "hodo_mapping.h"

namespace fs = std::filesystem;

//==============================================================================
// Constantes fijas
//==============================================================================
constexpr double kPre_ns = 2500.0;
constexpr int    kNCh = 64;
constexpr double kPitch = 2.0;

// Hot channels identificados en análisis previo (v8)

static const std::vector<int> kHotX={0,10};   // k indices (Ch53, Ch50)
static const std::vector<int> kHotY={3,10,15}; // k indices (Ch46, Ch34, Ch40)

//==============================================================================
// Extracción Q + tiempo de pulso (per strip)
//==============================================================================
struct QResult{double Q=0,pt=-1;};

static QResult get_Q(const std::vector<float>&wf,const std::vector<float>&t,
                      double tlo,double thi,int kNpost,bool use_loc){
    QResult r;
    const int n=(int)wf.size();
    if(n<4||(int)t.size()!=n)return r;
    double s=0,s2=0;int nb=0;
    for(int i=0;i<n;i++)if(t[i]<kPre_ns){s+=wf[i];s2+=wf[i]*wf[i];nb++;}
    if(nb<2)return r;
    const double ped=s/nb;
    const double var=s2/nb-ped*ped;
    const double noise=(var>0)?std::sqrt(var):0;
    double vmin=1e30;int imin=-1;
    for(int i=0;i<n;i++)if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if(imin<0)return r;
    r.pt=t[imin];
    // Baseline for Q integration
    double ped_Q = ped;
    if(use_loc && imin >= 12 && imin - 4 < n){
        double s_loc = 0;
        for(int j = -12; j <= -4; ++j) s_loc += wf[imin + j];
        ped_Q = s_loc / 9.0;
    }
    const double tilo=r.pt-kNpre*kDt_ns, tihi=r.pt+kNpost*kDt_ns;
    double Q=0;
    for(int i=0;i<n;i++)if(t[i]>=tilo&&t[i]<=tihi)Q+=(ped_Q-wf[i])*kDt_ns;
    r.Q=Q;
    return r;
}

//==============================================================================
// Resultado por run
//==============================================================================
struct RunResult {
    std::string fname, label;
    long long n_events=0;
    // Ventana de señal (fija en kMuonTlo..kMuonThi) y picos diagnósticos (ph > umbral)
    double t_lo_19=kMuonTlo, t_hi_19=kMuonThi, t_peak_19=0;
    double t_lo_18=kMuonTlo, t_hi_18=kMuonThi, t_peak_18=0;
    // Estadísticas de carga por plano (0=hod2018X,1=hod2018Y,2=hod2019X,3=hod2019Y)
    double ped_mean[4]={}, noise_mean[4]={}, sigmaQ[4]={};
    double mpvQ[4]={};   // MPV del ajuste Landau
    double snrQ[4]={};   // SNR_Q = MPV/sigma_Q
    double umbral_Q[4]={};  // valle en el espectro Q
    // Estimador robusto hod2019X (señal con ph > kPh19_mV en [kMuonTlo, kMuonThi])
    long long n_sig_19x = 0;
    double med_Q_19x = 0, q10_Q_19x = 0, q90_Q_19x = 0;
    // Landau fit hod2019X
    int lan_status = -1;
    double lan_chi2ndf = 0;
    bool lan_at_limit = false;
    bool lan_ok = false;
    int blo_19x = 0, bhi_19x = 0, min_b_19x = 0;
    double v_xlo_19x = 0, v_xhi_19x = 0;
    double min_v_19x = 0, c_lo_19x = 0, c_hi_19x = 0;
    std::string fail_reason_19x;
    bool has_valley_19x = false;
    // Tag-and-probe (hod2019X, all channels)
    long long n_tag_x=0, n_obs_x=0, n_acc_x=0;
    double eps_true_x=0, sigma_eps_x=0, sig_x=0;
    // Reproductibilidad (vs run anterior)
    double delta_mpv_rel=0;
};

//==============================================================================
// Clustering hod2018
//==============================================================================
struct Hit18{int k;double ph;};
struct Cl18{double pos,amp;};
static std::vector<Cl18> clust18(std::vector<Hit18> hits,double thr){
    std::vector<Hit18> act;
    for(const auto&h:hits)if(h.ph>=thr)act.push_back(h);
    if(act.empty())return{};
    std::sort(act.begin(),act.end(),[](const Hit18&a,const Hit18&b){return a.k<b.k;});
    std::vector<Cl18> out;double wp=0,wa=0;int hi=act[0].k;
    auto flush=[&](){out.push_back({wp/wa,wa});wp=wa=0;};
    auto add=[&](const Hit18&h){wp+=h.ph*(h.k+0.5-8.0)*kPitch;wa+=h.ph;hi=h.k;};
    hi=act[0].k;add(act[0]);
    for(int i=1;i<(int)act.size();i++){if(act[i].k>hi+1){flush();hi=act[i].k;}add(act[i]);}
    flush();return out;}

//==============================================================================
// Procesa un archivo (pasada única)
//==============================================================================
static RunResult process_run(const std::string& fname, const std::string& lbl,
                              TH1D* h_Qcomb_same, TH1D* h_Qcomb_mixed,
                              const std::vector<size_t>& shuf_global,
                              size_t global_offset) {
    RunResult rr; rr.fname=fname; rr.label=lbl;

    TFile* fin=TFile::Open(fname.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"WARNING: cannot open "<<fname<<"\n";return rr;}
    TTree* t=static_cast<TTree*>(fin->Get("data_dig"));
    if(!t){fin->Close();return rr;}

    std::array<std::vector<float>*,kNCh> cw{};
    std::vector<float>* tw=nullptr;
    t->SetBranchStatus("*",0);t->SetCacheSize(128LL<<20);
    for(int c=0;c<kNCh;c++){
        char bn[16];std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!t->GetBranch(bn))continue;
        t->SetBranchStatus(bn,1);t->SetBranchAddress(bn,&cw[c]);}
    t->SetBranchStatus("time",1);t->SetBranchAddress("time",&tw);

    const Long64_t nEnt=t->GetEntries();
    rr.n_events=nEnt;
    std::cout<<"  ["<<lbl<<"] "<<nEnt<<" eventos...\n";

    // Build inverse maps
    std::array<int,kNCh> inv19x{},inv19y{}; inv19x.fill(-1);inv19y.fill(-1);
    for(int k=0;k<16;k++){inv19x[kH19X[k]]=k;inv19y[kH19Y[k]]=k;}

    // PASO 0: t_peak diagnóstico por estación desde pulsos con ph > umbral
    // (scan de 3000 eventos con umbral de amplitud)
    TH1D h_pt19("_pt19","",100,2500,6000);
    TH1D h_pt18("_pt18","",100,2500,8000);
    const Long64_t scan_n=std::min(nEnt,(Long64_t)3000);
    bool checked_step = false;
    for(Long64_t iev=0;iev<scan_n;iev++){
        t->GetEntry(iev);
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;
        if (!checked_step && tv.size() > 1) {
            checked_step = true;
            if (std::fabs((tv[1] - tv[0]) - kDt_ns) > 1e-3) {
                std::fprintf(stderr, "Sampling step %.4f ns != kDt_ns %.4f ns; aborting\n", tv[1]-tv[0], kDt_ns);
                std::exit(1);
            }
        }
        for(int c=32;c<64;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const auto& wf=*cw[c]; const int n=(int)wf.size();
            double s=0;int nb=0;
            for(int i=0;i<n;i++)if(tv[i]<kPre_ns){s+=wf[i];nb++;}
            if(nb<2)continue;
            const double ped=s/nb;
            double vmin=1e30;int imin=-1;
            for(int i=0;i<n;i++)if(tv[i]>=2500&&tv[i]<=6000&&wf[i]<vmin){vmin=wf[i];imin=i;}
            if(imin>=0 && (ped-vmin) > kPh19_mV) h_pt19.Fill(tv[imin]);
        }
        for(int c=0;c<32;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const auto& wf=*cw[c]; const int n=(int)wf.size();
            double s=0;int nb=0;
            for(int i=0;i<n;i++)if(tv[i]<kPre_ns){s+=wf[i];nb++;}
            if(nb<2)continue;
            const double ped=s/nb;
            double vmin=1e30;int imin=-1;
            for(int i=0;i<n;i++)if(tv[i]>=2500&&tv[i]<=8000&&wf[i]<vmin){vmin=wf[i];imin=i;}
            if(imin>=0 && (ped-vmin) > kPh18_mV) h_pt18.Fill(tv[imin]);
        }
    }
    // t_peak diagnóstico por estación
    rr.t_peak_19 = (h_pt19.GetEntries()>5) ? h_pt19.GetBinCenter(h_pt19.GetMaximumBin()) : 0.0;
    rr.t_peak_18 = (h_pt18.GetEntries()>5) ? h_pt18.GetBinCenter(h_pt18.GetMaximumBin()) : 0.0;
    // Ventana de señal fija [kMuonTlo, kMuonThi] para ambos hodoscopios
    rr.t_lo_19 = kMuonTlo; rr.t_hi_19 = kMuonThi;
    rr.t_lo_18 = kMuonTlo; rr.t_hi_18 = kMuonThi;

    // Histogramas de carga por plano: [-1000, 40000] mV·ns con bins de 100 mV·ns (410 bins)
    const double kQlo = -1000.0, kQhi = 40000.0;
    const int kQnbins = 410;
    TH1D hQ19x(("hQ19x_"+lbl).c_str(),"",kQnbins,kQlo,kQhi);
    TH1D hQ19y(("hQ19y_"+lbl).c_str(),"",kQnbins,kQlo,kQhi);
    TH1D hQ18x(("hQ18x_"+lbl).c_str(),"",kQnbins,kQlo,kQhi);
    TH1D hQ18y(("hQ18y_"+lbl).c_str(),"",kQnbins,kQlo,kQhi);

    // Variables para pedestal/noise by plane
    double sum_noise[4]={};long long n_noise[4]={};

    // Tag-and-probe per event: store Q19x per strip, hod2018 cluster
    struct EvTP{
        double Q19x[16]={},pt19x[16]{};
        double Q19y[16]={},pt19y[16]{};
        int k18x=-1,k18y=-1;
        double amp18x=0,amp18y=0;
        int ncls18x=0,ncls18y=0;
    };
    std::vector<EvTP> evs; evs.reserve(nEnt);
    std::vector<double> sig_Q_19x;
    std::vector<double> nopulse_Q_pl[4];

    // Pasada única
    for(Long64_t iev=0;iev<nEnt;iev++){
        t->GetEntry(iev);
        if(iev%30000==0)std::cout<<"    "<<iev<<"/"<<nEnt<<"\r"<<std::flush;
        if(!tw||tw->empty())continue;
        const auto& tv=*tw;

        EvTP ev;
        std::fill(ev.pt19x,ev.pt19x+16,-1.0);
        std::fill(ev.pt19y,ev.pt19y+16,-1.0);

        double maxQ[4]={-1e9,-1e9,-1e9,-1e9};
        std::vector<Hit18> h18x,h18y;
        double sum_ped[4]={};double sum_noise_ev[4]={};int n_ch[4]={};

        for(int c=0;c<kNCh;c++){
            if(!cw[c]||cw[c]->empty())continue;
            const bool is19=(c>=32);
            const int npost = is19 ? kNpost19 : kNpost18;
            const double tlo=is19?rr.t_lo_19:rr.t_lo_18;
            const double thi=is19?rr.t_hi_19:rr.t_hi_18;
            const int pl = ch_to_plane(c);

            // Baseline
            double s=0,s2=0;int nb=0;
            for(int i=0;i<(int)tv.size();i++)if(tv[i]<kPre_ns){s+=(*cw[c])[i];s2+=(*cw[c])[i]*(*cw[c])[i];nb++;}
            if(nb<2)continue;
            const double ped=s/nb; const double var=s2/nb-ped*ped;
            const double noise=(var>0)?std::sqrt(var):0;
            sum_ped[pl]+=ped; sum_noise_ev[pl]+=noise; n_ch[pl]++;

            // Q
            const bool use_loc = is19 ? kQLocalPed19 : kQLocalPed18;
            QResult qr=get_Q(*cw[c],tv,tlo,thi,npost,use_loc);
            const double Q=qr.Q;

            if(Q>maxQ[pl])maxQ[pl]=Q;

            // Estimador de carga robusto: pulsos en hod2019X con ph > kPh19_mV en [kMuonTlo, kMuonThi]
            if(pl==2){
                double vmin=ped; int imin=-1;
                for(int i=0;i<(int)tv.size();i++)if(tv[i]>=tlo&&tv[i]<=thi&&(*cw[c])[i]<vmin){vmin=(*cw[c])[i];imin=i;}
                if(imin>=0 && (ped-vmin) > kPh19_mV && qr.pt >= kMuonTlo && qr.pt <= kMuonThi){
                    sig_Q_19x.push_back(Q);
                }
            }

            // F1: Q noise over channel-events with ph < 5 mV in [2500, 8000] ns
            double vmin_np = 1e30;
            for (int i = 0; i < (int)tv.size(); ++i) {
                if (tv[i] >= 2500.0 && tv[i] <= 8000.0 && (*cw[c])[i] < vmin_np)
                    vmin_np = (*cw[c])[i];
            }
            if ((ped - vmin_np) < 5.0) {
                const int i0 = 262; // sample at 4200 ns
                double ped_Q_np = ped;
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
                nopulse_Q_pl[pl].push_back(Q_noise);
            }

            // Fill histo planes
            if(pl==2){hQ19x.Fill(Q); int k=inv19x[c];if(k>=0){ev.Q19x[k]=Q;ev.pt19x[k]=qr.pt;}}
            else if(pl==3){hQ19y.Fill(Q); int k=inv19y[c];if(k>=0){ev.Q19y[k]=Q;ev.pt19y[k]=qr.pt;}}
            else if(pl==0){hQ18x.Fill(Q);
                if(c<=15){double vmin=ped;
                    for(int i=0;i<(int)tv.size();i++)if(tv[i]>=tlo&&tv[i]<=thi&&(*cw[c])[i]<vmin)vmin=(*cw[c])[i];
                    h18x.push_back({c,(double)(ped-vmin)});}}
            else{hQ18y.Fill(Q);
                if(c<=31){double vmin=ped;
                    for(int i=0;i<(int)tv.size();i++)if(tv[i]>=tlo&&tv[i]<=thi&&(*cw[c])[i]<vmin)vmin=(*cw[c])[i];
                    h18y.push_back({c-16,(double)(ped-vmin)});}}

            // Noise per plane
            sum_noise[pl]+=noise; n_noise[pl]++;
        }

        // hod2018 clusters
        {auto cls=clust18(h18x,10.0);ev.ncls18x=(int)cls.size();
         if(!cls.empty()){auto&b=*std::max_element(cls.begin(),cls.end(),[](const Cl18&a,const Cl18&b){return a.amp<b.amp;});
             ev.k18x=(int)std::round(b.pos/kPitch+7.5); ev.amp18x=b.amp;
             if(ev.k18x<0)ev.k18x=0;if(ev.k18x>15)ev.k18x=15;}}
        {auto cls=clust18(h18y,10.0);ev.ncls18y=(int)cls.size();
         if(!cls.empty()){auto&b=*std::max_element(cls.begin(),cls.end(),[](const Cl18&a,const Cl18&b){return a.amp<b.amp;});
             ev.k18y=(int)std::round(b.pos/kPitch+7.5); ev.amp18y=b.amp;
             if(ev.k18y<0)ev.k18y=0;if(ev.k18y>15)ev.k18y=15;}}

        evs.push_back(std::move(ev));
    }
    std::cout<<"    "<<nEnt<<"/"<<nEnt<<" ✓\n";
    fin->Close();

    // Estimador robusto hod2019X
    rr.n_sig_19x = (long long)sig_Q_19x.size();
    if(rr.n_sig_19x >= 50){
        std::sort(sig_Q_19x.begin(), sig_Q_19x.end());
        auto get_pct = [&](double pct) {
            double idx = pct * (sig_Q_19x.size() - 1);
            size_t i0 = static_cast<size_t>(std::floor(idx));
            size_t i1 = static_cast<size_t>(std::ceil(idx));
            double frac = idx - i0;
            return sig_Q_19x[i0] * (1.0 - frac) + sig_Q_19x[i1] * frac;
        };
        rr.med_Q_19x = get_pct(0.50);
        rr.q10_Q_19x = get_pct(0.10);
        rr.q90_Q_19x = get_pct(0.90);
    } else {
        rr.med_Q_19x = 0; rr.q10_Q_19x = 0; rr.q90_Q_19x = 0;
    }

    // Normaliza ruido por plano
    for(int pl=0;pl<4;pl++){
        if(n_noise[pl]>0){
            rr.noise_mean[pl]=sum_noise[pl]/n_noise[pl];
        }
        // F1: Pedestal noise of Q (RMS over channel-events with ph < 5 mV in [2500,8000] ns)
        if(!nopulse_Q_pl[pl].empty()){
            double sum=0, sum2=0;
            for(double q:nopulse_Q_pl[pl]) sum += q;
            double mean = sum / nopulse_Q_pl[pl].size();
            for(double q:nopulse_Q_pl[pl]) sum2 += (q - mean) * (q - mean);
            rr.sigmaQ[pl] = std::sqrt(sum2 / nopulse_Q_pl[pl].size());
        }
    }

    // Landau fit hod2019X (si hay datos)
    {
        TH1D* h = &hQ19x;
        if(h->Integral() > 10){
            // Step 4: Valle en rebinned histogram (500 mV·ns)
            if(rr.n_sig_19x < 50){
                rr.has_valley_19x = false;
                rr.umbral_Q[2] = 0;
                rr.fail_reason_19x = "N_sig < 50";
            } else {
                TH1D* h500 = static_cast<TH1D*>(h->Rebin(5, Form("%s_rebin500", h->GetName())));
                rr.v_xlo_19x = 5.0 * rr.sigmaQ[2];
                rr.v_xhi_19x = rr.med_Q_19x;
                int blo = h500->FindBin(rr.v_xlo_19x);
                int bhi = h500->FindBin(rr.v_xhi_19x);
                rr.blo_19x = blo; rr.bhi_19x = bhi;
                int min_b = -1;
                double min_v = 1e30;
                if(blo >= 1 && bhi <= h500->GetNbinsX() && bhi >= blo){
                    for(int b = blo; b <= bhi; b++){
                        double c = h500->GetBinContent(b);
                        if(c < min_v){
                            min_v = c;
                            min_b = b;
                        }
                    }
                }
                rr.min_b_19x = min_b;
                rr.min_v_19x = min_v;
                rr.c_lo_19x = (blo >= 1 && blo <= h500->GetNbinsX()) ? h500->GetBinContent(blo) : 0.0;
                rr.c_hi_19x = (bhi >= 1 && bhi <= h500->GetNbinsX()) ? h500->GetBinContent(bhi) : 0.0;

                const int n_bins = (bhi >= blo) ? (bhi - blo + 1) : 0;
                if(n_bins < 8){
                    rr.fail_reason_19x = "region span < 8 bins";
                } else if(rr.c_lo_19x < 20.0 || rr.c_hi_19x < 20.0){
                    rr.fail_reason_19x = "end content < 20";
                } else if(min_b < blo + 3 || min_b > bhi - 3){
                    rr.fail_reason_19x = "edge guard";
                } else if(min_v >= 0.5 * rr.c_lo_19x || min_v >= 0.5 * rr.c_hi_19x){
                    rr.fail_reason_19x = "min_v >= 0.5*c_end";
                } else {
                    rr.fail_reason_19x = "";
                }

                if(rr.fail_reason_19x.empty()){
                    rr.has_valley_19x = true;
                    rr.umbral_Q[2] = h500->GetBinCenter(min_b);
                } else {
                    rr.has_valley_19x = false;
                    rr.umbral_Q[2] = 0;
                }
                delete h500;
            }

        }
    }

    // Tag-and-probe (hod2019X, global shuffle + ventana temporal fija [kMuonTlo, kMuonThi])
    const size_t nev=evs.size();
    const double t_muon_lo=kMuonTlo;
    const double t_muon_hi=kMuonThi;

    long long n_obs=0,n_acc=0;
    const double Q_thr=kQcut19;  // from hodo_consts.h

    // Local shuffle (semilla fija por run)
    std::vector<size_t> local_shuf(nev);
    std::iota(local_shuf.begin(),local_shuf.end(),0);
    {std::mt19937 rng(42+std::hash<std::string>{}(lbl));
     std::shuffle(local_shuf.begin(),local_shuf.end(),rng);}

    for(size_t i=0;i<nev;i++){
        if(evs[i].ncls18x!=1||evs[i].amp18x<10)continue;
        const int kp=evs[i].k18x;
        if(kp<0||kp>15)continue;
        ++rr.n_tag_x;
        // Q_same
        const double Qs=evs[i].Q19x[kp];
        const bool ts_ok=(evs[i].pt19x[kp]>=t_muon_lo&&evs[i].pt19x[kp]<=t_muon_hi);
        const double Qsame=ts_ok?Qs:0.0;
        // Q_mixed
        const size_t j=local_shuf[i];
        const double Qm=evs[j].Q19x[kp];
        const bool tm_ok=(evs[j].pt19x[kp]>=t_muon_lo&&evs[j].pt19x[kp]<=t_muon_hi);
        const double Qmixed=tm_ok?Qm:0.0;
        if(Qsame>Q_thr)++n_obs;
        if(Qmixed>Q_thr)++n_acc;

        // Q_norm = Q / sigmaQ19x para histograma combinado
        if(rr.sigmaQ[2]>0){
            const double Qn_same=(ts_ok&&Qs>0)?Qs/rr.sigmaQ[2]:0.0;
            const double Qn_mixed=(tm_ok&&Qm>0)?Qm/rr.sigmaQ[2]:0.0;
            if(h_Qcomb_same)h_Qcomb_same->Fill(Qn_same);
            if(h_Qcomb_mixed)h_Qcomb_mixed->Fill(Qn_mixed);
        }
    }
    rr.n_obs_x = n_obs;
    rr.n_acc_x = n_acc;
    const double Nt=(double)rr.n_tag_x;
    if(Nt>=2){
        rr.eps_true_x=(n_obs-n_acc)/Nt;
        const double so=std::sqrt((double)n_obs/Nt*(1-(double)n_obs/Nt)/Nt);
        const double sa=std::sqrt((double)n_acc/Nt*(1-(double)n_acc/Nt)/Nt);
        rr.sigma_eps_x=std::sqrt(so*so+sa*sa);
        rr.sig_x=(rr.sigma_eps_x>0)?rr.eps_true_x/rr.sigma_eps_x:0;
    }
    return rr;
}

//==============================================================================
// main
//==============================================================================
int main(int argc,char**argv){
    std::time_t now=std::time(nullptr);char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("multirun_"+std::string(ds));
    fs::create_directories(outdir);
    std::cout<<"=== multirun_analysis ===\nOut: "<<outdir<<"\n\n";
    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    // Lista de archivos (v1=sin sufijo, v2..v10)
    std::vector<std::pair<std::string,std::string>> files={
        {"lectura_vx1740_ext_nim_1M.root","v1"},
        {"lectura_vx1740_ext_nim_1M_v2.root","v2"},
        {"lectura_vx1740_ext_nim_1M_v3.root","v3"},
        {"lectura_vx1740_ext_nim_1M_v4.root","v4"},
        {"lectura_vx1740_ext_nim_1M_v5.root","v5"},
        {"lectura_vx1740_ext_nim_1M_v6.root","v6"},
        {"lectura_vx1740_ext_nim_1M_v7.root","v7"},
        {"lectura_vx1740_ext_nim_1M_v8.root","v8"},
        {"lectura_vx1740_ext_nim_1M_v9.root","v9"},
        {"lectura_vx1740_ext_nim_1M_v10.root","v10"}
    };
    // Filtra archivos que no existen
    std::vector<std::pair<std::string,std::string>> valid_files;
    for(const auto& fl:files){
        if(fs::exists(fl.first)) valid_files.push_back(fl);
        else std::cout<<"  MISSING: "<<fl.first<<"\n";}
    std::cout<<"Archivos encontrados: "<<valid_files.size()<<"\n\n";

    // Histogramas combinados (Q/sigma, para acumular todos los runs)
    auto h_same_all=std::make_unique<TH1D>("h_same_all","H_{same} combinado;Q_{norm}=#sigma_Q;Events",100,-2,30);
    auto h_mixed_all=std::make_unique<TH1D>("h_mixed_all","H_{mixed} combinado;Q_{norm}=#sigma_Q;Events",100,-2,30);

    // Procesa cada run
    std::vector<RunResult> results;
    size_t global_off=0;
    for(const auto& fl:valid_files){
        std::cout<<"Processing "<<fl.second<<": "<<fl.first<<"\n";
        auto rr=process_run(fl.first,fl.second,h_same_all.get(),h_mixed_all.get(),{},global_off);
        global_off+=rr.n_events;
        results.push_back(rr);
    }

    // Reproducibilidad v8 vs v9 (si ambos existen)
    {const RunResult* rv8=nullptr,*rv9=nullptr,*rv10=nullptr;
     for(const auto&r:results){if(r.label=="v8")rv8=&r;if(r.label=="v9")rv9=&r;if(r.label=="v10")rv10=&r;}
     if(rv8&&rv9&&!rv8->lan_at_limit&&!rv9->lan_at_limit&&rv8->mpvQ[2]>0&&rv9->mpvQ[2]>0){
         const double delta=(rv9->mpvQ[2]-rv8->mpvQ[2])/rv8->mpvQ[2];
         std::cout<<"\nReproducibilidad v8 vs v9: ΔMPV_Q/MPV_Q = "<<100*delta<<"%";
         if(std::abs(delta)>0.20)std::cout<<" [FLAG: >20%!]";else std::cout<<" [OK]";std::cout<<"\n";
     }
     if(rv8&&rv9&&rv8->med_Q_19x>0&&rv9->med_Q_19x>0){
         const double delta=(rv9->med_Q_19x-rv8->med_Q_19x)/rv8->med_Q_19x;
         std::cout<<"Reproducibilidad v8 vs v9 (mediana Q): ΔMed_Q/Med_Q = "<<100*delta<<"%\n";
     }
     if(rv8&&rv10&&rv8->med_Q_19x>0&&rv10->med_Q_19x>0){
         const double delta=(rv10->med_Q_19x-rv8->med_Q_19x)/rv8->med_Q_19x;
         std::cout<<"Reproducibilidad v8 vs v10 (mediana Q): ΔMed_Q/Med_Q = "<<100*delta<<"%\n";}}

    // ── Figuras ───────────────────────────────────────────────────────────────
    const int nRuns=static_cast<int>(results.size());
    std::vector<double> run_nums,t_peaks,mpvQ_19x,sigQ_19x,snrQ_19x,eps_true,eps_err;
    for(int i=0;i<nRuns;i++){
        const auto&r=results[i];
        run_nums.push_back(i+1);
        t_peaks.push_back(r.t_peak_19);
        mpvQ_19x.push_back(r.med_Q_19x);
        sigQ_19x.push_back(r.sigmaQ[2]);
        snrQ_19x.push_back(r.snrQ[2]);
        eps_true.push_back(r.eps_true_x*100);
        eps_err.push_back(r.sigma_eps_x*100);
    }

    // SNR y Mediana Q vs run
    {TCanvas cc("c_snr","",1000,680);cc.SetLeftMargin(0.12);cc.SetBottomMargin(0.13);
     TGraphErrors gmpv(nRuns,run_nums.data(),mpvQ_19x.data(),
         std::vector<double>(nRuns,0.2).data(),std::vector<double>(nRuns,0).data());
     gmpv.SetTitle("Mediana Q hod2019X vs run (proxy de voltaje);Run;Mediana Q [mV#upointns]");
     gmpv.SetMarkerStyle(20);gmpv.SetMarkerColor(kBlue+1);gmpv.SetLineColor(kBlue+1);gmpv.SetLineWidth(2);
     gmpv.Draw("APL");
     cc.SaveAs((outdir/"mpv_and_pedestal_vs_run.png").string().c_str());}

    {TCanvas cc("c_ep","",1000,680);cc.SetLeftMargin(0.12);cc.SetBottomMargin(0.13);
     TGraphErrors gep(nRuns,run_nums.data(),eps_true.data(),
         std::vector<double>(nRuns,0.2).data(),eps_err.data());
     gep.SetTitle("ε_{true} (tag-and-probe) vs run;Run;ε_{true} [%]");
     gep.SetMarkerStyle(20);gep.SetMarkerColor(kBlue+1);gep.SetLineColor(kBlue+1);gep.SetLineWidth(2);
     gep.Draw("APL");
     TLine l0(0.5,0,nRuns+0.5,0);l0.SetLineColor(kRed);l0.SetLineStyle(2);l0.Draw();
     cc.SaveAs((outdir/"tagprobe_efficiency_vs_voltage.png").string().c_str());}

    // Histograma combinado same vs mixed
    {TCanvas cc("c_comb","",1000,700);cc.SetLogy();cc.SetLeftMargin(0.12);cc.SetBottomMargin(0.13);
     h_same_all->SetLineColor(kBlue+1);h_same_all->SetLineWidth(2);h_same_all->SetStats(0);
     h_mixed_all->SetLineColor(kOrange+7);h_mixed_all->SetLineWidth(2);
     h_same_all->SetTitle("Tag-and-probe COMBINADO (todos los runs);Q_{norm} = Q/#sigma_Q;Events");
     h_same_all->Draw("HIST");h_mixed_all->Draw("HIST SAME");
     TLegend lg(0.55,0.75,0.88,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
     lg.AddEntry(h_same_all.get(),"H_{same}","l");lg.AddEntry(h_mixed_all.get(),"H_{mixed}","l");lg.Draw();
     cc.SaveAs((outdir/"tagprobe_combined_same_vs_mixed.png").string().c_str());}

    // Significancia global combinada (eps_global = (sum N_same - sum N_mixed) / sum N_tag)
    long long sum_ntag = 0, sum_nsame = 0, sum_nmixed = 0;
    for(const auto& r : results){
        if(r.n_tag_x > 0){
            sum_ntag += r.n_tag_x;
            sum_nsame += r.n_obs_x;
            sum_nmixed += r.n_acc_x;
        }
    }
    const double eps_true_g = (sum_ntag > 0) ? (double)(sum_nsame - sum_nmixed) / sum_ntag : 0.0;
    const double sigma_eps_g = (sum_ntag > 0) ? std::sqrt((double)(sum_nsame + sum_nmixed)) / sum_ntag : 0.0;
    const double sig_global = (sigma_eps_g > 0) ? eps_true_g / sigma_eps_g : 0.0;

    // ── CSV por run ───────────────────────────────────────────────────────────
    {std::ofstream fc(outdir/"per_run_parameters.csv");
     fc<<"run,n_events,t_peak_18_ns,t_peak_19_ns,n_sig_19x,med_Q_19x,q10_Q_19x,q90_Q_19x,sigmaQ_19x_mVns,n_tag_x,n_same_x,n_mixed_x,eps_true_x_pct,sigma_eps_x_pct,sig_x\n";
     for(const auto&r:results){
         std::string med_str = (r.n_sig_19x >= 50) ? [&](){ char b[32]; std::snprintf(b, sizeof(b), "%.1f", r.med_Q_19x); return std::string(b); }() : "n/a";
         std::string q10_str = (r.n_sig_19x >= 50) ? [&](){ char b[32]; std::snprintf(b, sizeof(b), "%.1f", r.q10_Q_19x); return std::string(b); }() : "n/a";
         std::string q90_str = (r.n_sig_19x >= 50) ? [&](){ char b[32]; std::snprintf(b, sizeof(b), "%.1f", r.q90_Q_19x); return std::string(b); }() : "n/a";
         fc<<r.label<<','<<r.n_events<<','<<r.t_peak_18<<','<<r.t_peak_19<<','
           <<r.n_sig_19x<<','<<med_str<<','<<q10_str<<','<<q90_str<<','
           <<r.sigmaQ[2]<<','<<r.n_tag_x<<','<<r.n_obs_x<<','<<r.n_acc_x<<','
           <<r.eps_true_x*100<<','<<r.sigma_eps_x*100<<','<<r.sig_x<<'\n';
     }}

    // ── Summary ───────────────────────────────────────────────────────────────
    const double max_sig=[&]()->double{double m=0;for(const auto&r:results)m=std::max(m,std::abs(r.sig_x));return m;}();

    auto print_summ=[&](std::ostream&os){
        os<<"=== Análisis Multi-Run: barrido de voltaje MPPC ===\n\n";
        os<<"Parámetros de ventana muon: ["<<kMuonTlo<<", "<<kMuonThi<<"] ns (ambas estaciones)\n";
        os<<"Umbrales de pulso diagnóstico: hod2018 ph > "<<kPh18_mV<<" mV, hod2019 ph > "<<kPh19_mV<<" mV\n";
        os<<"Med_Q: (Q of pulses above the ph threshold; not a gain estimator)\n";
        os<<"sigmaQ19: σ_Q noise [mV·ns] (no-pulse channels)\n\n";
        os<<std::setw(6)<<std::left<<"Run"<<std::setw(9)<<"N"<<std::setw(8)<<"tpk18"<<std::setw(8)<<"tpk19"
          <<std::setw(8)<<"N_sig"<<std::setw(10)<<"Med_Q"<<std::setw(9)<<"Q_10"<<std::setw(9)<<"Q_90"
          <<std::setw(10)<<"sigmaQ19"
          <<std::setw(8)<<"N_tag"<<std::setw(8)<<"N_same"<<std::setw(8)<<"N_mixed"
          <<std::setw(12)<<"eps_true[%]"<<std::setw(8)<<"sig\n";
        os<<"-----------------------------------------------------------------------------------------------------------------\n";
        for(const auto&r:results){
            std::string med_str = (r.n_sig_19x >= 50) ? [&](){ char b[32]; std::snprintf(b, sizeof(b), "%.1f", r.med_Q_19x); return std::string(b); }() : "n/a";
            std::string q10_str = (r.n_sig_19x >= 50) ? [&](){ char b[32]; std::snprintf(b, sizeof(b), "%.1f", r.q10_Q_19x); return std::string(b); }() : "n/a";
            std::string q90_str = (r.n_sig_19x >= 50) ? [&](){ char b[32]; std::snprintf(b, sizeof(b), "%.1f", r.q90_Q_19x); return std::string(b); }() : "n/a";

            os<<std::setw(6)<<r.label<<std::setw(9)<<r.n_events
              <<std::setw(8)<<std::fixed<<std::setprecision(0)<<r.t_peak_18
              <<std::setw(8)<<r.t_peak_19
              <<std::setw(8)<<r.n_sig_19x
              <<std::setw(10)<<med_str
              <<std::setw(9)<<q10_str
              <<std::setw(9)<<q90_str
              <<std::setw(10)<<std::setprecision(1)<<r.sigmaQ[2]
              <<std::setw(8)<<r.n_tag_x
              <<std::setw(8)<<r.n_obs_x
              <<std::setw(8)<<r.n_acc_x
              <<std::setw(12)<<std::setprecision(4)<<r.eps_true_x*100
              <<std::setw(8)<<std::setprecision(2)<<r.sig_x<<"\n";
        }
        os<<"\nUmbral Q natural hod2019X (valle en [5*sigma_Q, median] mV·ns):\n";
        for(const auto& r : results){
            os<<"  "<<std::setw(6)<<std::left<<r.label<<": ";
            if(r.n_sig_19x < 50){
                os<<"no valley (N_sig="<<r.n_sig_19x<<" < 50)\n";
            } else {
                if(r.has_valley_19x && r.umbral_Q[2] > 0) os<<std::fixed<<std::setprecision(1)<<r.umbral_Q[2]<<" mV·ns";
                else os<<"no valley";
                os<<" (region=["<<std::fixed<<std::setprecision(1)<<r.v_xlo_19x<<", "<<r.v_xhi_19x<<"] mV·ns"
                  <<", bins=["<<r.blo_19x<<", "<<r.bhi_19x<<"], min_b="<<r.min_b_19x
                  <<", min_v="<<std::setprecision(0)<<r.min_v_19x
                  <<", c_lo="<<r.c_lo_19x<<", c_hi="<<r.c_hi_19x;
                if(!r.fail_reason_19x.empty()) os<<", "<<r.fail_reason_19x;
                os<<")\n";
            }
        }
        os<<"\nSignificancia global combinada (eps_global = (sum N_same - sum N_mixed) / sum N_tag, sigma = sqrt(sum N_same + sum N_mixed) / sum N_tag): "<<std::setprecision(2)<<sig_global<<"σ\n";
        os<<"ε_true_global = "<<std::setprecision(4)<<eps_true_g*100<<"% ± "<<sigma_eps_g*100<<"%\n\n";
        os<<"  cota superior ε_true < "<<std::setprecision(2)
          <<(std::abs(eps_true_g)+2*sigma_eps_g)*100<<"% (95% CL)\n";
    };
    print_summ(std::cout);
    {std::ofstream sf(outdir/"summary_multirun.txt");print_summ(sf);}

    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
