#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>
#include <Math/ProbFuncMathCore.h>

#include <algorithm>
#include <array>
#include <cmath>
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
constexpr int NC = 64, NR = 9;
constexpr double EXPOSURE_S = 3.2e-6;
constexpr std::array<int,16> X19={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
constexpr std::array<int,16> Y19={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};
const std::array<std::string,NR> FILES={
 "lectura_vx1740_ext_nim_1M.root","lectura_vx1740_ext_nim_1M_v2.root","lectura_vx1740_ext_nim_1M_v3.root",
 "lectura_vx1740_ext_nim_1M_v4.root","lectura_vx1740_ext_nim_1M_v5.root","lectura_vx1740_ext_nim_1M_v6.root",
 "lectura_vx1740_ext_nim_1M_v7.root","lectura_vx1740_ext_nim_1M_v8.root","lectura_vx1740_ext_nim_1M_v9.root"};
const std::array<std::string,NR> RUNS={"v1","v2","v3","v4","v5","v6","v7","v8","v9"};
struct Win { int a=0,b=-1; };
struct Dcr { std::array<long long,4> down{},up{}; long long entries=0; };
struct ChData { long long nin=0,nacc=0; std::vector<double> ph; };
struct Fit { double mpv=NAN,err=NAN,width=NAN; int status=-1; };

static Win window(const std::vector<float>& t,double lo,double hi) {
 auto a=std::lower_bound(t.begin(),t.end(),static_cast<float>(lo));
 auto b=std::upper_bound(t.begin(),t.end(),static_cast<float>(hi));
 if(a==t.end()||b==t.begin()||a>=b)return{};
 return{int(a-t.begin()),int(b-t.begin())-1};
}
static double median(std::vector<double>& v) {
 auto m=v.begin()+v.size()/2; std::nth_element(v.begin(),m,v.end()); double x=*m;
 if(v.size()%2==0)x=(x+*std::max_element(v.begin(),m))/2.; return x;
}
static std::array<long long,4> excursions(const std::vector<float>& wf,Win w,double ped,
                                          const std::array<double,4>& cuts,bool negative) {
 std::array<long long,4> counts{}; std::array<int,4> run{};
 for(int i=w.a;i<=w.b;++i) {
  double amp=negative?ped-wf[i]:wf[i]-ped;
  for(int j=0;j<4;++j) {
   if(amp>cuts[j])++run[j];
   else {if(run[j]>=2)++counts[j];run[j]=0;}
  }
 }
 for(int j=0;j<4;++j)if(run[j]>=2)++counts[j];
 return counts;
}
static double peak(const std::vector<float>& wf,Win w,double ped) {
 double p=0;for(int i=w.a;i<=w.b;++i)p=std::max(p,ped-double(wf[i]));return p;
}
static std::pair<std::string,int> plane_k(int ch) {
 if(ch<16)return{"hod2018_X",ch}; if(ch<32)return{"hod2018_Y",ch-16};
 for(int k=0;k<16;++k){if(X19[k]==ch)return{"hod2019_X",k};if(Y19[k]==ch)return{"hod2019_Y",k};}
 return{"unknown",-1};
}
static double poisson_excess_p(long long n,double lambda) {
 if(n<=0)return 1.; if(lambda<=0)return 0.;
 return ROOT::Math::poisson_cdf_c(static_cast<unsigned int>(n-1),lambda);
}
static Fit fit_full(const std::vector<double>& values,const std::string& name) {
 Fit out; if(values.size()<100)return out;
 double xmax=*std::max_element(values.begin(),values.end()); xmax=std::max(20.,xmax*1.02);
 TH1D h((name+"_h").c_str(),"",400,0,xmax); h.SetDirectory(nullptr);
 for(double x:values)h.Fill(x);
 TF1 f((name+"_f").c_str(),"gaus(0)+landau(3)",0,xmax);
 double peakx=h.GetXaxis()->GetBinCenter(h.GetMaximumBin());
 f.SetParameters(std::max(1.,h.GetMaximum()),std::clamp(peakx,0.,xmax),2.,
                 std::max(1.,h.GetMaximum()*.1),std::clamp(peakx+10.,0.,xmax),5.);
 f.SetParLimits(0,0.,h.GetEntries()*2.); f.SetParLimits(2,.1,30.);
 f.SetParLimits(3,0.,h.GetEntries()*2.); f.SetParLimits(5,.1,xmax);
 TFitResultPtr r=h.Fit(&f,"LQ0S"); out.status=int(r);
 if(r.Get()&&r->IsValid()){
  constexpr double landau_peak_shift=0.22278298;
  out.width=f.GetParameter(5);
  out.mpv=f.GetParameter(4)-landau_peak_shift*out.width;
  auto cov=r->GetCovarianceMatrix();
  double variance=cov(4,4)+landau_peak_shift*landau_peak_shift*cov(5,5)
                  -2.*landau_peak_shift*cov(4,5);
  out.err=variance>=0?std::sqrt(variance):NAN;
 }
 return out;
}

int main(int argc,char**argv) {
 fs::path outdir=argc>1?argv[1]:"outputs/channel_health_20261016";fs::create_directories(outdir);
 gROOT->SetBatch(kTRUE);
 std::array<std::array<Dcr,NC>,NR> dcr{};std::array<ChData,NC> chdata{};
 std::array<long long,NR> entries{};std::array<std::array<long long,2>,2> mcstrip{};
 constexpr int MCN=5000000; constexpr double zstart=160,z18=0,z19=98.4,zT1=124.1,zT2=134.1,half=16.5;
 // Deterministic cos^2(theta) ray MC; trigger acceptance and strip bins are counted together.
 std::mt19937_64 rng(42);std::uniform_real_distribution<double> u(0,1),xy(-180,180),phi(0,2*M_PI);
 long long ntrigger=0;std::array<std::array<long long,16>,4> strip_counts{};
 for(int n=0;n<MCN;++n){
  double U=u(rng),th=std::acos(std::cbrt(U)),ph=phi(rng),dx=std::sin(th)*std::cos(ph),dy=std::sin(th)*std::sin(ph),dz=-std::cos(th),x0=xy(rng),y0=xy(rng);
  auto prop=[&](double z){double t=(z-zstart)/dz;return std::pair<double,double>{x0+dx*t,y0+dy*t};};
  auto [x1,y1]=prop(zT1);auto[x2,y2]=prop(zT2);
  if(x1< -117.5||x1>32.5||y1< -117.5||y1>32.5||x2< -32.5||x2>117.5||y2< -32.5||y2>117.5)continue;
  ++ntrigger;
  for(int h=0;h<2;++h){auto [x,y]=prop(h?z19:z18);if(std::abs(x)>half||std::abs(y)>half)continue;
   int kx=std::clamp(int((x+half)/33.*16),0,15),ky=std::clamp(int((y+half)/33.*16),0,15);
   ++strip_counts[h*2][kx];++strip_counts[h*2+1][ky];
  }
 }
 if(!ntrigger){std::cerr<<"MC produced no triggers\n";return 3;}
 std::cout<<"MC generated="<<MCN<<" triggers="<<ntrigger<<"\n";
 std::ofstream mc(outdir/"cos2_strip_acceptance.csv");mc<<"hodoscope,axis,k,mc_trigger_count,n_mc_triggers,fraction_per_trigger\n"<<std::setprecision(10);
 for(int h=0;h<2;++h)for(int ax=0;ax<2;++ax)for(int k=0;k<16;++k){int p=h*2+ax;mc<<(h?"hod2019":"hod2018")<<','<<(ax?"Y":"X")<<','<<k<<','<<strip_counts[p][k]<<','<<ntrigger<<','<<double(strip_counts[p][k])/ntrigger<<'\n';}

 std::ofstream runsout(outdir/"dcr_thresholds_v1_v9.csv");
 runsout<<"run,ch,hodoscope,axis,k,entries,sigma4_down,sigma4_up,sigma4_down_Hz,sigma4_up_Hz,cut5_down,cut5_up,cut5_down_Hz,cut5_up_Hz,cut10_down,cut10_up,cut10_down_Hz,cut10_up_Hz,cut20_down,cut20_up,cut20_down_Hz,cut20_up_Hz\n";
 for(int r=0;r<NR;++r){
  std::string resolved = hodo::resolveRawFile(FILES[r]);
  std::unique_ptr<TFile> f(TFile::Open(resolved.c_str(),"READ"));if(!f||f->IsZombie()){std::cerr<<"Cannot open "<<resolved<<'\n';return 2;}
  auto*t=dynamic_cast<TTree*>(f->Get("data_dig"));if(!t){std::cerr<<"Missing data_dig\n";return 2;}entries[r]=t->GetEntries();
  t->SetBranchStatus("*",0);std::array<std::vector<float>*,NC>wf{};std::vector<float>*time=nullptr;
  for(int c=0;c<NC;++c){std::string b="Ch"+std::to_string(c);t->SetBranchStatus(b.c_str(),1);if(t->SetBranchAddress(b.c_str(),&wf[c])<0)return 2;}
  t->SetBranchStatus("time",1);if(t->SetBranchAddress("time",&time)<0)return 2;
  Win pre,mu,acc;
  std::vector<double> base,dev;
  for(long long e=0;e<entries[r];++e){t->GetEntry(e);if(!time||time->empty())continue;if(e==0){pre=window(*time,200,3400);mu=window(*time,4000,4400);acc=window(*time,2000,2400);}
   for(int c=0;c<NC;++c){if(!wf[c]||wf[c]->empty())continue;auto&v=*wf[c];Win w{pre.a,std::min(pre.b,int(v.size())-1)};if(w.b<w.a)continue;
    base.clear();for(int s=w.a;s<=w.b;++s)base.push_back(v[s]);double ped=median(base),mean=std::accumulate(base.begin(),base.end(),0.)/base.size();
    dev.clear();double ss=0;for(double q:base){dev.push_back(std::abs(q-ped));ss+=(q-mean)*(q-mean);}double sigmad=1.4826*median(dev),rms=std::sqrt(ss/std::max<size_t>(1,base.size()-1));
    double sig=std::max({sigmad,rms,1.0});std::array<double,4>cuts={4*sig,5,10,20};auto dn=excursions(v,w,ped,cuts,true),up=excursions(v,w,ped,cuts,false);
    for(int j=0;j<4;++j){dcr[r][c].down[j]+=dn[j];dcr[r][c].up[j]+=up[j];}++dcr[r][c].entries;
    if(r==7){double threshold=c<32?50.:35.;Win mw{mu.a,std::min(mu.b,int(v.size())-1)},aw{acc.a,std::min(acc.b,int(v.size())-1)};
     double p=peak(v,mw,ped),a=peak(v,aw,ped);if(p>threshold)++chdata[c].nin;if(a>threshold)++chdata[c].nacc;chdata[c].ph.push_back(p);}
   }
   if(e&&e%20000==0)std::cout<<RUNS[r]<<" "<<e<<"/"<<entries[r]<<"\r"<<std::flush;
  }
  std::cout<<RUNS[r]<<" entries="<<entries[r]<<"\n";
  for(int c=0;c<NC;++c){auto[plane,k]=plane_k(c);runsout<<RUNS[r]<<','<<c<<','<<plane.substr(0,7)<<','<<plane.back()<<','<<k<<','<<dcr[r][c].entries;
   for(int j=0;j<4;++j){double den=dcr[r][c].entries*EXPOSURE_S;double rd=den?dcr[r][c].down[j]/den:0,ru=den?dcr[r][c].up[j]/den:0;
    runsout<<','<<dcr[r][c].down[j]<<','<<dcr[r][c].up[j]<<','<<rd<<','<<ru;}
   runsout<<'\n';}
 }

 std::ofstream resp(outdir/"channel_response_v8.csv");
 resp<<"ch,hodoscope,axis,k,N_in,N_acc,p_value_poisson_upper,N_exp_k,eps_k,ph_entries,landau_mpv_full,landau_mpv_error,landau_width,fit_status\n"<<std::setprecision(10);
 std::array<std::vector<double>,2> pooled;std::array<Fit,NC> fits{};
 for(int c=0;c<NC;++c){auto[plane,k]=plane_k(c);int h=c<32?0:1,ax=plane.back()=='Y'?1:0;for(double p:chdata[c].ph)pooled[h].push_back(p);
  fits[c]=fit_full(chdata[c].ph,"full_ch"+std::to_string(c));double fraction=double(strip_counts[h*2+ax][k])/ntrigger;
  double nexp=fraction*entries[7],eps=nexp>0?(chdata[c].nin-chdata[c].nacc)/nexp:NAN;
  resp<<c<<','<<(h?"hod2019":"hod2018")<<','<<(ax?"Y":"X")<<','<<k<<','<<chdata[c].nin<<','<<chdata[c].nacc<<','<<poisson_excess_p(chdata[c].nin,chdata[c].nacc)<<','<<nexp<<','<<eps<<','<<chdata[c].ph.size()<<','<<fits[c].mpv<<','<<fits[c].err<<','<<fits[c].width<<','<<fits[c].status<<'\n';
 }
 Fit g18=fit_full(pooled[0],"full_hod2018"),g19=fit_full(pooled[1],"full_hod2019");
 std::ofstream summary(outdir/"summary_channel_health.txt",std::ios::trunc);summary<<std::fixed<<std::setprecision(6);
 summary<<"Channel health refit (2026-10-16)\n=================================\n\nInput files and entries (v10 excluded):\n";
 for(int r=0;r<NR;++r)summary<<"  "<<fs::absolute(FILES[r])<<"  "<<entries[r]<<"\n";
 summary<<"\nNoise: pedestal median per event/channel over 200-3400 ns; sigma=max(1.4826*MAD, sample RMS, 1 ADC).\n"
  <<"DCR pulse: local excursion runs of at least 2 samples. Rates use N_events*3.2e-6 s.\n"
  <<"For each run/channel the CSV reports negative and upward-polarity counts/rates for 4*sigma and fixed 5, 10, 20 ADC cuts.\n"
  <<"Upward-excursion control floor (all channels pooled):\n";
 for(int j=0;j<4;++j){long long up=0,ev=0;for(int r=0;r<NR;++r)for(int c=0;c<NC;++c){up+=dcr[r][c].up[j];ev+=dcr[r][c].entries;}
  summary<<"  cut "<<(j==0?"4*sigma":std::to_string(j==1?5:j==2?10:20)+" ADC")<<": "<<up<<" pulses, "<<up/(ev*EXPOSURE_S)<<" Hz/channel (aggregate equal-channel mean).\n";}
 summary<<"\nMuon response v8 uses one channel/one physical axis strip per count: N_in=max(pedestal-waveform)>plane threshold in [4000,4400] ns; N_acc same per-channel cut in [2000,2400] ns.\n"
  <<"p_value_poisson_upper=P(Poisson(lambda=N_acc)>=N_in); it is an excess-tail test, not the former 3x rule.\n"
  <<"N_exp_k = MC fraction per trigger crossing station strip k (also within its orthogonal active dimension) * v8 data_dig entries. MC uses 5,000,000 cos^2(theta) rays, geometry from geom_acceptance.py, seed 42. eps_k=(N_in-N_acc)/N_exp_k; no clipping is applied.\n"
  <<"N_in and N_acc are per axis/channel, not a coincidence of both axes.\n\n"
  <<"Full-spectrum Landau fits use all event peak heights (including below threshold), model gaus(0)+landau(3), with fit range from zero to observed maximum. Reported MPV is the ROOT Landau mode (location-0.22278298*width); its error propagates the fit covariance.\n"
  <<"Pooled full-spectrum grouped fits (all 32 channels per hodoscope, all v8 entries):\n"
  <<"  hod2018 MPV="<<g18.mpv<<" +/- "<<g18.err<<" ADC, Landau width="<<g18.width<<", fit_status="<<g18.status<<"\n"
  <<"  hod2019 MPV="<<g19.mpv<<" +/- "<<g19.err<<" ADC, Landau width="<<g19.width<<", fit_status="<<g19.status;
 if(g19.status!=0)summary<<" (fit not cleanly converged; nominal value/error are not reliable)";
 summary<<"\n\n"
  <<"Per-channel results and MPVs are in channel_response_v8.csv. DCR data are in dcr_thresholds_v1_v9.csv; strip acceptance in cos2_strip_acceptance.csv.\n"
  <<"Poisson consistency criterion for the requested list: p>=0.05 (one-sided upper-tail test). Consistent channels:\n";
 int consistent=0;for(int c=0;c<NC;++c)if(poisson_excess_p(chdata[c].nin,chdata[c].nacc)>=.05){++consistent;auto[p,k]=plane_k(c);summary<<"  Ch"<<c<<" "<<p<<" k="<<k<<" Nin="<<chdata[c].nin<<" Nacc="<<chdata[c].nacc<<" p="<<poisson_excess_p(chdata[c].nin,chdata[c].nacc)<<"\n";}
 summary<<"Count consistent="<<consistent<<"/64.\nNo causal or physics interpretation is made.\n";
 std::cout<<"Wrote refit outputs to "<<fs::absolute(outdir)<<"\n";
}
