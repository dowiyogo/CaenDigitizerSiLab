#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace fs=std::filesystem;
constexpr int NC=64,NR=9,NB=800;
constexpr double BIN=0.5, XMAX=NB*BIN;
constexpr std::array<int,16> X19={53,61,54,62,55,63,52,60,51,59,50,58,49,57,48,56};
constexpr std::array<int,16> Y19={37,45,38,46,39,47,36,44,35,43,34,42,33,41,32,40};
const std::array<std::string,NR> RUNS={"v1","v2","v3","v4","v5","v6","v7","v8","v9"};
const std::array<std::string,NR> FILES={
 "lectura_vx1740_ext_nim_1M.root","lectura_vx1740_ext_nim_1M_v2.root","lectura_vx1740_ext_nim_1M_v3.root",
 "lectura_vx1740_ext_nim_1M_v4.root","lectura_vx1740_ext_nim_1M_v5.root","lectura_vx1740_ext_nim_1M_v6.root",
 "lectura_vx1740_ext_nim_1M_v7.root","lectura_vx1740_ext_nim_1M_v8.root","lectura_vx1740_ext_nim_1M_v9.root"};
struct Win{int lo=0,hi=-1;};
struct PeakFit{double a=NAN,err=NAN,spacing12=NAN,spacing23=NAN,spacing34=NAN;int components=0,status=-1;bool resolved=false;};
struct ChannelHist{std::unique_ptr<TH1D> h;long long pulses=0;};
static Win get_window(const std::vector<float>&t,double lo,double hi){
 auto a=std::lower_bound(t.begin(),t.end(),float(lo)),b=std::upper_bound(t.begin(),t.end(),float(hi));
 if(a==t.end()||b==t.begin()||a>=b)return{};return{int(a-t.begin()),int(b-t.begin())-1};
}
static double median(std::vector<double>&v){
 size_t n=v.size();auto m=v.begin()+n/2;std::nth_element(v.begin(),m,v.end());double x=*m;
 if(!(n&1U))x=(x+*std::max_element(v.begin(),m))*.5;return x;
}
struct PlaneCh {std::string plane,axis;int k;};
static PlaneCh mapping(int ch){
 if(ch<16)return{"hod2018","X",ch};if(ch<32)return{"hod2018","Y",ch-16};
 for(int k=0;k<16;++k){if(X19[k]==ch)return{"hod2019","X",k};if(Y19[k]==ch)return{"hod2019","Y",k};}
 return{"unknown","?",-1};
}

static std::vector<std::pair<double,double>> peak_candidates(const TH1D&h){
 std::vector<double>s(h.GetNbinsX()+1,0);
 for(int b=3;b<h.GetNbinsX()-2;++b)s[b]=(h.GetBinContent(b-2)+2*h.GetBinContent(b-1)+3*h.GetBinContent(b)+2*h.GetBinContent(b+1)+h.GetBinContent(b+2))/9.;
 double maximum=*std::max_element(s.begin(),s.end());std::vector<std::pair<double,double>> out;
 if(maximum<8)return out;
 const double threshold=std::max(8.,maximum*.02);
 for(int b=4;b<h.GetNbinsX()-2;++b){
  if(s[b]>=threshold&&s[b]>s[b-1]&&s[b]>=s[b+1]){
   double left=std::min(s[b-1],s[b-2]),right=std::min(s[b+1],s[b+2]);
   if(s[b]-std::max(left,right)<std::max(5.,maximum*.015))continue;
   double x=h.GetXaxis()->GetBinCenter(b);
   if(!out.empty()&&x-out.back().first<1.5){if(s[b]>out.back().second)out.back()={x,s[b]};}
   else out.emplace_back(x,s[b]);
  }
 }
 return out;
}

static bool choose_spacing(const std::vector<std::pair<double,double>>&p,double&ped,double&spacing,int&components){
 if(p.size()<2)return false;
 double best=-1;bool found=false;int bestn=0;
 for(size_t i=0;i+1<p.size();++i){
  double a=p[i+1].first-p[i].first;if(a<2.0||a>100.)continue;
  int n=2,last=int(i+1);double score=std::min(p[i].second,p[i+1].second);
  for(int step=2;step<=3;++step){
   size_t chosen=p.size();double errbest=1e9;
   for(size_t j=last+1;j<p.size();++j){double e=std::abs(p[j].first-(p[i].first+step*a))/a;
    if(e<errbest){errbest=e;chosen=j;}}
   if(chosen<p.size()&&errbest<.22){++n;last=chosen;score+=p[chosen].second/(step+1.);}
   else break;
  }
  score+=n*std::log1p(p[i].second+p[i+1].second);
  if(!found||n>bestn||(n==bestn&&score>best)){best=score;bestn=n;ped=p[i].first;spacing=a;components=n;found=true;}
 }
 return found;
}

static PeakFit fit_spectrum(TH1D&h,const std::string&name){
 PeakFit r;
 if(h.GetEntries()<80)return r;
 auto peaks=peak_candidates(h);double ped=0,spacing=0;int ncomp=0;
 if(!choose_spacing(peaks,ped,spacing,ncomp))return r;
 ncomp=std::clamp(ncomp,2,4);r.components=ncomp;
 const int npar=2+2*ncomp;
 auto model=[ncomp](double*x,double*p){double y=0,xx=x[0];
  for(int n=0;n<ncomp;++n){double mean=p[0]+n*p[1],amp=p[2+2*n],sig=p[3+2*n];if(sig<=0)continue;
   double z=(xx-mean)/sig;y+=amp*std::exp(-.5*z*z);}
  return y;
 };
 TF1 fit(name.c_str(),model,0,XMAX,npar);
 fit.SetParName(0,"pedestal_mean");fit.SetParName(1,"spacing_1pe");
 fit.SetParameter(0,ped);fit.SetParameter(1,spacing);
 fit.SetParLimits(0,std::max(-2.,ped-3.),ped+3.);
 fit.SetParLimits(1,std::max(1.,spacing*.55),spacing*1.55);
 for(int n=0;n<ncomp;++n){
  int hi=2+2*n,wi=3+2*n;int b=h.GetXaxis()->FindBin(ped+n*spacing);
  double height=std::max(1.,h.GetBinContent(b));
  fit.SetParameter(hi,height);fit.SetParameter(wi,std::max(.5,std::min(8.,spacing*.22)));
  fit.SetParLimits(hi,0.,std::max(10.,h.GetMaximum()*2.));
  fit.SetParLimits(wi,.25,std::max(2.,spacing*.7));
 }
 TFitResultPtr fr=h.Fit(&fit,"LQ0S","",0,XMAX);r.status=int(fr);
 if(!fr.Get()||!fr->IsValid()||!std::isfinite(fit.GetParameter(1)))return r;
 r.a=fit.GetParameter(1);auto cov=fr->GetCovarianceMatrix();r.err=std::sqrt(std::max(0.,cov(1,1)));
 r.spacing12=r.a;r.spacing23=ncomp>=3?r.a:NAN;r.spacing34=ncomp>=4?r.a:NAN;
 // Require distinct resolved adjacent peaks and a successful fit, rather than reporting a forced seed.
 double sig0=fit.GetParameter(3),sig1=ncomp>1?fit.GetParameter(5):sig0;
 double low=std::max(1.,spacing*.55),high=spacing*1.55;
 bool away_from_bound=r.a>low+0.025*(high-low)&&r.a<high-0.025*(high-low);
 double one_pe_height=fit.GetParameter(4);
 r.resolved=r.status==0&&away_from_bound&&r.a>2.0&&r.a>1.5*std::max(sig0,sig1)&&one_pe_height>=std::max(5.,.02*h.GetMaximum());
 return r;
}

static double linear_zero(const std::array<double,NR>&y,const std::array<bool,NR>&ok){
 double sx=0,sy=0,sxx=0,sxy=0,n=0;
 for(int i=0;i<NR;++i)if(ok[i]&&std::isfinite(y[i])){double x=i+1;sx+=x;sy+=y[i];sxx+=x*x;sxy+=x*y[i];++n;}
 double d=n*sxx-sx*sx;if(n<4||d==0)return NAN;double b=(n*sxy-sx*sy)/d,a=(sy-b*sx)/n;
 double mean=sy/n,sst=0,sse=0;for(int i=0;i<NR;++i)if(ok[i]){double x=i+1;sst+=(y[i]-mean)*(y[i]-mean);sse+=(y[i]-(a+b*x))*(y[i]-(a+b*x));}
 double r2=sst>0?1.-sse/sst:0.;return r2>=.8&&std::abs(b)>1e-12?-a/b:NAN;
}

static void make_plots(const fs::path&out,const std::array<std::array<ChannelHist,NC>,NR>&ch,
 const std::array<std::array<PeakFit,NC>,NR>&fits,const std::array<std::array<double,4>,NR>&plane_med,
 const std::array<std::array<int,4>,NR>&plane_n){
 gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);gStyle->SetOptFit(0);
 TCanvas spectra("spe_examples","",1000,750);spectra.Divide(2,2);
 const std::array<int,4> rep={0,16,53,37};
 for(int p=0;p<4;++p){int best=-1;long long bestn=-1;
  for(int c=0;c<NC;++c){auto m=mapping(c);int gp=(m.plane=="hod2018"?0:2)+(m.axis=="Y"?1:0);if(gp!=p)continue;
   if(fits[7][c].resolved&&ch[7][c].pulses>bestn){best=c;bestn=ch[7][c].pulses;}}
  if(best<0)best=rep[p];spectra.cd(p+1);gPad->SetLogy();
  ch[7][best].h->SetTitle((mapping(best).plane+" "+mapping(best).axis+" v8 Ch"+std::to_string(best)+";Pulse amplitude [ADC];Pulses").c_str());
  ch[7][best].h->GetXaxis()->SetRangeUser(0,100);ch[7][best].h->SetLineColor(kBlue+1);ch[7][best].h->Draw("HIST");
 }
 spectra.SaveAs((out/"spe_spectra_examples.png").string().c_str());
 TCanvas c("a1pe_run","",1000,750);c.SetLeftMargin(.11);c.SetBottomMargin(.12);c.SetGrid();
 double ymax=0;for(int r=0;r<NR;++r)for(int p=0;p<4;++p)if(plane_n[r][p]&&std::isfinite(plane_med[r][p]))ymax=std::max(ymax,plane_med[r][p]);
 TH1D frame("spe_frame","Single-photoelectron peak spacing;Run index;A_{1pe} [ADC]",9,.5,9.5);frame.SetMinimum(0);frame.SetMaximum(std::max(10.,ymax*1.2));frame.SetStats(0);
 for(int i=0;i<NR;++i)frame.GetXaxis()->SetBinLabel(i+1, std::to_string(i+1).c_str());frame.GetXaxis()->SetTitle("Run (v1-v9)");frame.Draw();
 constexpr int color[]={kBlue+1,kRed+1};constexpr int marker[]={20,21};
 std::array<TGraphErrors,4> graphs;
 for(int p=0;p<4;++p){int h=p/2,ax=p%2,n=0;for(int r=0;r<NR;++r)if(plane_n[r][p]){graphs[p].SetPoint(n,r+1,plane_med[r][p]);graphs[p].SetPointError(n++,0,0);}
  graphs[p].SetMarkerColor(color[h]);graphs[p].SetLineColor(color[h]);graphs[p].SetMarkerStyle(marker[ax]);graphs[p].SetMarkerSize(1.2);graphs[p].Draw("PL SAME");}
 std::array<double,NR> y18{},y19{};std::array<bool,NR> ok18{},ok19{};
 for(int r=0;r<NR;++r){int a=0,b=1;int n=plane_n[r][a]+plane_n[r][b];if(n){y18[r]=plane_n[r][a]&&plane_n[r][b]?(plane_med[r][a]+plane_med[r][b])*.5:(plane_n[r][a]?plane_med[r][a]:plane_med[r][b]);ok18[r]=true;}
  a=2;b=3;n=plane_n[r][a]+plane_n[r][b];if(n){y19[r]=plane_n[r][a]&&plane_n[r][b]?(plane_med[r][a]+plane_med[r][b])*.5:(plane_n[r][a]?plane_med[r][a]:plane_med[r][b]);ok19[r]=true;}}
 for(int h=0;h<2;++h){const auto& yy=h?y19:y18;const auto& ok=h?ok19:ok18;double sx=0,sy=0,sxx=0,sxy=0,n=0;
  for(int r=0;r<NR;++r)if(ok[r]){double x=r+1;sx+=x;sy+=yy[r];sxx+=x*x;sxy+=x*yy[r];++n;}
  double d=n*sxx-sx*sx;if(n>=4&&d){double b=(n*sxy-sx*sy)/d,a=(sy-b*sx)/n,mean=sy/n,sst=0,sse=0;for(int r=0;r<NR;++r)if(ok[r]){double x=r+1;sst+=(yy[r]-mean)*(yy[r]-mean);sse+=(yy[r]-(a+b*x))*(yy[r]-(a+b*x));}double r2=sst>0?1.-sse/sst:0.;if(r2>=.8){TGraph line(81);for(int q=0;q<=80;++q){double x=1.+q*.1;line.SetPoint(q,x,a+b*x);}line.SetLineColor(color[h]);line.SetLineStyle(2);line.Draw("L SAME");double zero=-a/b;if(zero>0&&zero<9.5){TLine vline(zero,0,zero,frame.GetMaximum());vline.SetLineColor(color[h]);vline.SetLineStyle(3);vline.Draw();}}}
 }
 TLegend leg(.64,.68,.94,.93);leg.SetBorderSize(0);leg.SetFillStyle(0);
 leg.AddEntry(&graphs[0],"hod2018 X","lp");leg.AddEntry(&graphs[1],"hod2018 Y","lp");leg.AddEntry(&graphs[2],"hod2019 X","lp");leg.AddEntry(&graphs[3],"hod2019 Y","lp");leg.Draw();
 c.SaveAs((out/"a1pe_vs_run.png").string().c_str());
}

int main(int argc,char**argv){
 fs::path out=argc>1?argv[1]:"outputs/spe_gain_20261016";fs::create_directories(out);gROOT->SetBatch(kTRUE);
 std::array<std::array<ChannelHist,NC>,NR> ch{};std::array<std::array<std::unique_ptr<TH1D>,4>,NR> planes{};
 for(int r=0;r<NR;++r){
  for(int c=0;c<NC;++c){std::string n="spe_"+RUNS[r]+"_ch"+std::to_string(c);ch[r][c].h=std::make_unique<TH1D>(n.c_str(),"",NB,0,XMAX);ch[r][c].h->SetDirectory(nullptr);}
  for(int p=0;p<4;++p){std::string n="spe_plane_"+RUNS[r]+"_"+std::to_string(p);planes[r][p]=std::make_unique<TH1D>(n.c_str(),"",NB,0,XMAX);planes[r][p]->SetDirectory(nullptr);}
 }
 std::array<long long,NR> nentries{};std::vector<double> base,dev;
 for(int r=0;r<NR;++r){
  std::unique_ptr<TFile> file(TFile::Open(FILES[r].c_str(),"READ"));if(!file||file->IsZombie()){std::cerr<<"Cannot open "<<FILES[r]<<"\n";return 2;}
  auto*t=dynamic_cast<TTree*>(file->Get("data_dig"));if(!t){std::cerr<<"Missing data_dig in "<<FILES[r]<<"\n";return 2;}
  nentries[r]=t->GetEntries();t->SetBranchStatus("*",0);t->SetCacheSize(128LL*1024*1024);
  std::array<std::vector<float>*,NC>wf{};std::vector<float>*time=nullptr;
  for(int c=0;c<NC;++c){std::string b="Ch"+std::to_string(c);t->SetBranchStatus(b.c_str(),1);if(t->SetBranchAddress(b.c_str(),&wf[c])<0)return 2;}
  t->SetBranchStatus("time",1);if(t->SetBranchAddress("time",&time)<0)return 2;Win pre;
  for(long long e=0;e<nentries[r];++e){t->GetEntry(e);if(!time||time->empty())continue;if(e==0)pre=get_window(*time,200,3400);
   for(int c=0;c<NC;++c){if(!wf[c]||wf[c]->empty())continue;auto&v=*wf[c];int lo=pre.lo,hi=std::min(pre.hi,int(v.size())-1);if(hi<lo)continue;
    base.clear();base.reserve(hi-lo+1);for(int s=lo;s<=hi;++s)base.push_back(v[s]);double ped=median(base),mean=std::accumulate(base.begin(),base.end(),0.)/base.size();
    dev.clear();double ss=0;for(double q:base){dev.push_back(std::abs(q-ped));ss+=(q-mean)*(q-mean);}double sigma=std::max({1.4826*median(dev),std::sqrt(ss/std::max<size_t>(1,base.size()-1)),1.0});double threshold=2.5*sigma;
    int s=lo;while(s<=hi){double amp=ped-v[s];if(amp<=threshold){++s;continue;}int start=s;double peak=amp;
     while(s<=hi&&ped-v[s]>threshold){peak=std::max(peak,double(ped-v[s]));++s;}
     if(s-start>=2){ch[r][c].h->Fill(peak);++ch[r][c].pulses;auto m=mapping(c);int p=(m.plane=="hod2018"?0:2)+(m.axis=="Y"?1:0);planes[r][p]->Fill(peak);}
    }
   }
   if(e&&e%20000==0)std::cout<<RUNS[r]<<" "<<e<<"/"<<nentries[r]<<"\r"<<std::flush;
  }
  std::cout<<RUNS[r]<<" entries="<<nentries[r]<<"\n";
 }
 std::array<std::array<PeakFit,NC>,NR> fits{};std::array<std::array<double,4>,NR> plane_med{};std::array<std::array<int,4>,NR> plane_n{};
 std::ofstream csv(out/"spe_gain.csv");csv<<"record_type,run,ch,hodoscope,axis,k,pulses,components,a1pe_adc,a1pe_error_adc,spacing_1_2_adc,spacing_2_3_adc,spacing_3_4_adc,fit_status,resolved\n"<<std::setprecision(10);
 for(int r=0;r<NR;++r){std::array<std::vector<double>,4> vals;
  for(int c=0;c<NC;++c){auto m=mapping(c);int p=(m.plane=="hod2018"?0:2)+(m.axis=="Y"?1:0);fits[r][c]=fit_spectrum(*ch[r][c].h,"fit_"+RUNS[r]+"_"+std::to_string(c));auto&f=fits[r][c];
   if(f.resolved)vals[p].push_back(f.a);
   csv<<"channel,"<<RUNS[r]<<','<<c<<','<<m.plane<<','<<m.axis<<','<<m.k<<','<<ch[r][c].pulses<<','<<f.components<<',';
   if(f.resolved)csv<<f.a<<','<<f.err<<','<<f.spacing12<<','<<f.spacing23<<','<<f.spacing34;else csv<<",,,,";
   csv<<','<<f.status<<','<<(f.resolved?"y":"n")<<'\n';
  }
  for(int p=0;p<4;++p){if(vals[p].empty()){plane_med[r][p]=NAN;plane_n[r][p]=0;}else{std::sort(vals[p].begin(),vals[p].end());plane_n[r][p]=vals[p].size();size_t n=vals[p].size();plane_med[r][p]=n%2?vals[p][n/2]:(vals[p][n/2-1]+vals[p][n/2])*.5;}
   csv<<"plane_median,"<<RUNS[r]<<",,"<<(p<2?"hod2018":"hod2019")<<','<<(p%2?"Y":"X")<<",,"<<planes[r][p]->GetEntries()<<','<<plane_n[r][p]<<',';
   if(plane_n[r][p])csv<<plane_med[r][p];else csv<<"nan";csv<<",,,,,,"<<(plane_n[r][p]?"y":"n")<<'\n';
  }
 }
 // Ratios use each station's count-weighted median across its resolved X/Y strip fits.
 std::array<std::array<double,NR>,2> station{};std::array<std::array<bool,NR>,2> valid{};
 for(int r=0;r<NR;++r)for(int h=0;h<2;++h){int a=h?2:0,b=a+1;std::vector<double>x;if(plane_n[r][a])x.push_back(plane_med[r][a]);if(plane_n[r][b])x.push_back(plane_med[r][b]);if(!x.empty()){station[h][r]=x.size()==2?(x[0]+x[1])*.5:x[0];valid[h][r]=true;}}
 for(int r=0;r<NR;++r){double ratio=valid[0][r]&&valid[1][r]?station[1][r]/station[0][r]:NAN;
  csv<<"gain_ratio,"<<RUNS[r]<<",,hod2019_vs_hod2018,all,,,,";
  if(std::isfinite(ratio))csv<<ratio;else csv<<"nan";csv<<",,,,,0,y\n";}
 std::array<double,NR> ys18{},ys19{};std::array<bool,NR> ok18{},ok19{};for(int r=0;r<NR;++r){ys18[r]=station[0][r];ys19[r]=station[1][r];ok18[r]=valid[0][r];ok19[r]=valid[1][r];}
 double zero18=linear_zero(ys18,ok18),zero19=linear_zero(ys19,ok19);
 std::ofstream summary(out/"summary_spe_gain.txt");summary<<std::fixed<<std::setprecision(5)<<"SPE gain comparison; pulse threshold 2.5*sigma; pulse width >=2 samples. v10 excluded.\n"
  <<"Pedestal median and sigma=max(1.4826*MAD, sample RMS, 1 ADC) are calculated per event/channel over 200-3400 ns. Amplitudes are pedestal-minimum, binned at 0.5 ADC.\n"
  <<"A1pe is the fitted common spacing between Gaussian pedestal-tail and successive PE components; rows without supported peaks are unresolved, not extrapolated.\n"
  <<"Run-index linear zero crossings (run number treated as x=1..9; station median averaged over available X/Y plane medians; shown only for R^2>=0.8): ";
 if(std::isfinite(zero18))summary<<"hod2018="<<zero18<<"; ";else summary<<"hod2018=no stable linear estimate; ";
 if(std::isfinite(zero19))summary<<"hod2019="<<zero19<<".\n";else summary<<"hod2019=no stable linear estimate.\n";
 summary<<"No cause is inferred. Channel-level and plane-median values are in spe_gain.csv.\n";
 make_plots(out,ch,fits,plane_med,plane_n);
 std::cout<<"Output: "<<fs::absolute(out)<<"\n";
}
