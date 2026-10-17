// poisson_analysis.cxx
// Test estadístico de distribución de Poisson en tiempos de llegada de triggers.
// Lee SOLO la rama timestamp (Long_t, µs) — sin formas de onda, ejecución < 30 s.
//
// Física: muones cósmicos ~ proceso de Poisson estacionario con tasa λ.
//   Tests: (1) KS exponencial en ΔT, (2) χ² de Poisson en conteos/intervalo,
//          (3) índice de dispersión D=σ²/µ, (4) estabilidad de tasa.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o poisson_analysis poisson_analysis.cxx
// Usage: ./poisson_analysis [--v7 archivo] [--v8 archivo]

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// Flujo muones cósmicos a nivel del mar: ~1 muón/cm²/min.
// Área trigger: 65×65 mm = 42.25 cm² → λ_esperada ≈ 42.25/60 ≈ 0.704 Hz.
constexpr double kLambdaExp_Hz = 42.25 / 60.0;

//==============================================================================
// Datos por run
//==============================================================================
struct RunData {
    std::string label, filename;
    long long n = 0;
    std::vector<double> ts;           // timestamps [µs], relativos al primero
    std::vector<double> dt_all;       // todos los ΔT [µs]
    std::vector<double> dt_clean;     // ΔT sin gaps de DAQ
    std::vector<std::pair<int,double>> gaps; // (índice, duración µs)
    double T_total = 0;    // µs
    double lam_hat = 0;    // µs⁻¹
    double lam_hz  = 0;    // Hz
    double t_mean  = 0;    // µs = 1/lam_hat
    // Resultados del análisis
    double lam_fit = 0, lam_fit_err = 0;  // µs⁻¹
    double D_KS = 0, p_KS = 0;
    double tau_dead = 0;   // µs
    double D1ms = 0;       // índice de dispersión a Δτ=1ms
    double pchi2_1ms = 0;  // p-value chi² a Δτ=1ms
};

//==============================================================================
// Carga timestamps
//==============================================================================
static RunData load_run(const std::string& fname, const std::string& lbl) {
    RunData rd; rd.label = lbl; rd.filename = fname;
    TFile* f = TFile::Open(fname.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "ERROR: no se puede abrir " << fname << '\n'; return rd; }
    TTree* t = static_cast<TTree*>(f->Get("data_dig"));
    if (!t) { std::cerr << "ERROR: no data_dig en " << fname << '\n'; f->Close(); return rd; }

    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("timestamp", 1);
    Long_t ts_val = 0;
    t->SetBranchAddress("timestamp", &ts_val);
    const long long nEnt = t->GetEntries();
    rd.ts.resize(nEnt);
    for (long long i = 0; i < nEnt; i++) { t->GetEntry(i); rd.ts[i] = static_cast<double>(ts_val); }
    f->Close();

    // Ordena y calcula tiempos relativos
    std::sort(rd.ts.begin(), rd.ts.end());
    const double t0 = rd.ts[0];
    for (auto& v : rd.ts) v -= t0;

    rd.n       = nEnt;
    rd.T_total = rd.ts.back();               // µs
    rd.lam_hat = (nEnt - 1.0) / rd.T_total;  // µs⁻¹
    rd.lam_hz  = rd.lam_hat * 1e6;           // Hz
    rd.t_mean  = 1.0 / rd.lam_hat;           // µs

    // ΔT entre eventos consecutivos
    rd.dt_all.resize(nEnt - 1);
    for (long long i = 0; i < nEnt - 1; i++) rd.dt_all[i] = rd.ts[i+1] - rd.ts[i];

    // Gap de DAQ: ΔT > 50 × T_mean
    const double gap_thr = 50.0 * rd.t_mean;
    for (size_t i = 0; i < rd.dt_all.size(); i++) {
        if (rd.dt_all[i] > gap_thr) rd.gaps.push_back({(int)i, rd.dt_all[i]});
        else                         rd.dt_clean.push_back(rd.dt_all[i]);
    }
    // Dead time: ΔT mínimo observado
    if (!rd.dt_clean.empty())
        rd.tau_dead = *std::min_element(rd.dt_clean.begin(), rd.dt_clean.end());

    return rd;
}

//==============================================================================
// Ajuste exponencial + test KS
//==============================================================================
static void fit_exp_ks(RunData& rd) {
    if (rd.dt_clean.size() < 10) return;
    const int    nbins  = 200;
    const double range  = 10.0 * rd.t_mean;
    const double bin_w  = range / nbins;
    const double fit_lo = std::max(rd.tau_dead * 1.05, rd.t_mean * 0.02);
    const double fit_hi = 5.0 * rd.t_mean;

    TH1D* h = new TH1D("_hfit_", "", nbins, 0, range);
    for (double v : rd.dt_clean) if (v < range) h->Fill(v);

    TF1 fE("_fE_", "[0]*exp(-[1]*x)", fit_lo, fit_hi);
    fE.SetParameters(rd.dt_clean.size() * bin_w * rd.lam_hat, rd.lam_hat);
    fE.SetParLimits(1, rd.lam_hat * 0.05, rd.lam_hat * 20.0);
    if (h->Fit(&fE, "RQN") == 0) {
        rd.lam_fit     = fE.GetParameter(1);
        rd.lam_fit_err = fE.GetParError(1);
    } else {
        rd.lam_fit = rd.lam_hat; rd.lam_fit_err = rd.lam_hat / sqrt(rd.dt_clean.size());
    }
    delete h;

    // KS test vs Exp(lam_fit)
    std::vector<double> dts = rd.dt_clean;
    std::sort(dts.begin(), dts.end());
    const int n = (int)dts.size();
    rd.D_KS = 0;
    for (int i = 0; i < n; i++) {
        double F_theo  = 1.0 - std::exp(-rd.lam_fit * dts[i]);
        double d_plus  = std::abs((double)(i+1)/n - F_theo);
        double d_minus = std::abs((double)(i)/n   - F_theo);
        rd.D_KS = std::max(rd.D_KS, std::max(d_plus, d_minus));
    }
    rd.p_KS = TMath::KolmogorovProb(rd.D_KS * std::sqrt(n));
}

//==============================================================================
// Índice de dispersión D = σ²/µ para un Δτ dado [µs]
//==============================================================================
static std::pair<double,double> disp_idx(const std::vector<double>& ts,
                                          double T_total, double dtau) {
    const int k = (int)(T_total / dtau);
    if (k < 3) return {1.0, 1.0};
    std::vector<int> cnt(k, 0);
    for (double t : ts) { int b = (int)(t / dtau); if (b>=0&&b<k) cnt[b]++; }
    double mu=0, mu2=0;
    for (int c : cnt) { mu += c; mu2 += (double)c*c; }
    mu /= k; mu2 /= k;
    double D = mu > 0 ? (mu2 - mu*mu) / mu : 1.0;
    return {D, 2.0/std::sqrt(k)};   // D, 2σ error
}

//==============================================================================
// Test χ² vs Poisson(µ) en conteos por intervalo
//==============================================================================
static std::pair<double,double> chi2_poisson(const std::vector<double>& ts,
                                              double T_total, double dtau) {
    const int k = (int)(T_total / dtau);
    if (k < 5) return {0.0, 1.0};
    std::vector<int> cnt(k, 0);
    for (double t : ts) { int b = (int)(t/dtau); if(b>=0&&b<k) cnt[b]++; }
    double mu = 0; for (int c : cnt) mu += c; mu /= k;
    // Distribución de conteos
    std::map<int,int> freq; for (int c : cnt) freq[c]++;
    int n_max = cnt.empty() ? 0 : *std::max_element(cnt.begin(), cnt.end());

    // χ² agrupando bins con esperado < 5
    double chi2 = 0; int ndof = 0;
    int g_obs = 0; double g_exp = 0;
    for (int n = 0; n <= n_max + 5; n++) {
        double logP = (mu > 0 && n >= 0) ? -mu + n*std::log(mu) - TMath::LnGamma(n+1) : -1e30;
        double expN = std::exp(logP) * k;
        int obsN = freq.count(n) ? freq.at(n) : 0;
        g_obs += obsN; g_exp += expN;
        if (g_exp >= 5.0 || n == n_max + 5) {
            if (g_exp > 0) { chi2 += (g_obs-g_exp)*(g_obs-g_exp)/g_exp; ++ndof; }
            g_obs = 0; g_exp = 0;
        }
    }
    ndof = std::max(ndof - 2, 1);  // -1 param (µ estimado), mínimo 1
    return {chi2, TMath::Prob(chi2, ndof)};
}

//==============================================================================
// Figuras para un run
//==============================================================================
static void make_figs(RunData& rd, const fs::path& out) {
    const int colv = (rd.label=="v7") ? kBlue+1 : kOrange+7;
    const std::string L = rd.label;

    //── 1A: Histograma ΔT + fit exponencial ─────────────────────────────────
    {
        const int nb = 200; const double rng = 10.0*rd.t_mean;
        TH1D h(("hdt_"+L).c_str(),
            ("Inter-arrival time "+L+";#DeltaT [#mus];Entries").c_str(), nb, 0, rng);
        for (double v : rd.dt_clean) if (v<rng) h.Fill(v);
        h.SetLineColor(colv); h.SetLineWidth(2); h.SetStats(0);

        const double bin_w = rng/nb;
        const double fit_lo = std::max(rd.tau_dead*1.05, rd.t_mean*0.02);
        const double fit_hi = 5.0*rd.t_mean;
        TF1 fE(("fE_"+L).c_str(),"[0]*exp(-[1]*x)", fit_lo, fit_hi);
        fE.SetParameters(rd.dt_clean.size()*bin_w*rd.lam_fit, rd.lam_fit);
        fE.SetLineColor(kRed); fE.SetLineWidth(2);

        TCanvas c(("c1a_"+L).c_str(),"",1000,680);
        c.SetLogy(); c.SetLeftMargin(0.12); c.SetBottomMargin(0.13);
        h.Draw("HIST"); fE.Draw("SAME");
        TLatex lt; lt.SetNDC(); lt.SetTextSize(0.033);
        char buf[160];
        std::snprintf(buf,sizeof(buf),"#lambda_{hat}=%.3f Hz   #lambda_{fit}=%.3f#pm%.3f Hz",
                      rd.lam_hz, rd.lam_fit*1e6, rd.lam_fit_err*1e6);
        lt.DrawLatex(0.35,0.88,buf);
        std::snprintf(buf,sizeof(buf),"#tau_{dead}#approx%.0f #mus  (%.2f%%)",
                      rd.tau_dead, rd.tau_dead*rd.lam_hat*100);
        lt.DrawLatex(0.35,0.84,buf);
        std::snprintf(buf,sizeof(buf),"D_{KS}=%.5f  p_{KS}=%.4f",rd.D_KS,rd.p_KS);
        lt.DrawLatex(0.35,0.80,buf);
        if(!rd.gaps.empty()){
            std::snprintf(buf,sizeof(buf),"N_{gaps DAQ}=%d",(int)rd.gaps.size());
            lt.SetTextColor(kRed); lt.DrawLatex(0.35,0.76,buf); lt.SetTextColor(kBlack);}
        c.SaveAs((out/("interarrival_"+L+".png")).string().c_str());
    }

    //── 1B: Q-Q plot exponencial ─────────────────────────────────────────────
    {
        std::vector<double> dts = rd.dt_clean;
        std::sort(dts.begin(),dts.end());
        const int n = (int)dts.size();
        const int n_plot = std::min(n, 2000);
        std::vector<double> qemp, qtheo;
        for (int i = 0; i < n_plot; i++) {
            int idx = (int)((double)i/n_plot*n);
            if (idx >= n) idx = n-1;
            double F = (double)(idx+1)/n;
            F = std::min(F, 1.0 - 1e-9);
            qemp.push_back(dts[idx]);
            qtheo.push_back(-std::log(1.0-F)/rd.lam_fit);
        }
        TGraph g((int)qemp.size(), qtheo.data(), qemp.data());
        g.SetTitle(("Q-Q Exponencial "+L+";Cuantiles teóricos [#mus];Cuantiles empíricos [#mus]").c_str());
        g.SetMarkerStyle(20); g.SetMarkerSize(0.4); g.SetMarkerColor(colv);
        double xmax = *std::max_element(qtheo.begin(),qtheo.end());
        TF1 fd("_fd","x",0,xmax*1.1);fd.SetLineColor(kRed);fd.SetLineWidth(2);
        TCanvas c(("c_qq_"+L).c_str(),"",900,900);
        c.SetLeftMargin(0.13); c.SetBottomMargin(0.13);
        g.Draw("AP"); fd.Draw("SAME");
        c.SaveAs((out/("qq_interarrival_"+L+".png")).string().c_str());
    }

    //── 2A: Counts per interval (Δτ=1ms) ────────────────────────────────────
    {
        const double dtau = 1000.0;  // µs = 1 ms
        const int k = (int)(rd.T_total/dtau);
        if (k > 2) {
            std::vector<int> cnt(k,0);
            for (double t : rd.ts) { int b=(int)(t/dtau); if(b>=0&&b<k) cnt[b]++; }
            double mu=0; for(int c:cnt) mu+=c; mu/=k;
            int nmax=*std::max_element(cnt.begin(),cnt.end());
            TH1D hobs(("hcnt_"+L).c_str(),
                ("Counts per 1ms "+L+";N por bin;Bins").c_str(),
                nmax+4,-0.5,nmax+3.5);
            for (int c:cnt) hobs.Fill(c);
            hobs.SetLineColor(colv); hobs.SetLineWidth(2); hobs.SetStats(0);
            TH1D hpoi(("hpoi_"+L).c_str(),"",nmax+4,-0.5,nmax+3.5);
            for (int n=0;n<=nmax+3;n++){
                double logP=mu>0?-mu+n*std::log(mu)-TMath::LnGamma(n+1):-1e30;
                hpoi.SetBinContent(n+1,std::exp(logP)*k);}
            hpoi.SetLineColor(kRed); hpoi.SetLineStyle(2); hpoi.SetLineWidth(2);
            TCanvas c(("c2a_"+L).c_str(),"",1000,680);
            c.SetLogy(); c.SetLeftMargin(0.12); c.SetBottomMargin(0.13);
            double ymax=std::max(hobs.GetMaximum(),hpoi.GetMaximum())*2;
            hobs.GetYaxis()->SetRangeUser(0.5,ymax);
            hobs.Draw("HIST"); hpoi.Draw("HIST SAME");
            TLegend lg(0.60,0.80,0.92,0.92);lg.SetBorderSize(0);lg.SetFillStyle(0);
            lg.AddEntry(&hobs,"Observado","l");lg.AddEntry(&hpoi,"Poisson(#mu)","l");lg.Draw();
            TLatex lt;lt.SetNDC();lt.SetTextSize(0.032);
            char buf[80];
            std::snprintf(buf,sizeof(buf),"#mu=%.2f  D=%.3f  p_{#chi^2}=%.3f",
                          mu,rd.D1ms,rd.pchi2_1ms);
            lt.DrawLatex(0.14,0.93,buf);
            c.SaveAs((out/("counts_per_interval_1ms_"+L+".png")).string().c_str());
        }
    }

    //── 2B: Índice de dispersión D vs Δτ ─────────────────────────────────────
    {
        const std::vector<double> dtaus={100,500,1000,5000,10000,60000};
        int np=(int)dtaus.size();
        std::vector<double> xs(np),Ds(np),eDs(np);
        for (int i=0;i<np;i++){
            auto [D,eD]=disp_idx(rd.ts,rd.T_total,dtaus[i]);
            xs[i]=dtaus[i];Ds[i]=D;eDs[i]=eD;}
        TGraphErrors g(np,xs.data(),Ds.data(),std::vector<double>(np,0).data(),eDs.data());
        g.SetTitle(("Índice dispersión "+L+";#Delta#tau [#mus];D = #sigma^{2}/<N>").c_str());
        g.SetMarkerStyle(20);g.SetMarkerColor(colv);g.SetLineColor(colv);g.SetLineWidth(2);
        TCanvas c(("c2b_"+L).c_str(),"",1000,680);
        c.SetLogx();c.SetLeftMargin(0.12);c.SetBottomMargin(0.13);
        g.Draw("APL");
        double Dmax=*std::max_element(Ds.begin(),Ds.end());
        g.GetYaxis()->SetRangeUser(0,std::max(2.0,Dmax*1.3));
        TLine l1(dtaus.front(),1,dtaus.back(),1);
        l1.SetLineColor(kRed);l1.SetLineWidth(2);l1.SetLineStyle(2);l1.Draw();
        c.SaveAs((out/("dispersion_index_vs_binwidth_"+L+".png")).string().c_str());
    }

    //── 3: Rate vs tiempo + distribución de tasa ─────────────────────────────
    {
        const double win=60000.0;  // 60 s en µs
        // Limita a ≤3000 puntos
        double slide=std::max(win*0.1, rd.T_total/3000.0);
        std::vector<double> tcen,rates;
        for (double t=0;t+win<=rd.T_total;t+=slide){
            auto lo=std::lower_bound(rd.ts.begin(),rd.ts.end(),t);
            auto hi=std::lower_bound(rd.ts.begin(),rd.ts.end(),t+win);
            int cnt=(int)(hi-lo);
            tcen.push_back((t+win/2)/60e6);  // minutos
            rates.push_back(cnt/(win*1e-6));  // Hz
        }
        if (!tcen.empty()) {
            // Distribución de tasas
            double rmin=*std::min_element(rates.begin(),rates.end());
            double rmax=*std::max_element(rates.begin(),rates.end());
            TH1D hrate(("hrate_"+L).c_str(),
                ("Distribución de tasa "+L+";Tasa [Hz];Bins ventana").c_str(),
                40,std::max(0.0,rmin*0.8),rmax*1.2);
            for (double r:rates) hrate.Fill(r);
            hrate.SetLineColor(colv);hrate.SetLineWidth(2);
            hrate.Fit("gaus","Q");
            TCanvas cr(("cr_"+L).c_str(),"",1200,680);
            int np=(int)tcen.size();
            TGraph g(np,tcen.data(),rates.data());
            g.SetTitle(("Tasa vs tiempo "+L+";Tiempo [min];Tasa [Hz]").c_str());
            g.SetLineColor(colv);g.SetLineWidth(1);
            cr.SetLeftMargin(0.10);cr.SetBottomMargin(0.13);
            g.Draw("AL");
            g.GetYaxis()->SetRangeUser(0,rd.lam_hz*3+1);
            // Banda ±2σ_Poisson
            const double sig_band=2*std::sqrt(rd.lam_hz/(win*1e-6));
            double xlo=tcen.front(),xhi=tcen.back();
            TLine lmean(xlo,rd.lam_hz,xhi,rd.lam_hz);
            lmean.SetLineColor(kRed);lmean.SetLineWidth(2);lmean.Draw();
            TLine lup(xlo,rd.lam_hz+sig_band,xhi,rd.lam_hz+sig_band);
            lup.SetLineColor(kRed);lup.SetLineStyle(2);lup.SetLineWidth(1);lup.Draw();
            TLine ldn(xlo,rd.lam_hz-sig_band,xhi,rd.lam_hz-sig_band);
            ldn.SetLineColor(kRed);ldn.SetLineStyle(2);ldn.SetLineWidth(1);ldn.Draw();
            cr.SaveAs((out/("rate_vs_time_"+L+".png")).string().c_str());
        }
    }
}

//==============================================================================
// Figura de comparación (4 paneles)
//==============================================================================
static void make_comparison(const RunData& r0, const RunData& r1, const fs::path& out) {
    const RunData* rds[2]={&r0,&r1};
    const int cols[2]={kBlue+1,kOrange+7};
    TCanvas cc("c_comp","",1600,1200); cc.Divide(2,2);

    //── Panel 1: ΔT normalizados ──────────────────────────────────────────────
    cc.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    {
        TH1D* hs[2]={};
        for (int v=0;v<2;v++){
            const RunData&r=*rds[v];
            hs[v]=new TH1D(("_n"+r.label).c_str(),
                "ΔT normalizado;#DeltaT / <#DeltaT>;Densidad",100,0,10);
            for (double dt:r.dt_clean){double x=dt/r.t_mean;if(x<10)hs[v]->Fill(x);}
            if(hs[v]->Integral()>0)hs[v]->Scale(1.0/hs[v]->Integral()/(10.0/100));
            hs[v]->SetLineColor(cols[v]);hs[v]->SetLineWidth(2);hs[v]->SetStats(0);
        }
        hs[0]->Draw("HIST"); hs[1]->Draw("HIST SAME");
        TF1 fth("_fth","exp(-x)",0,10);fth.SetLineColor(kRed);fth.SetLineWidth(2);fth.Draw("SAME");
        TLegend lg(0.55,0.72,0.92,0.90);lg.SetBorderSize(0);lg.SetFillStyle(0);
        char lb0[60],lb1[60];
        std::snprintf(lb0,sizeof(lb0),"v7  #lambda=%.2f Hz",r0.lam_hz);
        std::snprintf(lb1,sizeof(lb1),"v8  #lambda=%.2f Hz",r1.lam_hz);
        lg.AddEntry(hs[0],lb0,"l");lg.AddEntry(hs[1],lb1,"l");
        lg.AddEntry(&fth,"Exp(1) teórico","l");lg.Draw();
    }

    //── Panel 2: Counts/intervalo v7 ─────────────────────────────────────────
    cc.cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    {
        for (int v=0;v<2;v++){
            const RunData&r=*rds[v];
            const double dtau=1000.0; const int k=(int)(r.T_total/dtau);
            if(k<2)continue;
            std::vector<int>cnt(k,0);
            for(double t:r.ts){int b=(int)(t/dtau);if(b>=0&&b<k)cnt[b]++;}
            double mu=0;for(int c:cnt)mu+=c;mu/=k;
            int nmax=*std::max_element(cnt.begin(),cnt.end());
            TH1D* hobs=new TH1D(("_co"+r.label).c_str(),
                "Counts per 1ms;N por bin;Bins",nmax+4,-0.5,nmax+3.5);
            for(int c:cnt)hobs->Fill(c);
            hobs->SetLineColor(cols[v]);hobs->SetLineWidth(2);hobs->SetStats(0);
            if(v==0)hobs->Draw("HIST");else hobs->Draw("HIST SAME");
        }
        // Poisson de v7
        {const RunData&r=r0;const double dtau=1000;const int k=(int)(r.T_total/dtau);
         if(k>2){std::vector<int>cnt(k,0);for(double t:r.ts){int b=(int)(t/dtau);if(b>=0&&b<k)cnt[b]++;}
         double mu=0;for(int c:cnt)mu+=c;mu/=k;
         int nmax=*std::max_element(cnt.begin(),cnt.end());
         TH1D*hpoi=new TH1D("_poi","",nmax+4,-0.5,nmax+3.5);
         for(int n=0;n<=nmax+3;n++){double lp=mu>0?-mu+n*std::log(mu)-TMath::LnGamma(n+1):-1e30;
             hpoi->SetBinContent(n+1,std::exp(lp)*k);}
         hpoi->SetLineColor(kRed);hpoi->SetLineStyle(2);hpoi->SetLineWidth(2);hpoi->Draw("HIST SAME");}}
        TLegend lg(0.55,0.75,0.92,0.90);lg.SetBorderSize(0);lg.SetFillStyle(0);
        for(int v=0;v<2;v++){char lb[20];std::snprintf(lb,sizeof(lb),"%s",rds[v]->label.c_str());
            lg.AddEntry((TObject*)nullptr, lb, "");}
        lg.AddEntry((TObject*)nullptr,"Poisson(#mu) v7","");lg.Draw();
    }

    //── Panel 3: D vs Δτ ─────────────────────────────────────────────────────
    cc.cd(3); gPad->SetLogx(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    {
        const std::vector<double> dtaus={100,500,1000,5000,10000,60000};
        int np=(int)dtaus.size();
        TGraphErrors* gs[2]={};
        double Dmax=1.5;
        for (int v=0;v<2;v++){
            const RunData&r=*rds[v];
            std::vector<double>xs(np),Ds(np),eDs(np);
            for(int i=0;i<np;i++){auto[D,eD]=disp_idx(r.ts,r.T_total,dtaus[i]);
                xs[i]=dtaus[i];Ds[i]=D;eDs[i]=eD;Dmax=std::max(Dmax,D);}
            gs[v]=new TGraphErrors(np,xs.data(),Ds.data(),
                std::vector<double>(np,0).data(),eDs.data());
            gs[v]->SetMarkerColor(cols[v]);gs[v]->SetLineColor(cols[v]);
            gs[v]->SetMarkerStyle(20+v);gs[v]->SetLineWidth(2);
        }
        gs[0]->SetTitle("Índice de dispersión;#Delta#tau [#mus];D = #sigma^{2}/<N>");
        gs[0]->Draw("APL"); gs[0]->GetYaxis()->SetRangeUser(0,Dmax*1.3);
        gs[1]->Draw("PL SAME");
        TLine l1(dtaus.front(),1,dtaus.back(),1);
        l1.SetLineColor(kRed);l1.SetLineWidth(2);l1.SetLineStyle(2);l1.Draw();
        TLegend lg(0.55,0.80,0.88,0.92);lg.SetBorderSize(0);lg.SetFillStyle(0);
        lg.AddEntry(gs[0],"v7","pl");lg.AddEntry(gs[1],"v8","pl");lg.Draw();
    }

    //── Panel 4: Rate vs tiempo ───────────────────────────────────────────────
    cc.cd(4); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    {
        const double win=60000.0;
        bool first=true;
        for (int v=0;v<2;v++){
            const RunData&r=*rds[v];
            double slide=std::max(win*0.1,r.T_total/500.0);
            std::vector<double>tx,ry;
            for(double t=0;t+win<=r.T_total;t+=slide){
                auto lo=std::lower_bound(r.ts.begin(),r.ts.end(),t);
                auto hi=std::lower_bound(r.ts.begin(),r.ts.end(),t+win);
                tx.push_back(t/r.T_total);  // fracción del run
                ry.push_back((int)(hi-lo)/(win*1e-6));}
            if(tx.empty())continue;
            TGraph*g=new TGraph((int)tx.size(),tx.data(),ry.data());
            g->SetTitle("Tasa vs tiempo;Fracción del run;Tasa [Hz]");
            g->SetLineColor(cols[v]);g->SetLineWidth(1);
            if(first){g->Draw("AL");g->GetYaxis()->SetRangeUser(0,3.0);first=false;}
            else g->Draw("L SAME");
        }
        TLegend lg(0.60,0.80,0.88,0.92);lg.SetBorderSize(0);lg.SetFillStyle(0);
        for(int v=0;v<2;v++){char lb[4];std::snprintf(lb,sizeof(lb),"v%d",(v==0?7:8));
            lg.AddEntry((TObject*)nullptr,"",""); }  // placeholder entries
        lg.Draw();
    }
    cc.SaveAs((out/"poisson_analysis_v7_vs_v8.png").string().c_str());
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
    const fs::path outdir=fs::path("outputs")/("poisson_analysis_"+std::string(ds));
    fs::create_directories(outdir);

    std::cout<<"=== poisson_analysis ===\nv7: "<<fv7<<"\nv8: "<<fv8<<"\nOut: "<<outdir<<"\n\n";
    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);

    // Carga y análisis básico
    RunData rd[2];
    rd[0]=load_run(fv7,"v7");
    rd[1]=load_run(fv8,"v8");

    // Ajuste exponencial + KS
    for(int v=0;v<2;v++) if(!rd[v].ts.empty()) fit_exp_ks(rd[v]);

    // Dispersión e índice χ² a Δτ=1ms
    for(int v=0;v<2;v++){
        if(rd[v].ts.empty())continue;
        auto[D,eD]=disp_idx(rd[v].ts,rd[v].T_total,1000.0);
        rd[v].D1ms=D;
        auto[chi2,pchi2]=chi2_poisson(rd[v].ts,rd[v].T_total,1000.0);
        rd[v].pchi2_1ms=pchi2;
    }

    // Imprime estadísticas básicas
    std::cout<<std::fixed<<std::setprecision(3);
    for(int v=0;v<2;v++){
        const RunData&r=rd[v];if(r.ts.empty())continue;
        std::cout<<"["<<r.label<<"]  N="<<r.n
                 <<"  T="<<r.T_total*1e-6/3600<<" h"
                 <<"  λ_hat="<<r.lam_hz<<" Hz"
                 <<"  T_mean="<<r.t_mean<<" µs"
                 <<"  τ_dead≈"<<r.tau_dead<<" µs"
                 <<"  N_gaps="<<r.gaps.size()<<"\n";
        for(const auto&[idx,dur]:r.gaps)
            std::cout<<"  Gap @event "<<idx<<"  dur="<<dur*1e-6<<"s\n";
    }
    std::cout<<"\n";

    // Genera figuras por run
    for(int v=0;v<2;v++) if(!rd[v].ts.empty()) make_figs(rd[v],outdir);
    if(!rd[0].ts.empty()&&!rd[1].ts.empty()) make_comparison(rd[0],rd[1],outdir);

    // Veredicto
    auto veredicto=[](const RunData&r)->std::string{
        if(r.ts.empty())return "SIN DATOS";
        bool ok_ks=r.p_KS>0.05,ok_chi2=r.pchi2_1ms>0.05,ok_D=std::abs(r.D1ms-1.0)<0.1;
        bool has_gaps=!r.gaps.empty();
        if(ok_ks&&ok_chi2&&ok_D&&!has_gaps)return "DISTRIBUCIÓN POISSON CONFIRMADA";
        if(!ok_ks&&r.tau_dead>r.t_mean*0.01)
            return "POISSON CON DEAD TIME ~"+std::to_string((int)r.tau_dead)+" µs";
        if(r.D1ms>1.5)return "CLUSTERING DETECTADO (D="+std::to_string(r.D1ms).substr(0,4)+")";
        if(has_gaps)return "GAPS DAQ: "+std::to_string(r.gaps.size())+" interrupción(es)";
        return "DESVIACIÓN MODERADA DEL PROCESO DE POISSON";
    };

    // Tabla summary
    auto print_table=[&](std::ostream&os){
        os<<"=== Test estadístico de Poisson — distribución temporal de triggers ===\n\n";
        os<<std::setw(33)<<std::left<<"Métrica"
          <<std::setw(18)<<"v7.root"<<std::setw(18)<<"v8.root\n";
        os<<"-------------------------------------------------------------------\n";
        auto R=[](double v,int p=3)->std::string{char b[32];std::snprintf(b,sizeof(b),"%.*f",p,v);return b;};
        auto row=[&](const std::string&m,const std::string&a,const std::string&b){
            os<<std::setw(33)<<m<<std::setw(18)<<a<<std::setw(18)<<b<<"\n";};
        row("N eventos",std::to_string(rd[0].n),std::to_string(rd[1].n));
        row("Duración [horas]",R(rd[0].T_total*1e-6/3600,2),R(rd[1].T_total*1e-6/3600,2));
        row("λ_hat [Hz]",R(rd[0].lam_hz),R(rd[1].lam_hz));
        row("λ_fit (exp) [Hz]",
            R(rd[0].lam_fit*1e6)+"±"+R(rd[0].lam_fit_err*1e6),
            R(rd[1].lam_fit*1e6)+"±"+R(rd[1].lam_fit_err*1e6));
        row("λ_esperada (geom.) [Hz]","~0.70","~0.70");
        row("λ_hat / λ_esperada",R(rd[0].lam_hz/kLambdaExp_Hz),R(rd[1].lam_hz/kLambdaExp_Hz));
        row("Dead time τ_dead [µs]",R(rd[0].tau_dead,0),R(rd[1].tau_dead,0));
        row("Fracción tiempo muerto [%]",
            R(rd[0].tau_dead*rd[0].lam_hat*100,3),
            R(rd[1].tau_dead*rd[1].lam_hat*100,3));
        row("D_KS",R(rd[0].D_KS,5),R(rd[1].D_KS,5));
        row("p-value KS",R(rd[0].p_KS),R(rd[1].p_KS));
        row("D (dispers. Δτ=1ms)",R(rd[0].D1ms),R(rd[1].D1ms));
        row("p-value χ² Poisson (Δτ=1ms)",R(rd[0].pchi2_1ms),R(rd[1].pchi2_1ms));
        row("N gaps DAQ",std::to_string(rd[0].gaps.size()),std::to_string(rd[1].gaps.size()));
        row("Veredicto",veredicto(rd[0]),veredicto(rd[1]));
        os<<"\nFiguras producidas:\n";
        os<<"  interarrival_v7/v8.png         — Histograma ΔT + fit exponencial\n";
        os<<"  qq_interarrival_v7/v8.png      — Q-Q plot vs exponencial teórica\n";
        os<<"  counts_per_interval_1ms_v7/v8  — Conteos/1ms vs Poisson esperado\n";
        os<<"  dispersion_index_vs_binwidth_*  — D vs Δτ\n";
        os<<"  rate_vs_time_v7/v8.png         — Tasa en ventanas de 60s\n";
        os<<"  poisson_analysis_v7_vs_v8.png  — Figura comparativa 4 paneles\n";
    };
    print_table(std::cout);
    {std::ofstream sf(outdir/"summary_poisson.txt");print_table(sf);}
    std::cout<<"\nOutput: "<<outdir<<"\n";
    return 0;
}
