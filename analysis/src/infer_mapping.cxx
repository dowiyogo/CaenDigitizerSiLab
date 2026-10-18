// infer_mapping.cxx
// Inferencia estadística del orden de canales hod2019 (X e Y) usando tres métodos
// independientes, con el hod2018 como plano de referencia confiable.
//
// NOTA HARDWARE: en hod2019, los cables del eje Y son bastante más largos que los del eje X.
//   Esto implica mayor atenuación de señal (menor amplitud) y posible desplazamiento temporal
//   en el eje Y respecto al X. Se analiza cada coordenada con la misma ventana temporal por
//   defecto, pero el offset Y se reporta si difiere.
//
// Método A — Calibración espacial inter-plano (principal):
//   Para cada canal c del banco hod2019X (48..63): en eventos donde c es el ÚNICO canal
//   activo en hod2019X, acumula la posición del canal más activo en hod2018X como estimador
//   x̂_c. Ordenando por x̂_c → mapeo inferido π_A.
//   Variante estricta: exactamente 1 activo en AMBOS planos (pocos eventos → reportado).
//
// Método B — Seriacón espectral por vector de Fiedler:
//   Matriz de co-activación A_ij entre pares de canales hod2019X. Vector de Fiedler del
//   Laplaciano → orden espacial.
//
// Método C — Suavidad de ocupancia:
//   R(π) = Σ(O_{k+1} - O_k)² para tres mapeos. Menor R = más físico.
//
// Build: g++ $(root-config --cflags --libs) -std=c++17 -O2 -o infer_mapping infer_mapping.cxx

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
// Parámetros configurables
//==============================================================================
struct Params {
    std::string input_file  = "lectura_vx1740_ext_nim_1M_v7.root";
    double pitch_mm         = 2.0;
    double zsep_mm          = 98.4;
    double n_sigma          = 5.0;
    double abs_thr_2018     = 20.0;  // ADC — umbral hod2018 (referencia)
    double abs_thr_2019     =  5.0;  // ADC — umbral hod2019 (señal muón ≈ 17 ADC)
    // Umbral aumentado para canal líder hod2018 en Método A relajado
    double abs_thr_2018_lead = 40.0; // ADC — requiere señal fuerte para referencia confiable
};

//==============================================================================
// Constantes del detector
//==============================================================================
constexpr double kPrePulse_ns  = 2500.0;
constexpr double kHod2019Tlo   = 3000.0,  kHod2019Thi = 5000.0;
constexpr double kHod2018Tlo   = 3500.0,  kHod2018Thi = 7000.0;
constexpr int    kNCh          = 64;

//==============================================================================
// Mapeo propuesto (π₀) — referencia con ~70% confianza de hardware
// Fuente: fotografías + etiquetado de cables. ORDEN INCIERTO dentro de cada banco.
//==============================================================================
static const std::array<int,16> kProp19X = {
    53, 61, 54, 62, 55, 63, 52, 60,
    51, 59, 50, 58, 49, 57, 48, 56
};
static const std::array<int,16> kProp19Y = {
    37, 45, 38, 46, 39, 47, 36, 44,
    35, 43, 34, 42, 33, 41, 32, 40
};

//==============================================================================
// Extracción de pulso (igual que reconstruct_v8)
//==============================================================================
struct PulseResult {
    bool   valid        = false;
    double pedestal     = 0.0;
    double noise        = 0.0;
    double pulse_height = 0.0;   // pedestal − min(ventana), positivo
    double pulse_time   = -1.0;
    bool   hit          = false;
};

static PulseResult extract_pulse(const std::vector<float>& wf,
                                  const std::vector<float>& t,
                                  double tlo, double thi,
                                  double nsig, double athr) {
    PulseResult r;
    const int n = static_cast<int>(wf.size());
    if (n < 4 || static_cast<int>(t.size()) != n) return r;
    double sum=0,sum2=0; int nb=0;
    for (int i=0;i<n;++i) if(t[i]<kPrePulse_ns){sum+=wf[i];sum2+=wf[i]*wf[i];++nb;}
    if (nb < 2) return r;
    r.valid    = true;
    r.pedestal = sum/nb;
    const double var=sum2/nb-r.pedestal*r.pedestal;
    r.noise    = (var>0)?std::sqrt(var):0.0;
    double vmin=std::numeric_limits<double>::infinity(); int imin=-1;
    for(int i=0;i<n;++i) if(t[i]>=tlo&&t[i]<=thi&&wf[i]<vmin){vmin=wf[i];imin=i;}
    if (imin<0) return r;
    r.pulse_height = r.pedestal - vmin;
    r.pulse_time   = t[imin];
    r.hit = r.pulse_height > std::max(nsig*r.noise, athr);
    return r;
}

//==============================================================================
// Utilidades de mapeo
//==============================================================================
// Inversa: canal c → índice en π₀ (-1 si no pertenece al banco)
static int prop_inv(const std::array<int,16>& prop, int ch) {
    for (int k=0;k<16;++k) if(prop[k]==ch) return k;
    return -1;
}

static double strip_coord(int k, double pitch) {
    return (k + 0.5 - 8.0) * pitch;
}

//==============================================================================
// Datos para Método A y B — coordenada X o Y
//==============================================================================
struct CoordData {
    // Banco de canales: hod2019X → 48..63 (offset 48); hod2019Y → 32..47 (offset 32)
    int base19;    // 48 o 32
    int base18;    // 0 (X) o 16 (Y)
    const std::array<int,16>* prop;  // puntero a kProp19X o kProp19Y

    // Método A — acumuladores por canal (índice = c - base19)
    long long nA[16]  = {};    // clean events (exactamente 1 en hod2019, 1 en hod2018)
    double sumX_A[16] = {};    // Σ x2018
    double sumX2_A[16]= {};    // Σ x2018²

    // Método A relajado — solo 1 activo en hod2019, canal líder hod2018 (ph > umbral_lead)
    long long nR[16]  = {};
    double sumX_R[16] = {};
    double sumX2_R[16]= {};

    // Scatter data (para figuras): (canal_idx_19, x18)
    std::vector<std::pair<int,double>> scatterA, scatterR;

    // Método B — co-activación entre pares de canales hod2019 (índice = c - base19)
    long long coact[16][16] = {};
    long long hits19[16]    = {};   // eventos donde cada canal hod2019 está activo

    // Ocupancia hod2019 (total hits por canal, cualquier evento)
    long long occ19[16] = {};

    // Ocupancia hod2018 (total hits, referencia)
    long long occ18[16] = {};

    // Tiempo de pulso medio por canal hod2019 (para diagnóstico de cables Y)
    double sum_pt19[16] = {};
    long long n_pt19[16]= {};
};

//==============================================================================
// Método B — vector de Fiedler del Laplaciano de la matriz de afinidad
//==============================================================================
static std::vector<double> compute_fiedler(const long long coact[16][16],
                                            const long long hits[16]) {
    const int n = 16;
    TMatrixDSym L(n);
    // Construye matriz de afinidad normalizada A_ij
    double A[16][16] = {};
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) {
        if (i==j) continue;
        const double denom = std::sqrt(static_cast<double>(hits[i]) *
                                        static_cast<double>(hits[j]));
        A[i][j] = (denom > 0.5) ? static_cast<double>(coact[i][j]) / denom : 0.0;
    }
    // Laplaciano L = D - A
    for (int i=0;i<n;++i) {
        double di = 0;
        for (int j=0;j<n;++j) di += A[i][j];
        for (int j=0;j<n;++j) L[i][j] = (i==j) ? di : -A[i][j];
    }
    TMatrixDSymEigen eigen(L);
    const TVectorD& evals  = eigen.GetEigenValues();
    const TMatrixD& evecs  = eigen.GetEigenVectors();

    // Encuentra el segundo autovalor más pequeño (saltar el trivial ≈ 0)
    int fiedler_col = -1;
    double second_smallest = std::numeric_limits<double>::infinity();
    double smallest = std::numeric_limits<double>::infinity();
    for (int k=0;k<n;++k) {
        if (evals(k) < smallest) { smallest = evals(k); }
    }
    for (int k=0;k<n;++k) {
        if (evals(k) > smallest + 1e-9 && evals(k) < second_smallest) {
            second_smallest = evals(k); fiedler_col = k;
        }
    }
    std::vector<double> fvec(n, 0.0);
    if (fiedler_col >= 0) {
        for (int i=0;i<n;++i) fvec[i] = evecs[i][fiedler_col];
    }
    return fvec;
}

//==============================================================================
// Método C — Rugosidad de ocupancia
//==============================================================================
// mapping[k] = canal_offset (0..15) para el índice k
static double roughness(const long long occ19[], const int mapping[16]) {
    double R = 0;
    for (int k=0;k<15;++k) {
        double diff = static_cast<double>(occ19[mapping[k+1]]) -
                      static_cast<double>(occ19[mapping[k]]);
        R += diff*diff;
    }
    return R;
}

//==============================================================================
// Kendall τ entre dos ordenaciones de 16 canales (ignorando canales sin datos)
// rank_a[offset] = índice asignado a canal (offset+base) por ordenación A
// rank_b[offset] = ídem por B
// Considera sólo los canales con datos en AMBAS ordenaciones (n_a[offset]>0 && n_b[offset]>0)
//==============================================================================
static double kendall_tau(const int rank_a[16], const int rank_b[16],
                          const long long valid[16]) {
    std::vector<int> chans;
    for (int i=0;i<16;++i) if(valid[i]>0 && rank_a[i]>=0 && rank_b[i]>=0) chans.push_back(i);
    const int m = static_cast<int>(chans.size());
    if (m < 2) return 0.0;
    long long C=0, D=0;
    for (int i=0;i<m;++i) for(int j=i+1;j<m;++j) {
        int da = rank_a[chans[i]] - rank_a[chans[j]];
        int db = rank_b[chans[i]] - rank_b[chans[j]];
        if (da*db > 0) ++C; else if (da*db < 0) ++D;
    }
    const long long total = static_cast<long long>(m)*(m-1)/2;
    return total > 0 ? static_cast<double>(C-D)/static_cast<double>(total) : 0.0;
}

//==============================================================================
// Gráficas auxiliares
//==============================================================================
static void style1(TH1D* h, const char* xt, const char* yt) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
    h->SetLineWidth(2); h->SetStats(0);
}
static void style2(TH2D* h, const char* xt, const char* yt) {
    h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
    h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
    h->SetStats(0);
}
static void png1(TH1D* h, const fs::path& p) {
    TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",1000,680);
    c.SetLeftMargin(0.12); c.SetBottomMargin(0.13); c.SetRightMargin(0.05);
    h->Draw("HIST"); c.SaveAs(p.string().c_str());
}
static void png2(TH2D* h, const fs::path& p) {
    TCanvas c(("c_"+std::string(h->GetName())).c_str(),"",900,820);
    c.SetLeftMargin(0.13); c.SetRightMargin(0.17); c.SetBottomMargin(0.12);
    h->Draw("COLZ"); c.SaveAs(p.string().c_str());
}

//==============================================================================
// Análisis completo para una coordenada (X o Y)
// Returns inferred ranking arrays and summary text
//==============================================================================
struct CoordResult {
    std::string coord;          // "X" o "Y"
    int   rank_A[16];          // rank_A[offset] = índice inferido (π_A) para canal offset+base
    int   rank_B[16];          // π_B
    int   rank_0[16];          // π₀ (propuesto)
    double x_hat[16];          // posición estimada [mm]
    double err[16];            // error estándar [mm]
    long long n_clean[16];     // eventos limpios por canal
    long long n_relax[16];     // eventos relajados por canal
    std::vector<double> fiedler;
    double tau_A0;             // Kendall τ(π_A, π₀)
    double tau_B0;             // Kendall τ(π_B, π₀)
    double tau_AB;             // Kendall τ(π_A, π_B)
    double roughness_0;
    double roughness_A;
    double roughness_B;
    double roughness_rand;
    double r_pearson_0;        // Pearson r con mapeo propuesto
    double r_pearson_A;        // Pearson r con mapeo inferido π_A
};

static CoordResult analyze_coord(CoordData& d, const std::string& coord,
                                  double pitch, const fs::path& outdir,
                                  TFile* rootout) {
    CoordResult res;
    res.coord = coord;

    // ── Construye inversas del mapeo propuesto ─────────────────────────────
    for (int c=0;c<16;++c) {
        int k0 = prop_inv(*d.prop, c + d.base19);
        res.rank_0[c] = k0;  // proposed index for channel (c + base19)
    }

    // ── Método A: ranking por x̂_c (relajado: mayor estadística) ──────────
    // Para cada canal c (offset 0..15): x̂_c = sumX_R[c] / nR[c]
    // Usamos datos relajados como principal (más estadística)
    for (int c=0;c<16;++c) {
        res.n_clean[c] = d.nA[c];
        res.n_relax[c] = d.nR[c];
        if (d.nR[c] > 0) {
            res.x_hat[c] = d.sumX_R[c] / d.nR[c];
            const double var = d.sumX2_R[c]/d.nR[c] - res.x_hat[c]*res.x_hat[c];
            const double sigma = (var>0) ? std::sqrt(var) : 0.0;
            res.err[c] = sigma / std::sqrt(static_cast<double>(d.nR[c]));
        } else {
            res.x_hat[c] = std::numeric_limits<double>::quiet_NaN();
            res.err[c]   = std::numeric_limits<double>::infinity();
        }
    }

    // Ordena los 16 canales por x̂_c para asignar índice π_A
    std::vector<int> order_A(16);
    std::iota(order_A.begin(), order_A.end(), 0);
    std::stable_sort(order_A.begin(), order_A.end(), [&](int a, int b) {
        if (d.nR[a] == 0 && d.nR[b] == 0) return a < b;
        if (d.nR[a] == 0) return false;
        if (d.nR[b] == 0) return true;
        return res.x_hat[a] < res.x_hat[b];
    });
    // order_A[k] = canal_offset con índice k en π_A
    // Inversa: rank_A[offset] = índice asignado
    for (int c=0;c<16;++c) res.rank_A[c] = -1;
    for (int k=0;k<16;++k) res.rank_A[order_A[k]] = k;

    // ── Método B: Fiedler ──────────────────────────────────────────────────
    res.fiedler = compute_fiedler(d.coact, d.hits19);

    // Ordena canales por Fiedler
    std::vector<int> order_B(16);
    std::iota(order_B.begin(), order_B.end(), 0);
    std::sort(order_B.begin(), order_B.end(), [&](int a, int b) {
        return res.fiedler[a] < res.fiedler[b];
    });

    // Elige orientación para que τ(π_B, π_A) > 0
    {
        // Compute ranks for both orientations and pick the one with τ > 0
        int tmpB[16], tmpB_rev[16];
        for (int k=0;k<16;++k) { tmpB[order_B[k]] = k; tmpB_rev[order_B[15-k]] = k; }
        // Compare with rank_A
        // Simple check: does forward orientation give positive sum?
        int score = 0;
        for (int c=0;c<16;++c) if(d.nR[c]>0 && res.rank_A[c]>=0)
            score += (tmpB[c] - tmpB_rev[c]) * (res.rank_A[c] - 7);
        if (score < 0) {
            // Reverse
            std::reverse(order_B.begin(), order_B.end());
        }
    }
    for (int c=0;c<16;++c) res.rank_B[c] = -1;
    for (int k=0;k<16;++k) res.rank_B[order_B[k]] = k;

    // ── Kendall τ ──────────────────────────────────────────────────────────
    res.tau_A0 = kendall_tau(res.rank_A, res.rank_0, d.nR);
    res.tau_B0 = kendall_tau(res.rank_B, res.rank_0, d.nR);
    res.tau_AB = kendall_tau(res.rank_A, res.rank_B, d.nR);

    // ── Método C: rugosidad ────────────────────────────────────────────────
    // Construye mappings como arrays de offsets
    int map0[16], mapA[16], mapB[16];
    for (int k=0;k<16;++k) map0[k] = prop_inv(*d.prop, (*d.prop)[k]) ; // identity by construction
    for (int k=0;k<16;++k) {
        // prop: kProp19X[k] = channel → offset = kProp19X[k] - base
        map0[k] = (*d.prop)[k] - d.base19;
        mapA[k] = order_A[k];
        mapB[k] = order_B[k];
    }
    res.roughness_0 = roughness(d.occ19, map0);
    res.roughness_A = roughness(d.occ19, mapA);
    res.roughness_B = roughness(d.occ19, mapB);

    // Rugosidad media de 1000 permutaciones aleatorias
    std::mt19937 rng(42);
    std::vector<int> perm(16); std::iota(perm.begin(), perm.end(), 0);
    double rand_sum = 0;
    for (int trial=0; trial<1000; ++trial) {
        std::shuffle(perm.begin(), perm.end(), rng);
        rand_sum += roughness(d.occ19, perm.data());
    }
    res.roughness_rand = rand_sum / 1000.0;

    // ── Pearson r para mapeo propuesto vs inferido ─────────────────────────
    // Acumula x2018 vs x2019(prop) y vs x2019(inf_A) para scatter
    {
        double sx18=0,sx18_2=0,sxp=0,sxp_2=0,sxp_18=0;
        double sxa=0,sxa_2=0,sxa_18=0;
        int nr=0;
        for (const auto& [idx,x18] : d.scatterR) {
            int k0 = res.rank_0[idx];  // proposed index
            int kA = res.rank_A[idx];  // inferred index A
            if (k0<0 || kA<0) continue;
            double x0 = strip_coord(k0, pitch);
            double xA = strip_coord(kA, pitch);
            sx18+=x18; sxp+=x0; sxa+=xA;
            sx18_2+=x18*x18; sxp_2+=x0*x0; sxa_2+=xA*xA;
            sxp_18+=x0*x18; sxa_18+=xA*x18;
            ++nr;
        }
        auto pearson = [&](double sxy, double sx, double sy, double sxx, double syy) -> double {
            if (nr<2) return 0;
            double num = nr*sxy - sx*sy;
            double den = std::sqrt((nr*sxx-sx*sx)*(nr*syy-sy*sy));
            return den>0 ? num/den : 0.0;
        };
        res.r_pearson_0 = pearson(sxp_18, sxp, sx18, sxp_2, sx18_2);
        res.r_pearson_A = pearson(sxa_18, sxa, sx18, sxa_2, sx18_2);
    }

    // ── Figuras ────────────────────────────────────────────────────────────
    const std::string cx = coord;
    constexpr double lo = -16.0, hi = 16.0;

    // Scatter plots (x2018 vs x2019) — propuesto y inferido
    auto h_scatter_0 = std::make_unique<TH2D>(
        ("h_scatter_"+cx+"_prop").c_str(),
        ("Scatter "+cx+" propuesto;x_{2018} [mm];x_{2019}(#pi_0) [mm];Eventos").c_str(),
        16,lo,hi,16,lo,hi);
    auto h_scatter_A = std::make_unique<TH2D>(
        ("h_scatter_"+cx+"_inf").c_str(),
        ("Scatter "+cx+" inferido;x_{2018} [mm];x_{2019}(#pi_A) [mm];Eventos").c_str(),
        16,lo,hi,16,lo,hi);
    for (const auto& [idx,x18] : d.scatterR) {
        int k0=res.rank_0[idx], kA=res.rank_A[idx];
        if(k0>=0) h_scatter_0->Fill(x18, strip_coord(k0,pitch));
        if(kA>=0) h_scatter_A->Fill(x18, strip_coord(kA,pitch));
    }
    style2(h_scatter_0.get(),"x_{2018} [mm]","x_{2019} [mm]");
    style2(h_scatter_A.get(),"x_{2018} [mm]","x_{2019} [mm]");
    png2(h_scatter_0.get(), outdir/("scatter_x2019_vs_x2018_propuesto_"+cx+".png"));
    {
        TCanvas c2(("cs_"+cx+"_inf").c_str(),"",900,820);
        c2.SetLeftMargin(0.13); c2.SetRightMargin(0.17); c2.SetBottomMargin(0.12);
        h_scatter_A->Draw("COLZ");
        TLatex lat; lat.SetTextSize(0.035); lat.SetNDC(true);
        char buf[80]; std::snprintf(buf,sizeof(buf),"r(inf)=%.3f  r(prop)=%.3f",
                                    res.r_pearson_A, res.r_pearson_0);
        lat.DrawLatex(0.14,0.93,buf);
        c2.SaveAs((outdir/("scatter_x2019_vs_x2018_inferido_"+cx+".png")).string().c_str());
    }

    // Posición estimada ± error por canal
    {
        std::vector<double> xs,ys,exs,eys;
        for (int k=0;k<16;++k) {
            int c_off = order_A[k];
            if (d.nR[c_off] < 1) continue;
            xs.push_back(k);
            ys.push_back(res.x_hat[c_off]);
            exs.push_back(0.3);
            eys.push_back(res.err[c_off]);
        }
        TGraphErrors g(static_cast<int>(xs.size()), xs.data(), ys.data(),
                       exs.data(), eys.data());
        g.SetTitle(("Posici\xC3\xB3n estimada hod2019 "+cx+";#\xC3\xADndice inferido k;x\xCC\x82_c [mm]").c_str());
        g.SetMarkerStyle(20); g.SetMarkerColor(kBlue+1); g.SetLineColor(kBlue+1);
        TCanvas cg(("cg_pos_"+cx).c_str(),"",900,680);
        cg.SetLeftMargin(0.12); cg.SetBottomMargin(0.12);
        g.GetXaxis()->SetTitle("Índice inferido k");
        g.GetYaxis()->SetTitle("#hat{x}_{c} [mm]");
        g.Draw("AP");
        // Línea de referencia uniforme
        TLine ref(-0.5, strip_coord(0,pitch), 15.5, strip_coord(15,pitch));
        ref.SetLineColor(kRed); ref.SetLineStyle(2); ref.SetLineWidth(2);
        ref.Draw("SAME");
        cg.SaveAs((outdir/("posicion_estimada_"+cx+".png")).string().c_str());
    }

    // Matriz de afinidad reordenada por π_A
    {
        auto h_aff = std::make_unique<TH2D>(
            ("h_aff_"+cx).c_str(),
            ("Afinidad co-activaci\xC3\xB3n hod2019 "+cx+" (reord. \xCF\x80_A)"
             ";Canal (en orden #pi_A);Canal (en orden #pi_A);A_{ij}").c_str(),
            16,-.5,15.5,16,-.5,15.5);
        for (int ki=0;ki<16;++ki) for(int kj=0;kj<16;++kj) {
            int ci=order_A[ki], cj=order_A[kj];
            double denom = std::sqrt(static_cast<double>(d.hits19[ci])*
                                     static_cast<double>(d.hits19[cj]));
            double aij = (denom>0) ? static_cast<double>(d.coact[ci][cj])/denom : 0;
            h_aff->SetBinContent(ki+1,kj+1,aij);
        }
        // Etiquetas: canal número
        for(int k=0;k<16;++k) {
            char l[8]; std::snprintf(l,sizeof(l),"Ch%d",order_A[k]+d.base19);
            h_aff->GetXaxis()->SetBinLabel(k+1,l);
            h_aff->GetYaxis()->SetBinLabel(k+1,l);
        }
        h_aff->SetStats(0);
        TCanvas caff(("caff_"+cx).c_str(),"",900,820);
        caff.SetLeftMargin(0.13); caff.SetRightMargin(0.17); caff.SetBottomMargin(0.17);
        h_aff->Draw("COLZ");
        caff.SaveAs((outdir/("afinidad_matrix_"+cx+".png")).string().c_str());
        if (rootout && !rootout->IsZombie()) h_aff->Write();
    }

    // Comparación de ocupancia: π₀ vs π_A
    {
        auto h_occ_0 = std::make_unique<TH1D>(("h_occ_"+cx+"_prop").c_str(),
            ("Ocupancia hod2019 "+cx+";Índice k;Hits").c_str(),16,-0.5,15.5);
        auto h_occ_A = std::make_unique<TH1D>(("h_occ_"+cx+"_inf").c_str(),
            ("Ocupancia hod2019 "+cx+";Índice k;Hits").c_str(),16,-0.5,15.5);
        for (int k=0;k<16;++k) {
            h_occ_0->SetBinContent(k+1, d.occ19[map0[k]]);
            h_occ_A->SetBinContent(k+1, d.occ19[mapA[k]]);
        }
        TCanvas cocc(("cocc_"+cx).c_str(),"",1000,680);
        cocc.SetLeftMargin(0.12); cocc.SetBottomMargin(0.12);
        h_occ_0->SetLineColor(kRed);  h_occ_0->SetLineWidth(2);
        h_occ_A->SetLineColor(kBlue+1); h_occ_A->SetLineWidth(2); h_occ_A->SetLineStyle(2);
        h_occ_0->SetStats(0);
        h_occ_0->SetTitle(("Ocupancia "+cx+";Índice k;Hits").c_str());
        h_occ_0->Draw("HIST"); h_occ_A->Draw("HIST SAME");
        TLatex lat; lat.SetNDC(true); lat.SetTextSize(0.032);
        lat.SetTextColor(kRed);   lat.DrawLatex(0.65,0.85,"#pi_0 (propuesto)");
        lat.SetTextColor(kBlue+1);lat.DrawLatex(0.65,0.80,"#pi_A (inferido)");
        cocc.SaveAs((outdir/("ocupancia_comparacion_"+cx+".png")).string().c_str());
        if(rootout&&!rootout->IsZombie()){h_occ_0->Write();h_occ_A->Write();}
    }

    // Fiedler vector
    {
        auto h_fied = std::make_unique<TH1D>(("h_fiedler_"+cx).c_str(),
            ("Vector de Fiedler hod2019 "+cx+";Canal offset;Valor Fiedler").c_str(),
            16,-0.5,15.5);
        for (int c=0;c<16;++c) {
            h_fied->SetBinContent(c+1, res.fiedler[c]);
            char l[8]; std::snprintf(l,sizeof(l),"Ch%d",c+d.base19);
            h_fied->GetXaxis()->SetBinLabel(c+1,l);
        }
        h_fied->SetLineWidth(2); h_fied->SetStats(0);
        TCanvas cfied(("cfied_"+cx).c_str(),"",1000,680);
        cfied.SetLeftMargin(0.12); cfied.SetBottomMargin(0.17);
        h_fied->Draw("HIST"); h_fied->LabelsOption("v","X");
        cfied.SaveAs((outdir/("fiedler_vector_"+cx+".png")).string().c_str());
        if(rootout&&!rootout->IsZombie()) h_fied->Write();
    }

    return res;
}

//==============================================================================
// Main
//==============================================================================
int main(int argc, char** argv) {
    Params p;
    for (int i=1;i<argc;++i) {
        std::string a=argv[i];
        if(a=="--pitch"   &&i+1<argc){p.pitch_mm         =std::stod(argv[++i]);continue;}
        if(a=="--nsigma"  &&i+1<argc){p.n_sigma           =std::stod(argv[++i]);continue;}
        if(a=="--abst18"  &&i+1<argc){p.abs_thr_2018      =std::stod(argv[++i]);continue;}
        if(a=="--abst19"  &&i+1<argc){p.abs_thr_2019      =std::stod(argv[++i]);continue;}
        if(a=="--abslead" &&i+1<argc){p.abs_thr_2018_lead =std::stod(argv[++i]);continue;}
        if(a[0]!='-'){p.input_file=a;continue;}
        std::cerr<<"Unknown option: "<<a<<'\n'; return 1;
    }

    std::time_t now=std::time(nullptr); char ds[32];
    std::strftime(ds,sizeof(ds),"%Y%m%d_%H%M%S",std::localtime(&now));
    const fs::path outdir=fs::path("outputs")/("mapping_inference_"+std::string(ds));
    fs::create_directories(outdir);

    std::cout<<"=== infer_mapping ===\n";
    std::cout<<"Input : "<<p.input_file<<"\n";
    std::cout<<"Output: "<<outdir<<"\n";
    std::cout<<"abs_thr_2018="<<p.abs_thr_2018<<" ADC  abs_thr_2019="<<p.abs_thr_2019<<" ADC"
             <<"  abs_lead="<<p.abs_thr_2018_lead<<" ADC  nsigma="<<p.n_sigma<<"\n\n";

    gROOT->SetBatch(kTRUE); TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0); gStyle->SetPalette(kViridis);

    const std::string resolved = hodo::resolveRawFile(p.input_file);
    TFile* fin=TFile::Open(resolved.c_str(),"READ");
    if(!fin||fin->IsZombie()){std::cerr<<"ERROR: cannot open "<<resolved<<'\n';return 1;}
    TTree* tree=dynamic_cast<TTree*>(fin->Get("data_dig"));
    if(!tree){std::cerr<<"ERROR: TTree 'data_dig' not found\n";return 1;}

    std::array<std::vector<float>*,kNCh> ch_wf{};
    std::vector<float>* time_wf=nullptr;
    tree->SetBranchStatus("*",0); tree->SetCacheSize(64LL*1024LL*1024LL);
    for(int c=0;c<kNCh;++c){
        char bn[16]; std::snprintf(bn,sizeof(bn),"Ch%d",c);
        if(!tree->GetBranch(bn)) continue;
        tree->SetBranchStatus(bn,1); tree->SetBranchAddress(bn,&ch_wf[c]);
    }
    tree->SetBranchStatus("time",1); tree->SetBranchAddress("time",&time_wf);

    const Long64_t nEnt=tree->GetEntries();
    std::cout<<"Total entries: "<<nEnt<<'\n';

    // Inicializar datos de análisis
    CoordData dX; dX.base19=48; dX.base18=0;  dX.prop=&kProp19X;
    CoordData dY; dY.base19=32; dY.base18=16; dY.prop=&kProp19Y;

    //==========================================================================
    // Bucle de eventos (único paso)
    //==========================================================================
    std::cout<<"Processing events...\n";
    for(Long64_t iev=0;iev<nEnt;++iev){
        tree->GetEntry(iev);
        if(iev%2000==0) std::cout<<"  "<<iev<<" / "<<nEnt<<"\r"<<std::flush;
        if(!time_wf||time_wf->empty()) continue;
        const auto& tv=*time_wf;

        // Extrae pulsos de los 64 canales
        std::array<PulseResult,kNCh> pr{};
        for(int c=0;c<kNCh;++c){
            if(!ch_wf[c]||ch_wf[c]->empty()) continue;
            const bool is19 = (c>=32);
            const double tlo = is19 ? kHod2019Tlo : kHod2018Tlo;
            const double thi = is19 ? kHod2019Thi : kHod2018Thi;
            const double thr = is19 ? p.abs_thr_2019 : p.abs_thr_2018;
            pr[c] = extract_pulse(*ch_wf[c],tv,tlo,thi,p.n_sigma,thr);
        }

        // Agrupa activos por plano (offset dentro del banco)
        std::vector<int> act19x, act19y, act18x, act18y;
        for(int c=0;c<kNCh;++c){
            if(!pr[c].hit) continue;
            if(c>=48&&c<=63) act19x.push_back(c-48);
            else if(c>=32&&c<=47) act19y.push_back(c-32);
            else if(c>=0&&c<=15) act18x.push_back(c);
            else if(c>=16&&c<=31) act18y.push_back(c-16);
        }

        // Función auxiliar para procesar una coordenada
        auto process_coord = [&](CoordData& d,
                                  const std::vector<int>& act19,
                                  const std::vector<int>& act18_ref) {
            const int n19 = static_cast<int>(act19.size());
            const int n18 = static_cast<int>(act18_ref.size());

            // Método B: co-activación hod2019 (cualquier event con ≥1 activo)
            for(int idx : act19) {
                d.hits19[idx]++;
                d.occ19[idx]++;
                if(pr[idx+d.base19].hit && pr[idx+d.base19].pulse_time > 0) {
                    d.sum_pt19[idx]  += pr[idx+d.base19].pulse_time;
                    d.n_pt19[idx]++;
                }
            }
            for(int i : act18_ref) d.occ18[i]++;
            for(int i=0;i<n19;++i)
                for(int j=i+1;j<n19;++j) {
                    d.coact[act19[i]][act19[j]]++;
                    d.coact[act19[j]][act19[i]]++;
                }

            // Método A — estricto: exactamente 1 en hod2019 Y 1 en hod2018
            if(n19==1 && n18==1) {
                const int c19 = act19[0];
                const int c18 = act18_ref[0];
                const double x18 = strip_coord(c18, p.pitch_mm);
                d.nA[c19]++;
                d.sumX_A[c19]  += x18;
                d.sumX2_A[c19] += x18*x18;
            }

            // Método A — relajado: exactamente 1 en hod2019, canal líder en hod2018 (ph > umbral_lead)
            if(n19==1) {
                const int c19 = act19[0];
                // Encuentra canal líder en hod2018 (máxima amplitud)
                int best18=-1; double best_ph=0;
                for(int c18 : act18_ref) {
                    const double ph = pr[c18+d.base18].pulse_height;
                    if(ph > best_ph){best_ph=ph; best18=c18;}
                }
                if(best18>=0 && best_ph >= p.abs_thr_2018_lead) {
                    const double x18 = strip_coord(best18, p.pitch_mm);
                    d.nR[c19]++;
                    d.sumX_R[c19]  += x18;
                    d.sumX2_R[c19] += x18*x18;
                    d.scatterR.push_back({c19, x18});
                }
            }
        };

        process_coord(dX, act19x, act18x);
        process_coord(dY, act19y, act18y);
    }
    std::cout<<"\nEvent loop done.\n\n";

    // Archivo ROOT de salida
    TFile* rootout = TFile::Open((outdir/"histos_mapping.root").string().c_str(),"RECREATE");

    // Analiza cada coordenada
    CoordResult rX = analyze_coord(dX, "X", p.pitch_mm, outdir, rootout);
    CoordResult rY = analyze_coord(dY, "Y", p.pitch_mm, outdir, rootout);

    if(rootout&&!rootout->IsZombie()) rootout->Close();

    //==========================================================================
    // Escribe CSVs de mapeo
    //==========================================================================
    // mapping_inferred.csv (π_A)
    {
        std::ofstream f(outdir/"mapping_inferred.csv");
        f<<"plane,coord,index,channel\n";
        // Para hod2019X: construir array inferido
        std::vector<int> ord_x(16,0), ord_y(16,0);
        for(int c=0;c<16;++c) if(rX.rank_A[c]>=0) ord_x[rX.rank_A[c]]=c+dX.base19;
        for(int c=0;c<16;++c) if(rY.rank_A[c]>=0) ord_y[rY.rank_A[c]]=c+dY.base19;
        for(int k=0;k<16;++k) f<<"hod2019,X,"<<k<<','<<ord_x[k]<<'\n';
        for(int k=0;k<16;++k) f<<"hod2019,Y,"<<k<<','<<ord_y[k]<<'\n';
    }
    // mapping_proposed.csv (π₀)
    {
        std::ofstream f(outdir/"mapping_proposed.csv");
        f<<"plane,coord,index,channel\n";
        for(int k=0;k<16;++k) f<<"hod2019,X,"<<k<<','<<kProp19X[k]<<'\n';
        for(int k=0;k<16;++k) f<<"hod2019,Y,"<<k<<','<<kProp19Y[k]<<'\n';
    }

    //==========================================================================
    // Summary a stdout y archivo
    //==========================================================================
    auto print_summary = [&](std::ostream& os) {
        os<<"=== INFERENCIA DE MAPEO hod2019 ===\n\n";
        os<<"Parámetros:\n";
        os<<"  input           = "<<p.input_file<<"\n";
        os<<"  pitch_mm        = "<<p.pitch_mm<<"\n";
        os<<"  n_sigma         = "<<p.n_sigma<<"\n";
        os<<"  abs_thr_2018    = "<<p.abs_thr_2018<<" ADC\n";
        os<<"  abs_thr_2019    = "<<p.abs_thr_2019<<" ADC\n";
        os<<"  abs_thr_lead18  = "<<p.abs_thr_2018_lead<<" ADC  (umbral canal líder hod2018 en Método A relajado)\n\n";

        os<<"NOTA HARDWARE: los cables del eje Y en hod2019 son bastante más largos que los del eje X.\n";
        os<<"  Esto implica mayor atenuación de señal (menor amplitud en hod2019Y vs hod2019X),\n";
        os<<"  potencial desplazamiento temporal (llegan más tarde al digitizador), y posiblemente\n";
        os<<"  mayor ruido. Las ventanas temporales se aplican igual a ambos ejes; si hod2019Y\n";
        os<<"  muestra pocas señales, ampliar la ventana Y o verificar el offset.\n\n";
        // Timing check
        for (auto* dd : {&dX, &dY}) {
            const char* cc = (dd==&dX) ? "X" : "Y";
            os<<"  hod2019"<<cc<<" tiempo medio de pulso por canal:\n";
            for(int c=0;c<16;++c) {
                if(dd->n_pt19[c]>5)
                    os<<"    Ch"<<(c+dd->base19)<<" (k_prop="<<prop_inv(*dd->prop,c+dd->base19)
                      <<"): t_mean="<<std::fixed<<std::setprecision(0)
                      <<dd->sum_pt19[c]/dd->n_pt19[c]<<" ns  (n="<<dd->n_pt19[c]<<")\n";
            }
        }
        os<<"\n";

        for (const CoordResult& res : {rX, rY}) {
            const std::string cx = res.coord;
            os<<"================================================================\n";
            os<<"HOD2019 COORDENADA "<<cx<<"\n";
            os<<"================================================================\n\n";

            os<<"Tabla Método A (relajado): x̂_c ± err_c vs índice propuesto vs inferido\n";
            os<<std::setw(6)<<"Canal"<<std::setw(8)<<"x̂_c"<<std::setw(8)<<"±err"
              <<std::setw(7)<<"n_rel"<<std::setw(7)<<"n_cl"<<std::setw(7)<<"π_A"
              <<std::setw(7)<<"π₀"<<"  ¿Ambiguo?\n";
            os<<"---------------------------------------------------------------\n";

            // Ordena por π_A para mostrar en orden físico
            std::vector<int> order_disp(16);
            for(int c=0;c<16;++c) order_disp[c]=c;
            std::sort(order_disp.begin(),order_disp.end(),[&](int a,int b){
                return res.rank_A[a] < res.rank_A[b];
            });
            auto& d = (cx=="X") ? dX : dY;

            for(int c_off : order_disp) {
                int ch = c_off + ((cx=="X") ? dX.base19 : dY.base19);
                long long nr = res.n_relax[c_off];
                long long nc = res.n_clean[c_off];
                int kA = res.rank_A[c_off];
                int k0 = res.rank_0[c_off];
                os<<std::setw(6)<<("Ch"+std::to_string(ch));
                if(nr>0) {
                    os<<std::fixed<<std::setprecision(2)
                      <<std::setw(8)<<res.x_hat[c_off]
                      <<std::setw(8)<<res.err[c_off];
                } else {
                    os<<std::setw(8)<<"N/A"<<std::setw(8)<<"---";
                }
                os<<std::setw(7)<<nr<<std::setw(7)<<nc
                  <<std::setw(7)<<kA<<std::setw(7)<<k0;
                // Ambigüedad: ¿algún otro canal tiene x̂ solapando dentro de 2σ?
                bool ambig=false;
                if(nr>0) {
                    for(int c2=0;c2<16;++c2) {
                        if(c2==c_off||res.n_relax[c2]==0) continue;
                        if(std::abs(res.x_hat[c_off]-res.x_hat[c2])
                           < 2*(res.err[c_off]+res.err[c2])) {ambig=true;break;}
                    }
                }
                os<<"  "<<(ambig?"AMB":"");
                os<<'\n';
            }
            os<<'\n';

            os<<"Concordancias (Kendall τ):\n";
            os<<std::fixed<<std::setprecision(4);
            os<<"  τ(π_A, π₀) = "<<res.tau_A0<<"  [Mét.A vs propuesto]\n";
            os<<"  τ(π_B, π₀) = "<<res.tau_B0<<"  [Mét.B vs propuesto]\n";
            os<<"  τ(π_A, π_B) = "<<res.tau_AB<<"  [Mét.A vs Mét.B]\n\n";

            os<<"Rugosidades R = Σ(O_{k+1} - O_k)²:\n";
            os<<"  R(π₀)        = "<<std::fixed<<std::setprecision(0)<<res.roughness_0<<"\n";
            os<<"  R(π_A)       = "<<res.roughness_A<<"\n";
            os<<"  R(π_B)       = "<<res.roughness_B<<"\n";
            os<<"  R(aleatorio) = "<<res.roughness_rand<<" (media 1000 permutaciones)\n\n";

            os<<"Correlación Pearson (eventos relajados Mét.A):\n";
            os<<std::fixed<<std::setprecision(4);
            os<<"  r(x2019(π₀), x2018) = "<<res.r_pearson_0<<"\n";
            os<<"  r(x2019(π_A), x2018) = "<<res.r_pearson_A<<"\n\n";
        }

        os<<"hod2019 X: tau_A0 = "<<rX.tau_A0<<"  tau_AB = "<<rX.tau_AB<<"\n";
        os<<"hod2019 Y: tau_A0 = "<<rY.tau_A0<<"  tau_AB = "<<rY.tau_AB<<"\n\n";

        os<<"Mapeo inferido π_A (recomendado si τ_AB > 0.9):\n";
        os<<"  hod2019 X: {";
        for(int k=0;k<16;++k){
            for(int c=0;c<16;++c) if(rX.rank_A[c]==k){os<<dX.base19+c;if(k<15)os<<",";break;}
        }
        os<<"}\n  hod2019 Y: {";
        for(int k=0;k<16;++k){
            for(int c=0;c<16;++c) if(rY.rank_A[c]==k){os<<dY.base19+c;if(k<15)os<<",";break;}
        }
        os<<"}\n";
    };

    print_summary(std::cout);
    {std::ofstream sf(outdir/"summary_mapping.txt"); print_summary(sf);}

    std::cout<<"\nOutput: "<<outdir<<"\n";
    fin->Close();
    return 0;
}
