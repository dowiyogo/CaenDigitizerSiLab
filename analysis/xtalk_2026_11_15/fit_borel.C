// fit_borel.C — seguimiento 2026-11-15: diafonía inmediata de hod2018 con un modelo de Borel, una definición para mayo y octubre.
// Entrada: <prefijo>_dark3.bin de dark3_extract.C (pulsos oscuros aislados, misma amplitud que aS), solo la región pre-trigger [16,37), y una
// máscara de eventos (uint8 por evento escrito; mayo = tiempo bueno de v8, octubre = todos).
// Modelo de la amplitud a (mV) por canal, histograma en bins de 0.2 mV en [0.55 G, 6.6 G] (picos de 1 a 6 p.e.; G = la del análisis, solo para
// la ventana y las semillas):
//   dN/da = A · Σ_{N=1}^{14} P(N) · Gauss(a; N·Gf + off (+ off1 si N = 1), σ_N² = s0² + N·s1²), integrada en cada bin (erf),
//   P(N) = (1 − q)·B(N; p) + q·B2(N; p), B(N; p) = e^{−pN} (pN)^{N−1} / N! (Borel: N p.e. por p.e. primario con diafonía inmediata),
//   B2 = Borel–Tanner con 2 primarios, 2/N · e^{−pN} (pN)^{N−2} / (N−2)! = B ⊗ B (apilamiento: dos pulsos oscuros independientes que el buscador
//   no separa).
//   Variantes E-G, término de fondo de centelleo: P(N) = (1 − b)·[(1 − q)·B + q·B2] + b·Σ_n (1 − r) r^{n−1} BT(N | n; p), con
//   BT(N | n; p) = n/N · e^{−pN} (pN)^{N−n} / (N−n)! (n primarios con distribución geométrica, cada uno con la misma cascada de Borel).
// Variantes: A q = 0; B q fijo al valor pasado por canal (R·τ, de la tasa medida y la longitud del pulso; ver xtalk_prep.py); C q libre en
// [0, 0.6]; D = B con off1 libre (pico de 1 p.e. con desplazamiento propio); E = off1 libre + fondo (b, r libres), q = 0; F = E con q libre;
// G = E con q fijo (R·τ); H = D ajustado solo en [0.55 G, 2.4 G] (picos de 1 y 2 p.e.: p del cociente 2/1 con el apilamiento R·τ restado).
// Ajuste: −ln L de Poisson por bin con Minuit2 (Migrad + Hesse); χ² de Baker–Cousins = 2 Σ [μ − n + n ln(n/μ)].
// Salida: <out>.csv (ch, variant, n_pulses, p, p_err, Gf_mV, Gf_err_mV, off_mV, off1_mV, s0_mV, s1_mV, q, q_err, b, r, chi2_BC, ndf, status).
// Uso: root -l -b -q 'analysis/xtalk_2026_11_15/fit_borel.C+("<prefijo>_dark3.bin", "<máscara.bin>", "<G 32 valores>", "<q 32 valores>",
//        "<out>")'
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {
#pragma pack(push, 1)
struct Rec {
  uint32_t ev;
  uint8_t ch, region;
  uint16_t k;
  float a, q;
};
#pragma pack(pop)
constexpr int kNmax = 16;
std::vector<double> split(const char* s) {
  std::vector<double> v;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ',')) v.push_back(std::stod(t));
  return v;
}
double lfact(int n) { return std::lgamma(n + 1.0); }
double bt(int N, int n, double p) { return N < n ? 0.0 : static_cast<double>(n) / N * std::exp(-p * N + (N - n) * std::log(p * N) - lfact(N - n)); }
void probs(double p, double q, double b, double r, double* P) {
  for (int N = 1; N <= kNmax; ++N) {
    double bg = 0;
    if (b > 0)
      for (int n = 1; n <= N; ++n) bg += (1 - r) * std::pow(r, n - 1) * bt(N, n, p);
    P[N] = (1 - b) * ((1 - q) * bt(N, 1, p) + q * bt(N, 2, p)) + b * bg;
  }
}
struct Fit {
  const std::vector<double>* edges;
  const std::vector<double>* cnt;
  // x: 0 A, 1 p, 2 G, 3 off, 4 s0, 5 s1, 6 q, 7 off1, 8 b, 9 r
  void mu(const double* x, std::vector<double>& m) const {
    double P[kNmax + 1];
    probs(x[1], x[6], x[8], x[9], P);
    const auto& e = *edges;
    m.assign(e.size() - 1, 0.0);
    for (int N = 1; N <= kNmax; ++N) {
      const double c = N * x[2] + x[3] + (N == 1 ? x[7] : 0.0);
      const double s = std::sqrt(x[4] * x[4] + N * x[5] * x[5]) * std::sqrt(2.0);
      for (std::size_t i = 0; i + 1 < e.size(); ++i) m[i] += x[0] * P[N] * 0.5 * (std::erf((e[i + 1] - c) / s) - std::erf((e[i] - c) / s));
    }
  }
  double operator()(const double* x) const {
    std::vector<double> m;
    mu(x, m);
    double l = 0;
    for (std::size_t i = 0; i < m.size(); ++i) {
      const double u = std::max(m[i], 1e-12);
      l += u - (*cnt)[i] * std::log(u);
    }
    return l;
  }
};
}  // namespace

void fit_borel(const char* bin, const char* mask, const char* gains, const char* qfix, const char* out) {
  const auto G = split(gains), Q = split(qfix);
  std::vector<uint8_t> mk;
  {
    FILE* f = fopen(mask, "rb");
    uint8_t b;
    while (fread(&b, 1, 1, f) == 1) mk.push_back(b);
    fclose(f);
  }
  std::vector<std::vector<float>> a(32);
  {
    FILE* f = fopen(bin, "rb");
    Rec r;
    while (fread(&r, sizeof r, 1, f) == 1)
      if (r.region == 0 && r.ev < mk.size() && mk[r.ev]) a[r.ch].push_back(r.a);
    fclose(f);
  }
  FILE* fo = fopen((std::string(out) + ".csv").c_str(), "w");
  fprintf(fo, "ch,variant,n_pulses,p,p_err,Gf_mV,Gf_err_mV,off_mV,off1_mV,s0_mV,s1_mV,q,q_err,b,r,chi2_BC,ndf,status\n");
  for (int c = 0; c < 32; ++c) {
    const double g = G[c], lo = 0.55 * g, hi = 6.6 * g;
    const int nb = static_cast<int>(std::floor((hi - lo) / 0.2));
    std::vector<double> edges(nb + 1), cnt(nb, 0.0);
    for (int i = 0; i <= nb; ++i) edges[i] = lo + 0.2 * i;
    double tot = 0;
    for (float v : a[c])
      if (v >= lo && v < edges[nb]) cnt[static_cast<std::size_t>((v - lo) / 0.2)] += 1, tot += 1;
    if (tot < 500) {
      for (const char* vn : {"A", "B", "C", "D", "E", "F", "G", "H"}) fprintf(fo, "%d,%s,%.0f,,,,,,,,,,,,,,,insufficient\n", c, vn, tot);
      continue;
    }
    // H: solo los picos de 1 y 2 p.e., [0.55 G, 2.4 G]
    const int nbH = static_cast<int>(std::floor((2.4 * g - lo) / 0.2));
    std::vector<double> edgesH(edges.begin(), edges.begin() + nbH + 1), cntH(cnt.begin(), cnt.begin() + nbH);
    Fit fit{&edges, &cnt}, fitH{&edgesH, &cntH};
    ROOT::Math::Functor fnAll(fit, 10), fnH(fitH, 10);
    for (const char* vn : {"A", "B", "C", "D", "E", "F", "G", "H"}) {
      const std::string v(vn);
      const bool bg = v == "E" || v == "F" || v == "G";
      std::unique_ptr<ROOT::Math::Minimizer> m(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
      m->SetMaxFunctionCalls(200000);
      m->SetTolerance(0.01);
      m->SetPrintLevel(0);
      m->SetFunction(v == "H" ? fnH : fnAll);
      m->SetLimitedVariable(0, "A", tot * 1.2, tot * 0.01, 0.1 * tot, 10 * tot);
      m->SetLimitedVariable(1, "p", 0.3, 0.01, 0.001, 0.95);
      m->SetLimitedVariable(2, "G", g, 0.01 * g, 0.7 * g, 1.4 * g);
      m->SetLimitedVariable(3, "off", -1.0, 0.05, -5.0, 3.0);
      m->SetLimitedVariable(4, "s0", 0.4, 0.02, 0.02, 3.0);
      m->SetLimitedVariable(5, "s1", 0.4, 0.02, 0.0, 3.0);
      if (v == "A" || v == "E")
        m->SetFixedVariable(6, "q", 0.0);
      else if (v == "C" || v == "F")
        m->SetLimitedVariable(6, "q", 0.5 * Q[c], 0.01, 0.0, 0.6);
      else
        m->SetFixedVariable(6, "q", Q[c]);
      if (v == "D" || bg)  // H: dos picos, G y off ya fijan las dos posiciones
        m->SetLimitedVariable(7, "off1", 0.0, 0.05, -3.0, 3.0);
      else
        m->SetFixedVariable(7, "off1", 0.0);
      if (bg) {
        m->SetLimitedVariable(8, "b", 0.1, 0.01, 0.0, 0.95);
        m->SetLimitedVariable(9, "r", 0.7, 0.01, 0.0, 0.97);
      } else {
        m->SetFixedVariable(8, "b", 0.0);
        m->SetFixedVariable(9, "r", 0.0);
      }
      m->Minimize();
      m->Hesse();
      const double* x = m->X();
      const double* ex = m->Errors();
      std::vector<double> mu;
      const Fit& fu = v == "H" ? fitH : fit;
      const std::vector<double>& cu = v == "H" ? cntH : cnt;
      fu.mu(x, mu);
      double chi2 = 0;
      for (std::size_t i = 0; i < mu.size(); ++i) chi2 += 2 * (mu[i] - cu[i] + (cu[i] > 0 ? cu[i] * std::log(cu[i] / std::max(mu[i], 1e-12)) : 0.0));
      const int nfree = static_cast<int>(m->NFree());
      fprintf(fo, "%d,%s,%.0f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%d,%d\n", c, vn, tot, x[1], ex[1], x[2], ex[2], x[3], x[7], x[4],
              x[5], x[6], (v == "C" || v == "F") ? ex[6] : 0.0, x[8], x[9], chi2, static_cast<int>(mu.size()) - nfree, m->Status());
    }
  }
  fclose(fo);
}
