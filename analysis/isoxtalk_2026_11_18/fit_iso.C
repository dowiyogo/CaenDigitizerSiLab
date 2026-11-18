// fit_iso.C — seguimiento 2026-11-18: ajuste de Borel del espectro de amplitud de los pulsos oscuros aislados de hod2018, sin término de apilamiento.
// Entrada: archivo de pares empaquetados {uint8 ch; float a_mV} (5 bytes) que escribe iso_xtalk.py para una campaña, región y ventana de aislamiento.
// Modelo por canal, histograma de a en bins de 0.2 mV en [0.55 G, 6.6 G] (picos de 1 a 6 p.e.; G = la del análisis, solo para la ventana y la semilla):
//   dN/da = A · Σ_{N=1}^{16} B(N; p) · Gauss(a; N·Gf + off (+ off1 si N = 1), σ_N² = s0² + N·s1²), integrada en el bin (erf),
//   B(N; p) = e^{−pN} (pN)^{N−1} / N!  (respuesta de 1 p.e. = gaussianas con la separación Gf y el desplazamiento off ajustados).
// Variantes: A = posiciones lineales (off1 = 0); D = off1 libre (el pico de 1 p.e. con desplazamiento propio).
// −ln L de Poisson por bin con Minuit2 (Migrad + Hesse); χ² de Baker–Cousins = 2 Σ [μ − n + n ln(n/μ)].
// Salida: <out>.csv (ch, variant, n_fit, p, p_err, Gf_mV, off_mV, off1_mV, s0_mV, s1_mV, chi2_BC, ndf, status).
// Uso: root -l -b -q 'analysis/isoxtalk_2026_11_18/fit_iso.C+("<pares.bin>", "<G 32 valores>", "<out>")'
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
constexpr int kNmax = 16;
struct Fit {
  const std::vector<double>* e;
  const std::vector<double>* n;
  void mu(const double* x, std::vector<double>& m) const {
    m.assign(e->size() - 1, 0.0);
    for (int N = 1; N <= kNmax; ++N) {
      const double b = std::exp(-x[1] * N + (N - 1) * std::log(x[1] * N) - std::lgamma(N + 1.0));
      const double c = N * x[2] + x[3] + (N == 1 ? x[6] : 0.0);
      const double s = std::sqrt(x[4] * x[4] + N * x[5] * x[5]) * std::sqrt(2.0);
      for (std::size_t i = 0; i + 1 < e->size(); ++i) m[i] += x[0] * b * 0.5 * (std::erf(((*e)[i + 1] - c) / s) - std::erf(((*e)[i] - c) / s));
    }
  }
  double operator()(const double* x) const {
    std::vector<double> m;
    mu(x, m);
    double l = 0;
    for (std::size_t i = 0; i < m.size(); ++i) {
      const double u = std::max(m[i], 1e-12);
      l += u - (*n)[i] * std::log(u);
    }
    return l;
  }
};
}  // namespace

void fit_iso(const char* bin, const char* gains, const char* out) {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string t;
    while (std::getline(ss, t, ',')) G.push_back(std::stod(t));
  }
  std::vector<std::vector<float>> a(32);
  {
    FILE* f = fopen(bin, "rb");
    uint8_t ch;
    float v;
    while (fread(&ch, 1, 1, f) == 1 && fread(&v, 4, 1, f) == 1) a[ch].push_back(v);
    fclose(f);
  }
  FILE* fo = fopen((std::string(out) + ".csv").c_str(), "w");
  fprintf(fo, "ch,variant,n_fit,p,p_err,Gf_mV,off_mV,off1_mV,s0_mV,s1_mV,chi2_BC,ndf,status\n");
  for (int c = 0; c < 32; ++c) {
    const double g = G[c], lo = 0.55 * g;
    const int nb = static_cast<int>(std::floor((6.6 * g - lo) / 0.2));
    std::vector<double> e(nb + 1), n(nb, 0.0);
    for (int i = 0; i <= nb; ++i) e[i] = lo + 0.2 * i;
    double tot = 0;
    for (float v : a[c])
      if (v >= lo && v < e[nb]) n[static_cast<std::size_t>((v - lo) / 0.2)] += 1, tot += 1;
    for (const char* vn : {"A", "D"}) {
      if (tot < 500) {
        fprintf(fo, "%d,%s,%.0f,,,,,,,,,,insufficient\n", c, vn, tot);
        continue;
      }
      Fit fit{&e, &n};
      ROOT::Math::Functor fn(fit, 7);
      std::unique_ptr<ROOT::Math::Minimizer> m(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
      m->SetMaxFunctionCalls(200000);
      m->SetTolerance(0.01);
      m->SetPrintLevel(0);
      m->SetFunction(fn);
      m->SetLimitedVariable(0, "A", tot * 1.2, tot * 0.01, 0.1 * tot, 10 * tot);
      m->SetLimitedVariable(1, "p", 0.3, 0.01, 0.001, 0.95);
      m->SetLimitedVariable(2, "G", g, 0.01 * g, 0.7 * g, 1.4 * g);
      m->SetLimitedVariable(3, "off", -1.0, 0.05, -5.0, 3.0);
      m->SetLimitedVariable(4, "s0", 0.4, 0.02, 0.02, 3.0);
      m->SetLimitedVariable(5, "s1", 0.4, 0.02, 0.0, 3.0);
      if (std::string(vn) == "D")
        m->SetLimitedVariable(6, "off1", 0.0, 0.05, -3.0, 3.0);
      else
        m->SetFixedVariable(6, "off1", 0.0);
      m->Minimize();
      m->Hesse();
      const double* x = m->X();
      const double* ex = m->Errors();
      std::vector<double> mu;
      fit.mu(x, mu);
      double chi2 = 0;
      for (int i = 0; i < nb; ++i) chi2 += 2 * (mu[i] - n[i] + (n[i] > 0 ? n[i] * std::log(n[i] / std::max(mu[i], 1e-12)) : 0.0));
      fprintf(fo, "%d,%s,%.0f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%d,%d\n", c, vn, tot, x[1], ex[1], x[2], x[3], x[6], x[4], x[5], chi2,
              nb - static_cast<int>(m->NFree()), m->Status());
    }
  }
  fclose(fo);
}
