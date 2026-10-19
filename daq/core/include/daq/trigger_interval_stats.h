// trigger_interval_stats.h
// Estadística de intervalos entre triggers (diseño 7.5). C++17 header-only, sin ROOT ni
// dependencias del build de CMake: analysis/ lo usa con -I../daq/core/include.
// Entrada: tiempos de trigger en ns (crecientes). Salida:
//  - tasa por máxima verosimilitud de un proceso de Poisson (ajuste exponencial de los Δt):
//    r = n / ΣΔt, error r/√n;
//  - prueba de exponencialidad: Kolmogorov-Smirnov (D y p asintótico) y χ² binado (p con la
//    gamma incompleta regularizada);
//  - huecos: Δt > T con exp(−r·T) = gap_probability;
//  - ráfagas: estadística de barrido; ventana W con μ = r·W = burst_window_expected; una ventana
//    que empieza en un trigger y contiene n triggers con P(N ≥ n | μ) < burst_probability es una
//    ráfaga (las ventanas marcadas no se solapan).
// Resolución declarada (resolution_ns): para datos con marcas cuantizadas (v7-v10: timestamp de
// la lectura, eventos leídos juntos comparten marca). Los Δt < resolución se cuentan como "no
// resueltos" y se excluyen de KS/χ²; la tasa usa todos.
// Funciones puras, sin estado: seguras entre hilos.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace daq::tis {

struct Options {
  double gapProbability = 1e-6;
  double burstProbability = 1e-6;
  double burstWindowExpected = 5.0;  // μ = r·W de la ventana de barrido
  int chi2Bins = 20;                 // bins equiprobables
  double resolutionNs = 0.0;
};

struct Summary {
  std::size_t nEvents = 0;
  std::size_t nIntervals = 0;
  std::size_t nUnresolved = 0;  // Δt < resolución
  double rateHz = 0, rateErrHz = 0;
  double meanDtS = 0, maxDtS = 0;
  double ksD = 0, ksP = 1;
  double chi2 = 0;
  int chi2Ndf = 0;
  double chi2P = 1;
  double gapThresholdS = 0;
  std::size_t nGaps = 0;
  double burstWindowS = 0;
  std::size_t nBursts = 0;
  std::size_t burstEvents = 0;  // triggers dentro de ventanas marcadas
  std::vector<std::size_t> gapIndices;    // índice i del intervalo (t[i+1] − t[i])
  std::vector<std::size_t> burstStarts;   // índice del primer trigger de cada ráfaga
};

namespace detail {

// ln Γ(x) (Lanczos).
inline double lnGamma(double x) { return std::lgamma(x); }

// Gamma incompleta regularizada (Numerical Recipes): serie para x < a + 1 (da P con precisión),
// fracción continua si no (da Q con precisión). gammaP + gammaQ = 1.
inline double gammaSeriesP(double a, double x) {
  double sum = 1.0 / a, del = sum, ap = a;
  for (int n = 0; n < 1000; ++n) {
    ap += 1;
    del *= x / ap;
    sum += del;
    if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
  }
  return sum * std::exp(-x + a * std::log(x) - lnGamma(a));
}
inline double gammaContinuedQ(double a, double x) {
  double b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
  for (int i = 1; i < 1000; ++i) {
    const double an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (std::fabs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (std::fabs(c) < 1e-300) c = 1e-300;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < 1e-16) break;
  }
  return std::exp(-x + a * std::log(x) - lnGamma(a)) * h;
}
inline double gammaP(double a, double x) {
  if (x <= 0) return 0.0;
  return x < a + 1 ? gammaSeriesP(a, x) : 1.0 - gammaContinuedQ(a, x);
}
inline double gammaQ(double a, double x) {
  if (x <= 0) return 1.0;
  return x < a + 1 ? 1.0 - gammaSeriesP(a, x) : gammaContinuedQ(a, x);
}

// P(N ≥ n | μ) de Poisson = P(n, μ) (gamma regularizada inferior).
inline double poissonUpperTail(std::size_t n, double mu) {
  if (n == 0) return 1.0;
  return gammaP(static_cast<double>(n), mu);
}

// p-valor asintótico de Kolmogorov para √n·D.
inline double kolmogorovP(double lambda) {
  if (lambda < 0.2) return 1.0;
  double s = 0;
  for (int k = 1; k <= 100; ++k) {
    const double term = 2.0 * ((k % 2) ? 1.0 : -1.0) * std::exp(-2.0 * k * k * lambda * lambda);
    s += term;
    if (std::fabs(term) < 1e-12) break;
  }
  return std::min(1.0, std::max(0.0, s));
}

}  // namespace detail

inline double poissonUpperTail(std::size_t n, double mu) { return detail::poissonUpperTail(n, mu); }
inline double chi2Pvalue(double chi2, int ndf) { return ndf > 0 ? detail::gammaQ(0.5 * ndf, 0.5 * chi2) : 1.0; }

// Tiempos en ns, crecientes.
inline Summary summarize(const std::vector<std::int64_t>& tNs, const Options& o = Options()) {
  Summary s;
  s.nEvents = tNs.size();
  if (tNs.size() < 2) return s;
  std::vector<double> dts;  // s
  dts.reserve(tNs.size() - 1);
  double sum = 0;
  for (std::size_t i = 1; i < tNs.size(); ++i) {
    const double dt = static_cast<double>(tNs[i] - tNs[i - 1]) * 1e-9;
    dts.push_back(dt);
    sum += dt;
    s.maxDtS = std::max(s.maxDtS, dt);
  }
  s.nIntervals = dts.size();
  s.meanDtS = sum / static_cast<double>(dts.size());
  s.rateHz = sum > 0 ? static_cast<double>(dts.size()) / sum : 0;
  s.rateErrHz = s.rateHz / std::sqrt(static_cast<double>(dts.size()));
  if (s.rateHz <= 0) return s;
  const double r = s.rateHz;
  const double res = o.resolutionNs * 1e-9;

  // KS y χ² con los Δt resueltos (distribución exponencial condicionada a Δt ≥ res).
  std::vector<double> resolved;
  for (double dt : dts) {
    if (dt < res)
      ++s.nUnresolved;
    else
      resolved.push_back(dt - res);
  }
  if (resolved.size() >= 5) {
    std::sort(resolved.begin(), resolved.end());
    const double n = static_cast<double>(resolved.size());
    double d = 0;
    for (std::size_t i = 0; i < resolved.size(); ++i) {
      const double f = 1.0 - std::exp(-r * resolved[i]);
      d = std::max(d, std::max(std::fabs(f - static_cast<double>(i) / n), std::fabs(static_cast<double>(i + 1) / n - f)));
    }
    s.ksD = d;
    const double sn = std::sqrt(n);
    s.ksP = detail::kolmogorovP((sn + 0.12 + 0.11 / sn) * d);
    const int nb = std::max(2, o.chi2Bins);
    std::vector<double> counts(static_cast<std::size_t>(nb), 0.0);
    for (double dt : resolved) {
      const double u = 1.0 - std::exp(-r * dt);  // uniforme en [0,1) si es exponencial
      const int b = std::min(nb - 1, static_cast<int>(u * nb));
      counts[static_cast<std::size_t>(b)] += 1;
    }
    const double expct = n / nb;
    for (double c : counts) s.chi2 += (c - expct) * (c - expct) / expct;
    s.chi2Ndf = nb - 2;  // un parámetro ajustado (r)
    s.chi2P = chi2Pvalue(s.chi2, s.chi2Ndf);
  }

  // Huecos.
  s.gapThresholdS = -std::log(o.gapProbability) / r;
  for (std::size_t i = 0; i < dts.size(); ++i)
    if (dts[i] > s.gapThresholdS) {
      ++s.nGaps;
      s.gapIndices.push_back(i);
    }

  // Ráfagas (barrido desde cada trigger, ventanas no solapadas).
  s.burstWindowS = o.burstWindowExpected / r;
  const auto wNs = static_cast<std::int64_t>(s.burstWindowS * 1e9);
  std::size_t j = 0;
  for (std::size_t i = 0; i < tNs.size();) {
    if (j < i) j = i;
    while (j < tNs.size() && tNs[j] - tNs[i] < wNs) ++j;
    const std::size_t n = j - i;
    if (detail::poissonUpperTail(n, o.burstWindowExpected) < o.burstProbability) {
      ++s.nBursts;
      s.burstEvents += n;
      s.burstStarts.push_back(i);
      i = j;
    } else {
      ++i;
    }
  }
  return s;
}

}  // namespace daq::tis
