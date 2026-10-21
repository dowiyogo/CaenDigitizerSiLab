// gain_fit.h
// Ganancia G (mV por p.e.) de un canal de hod2018 desde el espectro de amplitud de sus pulsos fuera
// de tiempo (diseño 7.3; handoff 8.11), por separación de picos. C++17 header-only, sin ROOT: lo usa
// el monitor de daqrun (una vez por hora) y lo puede usar el análisis (-I daq/core/include).
//  1. espectro suavizado (media de 3 bins) menos su tendencia (media móvil de ancho = semilla);
//  2. máximos locales (≥ sus vecinos a ±1 mV) entre lo_mV y hi_mV, separados al menos
//     min_separation_frac × semilla (0.3: con una semilla 2G se encuentran igual los picos a G, y la
//     comparación con la semilla lo marca ambiguo); se queda el más alto;
//  3. posición de cada pico = centroide del residuo positivo en ±1 mV;
//  4. G = pendiente de la recta posición contra índice consecutivo (mínimos cuadrados);
//  5. prueba de espaciado: cada diferencia entre picos consecutivos dentro de spacing_tolerance de G
//     (si no: "rejected");
//  6. comparación con la semilla (hora anterior o calibración de partida): G cerca de 2× o ½× la
//     semilla (dentro de spacing_tolerance) -> "ambiguous" (no se acepta; la semilla queda); cambio
//     mayor que max_relative_change -> "changed" (se acepta y se marca para revisión); si no "ok".
//  Menos de min_peaks picos: "insufficient". No determina el número de p.e. del primer pico (p0).
// En el monitor el resultado es PROVISIONAL mientras el run esté sin calibración (7.4).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace daq::gain {

struct Options {
  double loMv = 5.0, hiMv = 26.0;  // el bulto de 1-4 mV (mínimos de ruido y apilamiento) tapa el 1 p.e.
  double minSeparationFrac = 0.3;  // × semilla: resuelve picos a la mitad de la semilla (una semilla 2G no fuerza 2G)
  double spacingTolerance = 0.2;
  double maxRelativeChange = 0.2;
  int minPeaks = 3;
  double minPeakCounts = 20;  // cuentas del bin suavizado del pico
};

struct Result {
  std::string status = "insufficient";  // ok | changed | ambiguous | rejected | insufficient
  double gainMv = 0;                   // G medida (también si se rechaza o es ambigua)
  int peaks = 0;
  std::vector<double> positionsMv;
  bool accepted() const { return status == "ok" || status == "changed"; }
};

// spectrum: cuentas por bin de binMv desde 0 mV.
template <class T>
Result fitGain(const std::vector<T>& spectrum, double binMv, double seedMv, const Options& o = Options()) {
  Result r;
  const int n = static_cast<int>(spectrum.size());
  if (n < 5 || binMv <= 0 || seedMv <= 0) return r;
  std::vector<double> s(static_cast<std::size_t>(n), 0.0);
  for (int b = 1; b + 1 < n; ++b)
    s[static_cast<std::size_t>(b)] = (static_cast<double>(spectrum[static_cast<std::size_t>(b - 1)]) +
                                      static_cast<double>(spectrum[static_cast<std::size_t>(b)]) +
                                      static_cast<double>(spectrum[static_cast<std::size_t>(b + 1)])) / 3.0;
  const int half = std::max(1, static_cast<int>(std::lround(1.0 / binMv)));  // ±1 mV
  const int lo = std::max(half, static_cast<int>(std::ceil(o.loMv / binMv)));
  const int hi = std::min(n - 1 - half, static_cast<int>(std::floor(o.hiMv / binMv)));
  // Tendencia: media móvil de ancho = semilla (promedia la oscilación de período G); los picos se
  // buscan en el residuo r = s − tendencia (DAQ2-10: con registro completo y SW trigger el espectro
  // cae tan rápido que los picos son hombros, no máximos locales).
  std::vector<double> resid(static_cast<std::size_t>(n), 0.0);
  const int hw = std::max(1, static_cast<int>(std::lround(0.5 * seedMv / binMv)));
  for (int b = 0; b < n; ++b) {
    double sum = 0;
    int k = 0;
    for (int d = -hw; d <= hw; ++d)
      if (b + d >= 0 && b + d < n) {
        sum += s[static_cast<std::size_t>(b + d)];
        ++k;
      }
    resid[static_cast<std::size_t>(b)] = s[static_cast<std::size_t>(b)] - sum / k;
  }
  std::vector<int> maxima;
  for (int b = lo; b <= hi; ++b) {
    const double v = resid[static_cast<std::size_t>(b)];
    if (s[static_cast<std::size_t>(b)] < o.minPeakCounts || v <= 0) continue;
    bool isMax = true;
    for (int d = -half; d <= half && isMax; ++d)
      if (d != 0 && resid[static_cast<std::size_t>(b + d)] > v) isMax = false;
    if (!isMax) continue;
    const double minSep = o.minSeparationFrac * seedMv / binMv;
    if (!maxima.empty() && b - maxima.back() < minSep) {
      if (v > resid[static_cast<std::size_t>(maxima.back())]) maxima.back() = b;
      continue;
    }
    maxima.push_back(b);
  }
  for (int b : maxima) {
    double w = 0, x = 0;
    for (int d = -half; d <= half; ++d) {
      const double c = std::max(0.0, resid[static_cast<std::size_t>(b + d)]);  // centroide del residuo
      w += c;
      x += c * (b + d + 0.5) * binMv;
    }
    r.positionsMv.push_back(w > 0 ? x / w : (b + 0.5) * binMv);
  }
  r.peaks = static_cast<int>(r.positionsMv.size());
  if (r.peaks < o.minPeaks) return r;
  double si = 0, sp = 0, sii = 0, sip = 0;
  for (int i = 0; i < r.peaks; ++i) {
    si += i;
    sp += r.positionsMv[static_cast<std::size_t>(i)];
    sii += static_cast<double>(i) * i;
    sip += i * r.positionsMv[static_cast<std::size_t>(i)];
  }
  const double m = r.peaks;
  r.gainMv = (m * sip - si * sp) / (m * sii - si * si);
  for (int i = 1; i < r.peaks; ++i) {
    const double d = r.positionsMv[static_cast<std::size_t>(i)] - r.positionsMv[static_cast<std::size_t>(i - 1)];
    if (std::fabs(d - r.gainMv) > o.spacingTolerance * r.gainMv) {
      r.status = "rejected";
      return r;
    }
  }
  auto near = [&](double ref) { return std::fabs(r.gainMv - ref) <= o.spacingTolerance * ref; };
  if (near(2 * seedMv) || near(0.5 * seedMv))
    r.status = "ambiguous";
  else if (std::fabs(r.gainMv - seedMv) > o.maxRelativeChange * seedMv)
    r.status = "changed";
  else
    r.status = "ok";
  return r;
}

}  // namespace daq::gain
