// pulse_finder.h
// Búsqueda de pulsos (diseño 7.1), compartida por daqrun (muestras ADC, UShort_t) y los macros de
// análisis (muestras en mV de los archivos antiguos, float). C++17 header-only, sin ROOT.
//  - Línea base: media recortada de una región (promedio de las muestras a ±trim de la mediana).
//  - Pulso: mínimo local (muestra ≤ sus dos vecinas) con amplitud = (línea_base − x)·escala sobre
//    el umbral; dentro de min_separation muestras se queda el más grande.
//  - Umbral en mV o, si hay calibración, en p.e.: umbral_mV = p0_mV + n_pe·G_mV.
// La escala convierte unidades de muestra a mV: 1 para mV, lsbMv() para cuentas ADC (los pulsos
// del hodoscopio son negativos: la amplitud crece hacia abajo).
// Funciones puras: seguras entre hilos.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace daq::pulse {

// Paso del ADC en mV (VX1740: 12 bit, 2 Vpp) = 2000/4095 = 0.4884 mV.
inline constexpr double lsbMv(int nBit = 12, double maxVpp = 2.0) {
  return maxVpp * 1000.0 / (static_cast<double>(1 << nBit) - 1.0);
}

struct Region {
  int begin = 0;  // [begin, end)
  int end = 0;
  int size() const { return end - begin; }
};

struct Pulse {
  int sample = -1;
  double amplitudeMv = 0;
};

inline double thresholdFromPe(double nPe, double p0Mv, double gainMvPerPe) { return p0Mv + nPe * gainMvPerPe; }

template <class T>
double median(const T* x, Region r) {
  if (r.size() <= 0) return 0;
  std::vector<double> v(x + r.begin, x + r.end);
  const auto mid = v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2);
  std::nth_element(v.begin(), mid, v.end());
  double m = *mid;
  if (v.size() % 2 == 0) m = 0.5 * (m + *std::max_element(v.begin(), mid));
  return m;
}

// Media de las muestras a ±trim (unidades de x) de la mediana; mediana si ninguna cae dentro.
template <class T>
double trimmedMean(const T* x, Region r, double trim) {
  const double m = median(x, r);
  double s = 0;
  int n = 0;
  for (int i = r.begin; i < r.end; ++i)
    if (std::fabs(static_cast<double>(x[i]) - m) <= trim) {
      s += static_cast<double>(x[i]);
      ++n;
    }
  return n > 0 ? s / n : m;
}

// Pulsos en las regiones dadas de una forma de onda de n muestras. baseline en unidades de x;
// scale: unidades de x -> mV; threshold en mV. Los pulsos se agregan a out en orden de muestra.
template <class T>
void findPulses(const T* x, int n, double baseline, double scale, double thresholdMv, int minSeparation,
                const std::vector<Region>& regions, std::vector<Pulse>& out) {
  for (const Region& rg : regions) {
    const int b = std::max(rg.begin, 0), e = std::min(rg.end, n);
    bool havePending = false;
    Pulse pending;
    for (int i = b; i < e; ++i) {
      const double v = static_cast<double>(x[i]);
      const double left = i > 0 ? static_cast<double>(x[i - 1]) : v;
      const double right = i + 1 < n ? static_cast<double>(x[i + 1]) : v;
      if (!(v <= left && v <= right)) continue;
      const double amp = (baseline - v) * scale;
      if (amp <= thresholdMv) continue;
      if (havePending && i - pending.sample < minSeparation) {
        if (amp > pending.amplitudeMv) pending = Pulse{i, amp};
        continue;
      }
      if (havePending) out.push_back(pending);
      pending = Pulse{i, amp};
      havePending = true;
    }
    if (havePending) out.push_back(pending);
  }
}

// Atajo: cuenta de pulsos sobre el umbral.
template <class T>
int countPulses(const T* x, int n, double baseline, double scale, double thresholdMv, int minSeparation,
                const std::vector<Region>& regions) {
  std::vector<Pulse> p;
  findPulses(x, n, baseline, scale, thresholdMv, minSeparation, regions, p);
  return static_cast<int>(p.size());
}

}  // namespace daq::pulse
