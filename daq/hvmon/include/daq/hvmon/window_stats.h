// window_stats.h
// Estadística de una magnitud en una ventana deslizante de tiempo (p. ej. VMon e IMon de los últimos
// 10 min para la página /dash/ de hvmon). Guarda las muestras de la ventana y recalcula media y
// desviación desde ellas (sin sumas acumuladas que arrastren redondeo durante días).
// Seguridad entre hilos: pertenece al hilo de sondeo; el resumen se publica por PollerState.
#pragma once

#include <algorithm>
#include <cmath>
#include <deque>
#include <utility>

#include "daq/core/time.h"

namespace daq::hvmon {

struct WindowSummary {
  std::size_t n = 0;
  double mean = 0, std = 0, min = 0, max = 0;  // std: desviación poblacional
};

class WindowStats {
 public:
  explicit WindowStats(double windowS = 600.0) : windowNs_(static_cast<UtcNs>(windowS * 1e9)) {}
  void add(UtcNs t, double v) {
    samples_.emplace_back(t, v);
    prune(t);
  }
  // Resumen de las muestras con t > now − ventana.
  WindowSummary summary(UtcNs now) const {
    WindowSummary s;
    double sum = 0;
    for (const auto& [t, v] : samples_)
      if (t > now - windowNs_) {
        if (s.n == 0) s.min = s.max = v;
        ++s.n;
        sum += v;
        s.min = std::min(s.min, v);
        s.max = std::max(s.max, v);
      }
    if (s.n == 0) return s;
    s.mean = sum / static_cast<double>(s.n);
    double ss = 0;
    for (const auto& [t, v] : samples_)
      if (t > now - windowNs_) ss += (v - s.mean) * (v - s.mean);
    s.std = std::sqrt(ss / static_cast<double>(s.n));
    return s;
  }
  double windowS() const { return static_cast<double>(windowNs_) / 1e9; }

 private:
  void prune(UtcNs now) {
    while (!samples_.empty() && samples_.front().first <= now - windowNs_) samples_.pop_front();
  }
  UtcNs windowNs_;
  std::deque<std::pair<UtcNs, double>> samples_;
};

}  // namespace daq::hvmon
