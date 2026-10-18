// clock.h
// Reloj inyectable. Todo el código que mide o espera tiempo lo hace a través de Clock, para
// poder probar con tiempo acelerado.
// - SystemClock: system_clock para UTC, steady_clock para intervalos y esperas.
// - ManualClock: tiempo controlado por el test. sleep_until() avanza el reloj hasta el instante
//   pedido (si está en el futuro) y vuelve de inmediato; advance() lo mueve a mano. UTC y
//   steady avanzan juntos. Seguro entre hilos (mutex interno).
// Unidades: UTC en ns (UtcNs), instantes steady como std::chrono::steady_clock::time_point.
#pragma once

#include <chrono>
#include <mutex>

#include "daq/core/jthread.h"
#include "daq/core/time.h"

namespace daq {

class Clock {
 public:
  virtual ~Clock() = default;
  virtual UtcNs now_utc_ns() const = 0;
  virtual SteadyTime now_steady() const = 0;
  // Espera hasta t o hasta que se pida parada. Devuelve false si hubo parada.
  virtual bool sleep_until(SteadyTime t, const StopToken& token = {}) = 0;
  bool sleep_for(std::chrono::nanoseconds d, const StopToken& token = {}) {
    return sleep_until(now_steady() + d, token);
  }
};

class SystemClock final : public Clock {
 public:
  UtcNs now_utc_ns() const override { return utcNowNs(); }
  SteadyTime now_steady() const override { return steadyNow(); }
  bool sleep_until(SteadyTime t, const StopToken& token = {}) override;
};

class ManualClock final : public Clock {
 public:
  // Por defecto arranca en 2026-10-18T00:00:00Z.
  explicit ManualClock(UtcNs startUtc = 1792281600LL * 1000000000LL);
  UtcNs now_utc_ns() const override;
  SteadyTime now_steady() const override;
  bool sleep_until(SteadyTime t, const StopToken& token = {}) override;
  void advance(std::chrono::nanoseconds d);

 private:
  mutable std::mutex m_;
  UtcNs utc_;
  SteadyTime steady_;
};

}  // namespace daq
