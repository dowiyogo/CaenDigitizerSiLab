// ttt.h
// Desenrollado del TriggerTimeTag del VX1740, compartido por daqrun y el análisis. C++17
// header-only, sin ROOT (se usa con -I daq/core/include).
// Medido en DAQ2-6 (VX1740 serie 21): tick 8 ns, contador de 31 bit (vuelta cada 17.180 s); el
// bit 31 es una marca que queda en 1 tras la primera vuelta y NO forma parte del valor.
// Algoritmo: con la hora del PC de la lectura (t_read) y el origen del run (UTC antes del start),
// se elige el número de vueltas k que deja origen + ttt más cerca de t_read; nunca retrocede
// respecto del evento anterior.
// recompute(): corrección de los runs 1 y 2 de commissioning_2026_10, escritos antes de la
// máscara de 31 bit (ttt_ns y t_utc_ns guardados son incorrectos; ttt_raw y t_read_utc_ns no).
// Funciones puras: seguras entre hilos.
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

namespace daq::ttt {

struct Params {
  double tickNs = 8.0;
  int wrapBits = 31;
  std::uint32_t mask() const { return wrapBits >= 32 ? 0xFFFFFFFFu : static_cast<std::uint32_t>((1ull << wrapBits) - 1); }
  std::int64_t periodNs() const {
    return static_cast<std::int64_t>(std::llround(std::ldexp(1.0, wrapBits) * tickNs));
  }
};

// Un evento: devuelve ttt_ns; `last` es el ttt_ns del evento anterior (−1 si no hay) y se actualiza.
inline std::int64_t unwrap(std::uint32_t raw, std::int64_t tReadUtcNs, std::int64_t originUtcNs, const Params& p,
                           std::int64_t& last) {
  const std::int64_t period = p.periodNs();
  const auto base = static_cast<std::int64_t>(std::llround(static_cast<double>(raw & p.mask()) * p.tickNs));
  const std::int64_t elapsed = tReadUtcNs - originUtcNs;
  std::int64_t k = static_cast<std::int64_t>(std::floor(static_cast<double>(elapsed - base) / static_cast<double>(period) + 0.5));
  if (k < 0) k = 0;
  std::int64_t t = base + k * period;
  if (last >= 0 && t < last) t += ((last - t) / period + 1) * period;
  last = t;
  return t;
}

struct Recomputed {
  std::vector<std::int64_t> tttNs;   // ns desde el inicio de la adquisición
  std::vector<std::int64_t> tUtcNs;  // origen + ttt_ns
};

// Recalcula una secuencia de eventos en orden (ttt_raw y t_read_utc_ns del árbol `events`;
// origen = run_info_json.ttt_origin.utc_before_start_ns).
inline Recomputed recompute(const std::vector<std::uint32_t>& raw, const std::vector<std::int64_t>& tReadUtcNs,
                            std::int64_t originUtcNs, const Params& p = Params()) {
  Recomputed r;
  r.tttNs.reserve(raw.size());
  r.tUtcNs.reserve(raw.size());
  std::int64_t last = -1;
  for (std::size_t i = 0; i < raw.size() && i < tReadUtcNs.size(); ++i) {
    const std::int64_t t = unwrap(raw[i], tReadUtcNs[i], originUtcNs, p, last);
    r.tttNs.push_back(t);
    r.tUtcNs.push_back(originUtcNs + t);
  }
  return r;
}

}  // namespace daq::ttt
