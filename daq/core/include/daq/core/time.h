// time.h
// Utilidades de tiempo. Unidades: UTC en nanosegundos desde la época Unix (int64, de
// std::chrono::system_clock); intervalos con std::chrono::steady_clock. Nunca hora local.
// Todas las funciones son puras o leen el reloj del sistema; seguras entre hilos.
#pragma once

#include <chrono>
#include <cstdint>
#include <string>

namespace daq {

using UtcNs = std::int64_t;
using SteadyTime = std::chrono::steady_clock::time_point;

UtcNs utcNowNs();
SteadyTime steadyNow();

// "YYYYMMDDThhmmssZ" (para nombres de archivo); trunca a segundos.
std::string formatUtcCompact(UtcNs t);
// "YYYY-MM-DDThh:mm:ss.nnnnnnnnnZ".
std::string formatUtcIso(UtcNs t);

inline double toSeconds(std::chrono::nanoseconds d) {
  return std::chrono::duration<double>(d).count();
}
inline double toMilliseconds(std::chrono::nanoseconds d) {
  return std::chrono::duration<double, std::milli>(d).count();
}

}  // namespace daq
