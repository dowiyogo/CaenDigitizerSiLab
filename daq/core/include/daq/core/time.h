// time.h
// Utilidades de tiempo. Unidades: UTC en nanosegundos desde la época Unix (int64, de
// std::chrono::system_clock); intervalos con std::chrono::steady_clock. Nunca hora local.
// Todas las funciones son puras o leen el reloj del sistema; seguras entre hilos.
#pragma once

#include <chrono>
#include <cstdint>
#include <string>

#include "daq/core/result.h"

namespace daq {

using UtcNs = std::int64_t;
using SteadyTime = std::chrono::steady_clock::time_point;

UtcNs utcNowNs();
SteadyTime steadyNow();

// "YYYYMMDDThhmmssZ" (para nombres de archivo); trunca a segundos.
std::string formatUtcCompact(UtcNs t);
// "YYYY-MM-DDThh:mm:ss.nnnnnnnnnZ".
std::string formatUtcIso(UtcNs t);

// Hora local solo para textos para personas (correos), nunca para archivos de datos.
// setProcessTimeZone() fija TZ del proceso (p. ej. "America/Santiago") y llama tzset(); debe
// llamarse al arrancar, antes de crear hilos. Falla si la zona no existe en /usr/share/zoneinfo.
Status setProcessTimeZone(const std::string& tz);
// "YYYY-MM-DD hh:mm:ss -03" en la zona del proceso (localtime_r: seguro entre hilos).
std::string formatLocalShort(UtcNs t);

inline double toSeconds(std::chrono::nanoseconds d) {
  return std::chrono::duration<double>(d).count();
}
inline double toMilliseconds(std::chrono::nanoseconds d) {
  return std::chrono::duration<double, std::milli>(d).count();
}

}  // namespace daq
