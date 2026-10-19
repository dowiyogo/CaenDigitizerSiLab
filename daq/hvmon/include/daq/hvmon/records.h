// records.h
// Registros que van al archivo de cambios (hv_changes_YYYYMMDD.jsonl) y al anillo del estado:
// cambios de parámetros (type=change) y fases de los comandos de escritura (type=command).
// Separado de state.h para que retry.h y notifier.h no dependan de la configuración.
// Tipos de valor; serialización en state.cpp.
#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/time.h"

namespace daq::hvmon {

struct ChangeRecord {
  std::uint64_t id = 0;
  UtcNs tUtcNs = 0;          // t_resp de la lectura que mostró el valor nuevo
  UtcNs tPrevReadUtcNs = 0;  // t_req de la lectura anterior (el cambio ocurrió entre ambas)
  int slot = 0;
  int ch = 0;
  std::string param;
  std::string unit;          // unidad SI de old/new ("" para enteros)
  nlohmann::json oldValue;
  nlohmann::json newValue;
  std::string origin;        // "trip", "hvmon_command", "hvmon_retry" o "external"
  std::optional<std::uint64_t> commandId;
  bool duringDisconnect = false;  // la lectura anterior es de antes de una pérdida de conexión
  std::string readKind;      // "fast", "full" o "retry"
  nlohmann::json detail;     // hvmon_retry: {"k":..., "reason":...}; null en los demás
  nlohmann::json toJson() const;
};

// Registro de un comando de escritura (va al mismo JSONL que los cambios, con "type":"command").
struct CommandRecord {
  std::uint64_t id = 0;
  UtcNs tUtcNs = 0;
  int slot = 0;
  int ch = 0;
  std::string param = "Pw";
  int value = 0;
  std::string phase;      // request, rejected, aborted, noop, written, write_failed, verified, verify_failed
  std::string reason;     // texto de --authorize
  std::string requester;  // uid/pid del cliente o "hvmon"
  std::string detail;
  nlohmann::json toJson() const;
};

}  // namespace daq::hvmon
