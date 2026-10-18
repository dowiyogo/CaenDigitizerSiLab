// state.h
// Tipos del estado de hvmon y su composición en JSON (archivo /dev/shm/daqv2/hvmon_state.json y
// respuesta "status" del socket). Esquema documentado en docs/daq_v2_design.md, sección 4.5.
// - PollerState: lo produce el hilo de sondeo (lecturas, conexión, cambios recientes, comandos).
// - WriterHealth: lo produce el hilo escritor.
// - StateHub: reúne los Snapshot de cada productor; compose() es segura desde cualquier hilo.
// Ningún campo contiene secretos. Unidades SI; tiempos UTC ns.
#pragma once

#include <atomic>
#include <cstdint>
#include <deque>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/snapshot.h"
#include "daq/core/time.h"
#include "daq/hv/hv_types.h"
#include "daq/hvmon/config.h"

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
  std::string origin;        // "trip", "hvmon_command" o "external"
  std::optional<std::uint64_t> commandId;
  bool duringDisconnect = false;  // la lectura anterior es de antes de una pérdida de conexión
  std::string readKind;      // "fast" o "full"
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

struct ReadStats {
  std::uint64_t fastReads = 0, fullReads = 0, readErrors = 0;
  std::uint64_t connectAttempts = 0, connectFailures = 0, reconnections = 0;
  double sumFastMs = 0, maxFastMs = 0, sumFullMs = 0, maxFullMs = 0;
  nlohmann::json toJson() const;
};

struct PollerState {
  std::uint64_t readSeq = 0;
  bool connected = false;
  bool everConnected = false;
  UtcNs lastGoodReadUtcNs = 0;
  UtcNs lastErrorUtcNs = 0;
  std::string lastError;
  double nextConnectInS = 0;
  ReadStats stats;
  std::vector<hv::ParamInfo> params;                        // orden de la placa
  std::map<int, std::map<std::string, hv::ParamValue>> values;  // canal -> parámetro -> valor
  std::deque<ChangeRecord> recentChanges;
  nlohmann::json activeCommand;  // null si no hay
  nlohmann::json activeOffsets;
  nlohmann::json lastOffsets;
};

struct WriterHealth {
  bool ok = true;
  std::string lastError;
  UtcNs lastErrorUtcNs = 0;
  std::uint64_t errors = 0;
  std::uint64_t fsyncs = 0;
  std::uint64_t bytes = 0;
  std::uint64_t lines = 0;
  std::map<std::string, std::uint64_t> openFiles;  // ruta -> bytes escritos en esta sesión
  nlohmann::json toJson() const;
};

class StateHub {
 public:
  Snapshot<PollerState> poller;
  Snapshot<WriterHealth> writer;
  Snapshot<nlohmann::json> temperature{nlohmann::json()};
  std::atomic<bool> gecoRunning{false};
  std::atomic<std::uint64_t> writerDropped{0};
  UtcNs startedUtcNs = 0;

  // Estado completo. now: UTC de publicación; seq: número de publicación; final: último estado.
  nlohmann::json compose(const HvmonConfig& cfg, UtcNs now, std::uint64_t seq, bool final) const;
};

}  // namespace daq::hvmon
