// change_detector.h
// Detección de cambios entre dos lecturas consecutivas del mismo tipo (rápida o completa).
// Origen de cada cambio:
// - "trip": cambio de Pw o Status de un canal cuyo Status nuevo tiene un bit de trip
//   (internal_trip o external_trip, configurable);
// - "hvmon_command": cambio de Pw o Status de un canal con un comando de hvmon en curso
//   (desde la escritura hasta command_window_extra_s después de verificar), con su id;
// - "external": cualquier otro (GECO u otro cliente).
// La primera lectura tras arrancar no produce cambios. Tras una pérdida de conexión se compara con
// la última lectura buena (duringDisconnect = true): si nada cambió, no hay cambios espurios.
// Seguridad entre hilos: lo usa solo el hilo de sondeo.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "daq/core/time.h"
#include "daq/hv/hv_types.h"
#include "daq/hv/status_bits.h"
#include "daq/hvmon/state.h"

namespace daq::hvmon {

struct PendingCommand {
  std::uint64_t id = 0;
  int ch = 0;
  SteadyTime tStart{};
  SteadyTime tEnd{};
};

class ChangeDetector {
 public:
  ChangeDetector(hv::StatusBitMap bits, std::vector<std::string> tripFlags);

  std::vector<ChangeRecord> compare(const hv::ChannelReading& prev, const hv::ChannelReading& cur,
                                    const std::vector<std::string>& tracked, const std::string& readKind,
                                    bool duringDisconnect, SteadyTime now, std::uint64_t& nextId) const;

  void addPending(const PendingCommand& p);
  void setPendingEnd(std::uint64_t id, SteadyTime tEnd);
  void expire(SteadyTime now);
  const std::vector<PendingCommand>& pending() const { return pending_; }
  bool hasTrip(std::uint32_t status) const;

 private:
  hv::StatusBitMap bits_;
  std::vector<std::string> tripFlags_;
  std::vector<PendingCommand> pending_;
};

}  // namespace daq::hvmon
