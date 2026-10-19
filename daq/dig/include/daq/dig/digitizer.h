// digitizer.h
// Abstracción del digitizador (fase 3). DigitizerBackend es la interfaz que usa el bucle de
// lectura; en esta fase solo existen los fakes (sintético y replay de v8); la fase 4 agrega
// CaenDigitizerBackend (CAENDigitizer) sin tocar el bucle.
// Unidades: muestras ADC crudas (uint16, 12 bit en el VX1740), tiempos UTC en ns (CLOCK_REALTIME)
// y monotónicos en ns (CLOCK_MONOTONIC). TTT crudo en ticks de ttt_tick_ns.
// Seguridad entre hilos: una instancia de backend pertenece al hilo de lectura (único llamador,
// como el handle de CAENDigitizer). Los EventBlock se publican como shared_ptr<const ...>: después
// de la lectura son inmutables y los comparten el escritor y el monitor sin copia.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/core/time.h"

namespace daq::dig {

inline constexpr int kChannels = 64;
inline constexpr int kGroups = 8;

struct DigitizerConfig {
  std::string backend = "fake";  // "fake" (sintético), "replay"; "caen" en la fase 4
  int recordLength = 144;        // N, múltiplo de 3
  int postTriggerRegister = 93;  // 0x8114
  std::vector<int> groupDcOffset = std::vector<int>(kGroups, 32767);
  int maxEventsBlt = 16;
  double pollSleepS = 0.005;
  double tttTickNs = 8.0;
  int nBit = 12;
  double maxVpp = 2.0;
  bool softwareTrigger = false;  // modo calibración: el hilo de lectura envía SW triggers
  double swTriggerRateHz = 100;
  std::string scenarioFile;      // fake: escenario JSON; replay: configuración del replay
  nlohmann::json scenario;       // contenido (resuelto al cargar)

  static Result<DigitizerConfig> fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
};

// Conversión a mV del código antiguo (adc2mV de CaenDigitizerSiLab.h):
//   mV = (adc / (2^NBit − 1) − 1 + dcOffset / 0xFFFF) · MaxVpp · 1000
// y su inversa. Con polarización 0 el DAQ antiguo usó dcOffset = 0x7FFF.
struct AdcConversion {
  int nBit = 12;
  double maxVpp = 2.0;
  int dcOffset = 0x7FFF;
  double toMv(double adc) const {
    return (adc / ((1 << nBit) - 1.0) - 1.0 + static_cast<double>(dcOffset) / 65535.0) * maxVpp * 1000.0;
  }
  double toAdc(double mv) const {
    return (mv / (maxVpp * 1000.0) + 1.0 - static_cast<double>(dcOffset) / 65535.0) * ((1 << nBit) - 1.0);
  }
  double lsbMv() const { return maxVpp * 1000.0 / ((1 << nBit) - 1.0); }
  nlohmann::json toJson() const {
    return {{"formula", "mV = (adc/(2^NBit-1) - 1 + dc_offset/65535) * MaxVpp * 1000"},
            {"NBit", nBit}, {"MaxVpp_V", maxVpp}, {"dc_offset", dcOffset}, {"lsb_mV", lsbMv()}};
  }
};

struct Event {
  std::uint32_t counterHw = 0;   // contador de eventos de la placa
  std::uint32_t tttRaw = 0;      // TriggerTimeTag tal como llega
  std::uint8_t triggerKind = 0;  // 0 externo, 1 software
  std::vector<std::uint16_t> adc;  // [kChannels][N], canal mayor
};

struct EventBlock {
  std::vector<Event> events;
  UtcNs tReadUtcNs = 0;
  std::int64_t tReadMonoNs = 0;
  bool boardFull = false;  // la memoria de la placa se llenó desde la lectura anterior
  int recordLength = 0;
};
using EventBlockPtr = std::shared_ptr<const EventBlock>;

struct AcqOrigin {
  UtcNs utcBeforeStart = 0;  // UTC antes de SWStartAcquisition
  UtcNs utcAfterStart = 0;   // UTC después
};

class DigitizerBackend {
 public:
  virtual ~DigitizerBackend() = default;
  virtual Status configure(const DigitizerConfig& cfg) = 0;
  virtual Result<AcqOrigin> start() = 0;
  virtual Status stop() = 0;
  // Agrega a out hasta maxEvents eventos disponibles (0 si no hay nada); llena tRead*.
  virtual Result<std::size_t> readBlock(EventBlock& out, std::size_t maxEvents) = 0;
  virtual Status sendSoftwareTrigger() { return Status::failure(-1, "SW trigger no soportado"); }
  // Modelo, serie, firmware, versiones de biblioteca ("unknown" en los fakes).
  virtual nlohmann::json info() const = 0;
  virtual bool simulated() const = 0;
  // true cuando una fuente finita (replay) ya entregó todo.
  virtual bool exhausted() const { return false; }
};

// Monotónico del sistema en ns (CLOCK_MONOTONIC), para t_read_mono_ns.
std::int64_t monotonicNowNs();

}  // namespace daq::dig
