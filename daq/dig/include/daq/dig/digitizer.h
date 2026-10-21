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

#include <array>
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

// Auto-trigger por grupos del x740 (DAQ2-13; UM5483 sec. 1.5, 1.9, 1.12, 1.19; UM1935
// Set/GetGroupSelfTrigger, Set/GetChannelGroupMask, Set/GetGroupTriggerThreshold): cada grupo de 8
// canales pide trigger cuando al menos uno de sus canales habilitados pasa el umbral del grupo
// (0x1n80, cuentas ADC absolutas, común a los 8 canales; polaridad negativa con 0x8000 bit 6); el
// trigger global es el OR de los grupos habilitados en 0x810C[7:0] o, con nivel de mayoría m > 0
// (0x810C[26:24]), al menos m + 1 grupos dentro de la ventana 0x810C[23:20] (en relojes de 8 ns).
struct SelfTriggerConfig {
  bool enabled = false;
  std::uint32_t groupMask = 0;                 // grupos que participan (bits 0-7)
  std::vector<double> thresholdPe;             // por grupo (8 valores; los no usados se ignoran)
  std::vector<int> thresholdAdc;               // por grupo, cuentas ADC (de la configuración o calculado con la calibración)
  std::vector<int> maskedChannels;             // canales excluidos del OR de su grupo (0x1nA8)
  int majorityGroups = 1;                      // grupos requeridos dentro de la ventana (1 = OR)
  int windowClk = 1;                           // ventana de la mayoría, relojes de 8 ns (0-15)
  bool includeExternal = false;                // TRG IN en OR con el auto-trigger
  double targetRateHz = 50;                    // aviso si la tasa supera 2 × este valor
  nlohmann::json conversion = nlohmann::json(); // umbral p.e. -> ADC por grupo y umbral efectivo por canal (cabecera)
  std::uint32_t channelMask(int group) const;  // 8 bits: canales del grupo que participan
  int majorityLevel() const { return majorityGroups - 1; }  // valor del registro (m: dispara con m + 1)
  nlohmann::json toJson() const;
};

struct DigitizerConfig {
  std::string backend = "fake";  // "fake" (sintético), "replay", "caen" (VX1740 real, fase 4)
  int recordLength = 144;        // N pedido; la placa puede ajustarlo (625 -> 624): ver recordLength()
  int postTriggerRegister = 93;  // 0x8114
  std::vector<int> groupDcOffset = std::vector<int>(kGroups, 32767);
  int maxEventsBlt = 16;
  double pollSleepS = 0.005;
  double tttTickNs = 8.0;
  std::string tttTickSource = "nominal";  // de dónde sale ttt_tick_ns (p. ej. "medido en run008"); va a la cabecera
  int nBit = 12;
  double maxVpp = 2.0;
  bool softwareTrigger = false;  // modo calibración: el hilo de lectura envía SW triggers
  double swTriggerRateHz = 100;
  bool externalTrigger = true;   // TRG IN (ACQ_ONLY); false = deshabilitado
  // trigger.mode (DAQ2-13): "external" (por omisión), "software" o "self_hod2018"; fija
  // external_trigger y software_trigger. Vacío: se usan esas dos claves como antes.
  std::string triggerMode;
  SelfTriggerConfig self;
  std::string ioLevel = "NIM";   // nivel del panel frontal: NIM o TTL
  std::uint32_t groupEnableMask = 0xFF;
  // Enlace (backend caen), como el DAQ antiguo: CAEN_DGTZ_USB_A4818 = 5, linkNum 62578.
  int linkType = 5;
  int linkNum = 62578;
  int conetNode = 0;
  std::uint32_t vmeBaseAddress = 0;
  bool resetOnConfigure = true;
  // Captura de buffers crudos de ReadData (tests de decodificación): directorio y nº de bloques.
  std::string captureDir;
  int captureBlocks = 0;
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
  std::uint8_t triggerKind = 0;  // 0 externo, 1 software, 2 auto-trigger, 3 auto-trigger o externo
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
  // Longitud de registro real después de configure() (la placa puede ajustar la pedida).
  virtual int recordLength() const = 0;
};

// Monotónico del sistema en ns (CLOCK_MONOTONIC), para t_read_mono_ns.
std::int64_t monotonicNowNs();

}  // namespace daq::dig
