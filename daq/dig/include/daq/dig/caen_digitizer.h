// caen_digitizer.h
// CaenDigitizerBackend (fase 4): VX1740 por CAENDigitizer 2.17.3 (enlace A4818, como el DAQ
// antiguo: CAEN_DGTZ_USB_A4818, linkNum 62578). Enlazado con las copias versionadas de /usr/lib.
//  - open(): OpenDigitizer, GetInfo (modelo, serie, firmware ROC/AMC), versiones de CAENDigitizer,
//    CAENComm y CAENVME, información del A4818 por CAENComm_Info.
//  - configure(): Reset (opcional), longitud de registro, máscara de grupos, offsets DC por grupo,
//    self-trigger de grupos deshabilitado, TRG IN (ACQ_ONLY o deshabilitado), nivel NIM/TTL, SW
//    trigger (ACQ_ONLY en modo calibración), adquisición SW_CONTROLLED, eventos por bloque, registro
//    0x8114 de post-trigger. Luego RELEE todo y guarda lo escrito y lo leído (info()).
//  - readBlock(): estado de adquisición (0x8104) para memoria llena, ReadData MBLT, GetNumEvents,
//    GetEventInfo + DecodeEvent (UINT16) por evento; contador y TTT del encabezado.
//  - Cada llamada CAEN se comprueba; un error se registra con el código y el texto (x740_raw.h).
// Seguridad entre hilos: pertenece al hilo de lectura (único llamador del handle); open/configure
// los llama main antes de crear los hilos; info() después de configure o al terminar la lectura.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "daq/core/clock.h"
#include "daq/dig/digitizer.h"

namespace daq::dig {

class CaenDigitizerBackend final : public DigitizerBackend {
 public:
  explicit CaenDigitizerBackend(Clock& clock);
  ~CaenDigitizerBackend() override;

  Status open(const DigitizerConfig& link);
  Status configure(const DigitizerConfig& cfg) override;
  Result<AcqOrigin> start() override;
  Status stop() override;
  Result<std::size_t> readBlock(EventBlock& out, std::size_t maxEvents) override;
  Status sendSoftwareTrigger() override;
  nlohmann::json info() const override;
  bool simulated() const override { return false; }
  int recordLength() const override { return recordLength_; }

  // Lectura de registros (solo lectura); nombre -> valor.
  Result<std::uint32_t> readRegister(std::uint32_t addr);
  nlohmann::json dumpRegisters(const std::vector<std::pair<std::uint32_t, std::string>>& regs);
  static std::vector<std::pair<std::uint32_t, std::string>> standardRegisters();
  // Bit de memoria llena en 0x8104 (por omisión el 4, "Event Full"; ver DAQ2-6).
  void setMemoryFullMask(std::uint32_t m) { memFullMask_ = m; }
  void close();

 private:
  Status check(int code, const char* what);
  Clock& clock_;
  int handle_ = -1;
  bool open_ = false, running_ = false;
  DigitizerConfig cfg_;
  int recordLength_ = 0;
  char* buffer_ = nullptr;
  std::uint32_t bufferSize_ = 0;
  void* event_ = nullptr;
  nlohmann::json info_ = nlohmann::json::object();
  std::uint32_t memFullMask_ = 1u << 4;
  std::uint64_t memFullReads_ = 0, statusReads_ = 0, readDataCalls_ = 0, bytesRead_ = 0;
  std::uint32_t lastStatus_ = 0;
  int captured_ = 0;
  std::string lastError_;
};

}  // namespace daq::dig
