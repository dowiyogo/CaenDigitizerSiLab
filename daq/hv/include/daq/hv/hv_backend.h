// hv_backend.h
// Interfaz de acceso a una fuente de HV CAEN. En la fase 0 es de SOLO LECTURA: no existe
// ningún método que escriba en la fuente (las escrituras llegan en las fases 1-2 con su propia
// autorización, ver docs/daq_v2_design.md, sección 9.2).
// Seguridad entre hilos: una instancia NO es segura entre hilos; la usa un único hilo dueño
// (en hvmon, el hilo HV). Unidades: ver hv_types.h (SI por encima del backend).
#pragma once

#include <string>
#include <vector>

#include "daq/core/result.h"
#include "daq/hv/hv_types.h"

namespace daq::hv {

class HVBackend {
 public:
  virtual ~HVBackend() = default;

  virtual Status connect() = 0;
  virtual void disconnect() = 0;
  virtual bool isConnected() const = 0;

  virtual Result<CrateMap> crateMap() = 0;
  // Nombres, tipos, modo y unidades de los parámetros de un canal.
  virtual Result<std::vector<ParamInfo>> channelParamNames(int slot, int channel) = 0;
  // Una lectura de los parámetros pedidos en los canales pedidos (todos con el mismo instante
  // de referencia tReq/tResp). Falla completa si falla cualquier parámetro.
  virtual Result<ChannelReading> readChannels(int slot, const std::vector<int>& channels,
                                              const std::vector<std::string>& params) = 0;

  virtual std::string describe() const = 0;
};

}  // namespace daq::hv
