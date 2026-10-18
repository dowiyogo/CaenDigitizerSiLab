// hv_backend.h
// Interfaz de acceso a una fuente de HV CAEN. Solo hay UNA operación de escritura en todo el
// código: setChannelPw() (Pw de un canal a 0 o 1). No existe escritura de V0Set, I0Set, Trip ni de
// ningún otro parámetro. setChannelPw() solo la llama el hilo de sondeo de hvmon, tras validar
// allow_commands y la secuencia de docs/daq_v2_design.md, sección 9.2.
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

  // Escribe Pw de un canal. No verifica: la verificación por relectura la hace hvmon.
  virtual Status setChannelPw(int slot, int channel, bool on) = 0;

  virtual std::string describe() const = 0;
};

}  // namespace daq::hv
