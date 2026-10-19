// caen_hv_backend.h
// Backend sobre CAEN HV Wrapper 6.6 (libcaenhvwrapper). Usa CAENHV_InitSystem/DeinitSystem,
// funciones CAENHV_Get* (canal y placa) y dos escrituras: CAENHV_SetChParam("Pw") en
// setChannelPw() y CAENHV_SetBdParam("ClrAlarm") en clearBoardAlarm(); ninguna otra función escribe.
// - Es dueño del handle (RAII): el destructor llama a DeinitSystem.
// - Libera con CAENHV_Free cada buffer que reserva la Wrapper.
// - Convierte a SI en esta frontera: valor_SI = valor_placa * 10^Exp, con la unidad y el
//   exponente que informa la placa (GetChParamProp "Unit"/"Exp"); si la placa no los informa se
//   usa una tabla de respaldo por nombre (IMon/I0Set/I1Set en µA, resto V, V/s, s).
// Seguridad entre hilos: no es segura; un solo hilo dueño. El header de la Wrapper (que define
// macros uchar/ushort/ulong) solo se incluye en el .cpp.
#pragma once

#include <map>
#include <memory>
#include <string>
#include <utility>

#include "daq/core/clock.h"
#include "daq/hv/credentials.h"
#include "daq/hv/hv_backend.h"

namespace daq::hv {

struct ConnectionConfig {
  std::string system = "SY5527";  // tipo de sistema de la Wrapper
  std::string host;               // IP del mainframe
  Credentials credentials;
};

class CaenHVBackend final : public HVBackend {
 public:
  // El reloj debe vivir más que el backend; nullptr usa un SystemClock interno.
  explicit CaenHVBackend(ConnectionConfig config, const Clock* clock = nullptr);
  ~CaenHVBackend() override;
  CaenHVBackend(const CaenHVBackend&) = delete;
  CaenHVBackend& operator=(const CaenHVBackend&) = delete;

  Status connect() override;
  void disconnect() override;
  bool isConnected() const override { return handle_ >= 0; }

  Result<CrateMap> crateMap() override;
  Result<std::vector<ParamInfo>> channelParamNames(int slot, int channel) override;
  Result<ChannelReading> readChannels(int slot, const std::vector<int>& channels,
                                      const std::vector<std::string>& params) override;
  Result<std::vector<ParamInfo>> boardParamNames(int slot) override;
  Result<BoardReading> readBoard(int slot, const std::vector<std::string>& params) override;
  Status setChannelPw(int slot, int channel, bool on) override;
  Status clearBoardAlarm(int slot) override;
  std::string describe() const override;

  // Versión de la biblioteca (CAENHVLibSwRel).
  static std::string librarySoftwareRelease();

 private:
  Result<ParamInfo> paramInfo(int slot, int channel, const std::string& name);
  Result<ParamInfo> boardParamInfo(int slot, const std::string& name);
  Error wrapperError(int code, const std::string& where) const;

  ConnectionConfig config_;
  std::unique_ptr<Clock> ownClock_;
  const Clock* clock_;
  int handle_ = -1;
  std::map<std::pair<int, std::string>, ParamInfo> infoCache_;  // (slot, nombre)
  std::map<std::pair<int, std::string>, ParamInfo> bdInfoCache_;
};

}  // namespace daq::hv
