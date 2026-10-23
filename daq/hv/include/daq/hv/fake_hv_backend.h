// fake_hv_backend.h
// Backend HV simulado y guionable, para tests y para correr hvctl sin hardware (--fake).
// Reproduce el conjunto de 18 parámetros del A7030DP con los mismos nombres, tipos, modos y
// unidades que informa la placa real, de modo que la salida de hvctl se vea igual.
//
// Modelo por canal (unidades de la placa: V, µA, V/s, s):
// - Pw = 1 lleva la tensión verdadera hacia V0Set a RUp. Al apagar (a mano o por trip) la baja a
//   RDWn si PDwn = 1 ("Ramp") o la deja caer de inmediato si PDwn = 0 ("Kill", valor actual
//   de la placa). Esta interpretación de PDwn es una inferencia a verificar en la fase 2.
// - Corriente verdadera = load_uA·(V/V0Set) + excursiones guionadas·(V/V0Set).
// - Si la corriente supera I0Set con Pw = 1: bit overcurrent, IMon limitado a I0Set, y tras
//   Trip segundos seguidos el canal se apaga (Pw = 0, bit on = 0, bit internal_trip = 1). No
//   se vuelve a encender solo; internal_trip se borra al poner Pw = 1 (evento "set").
// - VMon/IMon = verdadero + offset + ruido gaussiano, cuantizados al paso de la placa.
// - "discharge_tau_s" (slot, DAQ2-14): con Pw = 1, una bajada de V0Set sigue V0Set + (V − V0Set)·e^(−t/τ)
//   (descarga de la capacidad del PCB del detector) y el bit ramp_down queda puesto mientras V supere
//   V0Set en más de 0.01 V. 0 (por omisión): bajada lineal a RDWn como antes.
// - Bits de Status según la convención del manual de la Wrapper (Rev. 13).
// El comportamiento exacto del A7030DP tras un trip es una inferencia hasta la fase 2.
//
// Guion (JSON): tiempos en segundos desde la construcción del backend, medidos con el Clock
// inyectado. Eventos: "set" (cambio externo de un parámetro, p. ej. hecho en GECO), "current"
// (excursión de corriente), "board" (parámetro de placa, p. ej. Temp), "connection_loss" (las llamadas fallan; después hay que llamar a
// connect()). Con "base": "default" el guion usa los slots del guion por defecto y solo
// aporta eventos y semilla. Ver daq/config/fake/README.md.
// Seguridad entre hilos: como HVBackend, un solo hilo dueño.
#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/hv/hv_backend.h"

namespace daq::hv {

class FakeHVBackend final : public HVBackend {
 public:
  // Lanza std::invalid_argument si el guion es inválido. El reloj debe vivir más que el backend.
  FakeHVBackend(const nlohmann::json& scenario, const Clock& clock);
  static Result<std::unique_ptr<FakeHVBackend>> fromFile(const std::string& path, const Clock& clock);

  // A7030DP en el slot 3 con los valores leídos el 2026-10-18 con HV apagada.
  static nlohmann::json defaultScenario();
  // Tabla de parámetros del A7030DP (nombres, tipos, modos, unidades, rangos).
  static std::vector<ParamInfo> a7030dpParams();
  // Parámetros de placa del A7030DP (BdStatus, HVMax, HIMax, Temp, ClrAlarm).
  static std::vector<ParamInfo> a7030dpBoardParams();

  Status connect() override;
  void disconnect() override;
  bool isConnected() const override { return connected_; }
  Result<CrateMap> crateMap() override;
  Result<std::vector<ParamInfo>> channelParamNames(int slot, int channel) override;
  Result<ChannelReading> readChannels(int slot, const std::vector<int>& channels,
                                      const std::vector<std::string>& params) override;
  // Como un "set" de Pw en el instante actual. Con "pw_write_ignored": true en el guion,
  // devuelve OK pero no cambia nada (para probar la verificación fallida).
  Result<std::vector<ParamInfo>> boardParamNames(int slot) override;
  Result<BoardReading> readBoard(int slot, const std::vector<std::string>& params) override;
  Status setChannelPw(int slot, int channel, bool on) override;
  Status setChannelV0Set(int slot, int channel, double volts) override;
  // Fake: borra los bits de trip retenidos del slot sin encender nada [I: a verificar con la placa].
  Status clearBoardAlarm(int slot) override;
  std::uint64_t clearAlarmWrites() const { return clearAlarmWrites_; }
  // Para tests: secuencia de escrituras ("Pw:<ch>=<v>", "ClrAlarm") en orden.
  const std::vector<std::string>& writeLog() const { return writeLog_; }
  std::string describe() const override { return "FakeHVBackend"; }

  // Solo para tests: tensión verdadera (sin offset ni ruido) en V, tras avanzar al instante actual.
  double trueVoltage(int slot, int channel);
  std::uint64_t callCount() const { return calls_; }
  std::uint64_t pwWrites() const { return pwWrites_; }
  std::uint64_t v0setWrites() const { return v0setWrites_; }

 private:
  struct Channel {
    std::map<std::string, double> num;          // parámetros numéricos (unidades de la placa)
    std::map<std::string, std::uint32_t> ints;  // ON/OFF, binarios
    double vTrue = 0.0;
    double ovcTime = 0.0;
    bool overcurrent = false;
    bool tripped = false;
    double loadUA = 0.0;
    double extraUA = 0.0;
    double imonOffsetUA = 0.0;
    double vmonOffsetV = 0.0;
    double dischargeTauS = 0.0;  // DAQ2-14: > 0: bajadas con Pw = 1 exponenciales (capacidad del PCB, τ)
  };
  struct Slot {
    SlotInfo info;
    std::vector<Channel> channels;
    double imonNoiseUA = 0.0;
    double vmonNoiseV = 0.0;
    double imonStepUA = 0.002;
    double vmonStepV = 0.01;
    std::map<std::string, double> board;  // parámetros de placa (unidades de la placa)
  };
  struct Event {
    double t = 0.0;
    std::string type;
    int slot = -1;
    int channel = -1;
    std::string param;
    double value = 0.0;
    double duration = 0.0;
  };

  void advanceToNow();
  void integrate(double dt);
  void apply(const Event& e);
  bool inConnectionLoss(double t) const;
  Status checkLink();
  Slot* findSlot(int slot);
  double seconds(SteadyTime t) const;
  std::uint32_t statusWord(const Channel& c) const;
  ParamValue value(const Slot& s, const Channel& c, const ParamInfo& info);

  const Clock& clock_;
  SteadyTime t0_;
  double simTime_ = 0.0;
  int crateSlots_ = 0;
  std::vector<Slot> slots_;
  std::vector<Event> events_;        // ordenados por tiempo, pendientes
  std::vector<std::pair<double, double>> losses_;  // [inicio, fin) de pérdidas de conexión
  std::vector<ParamInfo> params_;
  std::vector<ParamInfo> boardParams_;
  bool connected_ = false;
  std::uint64_t calls_ = 0;
  std::uint64_t pwWrites_ = 0;
  std::uint64_t v0setWrites_ = 0;
  std::uint64_t clearAlarmWrites_ = 0;
  std::vector<std::string> writeLog_;
  bool pwWriteIgnored_ = false;
  std::mt19937_64 rng_;
};

}  // namespace daq::hv
