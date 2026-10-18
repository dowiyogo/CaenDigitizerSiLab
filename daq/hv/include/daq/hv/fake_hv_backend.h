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
// - Bits de Status según la convención del manual de la Wrapper (Rev. 13).
// El comportamiento exacto del A7030DP tras un trip es una inferencia hasta la fase 2.
//
// Guion (JSON): tiempos en segundos desde la construcción del backend, medidos con el Clock
// inyectado. Eventos: "set" (cambio externo de un parámetro, p. ej. hecho en GECO), "current"
// (excursión de corriente), "connection_loss" (las llamadas fallan; después hay que llamar a
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

  Status connect() override;
  void disconnect() override;
  bool isConnected() const override { return connected_; }
  Result<CrateMap> crateMap() override;
  Result<std::vector<ParamInfo>> channelParamNames(int slot, int channel) override;
  Result<ChannelReading> readChannels(int slot, const std::vector<int>& channels,
                                      const std::vector<std::string>& params) override;
  // Como un "set" de Pw en el instante actual. Con "pw_write_ignored": true en el guion,
  // devuelve OK pero no cambia nada (para probar la verificación fallida).
  Status setChannelPw(int slot, int channel, bool on) override;
  std::string describe() const override { return "FakeHVBackend"; }

  // Solo para tests: tensión verdadera (sin offset ni ruido) en V, tras avanzar al instante actual.
  double trueVoltage(int slot, int channel);
  std::uint64_t callCount() const { return calls_; }
  std::uint64_t pwWrites() const { return pwWrites_; }

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
  };
  struct Slot {
    SlotInfo info;
    std::vector<Channel> channels;
    double imonNoiseUA = 0.0;
    double vmonNoiseV = 0.0;
    double imonStepUA = 0.002;
    double vmonStepV = 0.01;
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
  bool connected_ = false;
  std::uint64_t calls_ = 0;
  std::uint64_t pwWrites_ = 0;
  bool pwWriteIgnored_ = false;
  std::mt19937_64 rng_;
};

}  // namespace daq::hv
