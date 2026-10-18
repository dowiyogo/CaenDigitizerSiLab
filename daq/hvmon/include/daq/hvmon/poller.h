// poller.h
// Hilo de sondeo de hvmon: único dueño del backend HV (y por lo tanto de la sesión de la
// Wrapper). Cada fast_period_s lee el conjunto rápido (V0Set, I0Set, VMon, IMon, Status, Pw) de
// todos los canales; cada full_period_s todos los parámetros. Tras un error desconecta y
// reintenta con espera exponencial (1, 2, 4, ... hasta 60 s). Detecta cambios, alimenta la cola
// del escritor y publica PollerState en el StateHub.
//
// También ejecuta las peticiones que llegan por la cola de peticiones (desde el hilo de comandos):
// - Pw (única escritura): rechazada si allow_commands = false. Secuencia: relectura del canal
//   justo antes de escribir y aborto si difiere de la última lectura publicada (carrera con GECO
//   o con un trip); escritura; verificación en las lecturas rápidas siguientes hasta un timeout
//   V0Set/RUp (o V0Set/RDWn) + margen; registro de cada fase en el archivo de cambios. Un comando
//   a la vez; si la petición viene de hvmon mismo (internal), además ningún canal en rampa.
// - measure-offsets: promedia IMon y VMon de los canales con Pw = 0 durante N s de lecturas
//   rápidas; rechaza los canales con Pw = 1.
// step() hace una iteración y duerme hasta el próximo vencimiento (máximo 200 ms): con un
// ManualClock los tests recorren horas en milisegundos.
// Seguridad entre hilos: la instancia pertenece al hilo de sondeo; se comunica solo por colas y
// por el StateHub.
#pragma once

#include <deque>
#include <future>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/hv/hv_backend.h"
#include "daq/hvmon/change_detector.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/state.h"
#include "daq/hvmon/writer.h"

namespace daq::hvmon {

struct PollRequest {
  enum class Kind { Pw, Offsets };
  Kind kind = Kind::Pw;
  int ch = -1;
  bool on = false;
  std::string reason;
  std::string requester;
  bool internal = false;          // pedida por hvmon mismo (auto-retry, prueba de grupos)
  double seconds = 0;             // offsets
  std::vector<int> channels;      // offsets (vacío = todos)
  std::shared_ptr<std::promise<nlohmann::json>> reply = std::make_shared<std::promise<nlohmann::json>>();
};

class Poller {
 public:
  Poller(const HvmonConfig& cfg, hv::HVBackend& backend, Clock& clock, BoundedQueue<WriterItem>& writer,
         StateHub& hub, BoundedQueue<PollRequest>& requests);
  ~Poller();

  void run(const StopToken& st);
  void step(const StopToken& st = {});
  // Responde con error a todo lo pendiente (al apagar).
  void shutdown();

  // Para tests.
  const ReadStats& stats() const { return stats_; }
  const std::vector<SteadyTime>& connectAttemptTimes() const { return attemptTimes_; }
  bool connected() const { return connected_; }

 private:
  struct ActiveCommand {
    std::uint64_t id = 0;
    int ch = 0;
    bool on = false;
    SteadyTime deadline{};
    std::string reason, requester;
    std::shared_ptr<std::promise<nlohmann::json>> reply;
  };
  struct ActiveOffsets {
    std::vector<int> channels;
    std::vector<int> refused;
    double seconds = 0;
    SteadyTime tStart{};
    UtcNs tStartUtc = 0;
    UtcNs tEndUtc = 0;
    bool started = false;
    std::map<int, double> sumI, sumI2, sumV, sumV2;
    std::uint64_t n = 0;
    std::shared_ptr<std::promise<nlohmann::json>> reply;
  };

  void tryConnect(SteadyTime now);
  void onReadError(const Error& e, SteadyTime now);
  void fastRead(SteadyTime now);
  void fullRead(SteadyTime now);
  void handleRequest(PollRequest& r, SteadyTime now);
  void handlePw(PollRequest& r, SteadyTime now);
  void handleOffsets(PollRequest& r, SteadyTime now);
  void verifyCommand(const hv::ChannelReading& r, SteadyTime now);
  void accumulateOffsets(const hv::ChannelReading& r, SteadyTime now);
  void finishCommand(const std::string& phase, const std::string& detail, bool ok, SteadyTime now);
  void record(const CommandRecord& rec);
  void pushWriter(WriterItem item);
  void addChanges(std::vector<ChangeRecord>& changes);
  void publish(SteadyTime now);
  const hv::ParamValue* lastValue(int ch, const std::string& param) const;
  bool anyRamping() const;
  SteadyTime nextWake(SteadyTime now) const;
  std::chrono::nanoseconds secondsToNs(double s) const;

  const HvmonConfig& cfg_;
  hv::HVBackend& backend_;
  Clock& clock_;
  BoundedQueue<WriterItem>& writer_;
  StateHub& hub_;
  BoundedQueue<PollRequest>& requests_;
  ChangeDetector detector_;

  bool connected_ = false;
  bool everConnected_ = false;
  double backoffS_;
  SteadyTime nextConnect_{};
  SteadyTime nextFast_{};
  SteadyTime nextFull_{};
  bool fastGap_ = false;  // hubo desconexión desde la última lectura rápida buena
  bool fullGap_ = false;
  std::vector<std::string> allParams_;
  std::vector<std::string> fullTracked_;
  std::vector<std::string> boardParams_;  // vacío si la placa no informa parámetros
  std::optional<hv::ChannelReading> prevFast_;
  std::optional<hv::ChannelReading> prevFull_;
  std::uint64_t readSeq_ = 0;
  std::uint64_t nextChangeId_ = 1;
  std::uint64_t nextCommandId_ = 1;
  ReadStats stats_;
  std::vector<SteadyTime> attemptTimes_;
  PollerState state_;
  std::optional<ActiveCommand> active_;
  std::optional<ActiveOffsets> offsets_;
};

}  // namespace daq::hvmon
