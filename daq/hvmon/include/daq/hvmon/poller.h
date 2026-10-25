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
// - Barrido de tensión (DAQ2-14, bias_scan.h): bias-begin registra el nominal; v0set escribe V0Set
//   con los límites de BiasSession (relectura previa, verificación de la rampa por V0Set releído y
//   VMon − V0Set − desvío ≤ 0.05 V, aborto si hay trip o Status inesperado); bias-restore lleva todos
//   los canales al nominal en escalones de ≤ 0.5 V y Pw = 1, y cierra la sesión. La sesión se persiste:
//   si hvmon arranca con una sesión activa, restaura el nominal solo.
// - measure-offsets: promedia IMon y VMon de los canales con Pw = 0 durante N s de lecturas
//   rápidas; rechaza los canales con Pw = 1.
// step() hace una iteración y duerme hasta el próximo vencimiento (máximo 200 ms): con un
// ManualClock los tests recorren horas en milisegundos.
// Seguridad entre hilos: la instancia pertenece al hilo de sondeo; se comunica solo por colas y
// por el StateHub.
#pragma once

#include <deque>
#include <map>
#include <set>
#include <future>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/hv/hv_backend.h"
#include "daq/hvmon/bias_scan.h"
#include "daq/hvmon/change_detector.h"
#include "daq/core/file_system.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/notifier.h"
#include "daq/hvmon/retry.h"
#include "daq/hvmon/state.h"
#include "daq/hvmon/writer.h"

namespace daq::hvmon {

struct PollRequest {
  // NominalV0Set (XPLUS-1): V0Set al nominal registrado por hvmon, solo con el canal apagado (recuperación tras un corte).
  enum class Kind { Pw, Offsets, BiasBegin, V0Set, BiasRestore, BiasStatus, NominalV0Set };
  Kind kind = Kind::Pw;
  int ch = -1;
  bool on = false;
  std::string reason;
  std::string requester;
  bool internal = false;          // pedida por hvmon mismo (auto-retry, prueba de grupos)
  double seconds = 0;             // offsets
  std::vector<int> channels;      // offsets (vacío = todos)
  std::uint64_t retryToken = 0;   // != 0: escritura pedida por la máquina de reintento
  bool clearAlarmFirst = false;   // ClrAlarm de la placa antes de Pw = 1 (clear_alarm_before_retry)
  double value = 0;               // V0Set pedido (V)
  std::map<int, double> expected; // bias-begin: V0Set nominal esperado por canal
  double maxAboveV = 0.0;         // bias-begin (DAQ2-17): V0Set máximo sobre el nominal (0 a 1.0 V)
  double maxStepV = -1.0;         // bias-begin (DAQ2-17): escalón de la sesión (≤ 0.5 V; ≤ 0: el de siempre)
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
  // Fase 2: persistencia de la máquina de reintento (carga el estado anterior) y cola de correo.
  void enableRetryPersistence(FileSystem& fs);
  void setNotifyQueue(BoundedQueue<NotifyEvent>* q) { notifyQ_ = q; }
  const RetryMachine& retry() const { return retry_; }
  // DAQ2-14: persistencia de la sesión de barrido (si estaba activa, restaura el nominal al conectar).
  void enableBiasPersistence(FileSystem& fs, const std::string& path);
  const BiasSession& bias() const { return bias_; }

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
    std::uint64_t retryToken = 0;
    bool v0set = false;   // DAQ2-14: escritura de V0Set (target)
    bool offWrite = false;  // V0Set nominal con el canal apagado (restauración tras un trip): sin rampa
    double target = 0;
    double v0setBefore = 0, vmonBefore = 0;  // relectura justo antes de escribir V0Set
    int stableReads = 0;                     // lecturas seguidas sin rampa y dentro de tolerancia
    bool internal = false;                   // DAQ2-17: escritura de hvmon (restauración): sin el límite de deriva del desvío
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
  void handleBiasBegin(PollRequest& r);
  void handleV0Set(PollRequest& r, SteadyTime now);
  void handleNominalV0Set(PollRequest& r, SteadyTime now);
  void handleBiasRestore(PollRequest& r, SteadyTime now);
  void advanceRestore(SteadyTime now);
  void persistBias();
  void handleOffsets(PollRequest& r, SteadyTime now);
  void verifyCommand(const hv::ChannelReading& r, SteadyTime now);
  void accumulateOffsets(const hv::ChannelReading& r, SteadyTime now);
  void finishCommand(const std::string& phase, const std::string& detail, bool ok, SteadyTime now);
  void record(const CommandRecord& rec);
  void pushWriter(WriterItem item);
  void addChanges(std::vector<ChangeRecord>& changes);
  void publish(SteadyTime now);
  void runRetry(const hv::ChannelReading& rd, const std::vector<ChangeRecord>& changes, SteadyTime now);
  void applyRetry(RetryOutput out, SteadyTime now);
  nlohmann::json channelContext(int ch) const;
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
  RetryMachine retry_;
  std::map<int, std::pair<WindowStats, WindowStats>> window_;  // canal -> (VMon, IMon), 600 s
  FileSystem* fs_ = nullptr;
  BoundedQueue<NotifyEvent>* notifyQ_ = nullptr;
  BiasSession bias_;
  std::string biasPath_;
  // Restauración en curso: escalones pendientes (canal, V0Set; V0Set < 0 = Pw on) y quien la pidió.
  std::deque<std::pair<int, double>> restoreSteps_;
  bool restoring_ = false, restoreFailed_ = false, autoRestorePending_ = false;
  std::string restoreReason_, restoreDetail_;
  std::shared_ptr<std::promise<nlohmann::json>> restoreReply_;
  std::optional<std::shared_future<nlohmann::json>> restoreFut_;
  std::set<int> restorePwTried_;  // verificación del escalón en curso
};

}  // namespace daq::hvmon
