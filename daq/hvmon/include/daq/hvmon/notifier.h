// notifier.h
// Notificador por correo de hvmon (fase 2), en un hilo propio: nunca bloquea al hilo de sondeo.
// Eventos: gave_up (siempre), trip (notify_each_trip), recovery (notify_recovery), writer_error y
// connection_lost (> connection_lost_min, por defecto 5 min), test (hvctl notify-test, marcado
// PRUEBA). Límites: un correo por canal y episodio para trip/recovery/gave_up; uno cada
// problem_interval_min (30 min) para writer_error y connection_lost.
// Envío: curl lanzado con posix_spawn (Spawner inyectable). Opciones y credenciales van por stdin
// ("--config -"), nunca en argv; el cuerpo va en un archivo temporal con modo 600 que se borra al
// terminar (no contiene secretos). Encabezados Date (RFC 5322) y Message-ID
// (<utc-ns.aleatorio@host>); DAQ_NOTIFY_TO puede tener varias direcciones separadas por comas.
// Si curl falla se reintenta 3 veces (1, 5 y 15 min); luego se registra y se publica el fallo.
// Los secretos no aparecen en argv, logs, estado ni en el cuerpo del mensaje.
// Seguridad entre hilos: la instancia pertenece al hilo notificador; recibe eventos por
// BoundedQueue y publica su estado en StateHub::notifier.
#pragma once

#include <chrono>
#include <deque>
#include <future>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/core/jthread.h"
#include "daq/core/result.h"
#include "daq/core/secrets.h"
#include "daq/core/spawn.h"
#include "daq/hvmon/mapping.h"

namespace daq::hvmon {

class StateHub;

struct NotifierConfig {
  bool enabled = false;
  bool notifyEachTrip = false;
  bool notifyRecovery = false;
  bool notifyWriterErrors = true;
  bool notifyConnectionLost = true;
  double connectionLostMin = 5;
  double problemIntervalMin = 30;
  std::vector<double> retryDelaysMin{1, 5, 15};
  std::string curlPath = "curl";
  int maxTimeS = 30;
  std::string mappingFile = "daq/config/mapping/hodo_mapping_2026_10.json";
  std::string localTimezone = "America/Santiago";  // hora local entre paréntesis en los correos
  std::string pageSshUser = "xplus";                // línea "ssh -N -L ..." del correo
  std::string pageSshHost;                          // host Tailscale; vacío: la línea lo indica
  int pagePort = 8081;                              // = http.port de hvmon (no es clave de notifier)
  static Result<NotifierConfig> fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
};

struct SmtpSettings {
  std::string url, user, password, from;
  std::vector<std::string> to;
  // DAQ_SMTP_URL, DAQ_SMTP_USER, DAQ_SMTP_PASSWORD, DAQ_NOTIFY_FROM, DAQ_NOTIFY_TO.
  static Result<SmtpSettings> fromSecrets(Secrets& s);
};

struct NotifyEvent {
  std::string kind;  // gave_up, trip, recovery, writer_error, connection_lost
  int ch = -1;       // canal HV (-1 si no aplica)
  int slot = 3;
  std::uint64_t episode = 0;
  int k = 0;
  std::string reason;
  nlohmann::json timeline = nlohmann::json::array();
  nlohmann::json context = nlohmann::json::object();  // VMon_V, IMon_A, Status, status_flags, board_temperature_degC, geco_running
  UtcNs tUtcNs = 0;
  bool test = false;       // PRUEBA: sin límite por episodio y sin reintentos
  bool simulated = false;  // canal en test_mode (la fuente simulada se marca con setSimulatedBackend)
  std::shared_ptr<std::promise<nlohmann::json>> reply;  // opcional (notify-test)
};

class Notifier {
 public:
  Notifier(NotifierConfig cfg, SmtpSettings smtp, std::optional<Mapping> mapping, Spawner& spawner, Clock& clock,
           std::string tmpDir, StateHub& hub, std::string hostName = "");

  void run(BoundedQueue<NotifyEvent>& queue, const StopToken& st);
  // hvmon con FakeHVBackend: todo correo sale marcado [SIMULADO] (o [PRUEBA]).
  void setSimulatedBackend(bool v) { simulatedBackend_ = v; }

  // Para tests (y usados por run()).
  bool submit(NotifyEvent e);  // aplica habilitación y límites; false si se descarta
  void processDue();           // envía lo pendiente cuyo momento llegó
  void checkProblems();        // errores del escritor y conexión perdida, desde el StateHub
  std::string composeMessage(const NotifyEvent& e, std::string* subject = nullptr) const;
  std::string curlConfig(const std::string& bodyPath) const;
  nlohmann::json status() const;
  std::size_t pending() const { return pending_.size(); }

 private:
  struct Pending {
    NotifyEvent ev;
    std::string message;
    int attempts = 0;
    SteadyTime due{};
  };
  bool sendOnce(Pending& p, std::string& err, int& exitCode);
  void publish();
  std::string redact(std::string s) const;
  std::string rfc5322Date(UtcNs t) const;
  std::string groupText(int hvCh) const;

  NotifierConfig cfg_;
  SmtpSettings smtp_;
  std::optional<Mapping> mapping_;
  Spawner& spawner_;
  Clock& clock_;
  std::string tmpDir_;
  StateHub& hub_;
  std::string host_;
  std::deque<Pending> pending_;
  std::map<std::string, bool> sentKeys_;           // tipo:canal:episodio
  std::map<std::string, SteadyTime> lastProblem_;  // writer_error / connection_lost
  std::uint64_t lastWriterErrors_ = 0;
  bool connectionAlerted_ = false;
  std::uint64_t sent_ = 0, failed_ = 0, suppressed_ = 0;
  nlohmann::json lastResult_;
  std::uint64_t seq_ = 0;
  bool simulatedBackend_ = false;
};

}  // namespace daq::hvmon
