// retry.h
// Máquina de estados de reintento automático tras un trip, por canal HV (fase 2). Lógica pura:
// recibe las lecturas rápidas y los cambios detectados, y devuelve transiciones, acciones (pedir
// Pw = 1, opcionalmente ClrAlarm antes) y eventos para el notificador. La ejecución de las
// escrituras la hace el hilo de sondeo por la ruta con compuerta (allow_commands).
//
// Estados: NORMAL, WAIT(k), RETRY(k), STABILIZING(k), GAVE_UP, MANUAL, DISABLED. Tabla completa de
// transiciones en docs/daq_v2_design.md, sección 8. Todos los tiempos son UTC (ns) para que la
// espera sobreviva a un reinicio de hvmon; el estado se persiste en cada transición.
// Seguridad entre hilos: la usa solo el hilo de sondeo.
#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/core/time.h"
#include "daq/hv/status_bits.h"
#include "daq/hvmon/records.h"

namespace daq::hvmon {

enum class RetryState { Uninit, Disabled, Normal, Wait, Retry, Stabilizing, GaveUp, Manual };
const char* toString(RetryState s);

struct RetryConfig {
  std::vector<int> enabledChannels;            // vacío por defecto: ningún canal con reintento
  std::vector<double> waitMin{5, 30, 60};      // espera antes de cada intento k (min)
  double stableMin = 10;                       // ventana de estabilidad (min)
  double resetAfterH = 2;                      // NORMAL continuo para reiniciar el contador (h)
  double imonStableMaxA = 6e-6;                // IMon - offset máximo en la estabilidad (A)
  int maxRampingChannels = 1;
  bool inhibitWhenGecoRunning = false;
  bool clearAlarmBeforeRetry = false;
  bool testModeEnabled = false;
  std::vector<int> testChannels;
  std::vector<double> testWaitMin{1, 2, 3};
  double testStableMin = 1;
  bool testModeAllowDetectorChannels = false;  // test_mode en CH00-07 solo con esto en true
  std::string stateFile;                       // $XDG_STATE_HOME/daqv2/hvmon_retry_state.json

  // JSON con las claves de la tarea DAQ2-4 (imon_stable_max_uA en µA; internamente A).
  static Result<RetryConfig> fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
  bool enabled(int ch) const;
  bool inTestMode(int ch) const;
  const std::vector<double>& waits(int ch) const;
  double stableMinFor(int ch) const;
};

std::string defaultRetryStateFile();

// Lo que la máquina necesita de una lectura rápida de un canal (SI).
struct RetryChannelView {
  std::uint32_t pw = 0;
  std::uint32_t status = 0;
  double vmonV = 0, imonA = 0, v0setV = 0;
};

struct RetryInput {
  UtcNs nowUtc = 0;
  bool fresh = false;  // lectura buena recién hecha (conexión sana)
  bool geco = false;
  std::map<int, RetryChannelView> channels;
  std::vector<ChangeRecord> changes;   // cambios detectados en esta lectura
  std::map<int, double> imonOffsetA;   // offsets medidos (Pw = 0); 0 si no hay
};

struct RetryAction {
  int ch = 0;
  bool clearAlarmFirst = false;
  std::uint64_t token = 0;  // identifica el intento: los resultados de otro intento se ignoran
  int k = 0;
};

struct RetryTransition {
  int ch = 0;
  std::string from, to;  // con k, p. ej. "WAIT(2)"
  int k = 0;
  std::string reason;
  UtcNs tUtcNs = 0;
};

struct RetryEvent {
  std::string kind;  // "trip", "recovery", "gave_up"
  int ch = 0;
  std::uint64_t episode = 0;
  int k = 0;
  std::string reason;
  nlohmann::json timeline;
  UtcNs tUtcNs = 0;
};

struct RetryOutput {
  std::vector<RetryAction> actions;
  std::vector<RetryTransition> transitions;
  std::vector<RetryEvent> events;
  std::vector<std::string> warnings;
  bool persist = false;
  void append(RetryOutput&& o);
};

class RetryMachine {
 public:
  RetryMachine(RetryConfig cfg, hv::StatusBitMap bits, std::vector<std::string> tripFlags, std::vector<int> channels);

  // Estado persistido de una ejecución anterior; se aplica en la primera lectura de cada canal.
  void restore(const nlohmann::json& persisted);
  RetryOutput update(const RetryInput& in);
  // Resultado de la escritura pedida por una acción (fase del comando de Pw en el sondeo).
  RetryOutput onCommandResult(int ch, std::uint64_t token, const std::string& phase, const std::string& detail,
                              UtcNs nowUtc);

  nlohmann::json publish(UtcNs nowUtc) const;      // arreglo con un objeto por canal
  nlohmann::json persistJson(UtcNs nowUtc) const;  // archivo de persistencia
  RetryState state(int ch) const;
  int attempt(int ch) const;
  const RetryConfig& config() const { return cfg_; }

 private:
  struct Ch {
    RetryState st = RetryState::Uninit;
    int k = 0;
    UtcNs tripUtc = 0;
    UtcNs sinceUtc = 0;        // entrada al estado actual
    UtcNs stableSinceUtc = 0;  // 0 = ventana de estabilidad no iniciada
    UtcNs normalSinceUtc = 0;
    std::uint32_t waitStatus = 0;
    double waitV0 = 0;
    std::string lastOutcome;
    std::string blockReason;
    std::uint64_t episode = 0;
    std::uint64_t token = 0;
    bool gecoWarned = false;
    bool tripLatched = false;  // XPLUS-1: canal sin reintento con un trip ya avisado (se libera al volver encendido sin trip)
    nlohmann::json timeline = nlohmann::json::array();
    // Persistido (se aplica en la primera lectura)
    RetryState pState = RetryState::Uninit;
    int pK = 0;
    UtcNs pTripUtc = 0;
    std::uint64_t pEpisode = 0;
    nlohmann::json pTimeline;
  };

  static std::string label(RetryState s, int k);
  // fromLabel: etiqueta de origen si k ya cambió antes de la transición (vacío = estado actual).
  void go(Ch& c, int ch, RetryState to, const std::string& reason, UtcNs now, RetryOutput& out,
          const std::string& fromLabel = "");
  void failAttempt(Ch& c, int ch, const RetryChannelView& v, const std::string& reason, UtcNs now, RetryOutput& out,
                   bool trip);
  void enterWait(Ch& c, int ch, const RetryChannelView& v, int k, const std::string& reason, UtcNs now, RetryOutput& out);
  void initial(Ch& c, int ch, const RetryChannelView& v, UtcNs now, RetryOutput& out);
  bool hasTrip(std::uint32_t status) const;
  bool ramping(std::uint32_t status) const;
  bool stable(const RetryChannelView& v, double offsetA) const;
  std::uint32_t maskRamp(std::uint32_t status) const;
  void event(Ch& c, int ch, const std::string& kind, const std::string& reason, UtcNs now, RetryOutput& out);

  RetryConfig cfg_;
  hv::StatusBitMap bits_;
  std::vector<std::string> tripFlags_;
  std::vector<int> channels_;
  std::map<int, Ch> ch_;
  std::uint64_t nextToken_ = 1;
  std::uint64_t nextEpisode_ = 1;
};

}  // namespace daq::hvmon
