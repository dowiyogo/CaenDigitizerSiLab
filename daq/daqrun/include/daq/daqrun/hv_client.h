// hv_client.h
// Cliente HV de daqrun (diseño 2.1, 6.1, 6.3). Lee el estado publicado por hvmon
// (/dev/shm/daqv2/hvmon_state.json, o el del directorio simulado) cada poll_period_s (1 s) y
// manda cada muestra al hilo escritor por una cola. No abre sesión con el SY5527.
//  - alive = el archivo existe, no es "final", su t_utc_ns tiene menos de stale_after_s y el pid
//    de hvmon existe; stale = hvmon dice stale o alive = false.
//  - FlagEngine (en el hilo escritor) guarda el historial de lecturas (deduplicado por read_seq,
//    con t_req/t_resp de last_read) y calcula para un evento en t los bits por canal:
//      0 HV_OFF, 1 HV_OVC, 2 HV_RAMP, 3 HV_VDEV, 4 HV_UNKNOWN, 5 HV_AMBIGUOUS, 9 MAP_UNCONFIRMED.
//    Bit 4: sin lectura previa, lectura previa con más de stale_after_s, hvmon no vivo o stale en
//    t, sin lectura posterior al liberar el evento (línea de retardo flag_delay_s), o grupo HV del
//    canal desconocido. Bit 5: t en (t_req(a) − backdate, t_resp(b)] de un par consecutivo de
//    lecturas a, b con estado distinto del canal HV (Pw, Status, V0Set); backdate = trip_backdate_s
//    si el cambio es un trip; dentro del intervalo los bits 0-3 son el OR de ambos estados.
//    Bit 9: hv_group o face del canal sin confirmed = measurement (estático, del mapeo).
// Seguridad entre hilos: HvClient::run() en su hilo; latest() se lee desde cualquier hilo
// (Snapshot). FlagEngine pertenece al hilo escritor.
#pragma once

#include <array>
#include <cstdint>
#include <deque>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/core/jthread.h"
#include "daq/core/snapshot.h"
#include "daq/daqrun/config.h"
#include "daq/hvmon/mapping.h"
#include "daq/hvmon/notifier.h"
#include "daq/ttt.h"

namespace daq::daqrun {

inline constexpr int kHvChannels = 12;

namespace chflag {
inline constexpr std::uint16_t HvOff = 1u << 0, HvOvc = 1u << 1, HvRamp = 1u << 2, HvVdev = 1u << 3,
                               HvUnknown = 1u << 4, HvAmbiguous = 1u << 5, DataBlackout = 1u << 6,
                               LightLeak = 1u << 7, RetryStab = 1u << 8, MapUnconfirmed = 1u << 9,
                               PeProvisional = 1u << 10, StatInsufficient = 1u << 11;
}
namespace evtflag {
inline constexpr std::uint32_t Gap = 1u << 0, Burst = 1u << 1, HwCounterJump = 1u << 2, BoardFull = 1u << 3,
                               TttPcMismatch = 1u << 4, HvAnyFlag = 1u << 5, MonitorProvisional = 1u << 6,
                               NoCalibration = 1u << 7, HvGroupTest = 1u << 8;
}

struct HvChannelSample {
  bool present = false;
  std::uint32_t pw = 0, status = 0;
  double v0setV = 0, i0setA = 0, vmonV = 0, imonA = 0;
  std::string retryState;
};

struct HvSample {
  UtcNs receivedUtcNs = 0;  // hora de daqrun al leer el archivo
  bool fileOk = false;      // se pudo leer y parsear
  bool alive = false;
  bool stale = true;
  bool connected = false;
  std::string error;
  std::uint64_t stateSeq = 0, readSeq = 0;
  UtcNs stateUtcNs = 0, tReqUtcNs = 0, tRespUtcNs = 0;
  long pid = 0;
  std::optional<double> boardTempDegC;  // temperatura de la placa A7030DP (no del detector)
  std::array<HvChannelSample, kHvChannels> ch{};
  nlohmann::json newChanges = nlohmann::json::array();  // cambios con id > el último visto
  nlohmann::json raw;  // estado completo (solo en la primera muestra: cabecera del run)
};

// Parseo puro (tests). now: hora de lectura; pidAlive: comprobación del pid (inyectable).
HvSample parseHvState(const nlohmann::json& state, UtcNs now, double staleAfterS, bool pidAlive);

class HvClient {
 public:
  HvClient(const DaqrunConfig& cfg, Clock& clock, BoundedQueue<HvSample>& toWriter);
  // Se llama en el hilo del cliente con cada muestra (p. ej. el vigilante de hvmon).
  void setOnSample(std::function<void(const HvSample&)> f) { onSample_ = std::move(f); }
  void run(const StopToken& st);
  HvSample pollOnce();  // lee y parsea una vez (también la usa main para la cabecera)
  std::pair<HvSample, std::uint64_t> latest() const { return latest_.read(); }
  std::uint64_t pushedCount() const { return pushed_; }

 private:
  const DaqrunConfig& cfg_;
  Clock& clock_;
  BoundedQueue<HvSample>& out_;
  Snapshot<HvSample> latest_;
  std::int64_t lastChangeId_ = -1;
  std::uint64_t pushed_ = 0;
  std::function<void(const HvSample&)> onSample_;
};

// Vigilante de hvmon durante un run: si el estado no está vivo o está stale más de thresholdS
// seguidos, devuelve una vez (por episodio) un evento "hvmon_stale" para el notificador; el run
// sigue con el bit 4. Lógica pura, reloj de las muestras.
class HvmonWatchdog {
 public:
  HvmonWatchdog(double thresholdS, std::string origin) : thresholdS_(thresholdS), origin_(std::move(origin)) {}
  std::optional<hvmon::NotifyEvent> feed(const HvSample& s);
  bool alerted() const { return alerted_; }

 private:
  double thresholdS_;
  std::string origin_;
  UtcNs badSince_ = 0;
  bool alerted_ = false;
  std::string lastError_;
  std::uint64_t episode_ = 0;
};

struct ChannelFlags {
  std::array<std::uint16_t, 64> ch{};
  std::uint64_t hvSeq = 0;  // read_seq de la lectura usada
  bool anyHv = false;       // OR de los bits 0-5 de algún canal
};

class FlagEngine {
 public:
  FlagEngine(const DaqrunConfig& cfg, const hvmon::Mapping& mapping);
  void addSample(const HvSample& s);
  // Flags de un evento en t. Debe llamarse cuando ya pasó flag_delay_s (o al cierre).
  ChannelFlags flagsFor(UtcNs t) const;
  void prune(UtcNs keepAfter);  // descarta historial anterior
  std::uint16_t staticFlags(int ch) const { return static_[static_cast<std::size_t>(ch)]; }
  int hvChannelOf(int digCh) const { return hvOf_[static_cast<std::size_t>(digCh)]; }
  std::size_t historySize() const { return reads_.size(); }

 private:
  struct Read {
    std::uint64_t readSeq;
    UtcNs tReq, tResp;
    std::array<HvChannelSample, kHvChannels> ch;
  };
  struct Liveness {
    UtcNs t;
    bool ok;  // alive && !stale
  };
  std::uint16_t stateBits(const HvChannelSample& c, int h) const;
  bool differs(const HvChannelSample& a, const HvChannelSample& b) const;
  bool isTrip(const HvChannelSample& a, const HvChannelSample& b) const;
  double imonBaseline(int h) const;

  const DaqrunConfig& cfg_;
  std::array<int, 64> hvOf_{};  // canal HV del canal de digitizador, −1 si el grupo es desconocido
  std::array<std::uint16_t, 64> static_{};
  std::deque<Read> reads_;
  std::deque<Liveness> live_;
  std::array<std::deque<std::pair<UtcNs, double>>, kHvChannels> imonHist_;
};

// Vigilante de triggers (runs con trigger externo). Episodio "sin triggers" si (a) no llega ningún
// trigger durante thresholdS, o (b) la tasa en la ventana de windowS cae bajo minRateHz (triggers
// sueltos cada pocos minutos, como con el cable de trigger fallando en el run005, DAQ2-6). Devuelve "no_triggers" (hora
// del último trigger, tasa de la ventana, run, campaña, temperatura de la placa HV) y, cuando
// vuelve la normalidad, "triggers_resumed" con la duración del episodio. Uno de cada por episodio.
// Lógica pura: la hora la pone quien llama (hora de lectura de los bloques y del reloj).
class TriggerWatchdog {
 public:
  TriggerWatchdog(double thresholdS, std::string origin, UtcNs runStartUtc, double minRateHz = 0.0, double windowS = 600.0)
      : thresholdS_(thresholdS), origin_(std::move(origin)), last_(runStartUtc), runStart_(runStartUtc), minRate_(minRateHz),
        windowS_(windowS) {}
  // Un bloque con eventos leído a tReadUtc.
  std::optional<hvmon::NotifyEvent> onTriggers(UtcNs tReadUtc, std::uint64_t nEvents, std::optional<double> boardTempDegC);
  // Llamar periódicamente (p. ej. cada segundo).
  std::optional<hvmon::NotifyEvent> tick(UtcNs now, std::optional<double> boardTempDegC);
  bool alerted() const { return alerted_; }
  UtcNs lastTriggerUtc() const { return last_; }

 private:
  hvmon::NotifyEvent make(const std::string& kind, UtcNs t, std::optional<double> temp) const;
  double windowRateHz(UtcNs now) const;
  bool lowRate(UtcNs now) const;
  double thresholdS_;
  std::string origin_;
  UtcNs last_, runStart_;
  double minRate_, windowS_;
  std::deque<std::pair<UtcNs, std::uint64_t>> recent_;  // (hora de lectura, eventos) en la ventana
  bool any_ = false, alerted_ = false;
  UtcNs episodeStart_ = 0;
  std::uint64_t episode_ = 0;
};

// Flags de evento que no dependen de HV: contador de la placa, memoria llena, TTT contra PC,
// huecos y ráfagas con la tasa de la última hora (trigger_interval_stats.h).
class EventFlagger {
 public:
  EventFlagger(double gapProbability, double burstProbability, double tttPcToleranceS, int counterBits = 24);
  std::uint32_t flags(UtcNs tUtc, UtcNs tRead, std::uint32_t counterHw, bool boardFull, std::int64_t dtPrevNs);
  double rateHz() const;

 private:
  double gapP_, burstP_, tolS_;
  std::uint32_t counterMask_;
  bool haveCounter_ = false;
  std::uint32_t lastCounter_ = 0;
  std::deque<UtcNs> recent_;  // triggers de la última hora
};

// Desenrolla el TTT con la hora del PC de cada lectura: elige el número de vueltas que deja
// origen + ttt_ns más cerca de t_read, sin retroceder respecto del evento anterior. Solo se usan los
// wrap_bits bits bajos (VX1740: 31 bit de 8 ns; el bit 31 es una marca de vuelta, DAQ2-6).
class TttUnwrapper {
 public:
  TttUnwrapper(UtcNs originUtc, double tickNs, int wrapBits);
  std::int64_t tttNs(std::uint32_t raw, UtcNs tReadUtc);

 private:
  UtcNs origin_;
  ttt::Params p_;  // tick y bits del contador (daq/ttt.h, compartido con el análisis)
  std::int64_t last_ = -1;
};

// Reloj del TTT ajustado al PC (DAQ2-9). Con un tick fijo, t_utc se separa de la hora del PC a
// razón del error del tick (run 8: 8.000071 frente a 8 ns, ~0.77 s por día). Aquí:
//  - u = TTT desenrollado en ns de tick nominal; la vuelta se elige con la hora de lectura pasada a
//    la escala de u con el modelo vigente (no se pierde aunque el run dure días);
//  - ajuste por mínimos cuadrados y = c + s·u (y = t_read − origen) sobre las lecturas de los
//    últimos window_s; s = tick real / nominal; c absorbe la latencia media de lectura;
//  - re-anclaje cada reanchor_s (hora de lectura) y cuando lo pide el escritor (segmento nuevo):
//    el modelo t_utc = t_ancla + (u − u_ancla)·s sigue continuo en el ancla, salvo en el primer
//    ajuste válido, que ancla en el origen del run (t = origen + s·u) y deja registrado el salto
//    (corrige el tick provisional de los primeros minutos);
//  - con menos de min_points lecturas o un tramo menor que min_span_s se mantiene el tick vigente.
// Determinista y sin hilos: pertenece al hilo escritor.
class TttClock {
 public:
  struct Params {
    double nominalTickNs = 8.0;   // escala de u (unwrap)
    double initialTickNs = 8.0;   // tick provisional hasta el primer ajuste (digitizer.ttt_tick_ns)
    int wrapBits = 31;
    bool fit = true;              // false: tick fijo (comportamiento anterior)
    double windowS = 3600, reanchorS = 600, minSpanS = 300;
    std::size_t minPoints = 100;
  };
  struct Anchor {
    UtcNs tUtcNs = 0;       // t_utc del ancla
    std::int64_t uNs = 0;   // u del ancla
    double tickNs = 0;      // tick usado desde el ancla
    std::size_t nFit = 0;   // lecturas del ajuste (0: tick provisional)
    double spanS = 0, residualRmsMs = 0;
    std::int64_t jumpNs = 0;  // salto de t_utc en el ancla (solo el primer ajuste)
    std::string reason;     // "inicio", "periódico", "segmento"
    nlohmann::json toJson() const;
  };
  TttClock(UtcNs originUtc, const Params& p);
  // Un evento, en orden: ttt_ns (= t_utc − origen) y t_utc.
  std::pair<std::int64_t, UtcNs> map(std::uint32_t raw, UtcNs tReadUtc);
  void requestReanchor(const std::string& reason) { pendingReason_ = reason; }
  const std::vector<Anchor>& anchors() const { return anchors_; }
  double tickNs() const { return anchors_.back().tickNs; }
  nlohmann::json paramsJson() const;

 private:
  void reanchor(std::int64_t u, UtcNs tRead, const std::string& reason);
  UtcNs model(std::int64_t u) const;
  UtcNs origin_;
  Params p_;
  std::int64_t lastU_ = -1;
  std::deque<std::pair<std::int64_t, std::int64_t>> win_;  // (u, t_read − origen)
  std::vector<Anchor> anchors_;
  UtcNs lastAnchorRead_ = 0;
  bool fitted_ = false;
  std::string pendingReason_;
};

}  // namespace daq::daqrun
