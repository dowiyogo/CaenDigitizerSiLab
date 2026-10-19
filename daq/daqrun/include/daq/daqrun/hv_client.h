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
  double tick_;
  std::uint32_t mask_;  // bits del contador (los demás se ignoran)
  std::int64_t period_;
  std::int64_t last_ = -1;
};

}  // namespace daq::daqrun
