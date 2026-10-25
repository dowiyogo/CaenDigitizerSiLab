// retry.cpp
#include "daq/hvmon/retry.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>

namespace daq::hvmon {

namespace {

constexpr UtcNs kNsPerMin = 60LL * 1000000000LL;
UtcNs minutes(double m) { return static_cast<UtcNs>(m * static_cast<double>(kNsPerMin)); }

bool contains(const std::vector<int>& v, int x) { return std::find(v.begin(), v.end(), x) != v.end(); }

RetryState stateFromString(const std::string& s) {
  for (RetryState x : {RetryState::Disabled, RetryState::Normal, RetryState::Wait, RetryState::Retry,
                       RetryState::Stabilizing, RetryState::GaveUp, RetryState::Manual})
    if (s == toString(x)) return x;
  return RetryState::Uninit;
}

}  // namespace

const char* toString(RetryState s) {
  switch (s) {
    case RetryState::Uninit: return "UNINIT";
    case RetryState::Disabled: return "DISABLED";
    case RetryState::Normal: return "NORMAL";
    case RetryState::Wait: return "WAIT";
    case RetryState::Retry: return "RETRY";
    case RetryState::Stabilizing: return "STABILIZING";
    case RetryState::GaveUp: return "GAVE_UP";
    case RetryState::Manual: return "MANUAL";
  }
  return "?";
}

std::string defaultRetryStateFile() {
  const char* x = std::getenv("XDG_STATE_HOME");
  if (x && *x) return std::string(x) + "/daqv2/hvmon_retry_state.json";
  const char* home = std::getenv("HOME");
  return std::string(home ? home : "") + "/.local/state/daqv2/hvmon_retry_state.json";
}

void RetryOutput::append(RetryOutput&& o) {
  for (auto& a : o.actions) actions.push_back(a);
  for (auto& t : o.transitions) transitions.push_back(std::move(t));
  for (auto& e : o.events) events.push_back(std::move(e));
  for (auto& w : o.warnings) warnings.push_back(std::move(w));
  persist = persist || o.persist;
}

Result<RetryConfig> RetryConfig::fromJson(const nlohmann::json& j) {
  using R = Result<RetryConfig>;
  RetryConfig c;
  c.stateFile = defaultRetryStateFile();
  try {
    c.enabledChannels = j.value("enabled_channels", c.enabledChannels);
    c.waitMin = j.value("wait_min", c.waitMin);
    c.stableMin = j.value("stable_min", c.stableMin);
    c.resetAfterH = j.value("reset_after_h", c.resetAfterH);
    c.imonStableMaxA = j.value("imon_stable_max_uA", c.imonStableMaxA * 1e6) / 1e6;
    c.maxRampingChannels = j.value("max_ramping_channels", c.maxRampingChannels);
    c.inhibitWhenGecoRunning = j.value("inhibit_when_geco_running", c.inhibitWhenGecoRunning);
    c.clearAlarmBeforeRetry = j.value("clear_alarm_before_retry", c.clearAlarmBeforeRetry);
    c.testModeAllowDetectorChannels = j.value("test_mode_allow_detector_channels", c.testModeAllowDetectorChannels);
    c.stateFile = j.value("state_file", c.stateFile);
    if (j.contains("test_mode")) {
      const auto& t = j["test_mode"];
      c.testModeEnabled = t.value("enabled", c.testModeEnabled);
      c.testChannels = t.value("channels", c.testChannels);
      c.testWaitMin = t.value("wait_min", c.testWaitMin);
      c.testStableMin = t.value("stable_min", c.testStableMin);
    }
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("auto_retry inválido: ") + e.what());
  }
  if (c.waitMin.empty() || c.testWaitMin.empty()) return R::failure(-1, "auto_retry.wait_min vacío");
  for (double w : c.waitMin)
    if (w <= 0) return R::failure(-1, "auto_retry.wait_min debe ser positivo");
  for (double w : c.testWaitMin)
    if (w <= 0) return R::failure(-1, "auto_retry.test_mode.wait_min debe ser positivo");
  if (c.stableMin <= 0 || c.testStableMin <= 0 || c.resetAfterH <= 0 || c.imonStableMaxA <= 0 || c.maxRampingChannels < 1)
    return R::failure(-1, "auto_retry: stable_min, reset_after_h, imon_stable_max_uA y max_ramping_channels deben ser positivos");
  if (c.testModeEnabled) {
    for (int ch : c.testChannels) {
      if (ch >= 0 && ch <= 7 && !c.testModeAllowDetectorChannels)
        return R::failure(-1, "auto_retry.test_mode incluye CH0" + std::to_string(ch) +
                                  " (canal del detector) sin test_mode_allow_detector_channels = true");
      if (!contains(c.enabledChannels, ch))
        return R::failure(-1, "auto_retry.test_mode: CH" + std::to_string(ch) + " no está en enabled_channels");
    }
  }
  return c;
}

nlohmann::json RetryConfig::toJson() const {
  return nlohmann::json{{"enabled_channels", enabledChannels},
                        {"wait_min", waitMin},
                        {"stable_min", stableMin},
                        {"reset_after_h", resetAfterH},
                        {"imon_stable_max_uA", imonStableMaxA * 1e6},
                        {"max_ramping_channels", maxRampingChannels},
                        {"inhibit_when_geco_running", inhibitWhenGecoRunning},
                        {"clear_alarm_before_retry", clearAlarmBeforeRetry},
                        {"test_mode_allow_detector_channels", testModeAllowDetectorChannels},
                        {"state_file", stateFile},
                        {"test_mode",
                         {{"enabled", testModeEnabled}, {"channels", testChannels}, {"wait_min", testWaitMin},
                          {"stable_min", testStableMin}}}};
}

bool RetryConfig::enabled(int ch) const { return contains(enabledChannels, ch); }
bool RetryConfig::inTestMode(int ch) const { return testModeEnabled && contains(testChannels, ch); }
const std::vector<double>& RetryConfig::waits(int ch) const { return inTestMode(ch) ? testWaitMin : waitMin; }
double RetryConfig::stableMinFor(int ch) const { return inTestMode(ch) ? testStableMin : stableMin; }

RetryMachine::RetryMachine(RetryConfig cfg, hv::StatusBitMap bits, std::vector<std::string> tripFlags,
                           std::vector<int> channels)
    : cfg_(std::move(cfg)), bits_(std::move(bits)), tripFlags_(std::move(tripFlags)), channels_(std::move(channels)) {
  for (int ch : channels_) ch_[ch];
}

void RetryMachine::restore(const nlohmann::json& p) {
  if (!p.is_object() || !p.contains("channels")) return;
  for (auto it = p["channels"].begin(); it != p["channels"].end(); ++it) {
    const int ch = std::atoi(it.key().c_str());
    auto f = ch_.find(ch);
    if (f == ch_.end()) continue;
    const auto& j = it.value();
    f->second.pState = stateFromString(j.value("state", ""));
    f->second.pK = j.value("k", 0);
    f->second.pTripUtc = j.value("trip_utc_ns", static_cast<UtcNs>(0));
    f->second.pEpisode = j.value("episode", static_cast<std::uint64_t>(0));
    f->second.pTimeline = j.value("timeline", nlohmann::json::array());
    nextEpisode_ = std::max(nextEpisode_, f->second.pEpisode + 1);
  }
}

std::string RetryMachine::label(RetryState s, int k) {
  if (s == RetryState::Wait || s == RetryState::Retry || s == RetryState::Stabilizing)
    return std::string(toString(s)) + "(" + std::to_string(k) + ")";
  return toString(s);
}

bool RetryMachine::hasTrip(std::uint32_t status) const {
  return std::any_of(tripFlags_.begin(), tripFlags_.end(), [&](const std::string& f) { return bits_.has(status, f); });
}

bool RetryMachine::ramping(std::uint32_t status) const {
  return bits_.has(status, "ramp_up") || bits_.has(status, "ramp_down");
}

std::uint32_t RetryMachine::maskRamp(std::uint32_t status) const {
  std::uint32_t m = 0;
  for (const char* n : {"ramp_up", "ramp_down"}) {
    const int b = bits_.bitOf(n);
    if (b >= 0) m |= 1u << b;
  }
  return status & ~m;
}

bool RetryMachine::stable(const RetryChannelView& v, double offsetA) const {
  return v.pw == 1u && bits_.has(v.status, "on") && !hasTrip(v.status) && !ramping(v.status) &&
         (v.imonA - offsetA) < cfg_.imonStableMaxA;
}

void RetryMachine::go(Ch& c, int ch, RetryState to, const std::string& reason, UtcNs now, RetryOutput& out,
                      const std::string& fromLabel) {
  RetryTransition t;
  t.ch = ch;
  t.from = fromLabel.empty() ? label(c.st, c.k) : fromLabel;
  c.st = to;
  c.sinceUtc = now;
  c.blockReason.clear();
  t.to = label(c.st, c.k);
  t.k = c.k;
  t.reason = reason;
  t.tUtcNs = now;
  out.transitions.push_back(t);
  out.persist = true;
}

void RetryMachine::event(Ch& c, int ch, const std::string& kind, const std::string& reason, UtcNs now, RetryOutput& out) {
  RetryEvent e;
  e.kind = kind;
  e.ch = ch;
  e.episode = c.episode;
  e.k = c.k;
  e.reason = reason;
  e.timeline = c.timeline;
  e.tUtcNs = now;
  out.events.push_back(e);
}

void RetryMachine::enterWait(Ch& c, int ch, const RetryChannelView& v, int k, const std::string& reason, UtcNs now,
                             RetryOutput& out) {
  const std::string from = label(c.st, c.k);
  c.k = k;
  c.tripUtc = now;
  c.waitStatus = maskRamp(v.status);
  c.waitV0 = v.v0setV;
  c.gecoWarned = false;
  go(c, ch, RetryState::Wait, reason, now, out, from);
}

void RetryMachine::failAttempt(Ch& c, int ch, const RetryChannelView& v, const std::string& reason, UtcNs now,
                               RetryOutput& out, bool trip) {
  c.timeline.push_back({{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", reason}, {"k", c.k}});
  if (trip) event(c, ch, "trip", reason, now, out);
  const int n = static_cast<int>(cfg_.waits(ch).size());
  if (c.k >= n) {
    c.lastOutcome = "sin reintentos: " + reason;
    c.timeline.push_back({{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", "GAVE_UP"}, {"k", c.k}});
    go(c, ch, RetryState::GaveUp, reason + "; no quedan reintentos", now, out);
    event(c, ch, "gave_up", reason, now, out);
  } else {
    c.lastOutcome = "intento " + std::to_string(c.k) + " fallido: " + reason;
    enterWait(c, ch, v, c.k + 1, reason, now, out);
  }
}

void RetryMachine::initial(Ch& c, int ch, const RetryChannelView& v, UtcNs now, RetryOutput& out) {
  if (!cfg_.enabled(ch)) {
    // Canal sin reintento: queda DISABLED sin registro de cambio (solo se registran sus cambios
    // normales); no hay transición que persistir.
    c.st = RetryState::Disabled;
    c.sinceUtc = now;
    return;
  }
  const auto& waits = cfg_.waits(ch);
  if (c.pState == RetryState::Wait && v.pw == 0u && c.pK >= 1 && c.pK <= static_cast<int>(waits.size())) {
    const UtcNs limit = c.pTripUtc + minutes(waits[static_cast<std::size_t>(c.pK - 1)] + 10.0);
    if (now < limit) {
      c.episode = c.pEpisode;
      c.timeline = c.pTimeline.is_array() ? c.pTimeline : nlohmann::json::array();
      c.k = c.pK;
      c.tripUtc = c.pTripUtc;
      c.waitStatus = maskRamp(v.status);
      c.waitV0 = v.v0setV;
      go(c, ch, RetryState::Wait, "reanudado tras reiniciar hvmon (trip a las " + formatUtcIso(c.tripUtc) + ")", now, out);
      return;
    }
    go(c, ch, RetryState::Manual, "espera persistida vencida hace más de 10 min: sin acción automática", now, out);
    return;
  }
  if (c.pState == RetryState::GaveUp && v.pw == 0u) {
    c.episode = c.pEpisode;
    c.k = c.pK;
    go(c, ch, RetryState::GaveUp, "reanudado tras reiniciar hvmon (sigue apagado)", now, out);
    return;
  }
  if (v.pw == 1u) {
    c.k = (c.pState == RetryState::Normal || c.pState == RetryState::Stabilizing) ? c.pK : 0;
    c.normalSinceUtc = now;
    go(c, ch, RetryState::Normal, "encendido al arrancar", now, out);
  } else {
    go(c, ch, RetryState::Manual, "apagado al arrancar: sin acción automática", now, out);
  }
}

RetryOutput RetryMachine::update(const RetryInput& in) {
  RetryOutput out;
  const UtcNs now = in.nowUtc;
  // Canales en rampa (bits de Status) o en RETRY, para max_ramping_channels.
  auto rampingOthers = [&](int self) {
    int n = 0;
    for (const auto& kv : ch_) {
      if (kv.first == self) continue;
      auto v = in.channels.find(kv.first);
      const bool r = kv.second.st == RetryState::Retry || (v != in.channels.end() && ramping(v->second.status));
      if (r) ++n;
    }
    return n;
  };
  for (int ch : channels_) {
    auto vit = in.channels.find(ch);
    if (vit == in.channels.end()) continue;
    const RetryChannelView& v = vit->second;
    Ch& c = ch_[ch];
    if (c.st == RetryState::Uninit) {
      initial(c, ch, v, now, out);
      continue;
    }
    std::string manualWhy;
    bool tripChange = false;
    for (const auto& r : in.changes) {
      if (r.ch != ch) continue;
      if (r.origin == "external" && (r.param == "Pw" || r.param == "V0Set"))
        manualWhy = r.param + " " + r.oldValue.dump() + " -> " + r.newValue.dump();
      if (r.origin == "trip" && r.param == "Pw" && r.newValue == 0) tripChange = true;
    }
    if (c.st == RetryState::Disabled) {
      // XPLUS-1: el aviso de trip no depende de auto_retry.enabled_channels. Un canal sin reintento no escribe nada, pero su
      // trip genera un evento "trip" (correo si notify_each_trip), una vez por trip: se libera cuando vuelve a estar encendido
      // sin bit de trip.
      const bool newTripD = tripChange || (v.pw == 1u && hasTrip(v.status));
      if (newTripD && !c.tripLatched) {
        c.tripLatched = true;
        c.episode = nextEpisode_++;
        c.timeline = nlohmann::json::array({{{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", "trip"}, {"k", 0}}});
        event(c, ch, "trip", "trip (canal sin reintento automático)", now, out);
      } else if (c.tripLatched && v.pw == 1u && !hasTrip(v.status)) {
        c.tripLatched = false;
      }
      continue;
    }
    const bool tripNow = hasTrip(v.status) || tripChange;
    // En RETRY/STABILIZING solo cuenta un trip nuevo: Pw 1 -> 0 con origen trip, o bit de trip con
    // Pw = 1. Un bit retenido con el canal apagado (antes de que la relectura muestre Pw = 1) no.
    const bool newTrip = tripChange || (v.pw == 1u && hasTrip(v.status));
    const double offset = in.imonOffsetA.count(ch) ? in.imonOffsetA.at(ch) : 0.0;
    if (!manualWhy.empty() && c.st != RetryState::Manual) {
      c.stableSinceUtc = 0;
      c.lastOutcome = "acción manual";
      go(c, ch, RetryState::Manual, "acción manual (" + manualWhy + ")", now, out);
      continue;
    }
    switch (c.st) {
      case RetryState::Normal:
        if (tripNow) {
          if (c.k == 0) {
            c.episode = nextEpisode_++;
            c.timeline = nlohmann::json::array();
          }
          c.timeline.push_back({{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", "trip"}, {"k", c.k}});
          event(c, ch, "trip", "trip en NORMAL", now, out);
          if (c.k >= static_cast<int>(cfg_.waits(ch).size())) {
            c.lastOutcome = "trip sin reintentos disponibles";
            go(c, ch, RetryState::GaveUp, "trip y no quedan reintentos (contador " + std::to_string(c.k) + ")", now, out);
            event(c, ch, "gave_up", "trip y no quedan reintentos", now, out);
          } else {
            enterWait(c, ch, v, c.k + 1, "trip", now, out);
          }
        } else if (c.k > 0 && now - c.normalSinceUtc >= static_cast<UtcNs>(cfg_.resetAfterH * 60.0) * kNsPerMin) {
          c.k = 0;
          go(c, ch, RetryState::Normal, "contador de intentos reiniciado tras " + std::to_string(cfg_.resetAfterH) + " h en NORMAL",
             now, out);
        }
        break;
      case RetryState::Wait: {
        if (v.pw == 1u) {
          go(c, ch, RetryState::Manual, "el canal se encendió durante la espera", now, out);
          break;
        }
        if (maskRamp(v.status) != c.waitStatus || v.v0setV != c.waitV0) {
          go(c, ch, RetryState::Manual, "el estado del canal cambió durante la espera", now, out);
          break;
        }
        const auto& waits = cfg_.waits(ch);
        if (now - c.tripUtc < minutes(waits[static_cast<std::size_t>(c.k - 1)])) break;
        std::string block;
        if (!in.fresh)
          block = "sin lectura fresca (conexión)";
        else if (in.geco && cfg_.inhibitWhenGecoRunning)
          block = "GECO2020 en ejecución (inhibido por configuración)";
        else if (rampingOthers(ch) >= cfg_.maxRampingChannels)
          block = "otro canal en rampa (max_ramping_channels = " + std::to_string(cfg_.maxRampingChannels) + ")";
        if (!block.empty()) {
          c.blockReason = block;
          break;
        }
        if (in.geco && !c.gecoWarned) {
          out.warnings.push_back("GECO2020 en ejecución durante el reintento de CH" + std::to_string(ch) +
                                 " (no inhibe: inhibit_when_geco_running = false)");
          c.gecoWarned = true;
        }
        RetryAction a;
        a.ch = ch;
        a.clearAlarmFirst = cfg_.clearAlarmBeforeRetry;
        a.token = c.token = nextToken_++;
        a.k = c.k;
        out.actions.push_back(a);
        c.timeline.push_back({{"t_utc_ns", now},
                              {"t_utc", formatUtcIso(now)},
                              {"event", "retry"},
                              {"k", c.k},
                              {"wait_min", waits[static_cast<std::size_t>(c.k - 1)]}});
        go(c, ch, RetryState::Retry, "reintento " + std::to_string(c.k) + (a.clearAlarmFirst ? " (ClrAlarm y Pw = 1)" : " (Pw = 1)"),
           now, out);
        break;
      }
      case RetryState::Retry:
        if (newTrip) failAttempt(c, ch, v, "trip durante la rampa del intento " + std::to_string(c.k), now, out, true);
        break;
      case RetryState::Stabilizing:
        if (newTrip) {
          failAttempt(c, ch, v, "trip durante la estabilización del intento " + std::to_string(c.k), now, out, true);
        } else if (stable(v, offset)) {
          if (c.stableSinceUtc == 0) c.stableSinceUtc = now;
          if (now - c.stableSinceUtc >= minutes(cfg_.stableMinFor(ch))) {
            c.lastOutcome = "recuperado en el intento " + std::to_string(c.k);
            c.normalSinceUtc = now;
            c.timeline.push_back({{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", "recovered"}, {"k", c.k}});
            go(c, ch, RetryState::Normal, c.lastOutcome, now, out);
            event(c, ch, "recovery", c.lastOutcome, now, out);
          }
        } else {
          c.stableSinceUtc = 0;
        }
        break;
      case RetryState::Manual:
        if (stable(v, offset)) {
          if (c.stableSinceUtc == 0) c.stableSinceUtc = now;
          if (now - c.stableSinceUtc >= minutes(cfg_.stableMinFor(ch))) {
            c.k = 0;
            c.normalSinceUtc = now;
            c.lastOutcome = "estable tras acción manual";
            go(c, ch, RetryState::Normal, c.lastOutcome, now, out);
          }
        } else {
          c.stableSinceUtc = 0;
        }
        break;
      case RetryState::GaveUp:
      case RetryState::Disabled:
      case RetryState::Uninit:
        break;
    }
  }
  return out;
}

RetryOutput RetryMachine::onCommandResult(int ch, std::uint64_t token, const std::string& phase,
                                          const std::string& detail, UtcNs now) {
  RetryOutput out;
  auto it = ch_.find(ch);
  if (it == ch_.end()) return out;
  Ch& c = it->second;
  if (c.st != RetryState::Retry || token != c.token) return out;  // resultado de otro intento
  RetryChannelView v;
  v.status = c.waitStatus;
  v.v0setV = c.waitV0;
  if (phase == "verified" || phase == "noop") {
    c.stableSinceUtc = 0;
    c.timeline.push_back({{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", "ramp_done"}, {"k", c.k}});
    go(c, ch, RetryState::Stabilizing, "rampa completa (" + detail + ")", now, out);
  } else if (phase == "aborted") {
    c.lastOutcome = "intento abortado: " + detail;
    go(c, ch, RetryState::Manual, "lectura previa distinta de la publicada: " + detail, now, out);
  } else if (phase == "busy" || phase == "write_failed") {
    go(c, ch, RetryState::Wait, "escritura pospuesta (" + phase + "): " + detail, now, out);
  } else if (phase == "verify_failed") {
    failAttempt(c, ch, v, "la rampa no terminó a tiempo (" + detail + ")", now, out, false);
  } else {  // rejected u otro
    c.lastOutcome = "escritura rechazada: " + detail;
    go(c, ch, RetryState::Manual, "escritura rechazada (" + phase + "): " + detail, now, out);
  }
  return out;
}

nlohmann::json RetryMachine::publish(UtcNs now) const {
  nlohmann::json arr = nlohmann::json::array();
  for (int ch : channels_) {
    const Ch& c = ch_.at(ch);
    nlohmann::json j{{"ch", ch},
                     {"state", toString(c.st)},
                     {"k", c.k},
                     {"label", label(c.st, c.k)},
                     {"enabled", cfg_.enabled(ch)},
                     {"test_mode", cfg_.inTestMode(ch)},
                     {"episode", c.episode},
                     {"last_outcome", c.lastOutcome},
                     {"block_reason", c.blockReason},
                     {"since_utc_ns", c.sinceUtc}};
    j["trip_utc_ns"] = c.tripUtc ? nlohmann::json(c.tripUtc) : nlohmann::json();
    j["trip_utc"] = c.tripUtc ? nlohmann::json(formatUtcIso(c.tripUtc)) : nlohmann::json();
    if (c.st == RetryState::Wait && c.k >= 1) {
      const auto& waits = cfg_.waits(ch);
      const double left = static_cast<double>(c.tripUtc + minutes(waits[static_cast<std::size_t>(c.k - 1)]) - now) / 1e9;
      j["next_retry_in_s"] = std::max(0.0, left);
    } else {
      j["next_retry_in_s"] = nlohmann::json();
    }
    j["stable_for_s"] = c.stableSinceUtc ? nlohmann::json(static_cast<double>(now - c.stableSinceUtc) / 1e9) : nlohmann::json();
    arr.push_back(j);
  }
  return arr;
}

nlohmann::json RetryMachine::persistJson(UtcNs now) const {
  nlohmann::json chans = nlohmann::json::object();
  for (const auto& kv : ch_) {
    const Ch& c = kv.second;
    chans[std::to_string(kv.first)] = {{"state", toString(c.st)}, {"k", c.k}, {"trip_utc_ns", c.tripUtc},
                                       {"episode", c.episode}, {"last_outcome", c.lastOutcome},
                                       {"since_utc_ns", c.sinceUtc}, {"timeline", c.timeline}};
  }
  return nlohmann::json{{"schema", "daqv2.hvmon_retry_state/1"}, {"saved_utc_ns", now},
                        {"saved_utc", formatUtcIso(now)}, {"channels", chans}};
}

RetryState RetryMachine::state(int ch) const {
  auto it = ch_.find(ch);
  return it == ch_.end() ? RetryState::Uninit : it->second.st;
}

int RetryMachine::attempt(int ch) const {
  auto it = ch_.find(ch);
  return it == ch_.end() ? 0 : it->second.k;
}

}  // namespace daq::hvmon
