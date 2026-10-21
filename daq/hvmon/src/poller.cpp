// poller.cpp
#include "daq/hvmon/poller.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "daq/core/logger.h"

namespace daq::hvmon {

namespace {

std::string describeValue(const hv::ParamValue& v) {
  std::ostringstream os;
  if (v.isNumeric())
    os << v.si;
  else
    os << v.raw;
  return os.str();
}

bool contains(const std::vector<int>& v, int x) { return std::find(v.begin(), v.end(), x) != v.end(); }
bool containsStr(const std::vector<std::string>& v, const std::string& x) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

}  // namespace

Poller::Poller(const HvmonConfig& cfg, hv::HVBackend& backend, Clock& clock, BoundedQueue<WriterItem>& writer,
               StateHub& hub, BoundedQueue<PollRequest>& requests)
    : cfg_(cfg),
      backend_(backend),
      clock_(clock),
      writer_(writer),
      hub_(hub),
      requests_(requests),
      detector_(cfg.statusBits, cfg.tripFlags),
      backoffS_(cfg.backoffInitialS),
      nextConnect_(clock.now_steady()),
      retry_(cfg.retry, cfg.statusBits, cfg.tripFlags, cfg.channels) {
  state_.activeCommand = nlohmann::json();
  state_.activeOffsets = nlohmann::json();
  state_.lastOffsets = nlohmann::json();
}

Poller::~Poller() = default;

std::chrono::nanoseconds Poller::secondsToNs(double s) const {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(s));
}

void Poller::run(const StopToken& st) {
  Logger::setThreadName("poll");
  while (!st.stop_requested()) step(st);
  shutdown();
}

SteadyTime Poller::nextWake(SteadyTime now) const {
  SteadyTime t = connected_ ? std::min(nextFast_, nextFull_) : nextConnect_;
  if (active_) t = std::min(t, active_->deadline);
  const SteadyTime cap = now + std::chrono::milliseconds(200);
  return std::min(std::max(t, now), cap);
}

void Poller::step(const StopToken& st) {
  const SteadyTime now = clock_.now_steady();
  while (auto r = requests_.try_pop()) handleRequest(*r, now);
  if (!connected_ && now >= nextConnect_) tryConnect(now);
  if (connected_ && now >= nextFull_) {
    fullRead(now);
    nextFull_ += secondsToNs(cfg_.fullPeriodS);
    if (nextFull_ <= now) nextFull_ = now + secondsToNs(cfg_.fullPeriodS);
  }
  if (connected_ && now >= nextFast_) {
    fastRead(now);
    nextFast_ += secondsToNs(cfg_.fastPeriodS);
    if (nextFast_ <= now) nextFast_ = now + secondsToNs(cfg_.fastPeriodS);
  }
  if (active_ && now >= active_->deadline)
    finishCommand("verify_failed", "timeout de verificación sin alcanzar el estado pedido", false, now);
  if (autoRestorePending_ && connected_ && prevFast_ && !active_ && !restoring_) {
    autoRestorePending_ = false;
    if (!cfg_.allowCommands) {
      logError("sesión de barrido activa pero allow_commands = false: el V0Set NO se puede restaurar automáticamente");
    } else {
      PollRequest q;
      q.kind = PollRequest::Kind::BiasRestore;
      q.reason = "restauración automática del V0Set nominal al arrancar hvmon (sesión de barrido activa)";
      q.requester = "hvmon (arranque)";
      handleBiasRestore(q, now);
    }
  }
  if (restoring_) advanceRestore(now);
  detector_.expire(now);
  publish(now);
  SteadyTime wake = nextWake(now);
  if (wake <= now) wake = now + std::chrono::milliseconds(1);
  clock_.sleep_until(wake, st);
}

void Poller::tryConnect(SteadyTime now) {
  ++stats_.connectAttempts;
  if (attemptTimes_.size() < 10000) attemptTimes_.push_back(now);
  Status st = backend_.connect();
  if (st && allParams_.empty()) {
    auto p = backend_.channelParamNames(cfg_.slot, cfg_.channels.front());
    if (!p) {
      st = Status(p.error());
    } else {
      state_.params = p.value();
      for (const auto& info : p.value()) {
        allParams_.push_back(info.name);
        if (info.name != "VMon" && info.name != "IMon" && !containsStr(cfg_.fastTracked, info.name))
          fullTracked_.push_back(info.name);
      }
      for (const auto& f : cfg_.fastParams)
        if (!containsStr(allParams_, f)) st = Status::failure(-1, "la placa no tiene el parámetro " + f);
    }
    if (st) {
      // Parámetros de placa: opcionales (si la placa no los informa, sigue sin ellos).
      auto bp = backend_.boardParamNames(cfg_.slot);
      if (bp) {
        for (const auto& info : bp.value())
          if (!containsStr(cfg_.boardParamsExclude, info.name)) boardParams_.push_back(info.name);
      } else {
        logWarn("la placa no informa parámetros de placa: ", bp.error().message);
      }
    }
    if (!st) {
      allParams_.clear();
      fullTracked_.clear();
      boardParams_.clear();
      backend_.disconnect();
    }
  }
  if (!st) {
    ++stats_.connectFailures;
    if (st.error().message != state_.lastError) logWarn("conexión HV falló: ", st.error().message);
    state_.lastError = st.error().message;
    state_.lastErrorUtcNs = clock_.now_utc_ns();
    nextConnect_ = now + secondsToNs(backoffS_);
    backoffS_ = std::min(2 * backoffS_, cfg_.backoffMaxS);
    return;
  }
  connected_ = true;
  if (everConnected_) {
    ++stats_.reconnections;
    logInfo("reconectado a la fuente HV (reconexión ", stats_.reconnections, ")");
  } else {
    logInfo("conectado a la fuente HV: ", backend_.describe());
  }
  everConnected_ = true;
  backoffS_ = cfg_.backoffInitialS;
  nextFast_ = now;
  nextFull_ = now;
}

void Poller::onReadError(const Error& e, SteadyTime now) {
  ++stats_.readErrors;
  logWarn("lectura HV falló: ", e.message);
  state_.lastError = e.message;
  state_.lastErrorUtcNs = clock_.now_utc_ns();
  backend_.disconnect();
  connected_ = false;
  fastGap_ = fullGap_ = true;
  nextConnect_ = now + secondsToNs(backoffS_);
  backoffS_ = std::min(2 * backoffS_, cfg_.backoffMaxS);
  if (active_) finishCommand("verify_failed", "conexión perdida durante la verificación: " + e.message, false, now);
  if (offsets_) {
    offsets_->reply->set_value({{"ok", false}, {"error", "conexión perdida durante la medida de offsets"}});
    offsets_.reset();
  }
}

void Poller::pushWriter(WriterItem item) {
  if (!writer_.try_push(std::move(item), DropPolicy::DropNewest)) {
    const auto n = ++hub_.writerDropped;
    if (n == 1 || n % 1000 == 0) logError("cola del escritor llena: ", n, " datos descartados");
  }
}

void Poller::addChanges(std::vector<ChangeRecord>& changes) {
  for (auto& c : changes) {
    logInfo("cambio CH", c.ch, " ", c.param, ": ", c.oldValue.dump(), " -> ", c.newValue.dump(), " (", c.origin,
            c.duringDisconnect ? ", tras desconexión" : "", ")");
    pushWriter(ChangeItem{c});
    state_.recentChanges.push_back(c);
    while (static_cast<int>(state_.recentChanges.size()) > cfg_.changeRingSize) state_.recentChanges.pop_front();
  }
}

void Poller::fastRead(SteadyTime now) {
  auto r = backend_.readChannels(cfg_.slot, cfg_.channels, cfg_.fastParams);
  if (!r) {
    onReadError(r.error(), now);
    return;
  }
  const hv::ChannelReading& rd = r.value();
  ++readSeq_;
  ++stats_.fastReads;
  const double ms = static_cast<double>(rd.durationNs) / 1e6;
  stats_.sumFastMs += ms;
  stats_.maxFastMs = std::max(stats_.maxFastMs, ms);
  std::vector<ChangeRecord> changes;
  if (prevFast_) {
    changes = detector_.compare(*prevFast_, rd, cfg_.fastTracked, "fast", fastGap_, now, nextChangeId_);
    addChanges(changes);
  }
  fastGap_ = false;
  pushWriter(FastRows{rd, readSeq_});
  for (std::size_t p = 0; p < rd.params.size(); ++p)
    for (std::size_t c = 0; c < rd.channels.size(); ++c) state_.values[rd.channels[c]][rd.params[p].name] = rd.values[p][c];
  state_.lastGoodReadUtcNs = rd.tRespUtcNs;
  state_.lastReadReqUtcNs = rd.tReqUtcNs;
  state_.lastReadRespUtcNs = rd.tRespUtcNs;
  state_.lastReadSeq = readSeq_;
  state_.lastReadKind = "fast";
  for (std::size_t p = 0; p < rd.params.size(); ++p) {
    const bool v = rd.params[p].name == "VMon", i = rd.params[p].name == "IMon";
    if (!v && !i) continue;
    for (std::size_t c = 0; c < rd.channels.size(); ++c) {
      auto& w = window_.try_emplace(rd.channels[c], WindowStats(state_.windowS), WindowStats(state_.windowS)).first->second;
      (v ? w.first : w.second).add(rd.tRespUtcNs, rd.values[p][c].si);
    }
  }
  verifyCommand(rd, now);
  accumulateOffsets(rd, now);
  prevFast_ = rd;
  runRetry(rd, changes, now);
}

void Poller::fullRead(SteadyTime now) {
  auto r = backend_.readChannels(cfg_.slot, cfg_.channels, allParams_);
  if (!r) {
    onReadError(r.error(), now);
    return;
  }
  const hv::ChannelReading& rd = r.value();
  ++readSeq_;
  ++stats_.fullReads;
  const double ms = static_cast<double>(rd.durationNs) / 1e6;
  stats_.sumFullMs += ms;
  stats_.maxFullMs = std::max(stats_.maxFullMs, ms);
  if (prevFull_) {
    auto changes = detector_.compare(*prevFull_, rd, fullTracked_, "full", fullGap_, now, nextChangeId_);
    addChanges(changes);
  }
  fullGap_ = false;
  std::optional<hv::BoardReading> board;
  if (!boardParams_.empty()) {
    auto b = backend_.readBoard(cfg_.slot, boardParams_);
    if (!b) {
      onReadError(b.error(), now);
      return;
    }
    board = b.value();
    state_.board = board;
  }
  pushWriter(FullRows{rd, readSeq_, board});
  state_.params = rd.params;
  for (std::size_t p = 0; p < rd.params.size(); ++p)
    for (std::size_t c = 0; c < rd.channels.size(); ++c) state_.values[rd.channels[c]][rd.params[p].name] = rd.values[p][c];
  state_.lastGoodReadUtcNs = rd.tRespUtcNs;
  state_.lastReadReqUtcNs = rd.tReqUtcNs;
  state_.lastReadRespUtcNs = rd.tRespUtcNs;
  state_.lastReadSeq = readSeq_;
  state_.lastReadKind = "full";
  prevFull_ = rd;
}

const hv::ParamValue* Poller::lastValue(int ch, const std::string& param) const {
  auto it = state_.values.find(ch);
  if (it == state_.values.end()) return nullptr;
  auto v = it->second.find(param);
  return v == it->second.end() ? nullptr : &v->second;
}

bool Poller::anyRamping() const {
  for (const auto& kv : state_.values) {
    auto s = kv.second.find("Status");
    if (s == kv.second.end()) continue;
    if (cfg_.statusBits.has(s->second.raw, "ramp_up") || cfg_.statusBits.has(s->second.raw, "ramp_down")) return true;
  }
  return false;
}

void Poller::record(const CommandRecord& rec) { pushWriter(CommandItem{rec}); }

void Poller::handleRequest(PollRequest& r, SteadyTime now) {
  switch (r.kind) {
    case PollRequest::Kind::Pw: handlePw(r, now); break;
    case PollRequest::Kind::Offsets: handleOffsets(r, now); break;
    case PollRequest::Kind::BiasBegin: handleBiasBegin(r); break;
    case PollRequest::Kind::V0Set: handleV0Set(r, now); break;
    case PollRequest::Kind::BiasRestore: handleBiasRestore(r, now); break;
    case PollRequest::Kind::BiasStatus:
      r.reply->set_value(nlohmann::json{{"ok", true}, {"session", bias_.toJson(clock_.now_utc_ns())}, {"restoring", restoring_},
                                        {"pending_steps", restoreSteps_.size()}});
      break;
  }
}

// ---- Barrido de tensión (DAQ2-14) ----------------------------------------------------------------

void Poller::persistBias() {
  if (!fs_ || biasPath_.empty()) return;
  fs_->makeDirectories(parentDirectory(biasPath_), 0700);
  if (auto st = fs_->writeAtomic(biasPath_, bias_.toJson(clock_.now_utc_ns()).dump(1) + "\n", 0600); !st)
    logError("no se pudo guardar el estado del barrido de tensión: ", st.error().message);
}

void Poller::enableBiasPersistence(FileSystem& fs, const std::string& path) {
  fs_ = &fs;
  biasPath_ = path;
  if (!fs.exists(path)) return;
  auto content = fs.readFile(path);
  if (!content) {
    logError("no se pudo leer el estado del barrido de tensión ", path, ": ", content.error().message);
    return;
  }
  try {
    if (auto st = bias_.restore(nlohmann::json::parse(content.value())); !st) {
      logError(st.error().message, " (", path, ")");
      return;
    }
  } catch (const std::exception& e) {
    logError("estado del barrido de tensión ilegible en ", path, ": ", e.what());
    return;
  }
  if (bias_.active()) {
    logWarn("sesión de barrido de tensión activa al arrancar: se restaura el V0Set nominal al conectar");
    autoRestorePending_ = true;
  }
}

void Poller::handleBiasBegin(PollRequest& r) {
  auto reply = [&](bool ok, const std::string& detail) {
    if (!ok) logWarn("bias-begin rechazado: ", detail);
    r.reply->set_value(nlohmann::json{{"ok", ok}, {"result", ok ? "started" : "rejected"}, {"detail", detail},
                                      {"session", bias_.toJson(clock_.now_utc_ns())}});
  };
  if (!cfg_.allowCommands) return reply(false, "escrituras deshabilitadas (allow_commands = false)");
  if (r.reason.empty()) return reply(false, "falta authorize");
  if (!connected_ || !prevFast_) return reply(false, "sin conexión o sin lectura previa");
  if (active_ || restoring_) return reply(false, "hay un comando o una restauración en curso");
  std::map<int, double> v0, vm;
  for (int ch : bias_.limits().channels) {
    const auto* a = lastValue(ch, "V0Set");
    const auto* b = lastValue(ch, "VMon");
    const auto* pw = lastValue(ch, "Pw");
    const auto* st = lastValue(ch, "Status");
    if (!a || !b || !pw || !st) return reply(false, "sin lectura de CH" + std::to_string(ch));
    if (pw->raw != 1u || st->raw != 1u) return reply(false, "CH" + std::to_string(ch) + " no está encendido y estable (Status distinto de 0x1)");
    v0[ch] = a->si;
    vm[ch] = b->si;
  }
  if (auto st = bias_.begin(r.expected, v0, vm, r.reason, clock_.now_utc_ns()); !st) return reply(false, st.error().message);
  persistBias();
  CommandRecord rec;
  rec.id = nextCommandId_++;
  rec.slot = cfg_.slot;
  rec.ch = -1;
  rec.param = "bias_begin";
  rec.phase = "verified";
  rec.reason = r.reason;
  rec.requester = r.requester;
  rec.detail = bias_.toJson(clock_.now_utc_ns()).dump();
  rec.tUtcNs = clock_.now_utc_ns();
  record(rec);
  logInfo("sesión de barrido de tensión iniciada (", r.requester, "): ", r.reason);
  reply(true, "nominal registrado");
}

void Poller::handleV0Set(PollRequest& r, SteadyTime now) {
  CommandRecord rec;
  rec.id = nextCommandId_++;
  rec.slot = cfg_.slot;
  rec.ch = r.ch;
  rec.param = "V0Set";
  rec.valueV = r.value;
  rec.reason = r.reason;
  rec.requester = r.requester;
  rec.tUtcNs = clock_.now_utc_ns();
  rec.phase = "request";
  record(rec);
  auto finish = [&](const std::string& phase, const std::string& detail, bool ok) {
    rec.phase = phase;
    rec.detail = detail;
    rec.tUtcNs = clock_.now_utc_ns();
    record(rec);
    if (!ok) logWarn("comando ", rec.id, " V0Set CH", r.ch, "=", r.value, ": ", phase, ": ", detail);
    r.reply->set_value(nlohmann::json{{"ok", ok}, {"command_id", rec.id}, {"result", phase}, {"detail", detail}, {"ch", r.ch}, {"value_V", r.value}});
  };
  if (!cfg_.allowCommands) return finish("rejected", "escrituras deshabilitadas (allow_commands = false)", false);
  if (r.reason.empty()) return finish("rejected", "falta authorize", false);
  if (!contains(cfg_.channels, r.ch)) return finish("rejected", "canal fuera de la configuración", false);
  if (!connected_ || !prevFast_) return finish("rejected", "sin conexión o sin lectura previa", false);
  if (active_) return finish(r.internal ? "busy" : "rejected", "otro comando en curso (id " + std::to_string(active_->id) + ")", false);
  auto pre = backend_.readChannels(cfg_.slot, {r.ch}, cfg_.fastParams);
  if (!pre) {
    finish("aborted", "la relectura previa falló: " + pre.error().message, false);
    onReadError(pre.error(), now);
    return;
  }
  for (const char* p : {"Pw", "Status", "V0Set"}) {
    const hv::ParamValue* last = lastValue(r.ch, p);
    const hv::ParamValue* cur = pre.value().get(p, r.ch);
    const bool same = last && cur && (last->isNumeric() ? last->boardValue == cur->boardValue : last->raw == cur->raw);
    if (!same)
      return finish("aborted", std::string("el canal cambió entre la última lectura publicada y la relectura (") + p + ")", false);
  }
  const auto* pw = pre.value().get("Pw", r.ch);
  const auto* st = pre.value().get("Status", r.ch);
  const auto* v0 = pre.value().get("V0Set", r.ch);
  // Con el canal apagado solo se admite, durante una restauración, escribir el nominal (no hay rampa: al
  // encender, la placa sube de 0 al V0Set a RUp como siempre).
  const bool offWrite = pw->raw == 0u && r.internal && restoring_ && bias_.active() && bias_.inSession(r.ch) &&
                        std::fabs(r.value - bias_.nominal(r.ch)) <= 1e-6;
  if (!offWrite && (pw->raw != 1u || st->raw != 1u)) {
    char b[96];
    std::snprintf(b, sizeof b, "el canal no está encendido y estable (Pw %u, Status 0x%X)", pw->raw, st->raw);
    return finish("rejected", b, false);
  }
  if (!offWrite)
    if (auto ck = bias_.checkStep(r.ch, r.value, v0->si); !ck) return finish("rejected", ck.error().message, false);
  if (std::fabs(r.value - v0->si) <= bias_.limits().matchTolV) return finish("noop", "V0Set ya es el pedido", true);
  const hv::ParamValue* rate = lastValue(r.ch, r.value > v0->si ? "RUp" : "RDWn");
  const double rateVps = (rate && rate->si > 0) ? rate->si : 1.0;
  const double timeoutS = std::fabs(r.value - v0->si) / rateVps + cfg_.verifyMarginS + 30.0;
  detector_.addPending({rec.id, r.ch, now, now + secondsToNs(timeoutS + cfg_.commandWindowExtraS)});
  if (auto w = backend_.setChannelV0Set(cfg_.slot, r.ch, r.value); !w) {
    detector_.setPendingEnd(rec.id, now);
    return finish("write_failed", w.error().message, false);
  }
  rec.phase = "written";
  rec.detail = "timeout de verificación " + std::to_string(timeoutS) + " s";
  rec.tUtcNs = clock_.now_utc_ns();
  record(rec);
  logInfo("comando ", rec.id, ": V0Set CH", r.ch, "=", r.value, " V escrito (", r.requester, "): ", r.reason);
  ActiveCommand a{rec.id, r.ch, true, now + secondsToNs(timeoutS), r.reason, r.requester, r.reply, 0};
  a.v0set = true;
  a.offWrite = offWrite;
  a.target = r.value;
  a.v0setBefore = v0->si;
  if (const auto* vm = pre.value().get("VMon", r.ch)) a.vmonBefore = vm->si;
  active_ = a;
}

void Poller::handleBiasRestore(PollRequest& r, SteadyTime now) {
  if (!bias_.active()) {
    r.reply->set_value(nlohmann::json{{"ok", true}, {"result", "noop"}, {"detail", "sin sesión de barrido activa"}});
    return;
  }
  if (restoring_) {
    r.reply->set_value(nlohmann::json{{"ok", false}, {"result", "busy"}, {"detail", "ya hay una restauración en curso"}});
    return;
  }
  if (!cfg_.allowCommands) {
    r.reply->set_value(nlohmann::json{{"ok", false}, {"result", "rejected"}, {"detail", "escrituras deshabilitadas (allow_commands = false)"}});
    return;
  }
  restoreReason_ = r.reason.empty() ? std::string("restauración del nominal") : r.reason;
  restoreReply_ = r.reply;
  restoring_ = true;
  restoreFailed_ = false;
  restoreDetail_.clear();
  restoreSteps_.clear();
  restorePwTried_.clear();
  logWarn("restauración del V0Set nominal pedida (", r.requester, "): ", restoreReason_);
  advanceRestore(now);
}

void Poller::advanceRestore(SteadyTime now) {
  if (!restoring_ || active_ || !connected_ || !prevFast_) return;
  if (restoreFut_ && restoreFut_->wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
    const auto res = restoreFut_->get();
    restoreFut_.reset();
    if (!res.value("ok", false)) {
      restoreFailed_ = true;  // se sigue con los escalones (lo escrito ya está en la placa) y se informa al final
      restoreDetail_ += "CH" + std::to_string(res.value("ch", -1)) + " " + res.value("result", "") + ": " + res.value("detail", "") + "; ";
    }
  }
  // Plan fresco desde la última lectura: escalones de ≤ max_step_V hacia el nominal; luego Pw = 1.
  std::map<int, double> v0;
  for (int ch : bias_.limits().channels)
    if (const auto* a = lastValue(ch, "V0Set")) v0[ch] = a->si;
  for (int ch : bias_.limits().channels) {
    const auto* pw = lastValue(ch, "Pw");
    if (pw && pw->raw != 1u && v0.count(ch) && std::fabs(v0[ch] - bias_.nominal(ch)) > bias_.limits().matchTolV) {
      PollRequest q;
      q.kind = PollRequest::Kind::V0Set;
      q.ch = ch;
      q.value = bias_.nominal(ch);
      q.internal = true;
      q.reason = restoreReason_ + " (canal apagado: V0Set nominal antes de encender)";
      q.requester = "hvmon (restauración del barrido)";
      auto fut = q.reply->get_future().share();
      handleV0Set(q, now);
      if (fut.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
        restoreFut_ = fut;
        return;
      }
      const auto res = fut.get();
      restoreFailed_ = true;
      restoreDetail_ += "CH" + std::to_string(ch) + ": " + res.value("detail", "") + "; ";
      v0.erase(ch);  // no se insiste; se sigue con los demás
    }
  }
  auto plan = bias_.restorePlan(v0);
  for (auto it = plan.begin(); it != plan.end();) {  // los canales apagados no van por escalones
    const auto* pw = lastValue(it->first, "Pw");
    it = (pw && pw->raw != 1u) ? plan.erase(it) : it + 1;
  }
  if (!plan.empty()) {
    PollRequest q;
    q.kind = PollRequest::Kind::V0Set;
    q.ch = plan.front().first;
    q.value = plan.front().second;
    q.internal = true;
    q.reason = restoreReason_;
    q.requester = "hvmon (restauración del barrido)";
    auto fut = q.reply->get_future().share();
    handleV0Set(q, now);
    if (fut.wait_for(std::chrono::seconds(0)) != std::future_status::ready) restoreFut_ = fut;
    if (fut.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
      const auto res = fut.get();
      if (!res.value("ok", false) && res.value("result", "") != "busy") {
        restoreFailed_ = true;
        restoreDetail_ += "CH" + std::to_string(q.ch) + ": " + res.value("detail", "") + "; ";
        // No se repite un escalón rechazado: se termina con error y la sesión queda activa.
        restoreSteps_.clear();
        plan.clear();
      } else {
        return;  // noop o busy: se reintenta en el próximo paso
      }
    } else {
      return;  // escrito: se verifica en las lecturas siguientes
    }
  }
  {
    for (int ch : bias_.limits().channels) {
      const auto* pw = lastValue(ch, "Pw");
      if (pw && pw->raw != 1u && !restorePwTried_.count(ch)) {
        restorePwTried_.insert(ch);  // un intento por canal en esta restauración
        PollRequest q;
        q.kind = PollRequest::Kind::Pw;
        q.ch = ch;
        q.on = true;
        q.internal = true;
        q.reason = restoreReason_;
        q.requester = "hvmon (restauración del barrido)";
        handlePw(q, now);
        return;
      }
    }
  }
  // Fin: todo en el nominal y Pw = 1 (si no falló nada) -> la sesión se cierra.
  bool ok = !restoreFailed_;
  for (int ch : bias_.limits().channels) {
    const auto* a = lastValue(ch, "V0Set");
    const auto* pw = lastValue(ch, "Pw");
    if (!a || !pw || std::fabs(a->si - bias_.nominal(ch)) > bias_.limits().matchTolV || pw->raw != 1u) ok = false;
  }
  restoring_ = false;
  if (ok) {
    bias_.end();
    persistBias();
    logInfo("V0Set nominal restaurado en todos los canales; sesión de barrido cerrada");
  } else {
    logError("restauración del V0Set nominal INCOMPLETA: ", restoreDetail_, " (la sesión queda activa)");
  }
  if (restoreReply_) {
    restoreReply_->set_value(nlohmann::json{{"ok", ok}, {"result", ok ? "restored" : "restore_failed"}, {"detail", restoreDetail_},
                                            {"session", bias_.toJson(clock_.now_utc_ns())}});
    restoreReply_.reset();
  }
}

void Poller::handlePw(PollRequest& r, SteadyTime now) {
  CommandRecord rec;
  rec.id = nextCommandId_++;
  rec.slot = cfg_.slot;
  rec.ch = r.ch;
  rec.value = r.on ? 1 : 0;
  rec.reason = r.reason;
  rec.requester = r.requester;
  rec.tUtcNs = clock_.now_utc_ns();
  rec.phase = "request";
  record(rec);
  auto finish = [&](const std::string& phase, const std::string& detail, bool ok) {
    rec.phase = phase;
    rec.detail = detail;
    rec.tUtcNs = clock_.now_utc_ns();
    record(rec);
    if (!ok) logWarn("comando ", rec.id, " Pw CH", r.ch, "=", rec.value, ": ", phase, ": ", detail);
    r.reply->set_value(nlohmann::json{{"ok", ok}, {"command_id", rec.id}, {"result", phase}, {"detail", detail},
                                      {"ch", r.ch}, {"on", r.on}});
    if (r.retryToken) applyRetry(retry_.onCommandResult(r.ch, r.retryToken, phase, detail, clock_.now_utc_ns()), now);
  };
  if (!cfg_.allowCommands) return finish("rejected", "escrituras deshabilitadas (allow_commands = false)", false);
  if (!contains(cfg_.channels, r.ch)) return finish("rejected", "canal fuera de la configuración", false);
  if (!connected_ || !prevFast_) return finish("rejected", "sin conexión o sin lectura previa", false);
  if (active_)
    return finish(r.internal ? "busy" : "rejected", "otro comando en curso (id " + std::to_string(active_->id) + ")",
                  false);
  if (r.internal && anyRamping()) return finish("busy", "hay un canal en rampa (petición interna)", false);

  // Relectura inmediata antes de escribir: aborta si el canal cambió desde la última lectura publicada.
  auto pre = backend_.readChannels(cfg_.slot, {r.ch}, cfg_.fastParams);
  if (!pre) {
    finish("aborted", "la relectura previa falló: " + pre.error().message, false);
    onReadError(pre.error(), now);
    return;
  }
  for (const char* p : {"Pw", "Status", "V0Set"}) {
    const hv::ParamValue* last = lastValue(r.ch, p);
    const hv::ParamValue* cur = pre.value().get(p, r.ch);
    const bool same = last && cur && (last->isNumeric() ? last->boardValue == cur->boardValue : last->raw == cur->raw);
    if (!same)
      return finish("aborted",
                    std::string("el canal cambió entre la última lectura publicada y la relectura (") + p + ": " +
                        (last ? describeValue(*last) : "?") + " -> " + (cur ? describeValue(*cur) : "?") + ")",
                    false);
  }
  if (pre.value().get("Pw", r.ch)->raw == static_cast<std::uint32_t>(rec.value))
    return finish("noop", "el canal ya está en el estado pedido", true);

  const hv::ParamValue* v0 = lastValue(r.ch, "V0Set");
  const hv::ParamValue* rate = lastValue(r.ch, r.on ? "RUp" : "RDWn");
  const double rateVps = (rate && rate->si > 0) ? rate->si : 1.0;
  const double timeoutS = (v0 ? std::fabs(v0->si) : 0.0) / rateVps + cfg_.verifyMarginS;
  detector_.addPending({rec.id, r.ch, now, now + secondsToNs(timeoutS + cfg_.commandWindowExtraS)});

  if (r.clearAlarmFirst) {
    // ClrAlarm de la placa antes de Pw = 1 (clear_alarm_before_retry); misma compuerta allow_commands.
    CommandRecord ca = rec;
    ca.param = "ClrAlarm";
    ca.value = 1;
    ca.phase = "request";
    ca.tUtcNs = clock_.now_utc_ns();
    record(ca);
    if (auto st = backend_.clearBoardAlarm(cfg_.slot); !st) {
      ca.phase = "write_failed";
      ca.detail = st.error().message;
      record(ca);
      detector_.setPendingEnd(rec.id, now);
      return finish("write_failed", "ClrAlarm: " + st.error().message, false);
    }
    ca.phase = "written";
    ca.tUtcNs = clock_.now_utc_ns();
    record(ca);
  }
  if (auto st = backend_.setChannelPw(cfg_.slot, r.ch, r.on); !st) {
    detector_.setPendingEnd(rec.id, now);
    return finish("write_failed", st.error().message, false);
  }
  rec.phase = "written";
  rec.detail = "timeout de verificación " + std::to_string(timeoutS) + " s";
  rec.tUtcNs = clock_.now_utc_ns();
  record(rec);
  logInfo("comando ", rec.id, ": Pw CH", r.ch, "=", rec.value, " escrito (", r.requester, "): ", r.reason);
  active_ = ActiveCommand{rec.id, r.ch, r.on, now + secondsToNs(timeoutS), r.reason, r.requester, r.reply, r.retryToken};
}

void Poller::verifyCommand(const hv::ChannelReading& rd, SteadyTime now) {
  if (!active_) return;
  const int ch = active_->ch;
  const hv::ParamValue* pw = rd.get("Pw", ch);
  const hv::ParamValue* st = rd.get("Status", ch);
  const hv::ParamValue* vmon = rd.get("VMon", ch);
  const hv::ParamValue* v0 = rd.get("V0Set", ch);
  if (!pw || !st || !vmon || !v0) return;
  bool done;
  if (active_->v0set && active_->offWrite) {
    done = pw->raw == 0u && std::fabs(v0->si - active_->target) <= bias_.limits().matchTolV;
  } else if (active_->v0set) {
    // Cualquier bit fuera de on / ramp_up / ramp_down durante la rampa: abortar (trip, sobrecorriente, ...).
    const std::uint32_t okMask = (1u << 0) | (1u << 1) | (1u << 2);
    if (pw->raw != 1u || (st->raw & ~okMask) != 0u || !cfg_.statusBits.has(st->raw, "on")) {
      std::ostringstream os;
      os << "Status inesperado durante la rampa de V0Set: Pw " << pw->raw << ", Status 0x" << std::hex << st->raw;
      finishCommand("verify_failed", os.str(), false, now);
      return;
    }
    done = !cfg_.statusBits.has(st->raw, "ramp_up") && !cfg_.statusBits.has(st->raw, "ramp_down") &&
           bias_.reached(ch, active_->target, v0->si, vmon->si, active_->v0setBefore, active_->vmonBefore);
  } else if (active_->on)
    done = pw->raw == 1u && cfg_.statusBits.has(st->raw, "on") && !cfg_.statusBits.has(st->raw, "ramp_up") &&
           std::fabs(vmon->si - v0->si) <= cfg_.vmonOnToleranceV;
  else
    // DAQ2-10: apagado verificado por el bit on de Status en 0 (y Pw = 0), no por VMon: la
    // capacidad del PCB del detector descarga con τ ≈ 115 s (CH07), y CH08 (V0Set 0) se daba por
    // apagado con Status = 0x1 por su VMon de 1.5 V de offset.
    done = pw->raw == 0u && !cfg_.statusBits.has(st->raw, "on");
  if (done) {
    if (active_->v0set && !active_->offWrite) {
      bias_.updateOffset(ch, vmon->si - v0->si);
      persistBias();
    }
    std::ostringstream os;
    os << "VMon=" << vmon->si << " V, Status=0x" << std::hex << st->raw;
    finishCommand("verified", os.str(), true, now);
  }
}

void Poller::finishCommand(const std::string& phase, const std::string& detail, bool ok, SteadyTime now) {
  if (!active_) return;
  CommandRecord rec;
  rec.id = active_->id;
  rec.slot = cfg_.slot;
  rec.ch = active_->ch;
  rec.value = active_->on ? 1 : 0;
  if (active_->v0set) {
    rec.param = "V0Set";
    rec.valueV = active_->target;
  }
  rec.phase = phase;
  rec.detail = detail;
  rec.reason = active_->reason;
  rec.requester = active_->requester;
  rec.tUtcNs = clock_.now_utc_ns();
  record(rec);
  if (ok)
    logInfo("comando ", rec.id, " verificado: ", detail);
  else
    logWarn("comando ", rec.id, " ", phase, ": ", detail);
  auto res = nlohmann::json{{"ok", ok}, {"command_id", rec.id}, {"result", phase}, {"detail", detail}, {"ch", rec.ch}, {"on", active_->on}};
  if (active_->v0set) res["value_V"] = active_->target;
  active_->reply->set_value(res);
  detector_.setPendingEnd(rec.id, now + secondsToNs(cfg_.commandWindowExtraS));
  const std::uint64_t token = active_->retryToken;
  const int ch = active_->ch;
  active_.reset();
  if (token) applyRetry(retry_.onCommandResult(ch, token, phase, detail, clock_.now_utc_ns()), now);
}

void Poller::handleOffsets(PollRequest& r, SteadyTime now) {
  (void)now;
  auto fail = [&](const std::string& msg, nlohmann::json extra = nlohmann::json::object()) {
    extra["ok"] = false;
    extra["error"] = msg;
    r.reply->set_value(extra);
  };
  if (offsets_) return fail("ya hay una medida de offsets en curso");
  if (!connected_ || !prevFast_) return fail("sin conexión o sin lectura previa");
  if (r.seconds < 1 || r.seconds > 3600) return fail("seconds debe estar entre 1 y 3600");
  ActiveOffsets o;
  o.seconds = r.seconds;
  o.reply = r.reply;
  const std::vector<int> wanted = r.channels.empty() ? cfg_.channels : r.channels;
  for (int ch : wanted) {
    if (!contains(cfg_.channels, ch)) return fail("canal fuera de la configuración: " + std::to_string(ch));
    const hv::ParamValue* pw = lastValue(ch, "Pw");
    if (!pw || pw->raw != 0u)
      o.refused.push_back(ch);
    else
      o.channels.push_back(ch);
  }
  if (o.channels.empty()) return fail("todos los canales pedidos tienen Pw = 1", {{"refused", o.refused}});
  logInfo("medida de offsets: ", o.channels.size(), " canales, ", r.seconds, " s; rechazados (Pw = 1): ", o.refused.size());
  offsets_ = std::move(o);
}

void Poller::accumulateOffsets(const hv::ChannelReading& rd, SteadyTime now) {
  if (!offsets_) return;
  ActiveOffsets& o = *offsets_;
  if (!o.started) {
    o.started = true;
    o.tStart = now;
    o.tStartUtc = rd.tReqUtcNs;
  }
  for (int ch : o.channels) {
    const hv::ParamValue* pw = rd.get("Pw", ch);
    const hv::ParamValue* i = rd.get("IMon", ch);
    const hv::ParamValue* v = rd.get("VMon", ch);
    if (!pw || !i || !v) return;
    if (pw->raw != 0u) {
      o.reply->set_value({{"ok", false}, {"error", "CH" + std::to_string(ch) + " se encendió durante la medida"}});
      offsets_.reset();
      return;
    }
    o.sumI[ch] += i->si;
    o.sumI2[ch] += i->si * i->si;
    o.sumV[ch] += v->si;
    o.sumV2[ch] += v->si * v->si;
  }
  ++o.n;
  o.tEndUtc = rd.tRespUtcNs;
  if (now - o.tStart < secondsToNs(o.seconds)) return;
  const double n = static_cast<double>(o.n);
  nlohmann::json chans = nlohmann::json::array();
  for (int ch : o.channels) {
    const double mi = o.sumI[ch] / n, mv = o.sumV[ch] / n;
    const double si = n > 1 ? std::sqrt(std::max(0.0, (o.sumI2[ch] - n * mi * mi) / (n - 1))) : 0.0;
    const double sv = n > 1 ? std::sqrt(std::max(0.0, (o.sumV2[ch] - n * mv * mv) / (n - 1))) : 0.0;
    chans.push_back({{"ch", ch}, {"imon_offset_A", mi}, {"imon_std_A", si}, {"vmon_offset_V", mv}, {"vmon_std_V", sv}});
  }
  nlohmann::json j{{"type", "offsets"},        {"t_start_utc_ns", o.tStartUtc}, {"t_end_utc_ns", o.tEndUtc},
                   {"t_start_utc", formatUtcIso(o.tStartUtc)}, {"seconds", o.seconds}, {"n_reads", o.n},
                   {"slot", cfg_.slot},       {"channels", chans},             {"refused_pw_on", o.refused}};
  pushWriter(OffsetItem{j});
  state_.lastOffsets = j;
  logInfo("offsets medidos: ", o.n, " lecturas en ", o.seconds, " s");
  o.reply->set_value({{"ok", true}, {"offsets", j}});
  offsets_.reset();
}

void Poller::enableRetryPersistence(FileSystem& fs) {
  fs_ = &fs;
  const std::string path = cfg_.retry.stateFile;
  if (!fs.exists(path)) return;
  auto content = fs.readFile(path);
  if (!content) {
    logWarn("no se pudo leer el estado de reintento ", path, ": ", content.error().message);
    return;
  }
  try {
    const auto j = nlohmann::json::parse(content.value());
    retry_.restore(j);
    logInfo("estado de reintento cargado de ", path, " (guardado ", j.value("saved_utc", "?"), ")");
  } catch (const std::exception& e) {
    logWarn("estado de reintento ilegible en ", path, ": ", e.what(), " (se ignora)");
  }
}

nlohmann::json Poller::channelContext(int ch) const {
  nlohmann::json c = nlohmann::json::object();
  auto it = state_.values.find(ch);
  if (it != state_.values.end()) {
    for (const auto& info : state_.params) {
      auto v = it->second.find(info.name);
      if (v == it->second.end()) continue;
      if (info.name == "VMon") c["VMon_V"] = v->second.si;
      if (info.name == "IMon") c["IMon_A"] = v->second.si;
      if (info.name == "V0Set") c["V0Set_V"] = v->second.si;
      if (info.name == "Status") {
        c["Status"] = v->second.raw;
        c["status_flags"] = cfg_.statusBits.decode(v->second.raw);
      }
    }
  }
  c["board_temperature_degC"] = nlohmann::json();
  if (state_.board)
    if (const auto* t = state_.board->get("Temp")) c["board_temperature_degC"] = t->si;
  c["geco_running"] = hub_.gecoRunning.load();
  return c;
}

void Poller::runRetry(const hv::ChannelReading& rd, const std::vector<ChangeRecord>& changes, SteadyTime now) {
  RetryInput in;
  in.nowUtc = clock_.now_utc_ns();
  in.fresh = connected_;
  in.geco = hub_.gecoRunning.load();
  in.changes = changes;
  for (int ch : rd.channels) {
    RetryChannelView v;
    if (const auto* x = rd.get("Pw", ch)) v.pw = x->raw;
    if (const auto* x = rd.get("Status", ch)) v.status = x->raw;
    if (const auto* x = rd.get("VMon", ch)) v.vmonV = x->si;
    if (const auto* x = rd.get("IMon", ch)) v.imonA = x->si;
    if (const auto* x = rd.get("V0Set", ch)) v.v0setV = x->si;
    in.channels[ch] = v;
  }
  if (state_.lastOffsets.is_object() && state_.lastOffsets.contains("channels"))
    for (const auto& c : state_.lastOffsets["channels"]) in.imonOffsetA[c.value("ch", -1)] = c.value("imon_offset_A", 0.0);
  applyRetry(retry_.update(in), now);
}

void Poller::applyRetry(RetryOutput out, SteadyTime now) {
  std::vector<ChangeRecord> recs;
  for (const auto& t : out.transitions) {
    ChangeRecord r;
    r.id = nextChangeId_++;
    r.tUtcNs = t.tUtcNs;
    r.tPrevReadUtcNs = t.tUtcNs;
    r.slot = cfg_.slot;
    r.ch = t.ch;
    r.param = "auto_retry";
    r.oldValue = t.from;
    r.newValue = t.to;
    r.origin = "hvmon_retry";
    r.readKind = "retry";
    r.detail = {{"k", t.k}, {"reason", t.reason}};
    recs.push_back(r);
  }
  addChanges(recs);
  for (const auto& w : out.warnings) logWarn(w);
  for (const auto& e : out.events) {
    NotifyEvent n;
    n.kind = e.kind;
    n.ch = e.ch;
    n.slot = cfg_.slot;
    n.episode = e.episode;
    n.k = e.k;
    n.reason = e.reason;
    n.timeline = e.timeline;
    n.tUtcNs = e.tUtcNs;
    n.context = channelContext(e.ch);
    n.simulated = cfg_.retry.inTestMode(e.ch);
    if (notifyQ_ && !notifyQ_->try_push(std::move(n), DropPolicy::DropNewest)) logError("cola del notificador llena");
  }
  if (out.persist && fs_) {
    const std::string path = cfg_.retry.stateFile;
    fs_->makeDirectories(parentDirectory(path), 0700);
    if (auto st = fs_->writeAtomic(path, retry_.persistJson(clock_.now_utc_ns()).dump(1) + "\n", 0600); !st)
      logError("no se pudo guardar el estado de reintento: ", st.error().message);
  }
  for (const auto& a : out.actions) {
    PollRequest r;
    r.kind = PollRequest::Kind::Pw;
    r.ch = a.ch;
    r.on = true;
    r.internal = true;
    r.retryToken = a.token;
    r.clearAlarmFirst = a.clearAlarmFirst;
    r.reason = "auto-retry intento " + std::to_string(a.k);
    r.requester = "hvmon (auto-retry)";
    handlePw(r, now);
  }
}

void Poller::publish(SteadyTime now) {
  state_.retry = retry_.publish(clock_.now_utc_ns());
  state_.readSeq = readSeq_;
  state_.connected = connected_;
  state_.everConnected = everConnected_;
  state_.stats = stats_;
  for (const auto& [ch, w] : window_) state_.window[ch] = {w.first.summary(state_.lastGoodReadUtcNs), w.second.summary(state_.lastGoodReadUtcNs)};
  state_.nextConnectInS =
      connected_ ? 0.0 : std::max(0.0, std::chrono::duration<double>(nextConnect_ - now).count());
  if (active_)
    state_.activeCommand = {{"command_id", active_->id}, {"ch", active_->ch}, {"on", active_->on},
                            {"deadline_in_s", std::chrono::duration<double>(active_->deadline - now).count()},
                            {"requester", active_->requester}, {"reason", active_->reason}};
  else
    state_.activeCommand = nlohmann::json();
  if (offsets_)
    state_.activeOffsets = {{"channels", offsets_->channels}, {"seconds", offsets_->seconds}, {"n_reads", offsets_->n}};
  else
    state_.activeOffsets = nlohmann::json();
  hub_.poller.publish(state_);
}

void Poller::shutdown() {
  const nlohmann::json bye{{"ok", false}, {"error", "hvmon se está apagando"}};
  if (active_) {
    active_->reply->set_value(bye);
    active_.reset();
  }
  if (offsets_) {
    offsets_->reply->set_value(bye);
    offsets_.reset();
  }
  while (auto r = requests_.try_pop()) r->reply->set_value(bye);
  backend_.disconnect();
  connected_ = false;
  publish(clock_.now_steady());
}

}  // namespace daq::hvmon
