// hv_client.cpp
#include "daq/daqrun/hv_client.h"

#include <signal.h>

#include <algorithm>
#include <cerrno>
#include <cmath>

#include "daq/core/file_system.h"
#include "daq/core/logger.h"
#include "daq/trigger_interval_stats.h"

namespace daq::daqrun {

namespace {

constexpr double kNs = 1e9;

std::uint32_t uintValue(const nlohmann::json& v, const char* k) {
  if (!v.contains(k) || !v[k].is_number()) return 0;
  return static_cast<std::uint32_t>(v[k].get<double>());
}
double numValue(const nlohmann::json& v, const char* k) {
  if (!v.contains(k) || !v[k].is_number()) return 0.0;
  return v[k].get<double>();
}

bool pidExists(long pid) {
  if (pid <= 0) return false;
  return ::kill(static_cast<pid_t>(pid), 0) == 0 || errno == EPERM;
}

}  // namespace

HvSample parseHvState(const nlohmann::json& s, UtcNs now, double staleAfterS, bool pidAlive) {
  HvSample h;
  h.receivedUtcNs = now;
  try {
    h.fileOk = s.value("schema", "") == "daqv2.hvmon_state/1";
    if (!h.fileOk) {
      h.error = "esquema desconocido";
      return h;
    }
    h.stateSeq = s.value("seq", std::uint64_t{0});
    h.stateUtcNs = s.value("t_utc_ns", UtcNs{0});
    h.pid = s.value("pid", 0L);
    const bool final = s.value("final", false);
    const auto conn = s.value("connection", nlohmann::json::object());
    h.connected = conn.value("connected", false);
    const bool hvStale = conn.value("stale", true);
    if (s.contains("last_read")) {
      const auto& lr = s["last_read"];
      h.readSeq = lr.value("read_seq", std::uint64_t{0});
      h.tReqUtcNs = lr.value("t_req_utc_ns", UtcNs{0});
      h.tRespUtcNs = lr.value("t_resp_utc_ns", UtcNs{0});
    } else {  // hvmon anterior a DAQ2-5: sin t_req; se usa la última lectura buena para ambos
      h.readSeq = s.value("read_seq", std::uint64_t{0});
      h.tReqUtcNs = h.tRespUtcNs = conn.value("last_good_read_utc_ns", UtcNs{0});
    }
    const double age = static_cast<double>(now - h.stateUtcNs) / kNs;
    h.alive = !final && pidAlive && age <= staleAfterS;
    h.stale = !h.alive || hvStale;
    if (final) h.error = "hvmon terminó (estado final)";
    else if (!pidAlive) h.error = "el proceso de hvmon no existe";
    else if (age > staleAfterS) h.error = "estado de hvmon viejo (" + std::to_string(age) + " s)";
    for (const auto& c : s.value("channels", nlohmann::json::array())) {
      const int ch = c.value("ch", -1);
      if (ch < 0 || ch >= kHvChannels) continue;
      HvChannelSample& x = h.ch[static_cast<std::size_t>(ch)];
      const auto v = c.value("values", nlohmann::json::object());
      x.present = !v.empty();
      x.pw = uintValue(v, "Pw");
      x.status = uintValue(v, "Status");
      x.v0setV = numValue(v, "V0Set_V");
      x.i0setA = numValue(v, "I0Set_A");
      x.vmonV = numValue(v, "VMon_V");
      x.imonA = numValue(v, "IMon_A");
      if (c.contains("auto_retry") && c["auto_retry"].is_object()) x.retryState = c["auto_retry"].value("label", "");
      else if (c.contains("auto_retry") && c["auto_retry"].is_string()) x.retryState = c["auto_retry"].get<std::string>();
    }
    h.newChanges = s.value("changes", nlohmann::json::array());
  } catch (const std::exception& e) {
    h.fileOk = false;
    h.alive = false;
    h.stale = true;
    h.error = std::string("estado de hvmon inválido: ") + e.what();
  }
  return h;
}

HvClient::HvClient(const DaqrunConfig& cfg, Clock& clock, BoundedQueue<HvSample>& toWriter)
    : cfg_(cfg), clock_(clock), out_(toWriter) {}

HvSample HvClient::pollOnce() {
  PosixFileSystem fs;
  const UtcNs now = clock_.now_utc_ns();
  HvSample s;
  auto text = fs.readFile(cfg_.hvStateFile);
  if (!text) {
    s.receivedUtcNs = now;
    s.error = "hvmon no corre o no publica estado: " + text.error().message;
    return s;
  }
  nlohmann::json j;
  try {
    j = nlohmann::json::parse(text.value());
  } catch (const std::exception& e) {
    s.receivedUtcNs = now;
    s.error = std::string("estado de hvmon ilegible: ") + e.what();
    return s;
  }
  s = parseHvState(j, now, cfg_.hvStaleAfterS, pidExists(j.value("pid", 0L)));
  s.raw = std::move(j);
  // Deduplicación de cambios por id.
  nlohmann::json fresh = nlohmann::json::array();
  for (const auto& c : s.newChanges) {
    const std::int64_t id = c.value("id", std::int64_t{-1});
    if (id > lastChangeId_) {
      fresh.push_back(c);
      lastChangeId_ = std::max(lastChangeId_, id);
    }
  }
  s.newChanges = fresh;
  return s;
}

void HvClient::run(const StopToken& st) {
  Logger::setThreadName("hvclient");
  const auto period = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(cfg_.hvPollPeriodS));
  SteadyTime next = clock_.now_steady();
  bool lastAlive = true;
  while (!st.stop_requested()) {
    HvSample s = pollOnce();
    if (s.alive != lastAlive) {
      if (s.alive) logInfo("hvmon vivo de nuevo (read_seq ", s.readSeq, ")");
      else logWarn("hvmon no disponible: ", s.error, " (bit 4 en todos los canales)");
      lastAlive = s.alive;
    }
    HvSample forPage = s;
    forPage.raw = nlohmann::json();  // la página no necesita el estado completo
    latest_.publish(forPage);
    s.raw = nlohmann::json();
    if (!out_.push(std::move(s))) break;  // cola cerrada: el escritor terminó
    ++pushed_;
    next += period;
    if (next < clock_.now_steady()) next = clock_.now_steady() + period;
    clock_.sleep_until(next, st);
  }
}

// ---------------------------------------------------------------------------------------------

FlagEngine::FlagEngine(const DaqrunConfig& cfg, const hvmon::Mapping& mapping) : cfg_(cfg) {
  hvOf_.fill(-1);
  static_.fill(chflag::MapUnconfirmed);  // sin entrada en el mapeo: no confirmado
  for (const auto& e : mapping.channels) {
    if (e.digCh < 0 || e.digCh >= 64) continue;
    const auto i = static_cast<std::size_t>(e.digCh);
    if (e.hvGroup) {
      auto it = mapping.groupToHvCh.find(*e.hvGroup);
      hvOf_[i] = e.hvCh ? *e.hvCh : (it != mapping.groupToHvCh.end() ? it->second : *e.hvGroup - 1);
    }
    static_[i] = (e.groupConfirmed() && e.faceIsConfirmed()) ? 0 : chflag::MapUnconfirmed;
  }
}

void FlagEngine::addSample(const HvSample& s) {
  live_.push_back(Liveness{s.receivedUtcNs, s.fileOk && s.alive && !s.stale});
  if (!s.fileOk || !s.connected || s.tRespUtcNs <= 0) return;
  if (!reads_.empty() && reads_.back().readSeq == s.readSeq) return;
  if (!reads_.empty() && s.tRespUtcNs <= reads_.back().tResp) return;
  reads_.push_back(Read{s.readSeq, s.tReqUtcNs, s.tRespUtcNs, s.ch});
  const auto onBit = cfg_.statusBits.bitOf("on");
  for (int h = 0; h < kHvChannels; ++h) {
    const auto& c = s.ch[static_cast<std::size_t>(h)];
    if (!c.present || c.pw != 1 || (onBit >= 0 && !(c.status & (1u << onBit)))) continue;
    if (stateBits(c, h) & (chflag::HvRamp | chflag::HvOff)) continue;
    auto& d = imonHist_[static_cast<std::size_t>(h)];
    d.emplace_back(s.tRespUtcNs, c.imonA);
    const UtcNs cut = s.tRespUtcNs - static_cast<UtcNs>(cfg_.imonBaselineWindowS * kNs);
    while (!d.empty() && d.front().first < cut) d.pop_front();
  }
}

double FlagEngine::imonBaseline(int h) const {
  const auto& d = imonHist_[static_cast<std::size_t>(h)];
  if (d.size() < 60) return 0.0;  // menos de un minuto de NORMAL: sin línea base
  std::vector<double> v;
  v.reserve(d.size());
  for (const auto& p : d) v.push_back(p.second);
  auto mid = v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2);
  std::nth_element(v.begin(), mid, v.end());
  return *mid;
}

std::uint16_t FlagEngine::stateBits(const HvChannelSample& c, int h) const {
  if (!c.present) return chflag::HvUnknown;
  const bool expected = std::find(cfg_.expectedOn.begin(), cfg_.expectedOn.end(), h) != cfg_.expectedOn.end();
  const auto& sb = cfg_.statusBits;
  const bool on = sb.has(c.status, "on");
  bool trip = false;
  for (const auto& f : cfg_.tripFlags) trip = trip || sb.has(c.status, f);
  const bool ramp = sb.has(c.status, "ramp_up") || sb.has(c.status, "ramp_down");
  std::uint16_t f = 0;
  if (expected && (c.pw == 0 || !on || trip)) f |= chflag::HvOff;
  if (sb.has(c.status, "overcurrent")) f |= chflag::HvOvc;
  if (ramp) f |= chflag::HvRamp;
  if (c.pw == 1 && !ramp && std::fabs(c.vmonV - c.v0setV) > cfg_.vmonToleranceV) f |= chflag::HvVdev;
  return f;
}

bool FlagEngine::differs(const HvChannelSample& a, const HvChannelSample& b) const {
  return a.pw != b.pw || a.status != b.status || a.v0setV != b.v0setV;
}

bool FlagEngine::isTrip(const HvChannelSample& a, const HvChannelSample& b) const {
  bool trip = false;
  for (const auto& f : cfg_.tripFlags) trip = trip || cfg_.statusBits.has(b.status, f);
  return a.pw == 1 && b.pw == 0 && trip;
}

ChannelFlags FlagEngine::flagsFor(UtcNs t) const {
  ChannelFlags out;
  // Estado de hvmon en t (por la hora de lectura de daqrun).
  bool liveOk = false;
  for (auto it = live_.rbegin(); it != live_.rend(); ++it)
    if (it->t <= t) {
      liveOk = it->ok && static_cast<double>(t - it->t) / kNs <= cfg_.hvStaleAfterS;
      break;
    }
  const Read* prev = nullptr;
  const Read* next = nullptr;
  for (const auto& r : reads_) {
    if (r.tResp <= t) prev = &r;
    else {
      next = &r;
      break;
    }
  }
  const bool unknownAll = !liveOk || !prev || !next || static_cast<double>(t - prev->tResp) / kNs > cfg_.hvStaleAfterS;
  if (prev) out.hvSeq = prev->readSeq;
  // Intervalos de ambigüedad por canal HV.
  std::array<std::uint16_t, kHvChannels> ambBits{};
  std::array<bool, kHvChannels> amb{};
  const UtcNs back = static_cast<UtcNs>(cfg_.tripBackdateS * kNs);
  for (std::size_t i = 1; i < reads_.size(); ++i) {
    const Read& a = reads_[i - 1];
    const Read& b = reads_[i];
    if (b.tResp < t) continue;
    if (a.tReq - back > t) break;
    for (int h = 0; h < kHvChannels; ++h) {
      const auto& ca = a.ch[static_cast<std::size_t>(h)];
      const auto& cb = b.ch[static_cast<std::size_t>(h)];
      if (!differs(ca, cb)) continue;
      const UtcNs start = a.tReq - (isTrip(ca, cb) ? back : 0);
      if (t > start && t <= b.tResp) {
        amb[static_cast<std::size_t>(h)] = true;
        ambBits[static_cast<std::size_t>(h)] |= static_cast<std::uint16_t>(stateBits(ca, h) | stateBits(cb, h));
      }
    }
  }
  std::array<std::uint16_t, kHvChannels> base{};
  if (!unknownAll)
    for (int h = 0; h < kHvChannels; ++h) {
      const auto& c = prev->ch[static_cast<std::size_t>(h)];
      std::uint16_t f = stateBits(c, h);
      if (c.pw == 1 && !(f & chflag::HvOvc)) {
        const double bl = imonBaseline(h);
        if (bl > 0 && c.imonA > cfg_.ovcKOverBaseline * bl) f |= chflag::HvOvc;
      }
      base[static_cast<std::size_t>(h)] = f;
    }
  for (int d = 0; d < 64; ++d) {
    std::uint16_t f = static_[static_cast<std::size_t>(d)];
    const int h = hvOf_[static_cast<std::size_t>(d)];
    if (h < 0 || h >= kHvChannels || unknownAll) {
      f |= chflag::HvUnknown;
    } else {
      f |= base[static_cast<std::size_t>(h)];
      if (amb[static_cast<std::size_t>(h)]) f |= static_cast<std::uint16_t>(chflag::HvAmbiguous | ambBits[static_cast<std::size_t>(h)]);
    }
    out.ch[static_cast<std::size_t>(d)] = f;
    if (f & 0x3F) out.anyHv = true;
  }
  return out;
}

void FlagEngine::prune(UtcNs keepAfter) {
  while (reads_.size() > 2 && reads_[1].tResp < keepAfter) reads_.pop_front();
  while (live_.size() > 2 && live_[1].t < keepAfter) live_.pop_front();
}

// ---------------------------------------------------------------------------------------------

EventFlagger::EventFlagger(double gapProbability, double burstProbability, double tttPcToleranceS, int counterBits)
    : gapP_(gapProbability),
      burstP_(burstProbability),
      tolS_(tttPcToleranceS),
      counterMask_(counterBits >= 32 ? 0xFFFFFFFFu : static_cast<std::uint32_t>((1ull << counterBits) - 1)) {}

double EventFlagger::rateHz() const {
  if (recent_.size() < 20) return 0.0;
  const double span = static_cast<double>(recent_.back() - recent_.front()) / kNs;
  return span > 0 ? static_cast<double>(recent_.size() - 1) / span : 0.0;
}

std::uint32_t EventFlagger::flags(UtcNs tUtc, UtcNs tRead, std::uint32_t counterHw, bool boardFull, std::int64_t dtPrevNs) {
  std::uint32_t f = 0;
  if (haveCounter_ && ((counterHw - lastCounter_) & counterMask_) > 1) f |= evtflag::HwCounterJump;
  haveCounter_ = true;
  lastCounter_ = counterHw;
  if (boardFull) f |= evtflag::BoardFull;
  if (std::fabs(static_cast<double>(tRead - tUtc)) / kNs > tolS_) f |= evtflag::TttPcMismatch;
  const double r = rateHz();  // tasa antes de este evento
  if (r > 0 && dtPrevNs > 0 && static_cast<double>(dtPrevNs) / kNs > -std::log(gapP_) / r) f |= evtflag::Gap;
  recent_.push_back(tUtc);
  while (!recent_.empty() && recent_.front() < tUtc - static_cast<UtcNs>(3600 * kNs)) recent_.pop_front();
  if (r > 0) {
    const double muW = 5.0;
    const UtcNs w = static_cast<UtcNs>(muW / r * kNs);
    std::size_t n = 0;
    for (auto it = recent_.rbegin(); it != recent_.rend() && *it > tUtc - w; ++it) ++n;
    if (tis::poissonUpperTail(n, muW) < burstP_) f |= evtflag::Burst;
  }
  return f;
}

TttUnwrapper::TttUnwrapper(UtcNs originUtc, double tickNs, int wrapBits)
    : origin_(originUtc),
      tick_(tickNs),
      period_(static_cast<std::int64_t>(std::llround(static_cast<double>(wrapBits >= 32 ? (1ull << 32) : (1ull << wrapBits)) * tickNs))) {}

std::int64_t TttUnwrapper::tttNs(std::uint32_t raw, UtcNs tReadUtc) {
  const auto base = static_cast<std::int64_t>(std::llround(static_cast<double>(raw) * tick_));
  const std::int64_t elapsed = tReadUtc - origin_;
  std::int64_t k = static_cast<std::int64_t>(std::floor(static_cast<double>(elapsed - base) / static_cast<double>(period_) + 0.5));
  if (k < 0) k = 0;
  std::int64_t t = base + k * period_;
  if (last_ >= 0 && t < last_) t += ((last_ - t) / period_ + 1) * period_;
  last_ = t;
  return t;
}

}  // namespace daq::daqrun
