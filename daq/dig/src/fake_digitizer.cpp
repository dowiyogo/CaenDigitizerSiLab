// fake_digitizer.cpp
#include "daq/dig/fake_digitizer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace daq::dig {

namespace {

constexpr double kSampleNs = 16.0;
constexpr double kNsPerS = 1e9;

int stationIndex(const std::string& s) { return s == "hod2019" ? 1 : 0; }

UtcNs secondsToNs(double s) { return static_cast<UtcNs>(std::llround(s * kNsPerS)); }

void readPair(const nlohmann::json& j, const char* key, double out[2]) {
  if (!j.contains(key)) return;
  const auto& v = j[key];
  if (v.is_number()) {
    out[0] = out[1] = v.get<double>();
  } else {
    out[0] = v.value("hod2018", out[0]);
    out[1] = v.value("hod2019", out[1]);
  }
}

nlohmann::json pairJson(const double v[2]) { return {{"hod2018", v[0]}, {"hod2019", v[1]}}; }

}  // namespace

std::vector<ChannelInfo> defaultChannelInfo() {
  std::vector<ChannelInfo> v(kChannels);
  for (int ch = 0; ch < kChannels; ++ch) {
    v[static_cast<std::size_t>(ch)].station = ch < 32 ? 0 : 1;
    v[static_cast<std::size_t>(ch)].plane = (ch % 32) < 16 ? 'X' : 'Y';
    v[static_cast<std::size_t>(ch)].strip = ch % 16;
  }
  return v;
}

Result<SyntheticParams> SyntheticParams::fromJson(const nlohmann::json& j) {
  using R = Result<SyntheticParams>;
  SyntheticParams p;
  try {
    p.seed = j.value("seed", p.seed);
    p.triggerRateHz = j.value("trigger_rate_Hz", p.triggerRateHz);
    p.signalSample = j.value("signal_sample", p.signalSample);
    readPair(j, "baseline_adc", p.baselineAdc);
    p.baselineSpreadAdc = j.value("baseline_spread_adc", p.baselineSpreadAdc);
    readPair(j, "noise_mV", p.noiseMv);
    readPair(j, "dark_rate_Hz", p.darkRateHz);
    if (j.contains("gain_mV_per_pe")) {
      const auto& g = j["gain_mV_per_pe"];
      for (int s = 0; s < 2; ++s) {
        const auto& v = g.at(s == 0 ? "hod2018" : "hod2019");
        if (v.is_array()) {
          p.gainMinMv[s] = v.at(0).get<double>();
          p.gainMaxMv[s] = v.at(1).get<double>();
        } else {
          p.gainMinMv[s] = p.gainMaxMv[s] = v.get<double>();
        }
      }
    }
    readPair(j, "crosstalk_prob", p.crosstalkProb);
    p.decaySamples = j.value("decay_samples", p.decaySamples);
    readPair(j, "muon_pe_mean", p.muonPeMean);
    p.muonEfficiency = j.value("muon_efficiency", p.muonEfficiency);
    p.slopeSigmaStrips = j.value("slope_sigma_strips", p.slopeSigmaStrips);
    p.tttWrapBits = j.value("ttt_wrap_bits", p.tttWrapBits);
    p.counterBits = j.value("counter_bits", p.counterBits);
    p.boardMemoryEvents = j.value("board_memory_events", p.boardMemoryEvents);
    if (j.contains("events")) p.events = j["events"];
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("escenario del FakeDigitizer inválido: ") + e.what());
  }
  if (p.triggerRateHz < 0 || p.decaySamples <= 0 || p.tttWrapBits < 8 || p.tttWrapBits > 32 || p.counterBits < 8 ||
      p.counterBits > 32 || p.boardMemoryEvents == 0)
    return R::failure(-1, "escenario del FakeDigitizer: parámetros fuera de rango");
  for (const auto& e : p.events) {
    const std::string t = e.value("type", "");
    if (t != "hv_group_off" && t != "light_leak" && t != "gain_drift" && t != "gain_step" && t != "trigger_gap" &&
        t != "trigger_burst")
      return R::failure(-1, "escenario del FakeDigitizer: evento desconocido '" + t + "'");
    if (!e.contains("t_s")) return R::failure(-1, "escenario del FakeDigitizer: evento sin t_s");
  }
  return p;
}

nlohmann::json SyntheticParams::toJson() const {
  return {{"seed", seed},
          {"trigger_rate_Hz", triggerRateHz},
          {"signal_sample", signalSample},
          {"baseline_adc", pairJson(baselineAdc)},
          {"baseline_spread_adc", baselineSpreadAdc},
          {"noise_mV", pairJson(noiseMv)},
          {"dark_rate_Hz", pairJson(darkRateHz)},
          {"gain_mV_per_pe",
           {{"hod2018", {gainMinMv[0], gainMaxMv[0]}}, {"hod2019", {gainMinMv[1], gainMaxMv[1]}}}},
          {"crosstalk_prob", pairJson(crosstalkProb)},
          {"decay_samples", decaySamples},
          {"muon_pe_mean", pairJson(muonPeMean)},
          {"muon_efficiency", muonEfficiency},
          {"slope_sigma_strips", slopeSigmaStrips},
          {"ttt_wrap_bits", tttWrapBits},
          {"counter_bits", counterBits},
          {"board_memory_events", boardMemoryEvents},
          {"events", events}};
}

FakeDigitizer::FakeDigitizer(SyntheticParams p, Clock& clock, std::vector<ChannelInfo> channels)
    : p_(std::move(p)), clock_(clock), ch_(std::move(channels)), rng_(p_.seed) {
  if (ch_.size() != static_cast<std::size_t>(kChannels)) ch_ = defaultChannelInfo();
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  baseline_.resize(kChannels);
  gainMv_.resize(kChannels);
  for (int c = 0; c < kChannels; ++c) {
    const int s = ch_[static_cast<std::size_t>(c)].station;
    baseline_[static_cast<std::size_t>(c)] = p_.baselineAdc[s] + p_.baselineSpreadAdc * u(rng_);
    gainMv_[static_cast<std::size_t>(c)] = p_.gainMinMv[s] + (p_.gainMaxMv[s] - p_.gainMinMv[s]) * u01(rng_);
  }
}

Status FakeDigitizer::configure(const DigitizerConfig& cfg) {
  cfg_ = cfg;
  configured_ = true;
  return okStatus();
}

Result<AcqOrigin> FakeDigitizer::start() {
  if (!configured_) return Result<AcqOrigin>::failure(-1, "FakeDigitizer: start() sin configure()");
  AcqOrigin o;
  o.utcBeforeStart = clock_.now_utc_ns();
  t0_ = o.utcBeforeStart;
  o.utcAfterStart = clock_.now_utc_ns();
  running_ = true;
  generatedUntil_ = t0_;
  injected_.clear();
  memory_.clear();
  // Ráfagas del guion: n triggers uniformes en [t, t + duración).
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  for (const auto& e : p_.events) {
    if (e.value("type", "") != "trigger_burst") continue;
    const UtcNs t = t0_ + secondsToNs(e["t_s"].get<double>());
    const double d = e.value("duration_s", 1.0);
    const int n = e.value("n", 20);
    for (int i = 0; i < n; ++i) injected_.push_back(t + secondsToNs(d * u01(rng_)));
  }
  std::sort(injected_.begin(), injected_.end());
  nextPoisson_ = cfg_.softwareTrigger ? INT64_MAX : drawNextPoisson(t0_);
  return o;
}

Status FakeDigitizer::stop() {
  running_ = false;
  return okStatus();
}

bool FakeDigitizer::inGap(UtcNs t, UtcNs& gapEnd) const {
  for (const auto& e : p_.events) {
    if (e.value("type", "") != "trigger_gap") continue;
    const UtcNs a = t0_ + secondsToNs(e["t_s"].get<double>());
    const UtcNs b = a + secondsToNs(e.value("duration_s", 60.0));
    if (t >= a && t < b) {
      gapEnd = b;
      return true;
    }
  }
  return false;
}

UtcNs FakeDigitizer::drawNextPoisson(UtcNs from) {
  if (p_.triggerRateHz <= 0) return INT64_MAX;
  std::exponential_distribution<double> ex(p_.triggerRateHz);
  UtcNs t = from + secondsToNs(ex(rng_));
  UtcNs gapEnd = 0;
  while (inGap(t, gapEnd)) t = gapEnd + secondsToNs(ex(rng_));  // sin memoria: se redibuja desde el fin
  return t;
}

void FakeDigitizer::generateUntil(UtcNs now) {
  const std::uint32_t counterMask =
      p_.counterBits >= 32 ? 0xFFFFFFFFu : static_cast<std::uint32_t>((1ull << p_.counterBits) - 1);
  for (;;) {
    const UtcNs nextInj = injected_.empty() ? INT64_MAX : injected_.front();
    const UtcNs t = std::min(nextPoisson_, nextInj);
    if (t > now || t == INT64_MAX) break;
    std::uint8_t kind = 0;
    if (t == nextInj) {
      kind = static_cast<std::uint8_t>(cfg_.softwareTrigger ? 1 : 0);
      injected_.erase(injected_.begin());
    } else {
      nextPoisson_ = drawNextPoisson(t);
    }
    ++generated_;
    counter_ = (counter_ + 1) & counterMask;
    if (memory_.size() >= p_.boardMemoryEvents) {  // memoria llena: el trigger se pierde
      ++lostFull_;
      fullSinceLastRead_ = true;
      continue;
    }
    // El contador cuenta todos los triggers: los perdidos por memoria llena dejan un salto.
    memory_.push_back(Pending{t, kind, counter_});
  }
  generatedUntil_ = now;
}

Status FakeDigitizer::sendSoftwareTrigger() {
  if (!running_) return Status::failure(-1, "FakeDigitizer: adquisición detenida");
  const UtcNs now = clock_.now_utc_ns();
  injected_.insert(std::upper_bound(injected_.begin(), injected_.end(), now), now);
  return okStatus();
}

double FakeDigitizer::gainMv(int ch, UtcNs t) const {
  double g = gainMv_[static_cast<std::size_t>(ch)];
  const ChannelInfo& ci = ch_[static_cast<std::size_t>(ch)];
  for (const auto& e : p_.events) {
    const std::string type = e.value("type", "");
    if (type != "gain_drift" && type != "gain_step") continue;
    const UtcNs ts = t0_ + secondsToNs(e["t_s"].get<double>());
    if (t < ts) continue;
    bool applies = true;
    if (e.contains("station")) applies = stationIndex(e["station"].get<std::string>()) == ci.station;
    if (e.contains("channels")) {
      applies = false;
      for (const auto& c : e["channels"]) applies = applies || c.get<int>() == ch;
    }
    if (!applies) continue;
    if (type == "gain_drift")
      g *= 1.0 + e.value("per_hour", 0.0) * static_cast<double>(t - ts) / (3600.0 * kNsPerS);
    else
      g *= e.value("factor", 1.0);
  }
  return g;
}

bool FakeDigitizer::channelOff(int ch, UtcNs t) const {
  const int g = ch_[static_cast<std::size_t>(ch)].hvGroup;
  if (g <= 0) return false;
  for (const auto& e : p_.events) {
    if (e.value("type", "") != "hv_group_off" || e.value("group", 0) != g) continue;
    const UtcNs a = t0_ + secondsToNs(e["t_s"].get<double>());
    const UtcNs b = a + secondsToNs(e.value("duration_s", 300.0));
    if (t >= a && t < b) return true;
  }
  return false;
}

double FakeDigitizer::lightLeakExtraHz(int ch, UtcNs t) const {
  double extra = 0;
  const std::string& face = ch_[static_cast<std::size_t>(ch)].face;
  if (face.empty()) return 0;
  for (const auto& e : p_.events) {
    if (e.value("type", "") != "light_leak" || e.value("face", "") != face) continue;
    const UtcNs a = t0_ + secondsToNs(e["t_s"].get<double>());
    const UtcNs b = a + secondsToNs(e.value("duration_s", 120.0));
    if (t >= a && t < b) extra += e.value("extra_rate_Hz", 1e5);
  }
  return extra;
}

Event FakeDigitizer::synthesize(const Pending& pd) {
  const int n = cfg_.recordLength;
  const double lsb = AdcConversion{cfg_.nBit, cfg_.maxVpp, 0x7FFF}.lsbMv();
  const double maxAdc = static_cast<double>((1 << cfg_.nBit) - 1);
  Event ev;
  ev.counterHw = pd.counter;
  const double tttNs = static_cast<double>(pd.t - t0_);
  const std::uint64_t ticks = static_cast<std::uint64_t>(tttNs / cfg_.tttTickNs);
  const std::uint64_t wrap = p_.tttWrapBits >= 32 ? (1ull << 32) : (1ull << p_.tttWrapBits);
  ev.tttRaw = static_cast<std::uint32_t>(ticks % wrap);
  ev.triggerKind = pd.kind;
  ev.adc.assign(static_cast<std::size_t>(kChannels * n), 0);

  // Traza del muón (solo triggers externos): posición en unidades de tira y pendiente por estación.
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::normal_distribution<double> slope(0.0, p_.slopeSigmaStrips);
  double pos[2][2];  // [estación][plano X=0, Y=1]
  const bool muon = pd.kind == 0;
  if (muon) {
    pos[0][0] = 16.0 * u01(rng_);
    pos[0][1] = 16.0 * u01(rng_);
    pos[1][0] = pos[0][0] + slope(rng_);
    pos[1][1] = pos[0][1] + slope(rng_);
  }
  const double tail = 8.0 * p_.decaySamples;
  std::vector<double> wf(static_cast<std::size_t>(n));
  auto addPulse = [&](double s, double ampMv) {
    const double a = ampMv / lsb;
    const int k0 = std::max(0, static_cast<int>(std::ceil(s)));
    const int k1 = std::min(n, static_cast<int>(s + tail) + 1);
    for (int k = k0; k < k1; ++k) wf[static_cast<std::size_t>(k)] -= a * std::exp(-(k - s) / p_.decaySamples);
  };
  auto drawPe = [&](int st) {
    int pe = 1;
    while (u01(rng_) < p_.crosstalkProb[st]) ++pe;
    return pe;
  };
  for (int c = 0; c < kChannels; ++c) {
    const ChannelInfo& ci = ch_[static_cast<std::size_t>(c)];
    const int st = ci.station;
    std::fill(wf.begin(), wf.end(), baseline_[static_cast<std::size_t>(c)]);
    const bool off = channelOff(c, pd.t);
    if (!off) {
      const double g = gainMv(c, pd.t);
      const double rate = p_.darkRateHz[st] + lightLeakExtraHz(c, pd.t);
      const double span = (n + tail) * kSampleNs * 1e-9;
      std::poisson_distribution<int> nd(rate * span);
      const int nDark = nd(rng_);
      std::normal_distribution<double> spread(1.0, 0.1);
      for (int i = 0; i < nDark; ++i) addPulse(-tail + (n + tail) * u01(rng_), drawPe(st) * g * spread(rng_));
      if (muon) {
        const double x = pos[st][ci.plane == 'X' ? 0 : 1];
        if (x >= 0 && x < 16 && static_cast<int>(x) == ci.strip && u01(rng_) < p_.muonEfficiency) {
          std::poisson_distribution<int> np(p_.muonPeMean[st]);
          const int pe = std::max(1, np(rng_));
          addPulse(p_.signalSample + (u01(rng_) - 0.5), pe * g);
        }
      }
    }
    std::normal_distribution<double> noise(0.0, p_.noiseMv[st] / lsb);
    std::uint16_t* out = ev.adc.data() + static_cast<std::size_t>(c * n);
    for (int k = 0; k < n; ++k) {
      const double v = std::round(wf[static_cast<std::size_t>(k)] + noise(rng_));
      out[k] = static_cast<std::uint16_t>(std::min(maxAdc, std::max(0.0, v)));
    }
  }
  return ev;
}

Result<std::size_t> FakeDigitizer::readBlock(EventBlock& out, std::size_t maxEvents) {
  // Tras stop() solo se vacía la memoria (como ReadData después de SWStopAcquisition).
  if (!configured_) return Result<std::size_t>::failure(-1, "FakeDigitizer: readBlock() sin configure()");
  const UtcNs now = clock_.now_utc_ns();
  if (running_) generateUntil(now);
  out.tReadUtcNs = now;
  out.tReadMonoNs = monotonicNowNs();
  out.recordLength = cfg_.recordLength;
  out.boardFull = out.boardFull || fullSinceLastRead_;
  fullSinceLastRead_ = false;
  std::size_t n = 0;
  while (n < maxEvents && !memory_.empty()) {
    const Pending pd = memory_.front();
    memory_.pop_front();
    out.events.push_back(synthesize(pd));
    ++n;
  }
  return n;
}

nlohmann::json FakeDigitizer::info() const {
  return {{"model", "FakeDigitizer (sintético)"},
          {"serial", "none"},
          {"roc_firmware", "unknown"},
          {"amc_firmware", "unknown"},
          {"caen_digitizer_version", "unknown"},
          {"caen_comm_version", "unknown"},
          {"caen_vme_version", "unknown"},
          {"scenario", p_.toJson()}};
}

}  // namespace daq::dig
