// data_flags.cpp
#include "daq/daqrun/data_flags.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "daq/core/time.h"
#include "daq/trigger_interval_stats.h"

namespace daq::daqrun {

namespace {
constexpr double kNs = 1e9;
// P(N ≤ n | μ) = Q(n + 1, μ) (gamma regularizada superior).
double poissonLowerTail(std::uint64_t n, double mu) { return tis::detail::gammaQ(static_cast<double>(n) + 1.0, mu); }
double median(std::vector<double> v) {
  if (v.empty()) return -1;
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
std::array<double, 2> pair2(const nlohmann::json& j, const char* key, std::array<double, 2> def) {
  if (!j.contains(key)) return def;
  const auto& v = j[key];
  return {v.at(0).get<double>(), v.at(1).get<double>()};
}
}  // namespace

DataFlagConfig DataFlagConfig::fromJson(const nlohmann::json& j) {
  DataFlagConfig c;
  if (!j.is_object()) return c;
  c.enabled = j.value("enabled", c.enabled);
  c.windowS = j.value("window_s", c.windowS);
  c.blackoutWindowS = j.value("blackout_window_s", c.blackoutWindowS);
  c.imonWindowS = j.value("imon_window_s", c.imonWindowS);
  c.baselineHours = j.value("baseline_hours", c.baselineHours);
  c.baselineMinS = j.value("baseline_min_s", c.baselineMinS);
  c.minExpected = j.value("min_expected_counts", c.minExpected);
  c.minBaselineCounts = j.value("min_baseline_counts", c.minBaselineCounts);
  c.significanceP = j.value("significance_p", c.significanceP);
  c.blackoutRatio = j.value("blackout_ratio", c.blackoutRatio);
  c.lightLeakRatio = j.value("light_leak_ratio", c.lightLeakRatio);
  c.blackoutThresholdMv = pair2(j, "blackout_threshold_mV", c.blackoutThresholdMv);
  c.lightLeakThresholdMv = pair2(j, "light_leak_threshold_mV", c.lightLeakThresholdMv);
  c.clearEvaluations = j.value("clear_evaluations", c.clearEvaluations);
  c.imonK = j.value("imon_k_over_baseline", c.imonK);
  c.imonFracI0Set = j.value("imon_frac_i0set", c.imonFracI0Set);
  if (j.contains("baselines"))
    for (const auto& [k, v] : j["baselines"].items()) c.baselines[k] = v.get<double>();
  return c;
}

nlohmann::json DataFlagConfig::toJson() const {
  return {{"enabled", enabled},
          {"window_s", windowS},
          {"blackout_window_s", blackoutWindowS},
          {"imon_window_s", imonWindowS},
          {"baseline_hours", baselineHours},
          {"baseline_min_s", baselineMinS},
          {"min_expected_counts", minExpected},
          {"min_baseline_counts", minBaselineCounts},
          {"significance_p", significanceP},
          {"blackout_ratio", blackoutRatio},
          {"light_leak_ratio", lightLeakRatio},
          {"blackout_threshold_mV", blackoutThresholdMv},
          {"light_leak_threshold_mV", lightLeakThresholdMv},
          {"clear_evaluations", clearEvaluations},
          {"imon_k_over_baseline", imonK},
          {"imon_frac_i0set", imonFracI0Set},
          {"baselines", baselines},
          {"baseline_rule", "mediana de las ventanas de la unidad en las baseline_hours anteriores, terminadas antes de la "
                            "ventana actual y sin episodio activo; base explícita si está en baselines"}};
}

nlohmann::json UnitStatus::toJson() const {
  return {{"key", key},   {"name", name},     {"kind", kind},   {"channels", channels}, {"rate", rateHz},
          {"baseline", baselineHz >= 0 ? nlohmann::json(baselineHz) : nlohmann::json()},
          {"mu", mu},     {"n", n},           {"ratio", ratio}, {"p", p},               {"live", live},
          {"decided", decided}, {"active", active}, {"episode", episode}, {"since_utc_ns", since}};
}

nlohmann::json DataFlagView::toJson() const {
  nlohmann::json u = nlohmann::json::array();
  for (const auto& x : units) u.push_back(x.toJson());
  return {{"t_utc_ns", tUtcNs}, {"units", u}};
}

DataFlagEngine::DataFlagEngine(const DataFlagConfig& cfg, const hvmon::Mapping& mapping, std::string origin)
    : cfg_(cfg), origin_(std::move(origin)) {
  std::map<std::string, std::vector<int>> faces;
  std::map<int, std::vector<int>> groups;
  for (const auto& e : mapping.channels) {
    if (e.digCh < 0 || e.digCh >= 64) continue;
    if (!e.face.empty()) faces[e.station + " " + e.face].push_back(e.digCh);
    if (e.hvGroup) groups[*e.hvGroup].push_back(e.digCh);
  }
  for (auto& [name, chs] : faces) {
    Unit u;
    u.st = UnitStatus{"face:" + name, name, "face", chs};
    faces_.push_back(u);
  }
  for (auto& [g, chs] : groups) {
    Unit u;
    u.st = UnitStatus{"group:" + std::to_string(g), "grupo " + std::to_string(g), "group", chs};
    groups_.push_back(u);
  }
  for (int h = 0; h < 8; ++h) {
    Unit u;
    char name[16];
    std::snprintf(name, sizeof name, "CH%02d", h);
    u.st = UnitStatus{std::string("imon:") + name, name, "imon", {h}};
    imon_.push_back(u);
  }
}

double DataFlagEngine::baselineOf(Unit& u, UtcNs windowStart, double windowS, double* spanS) const {
  if (auto it = cfg_.baselines.find(u.st.key); it != cfg_.baselines.end()) {
    if (spanS) *spanS = 1e12;  // base explícita: sin límite de cuentas
    return it->second;
  }
  const UtcNs keepAfter = windowStart - static_cast<UtcNs>(cfg_.baselineHours * 3600 * kNs);
  while (!u.windows.empty() && u.windows.front().first < keepAfter) u.windows.pop_front();
  std::vector<double> v;
  UtcNs first = 0, last = 0;
  for (const auto& [t, r] : u.windows)
    if (t <= windowStart) {
      if (v.empty()) first = t;
      last = t;
      v.push_back(r);
    }
  // Historia mínima: el tramo cubierto por las ventanas usables (inicio de la primera a fin de la última).
  if (v.empty() || static_cast<double>(last - first) / kNs + windowS < cfg_.baselineMinS) return -1;
  if (spanS) *spanS = static_cast<double>(last - first) / kNs + windowS;
  return median(v);
}

hvmon::NotifyEvent DataFlagEngine::event(const std::string& kind, const Unit& u, UtcNs t, const std::string& reason) const {
  hvmon::NotifyEvent e;
  e.kind = kind;
  // Canal HV del asunto: grupo g -> CH(g−1); IMon: su canal; cara: ninguno.
  e.ch = u.st.kind == "face" ? -1 : (u.st.kind == "group" ? std::stoi(u.st.key.substr(6)) - 1 : u.st.channels.front());
  e.episode = u.st.episode;
  e.tUtcNs = t;
  e.reason = reason;
  e.origin = origin_;
  e.context = {{"unit", u.st.key}, {"unit_name", u.st.name}, {"rate", u.st.rateHz}, {"baseline", u.st.baselineHz},
               {"n", u.st.n}, {"mu", u.st.mu}, {"ratio", u.st.ratio}, {"p", u.st.p}};
  return e;
}

std::vector<hvmon::NotifyEvent> DataFlagEngine::addPeriod(const MonitorSummary& s, const HvSample* hv) {
  std::vector<hvmon::NotifyEvent> out;
  Period p;
  p.t0 = s.tStartUtcNs;
  p.t1 = s.tEndUtcNs;
  p.blk = s.blackoutCounts;
  p.leak = s.leakCounts;
  p.low = s.lowCounts;
  p.time = s.ootTimeS;
  if (hv && hv->alive && !hv->stale)
    for (int h = 0; h < 12; ++h) {
      const auto& c = hv->ch[static_cast<std::size_t>(h)];
      p.on[static_cast<std::size_t>(h)] = c.present && c.pw == 1;
      p.imon[static_cast<std::size_t>(h)] = c.imonA;
      p.i0set[static_cast<std::size_t>(h)] = c.i0setA;
    }
  periods_.push_back(p);
  // Se guardan los períodos de la ventana más larga y uno más (los períodos reales duran algo más
  // de summary_period_s: cierran con la lectura siguiente). Un período está en la ventana si al
  // menos la mitad de él cae dentro; la ventana está completa si esos períodos cubren 95 % de w.
  const double longest = std::max({cfg_.windowS, cfg_.blackoutWindowS, cfg_.imonWindowS});
  const UtcNs plen = p.t1 - p.t0;
  while (periods_.size() > 1 && periods_[1].t0 < p.t1 - static_cast<UtcNs>(longest * kNs) - plen) periods_.pop_front();
  auto inWindow = [](const Period& q, UtcNs ws) { return q.t1 - (q.t1 - q.t0) / 2 > ws; };
  auto windowOf = [&](double w) {
    const UtcNs ws = p.t1 - static_cast<UtcNs>(w * kNs);
    UtcNs first = p.t0;
    for (const auto& q : periods_)
      if (inWindow(q, ws)) {
        first = q.t0;
        break;
      }
    return std::make_pair(ws, static_cast<double>(p.t1 - first) / kNs >= 0.95 * w);
  };

  bits_ = DataFlagBits{};
  bits_.provisional = false;
  view_.tUtcNs = p.t1;
  view_.units.clear();

  auto evaluate = [&](Unit& u, bool leak) {
    const double wLen = leak ? cfg_.windowS : cfg_.blackoutWindowS;
    const auto [wStart, fullWindow] = windowOf(wLen);
    std::uint64_t n = 0;
    double T = 0;
    int live = 0;
    for (int c : u.st.channels) {
      std::uint64_t low = 0;
      for (const auto& q : periods_) {
        if (!inWindow(q, wStart)) continue;
        n += (leak ? q.leak : q.blk)[static_cast<std::size_t>(c)];
        T += q.time[static_cast<std::size_t>(c)];
        low += q.low[static_cast<std::size_t>(c)];
      }
      live += low > 0;
    }
    u.st.n = n;
    u.st.live = live;
    u.st.rateHz = T > 0 ? static_cast<double>(n) / T : 0;
    double span = 0;
    u.st.baselineHz = fullWindow ? baselineOf(u, wStart, wLen, &span) : -1;
    u.st.mu = u.st.baselineHz > 0 ? u.st.baselineHz * T : 0;
    u.st.ratio = u.st.mu > 0 ? static_cast<double>(n) / u.st.mu : 0;
    // μ ≥ mínimo; para el exceso (bit 7) también n ≥ mínimo con una base de suficientes cuentas.
    const bool baseCounts = u.st.mu * span / wLen >= cfg_.minBaselineCounts;
    u.st.decided = fullWindow && u.st.baselineHz >= 0 &&
                   (u.st.mu >= cfg_.minExpected || (leak && u.st.mu > 0 && static_cast<double>(n) >= cfg_.minExpected && baseCounts));
    bool cond = false;
    if (u.st.decided) {
      u.st.p = leak ? tis::poissonUpperTail(static_cast<std::size_t>(n), u.st.mu) : poissonLowerTail(n, u.st.mu);
      cond = leak ? (u.st.ratio > cfg_.lightLeakRatio && u.st.p < cfg_.significanceP)
                  : (u.st.ratio < cfg_.blackoutRatio && u.st.p < cfg_.significanceP);
    } else {
      u.st.p = 1;
    }
    if (cond && !u.st.active) {
      u.st.active = true;
      ++u.st.episode;
      u.st.since = p.t1;
      char why[256];
      std::snprintf(why, sizeof why, "%s: %llu pulsos en %.0f s frente a %.1f esperados (razón %.2f, p = %.2g), tasa %.1f Hz por canal, base %.1f Hz",
                    u.st.name.c_str(), static_cast<unsigned long long>(n), wLen, u.st.mu, u.st.ratio, u.st.p, u.st.rateHz,
                    u.st.baselineHz);
      out.push_back(event(leak ? "light_leak" : "data_blackout", u, p.t1, why));
    }
    // Sin decisión (p. ej. tras un hueco de trigger) el episodio sigue como estaba.
    if (cond) {
      u.falseRun = 0;
    } else if (u.st.decided && u.st.active && ++u.falseRun >= cfg_.clearEvaluations) {
      u.st.active = false;
      u.falseRun = 0;
    }
    // Historia para la base: solo ventanas completas y sin episodio activo.
    if (fullWindow && !u.st.active) u.windows.emplace_back(p.t1, u.st.rateHz);
    const std::uint16_t bit = leak ? chflag::LightLeak : chflag::DataBlackout;
    for (int c : u.st.channels) {
      auto& b = bits_.ch[static_cast<std::size_t>(c)];
      if (u.st.active) b |= bit;
      if (!u.st.decided) b |= chflag::StatInsufficient;
    }
    if (u.st.baselineHz < 0) bits_.provisional = true;
    auto& h = view_.history[u.st.key];
    h.push_back({p.t1, u.st.rateHz, u.st.baselineHz, u.st.active});
    if (h.size() > 1440) h.erase(h.begin());
    view_.units.push_back(u.st);
  };
  for (auto& u : faces_) evaluate(u, true);
  for (auto& u : groups_) evaluate(u, false);

  // IMon por canal HV.
  const auto [iStart, iFull] = windowOf(cfg_.imonWindowS);
  for (auto& u : imon_) {
    const auto h = static_cast<std::size_t>(u.st.channels.front());
    double sum = 0, i0 = 0;
    int n = 0;
    for (const auto& q : periods_)
      if (inWindow(q, iStart) && q.on[h]) {
        sum += q.imon[h];
        i0 = q.i0set[h];
        ++n;
      }
    const bool on = n > 0 && p.on[h];
    u.st.rateHz = n > 0 ? sum / n : 0;
    u.st.n = static_cast<std::uint64_t>(n);
    u.st.baselineHz = iFull && on ? baselineOf(u, iStart, cfg_.imonWindowS, nullptr) : -1;
    u.st.ratio = u.st.baselineHz > 0 ? u.st.rateHz / u.st.baselineHz : 0;
    u.st.mu = i0;  // I0Set (A), para el panel
    const bool overK = u.st.baselineHz > 0 && u.st.ratio > cfg_.imonK;
    const bool overI0 = on && i0 > 0 && u.st.rateHz > cfg_.imonFracI0Set * i0;
    u.st.decided = on && (u.st.baselineHz > 0 || i0 > 0);
    const bool cond = on && (overK || overI0);
    if (cond && !u.st.active) {
      u.st.active = true;
      ++u.st.episode;
      u.st.since = p.t1;
      char why[256];
      std::snprintf(why, sizeof why, "%s: IMon medio %.3f µA en %.0f s, base %.3f µA (razón %.2f), I0Set %.1f µA (%.0f %%)",
                    u.st.name.c_str(), u.st.rateHz * 1e6, cfg_.imonWindowS, u.st.baselineHz * 1e6, u.st.ratio, i0 * 1e6,
                    i0 > 0 ? 100.0 * u.st.rateHz / i0 : 0.0);
      out.push_back(event("imon_warning", u, p.t1, why));
    }
    if (cond) {
      u.falseRun = 0;
    } else if (u.st.active && ++u.falseRun >= cfg_.clearEvaluations) {
      u.st.active = false;
      u.falseRun = 0;
    }
    if (iFull && on && !u.st.active) u.windows.emplace_back(p.t1, u.st.rateHz);
    auto& hist = view_.history[u.st.key];
    hist.push_back({p.t1, u.st.rateHz, u.st.baselineHz, u.st.active});
    if (hist.size() > 1440) hist.erase(hist.begin());
    view_.units.push_back(u.st);
  }
  if (!cfg_.enabled) bits_ = DataFlagBits{std::array<std::uint16_t, 64>{}, false};
  return cfg_.enabled ? out : std::vector<hvmon::NotifyEvent>{};
}

}  // namespace daq::daqrun
