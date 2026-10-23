// bias_scan.cpp
#include "daq/hvmon/bias_scan.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace daq::hvmon {

namespace {
std::string fmt(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%.3f", v);
  return b;
}
}  // namespace

bool BiasSession::inSession(int ch) const { return std::find(lim_.channels.begin(), lim_.channels.end(), ch) != lim_.channels.end(); }

Status BiasSession::begin(const std::map<int, double>& expected, const std::map<int, double>& v0set, const std::map<int, double>& vmon,
                          const std::string& reason, UtcNs now, double maxAboveV, double maxStepV) {
  if (active_) return Status::failure(-1, "ya hay una sesión de barrido activa (restaurar primero)");
  if (!std::isfinite(maxAboveV) || maxAboveV < 0 || maxAboveV > kMaxAboveCapV + 1e-9)
    return Status::failure(-1, "max_above_V " + fmt(maxAboveV) + " V fuera de 0-" + fmt(kMaxAboveCapV) + " V");
  if (maxStepV > 0 && maxStepV > kMaxStepCapV + 1e-9)
    return Status::failure(-1, "max_step_V " + fmt(maxStepV) + " V mayor que " + fmt(kMaxStepCapV) + " V");
  std::map<int, double> nom, off;
  for (int ch : lim_.channels) {
    auto e = expected.find(ch), v = v0set.find(ch), m = vmon.find(ch);
    if (e == expected.end()) return Status::failure(-1, "falta el V0Set nominal esperado de CH" + std::to_string(ch));
    if (v == v0set.end() || m == vmon.end()) return Status::failure(-1, "sin lectura de V0Set/VMon de CH" + std::to_string(ch));
    if (std::fabs(v->second - e->second) > lim_.matchTolV)
      return Status::failure(-1, "CH" + std::to_string(ch) + ": V0Set leído " + fmt(v->second) + " V distinto del nominal esperado " +
                                     fmt(e->second) + " V: no se empieza");
    nom[ch] = v->second;
    off[ch] = m->second - v->second;
  }
  nominal_ = nom;
  offset_ = off;
  offsetNominal_ = off;
  lim_.maxAboveV = maxAboveV;
  lim_.maxStepV = maxStepV > 0 ? maxStepV : baseMaxStepV_;
  startUtc_ = now;
  reason_ = reason;
  active_ = true;
  return okStatus();
}

Status BiasSession::checkStep(int ch, double target, double currentV0Set) const {
  if (!active_) return Status::failure(-1, "sin sesión de barrido activa: V0Set no se escribe");
  if (!inSession(ch)) return Status::failure(-1, "CH" + std::to_string(ch) + " no es parte de la sesión de barrido");
  if (!std::isfinite(target)) return Status::failure(-1, "V0Set pedido inválido");
  const double nom = nominal_.at(ch);
  if (target > nom + lim_.maxAboveV + 1e-6)
    return Status::failure(-1, "V0Set " + fmt(target) + " V sobre el nominal " + fmt(nom) + " V + " + fmt(lim_.maxAboveV) + " V (límite de la sesión)");
  if (target < nom - lim_.maxBelowV - 1e-6)
    return Status::failure(-1, "V0Set " + fmt(target) + " V más de " + fmt(lim_.maxBelowV) + " V bajo el nominal " + fmt(nom) + " V");
  if (std::fabs(target - currentV0Set) > lim_.maxStepV + 1e-6)
    return Status::failure(-1, "escalón " + fmt(target - currentV0Set) + " V mayor que " + fmt(lim_.maxStepV) + " V");
  if (std::max(target, currentV0Set) > nom + 1e-6 && std::fabs(target - currentV0Set) > kMaxStepAboveV + 1e-6)
    return Status::failure(-1, "escalón " + fmt(target - currentV0Set) + " V sobre el nominal mayor que " + fmt(kMaxStepAboveV) + " V");
  return okStatus();
}

bool BiasSession::reached(int ch, double target, double v0setRead, double vmon, double v0setBefore, double vmonBefore) const {
  if (!active_ || !inSession(ch)) return false;
  if (std::fabs(v0setRead - target) > lim_.matchTolV) return false;
  if (std::fabs(vmon - v0setRead - offset_.at(ch)) > lim_.vmonTolV) return false;
  const double step = target - v0setBefore, moved = vmon - vmonBefore;
  return step * moved > 0 && std::fabs(moved) >= lim_.minMoveFrac * std::fabs(step);
}

bool BiasSession::offsetDriftOk(int ch, double v0setRead, double vmon) const {
  return active_ && inSession(ch) && offsetNominal_.count(ch) && std::fabs(vmon - v0setRead - offsetNominal_.at(ch)) <= lim_.maxOffsetDriftV + 1e-9;
}

bool BiasSession::withinTolerance(int ch, double v0setRead, double vmon) const {
  return active_ && inSession(ch) && std::fabs(vmon - v0setRead - offset_.at(ch)) <= lim_.vmonTolV;
}

std::vector<std::pair<int, double>> BiasSession::restorePlan(const std::map<int, double>& v0set) const {
  std::vector<std::pair<int, double>> plan;
  for (int ch : lim_.channels) {
    auto it = v0set.find(ch);
    if (it == v0set.end() || !nominal_.count(ch)) continue;
    double v = it->second;
    const double nom = nominal_.at(ch);
    while (nom - v > lim_.matchTolV) {  // desde abajo
      v = std::min(nom, v + lim_.maxStepV);
      plan.push_back({ch, v});
    }
    while (v - nom > lim_.matchTolV) {  // desde arriba (DAQ2-17): escalones ≤ kMaxStepAboveV
      v = std::max(nom, v - std::min(lim_.maxStepV, kMaxStepAboveV));
      plan.push_back({ch, v});
    }
  }
  return plan;
}

nlohmann::json BiasSession::toJson(UtcNs now) const {
  nlohmann::json nom = nlohmann::json::object(), off = nlohmann::json::object(), off0 = nlohmann::json::object();
  for (const auto& [ch, v] : nominal_) nom[std::to_string(ch)] = v;
  for (const auto& [ch, v] : offset_) off[std::to_string(ch)] = v;
  for (const auto& [ch, v] : offsetNominal_) off0[std::to_string(ch)] = v;
  return {{"schema", "daqv2.hvmon_bias_session/1"}, {"active", active_}, {"nominal_V0Set_V", nom}, {"vmon_offset_V", off}, {"vmon_offset_nominal_V", off0},
          {"start_utc_ns", startUtc_}, {"start_utc", startUtc_ ? formatUtcIso(startUtc_) : ""}, {"reason", reason_},
          {"saved_utc", formatUtcIso(now)},
          {"limits", {{"max_below_V", lim_.maxBelowV}, {"max_step_V", lim_.maxStepV}, {"max_above_V", lim_.maxAboveV}, {"vmon_tol_V", lim_.vmonTolV},
                      {"match_tol_V", lim_.matchTolV},
                      {"min_move_frac", lim_.minMoveFrac}, {"ramp_tau_s", lim_.rampTauS}, {"opposite_V", lim_.oppositeV},
                      {"max_offset_drift_V", lim_.maxOffsetDriftV}}}};
}

Status BiasSession::restore(const nlohmann::json& j) {
  try {
    if (j.value("schema", "") != "daqv2.hvmon_bias_session/1") return Status::failure(-1, "estado de barrido: schema distinto");
    nominal_.clear();
    offset_.clear();
    for (const auto& [k, v] : j.at("nominal_V0Set_V").items()) nominal_[std::stoi(k)] = v.get<double>();
    for (const auto& [k, v] : j.at("vmon_offset_V").items()) offset_[std::stoi(k)] = v.get<double>();
    offsetNominal_ = offset_;  // estado anterior a DAQ2-17 sin el desvío al nominal: el último verificado
    if (j.contains("vmon_offset_nominal_V"))
      for (const auto& [k, v] : j.at("vmon_offset_nominal_V").items()) offsetNominal_[std::stoi(k)] = v.get<double>();
    startUtc_ = j.value("start_utc_ns", static_cast<UtcNs>(0));
    reason_ = j.value("reason", "");
    active_ = j.value("active", false);
    // DAQ2-17: límites de la sesión persistida (una sesión sobre el nominal se restaura con su escalón), con los topes fijos.
    const auto lim = j.value("limits", nlohmann::json::object());
    lim_.maxAboveV = std::clamp(lim.value("max_above_V", 0.0), 0.0, kMaxAboveCapV);
    lim_.maxStepV = std::clamp(lim.value("max_step_V", baseMaxStepV_), 0.01, kMaxStepCapV);
    for (int ch : lim_.channels)
      if (active_ && (!nominal_.count(ch) || !offset_.count(ch))) return Status::failure(-1, "estado de barrido incompleto");
  } catch (const std::exception& e) {
    return Status::failure(-1, std::string("estado de barrido ilegible: ") + e.what());
  }
  return okStatus();
}

}  // namespace daq::hvmon
