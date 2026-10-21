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
                          const std::string& reason, UtcNs now) {
  if (active_) return Status::failure(-1, "ya hay una sesión de barrido activa (restaurar primero)");
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
  if (target > nom + 1e-6) return Status::failure(-1, "V0Set " + fmt(target) + " V sobre el nominal " + fmt(nom) + " V");
  if (target < nom - lim_.maxBelowV - 1e-6)
    return Status::failure(-1, "V0Set " + fmt(target) + " V más de " + fmt(lim_.maxBelowV) + " V bajo el nominal " + fmt(nom) + " V");
  if (std::fabs(target - currentV0Set) > lim_.maxStepV + 1e-6)
    return Status::failure(-1, "escalón " + fmt(target - currentV0Set) + " V mayor que " + fmt(lim_.maxStepV) + " V");
  return okStatus();
}

bool BiasSession::reached(int ch, double target, double v0setRead, double vmon) const {
  if (!active_ || !inSession(ch)) return false;
  return std::fabs(v0setRead - target) <= lim_.matchTolV && std::fabs(vmon - v0setRead - offset_.at(ch)) <= lim_.vmonTolV;
}

std::vector<std::pair<int, double>> BiasSession::restorePlan(const std::map<int, double>& v0set) const {
  std::vector<std::pair<int, double>> plan;
  for (int ch : lim_.channels) {
    auto it = v0set.find(ch);
    if (it == v0set.end() || !nominal_.count(ch)) continue;
    double v = it->second;
    const double nom = nominal_.at(ch);
    while (nom - v > lim_.matchTolV) {
      v = std::min(nom, v + lim_.maxStepV);
      plan.push_back({ch, v});
    }
  }
  return plan;
}

nlohmann::json BiasSession::toJson(UtcNs now) const {
  nlohmann::json nom = nlohmann::json::object(), off = nlohmann::json::object();
  for (const auto& [ch, v] : nominal_) nom[std::to_string(ch)] = v;
  for (const auto& [ch, v] : offset_) off[std::to_string(ch)] = v;
  return {{"schema", "daqv2.hvmon_bias_session/1"}, {"active", active_}, {"nominal_V0Set_V", nom}, {"vmon_offset_V", off},
          {"start_utc_ns", startUtc_}, {"start_utc", startUtc_ ? formatUtcIso(startUtc_) : ""}, {"reason", reason_},
          {"saved_utc", formatUtcIso(now)},
          {"limits", {{"max_below_V", lim_.maxBelowV}, {"max_step_V", lim_.maxStepV}, {"vmon_tol_V", lim_.vmonTolV}, {"match_tol_V", lim_.matchTolV}}}};
}

Status BiasSession::restore(const nlohmann::json& j) {
  try {
    if (j.value("schema", "") != "daqv2.hvmon_bias_session/1") return Status::failure(-1, "estado de barrido: schema distinto");
    nominal_.clear();
    offset_.clear();
    for (const auto& [k, v] : j.at("nominal_V0Set_V").items()) nominal_[std::stoi(k)] = v.get<double>();
    for (const auto& [k, v] : j.at("vmon_offset_V").items()) offset_[std::stoi(k)] = v.get<double>();
    startUtc_ = j.value("start_utc_ns", static_cast<UtcNs>(0));
    reason_ = j.value("reason", "");
    active_ = j.value("active", false);
    for (int ch : lim_.channels)
      if (active_ && (!nominal_.count(ch) || !offset_.count(ch))) return Status::failure(-1, "estado de barrido incompleto");
  } catch (const std::exception& e) {
    return Status::failure(-1, std::string("estado de barrido ilegible: ") + e.what());
  }
  return okStatus();
}

}  // namespace daq::hvmon
