// bias_scan.h
// Sesión de barrido de tensión de hvmon (DAQ2-14): la única vía para escribir V0Set.
//  - begin(): registra el V0Set nominal de los canales de la sesión (CH00-07) leído de la placa y lo
//    compara con el esperado (si difiere en más de match_tol_V, no empieza); guarda también el
//    desvío VMon − V0Set de cada canal al nominal (en la placa real −0.14 a −0.38 V, DAQ2-6), con el
//    que se verifica la rampa: el criterio "VMon a 0.05 V del V0Set nuevo" se aplica al desvío.
//  - checkStep(): límites duros de cada escritura: solo canales de la sesión; V0Set ≤ nominal;
//    V0Set ≥ nominal − max_below_V (2.0 V); |V0Set nuevo − V0Set actual| ≤ max_step_V (0.5 V).
//  - reached(): V0Set releído = pedido (± match_tol_V), |VMon − V0Set − desvío| ≤ vmon_tol_V (0.10 V) y
//    VMon se movió en el sentido del escalón al menos la mitad del escalón (VMon − VMon antes de
//    escribir, con el signo de V0Set pedido − V0Set antes). El desvío de referencia es el del último
//    escalón verificado (updateOffset): el desvío de la placa cambia ~±0.045 V en un día y ≤ 0.02 V/µA
//    con la corriente (2026-10-20). Criterio aceptado por René (2026-10-21).
//  - restorePlan(): escalones de ≤ max_step_V hacia el nominal, canal por canal.
//  - El estado se persiste (JSON) para que un hvmon que arranca con una sesión activa restaure el
//    nominal solo (Poller::enableBiasPersistence).
// Lógica pura, sin hardware ni hilos (la usa el hilo de sondeo; la prueban los tests).
#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/core/time.h"

namespace daq::hvmon {

struct BiasLimits {
  std::vector<int> channels{0, 1, 2, 3, 4, 5, 6, 7};
  double maxBelowV = 2.0;
  double maxStepV = 0.5;
  double vmonTolV = 0.10;
  double minMoveFrac = 0.5;  // VMon debe moverse al menos esta fracción del escalón, en su sentido
  double matchTolV = 0.005;
};

class BiasSession {
 public:
  explicit BiasSession(BiasLimits lim = {}) : lim_(std::move(lim)) {}
  // expected: V0Set nominal esperado por canal; v0set/vmon: última lectura.
  Status begin(const std::map<int, double>& expected, const std::map<int, double>& v0set, const std::map<int, double>& vmon,
               const std::string& reason, UtcNs now);
  Status checkStep(int ch, double target, double currentV0Set) const;
  bool reached(int ch, double target, double v0setRead, double vmon, double v0setBefore, double vmonBefore) const;
  std::vector<std::pair<int, double>> restorePlan(const std::map<int, double>& v0set) const;
  void end() { active_ = false; }
  void updateOffset(int ch, double off) {
    if (offset_.count(ch)) offset_[ch] = off;
  }
  bool active() const { return active_; }
  bool inSession(int ch) const;
  double nominal(int ch) const { return nominal_.at(ch); }
  double offset(int ch) const { return offset_.at(ch); }
  const BiasLimits& limits() const { return lim_; }
  nlohmann::json toJson(UtcNs now) const;
  Status restore(const nlohmann::json& j);

 private:
  BiasLimits lim_;
  bool active_ = false;
  std::map<int, double> nominal_, offset_;
  UtcNs startUtc_ = 0;
  std::string reason_;
};

}  // namespace daq::hvmon
