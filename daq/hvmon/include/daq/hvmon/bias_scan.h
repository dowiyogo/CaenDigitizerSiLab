// bias_scan.h
// Sesión de barrido de tensión de hvmon (DAQ2-14): la única vía para escribir V0Set.
//  - begin(): registra el V0Set nominal de los canales de la sesión (CH00-07) leído de la placa y lo
//    compara con el esperado (si difiere en más de match_tol_V, no empieza); guarda también el
//    desvío VMon − V0Set de cada canal al nominal (en la placa real −0.14 a −0.38 V, DAQ2-6), con el
//    que se verifica la rampa: el criterio "VMon a 0.05 V del V0Set nuevo" se aplica al desvío.
//  - checkStep(): límites duros de cada escritura: solo canales de la sesión; V0Set ≤ nominal + max_above_V;
//    V0Set ≥ nominal − max_below_V (2.0 V); |V0Set nuevo − V0Set actual| ≤ max_step_V (0.5 V).
//  - DAQ2-17 (2026-10-23, barrido hacia arriba autorizado por René): max_above_V y max_step_V se fijan por sesión en
//    begin(); max_above_V va de 0 (por omisión: nunca sobre el nominal) a kMaxAboveCapV = 1.0 V, tope fijo en el código;
//    toda escritura con V0Set nuevo o actual sobre el nominal tiene además un escalón ≤ kMaxStepAboveV = 0.25 V.
//  - reached(): V0Set releído = pedido (± match_tol_V), |VMon − V0Set − desvío| ≤ vmon_tol_V (0.10 V) y
//    VMon se movió en el sentido del escalón al menos la mitad del escalón (VMon − VMon antes de
//    escribir, con el signo de V0Set pedido − V0Set antes). El desvío de referencia es el del último
//    escalón verificado (updateOffset): el desvío de la placa cambia ~±0.045 V en un día y ≤ 0.02 V/µA
//    con la corriente (2026-10-20). Criterio aceptado por René (2026-10-21).
//  - restorePlan(): escalones de ≤ max_step_V hacia el nominal, canal por canal, desde abajo o desde arriba
//    (sobre el nominal, ≤ kMaxStepAboveV).
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

inline constexpr double kMaxAboveCapV = 1.0;    // tope fijo de V0Set sobre el nominal (DAQ2-17)
inline constexpr double kMaxStepAboveV = 0.25;  // escalón máximo con V0Set sobre el nominal (DAQ2-17)
inline constexpr double kMaxStepCapV = 0.5;     // escalón máximo en cualquier sesión

struct BiasLimits {
  std::vector<int> channels{0, 1, 2, 3, 4, 5, 6, 7};
  double maxBelowV = 2.0;
  double maxStepV = 0.5;
  double maxAboveV = 0.0;  // por sesión (begin); 0 = solo en o bajo el nominal
  double vmonTolV = 0.10;
  double minMoveFrac = 0.5;  // VMon debe moverse al menos esta fracción del escalón, en su sentido
  // DAQ2-14 (2026-10-23): la rampa se espera hasta 6τ (τ ≈ 115 s, descarga del PCB medida en DAQ2-10);
  // VMon moviéndose en contra del escalón más de opposite_V es un riesgo (aborto).
  double rampTauS = 115.0;
  double oppositeV = 0.10;
  double matchTolV = 0.005;
};

class BiasSession {
 public:
  explicit BiasSession(BiasLimits lim = {}) : lim_(std::move(lim)) {}
  // expected: V0Set nominal esperado por canal; v0set/vmon: última lectura.
  // maxAboveV: V0Set máximo sobre el nominal en esta sesión (0 a kMaxAboveCapV); maxStepV > 0: escalón de la sesión
  // (≤ kMaxStepCapV); ≤ 0: el de la construcción.
  Status begin(const std::map<int, double>& expected, const std::map<int, double>& v0set, const std::map<int, double>& vmon,
               const std::string& reason, UtcNs now, double maxAboveV = 0.0, double maxStepV = -1.0);
  Status checkStep(int ch, double target, double currentV0Set) const;
  bool reached(int ch, double target, double v0setRead, double vmon, double v0setBefore, double vmonBefore) const;
  bool withinTolerance(int ch, double v0setRead, double vmon) const;  // solo |VMon − V0Set − desvío| ≤ vmon_tol_V
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
  double baseMaxStepV_ = lim_.maxStepV;
  bool active_ = false;
  std::map<int, double> nominal_, offset_;
  UtcNs startUtc_ = 0;
  std::string reason_;
};

}  // namespace daq::hvmon
