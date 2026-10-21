// data_flags.h
// Bits de datos de daqrun (fase 5, DAQ2-9; diseño 6.1 bits 6, 7 y 11, 6.2 bit 6, 6.4) y aviso de IMon.
// Lógica pura, sin ROOT ni hilos: la alimenta el hilo monitor con un resumen por período
// (MonitorSummary, 60 s) y las muestras de hvmon; el escritor lee los bits vigentes (Snapshot).
//
// Unidades (del mapeo JSON):
//  - unidades de luz del mapeo (light_units, r3: canales expuestos a la luz de la tarde) o, si no hay,
//    caras físicas (estación + cara) para el bit 7, LIGHT_LEAK; hod2018 no tiene caras en el mapeo
//    de octubre, así que el bit 7 cubre hod2019 (la luz de mayo entró por la cara Trasera);
//  - grupos HV (1-8) para el bit 6, DATA_BLACKOUT; los canales de grupo desconocido no entran.
// Conteo: todos los pulsos de la zona sin señal (pedestal + control + post, más tail sin pulso en
// tiempo; diseño 7.1, nunca uno por registro) sobre un umbral por estación:
//  - bit 6: blackout_threshold_mV (3 mV hod2018, 1.5 mV hod2019: con HV apagada la tasa es 0 en la
//    mediana de los canales, run004; con HV encendida 1.6 MHz y 14 kHz de mediana);
//  - bit 7: light_leak_threshold_mV (20 mV hod2018, 3 mV hod2019, los provisionales).
// Ventanas: window_s (600 s) para las caras, blackout_window_s (300 s) para los grupos e
// imon_window_s (120 s) para IMon (una ventana de 10 min tarda ~9 min en bajar de 0.1 o en pasar
// 1.5×; ver los tests). Tasa de la unidad = Σ cuentas /
// Σ tiempo expuesto de sus canales (Hz por canal); μ = tasa_base × Σ tiempo expuesto.
// Tasa base (documentada en 15.7): la mediana de las tasas de ventana de la unidad en las
// baseline_hours (3 h) anteriores, solo ventanas que terminan antes de que empiece la actual y
// sin episodio activo; hace falta historia de baseline_min_s (30 min). Una base explícita en
// `baselines` (p. ej. la de la noche) la reemplaza.
//  - bit 6: μ ≥ min_expected, P(N ≤ n | μ) < significance_p y n/μ < blackout_ratio (0.1);
//  - bit 7: (μ ≥ min_expected, o n ≥ min_expected con una base de min_baseline_counts cuentas),
//    P(N ≥ n | μ) < significance_p y n/μ > light_leak_ratio (2);
//  - bit 11: la prueba de la unidad no se pudo hacer (sin base o μ < mínimo);
//  - bit 6 de evento (MONITOR_PROVISIONAL): alguna unidad sin base todavía.
// Episodio: empieza cuando se cumple la condición y termina tras clear_evaluations evaluaciones
// decididas seguidas sin cumplirla (una ventana sin decisión no lo cierra); un correo por unidad y episodio (notifier: light_leak, data_blackout,
// imon_warning; dedup por context.unit y episodio).
// IMon (por canal HV con Pw = 1): media de la ventana > imon_k × base (mediana de las medias por
// período de las horas anteriores, como arriba) o > imon_frac_i0set × I0Set.
#pragma once

#include <array>
#include <cstdint>
#include <deque>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/time.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/monitor.h"
#include "daq/hvmon/mapping.h"
#include "daq/hvmon/notifier.h"

namespace daq::daqrun {

// DataFlagConfig: en config.h (clave data_flags de la configuración de daqrun).

// Bits vigentes para el escritor.
struct DataFlagBits {
  std::array<std::uint16_t, 64> ch{};  // bits 6, 7 y 11 de ch_flags
  bool provisional = true;            // bit 6 de evt_flags
};

struct UnitStatus {
  std::string key, name, kind;  // kind: "face" | "group" | "imon"
  std::vector<int> channels;
  double rateHz = 0, baselineHz = -1, mu = 0, ratio = 0, p = 1;  // imon: rate = IMon medio (A), base en A
  std::uint64_t n = 0;
  int live = 0;  // canales con algún pulso > 2 LSB en la ventana
  bool decided = false, active = false;
  std::uint64_t episode = 0;
  UtcNs since = 0;  // inicio del episodio activo
  nlohmann::json toJson() const;
};

struct UnitPoint {
  UtcNs t;
  double rate, baseline;
  bool active;
};

struct DataFlagView {
  UtcNs tUtcNs = 0;
  std::vector<UnitStatus> units;                       // caras, grupos e IMon
  std::map<std::string, std::vector<UnitPoint>> history;  // por unidad, últimas 24 h
  nlohmann::json toJson() const;                        // sin la historia
};

class DataFlagEngine {
 public:
  DataFlagEngine(const DataFlagConfig& cfg, const hvmon::Mapping& mapping, std::string origin);
  // Un período cerrado del monitor y la última muestra de hvmon. Devuelve los correos a enviar.
  std::vector<hvmon::NotifyEvent> addPeriod(const MonitorSummary& s, const HvSample* hv);
  const DataFlagBits& bits() const { return bits_; }
  const DataFlagView& view() const { return view_; }
  const DataFlagConfig& config() const { return cfg_; }

 private:
  struct Unit {
    UnitStatus st;
    std::deque<std::pair<UtcNs, double>> windows;  // (fin de ventana, tasa) para la base
    int falseRun = 0;
  };
  struct Period {
    UtcNs t0, t1;
    std::array<std::uint32_t, 64> blk{}, leak{}, low{};
    std::array<double, 64> time{};
    std::array<double, 12> imon{}, i0set{};
    std::array<bool, 12> on{};
  };
  double baselineOf(Unit& u, UtcNs windowStart, double windowS, double* spanS) const;
  hvmon::NotifyEvent event(const std::string& kind, const Unit& u, UtcNs t, const std::string& reason) const;

  DataFlagConfig cfg_;
  std::string origin_;
  std::vector<Unit> faces_, groups_, imon_;
  std::deque<Period> periods_;
  DataFlagBits bits_;
  DataFlagView view_;
};

}  // namespace daq::daqrun
