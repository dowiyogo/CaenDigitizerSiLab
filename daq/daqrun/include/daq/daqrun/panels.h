// panels.h
// Paneles de monitoreo en tiempo de daqrun (fase 5, DAQ2-9; diseño 7.2). Datos planos, sin ROOT: los
// llena el hilo monitor (Panels, dueño Monitor) y el hilo principal los convierte en histogramas
// (daq/io/panels_root.h). Cada histograma de eventos tiene dos versiones: acumulada del run y la
// última hora (60 rebanadas de summary_period_s).
//  - amplitud en tiempo por canal: línea base − mínimo en la región signal (mV; p.e. cuando haya
//    calibración);
//  - integral en tiempo por canal: Σ (línea base − x)·16 ns en signal (mV·ns, equivalente QDC);
//  - muestra del mínimo del registro por canal, en eventos con amplitud > hit_threshold_mV;
//  - mapas X × Y por estación (tiras en orden de coordenada del mapeo) con amplitud en tiempo >
//    hit_threshold_mV; matriz de coincidencias entre tiras de un plano (16 × 16); multiplicidad por
//    plano en dos umbrales;
//  - persistencia por grupo del digitizador (muestra × ADC de sus 8 canales, últimos ~1000 eventos:
//    bloque anterior + bloque en curso de persistence_events/2) y el último evento (canal × muestra);
//  - mapa canal × tiempo de la tasa fuera de tiempo (umbral provisional) y canales vivos por plano
//    contra el tiempo (> 2 LSB en el período), últimas 24 h;
//  - G_c por hora de hod2018 (gain_fit.h: semilla de la hora anterior, rechazo de G/2 y 2G por la
//    separación de picos, ambiguos marcados) y tasa oscura > 3 mV por hora de hod2019 (indicador de
//    ganancia). PROVISIONAL mientras el run esté sin calibración.
#pragma once

#include <array>
#include <cstdint>
#include <deque>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/time.h"
#include "daq/dig/digitizer.h"
#include "daq/hvmon/mapping.h"

namespace daq::daqrun {

struct DaqrunConfig;
struct MonitorSummary;

// Histograma plano (x más rápido) con versión acumulada y de la última hora.
struct PanelHist {
  std::string name, title;
  int nx = 0, ny = 1;
  double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  std::vector<std::uint64_t> cum;
  std::vector<std::uint32_t> hour;
};

class RollingHist {
 public:
  void init(std::string name, std::string title, int nx, double x0, double x1, int ny = 1, double y0 = 0, double y1 = 1,
            std::size_t slices = 60);
  void fill(double x, double y = 0.5);
  void fillBin(int ix, int iy = 0);
  void closeSlice();  // fin de período: la rebanada más vieja sale de la última hora
  const PanelHist& view() const { return h_; }

 private:
  PanelHist h_;
  std::size_t slices_ = 60;
  std::vector<std::uint32_t> cur_;
  std::deque<std::vector<std::uint32_t>> old_;
};

struct GainPoint {
  UtcNs hourStartUtcNs = 0;
  int ch = 0;
  double value = 0;       // hod2018: G (mV/p.e.); hod2019: tasa > 3 mV (Hz)
  std::string status;     // hod2018: ok | changed | ambiguous | rejected | insufficient; hod2019: "dark_rate"
  int peaks = 0;
};

struct PanelView {
  UtcNs tUtcNs = 0;
  std::uint64_t events = 0;
  int recordLength = 0;
  std::vector<PanelHist> hists;                      // amplitud, integral, mínimo, mapas, coincidencias, multiplicidad
  int persistAdcBins = 256;
  std::array<std::vector<std::uint32_t>, 8> persistence;  // [muestra][ADC/16]
  std::vector<std::uint16_t> lastEvent;              // [canal][muestra]
  UtcNs lastEventUtcNs = 0;
  std::vector<UtcNs> rateMapT;                       // fin de cada período (últimas 24 h)
  std::vector<std::array<float, 64>> rateMap;        // Hz por canal y período
  std::vector<std::array<int, 4>> livePerPlane;      // por período: hod2018 X, Y, hod2019 X, Y
  std::vector<GainPoint> gains;
  nlohmann::json gainSummary() const;                // última hora: G por canal y conteo por estado
  // Fracción de las coincidencias entre tiras (i ≠ j, acumulado del run) que son vecinas |i − j| = 1,
  // por orden (0 medido r3, 1 hoja de octubre) y plano (hod2018 X, Y, hod2019 X, Y); −1 sin pares.
  std::array<std::array<double, 4>, 2> neighbourFrac{{{-1, -1, -1, -1}, {-1, -1, -1, -1}}};
  std::array<std::array<double, 4>, 2> neighbourPairs{};
};

struct PanelConfig {
  double inTimeAmpMaxMv = 1000, inTimeAmpBinMv = 2;
  double integralMaxMvNs = 40000, integralBinMvNs = 100;
  std::array<double, 2> hitThresholdMv{50.0, 5.0};
  std::array<std::array<double, 2>, 2> multThresholdMv{{{20.0, 50.0}, {3.0, 5.0}}};
  int persistenceEvents = 1000;
  double gainPeriodS = 3600;
  double gainSeedMv = 4.0;  // vista previa del run 6 (mediana 3.978 mV); la hora siguiente usa la anterior
  static PanelConfig fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
};

class Panels {
 public:
  Panels(const DaqrunConfig& cfg, const hvmon::Mapping* mapping, int recordLength, double lsbMv);
  // Un evento: x de los 64 canales y sus líneas base (cuentas ADC).
  void event(const dig::Event& ev, const std::array<double, 64>& baselineAdc, UtcNs tUtc);
  // Fin de período: rebanadas, mapa de tasa, vivos por plano; espectros fuera de tiempo acumulados
  // (para la ganancia por hora).
  void period(const MonitorSummary& s, const std::vector<std::uint32_t>& ampCum, int ampBins, double ampBinMv);
  const PanelView& view() const { return v_; }
  void refresh() { publishHists(); }  // copia los histogramas a la vista (entre períodos)

 private:
  const DaqrunConfig& cfg_;
  int N_;
  double lsb_;
  std::array<int, 64> station_{}, plane_{}, strip_{}, stripSheet_{};  // plano 0-3 (hod2018 X, Y, hod2019 X, Y); −1 sin mapeo
  RollingHist amp_, integ_, minS_;
  // [orden][...]: 0 = coordenada medida (strip del mapeo, r3), 1 = la de la hoja de octubre (strip_sheet)
  std::array<std::array<RollingHist, 2>, 2> hit_;
  std::array<std::array<RollingHist, 4>, 2> coinc_;
  std::array<std::array<RollingHist, 2>, 4> mult_;
  std::array<std::vector<std::uint32_t>, 8> persPrev_, persCur_;
  int persCount_ = 0;
  PanelView v_;
  // Ganancia por hora.
  UtcNs hourStart_ = 0;
  std::vector<std::uint32_t> ampAtHour_;
  std::array<double, 64> seed_{};
  std::array<double, 64> darkCounts_{}, darkTime_{};
  void publishHists();
};

}  // namespace daq::daqrun
