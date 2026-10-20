// monitor.h
// Hilo monitor de daqrun, fase 3 (diseño 7.1). Datos planos, sin ROOT: el hilo principal convierte
// la vista en histogramas de la página y el escritor guarda los resúmenes en el árbol `monitor`.
// Calcula solo:
//  - tasa de trigger (por período de resumen) y distribución de Δt (log10 Δt/s);
//  - por canal, pulsos fuera de tiempo sobre el umbral provisional (hod2018 20 mV, hod2019 3 mV:
//    estado sin calibración) en pedestal + control + post (+ tail si el canal no tiene pulso en
//    tiempo), su tiempo expuesto y su tasa; pulsos en tiempo (región signal);
//  - espectro de amplitud de los pulsos fuera de tiempo, en mV, paso de 1 LSB (0.4884 mV);
//  - niveles de cola (del proveedor que pasa main) y eventos descartados por el monitor.
// Línea base por evento y canal: media recortada (±trim_mV de la mediana) del pedestal.
// Seguridad entre hilos: process()/takeSummary() en el hilo monitor; view() se lee desde
// cualquier hilo (Snapshot).
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/snapshot.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/hv_client.h"
#include "daq/dig/digitizer.h"

namespace daq::daqrun {

struct MonitorSummary {
  UtcNs tStartUtcNs = 0, tEndUtcNs = 0;
  std::uint64_t nEvents = 0;
  double triggerRateHz = 0;
  std::array<std::uint32_t, 64> ootCounts{}, inTimeCounts{};
  std::array<double, 64> ootTimeS{}, ootRateHz{};
  std::array<std::uint32_t, 64> lowCounts{};  // pulsos fuera de tiempo > 2 LSB (canal vivo; no va al árbol)
  std::array<std::uint32_t, 64> blackoutCounts{}, leakCounts{};  // sobre los umbrales de data_flags (bits 6 y 7)
  nlohmann::json queues;  // niveles de cola al cerrar el período
  std::uint64_t monitorDropped = 0;
};

struct MonitorView {
  std::uint64_t events = 0;
  UtcNs firstEventUtcNs = 0, lastEventUtcNs = 0;
  std::array<std::uint64_t, 64> ootCounts{}, inTimeCounts{};
  std::array<double, 64> ootTimeS{};
  int ampBins = 0;
  double ampBinMv = 0;
  std::vector<std::uint32_t> amp;  // [64][ampBins]
  static constexpr int kDtBins = 70;
  static constexpr double kDtLogMin = -4.0, kDtLogMax = 3.0;  // log10(Δt/s)
  std::array<std::uint32_t, kDtBins> dtHist{};
  std::vector<std::pair<UtcNs, double>> rateHistory;  // (fin del período, tasa Hz)
  MonitorSummary last;
  std::vector<MonitorSummary> recent;  // últimos kRecent períodos (10 × 60 s): página /dash/
  static constexpr std::size_t kRecent = 10;
};

class Monitor {
 public:
  Monitor(const DaqrunConfig& cfg, UtcNs originUtc, std::function<nlohmann::json()> queueLevels = {});
  void process(const dig::EventBlock& blk);
  // Si pasó summary_period_s desde el inicio del período, cierra el período y lo devuelve.
  bool takeSummary(UtcNs now, MonitorSummary& out, std::uint64_t monitorDropped);
  void publishView();
  std::pair<MonitorView, std::uint64_t> view() const { return view_.read(); }

 private:
  const DaqrunConfig& cfg_;
  TttUnwrapper ttt_;
  std::function<nlohmann::json()> queues_;
  double lsb_;
  std::vector<int> regionsPed_;
  std::int64_t lastTttNs_ = -1;
  UtcNs origin_;
  MonitorSummary cur_;
  MonitorView v_;
  Snapshot<MonitorView> view_;
  std::vector<std::uint16_t> scratch_;
};

}  // namespace daq::daqrun
