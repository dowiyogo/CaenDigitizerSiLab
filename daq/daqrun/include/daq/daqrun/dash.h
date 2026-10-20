// dash.h
// Números de la página /dash/ de daqrun (objeto "dash" del JSON de estado). Función pura: la
// llama el hilo principal con la vista del monitor, la última muestra de hvmon y el disco libre.
//  - tasa de trigger del último período (1 min) y de los últimos MonitorView::kRecent (10 min);
//  - segundos desde el último trigger (t_utc del último evento);
//  - por estación: canales vivos = con al menos un pulso fuera de tiempo > 2 LSB en los últimos
//    10 min (hod2019 tiene canales con tasa oscura muy baja, por eso no se usa el umbral de 3 mV),
//    y mediana entre los vivos de la tasa fuera de tiempo sobre el umbral PROVISIONAL (hod2018
//    20 mV, hod2019 3 mV) en esos 10 min;
//  - avisos activos (texto) que la página muestra en rojo.
#pragma once

#include <array>
#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/time.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/monitor.h"

namespace daq::daqrun {

struct DashInputs {
  std::string campaign, runName;
  bool simulated = false;
  UtcNs runStartUtcNs = 0, nowUtcNs = 0;
  const MonitorView* view = nullptr;
  const HvSample* hv = nullptr;
  double diskFreeBytes = -1;            // < 0: desconocido
  double diskLowBytes = 20e9;           // aviso bajo este valor
  bool watchTriggers = false;           // run con trigger externo
  bool noTriggerAlert = false;          // TriggerWatchdog en episodio
  double noTriggerAlertS = 300;         // aviso directo si no hay trigger en este tiempo
  std::array<double, 2> thresholdMv{20.0, 3.0};
  nlohmann::json readout = nlohmann::json::object();  // ReadoutStats::toJson()
  nlohmann::json dataFlags = nlohmann::json();         // DataFlagView::toJson() (caras, grupos, IMon)
};

nlohmann::json dashJson(const DashInputs& in);

}  // namespace daq::daqrun
