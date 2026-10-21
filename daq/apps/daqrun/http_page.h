// http_page.h
// Página de solo lectura de daqrun (THttpServer en 127.0.0.1:8080 por omisión; remoto por túnel
// SSH). Fase 3: tasa de trigger, Δt, tasa de pulsos fuera de tiempo por canal, espectro de
// amplitud (mV, umbrales PROVISIONALES: estado sin calibración), niveles de cola, resumen de flags
// y estado de hvmon (aviso si no corre). /dash/: página numérica (dash.html embebido) que lee
// /status/state_json/root.json.gz (objeto "dash") cada 2 s. Todo objeto ROOT pertenece al hilo principal; update() y
// processRequests() se llaman solo desde él.
#pragma once

#include <memory>
#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/monitor.h"
#include "daq/io/panels_root.h"

class THttpServer;
class TH1D;
class TH2F;
class TGraph;
class TNamed;

namespace daq::daqrun {

class DaqrunPage {
 public:
  static Result<std::unique_ptr<DaqrunPage>> create(const DaqrunConfig& cfg, const std::string& bindAddress,
                                                     const std::string& mappingText = "{}");
  ~DaqrunPage();
  void update(const MonitorView& v, const HvSample& hv, const nlohmann::json& status);
  // Paneles en tiempo y bits de datos (fase 5); solo cuando cambió la vista del monitor.
  void updatePanels(const PanelView& p, const DataFlagView& d, double leakRatio);
  void processRequests();
  const std::string& url() const { return url_; }

 private:
  DaqrunPage() = default;
  std::unique_ptr<THttpServer> server_;
  std::unique_ptr<TH1D> dt_, oot_, chFlags_, evtFlags_;
  std::unique_ptr<TH2F> amp_;
  std::unique_ptr<TGraph> rate_;
  std::unique_ptr<TNamed> table_, json_, mapping_;
  io::PanelRoot panels_;
  std::string url_;
};

}  // namespace daq::daqrun
