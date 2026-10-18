// http_page.h
// Página web de solo lectura de hvmon (THttpServer, puerto 8081, enlazado a 127.0.0.1 o a la
// dirección de tailscale0). Objetos servidos:
//   /CHxx/VMon_CHxx, /CHxx/IMon_CHxx   TGraph de las últimas 24 h (anillo de tamaño fijo, un punto
//                                      cada graph_period_s), eje x en tiempo UTC, unidades SI
//   /board/BoardTemp_degC              TGraph de la temperatura de la placa A7030DP (no del
//                                      detector; sigue la del rack), 24 h
//   /status/table                      TNamed; el título es la tabla de estado (texto)
//   /status/state_json                 TNamed; el título es el JSON completo del estado
// SetReadOnly(kTRUE) y ningún comando registrado: no hay acciones por HTTP. SetTimer(0): las
// peticiones solo se atienden cuando el hilo principal llama processRequests().
// Seguridad entre hilos: TODO objeto ROOT de esta clase pertenece al hilo principal.
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/state.h"

class THttpServer;
class TGraph;
class TNamed;

namespace daq::hvmon {

class HttpPage {
 public:
  static Result<std::unique_ptr<HttpPage>> create(const HvmonConfig& cfg, const std::string& bindAddress);
  ~HttpPage();

  // Agrega un punto por canal si pasó graph_period_s y actualiza la tabla y el JSON.
  void update(const PollerState& p, const nlohmann::json& state, UtcNs now);
  void processRequests();
  std::string url() const { return url_; }

 private:
  HttpPage() = default;
  struct Ring {
    std::vector<double> t, y;
    std::size_t head = 0, n = 0;
    void push(double tt, double yy);
  };
  struct Channel {
    TGraph* vmon = nullptr;  // dueños: graphs_
    TGraph* imon = nullptr;
    Ring vRing, iRing;
  };
  void refill(TGraph* g, const Ring& r);

  const HvmonConfig* cfg_ = nullptr;
  std::unique_ptr<THttpServer> server_;
  std::vector<std::unique_ptr<TGraph>> graphs_;
  std::unique_ptr<TNamed> table_;
  std::unique_ptr<TNamed> json_;
  std::map<int, Channel> channels_;
  TGraph* boardTemp_ = nullptr;
  Ring boardRing_;
  UtcNs lastBoardUtc_ = 0;
  UtcNs lastGraphUtc_ = 0;
  std::string url_;
};

}  // namespace daq::hvmon
