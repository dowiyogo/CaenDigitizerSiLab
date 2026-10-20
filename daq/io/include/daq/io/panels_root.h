// panels_root.h
// Objetos ROOT de los paneles en tiempo de daqrun (fase 5): los mismos para la página
// (THttpServer, hilo principal) y para el modo offline (archivo ROOT). Carpetas:
//   overview/  tasa de trigger (la registra la página), mapa canal × tiempo de la tasa fuera de
//              tiempo, tasa por cara con su banda de base, canales vivos por plano
//   intime/    amplitud, integral y muestra del mínimo por canal (acumulado y _1h)
//   hits/      mapas X × Y por estación, coincidencias por plano, multiplicidad (acumulado y _1h)
//   scope/     persistencia por grupo del digitizador y último evento
//   gain/      G por hora de hod2018 (valor y estado) y tasa > 3 mV de hod2019, PROVISIONAL
//   dataflags/ tasa por cara y por grupo, base, IMon por canal HV
// Banda de base de una cara: de la base a base × light_leak_ratio (umbral del bit 7).
// Pertenece a un solo hilo (el dueño de los objetos ROOT).
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "daq/daqrun/data_flags.h"
#include "daq/daqrun/panels.h"

class TH1;
class TH2F;
class TGraph;
class TGraphAsymmErrors;
class THttpServer;
class TDirectory;
class TObject;

namespace daq::io {

class PanelRoot {
 public:
  PanelRoot();
  ~PanelRoot();
  // Crea o actualiza los objetos. leakRatio: umbral de razón del bit 7 (banda de base).
  void update(const daqrun::PanelView& p, const daqrun::DataFlagView& d, double leakRatio);
  void registerIn(THttpServer& s);  // en sus carpetas; se puede llamar tras cada update (solo agrega los nuevos)
  void unregisterFrom(THttpServer& s);
  void write(TDirectory& dir) const;  // offline: subdirectorios por carpeta
  // Ítems por omisión de la vista general (rutas de la página).
  static std::vector<std::string> defaultDrawItems();

 private:
  TH2F* h2(const std::string& folder, const std::string& name, const std::string& title, int nx, double x0, double x1, int ny,
           double y0, double y1);
  TH1* h1(const std::string& folder, const std::string& name, const std::string& title, int nx, double x0, double x1);
  TGraph* graph(const std::string& folder, const std::string& name, const std::string& title);
  std::map<std::string, std::pair<std::string, std::unique_ptr<TObject>>> obj_;  // nombre -> (carpeta, objeto)
  std::vector<TObject*> registered_;  // ya registrados (registerIn agrega solo los nuevos)
};

}  // namespace daq::io
