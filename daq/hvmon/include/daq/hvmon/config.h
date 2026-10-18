// config.h
// Configuración de hvmon (JSON). Valores por defecto = configuración real de la fase 1 (solo
// lectura, allow_commands = false). Unidades: s, V, A. Ningún secreto vive aquí: usuario y
// contraseña del SY5527 vienen de secrets.env (daq/core/secrets.h).
// Tipo de valor; seguro entre hilos para lectura concurrente una vez construido.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/hv/status_bits.h"

namespace daq::hvmon {

struct HvmonConfig {
  // HV
  std::string system = "SY5527";
  std::string host = "192.168.1.118";  // DAQ_HV_HOST tiene prioridad
  int slot = 3;
  std::vector<int> channels{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  std::vector<std::string> fastParams{"V0Set", "I0Set", "VMon", "IMon", "Status", "Pw"};
  std::vector<std::string> fastTracked{"V0Set", "I0Set", "Pw", "Status"};  // cambios en lecturas rápidas
  double fastPeriodS = 1.0;
  double fullPeriodS = 60.0;
  double backoffInitialS = 1.0;
  double backoffMaxS = 60.0;
  double staleAfterS = 3.0;
  hv::StatusBitMap statusBits = hv::StatusBitMap::wrapperRev13();
  std::vector<std::string> tripFlags{"internal_trip", "external_trip"};
  std::string legacyCredentials = "~/.config/daqv2/hv_credentials.json";

  // Escritura (única operación: Pw de un canal)
  bool allowCommands = false;
  double verifyMarginS = 15.0;    // timeout = V0Set/RUp (o V0Set/RDWn) + margen
  // La verificación usa Status (on, sin bits de rampa) y VMon como control grueso: con HV apagada
  // el A7030DP lee VMon de 0 a 1.47 V (offset), así que las tolerancias no pueden ser menores.
  double vmonOnToleranceV = 2.0;  // |VMon - V0Set| máximo al terminar una rampa
  double vmonOffMaxV = 2.0;       // VMon máximo para dar por apagado
  double commandWindowExtraS = 5.0;  // ventana de origen hvmon_command tras la verificación

  // Archivos
  std::string slowcontrolDir = "/data/daq/slowcontrol";
  double fsyncPeriodS = 10.0;
  std::string runtimeDir;  // por defecto $XDG_RUNTIME_DIR/daqv2
  std::string stateDir = "/dev/shm/daqv2";
  std::string logFile;     // vacío: solo stderr (journal de systemd)
  int changeRingSize = 50;
  std::size_t writerQueueCapacity = 100000;

  // Otros
  std::string gecoPattern = "bin/CAENGECO2020";
  double gecoPeriodS = 10.0;
  std::string httpBind = "127.0.0.1";  // "127.0.0.1", "localhost" o "tailscale"
  int httpPort = 8081;
  bool httpEnabled = true;
  std::vector<int> httpChannels{0, 1, 2, 3, 4, 5, 6, 7};
  double graphPeriodS = 10.0;
  double graphHours = 24.0;
  std::string tempBackend = "none";
  double tempPeriodS = 30.0;
  nlohmann::json tempConfig = nlohmann::json::object();

  std::string lockPath() const { return runtimeDir + "/hvmon.lock"; }
  std::string socketPath() const { return runtimeDir + "/hvmon.sock"; }
  std::string statePath() const { return stateDir + "/hvmon_state.json"; }

  static Result<HvmonConfig> fromJson(const nlohmann::json& j);
  static Result<HvmonConfig> fromFile(const std::string& path);
  nlohmann::json toJson() const;
};

std::string defaultRuntimeDir();  // $XDG_RUNTIME_DIR/daqv2 (o /run/user/<uid>/daqv2)

}  // namespace daq::hvmon
