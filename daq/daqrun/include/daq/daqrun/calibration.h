// calibration.h
// Calibración de p.e. de una campaña (diseño 7.4; DAQ2-10). Archivo JSON versionado
// (schema daqv2.pe_calibration/1) con el run del que sale, el estado de HV y la temperatura de la
// placa, y por canal: G (mV/p.e.), pedestal (mV), estado del ajuste y si es provisional.
//  - hod2018: G por separación de picos (gain_fit.h); estados ok | changed | ambiguous | rejected |
//    insufficient; los que no son ok quedan provisionales (bit 10 de ch_flags) con el valor que
//    tengan (o sin G).
//  - hod2019: sin picos resueltos: tasa oscura > 3 mV como indicador; sin G, provisional.
// daqrun lo carga si calibration.file está en la configuración: el run sale del estado "sin
// calibración" (bit 7 de evt_flags en 0), la cabecera guarda el JSON completo y su versión.
#pragma once

#include <array>
#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"

namespace daq::daqrun {

struct PeCalibration {
  std::string version, path;
  nlohmann::json raw;
  std::array<double, 64> gainMv{}, pedestalMv{};
  std::array<bool, 64> provisional{};
  static Result<PeCalibration> fromJson(const nlohmann::json& j);
  static Result<PeCalibration> fromFile(const std::string& path);
};

}  // namespace daq::daqrun
