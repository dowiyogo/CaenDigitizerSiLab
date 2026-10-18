// status_bits.h
// Decodificación de la palabra Status de un canal HV en banderas con nombre. Las posiciones
// vienen de la configuración; los valores por defecto son los del manual del HV Wrapper
// (Rev. 13, 2016, PARAM_TYPE_CHSTATUS), según docs/daq_v2_design.md, sección 6.1. El
// significado para el A7030DP se confirma en la fase 2.
// Tipos de valor; seguros entre hilos para lectura concurrente.
#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <nlohmann/json_fwd.hpp>

#include "daq/core/result.h"

namespace daq::hv {

class StatusBitMap {
 public:
  // on 0, ramp_up 1, ramp_down 2, overcurrent 3, overvoltage 4, undervoltage 5,
  // external_trip 6, max_v 7, external_disable 8, internal_trip 9, calibration_error 10,
  // unplugged 11.
  static StatusBitMap wrapperRev13();
  // Objeto JSON {"nombre": bit, ...}; valida bits en [0, 31] distintos y nombres únicos.
  static Result<StatusBitMap> fromJson(const nlohmann::json& j);

  // Banderas activas en el orden de bit; los bits activos sin nombre aparecen como "bitN".
  std::vector<std::string> decode(std::uint32_t word) const;
  bool has(std::uint32_t word, std::string_view name) const;  // false si el nombre no existe
  int bitOf(std::string_view name) const;                       // -1 si no existe
  const std::vector<std::pair<std::string, int>>& bits() const { return bits_; }

 private:
  std::vector<std::pair<std::string, int>> bits_;  // ordenado por bit
};

}  // namespace daq::hv
