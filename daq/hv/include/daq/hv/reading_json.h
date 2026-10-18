// reading_json.h
// Serialización de una lectura HV a JSON y vuelta, compartida por hvctl (modo directo y cliente)
// y por hvmon (respuesta "read" del socket). Los parámetros van en un arreglo ordenado (el orden
// de la placa), cada uno con tipo, unidades y valores por canal:
//   {"slot":3, "channels":[0..], "t_req_utc_ns":..., "t_resp_utc_ns":..., "duration_ns":...,
//    "params":[{"name":"VMon","type":"numeric","type_code":0,"unit_si":"V","unit_board":"V",
//               "exp":0,"values_si":[...],"values_board":[...]},
//              {"name":"Status","type":"chstatus","type_code":2,"values":[...],"decoded":[[...]]}]}
// Funciones puras, seguras entre hilos.
#pragma once

#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/hv/hv_types.h"
#include "daq/hv/status_bits.h"

namespace daq::hv {

nlohmann::json readingToJson(const ChannelReading& r, const StatusBitMap& bits);
Result<ChannelReading> readingFromJson(const nlohmann::json& j);

// Nombre de columna con unidad SI: "VMon_V", "IMon_A", "RUp_V_per_s", "Trip_s", "Status".
std::string columnName(const ParamInfo& p);

}  // namespace daq::hv
