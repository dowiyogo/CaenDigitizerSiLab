// hv_types.h
// Tipos de datos del control lento de HV, independientes del backend (Wrapper real o fake).
// Unidades: los valores numéricos de ParamValue::si están en SI (V, A, V/s, s); boardValue
// guarda el número tal como lo entrega la placa (p. ej. µA) para mostrarlo como hv_readall.
// Tiempos en UTC ns; duraciones en ns de steady_clock. Tipos de valor, sin estado compartido.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "daq/core/time.h"

namespace daq::hv {

// Valores iguales a PARAM_TYPE_* de CAENHVWrapper.h.
enum class ParamType : int {
  Numeric = 0,
  OnOff = 1,
  ChStatus = 2,
  BdStatus = 3,
  Binary = 4,
  String = 5,
  Enum = 6,
  Cmd = 7,
  Unknown = -1
};
// Valores iguales a PARAM_MODE_* de CAENHVWrapper.h.
enum class ParamMode : int { ReadOnly = 0, WriteOnly = 1, ReadWrite = 2, Unknown = -1 };

const char* toString(ParamType t);
const char* toString(ParamMode m);
ParamType paramTypeFromInt(int v);
ParamMode paramModeFromInt(int v);
std::optional<ParamType> paramTypeFromString(const std::string& s);

struct ParamInfo {
  std::string name;
  ParamType type = ParamType::Unknown;
  ParamMode mode = ParamMode::Unknown;
  // Solo numéricos: unidad PARAM_UN_* y exponente decimal de la placa.
  int unitCode = -1;
  int exponent = 0;
  bool unitFromBoard = false;    // false: unidad tomada de la tabla de respaldo
  std::optional<double> minBoard;  // en unidades de la placa
  std::optional<double> maxBoard;
  // Solo ON/OFF: textos de estado.
  std::string onState;
  std::string offState;

  std::string siUnit() const;     // "V", "A", "V/s", "s" o "" si no es numérico
  std::string boardUnit() const;  // "V", "uA", ...
};

struct SlotInfo {
  int slot = 0;
  std::string model;        // vacío si el slot está vacío
  std::string description;
  int channels = 0;
  unsigned serial = 0;
  int fwMajor = 0;
  int fwMinor = 0;
  bool empty() const { return model.empty(); }
};

struct CrateMap {
  std::vector<SlotInfo> slots;  // un elemento por slot del crate, incluidos los vacíos
};

struct ParamValue {
  ParamType type = ParamType::Unknown;
  double si = 0.0;          // numéricos: valor en SI
  double boardValue = 0.0;  // numéricos: valor en unidades de la placa
  std::uint32_t raw = 0;    // no numéricos: entero tal como lo entrega la Wrapper
  bool isNumeric() const { return type == ParamType::Numeric; }
};

struct ChannelReading {
  int slot = 0;
  std::vector<int> channels;
  std::vector<ParamInfo> params;
  std::vector<std::vector<ParamValue>> values;  // values[iParam][iChannel]
  UtcNs tReqUtcNs = 0;    // UTC justo antes de la primera llamada
  UtcNs tRespUtcNs = 0;   // UTC justo después de la última llamada
  std::int64_t durationNs = 0;  // duración medida con steady_clock

  // nullptr si el parámetro o el canal no están en la lectura.
  const ParamValue* get(const std::string& param, int channel) const;
};

// Construye el valor de un parámetro numérico desde el número de la placa.
ParamValue numericValue(const ParamInfo& info, double boardValue);
// Construye el valor de un parámetro entero (ON/OFF, Status, binario).
ParamValue integerValue(const ParamInfo& info, std::uint32_t raw);

}  // namespace daq::hv
