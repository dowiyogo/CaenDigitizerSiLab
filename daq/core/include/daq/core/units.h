// units.h
// Conversión de unidades en la frontera con el hardware. Por encima de los backends todo está
// en SI: V, A, V/s, s. El A7030DP informa corrientes en µA; la Wrapper describe cada parámetro
// numérico con una unidad (PARAM_UN_*) y un exponente decimal (Exp: -6 = micro).
// Funciones puras, seguras entre hilos.
#pragma once

#include <cmath>
#include <string>

namespace daq::units {

constexpr double kMicro = 1e-6;

constexpr double microampereToAmpere(double uA) noexcept { return uA / 1e6; }
constexpr double ampereToMicroampere(double A) noexcept { return A * 1e6; }
constexpr double volt(double v) noexcept { return v; }  // V se mantiene V

// valor_SI = valor_placa * 10^exp. Con exp < 0 se divide por 10^-exp (exacto en potencias de 10),
// para que 15 µA dé 1.5e-05 A y no 1.4999999999999999e-05.
inline double applyDecimalExponent(double value, int exp) {
  return exp < 0 ? value / std::pow(10.0, -exp) : value * std::pow(10.0, exp);
}
inline double removeDecimalExponent(double valueSI, int exp) { return applyDecimalExponent(valueSI, -exp); }

// Símbolo SI de un código PARAM_UN_* de CAENHVWrapper.h (0 NONE ... 12 APS).
inline std::string unitSymbol(int unitCode) {
  switch (unitCode) {
    case 0: return "";
    case 1: return "A";
    case 2: return "V";
    case 3: return "W";
    case 4: return "degC";
    case 5: return "Hz";
    case 6: return "bar";
    case 7: return "V/s";
    case 8: return "s";
    case 9: return "rpm";
    case 10: return "count";
    case 11: return "bit";
    case 12: return "A/s";
    default: return "?";
  }
}

// Prefijo para el exponente decimal (-6 -> "u", -3 -> "m", 0 -> "", 3 -> "k", 6 -> "M").
inline std::string decimalPrefix(int exp) {
  switch (exp) {
    case -9: return "n";
    case -6: return "u";
    case -3: return "m";
    case 0: return "";
    case 3: return "k";
    case 6: return "M";
    default: return "e" + std::to_string(exp);
  }
}

}  // namespace daq::units
