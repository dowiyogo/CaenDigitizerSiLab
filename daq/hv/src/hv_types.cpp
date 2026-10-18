// hv_types.cpp
#include "daq/hv/hv_types.h"

#include <algorithm>

#include "daq/core/units.h"

namespace daq::hv {

const char* toString(ParamType t) {
  switch (t) {
    case ParamType::Numeric: return "numeric";
    case ParamType::OnOff: return "onoff";
    case ParamType::ChStatus: return "chstatus";
    case ParamType::BdStatus: return "bdstatus";
    case ParamType::Binary: return "binary";
    case ParamType::String: return "string";
    case ParamType::Enum: return "enum";
    case ParamType::Cmd: return "cmd";
    case ParamType::Unknown: return "unknown";
  }
  return "unknown";
}

const char* toString(ParamMode m) {
  switch (m) {
    case ParamMode::ReadOnly: return "rdonly";
    case ParamMode::WriteOnly: return "wronly";
    case ParamMode::ReadWrite: return "rdwr";
    case ParamMode::Unknown: return "unknown";
  }
  return "unknown";
}

ParamType paramTypeFromInt(int v) {
  if (v >= 0 && v <= 7) return static_cast<ParamType>(v);
  return ParamType::Unknown;
}

ParamMode paramModeFromInt(int v) {
  if (v >= 0 && v <= 2) return static_cast<ParamMode>(v);
  return ParamMode::Unknown;
}

std::optional<ParamType> paramTypeFromString(const std::string& s) {
  for (int v = 0; v <= 7; ++v) {
    if (s == toString(static_cast<ParamType>(v))) return static_cast<ParamType>(v);
  }
  return std::nullopt;
}

std::string ParamInfo::siUnit() const {
  if (type != ParamType::Numeric) return "";
  return units::unitSymbol(unitCode);
}

std::string ParamInfo::boardUnit() const {
  if (type != ParamType::Numeric) return "";
  return units::decimalPrefix(exponent) + units::unitSymbol(unitCode);
}

const ParamValue* ChannelReading::get(const std::string& param, int channel) const {
  auto ip = std::find_if(params.begin(), params.end(), [&](const ParamInfo& p) { return p.name == param; });
  auto ic = std::find(channels.begin(), channels.end(), channel);
  if (ip == params.end() || ic == channels.end()) return nullptr;
  const auto pi = static_cast<std::size_t>(ip - params.begin());
  const auto ci = static_cast<std::size_t>(ic - channels.begin());
  return &values[pi][ci];
}

const ParamValue* BoardReading::get(const std::string& param) const {
  for (std::size_t i = 0; i < params.size(); ++i)
    if (params[i].name == param) return &values[i];
  return nullptr;
}

std::vector<std::string> decodeBoardStatus(std::uint32_t word) {
  static const char* names[] = {"power_fail", "firmware_checksum_error", "hv_calibration_error",
                                "temperature_calibration_error", "under_temperature", "over_temperature"};
  std::vector<std::string> out;
  for (int b = 0; b < 32; ++b) {
    if (((word >> b) & 1u) == 0u) continue;
    out.push_back(b < 6 ? std::string(names[b]) : "bit" + std::to_string(b));
  }
  return out;
}

ParamValue numericValue(const ParamInfo& info, double boardValue) {
  ParamValue v;
  v.type = ParamType::Numeric;
  v.boardValue = boardValue;
  v.si = units::applyDecimalExponent(boardValue, info.exponent);
  return v;
}

ParamValue integerValue(const ParamInfo& info, std::uint32_t raw) {
  ParamValue v;
  v.type = info.type;
  v.raw = raw;
  return v;
}

}  // namespace daq::hv
