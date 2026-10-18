// status_bits.cpp
#include "daq/hv/status_bits.h"

#include <algorithm>
#include <set>

#include <nlohmann/json.hpp>

namespace daq::hv {

StatusBitMap StatusBitMap::wrapperRev13() {
  StatusBitMap m;
  m.bits_ = {{"on", 0},           {"ramp_up", 1},          {"ramp_down", 2},
             {"overcurrent", 3},  {"overvoltage", 4},      {"undervoltage", 5},
             {"external_trip", 6}, {"max_v", 7},           {"external_disable", 8},
             {"internal_trip", 9}, {"calibration_error", 10}, {"unplugged", 11}};
  return m;
}

Result<StatusBitMap> StatusBitMap::fromJson(const nlohmann::json& j) {
  if (!j.is_object()) return Result<StatusBitMap>::failure(-1, "status_bits: se esperaba un objeto");
  StatusBitMap m;
  std::set<int> used;
  for (auto it = j.begin(); it != j.end(); ++it) {
    if (!it.value().is_number_integer())
      return Result<StatusBitMap>::failure(-1, "status_bits." + it.key() + ": se esperaba un entero");
    const int bit = it.value().get<int>();
    if (bit < 0 || bit > 31)
      return Result<StatusBitMap>::failure(-1, "status_bits." + it.key() + ": bit fuera de [0, 31]");
    if (!used.insert(bit).second)
      return Result<StatusBitMap>::failure(-1, "status_bits: bit " + std::to_string(bit) + " repetido");
    m.bits_.emplace_back(it.key(), bit);
  }
  std::sort(m.bits_.begin(), m.bits_.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
  return m;
}

std::vector<std::string> StatusBitMap::decode(std::uint32_t word) const {
  std::vector<std::string> out;
  for (int b = 0; b < 32; ++b) {
    if (((word >> b) & 1u) == 0u) continue;
    auto it = std::find_if(bits_.begin(), bits_.end(), [b](const auto& p) { return p.second == b; });
    out.push_back(it != bits_.end() ? it->first : "bit" + std::to_string(b));
  }
  return out;
}

int StatusBitMap::bitOf(std::string_view name) const {
  auto it = std::find_if(bits_.begin(), bits_.end(), [&](const auto& p) { return p.first == name; });
  return it == bits_.end() ? -1 : it->second;
}

bool StatusBitMap::has(std::uint32_t word, std::string_view name) const {
  const int b = bitOf(name);
  return b >= 0 && ((word >> b) & 1u) != 0u;
}

}  // namespace daq::hv
