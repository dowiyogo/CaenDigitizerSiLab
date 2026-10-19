// mapping.cpp
#include "daq/hvmon/mapping.h"

#include <cstdlib>
#include <fstream>

#include <nlohmann/json.hpp>

namespace daq::hvmon {

Result<Mapping> Mapping::fromFile(const std::string& path) {
  using R = Result<Mapping>;
  std::ifstream in(path);
  if (!in) return R::failure(-1, "no se pudo abrir el mapeo " + path);
  try {
    nlohmann::json j;
    in >> j;
    Mapping m;
    m.version = j.at("mapping_version").get<std::string>();
    for (auto it = j.at("hv_group_to_sy5527_channel").begin(); it != j["hv_group_to_sy5527_channel"].end(); ++it)
      m.groupToHvCh[std::atoi(it.key().c_str())] = it.value().get<int>();
    for (const auto& e : j.at("channels")) {
      MappingEntry x;
      x.digCh = e.at("dig_ch").get<int>();
      x.station = e.value("station", "");
      x.plane = e.value("plane", "");
      x.strip = e.value("strip", -1);
      if (e.contains("hv_group") && e["hv_group"].is_number_integer()) x.hvGroup = e["hv_group"].get<int>();
      if (e.contains("face") && e["face"].is_string()) x.face = e["face"].get<std::string>();
      if (e.contains("confirmed")) x.hvGroupConfirmed = e["confirmed"].value("hv_group", "unknown");
      m.channels.push_back(x);
    }
    return m;
  } catch (const std::exception& ex) {
    return R::failure(-1, "mapeo inválido " + path + ": " + ex.what());
  }
}

std::optional<int> Mapping::groupOfHvChannel(int hvCh) const {
  for (const auto& kv : groupToHvCh)
    if (kv.second == hvCh) return kv.first;
  return std::nullopt;
}

std::vector<MappingEntry> Mapping::forGroup(int group) const {
  std::vector<MappingEntry> out;
  for (const auto& e : channels)
    if (e.hvGroup && *e.hvGroup == group) out.push_back(e);
  return out;
}

int Mapping::unknownGroupCount(const std::string& station) const {
  int n = 0;
  for (const auto& e : channels)
    if (e.station == station && !e.hvGroup) ++n;
  return n;
}

}  // namespace daq::hvmon
