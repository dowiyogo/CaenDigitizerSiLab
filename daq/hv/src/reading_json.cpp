// reading_json.cpp
#include "daq/hv/reading_json.h"

namespace daq::hv {

std::string columnName(const ParamInfo& p) {
  if (p.type != ParamType::Numeric) return p.name;
  std::string u = p.siUnit();
  if (u.empty()) return p.name;
  std::string tag;
  for (char c : u) tag += (c == '/') ? std::string("_per_") : std::string(1, c);
  return p.name + "_" + tag;
}

nlohmann::json readingToJson(const ChannelReading& r, const StatusBitMap& bits) {
  nlohmann::json params = nlohmann::json::array();
  for (std::size_t p = 0; p < r.params.size(); ++p) {
    const auto& info = r.params[p];
    nlohmann::json e{{"name", info.name}, {"type", toString(info.type)}, {"type_code", static_cast<int>(info.type)}};
    if (info.type == ParamType::Numeric) {
      nlohmann::json si = nlohmann::json::array(), board = nlohmann::json::array();
      for (const auto& v : r.values[p]) {
        si.push_back(v.si);
        board.push_back(v.boardValue);
      }
      e["unit_si"] = info.siUnit();
      e["unit_board"] = info.boardUnit();
      e["unit_code"] = info.unitCode;
      e["exp"] = info.exponent;
      e["unit_from_board"] = info.unitFromBoard;
      e["values_si"] = si;
      e["values_board"] = board;
    } else {
      nlohmann::json vals = nlohmann::json::array();
      for (const auto& v : r.values[p]) vals.push_back(v.raw);
      e["values"] = vals;
      if (info.type == ParamType::ChStatus) {
        nlohmann::json dec = nlohmann::json::array();
        for (const auto& v : r.values[p]) dec.push_back(bits.decode(v.raw));
        e["decoded"] = dec;
      }
    }
    params.push_back(e);
  }
  return nlohmann::json{{"slot", r.slot},
                        {"channels", r.channels},
                        {"t_req_utc_ns", r.tReqUtcNs},
                        {"t_resp_utc_ns", r.tRespUtcNs},
                        {"t_req_utc", formatUtcIso(r.tReqUtcNs)},
                        {"duration_ns", r.durationNs},
                        {"params", params}};
}

Result<ChannelReading> readingFromJson(const nlohmann::json& j) {
  using R = Result<ChannelReading>;
  try {
    ChannelReading r;
    r.slot = j.at("slot").get<int>();
    r.channels = j.at("channels").get<std::vector<int>>();
    r.tReqUtcNs = j.at("t_req_utc_ns").get<std::int64_t>();
    r.tRespUtcNs = j.at("t_resp_utc_ns").get<std::int64_t>();
    r.durationNs = j.at("duration_ns").get<std::int64_t>();
    for (const auto& e : j.at("params")) {
      ParamInfo info;
      info.name = e.at("name").get<std::string>();
      info.type = paramTypeFromInt(e.at("type_code").get<int>());
      std::vector<ParamValue> row;
      if (info.type == ParamType::Numeric) {
        info.unitCode = e.value("unit_code", -1);
        info.exponent = e.value("exp", 0);
        info.unitFromBoard = e.value("unit_from_board", false);
        const auto& si = e.at("values_si");
        const auto& board = e.at("values_board");
        if (si.size() != r.channels.size() || board.size() != r.channels.size())
          return R::failure(-1, "lectura JSON: tamaño de valores incorrecto en " + info.name);
        for (std::size_t c = 0; c < r.channels.size(); ++c) {
          ParamValue v;
          v.type = ParamType::Numeric;
          v.si = si[c].get<double>();
          v.boardValue = board[c].get<double>();
          row.push_back(v);
        }
      } else {
        const auto& vals = e.at("values");
        if (vals.size() != r.channels.size()) return R::failure(-1, "lectura JSON: tamaño incorrecto en " + info.name);
        for (const auto& x : vals) row.push_back(integerValue(info, x.get<std::uint32_t>()));
      }
      r.params.push_back(info);
      r.values.push_back(std::move(row));
    }
    return r;
  } catch (const std::exception& ex) {
    return R::failure(-1, std::string("lectura JSON inválida: ") + ex.what());
  }
}

}  // namespace daq::hv
