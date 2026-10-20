// state.cpp
#include "daq/hvmon/state.h"

#include <unistd.h>

#include "daq/hv/reading_json.h"

namespace daq::hvmon {

nlohmann::json ChangeRecord::toJson() const {
  nlohmann::json j{{"type", "change"},     {"id", id},           {"t_utc_ns", tUtcNs},
                   {"t_utc", formatUtcIso(tUtcNs)}, {"t_prev_read_utc_ns", tPrevReadUtcNs},
                   {"slot", slot},         {"ch", ch},           {"param", param},
                   {"unit", unit},         {"old", oldValue},    {"new", newValue},
                   {"origin", origin},     {"during_disconnect", duringDisconnect},
                   {"read_kind", readKind}};
  j["command_id"] = commandId ? nlohmann::json(*commandId) : nlohmann::json();
  if (!detail.is_null()) j["detail"] = detail;
  return j;
}

nlohmann::json CommandRecord::toJson() const {
  return nlohmann::json{{"type", "command"}, {"command_id", id}, {"t_utc_ns", tUtcNs}, {"t_utc", formatUtcIso(tUtcNs)},
                        {"slot", slot},      {"ch", ch},         {"param", param},    {"value", value},
                        {"phase", phase},    {"reason", reason}, {"requester", requester}, {"detail", detail}};
}

nlohmann::json ReadStats::toJson() const {
  return nlohmann::json{{"fast_reads", fastReads},
                        {"full_reads", fullReads},
                        {"read_errors", readErrors},
                        {"connect_attempts", connectAttempts},
                        {"connect_failures", connectFailures},
                        {"reconnections", reconnections},
                        {"fast_read_mean_ms", fastReads ? sumFastMs / static_cast<double>(fastReads) : 0.0},
                        {"fast_read_max_ms", maxFastMs},
                        {"full_read_mean_ms", fullReads ? sumFullMs / static_cast<double>(fullReads) : 0.0},
                        {"full_read_max_ms", maxFullMs}};
}

nlohmann::json WriterHealth::toJson() const {
  return nlohmann::json{{"ok", ok},         {"last_error", lastError}, {"last_error_utc_ns", lastErrorUtcNs},
                        {"errors", errors}, {"fsyncs", fsyncs},        {"bytes", bytes},
                        {"lines", lines},   {"open_files", openFiles}};
}

namespace {
nlohmann::json hub_notifier(const StateHub& h) { return h.notifier.read().first; }
}  // namespace

nlohmann::json StateHub::compose(const HvmonConfig& cfg, UtcNs now, std::uint64_t seq, bool final) const {
  const auto [p, pseq] = poller.read();
  const auto [w, wseq] = writer.read();
  const auto [temp, tseq] = temperature.read();
  (void)wseq;
  (void)tseq;
  const double age = p.lastGoodReadUtcNs > 0 ? static_cast<double>(now - p.lastGoodReadUtcNs) / 1e9 : -1.0;
  const bool stale = !p.connected || p.lastGoodReadUtcNs == 0 || age > cfg.staleAfterS;

  nlohmann::json channels = nlohmann::json::array();
  for (int ch : cfg.channels) {
    nlohmann::json values = nlohmann::json::object();
    nlohmann::json flags = nlohmann::json::array();
    auto it = p.values.find(ch);
    if (it != p.values.end()) {
      for (const auto& info : p.params) {
        auto v = it->second.find(info.name);
        if (v == it->second.end()) continue;
        if (info.type == hv::ParamType::Numeric)
          values[hv::columnName(info)] = v->second.si;
        else
          values[info.name] = v->second.raw;
        if (info.type == hv::ParamType::ChStatus) flags = cfg.statusBits.decode(v->second.raw);
      }
    }
    nlohmann::json retry{{"state", "UNINIT"}};
    for (const auto& r : p.retry)
      if (r.value("ch", -1) == ch) retry = r;
    nlohmann::json win;  // null sin lecturas rápidas todavía
    if (auto wi = p.window.find(ch); wi != p.window.end()) {
      auto js = [](const WindowSummary& x) {
        return nlohmann::json{{"n", x.n}, {"mean", x.mean}, {"std", x.std}, {"min", x.min}, {"max", x.max}};
      };
      win = {{"window_s", p.windowS}, {"VMon_V", js(wi->second.first)}, {"IMon_A", js(wi->second.second)}};
    }
    channels.push_back({{"ch", ch}, {"values", values}, {"status_flags", flags}, {"auto_retry", retry}, {"stats_10min", win}});
  }
  nlohmann::json changes = nlohmann::json::array();
  for (const auto& c : p.recentChanges) changes.push_back(c.toJson());
  nlohmann::json board;  // null si la placa no informa parámetros
  nlohmann::json boardTemp;
  if (p.board) {
    nlohmann::json values = nlohmann::json::object();
    nlohmann::json flags = nlohmann::json::array();
    for (std::size_t i = 0; i < p.board->params.size(); ++i) {
      const auto& info = p.board->params[i];
      const auto& v = p.board->values[i];
      if (info.type == hv::ParamType::Numeric)
        values[hv::columnName(info)] = v.si;
      else
        values[info.name] = v.raw;
      if (info.type == hv::ParamType::BdStatus) flags = hv::decodeBoardStatus(v.raw);
      if (info.name == "Temp") boardTemp = v.si;
    }
    board = {{"slot", p.board->slot}, {"t_utc_ns", p.board->tRespUtcNs}, {"values", values}, {"status_flags", flags}};
  }

  return nlohmann::json{
      {"schema", "daqv2.hvmon_state/1"},
      {"seq", seq},
      {"t_utc_ns", now},
      {"t_utc", formatUtcIso(now)},
      {"final", final},
      {"pid", static_cast<long>(::getpid())},
      {"started_utc_ns", startedUtcNs},
      {"system", cfg.system},
      {"slot", cfg.slot},
      {"allow_commands", cfg.allowCommands},
      {"connection",
       {{"connected", p.connected},
        {"ever_connected", p.everConnected},
        {"stale", stale},
        {"last_good_read_utc_ns", p.lastGoodReadUtcNs},
        {"age_s", age},
        {"last_error", p.lastError},
        {"last_error_utc_ns", p.lastErrorUtcNs},
        {"next_connect_in_s", p.nextConnectInS}}},
      {"read_seq", p.readSeq},
      {"last_read",
       {{"read_seq", p.lastReadSeq},
        {"kind", p.lastReadKind},
        {"t_req_utc_ns", p.lastReadReqUtcNs},
        {"t_resp_utc_ns", p.lastReadRespUtcNs}}},
      {"poller_seq", pseq},
      {"stats", p.stats.toJson()},
      {"channels", channels},
      {"changes", changes},
      {"geco_running", gecoRunning.load()},
      {"writer", [&] {
         auto j = w.toJson();
         j["dropped"] = writerDropped.load();
         return j;
       }()},
      {"active_command", p.activeCommand},
      {"active_offsets", p.activeOffsets},
      {"last_offsets", p.lastOffsets},
      {"temperature", temp},
      {"notifier", hub_notifier(*this)},
      {"compressor", compressor.read().first},
      {"board", board},
      {"board_temperature_degC", boardTemp},
      {"board_temperature_note", "temperatura de la placa A7030DP (sigue la del rack), no la del detector"}};
}

}  // namespace daq::hvmon
