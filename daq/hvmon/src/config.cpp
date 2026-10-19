// config.cpp
#include "daq/hvmon/config.h"

#include <unistd.h>

#include <cstdlib>
#include <fstream>

#include "daq/hv/credentials.h"

namespace daq::hvmon {

std::string defaultRuntimeDir() {
  const char* x = std::getenv("XDG_RUNTIME_DIR");
  if (x && *x) return std::string(x) + "/daqv2";
  return "/run/user/" + std::to_string(::getuid()) + "/daqv2";
}

Result<HvmonConfig> HvmonConfig::fromJson(const nlohmann::json& j) {
  using R = Result<HvmonConfig>;
  HvmonConfig c;
  c.runtimeDir = defaultRuntimeDir();
  try {
    c.system = j.value("system", c.system);
    c.host = j.value("host", c.host);
    c.slot = j.value("slot", c.slot);
    c.channels = j.value("channels", c.channels);
    c.fastParams = j.value("fast_params", c.fastParams);
    c.fastTracked = j.value("fast_tracked", c.fastTracked);
    c.fastPeriodS = j.value("fast_period_s", c.fastPeriodS);
    c.fullPeriodS = j.value("full_period_s", c.fullPeriodS);
    c.backoffInitialS = j.value("backoff_initial_s", c.backoffInitialS);
    c.backoffMaxS = j.value("backoff_max_s", c.backoffMaxS);
    c.staleAfterS = j.value("stale_after_s", c.staleAfterS);
    if (j.contains("status_bits")) {
      auto sb = hv::StatusBitMap::fromJson(j["status_bits"]);
      if (!sb) return R(sb.error());
      c.statusBits = sb.value();
    }
    c.tripFlags = j.value("trip_flags", c.tripFlags);
    c.legacyCredentials = j.value("legacy_credentials", c.legacyCredentials);
    c.boardParamsExclude = j.value("board_params_exclude", c.boardParamsExclude);
    c.allowCommands = j.value("allow_commands", c.allowCommands);
    c.verifyMarginS = j.value("verify_margin_s", c.verifyMarginS);
    c.vmonOnToleranceV = j.value("vmon_on_tolerance_V", c.vmonOnToleranceV);
    c.vmonOffMaxV = j.value("vmon_off_max_V", c.vmonOffMaxV);
    c.commandWindowExtraS = j.value("command_window_extra_s", c.commandWindowExtraS);
    c.slowcontrolDir = j.value("slowcontrol_dir", c.slowcontrolDir);
    c.fsyncPeriodS = j.value("fsync_period_s", c.fsyncPeriodS);
    c.runtimeDir = hv::expandHome(j.value("runtime_dir", c.runtimeDir));
    c.stateDir = j.value("state_dir", c.stateDir);
    c.logFile = j.value("log_file", c.logFile);
    c.changeRingSize = j.value("change_ring_size", c.changeRingSize);
    c.writerQueueCapacity = j.value("writer_queue_capacity", c.writerQueueCapacity);
    c.gecoPattern = j.value("geco_pattern", c.gecoPattern);
    c.gecoPeriodS = j.value("geco_period_s", c.gecoPeriodS);
    if (j.contains("http")) {
      const auto& h = j["http"];
      c.httpEnabled = h.value("enabled", c.httpEnabled);
      c.httpBind = h.value("bind", c.httpBind);
      c.httpPort = h.value("port", c.httpPort);
      c.httpChannels = h.value("channels", c.httpChannels);
      c.graphPeriodS = h.value("graph_period_s", c.graphPeriodS);
      c.graphHours = h.value("graph_hours", c.graphHours);
    }
    if (j.contains("temperature")) {
      const auto& t = j["temperature"];
      c.tempBackend = t.value("backend", c.tempBackend);
      c.tempPeriodS = t.value("period_s", c.tempPeriodS);
      c.tempConfig = t.value("config", nlohmann::json::object());
    }
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("configuración de hvmon inválida: ") + e.what());
  }
  if (c.channels.empty()) return R::failure(-1, "channels vacío");
  if (c.fastPeriodS <= 0 || c.fullPeriodS <= 0 || c.fsyncPeriodS <= 0 || c.backoffInitialS <= 0 ||
      c.backoffMaxS < c.backoffInitialS || c.gecoPeriodS <= 0 || c.tempPeriodS <= 0 || c.graphPeriodS <= 0)
    return R::failure(-1, "períodos inválidos en la configuración");
  for (const char* p : {"VMon", "IMon", "Status", "Pw", "V0Set"}) {
    bool found = false;
    for (const auto& f : c.fastParams) found = found || f == p;
    if (!found) return R::failure(-1, std::string("fast_params debe incluir ") + p);
  }
  if (c.httpBind != "127.0.0.1" && c.httpBind != "localhost" && c.httpBind != "tailscale")
    return R::failure(-1, "http.bind debe ser 127.0.0.1, localhost o tailscale (nunca la red de la universidad)");
  for (const auto& f : c.tripFlags)
    if (c.statusBits.bitOf(f) < 0) return R::failure(-1, "trip_flags: bit desconocido " + f);
  return c;
}

Result<HvmonConfig> HvmonConfig::fromFile(const std::string& path) {
  std::ifstream in(hv::expandHome(path));
  if (!in) return Result<HvmonConfig>::failure(-1, "no se pudo abrir la configuración " + path);
  nlohmann::json j;
  try {
    in >> j;
  } catch (const std::exception& e) {
    return Result<HvmonConfig>::failure(-1, "JSON inválido en " + path + ": " + e.what());
  }
  return fromJson(j);
}

nlohmann::json HvmonConfig::toJson() const {
  nlohmann::json bits = nlohmann::json::object();
  for (const auto& b : statusBits.bits()) bits[b.first] = b.second;
  return nlohmann::json{{"system", system}, {"host", host}, {"slot", slot}, {"channels", channels},
                        {"fast_params", fastParams}, {"fast_tracked", fastTracked}, {"fast_period_s", fastPeriodS},
                        {"full_period_s", fullPeriodS}, {"backoff_initial_s", backoffInitialS},
                        {"backoff_max_s", backoffMaxS}, {"stale_after_s", staleAfterS}, {"status_bits", bits},
                        {"trip_flags", tripFlags}, {"board_params_exclude", boardParamsExclude}, {"allow_commands", allowCommands},
                        {"verify_margin_s", verifyMarginS}, {"vmon_on_tolerance_V", vmonOnToleranceV},
                        {"vmon_off_max_V", vmonOffMaxV}, {"slowcontrol_dir", slowcontrolDir},
                        {"fsync_period_s", fsyncPeriodS}, {"runtime_dir", runtimeDir}, {"state_dir", stateDir},
                        {"log_file", logFile}, {"geco_pattern", gecoPattern}, {"geco_period_s", gecoPeriodS},
                        {"http", {{"enabled", httpEnabled}, {"bind", httpBind}, {"port", httpPort},
                                  {"channels", httpChannels}, {"graph_period_s", graphPeriodS},
                                  {"graph_hours", graphHours}}},
                        {"temperature", {{"backend", tempBackend}, {"period_s", tempPeriodS}, {"config", tempConfig}}}};
}

}  // namespace daq::hvmon
