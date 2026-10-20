// dash.cpp
#include "daq/daqrun/dash.h"

#include <algorithm>
#include <cmath>

namespace daq::daqrun {

nlohmann::json dashJson(const DashInputs& in) {
  const MonitorView& v = *in.view;
  nlohmann::json j{{"campaign", in.campaign}, {"run", in.runName}, {"simulated", in.simulated}, {"t_utc_ns", in.nowUtcNs},
                   {"elapsed_s", static_cast<double>(in.nowUtcNs - in.runStartUtcNs) / 1e9}, {"events", v.events}};
  j["rate_1min_Hz"] = v.recent.empty() ? nlohmann::json() : nlohmann::json(v.last.triggerRateHz);
  std::uint64_t nEv = 0;
  double span = 0;
  for (const auto& s : v.recent) {
    nEv += s.nEvents;
    span += static_cast<double>(s.tEndUtcNs - s.tStartUtcNs) / 1e9;
  }
  j["rate_10min_Hz"] = span > 0 ? nlohmann::json(static_cast<double>(nEv) / span) : nlohmann::json();
  j["window_s"] = span;
  const double since = v.events > 0 ? static_cast<double>(in.nowUtcNs - v.lastEventUtcNs) / 1e9 : -1.0;
  j["since_last_trigger_s"] = v.events > 0 ? nlohmann::json(since) : nlohmann::json();
  j["disk_free_GB"] = in.diskFreeBytes >= 0 ? nlohmann::json(in.diskFreeBytes / 1e9) : nlohmann::json();
  const HvSample* hv = in.hv;
  const double hvAge = hv && hv->stateUtcNs > 0 ? static_cast<double>(in.nowUtcNs - hv->stateUtcNs) / 1e9 : -1.0;
  j["hvmon"] = {{"alive", hv && hv->alive}, {"stale", !hv || hv->stale}, {"age_s", hvAge}};
  nlohmann::json stations = nlohmann::json::array();
  for (int st = 0; st < 2; ++st) {
    std::vector<double> rates;
    int live = 0;
    for (int c = 32 * st; c < 32 * st + 32; ++c) {
      const auto ci = static_cast<std::size_t>(c);
      std::uint64_t low = 0, cnt = 0;
      double t = 0;
      for (const auto& s : v.recent) {
        low += s.lowCounts[ci];
        cnt += s.ootCounts[ci];
        t += s.ootTimeS[ci];
      }
      if (low == 0) continue;
      ++live;
      if (t > 0) rates.push_back(static_cast<double>(cnt) / t);
    }
    nlohmann::json med;
    if (!rates.empty()) {
      std::sort(rates.begin(), rates.end());
      const std::size_t n = rates.size();
      med = n % 2 ? rates[n / 2] : 0.5 * (rates[n / 2 - 1] + rates[n / 2]);
    }
    stations.push_back({{"name", st == 0 ? "hod2018" : "hod2019"}, {"live", live}, {"channels", 32},
                        {"median_oot_Hz", med}, {"threshold_mV", in.thresholdMv[static_cast<std::size_t>(st)]}});
  }
  j["stations"] = stations;
  nlohmann::json alerts = nlohmann::json::array();
  if (!hv || !hv->alive || hv->stale) alerts.push_back("hvmon no vivo o stale: bit 4 (HV_UNKNOWN) en todos los canales");
  if (in.watchTriggers && (in.noTriggerAlert || (v.events > 0 && since > in.noTriggerAlertS) ||
                           (v.events == 0 && static_cast<double>(in.nowUtcNs - in.runStartUtcNs) / 1e9 > in.noTriggerAlertS)))
    alerts.push_back("sin triggers externos (o tasa bajo el mínimo)");
  const auto& r = in.readout;
  if (r.value("read_errors", 0) > 0) alerts.push_back("errores de lectura del digitizador: " + std::to_string(r.value("read_errors", 0)));
  if (r.value("board_full_blocks", 0) > 0) alerts.push_back("memoria de la placa llena en " + std::to_string(r.value("board_full_blocks", 0)) + " bloques");
  if (r.value("writer_blocked", 0) > 0) alerts.push_back("cola de escritura llena " + std::to_string(r.value("writer_blocked", 0)) + " veces");
  if (!r.value("last_error", std::string()).empty()) alerts.push_back("último error: " + r.value("last_error", std::string()));
  if (in.diskFreeBytes >= 0 && in.diskFreeBytes < in.diskLowBytes) alerts.push_back("poco disco libre");
  j["alerts"] = alerts;
  return j;
}

}  // namespace daq::daqrun
