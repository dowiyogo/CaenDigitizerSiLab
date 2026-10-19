// digitizer.cpp
#include "daq/dig/digitizer.h"

#include <time.h>

namespace daq::dig {

std::int64_t monotonicNowNs() {
  timespec ts{};
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<std::int64_t>(ts.tv_sec) * 1000000000LL + ts.tv_nsec;
}

Result<DigitizerConfig> DigitizerConfig::fromJson(const nlohmann::json& j) {
  using R = Result<DigitizerConfig>;
  DigitizerConfig c;
  try {
    c.backend = j.value("backend", c.backend);
    c.recordLength = j.value("record_length_samples", c.recordLength);
    c.postTriggerRegister = j.value("post_trigger_register", c.postTriggerRegister);
    c.groupDcOffset = j.value("group_dc_offset", c.groupDcOffset);
    c.maxEventsBlt = j.value("max_events_blt", c.maxEventsBlt);
    c.pollSleepS = j.value("poll_sleep_s", c.pollSleepS);
    c.tttTickNs = j.value("ttt_tick_ns", c.tttTickNs);
    c.nBit = j.value("n_bit", c.nBit);
    c.maxVpp = j.value("max_vpp_V", c.maxVpp);
    c.softwareTrigger = j.value("software_trigger", c.softwareTrigger);
    c.swTriggerRateHz = j.value("sw_trigger_rate_Hz", c.swTriggerRateHz);
    c.scenarioFile = j.value("scenario_file", c.scenarioFile);
    if (j.contains("scenario")) c.scenario = j["scenario"];
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("digitizer inválido: ") + e.what());
  }
  if (c.recordLength <= 0 || c.recordLength % 3 != 0) return R::failure(-1, "record_length_samples debe ser múltiplo de 3 y > 0");
  if (static_cast<int>(c.groupDcOffset.size()) != kGroups) return R::failure(-1, "group_dc_offset debe tener 8 valores");
  if (c.maxEventsBlt <= 0 || c.pollSleepS <= 0 || c.tttTickNs <= 0) return R::failure(-1, "digitizer: parámetros de lectura inválidos");
  if (c.backend != "fake" && c.backend != "replay")
    return R::failure(-1, "digitizer.backend = " + c.backend + ": en la fase 3 solo hay 'fake' y 'replay'");
  return c;
}

nlohmann::json DigitizerConfig::toJson() const {
  return {{"backend", backend},
          {"record_length_samples", recordLength},
          {"post_trigger_register", postTriggerRegister},
          {"group_dc_offset", groupDcOffset},
          {"max_events_blt", maxEventsBlt},
          {"poll_sleep_s", pollSleepS},
          {"ttt_tick_ns", tttTickNs},
          {"n_bit", nBit},
          {"max_vpp_V", maxVpp},
          {"software_trigger", softwareTrigger},
          {"sw_trigger_rate_Hz", swTriggerRateHz},
          {"scenario_file", scenarioFile},
          {"scenario", scenario}};
}

}  // namespace daq::dig
