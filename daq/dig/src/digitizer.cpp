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
    c.tttTickSource = j.value("ttt_tick_source", c.tttTickSource);
    c.nBit = j.value("n_bit", c.nBit);
    c.maxVpp = j.value("max_vpp_V", c.maxVpp);
    c.softwareTrigger = j.value("software_trigger", c.softwareTrigger);
    c.swTriggerRateHz = j.value("sw_trigger_rate_Hz", c.swTriggerRateHz);
    c.scenarioFile = j.value("scenario_file", c.scenarioFile);
    c.externalTrigger = j.value("external_trigger", c.externalTrigger);
    c.ioLevel = j.value("io_level", c.ioLevel);
    if (j.contains("group_enable_mask"))
      c.groupEnableMask = j["group_enable_mask"].is_string()
                              ? static_cast<std::uint32_t>(std::stoul(j["group_enable_mask"].get<std::string>(), nullptr, 0))
                              : j["group_enable_mask"].get<std::uint32_t>();
    if (j.contains("link")) {
      const auto& l = j["link"];
      c.linkType = l.value("type", c.linkType);
      c.linkNum = l.value("link_num", c.linkNum);
      c.conetNode = l.value("conet_node", c.conetNode);
      if (l.contains("vme_base_address"))
        c.vmeBaseAddress = l["vme_base_address"].is_string()
                               ? static_cast<std::uint32_t>(std::stoul(l["vme_base_address"].get<std::string>(), nullptr, 0))
                               : l["vme_base_address"].get<std::uint32_t>();
    }
    c.resetOnConfigure = j.value("reset_on_configure", c.resetOnConfigure);
    c.captureDir = j.value("capture_dir", c.captureDir);
    c.captureBlocks = j.value("capture_blocks", c.captureBlocks);
    if (j.contains("scenario")) c.scenario = j["scenario"];
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("digitizer inválido: ") + e.what());
  }
  if (c.recordLength <= 0 || c.recordLength > 100000) return R::failure(-1, "record_length_samples fuera de rango");
  if (c.recordLength % 3 != 0 && c.backend != "caen")
    return R::failure(-1, "record_length_samples debe ser múltiplo de 3 (la placa redondea hacia abajo)");
  if (c.ioLevel != "NIM" && c.ioLevel != "TTL") return R::failure(-1, "io_level debe ser NIM o TTL");
  if (c.groupEnableMask == 0 || c.groupEnableMask > 0xFF) return R::failure(-1, "group_enable_mask fuera de rango");
  if (static_cast<int>(c.groupDcOffset.size()) != kGroups) return R::failure(-1, "group_dc_offset debe tener 8 valores");
  if (c.maxEventsBlt <= 0 || c.pollSleepS <= 0 || c.tttTickNs <= 0) return R::failure(-1, "digitizer: parámetros de lectura inválidos");
  if (c.backend != "fake" && c.backend != "replay" && c.backend != "caen")
    return R::failure(-1, "digitizer.backend = " + c.backend + ": debe ser fake, replay o caen");
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
          {"ttt_tick_source", tttTickSource},
          {"n_bit", nBit},
          {"max_vpp_V", maxVpp},
          {"software_trigger", softwareTrigger},
          {"sw_trigger_rate_Hz", swTriggerRateHz},
          {"scenario_file", scenarioFile},
          {"external_trigger", externalTrigger},
          {"io_level", ioLevel},
          {"group_enable_mask", groupEnableMask},
          {"link", {{"type", linkType}, {"link_num", linkNum}, {"conet_node", conetNode}, {"vme_base_address", vmeBaseAddress}}},
          {"reset_on_configure", resetOnConfigure},
          {"capture_dir", captureDir},
          {"capture_blocks", captureBlocks},
          {"scenario", scenario}};
}

}  // namespace daq::dig
