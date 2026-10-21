// digitizer.cpp
#include "daq/dig/digitizer.h"

#include <time.h>

namespace daq::dig {

std::int64_t monotonicNowNs() {
  timespec ts{};
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<std::int64_t>(ts.tv_sec) * 1000000000LL + ts.tv_nsec;
}

std::uint32_t SelfTriggerConfig::channelMask(int group) const {
  std::uint32_t m = 0xFF;
  for (int c : maskedChannels)
    if (c / 8 == group) m &= ~(1u << (c % 8));
  return m;
}

nlohmann::json SelfTriggerConfig::toJson() const {
  nlohmann::json masks = nlohmann::json::array(), groups = nlohmann::json::array();
  for (int g = 0; g < kGroups; ++g) {
    masks.push_back(channelMask(g));
    if ((groupMask >> g) & 1u) groups.push_back(g);
  }
  return {{"enabled", enabled},           {"groups", groups},                  {"group_mask", groupMask},          {"threshold_pe", thresholdPe},
          {"threshold_adc", thresholdAdc}, {"masked_channels", maskedChannels}, {"channel_masks", masks},
          {"majority_groups", majorityGroups}, {"majority_level_register", majorityLevel()},
          {"window_clk", windowClk},       {"window_ns", windowClk * 8},       {"include_external", includeExternal},
          {"target_rate_Hz", targetRateHz}, {"conversion", conversion}};
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
    if (j.contains("trigger")) {
      const auto& t = j["trigger"];
      c.triggerMode = t.value("mode", std::string("external"));
      if (c.triggerMode == "external") {
        c.externalTrigger = true;
        c.softwareTrigger = false;
      } else if (c.triggerMode == "software") {
        c.externalTrigger = false;
        c.softwareTrigger = true;
      } else if (c.triggerMode == "self_hod2018") {
        const auto s = t.value("self", nlohmann::json::object());
        auto& st = c.self;
        st.enabled = true;
        for (int g : s.value("groups", std::vector<int>{0, 1, 2, 3})) {
          if (g < 0 || g >= kGroups) return R::failure(-1, "trigger.self.groups fuera de 0-7");
          st.groupMask |= 1u << g;
        }
        if (s.contains("threshold_pe")) {
          if (s["threshold_pe"].is_number())
            st.thresholdPe.assign(kGroups, s["threshold_pe"].get<double>());
          else
            st.thresholdPe = s["threshold_pe"].get<std::vector<double>>();
        }
        st.thresholdAdc = s.value("threshold_adc", st.thresholdAdc);
        st.maskedChannels = s.value("masked_channels", st.maskedChannels);
        st.majorityGroups = s.value("majority_groups", st.majorityGroups);
        st.windowClk = s.value("window_clk", st.windowClk);
        st.includeExternal = s.value("include_external", st.includeExternal);
        st.targetRateHz = s.value("target_rate_Hz", st.targetRateHz);
        c.externalTrigger = st.includeExternal;
        c.softwareTrigger = false;
        int nGroups = 0;
        for (int g = 0; g < kGroups; ++g) nGroups += (st.groupMask >> g) & 1u;
        if (nGroups == 0) return R::failure(-1, "trigger.self: ningún grupo");
        if (st.majorityGroups < 1 || st.majorityGroups > nGroups || st.majorityGroups > 8)
          return R::failure(-1, "trigger.self.majority_groups debe estar entre 1 y el número de grupos");
        if (st.windowClk < 0 || st.windowClk > 15) return R::failure(-1, "trigger.self.window_clk fuera de 0-15 (0x810C[23:20])");
        if (st.thresholdPe.empty() && st.thresholdAdc.empty()) return R::failure(-1, "trigger.self: falta threshold_pe o threshold_adc");
        if (!st.thresholdPe.empty() && static_cast<int>(st.thresholdPe.size()) != kGroups)
          return R::failure(-1, "trigger.self.threshold_pe: un número o 8 valores");
        if (!st.thresholdAdc.empty() && static_cast<int>(st.thresholdAdc.size()) != kGroups)
          return R::failure(-1, "trigger.self.threshold_adc: 8 valores");
        for (int a : st.thresholdAdc)
          if (a < 0 || a > 4095) return R::failure(-1, "trigger.self.threshold_adc fuera de 0-4095 (0x1n80[11:0])");
        for (int ch : st.maskedChannels)
          if (ch < 0 || ch >= kChannels) return R::failure(-1, "trigger.self.masked_channels fuera de 0-63");
        if (st.targetRateHz <= 0) return R::failure(-1, "trigger.self.target_rate_Hz debe ser > 0");
      } else {
        return R::failure(-1, "trigger.mode = " + c.triggerMode + ": debe ser external, software o self_hod2018");
      }
    }
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
          {"trigger", {{"mode", triggerMode.empty() ? (softwareTrigger ? "software" : "external") : triggerMode}, {"self", self.toJson()}}},
          {"io_level", ioLevel},
          {"group_enable_mask", groupEnableMask},
          {"link", {{"type", linkType}, {"link_num", linkNum}, {"conet_node", conetNode}, {"vme_base_address", vmeBaseAddress}}},
          {"reset_on_configure", resetOnConfigure},
          {"capture_dir", captureDir},
          {"capture_blocks", captureBlocks},
          {"scenario", scenario}};
}

}  // namespace daq::dig
