// test_self_trigger.cpp — DAQ2-13: modo de trigger, conversión de umbrales p.e. -> ADC y aviso de tasa alta.
#include <gtest/gtest.h>

#include <cmath>

#include "daq/daqrun/calibration.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/self_trigger.h"
#include "daq/dig/digitizer.h"

using namespace daq;
using namespace daq::daqrun;

namespace {
nlohmann::json selfCfg() {
  return {{"backend", "fake"},
          {"trigger",
           {{"mode", "self_hod2018"},
            {"self",
             {{"groups", {0, 1, 2, 3}}, {"threshold_pe", 10.0}, {"masked_channels", {2, 9, 15}}, {"majority_groups", 2}, {"window_clk", 1},
              {"target_rate_Hz", 50.0}}}}}};
}
}  // namespace

TEST(TriggerMode, DefaultExternalAndModesSetTheTriggerSources) {
  auto d = dig::DigitizerConfig::fromJson({{"backend", "fake"}}).value();
  EXPECT_TRUE(d.externalTrigger);
  EXPECT_FALSE(d.self.enabled);
  auto s = dig::DigitizerConfig::fromJson({{"backend", "fake"}, {"trigger", {{"mode", "software"}}}}).value();
  EXPECT_TRUE(s.softwareTrigger);
  EXPECT_FALSE(s.externalTrigger);
  auto a = dig::DigitizerConfig::fromJson(selfCfg()).value();
  EXPECT_TRUE(a.self.enabled);
  EXPECT_FALSE(a.externalTrigger);
  EXPECT_FALSE(a.softwareTrigger);
  EXPECT_EQ(a.self.groupMask, 0x0Fu);
  EXPECT_EQ(a.self.majorityLevel(), 1);
  EXPECT_EQ(a.self.channelMask(0), 0xFBu);  // ch2
  EXPECT_EQ(a.self.channelMask(1), 0x7Du);  // ch9, ch15
  EXPECT_EQ(a.self.channelMask(2), 0xFFu);
  // ida y vuelta por toJson
  auto b = dig::DigitizerConfig::fromJson(a.toJson()).value();
  EXPECT_EQ(b.self.groupMask, a.self.groupMask);
  EXPECT_EQ(b.self.maskedChannels, a.self.maskedChannels);
  EXPECT_FALSE(dig::DigitizerConfig::fromJson({{"backend", "fake"}, {"trigger", {{"mode", "otro"}}}}));
  auto bad = selfCfg();
  bad["trigger"]["self"]["majority_groups"] = 5;
  EXPECT_FALSE(dig::DigitizerConfig::fromJson(bad));
  bad = selfCfg();
  bad["trigger"]["self"]["window_clk"] = 16;
  EXPECT_FALSE(dig::DigitizerConfig::fromJson(bad));
}

TEST(TriggerMode, ThresholdFromCalibrationUsesMedianBaselineAndGain) {
  auto d = dig::DigitizerConfig::fromJson(selfCfg()).value();
  PeCalibration cal;
  cal.version = "test";
  for (int c = 0; c < 64; ++c) {
    cal.gainMv[static_cast<std::size_t>(c)] = c < 32 ? 4.0 : 0.0;
    cal.pedestalMv[static_cast<std::size_t>(c)] = 25.0;
  }
  cal.pedestalMv[1] = 35.0;  // un canal con la línea base 10 mV más arriba
  ASSERT_TRUE(applySelfTriggerThresholds(d, cal));
  dig::AdcConversion ac{12, 2.0, 32767};
  const int expected = static_cast<int>(std::lround(ac.toAdc(25.0) - 10.0 * 4.0 / ac.lsbMv()));
  EXPECT_EQ(d.self.thresholdAdc[0], expected);
  const auto& eff = d.self.conversion["0"]["effective_threshold_pe"];
  EXPECT_NEAR(eff["0"].get<double>(), 10.0, 0.2);
  EXPECT_NEAR(eff["1"].get<double>(), 12.5, 0.2);  // 10 mV más de línea base = 2.5 p.e. más
  EXPECT_FALSE(eff.contains("2"));                  // enmascarado
  EXPECT_EQ(d.self.thresholdAdc[4], 0);             // grupo que no participa
}

TEST(TriggerMode, RateHighAlertOncePerEpisode) {
  const UtcNs t0 = 1'000'000'000'000'000'000LL;
  TriggerWatchdog w(300, "test", t0, 0.1, 60);
  w.setLabel("del auto-trigger de hod2018");
  w.setMaxRate(100);
  UtcNs t = t0;
  for (int i = 0; i < 59; ++i) {
    t += 1'000'000'000LL;
    w.onTriggers(t, 150, std::nullopt);  // 150 Hz
    EXPECT_FALSE(w.tickHigh(t, std::nullopt));  // ventana incompleta
  }
  t += 2'000'000'000LL;
  w.onTriggers(t, 300, std::nullopt);
  auto e = w.tickHigh(t, std::nullopt);
  ASSERT_TRUE(e);
  EXPECT_EQ(e->kind, "trigger_rate_high");
  EXPECT_FALSE(w.tickHigh(t, std::nullopt));
  EXPECT_TRUE(w.highAlerted());
}
