// Tests de utilidades de tiempo y de los relojes.
#include <gtest/gtest.h>

#include <chrono>

#include "daq/core/clock.h"
#include "daq/core/time.h"

using namespace std::chrono_literals;

TEST(Time, FormatCompactAndIso) {
  // 2026-10-18T21:56:18.617305369Z
  const daq::UtcNs t = 1792360578617305369LL;
  EXPECT_EQ(daq::formatUtcCompact(t), "20261018T215618Z");
  EXPECT_EQ(daq::formatUtcIso(t), "2026-10-18T21:56:18.617305369Z");
  EXPECT_EQ(daq::formatUtcIso(0), "1970-01-01T00:00:00.000000000Z");
  EXPECT_EQ(daq::formatUtcIso(-1), "1969-12-31T23:59:59.999999999Z");
}

TEST(Time, UtcNowIsPlausible) {
  const daq::UtcNs t = daq::utcNowNs();
  EXPECT_GT(t, 1792281600LL * 1000000000LL);  // posterior a 2026-10-18
}

TEST(ManualClock, SleepUntilAdvancesDeterministically) {
  daq::ManualClock c;
  const auto s0 = c.now_steady();
  const daq::UtcNs u0 = c.now_utc_ns();
  EXPECT_EQ(daq::formatUtcIso(u0), "2026-10-18T00:00:00.000000000Z");
  EXPECT_TRUE(c.sleep_until(s0 + 65750ms));
  EXPECT_EQ(c.now_steady() - s0, 65750ms);
  EXPECT_EQ(c.now_utc_ns() - u0, 65750LL * 1000000);
  // Un instante pasado no retrocede el reloj.
  EXPECT_TRUE(c.sleep_until(s0 + 1s));
  EXPECT_EQ(c.now_steady() - s0, 65750ms);
  c.advance(250ms);
  EXPECT_EQ(c.now_steady() - s0, 66s);
  EXPECT_TRUE(c.sleep_for(4s));
  EXPECT_EQ(c.now_steady() - s0, 70s);
}

TEST(ManualClock, SleepReturnsFalseWhenStopped) {
  daq::ManualClock c;
  daq::StopSource s;
  s.request_stop();
  const auto s0 = c.now_steady();
  EXPECT_FALSE(c.sleep_until(s0 + 10s, s.get_token()));
  EXPECT_EQ(c.now_steady(), s0);
}

TEST(SystemClock, SleepForWaitsAtLeastTheDuration) {
  daq::SystemClock c;
  const auto t0 = c.now_steady();
  EXPECT_TRUE(c.sleep_for(20ms));
  EXPECT_GE(c.now_steady() - t0, 20ms);
}
