// Tests del FakeDigitizer sintético y del bucle de lectura (reloj manual, sin hardware).
#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <thread>

#include "daq/core/clock.h"
#include "daq/dig/fake_digitizer.h"
#include "daq/dig/readout.h"
#include "daq/pulse_finder.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::dig;

namespace {

DigitizerConfig digCfg(bool sw = false) {
  DigitizerConfig c;
  c.softwareTrigger = sw;
  return c;
}

SyntheticParams params(nlohmann::json extra = nlohmann::json::object()) {
  nlohmann::json j{{"seed", 3}, {"trigger_rate_Hz", 5.0}};
  for (auto it = extra.begin(); it != extra.end(); ++it) j[it.key()] = it.value();
  auto p = SyntheticParams::fromJson(j);
  EXPECT_TRUE(p.ok()) << (p.ok() ? "" : p.error().message);
  return p.value();
}

std::vector<ChannelInfo> withGroups() {
  auto v = defaultChannelInfo();
  for (int c = 0; c < 64; ++c) v[static_cast<std::size_t>(c)].hvGroup = 1 + (c % 8);
  for (int c = 32; c < 48; c += 2) v[static_cast<std::size_t>(c)].face = "Delantera";
  return v;
}

}  // namespace

TEST(FakeDigitizer, PoissonRateCounterTttAndNoise) {
  ManualClock clock;
  FakeDigitizer fd(params(), clock);
  ASSERT_TRUE(fd.configure(digCfg()).ok());
  auto o = fd.start();
  ASSERT_TRUE(o.ok());
  std::vector<Event> all;
  for (int i = 0; i < 2000; ++i) {
    clock.advance(1s);
    EventBlock b;
    ASSERT_TRUE(fd.readBlock(b, 1000).ok());
    for (auto& e : b.events) all.push_back(std::move(e));
  }
  // 2000 s a 5 Hz: 10000 ± 300 (3 σ).
  EXPECT_NEAR(static_cast<double>(all.size()), 10000.0, 300.0);
  for (std::size_t i = 1; i < all.size(); ++i) EXPECT_EQ(all[i].counterHw, all[i - 1].counterHw + 1);
  // TTT en ticks de 8 ns, vuelta cada 2^31 ticks (17.18 s): hay vueltas en 2000 s.
  int wraps = 0;
  for (std::size_t i = 1; i < all.size(); ++i) wraps += all[i].tttRaw < all[i - 1].tttRaw;
  EXPECT_GE(wraps, 100);
  // Ruido de hod2019 (0.25 mV): desviación de las muestras del pedestal, sin los pulsos > 3 mV.
  const double lsb = pulse::lsbMv();
  double s = 0, s2 = 0;
  int n = 0;
  for (std::size_t i = 0; i < 200; ++i)
    for (int c = 40; c < 44; ++c) {
      const std::uint16_t* x = all[i].adc.data() + c * 144;
      const double m = pulse::median(x, {0, 40});
      for (int k = 0; k < 40; ++k)
        if (std::fabs(x[k] - m) * lsb < 2.0) {
          s += (x[k] - m) * lsb;
          s2 += (x[k] - m) * (x[k] - m) * lsb * lsb;
          ++n;
        }
    }
  const double sigma = std::sqrt(s2 / n - (s / n) * (s / n));
  EXPECT_NEAR(sigma, 0.25, 0.08);  // incluye la cuantización (0.488/√12 = 0.14 mV)
}

TEST(FakeDigitizer, HvGroupOffLightLeakGapBurstAndMemoryFull) {
  ManualClock clock;
  auto p = params({{"trigger_rate_Hz", 2.0},
                   {"board_memory_events", 50},
                   {"events",
                    {{{"t_s", 100}, {"type", "hv_group_off"}, {"group", 3}, {"duration_s", 100}},
                     {{"t_s", 300}, {"type", "trigger_gap"}, {"duration_s", 60}},
                     {{"t_s", 500}, {"type", "trigger_burst"}, {"n", 40}, {"duration_s", 1.0}},
                     {{"t_s", 600}, {"type", "light_leak"}, {"face", "Delantera"}, {"duration_s", 100}, {"extra_rate_Hz", 2e6}}}}});
  FakeDigitizer fd(p, clock, withGroups());
  ASSERT_TRUE(fd.configure(digCfg()).ok());
  const UtcNs t0 = fd.start().value().utcBeforeStart;
  std::vector<std::pair<UtcNs, Event>> evs;
  for (int i = 0; i < 800; ++i) {
    clock.advance(1s);
    EventBlock b;
    ASSERT_TRUE(fd.readBlock(b, 1000).ok());
    for (auto& e : b.events) evs.emplace_back(b.tReadUtcNs, std::move(e));
  }
  // Grupo 3 apagado entre 100 y 200 s: sus canales sin pulsos (> 10 mV) en ese intervalo.
  const double lsb = pulse::lsbMv();
  auto pulsesOn = [&](const Event& e, int c, double thr) {
    const std::uint16_t* x = e.adc.data() + c * 144;
    return pulse::countPulses(x, 144, pulse::trimmedMean(x, {0, 40}, 2.0 / lsb), lsb, thr, 4, {{0, 144}});
  };
  int offPulses = 0, onPulses = 0, offEvents = 0;
  for (const auto& [t, e] : evs) {
    const double s = static_cast<double>(t - t0) / 1e9;
    if (s > 101 && s < 199) {
      ++offEvents;
      for (int c = 2; c < 32; c += 8) offPulses += pulsesOn(e, c, 10.0);  // canales del grupo 3 (c % 8 == 2)
      for (int c = 3; c < 32; c += 8) onPulses += pulsesOn(e, c, 10.0);   // grupo 4, encendido
    }
  }
  EXPECT_GT(offEvents, 100);
  EXPECT_EQ(offPulses, 0);
  EXPECT_GT(onPulses, 100);
  EXPECT_TRUE(fd.channelOff(2, t0 + 150 * 1000000000LL));
  EXPECT_FALSE(fd.channelOff(2, t0 + 250 * 1000000000LL));
  // Hueco 300-360 s: ningún trigger leído con tRead en (301, 360].
  int inGap = 0, inBurst = 0;
  for (const auto& [t, e] : evs) {
    const double s = static_cast<double>(t - t0) / 1e9;
    inGap += s > 301 && s <= 360;
    inBurst += s > 500 && s <= 502;
  }
  EXPECT_EQ(inGap, 0);
  EXPECT_GE(inBurst, 40);
  // Fuga de luz en la cara Delantera (ch 32, 34, ...): más pulsos de hod2019 > 3 mV que en la Trasera.
  int leak = 0, other = 0;
  for (const auto& [t, e] : evs) {
    const double s = static_cast<double>(t - t0) / 1e9;
    if (s > 601 && s < 699) {
      leak += pulsesOn(e, 34, 3.0);
      other += pulsesOn(e, 35, 3.0);
    }
  }
  EXPECT_GT(leak, 5 * std::max(other, 1));
  // Memoria de 50 eventos: sin leer durante 100 s a 2 Hz se pierden triggers y se marca board_full.
  clock.advance(100s);
  EventBlock b;
  ASSERT_TRUE(fd.readBlock(b, 1000).ok());
  EXPECT_TRUE(b.boardFull);
  EXPECT_EQ(b.events.size(), 50u);
  EXPECT_GT(fd.eventsLostBoardFull(), 0u);
  bool jump = false;
  for (std::size_t i = 1; i < b.events.size(); ++i) jump = jump || b.events[i].counterHw != b.events[i - 1].counterHw + 1;
  // Los perdidos son los más nuevos (memoria llena): el salto aparece con el bloque siguiente.
  clock.advance(10s);  // varios triggers nuevos (a 2 Hz, 1 s puede no traer ninguno)
  EventBlock b2;
  ASSERT_TRUE(fd.readBlock(b2, 1000).ok());
  if (!b2.events.empty()) jump = jump || b2.events.front().counterHw != b.events.back().counterHw + 1;
  EXPECT_TRUE(jump);
}

TEST(FakeDigitizer, SoftwareTriggerModeAndGainStep) {
  ManualClock clock;
  auto p = params({{"events", {{{"t_s", 10}, {"type", "gain_step"}, {"station", "hod2018"}, {"factor", 0.5}}}}});
  FakeDigitizer fd(p, clock);
  ASSERT_TRUE(fd.configure(digCfg(true)).ok());
  const UtcNs t0 = fd.start().value().utcBeforeStart;
  clock.advance(5s);
  EventBlock b;
  ASSERT_TRUE(fd.readBlock(b, 100).ok());
  EXPECT_TRUE(b.events.empty());  // sin triggers externos en modo calibración
  for (int i = 0; i < 10; ++i) ASSERT_TRUE(fd.sendSoftwareTrigger().ok());
  ASSERT_TRUE(fd.readBlock(b, 100).ok());
  ASSERT_EQ(b.events.size(), 10u);
  for (const auto& e : b.events) EXPECT_EQ(e.triggerKind, 1);
  EXPECT_NEAR(fd.gainMv(0, t0 + 20 * 1000000000LL) / fd.gainMv(0, t0), 0.5, 1e-12);
  EXPECT_NEAR(fd.gainMv(40, t0 + 20 * 1000000000LL) / fd.gainMv(40, t0), 1.0, 1e-12);
  EXPECT_FALSE(SyntheticParams::fromJson({{"events", {{{"t_s", 1}, {"type", "otra_cosa"}}}}}).ok());
}

TEST(Readout, FullWriterQueueBlocksReadoutWithoutLossMonitorDrops) {
  ManualClock clock;
  FakeDigitizer fd(params({{"trigger_rate_Hz", 50.0}}), clock);
  ASSERT_TRUE(fd.configure(digCfg()).ok());
  ASSERT_TRUE(fd.start().ok());
  BoundedQueue<EventBlockPtr> writerQ(4), monitorQ(2);
  ReadoutStats st;
  ReadoutOptions o;
  o.maxEventsBlt = 8;
  o.maxEvents = 2000;
  ReadoutLoop loop(fd, clock, writerQ, monitorQ, o, st);
  std::atomic<std::uint64_t> received{0};
  JThread rt([&](StopToken s) { loop.run(s); });
  // El consumidor arranca tarde: la cola de escritura se llena y la lectura espera.
  std::this_thread::sleep_for(200ms);
  std::uint32_t last = 0;
  bool contiguous = true, first = true;
  while (auto b = writerQ.pop()) {
    for (const auto& e : (*b)->events) {
      if (!first && e.counterHw != last + 1) contiguous = false;
      last = e.counterHw;
      first = false;
      ++received;
    }
  }
  rt.join();
  EXPECT_TRUE(st.finished.load());
  EXPECT_EQ(st.finishReason(), "max_events");
  EXPECT_EQ(received.load(), st.events.load());  // nada se perdió en el camino del escritor
  EXPECT_GE(received.load(), 2000u);
  EXPECT_TRUE(contiguous);
  EXPECT_GT(st.writerBlocked.load(), 0u);
  EXPECT_EQ(st.writerMaxOccupancy.load(), 4u);
  EXPECT_GT(st.monitorDropped.load(), 0u);  // el monitor nunca consumió
  EXPECT_EQ(fd.eventsLostBoardFull(), 0u);
  EXPECT_TRUE(writerQ.closed());
  EXPECT_TRUE(monitorQ.closed());
}
