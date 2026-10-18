// Tests de FakeHVBackend con tiempo acelerado (ManualClock).
#include <gtest/gtest.h>

#include <chrono>
#include <cstdio>
#include <fstream>

#include "daq/core/clock.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hv/status_bits.h"

using namespace std::chrono_literals;
using daq::ManualClock;
using daq::hv::FakeHVBackend;
using daq::hv::StatusBitMap;

namespace {

// Guion de un canal sin ruido ni offsets, para comprobar tiempos exactos.
nlohmann::json quietScenario(nlohmann::json events) {
  auto sc = FakeHVBackend::defaultScenario();
  auto& s = sc["slots"][0];
  s["imon_noise_uA"] = 0.0;
  s["vmon_noise_V"] = 0.0;
  s["imon_offset_uA"] = 0.0;
  s["vmon_offset_V"] = 0.0;
  sc["events"] = std::move(events);
  return sc;
}

struct Fixture {
  ManualClock clock;
  FakeHVBackend be;
  StatusBitMap bits = StatusBitMap::wrapperRev13();
  explicit Fixture(const nlohmann::json& sc) : be(sc, clock) { EXPECT_TRUE(be.connect().ok()); }
  void at(std::chrono::milliseconds t) { clock.sleep_until(daq::SteadyTime{} + t); }
  std::uint32_t status(int ch) { return be.readChannels(3, {ch}, {"Status"}).value().values[0][0].raw; }
  double param(int ch, const char* p) { return be.readChannels(3, {ch}, {p}).value().values[0][0].si; }
  std::uint32_t pw(int ch) { return be.readChannels(3, {ch}, {"Pw"}).value().values[0][0].raw; }
};

}  // namespace

TEST(FakeHV, ReproducesA7030DPParameterSet) {
  ManualClock clock;
  FakeHVBackend be(FakeHVBackend::defaultScenario(), clock);
  ASSERT_TRUE(be.connect().ok());
  auto p = be.channelParamNames(3, 0);
  ASSERT_TRUE(p.ok());
  const std::vector<std::pair<std::string, int>> expected = {
      {"V0Set", 0}, {"I0Set", 0}, {"V1Set", 0}, {"I1Set", 0}, {"RUp", 0}, {"RDWn", 0},
      {"Trip", 0}, {"SVMax", 0}, {"VMon", 0}, {"IMon", 0}, {"Status", 2}, {"Pw", 1},
      {"POn", 1}, {"PDwn", 1}, {"TripInt", 4}, {"TripExt", 4}, {"ZCDetect", 1}, {"ZCAdjust", 1}};
  ASSERT_EQ(p.value().size(), expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    EXPECT_EQ(p.value()[i].name, expected[i].first);
    EXPECT_EQ(static_cast<int>(p.value()[i].type), expected[i].second) << expected[i].first;
  }
  auto m = be.crateMap();
  ASSERT_TRUE(m.ok());
  EXPECT_EQ(m.value().slots.at(3).model, "A7030DP");
  EXPECT_EQ(m.value().slots.at(3).channels, 12);
  EXPECT_EQ(m.value().slots.size(), 4u);
  EXPECT_EQ(m.value().slots.at(3).serial, 58015u);
}

TEST(FakeHV, DefaultScenarioMatchesReadOf20261018) {
  ManualClock clock;
  FakeHVBackend be(FakeHVBackend::defaultScenario(), clock);
  ASSERT_TRUE(be.connect().ok());
  auto r = be.readChannels(3, {0, 6, 8}, {"V0Set", "I0Set", "Trip", "Pw", "TripExt"});
  ASSERT_TRUE(r.ok());
  EXPECT_DOUBLE_EQ(r.value().get("V0Set", 0)->boardValue, 65.75);
  EXPECT_DOUBLE_EQ(r.value().get("V0Set", 6)->boardValue, 66.15);
  EXPECT_NEAR(r.value().get("I0Set", 0)->si, 15e-6, 1e-15);  // SI: A
  EXPECT_DOUBLE_EQ(r.value().get("I0Set", 8)->boardValue, 5.0);
  EXPECT_DOUBLE_EQ(r.value().get("Trip", 8)->si, 10.0);
  EXPECT_EQ(r.value().get("Pw", 0)->raw, 0u);
  EXPECT_EQ(r.value().get("TripExt", 6)->raw, 16u);
}

TEST(FakeHV, RampFromZeroTo6575VAt1VPerSecondTakes6575s) {
  Fixture f(quietScenario({{{"t_s", 10.0}, {"type", "set"}, {"ch", 0}, {"param", "Pw"}, {"value", 1}}}));
  f.at(9s);
  EXPECT_EQ(f.pw(0), 0u);
  f.at(10s);
  EXPECT_EQ(f.pw(0), 1u);
  EXPECT_TRUE(f.bits.has(f.status(0), "ramp_up"));
  f.at(10s + 32875ms);  // mitad de la rampa
  EXPECT_NEAR(f.be.trueVoltage(3, 0), 32.875, 1e-6);
  f.at(10s + 65700ms);
  EXPECT_TRUE(f.bits.has(f.status(0), "ramp_up"));
  EXPECT_NEAR(f.be.trueVoltage(3, 0), 65.70, 1e-6);
  f.at(10s + 65750ms);
  EXPECT_NEAR(f.be.trueVoltage(3, 0), 65.75, 1e-9);
  EXPECT_FALSE(f.bits.has(f.status(0), "ramp_up"));
  EXPECT_TRUE(f.bits.has(f.status(0), "on"));
  EXPECT_NEAR(f.param(0, "VMon"), 65.75, 0.005);
  // Corriente de la carga a V0Set: 2 µA (en SI, 2e-6 A).
  EXPECT_NEAR(f.param(0, "IMon"), 2e-6, 1.1e-9);
}

TEST(FakeHV, RampDownAt50VPerSecond) {
  auto sc = quietScenario({{{"t_s", 1.0}, {"type", "set"}, {"ch", 1}, {"param", "Pw"}, {"value", 0}}});
  sc["slots"][0]["set"]["Pw"] = nlohmann::json::array({1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0});
  sc["slots"][0]["set"]["PDwn"] = 1;  // Ramp
  Fixture f(sc);
  EXPECT_NEAR(f.be.trueVoltage(3, 1), 65.75, 1e-9);
  f.at(1s + 500ms);
  EXPECT_TRUE(f.bits.has(f.status(1), "ramp_down"));
  EXPECT_FALSE(f.bits.has(f.status(1), "on"));
  EXPECT_NEAR(f.be.trueVoltage(3, 1), 40.75, 1e-6);
  f.at(1s + 1315ms);
  EXPECT_NEAR(f.be.trueVoltage(3, 1), 0.0, 1e-9);
  EXPECT_EQ(f.status(1), 0u);
}

TEST(FakeHV, KillDropsVoltageImmediately) {
  auto sc = quietScenario({{{"t_s", 1.0}, {"type", "set"}, {"ch", 1}, {"param", "Pw"}, {"value", 0}}});
  sc["slots"][0]["set"]["Pw"] = 1;  // PDwn = 0 (Kill) en el guion por defecto
  Fixture f(sc);
  f.at(1s + 10ms);
  EXPECT_NEAR(f.be.trueVoltage(3, 1), 0.0, 1e-9);
  EXPECT_EQ(f.status(1), 0u);
}

TEST(FakeHV, OvercurrentLongerThanTripTripsAndStaysOff) {
  // Canal 2 encendido; a t = 100 s una excursión de +20 µA (total 22 > I0Set = 15) durante 30 s.
  auto sc = quietScenario({{{"t_s", 100.0}, {"type", "current"}, {"ch", 2}, {"extra_uA", 20.0}, {"duration_s", 30.0}}});
  sc["slots"][0]["set"]["Pw"] = nlohmann::json::array({1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0});
  Fixture f(sc);
  f.at(99s);
  EXPECT_EQ(f.status(2), 0x0001u);
  f.at(100s + 100ms);
  const auto s1 = f.status(2);
  EXPECT_TRUE(f.bits.has(s1, "on"));
  EXPECT_TRUE(f.bits.has(s1, "overcurrent"));
  EXPECT_NEAR(f.param(2, "IMon"), 15e-6, 1e-12);  // limitado a I0Set
  f.at(104s + 900ms);  // 4.9 s en sobrecorriente: todavía encendido
  EXPECT_EQ(f.pw(2), 1u);
  f.at(105s + 100ms);  // pasó Trip = 5 s
  EXPECT_EQ(f.pw(2), 0u);
  const auto s2 = f.status(2);
  EXPECT_FALSE(f.bits.has(s2, "on"));
  EXPECT_TRUE(f.bits.has(s2, "internal_trip"));
  EXPECT_FALSE(f.bits.has(s2, "overcurrent"));
  // No se vuelve a encender solo, ni siquiera cuando termina la excursión.
  f.at(200s);
  EXPECT_EQ(f.pw(2), 0u);
  EXPECT_TRUE(f.bits.has(f.status(2), "internal_trip"));
  EXPECT_NEAR(f.be.trueVoltage(3, 2), 0.0, 1e-9);
  // Los demás canales siguen encendidos.
  EXPECT_EQ(f.pw(3), 1u);
}

TEST(FakeHV, ShortOvercurrentDoesNotTrip) {
  auto sc = quietScenario({{{"t_s", 10.0}, {"type", "current"}, {"ch", 0}, {"extra_uA", 20.0}, {"duration_s", 3.0}}});
  sc["slots"][0]["set"]["Pw"] = 1;
  Fixture f(sc);
  f.at(12s);
  EXPECT_TRUE(f.bits.has(f.status(0), "overcurrent"));
  f.at(20s);
  EXPECT_EQ(f.pw(0), 1u);
  EXPECT_EQ(f.status(0), 0x0001u);
}

TEST(FakeHV, ManualPwOnClearsTripAndRampsUp) {
  auto sc = quietScenario({{{"t_s", 1.0}, {"type", "current"}, {"ch", 4}, {"extra_uA", 50.0}, {"duration_s", 10.0}},
                           {{"t_s", 60.0}, {"type", "set"}, {"ch", 4}, {"param", "Pw"}, {"value", 1}}});
  sc["slots"][0]["set"]["Pw"] = 1;
  Fixture f(sc);
  f.at(30s);
  EXPECT_TRUE(f.bits.has(f.status(4), "internal_trip"));
  f.at(61s);
  const auto s = f.status(4);
  EXPECT_FALSE(f.bits.has(s, "internal_trip"));
  EXPECT_TRUE(f.bits.has(s, "ramp_up"));
  EXPECT_NEAR(f.be.trueVoltage(3, 4), 1.0, 1e-6);
}

TEST(FakeHV, ConnectionLossReturnsErrorsAndRecoveryNeedsConnect) {
  auto sc = quietScenario({{{"t_s", 100.0}, {"type", "connection_loss"}, {"duration_s", 30.0}}});
  ManualClock clock;
  FakeHVBackend be(sc, clock);
  ASSERT_TRUE(be.connect().ok());
  clock.sleep_until(daq::SteadyTime{} + 50s);
  EXPECT_TRUE(be.readChannels(3, {0}, {"VMon"}).ok());
  clock.sleep_until(daq::SteadyTime{} + 101s);
  auto r = be.readChannels(3, {0}, {"VMon"});
  ASSERT_FALSE(r.ok());
  EXPECT_FALSE(r.error().message.empty());
  EXPECT_FALSE(be.isConnected());
  EXPECT_FALSE(be.connect().ok());  // el mainframe sigue inalcanzable
  EXPECT_FALSE(be.crateMap().ok());
  clock.sleep_until(daq::SteadyTime{} + 131s);
  EXPECT_FALSE(be.readChannels(3, {0}, {"VMon"}).ok());  // hay que reconectar
  ASSERT_TRUE(be.connect().ok());
  EXPECT_TRUE(be.isConnected());
  EXPECT_TRUE(be.readChannels(3, {0}, {"VMon"}).ok());
}

TEST(FakeHV, ErrorsAreExplicitNotZeros) {
  ManualClock clock;
  FakeHVBackend be(FakeHVBackend::defaultScenario(), clock);
  EXPECT_FALSE(be.readChannels(3, {0}, {"VMon"}).ok());  // no conectado
  ASSERT_TRUE(be.connect().ok());
  EXPECT_FALSE(be.readChannels(3, {0}, {"NoExiste"}).ok());
  EXPECT_FALSE(be.readChannels(3, {12}, {"VMon"}).ok());
  EXPECT_FALSE(be.readChannels(5, {0}, {"VMon"}).ok());
  EXPECT_FALSE(be.channelParamNames(4, 0).ok());
}

TEST(FakeHV, ReadingCarriesUtcAndDuration) {
  ManualClock clock;
  FakeHVBackend be(FakeHVBackend::defaultScenario(), clock);
  ASSERT_TRUE(be.connect().ok());
  clock.advance(5s);
  auto r = be.readChannels(3, {0, 1}, {"VMon", "IMon"});
  ASSERT_TRUE(r.ok());
  EXPECT_EQ(r.value().tReqUtcNs, 1792281600LL * 1000000000LL + 5000000000LL);
  EXPECT_GE(r.value().tRespUtcNs, r.value().tReqUtcNs);
  EXPECT_GE(r.value().durationNs, 0);
}

TEST(FakeHV, InvalidScenarioThrows) {
  ManualClock clock;
  EXPECT_THROW(FakeHVBackend(nlohmann::json::object(), clock), std::invalid_argument);
  auto sc = FakeHVBackend::defaultScenario();
  sc["events"] = {{{"t_s", 1.0}, {"type", "boom"}}};
  EXPECT_THROW(FakeHVBackend(sc, clock), std::invalid_argument);
  sc["events"] = {{{"t_s", 1.0}, {"type", "set"}, {"ch", 99}, {"param", "Pw"}, {"value", 1}}};
  EXPECT_THROW(FakeHVBackend(sc, clock), std::invalid_argument);
}

TEST(FakeHV, FromFileWithBaseDefaultUsesDefaultSlotsAndFileEvents) {
  const std::string path = ::testing::TempDir() + "fake_base_default.json";
  {
    std::ofstream o(path);
    o << R"({"base": "default", "events": [{"t_s": 2, "type": "set", "ch": 5, "param": "Pw", "value": 1}]})";
  }
  ManualClock clock;
  auto be = FakeHVBackend::fromFile(path, clock);
  ASSERT_TRUE(be.ok()) << be.error().message;
  ASSERT_TRUE(be.value()->connect().ok());
  clock.sleep_until(daq::SteadyTime{} + 3s);
  auto r = be.value()->readChannels(3, {5}, {"Pw", "V0Set"});
  ASSERT_TRUE(r.ok());
  EXPECT_EQ(r.value().get("Pw", 5)->raw, 1u);
  EXPECT_DOUBLE_EQ(r.value().get("V0Set", 5)->boardValue, 65.70);
  EXPECT_FALSE(FakeHVBackend::fromFile(::testing::TempDir() + "no_existe.json", clock).ok());
  std::remove(path.c_str());
}
