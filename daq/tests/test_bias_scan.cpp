// Tests del barrido de tensión de hvmon (DAQ2-14) con FakeHVBackend y ManualClock: límites de cada
// escritura de V0Set, verificación de la rampa, abortos y restauración del nominal (pedida, tras un
// trip y al arrancar hvmon con una sesión persistida).
#include <gtest/gtest.h>

#include <chrono>

#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hvmon/bias_scan.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/poller.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::hvmon;

namespace {

const std::map<int, double> kNominal{{0, 65.75}, {1, 65.75}, {2, 65.90}, {3, 65.85}, {4, 65.75}, {5, 65.70}, {6, 66.15}, {7, 65.90}};

nlohmann::json scenario(nlohmann::json events = nlohmann::json::array()) {
  auto sc = hv::FakeHVBackend::defaultScenario();
  auto& s = sc["slots"][0];
  s["imon_noise_uA"] = 0.0;
  s["vmon_noise_V"] = 0.0;
  s["set"]["Pw"] = nlohmann::json{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
  s["set"]["RDWn"] = nlohmann::json{1, 1, 1, 1, 1, 1, 1, 1, 50, 50, 50, 50};
  sc["events"] = std::move(events);
  return sc;
}

struct Rig {
  ManualClock clock;
  HvmonConfig cfg;
  hv::FakeHVBackend backend;
  BoundedQueue<WriterItem> writer{1000000};
  BoundedQueue<PollRequest> requests{16};
  StateHub hub;
  MemoryFileSystem* fs;
  Poller poller;
  Rig(const nlohmann::json& sc, MemoryFileSystem& f, bool allow = true)
      : cfg(makeCfg(allow)), backend(sc, clock), fs(&f), poller(cfg, backend, clock, writer, hub, requests) {
    poller.enableBiasPersistence(f, "/state/hvmon_bias_state.json");
  }
  static HvmonConfig makeCfg(bool allow) {
    auto c = HvmonConfig::fromJson(nlohmann::json::object()).value();
    c.allowCommands = allow;
    return c;
  }
  void run(std::chrono::milliseconds d) {
    const auto end = clock.now_steady() + d;
    while (clock.now_steady() < end) poller.step();
    while (writer.try_pop()) {
    }
  }
  nlohmann::json send(PollRequest r, std::chrono::milliseconds maxWait = 600s) {
    r.requester = "gtest";
    if (r.reason.empty()) r.reason = "test";
    auto f = r.reply->get_future();
    EXPECT_TRUE(requests.try_push(std::move(r), DropPolicy::DropNewest));
    const auto end = clock.now_steady() + maxWait;
    while (f.wait_for(0s) != std::future_status::ready && clock.now_steady() < end) poller.step();
    while (writer.try_pop()) {
    }
    return f.wait_for(0s) == std::future_status::ready ? f.get() : nlohmann::json{{"ok", false}, {"result", "test_timeout"}};
  }
  nlohmann::json begin(std::map<int, double> exp = kNominal) {
    PollRequest r;
    r.kind = PollRequest::Kind::BiasBegin;
    r.expected = std::move(exp);
    return send(std::move(r));
  }
  nlohmann::json v0set(int ch, double v) {
    PollRequest r;
    r.kind = PollRequest::Kind::V0Set;
    r.ch = ch;
    r.value = v;
    return send(std::move(r));
  }
  nlohmann::json restore() {
    PollRequest r;
    r.kind = PollRequest::Kind::BiasRestore;
    return send(std::move(r), 3600s);
  }
  double v0(int ch) { return backend.readChannels(3, {ch}, {"V0Set"}).value().get("V0Set", ch)->si; }
};

}  // namespace

TEST(BiasSession, LimitsAndRestorePlan) {
  BiasSession s;
  std::map<int, double> vm;
  for (const auto& [ch, v] : kNominal) vm[ch] = v - 0.3;
  EXPECT_FALSE(s.checkStep(0, 65.5, 65.75));  // sin sesión
  auto wrong = kNominal;
  wrong[3] = 65.80;
  EXPECT_FALSE(s.begin(wrong, kNominal, vm, "t", 1));  // el nominal leído no es el esperado
  ASSERT_TRUE(s.begin(kNominal, kNominal, vm, "t", 1));
  EXPECT_TRUE(s.checkStep(0, 65.25, 65.75));
  EXPECT_FALSE(s.checkStep(0, 65.80, 65.75));   // sobre el nominal
  EXPECT_FALSE(s.checkStep(0, 65.15, 65.75));   // escalón de 0.6 V
  EXPECT_FALSE(s.checkStep(0, 63.70, 63.75));   // más de 2.0 V bajo el nominal
  EXPECT_TRUE(s.checkStep(0, 63.75, 64.25));    // justo 2.0 V bajo el nominal
  EXPECT_FALSE(s.checkStep(8, 0, 0));           // canal fuera de la sesión
  // escalón 65.75 -> 65.25 V con VMon 65.45 antes (desvío −0.30 V)
  EXPECT_TRUE(s.reached(0, 65.25, 65.25, 64.95, 65.75, 65.45));   // desvío −0.30 V como al nominal, movió −0.50 V
  EXPECT_TRUE(s.reached(0, 65.25, 65.25, 65.04, 65.75, 65.45));   // desvío −0.21 V: dentro de 0.10 V
  EXPECT_FALSE(s.reached(0, 65.25, 65.25, 65.25, 65.75, 65.45));  // VMon = V0Set: 0.30 V fuera del desvío
  EXPECT_FALSE(s.reached(0, 65.25, 65.24, 64.95, 65.75, 65.45));  // V0Set releído distinto del pedido
  auto cur = kNominal;
  cur[0] = 64.25;
  const auto plan = s.restorePlan(cur);
  ASSERT_EQ(plan.size(), 3u);
  EXPECT_DOUBLE_EQ(plan[0].second, 64.75);
  EXPECT_DOUBLE_EQ(plan[2].second, 65.75);
  BiasSession t;
  ASSERT_TRUE(t.restore(s.toJson(2)));
  EXPECT_TRUE(t.active());
  EXPECT_DOUBLE_EQ(t.nominal(6), 66.15);
}

TEST(BiasScanPoller, RefusedWithoutAllowCommandsOrSession) {
  MemoryFileSystem fs;
  Rig r(scenario(), fs, false);
  r.run(5s);
  EXPECT_FALSE(r.begin().value("ok", true));
  Rig q(scenario(), fs, true);
  q.run(5s);
  EXPECT_FALSE(q.v0set(0, 65.25).value("ok", true));  // sin sesión
  EXPECT_EQ(q.backend.v0setWrites(), 0u);
}

TEST(BiasScanPoller, StepsAreLimitedRampIsVerifiedAndRestoreReturnsToNominal) {
  MemoryFileSystem fs;
  Rig r(scenario(), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  EXPECT_TRUE(fs.exists("/state/hvmon_bias_state.json"));
  EXPECT_FALSE(r.v0set(0, 65.80).value("ok", true));  // sobre el nominal
  EXPECT_FALSE(r.v0set(0, 65.15).value("ok", true));  // escalón > 0.5 V
  EXPECT_EQ(r.backend.v0setWrites(), 0u);
  for (double v : {65.25, 64.75, 64.25, 63.75}) {
    const auto res = r.v0set(0, v);
    ASSERT_TRUE(res.value("ok", false)) << res.dump();
    EXPECT_EQ(res.value("result", ""), "verified");
    EXPECT_NEAR(r.v0(0), v, 1e-6);
  }
  EXPECT_FALSE(r.v0set(0, 63.25).value("ok", true));  // más de 2.0 V bajo el nominal
  EXPECT_EQ(r.backend.v0setWrites(), 4u);
  const auto res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_NEAR(r.v0(0), 65.75, 1e-6);
  EXPECT_EQ(r.backend.v0setWrites(), 8u);  // cuatro escalones de vuelta
  EXPECT_FALSE(r.poller.bias().active());
  EXPECT_FALSE(nlohmann::json::parse(fs.content("/state/hvmon_bias_state.json")).value("active", true));
}

TEST(BiasScanPoller, BeginRefusedWhenNominalDiffers) {
  MemoryFileSystem fs;
  Rig r(scenario(), fs);
  r.run(5s);
  auto exp = kNominal;
  exp[2] = 65.80;
  EXPECT_FALSE(r.begin(exp).value("ok", true));
  EXPECT_FALSE(r.poller.bias().active());
}

TEST(BiasScanPoller, TripDuringRampAbortsAndRestoreSwitchesTheChannelBackOn) {
  MemoryFileSystem fs;
  // de 20 a 40 s CH01 recibe 30 µA extra: sobrecorriente y trip tras 5 s (Trip = 5 s)
  Rig r(scenario(nlohmann::json::array({{{"t_s", 20.0}, {"type", "current"}, {"slot", 3}, {"ch", 1}, {"extra_uA", 30.0}, {"duration_s", 20.0}}})), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  ASSERT_TRUE(r.v0set(1, 65.25).value("ok", false));
  r.run(30s);  // CH01 se dispara
  EXPECT_FALSE(r.v0set(1, 64.75).value("ok", true));  // el canal ya no está encendido y estable
  r.run(20s);   // vuelve la corriente normal
  const auto res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_NEAR(r.v0(1), 65.75, 1e-6);
  EXPECT_EQ(r.backend.readChannels(3, {1}, {"Pw"}).value().get("Pw", 1)->raw, 1u);
}

TEST(BiasScanPoller, HvmonStartingWithActiveSessionRestoresNominalByItself) {
  MemoryFileSystem fs;
  {
    Rig r(scenario(), fs);
    r.run(5s);
    ASSERT_TRUE(r.begin().value("ok", false));
    for (int ch = 0; ch < 8; ++ch) ASSERT_TRUE(r.v0set(ch, kNominal.at(ch) - 0.5).value("ok", false));
  }
  // "Reinicio" de hvmon: la placa queda con V0Set bajo el nominal y el estado persistido activo.
  auto sc = scenario();
  for (int ch = 0; ch < 8; ++ch) sc["slots"][0]["set"]["V0Set"][ch] = kNominal.at(ch) - 0.5;
  Rig r2(sc, fs);
  EXPECT_TRUE(r2.poller.bias().active());
  r2.run(300s);
  for (int ch = 0; ch < 8; ++ch) EXPECT_NEAR(r2.v0(ch), kNominal.at(ch), 1e-6);
  EXPECT_FALSE(r2.poller.bias().active());
}

TEST(BiasScanPoller, OvercurrentDuringRampFailsVerification) {
  MemoryFileSystem fs;
  auto sc = scenario(nlohmann::json::array({{{"t_s", 8.0}, {"type", "current"}, {"slot", 3}, {"ch", 2}, {"extra_uA", 30.0}, {"duration_s", 3.0}}}));
  sc["slots"][0]["set"]["RDWn"][2] = 0.05;  // rampa lenta: 0.5 V en 10 s
  Rig r(sc, fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  const auto res = r.v0set(2, 65.40);
  EXPECT_FALSE(res.value("ok", true)) << res.dump();
  EXPECT_EQ(res.value("result", ""), "verify_failed");
  EXPECT_NE(res.value("detail", "").find("Status inesperado"), std::string::npos) << res.dump();
}

TEST(BiasSession, OffsetReferenceFollowsLastVerifiedStep) {
  BiasSession s;
  std::map<int, double> vm;
  for (const auto& [ch, v] : kNominal) vm[ch] = v - 0.30;
  ASSERT_TRUE(s.begin(kNominal, kNominal, vm, "t", 1));
  // el desvío de la placa se corrió 0.04 V en el escalón anterior: se verifica contra ese
  s.updateOffset(0, -0.34);
  EXPECT_TRUE(s.reached(0, 65.25, 65.25, 64.90, 65.75, 65.41));
  EXPECT_FALSE(s.reached(0, 65.25, 65.25, 65.03, 65.75, 65.41));  // 0.12 V fuera del desvío vigente
}

TEST(BiasSession, VMonMustMoveWithTheStep) {
  BiasSession s;
  std::map<int, double> vm;
  for (const auto& [ch, v] : kNominal) vm[ch] = v - 0.30;
  ASSERT_TRUE(s.begin(kNominal, kNominal, vm, "t", 1));
  s.updateOffset(0, -0.30);
  // V0Set 65.25 -> 64.75 V (−0.50 V); VMon antes 64.95
  EXPECT_TRUE(s.reached(0, 64.75, 64.75, 64.45, 65.25, 64.95));   // movió −0.50 V
  EXPECT_TRUE(s.reached(0, 64.75, 64.75, 64.53, 65.25, 64.95));   // movió −0.42 V (≥ la mitad), desvío −0.22 V
  EXPECT_FALSE(s.reached(0, 64.75, 64.75, 64.73, 65.25, 64.83));  // movió −0.10 V: menos de la mitad
  // subida (restauración) 64.75 -> 65.25 V
  EXPECT_TRUE(s.reached(0, 65.25, 65.25, 64.95, 64.75, 64.45));
  EXPECT_FALSE(s.reached(0, 65.25, 65.25, 64.95, 64.75, 65.05));  // VMon bajó en un escalón de subida
}
