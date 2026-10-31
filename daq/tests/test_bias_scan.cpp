// Tests del barrido de tensión de hvmon (DAQ2-14) con FakeHVBackend y ManualClock: límites de cada
// escritura de V0Set, verificación de la rampa, abortos y restauración del nominal (pedida, tras un
// trip y al arrancar hvmon con una sesión persistida). Seguimiento 2026-10-31 (7): RUp dentro de la sesión (fijar, leer,
// restaurar el previo con bias-restore, abortar a mitad, fuera de rango, reinicio de hvmon con la sesión activa).
#include <gtest/gtest.h>

#include <chrono>
#include <variant>

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
  std::vector<CommandRecord> cmds;  // registros de comando que irían a hv_changes
  void drain() {
    while (auto it = writer.try_pop())
      if (const auto* k = std::get_if<CommandItem>(&*it)) cmds.push_back(k->rec);
  }
  void run(std::chrono::milliseconds d) {
    const auto end = clock.now_steady() + d;
    while (clock.now_steady() < end) poller.step();
    drain();
  }
  nlohmann::json send(PollRequest r, std::chrono::milliseconds maxWait = 600s) {
    r.requester = "gtest";
    if (r.reason.empty()) r.reason = "test";
    auto f = r.reply->get_future();
    EXPECT_TRUE(requests.try_push(std::move(r), DropPolicy::DropNewest));
    const auto end = clock.now_steady() + maxWait;
    while (f.wait_for(0s) != std::future_status::ready && clock.now_steady() < end) poller.step();
    drain();
    return f.wait_for(0s) == std::future_status::ready ? f.get() : nlohmann::json{{"ok", false}, {"result", "test_timeout"}};
  }
  nlohmann::json begin(std::map<int, double> exp = kNominal, double maxAboveV = 0.0, double maxStepV = -1.0) {
    PollRequest r;
    r.kind = PollRequest::Kind::BiasBegin;
    r.expected = std::move(exp);
    r.maxAboveV = maxAboveV;
    r.maxStepV = maxStepV;
    return send(std::move(r));
  }
  nlohmann::json v0set(int ch, double v, std::chrono::milliseconds maxWait = 900s) {
    PollRequest r;
    r.kind = PollRequest::Kind::V0Set;
    r.ch = ch;
    r.value = v;
    return send(std::move(r), maxWait);
  }
  nlohmann::json restore() {
    PollRequest r;
    r.kind = PollRequest::Kind::BiasRestore;
    return send(std::move(r), 3600s);
  }
  double v0(int ch) { return backend.readChannels(3, {ch}, {"V0Set"}).value().get("V0Set", ch)->si; }
  nlohmann::json rup(int ch, double vps) {
    PollRequest r;
    r.kind = PollRequest::Kind::RUp;
    r.ch = ch;
    r.value = vps;
    return send(std::move(r), 60s);
  }
  double rupRead(int ch) { return backend.readChannels(3, {ch}, {"RUp"}).value().get("RUp", ch)->si; }
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
  EXPECT_EQ(res.value("result", ""), "hazard");
  EXPECT_NE(res.value("detail", "").find("riesgo"), std::string::npos) << res.dump();
}

// DAQ2-14 (2026-10-23): bajada con la descarga del PCB (τ = 115 s): Status on + ramp_down durante
// minutos es un estado esperado; se espera hasta que la rampa termina (timeout 6τ) y se verifica.
TEST(BiasScanPoller, SlowDischargeWithTau115IsWaitedNotAborted) {
  MemoryFileSystem fs;
  auto sc = scenario();
  sc["slots"][0]["discharge_tau_s"] = 115.0;
  Rig r(sc, fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  const auto t0 = r.clock.now_steady();
  const auto st0 = r.backend.readChannels(3, {0}, {"Status"});
  const auto res = r.v0set(0, 65.25);
  const double waited = std::chrono::duration<double>(r.clock.now_steady() - t0).count();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_EQ(res.value("result", ""), "verified");
  // V − V0Set baja de 0.5 a 0.01 V en τ·ln(50) ≈ 450 s: el bit ramp_down estuvo puesto todo ese tiempo
  EXPECT_GT(waited, 400.0);
  EXPECT_LT(waited, 6 * 115.0);
  (void)st0;
}

TEST(BiasScanPoller, RampLongerThanSixTauTimesOutAsRecoverable) {
  MemoryFileSystem fs;
  auto sc = scenario();
  sc["slots"][0]["discharge_tau_s"] = 400.0;  // ln(50)·400 s ≫ 6·115 s
  Rig r(sc, fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  const auto res = r.v0set(0, 65.25);
  EXPECT_EQ(res.value("result", ""), "verify_failed") << res.dump();  // fallo recuperable, no riesgo
}

TEST(BiasScanPoller, VMonMovingAgainstTheStepIsAHazard) {
  MemoryFileSystem fs;
  // durante la bajada lenta de CH04, alguien (GECO) sube su V0Set a 66.30 V: VMon sube en contra del escalón
  auto sc = scenario(nlohmann::json::array({{{"t_s", 30.0}, {"type", "set"}, {"slot", 3}, {"ch", 4}, {"param", "V0Set"}, {"value", 66.30}}}));
  sc["slots"][0]["discharge_tau_s"] = 115.0;
  Rig r(sc, fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  const auto res = r.v0set(4, 65.25);
  EXPECT_EQ(res.value("result", ""), "hazard") << res.dump();
  EXPECT_NE(res.value("detail", "").find("en contra del escalón"), std::string::npos) << res.dump();
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

// DAQ2-17 (2026-10-23): barrido hacia arriba autorizado por René: hasta nominal + 1.0 V (tope en el código), escalones de
// 0.25 V; por omisión una sesión sigue sin poder pasar el nominal; la restauración baja desde arriba en escalones ≤ 0.25 V.
TEST(BiasSession, UpwardBoundIsPerSessionWithHardCapAndSmallSteps) {
  std::map<int, double> vm;
  for (const auto& [ch, v] : kNominal) vm[ch] = v - 0.3;
  BiasSession d;
  ASSERT_TRUE(d.begin(kNominal, kNominal, vm, "t", 1));  // sin max_above: como antes
  EXPECT_FALSE(d.checkStep(0, 66.00, 65.75));
  BiasSession bad;
  EXPECT_FALSE(bad.begin(kNominal, kNominal, vm, "t", 1, 1.5));        // más que el tope de 1.0 V
  EXPECT_FALSE(bad.begin(kNominal, kNominal, vm, "t", 1, 1.0, 0.75));  // escalón de sesión > 0.5 V
  BiasSession s;
  ASSERT_TRUE(s.begin(kNominal, kNominal, vm, "t", 1, 1.0, 0.25));
  EXPECT_TRUE(s.checkStep(0, 66.00, 65.75));   // +0.25 V
  EXPECT_FALSE(s.checkStep(0, 66.10, 65.75));  // escalón de 0.35 V
  EXPECT_TRUE(s.checkStep(0, 66.75, 66.50));   // justo nominal + 1.0 V
  EXPECT_FALSE(s.checkStep(0, 66.80, 66.75));  // sobre nominal + 1.0 V
  EXPECT_TRUE(s.checkStep(0, 65.50, 65.75));   // bajo el nominal, escalón de la sesión
  EXPECT_FALSE(s.checkStep(0, 65.25, 65.75));  // 0.5 V > escalón de la sesión (0.25 V)
  BiasSession w;  // escalón de sesión de 0.5 V: bajo el nominal vale, sobre el nominal sigue el tope de 0.25 V
  ASSERT_TRUE(w.begin(kNominal, kNominal, vm, "t", 1, 1.0, 0.5));
  EXPECT_TRUE(w.checkStep(0, 65.25, 65.75));
  EXPECT_FALSE(w.checkStep(0, 66.25, 65.75));
  EXPECT_FALSE(w.checkStep(0, 65.75, 66.25));  // también al bajar desde arriba
  auto cur = kNominal;
  cur[0] = 66.75;
  const auto plan = s.restorePlan(cur);
  ASSERT_EQ(plan.size(), 4u);
  EXPECT_DOUBLE_EQ(plan[0].second, 66.50);
  EXPECT_DOUBLE_EQ(plan[3].second, 65.75);
  const auto planW = w.restorePlan(cur);  // con escalón de sesión 0.5 V, desde arriba igual 0.25 V
  EXPECT_EQ(planW.size(), 4u);
  BiasSession t;  // la sesión persistida conserva sus límites
  ASSERT_TRUE(t.restore(s.toJson(2)));
  EXPECT_DOUBLE_EQ(t.limits().maxAboveV, 1.0);
  EXPECT_DOUBLE_EQ(t.limits().maxStepV, 0.25);
  EXPECT_TRUE(t.checkStep(0, 66.50, 66.75));
}

TEST(BiasScanPoller, UpwardScanToPlusOneVoltAndRestoreDownToNominal) {
  MemoryFileSystem fs;
  Rig r(scenario(), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin(kNominal, 1.0, 0.25).value("ok", false));
  EXPECT_FALSE(r.v0set(0, 66.25).value("ok", true));  // escalón de 0.5 V sobre el nominal
  EXPECT_EQ(r.backend.v0setWrites(), 0u);
  for (double v : {66.00, 66.25, 66.50, 66.75}) {
    const auto res = r.v0set(0, v);
    ASSERT_TRUE(res.value("ok", false)) << res.dump();
    EXPECT_EQ(res.value("result", ""), "verified");
    EXPECT_NEAR(r.v0(0), v, 1e-6);
  }
  EXPECT_FALSE(r.v0set(0, 67.00).value("ok", true));  // sobre nominal + 1.0 V
  EXPECT_EQ(r.backend.v0setWrites(), 4u);
  const auto res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_NEAR(r.v0(0), 65.75, 1e-6);
  EXPECT_EQ(r.backend.v0setWrites(), 8u);  // cuatro escalones de 0.25 V de vuelta
  EXPECT_FALSE(r.poller.bias().active());
}

TEST(BiasScanPoller, DefaultSessionStillRefusesAboveNominal) {
  MemoryFileSystem fs;
  Rig r(scenario(), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  EXPECT_FALSE(r.v0set(0, 66.00).value("ok", true));
  EXPECT_EQ(r.backend.v0setWrites(), 0u);
}

TEST(BiasScanPoller, TripAboveNominalRestoresNominalAndSwitchesOn) {
  MemoryFileSystem fs;
  Rig r(scenario(nlohmann::json::array({{{"t_s", 60.0}, {"type", "current"}, {"slot", 3}, {"ch", 1}, {"extra_uA", 30.0}, {"duration_s", 20.0}}})), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin(kNominal, 1.0, 0.25).value("ok", false));
  for (double v : {66.00, 66.25}) ASSERT_TRUE(r.v0set(1, v).value("ok", false));
  r.run(60s);  // CH01 se dispara a nominal + 0.5 V
  EXPECT_EQ(r.backend.readChannels(3, {1}, {"Pw"}).value().get("Pw", 1)->raw, 0u);
  r.run(30s);
  const auto res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_NEAR(r.v0(1), 65.75, 1e-6);
  EXPECT_EQ(r.backend.readChannels(3, {1}, {"Pw"}).value().get("Pw", 1)->raw, 1u);
  EXPECT_FALSE(r.poller.bias().active());
}

TEST(BiasScanPoller, HvmonStartingWithActiveUpwardSessionRestoresDownByItself) {
  MemoryFileSystem fs;
  {
    Rig r(scenario(), fs);
    r.run(5s);
    ASSERT_TRUE(r.begin(kNominal, 1.0, 0.25).value("ok", false));
    for (int ch = 0; ch < 8; ++ch)
      for (double d : {0.25, 0.50}) ASSERT_TRUE(r.v0set(ch, kNominal.at(ch) + d).value("ok", false));
  }
  auto sc = scenario();
  for (int ch = 0; ch < 8; ++ch) sc["slots"][0]["set"]["V0Set"][ch] = kNominal.at(ch) + 0.5;
  Rig r2(sc, fs);
  EXPECT_TRUE(r2.poller.bias().active());
  r2.run(600s);
  for (int ch = 0; ch < 8; ++ch) EXPECT_NEAR(r2.v0(ch), kNominal.at(ch), 1e-6);
  EXPECT_FALSE(r2.poller.bias().active());
  EXPECT_EQ(r2.backend.v0setWrites(), 16u);  // dos escalones de 0.25 V por canal
}

// DAQ2-17 (René): el desvío de referencia es el del último escalón verificado, pero un desvío que se aleja más de 0.15 V
// del desvío al nominal es un riesgo (aborto); la restauración no tiene ese límite.
TEST(BiasSession, OffsetDriftFromNominalIsBounded) {
  BiasSession s;
  std::map<int, double> vm;
  for (const auto& [ch, v] : kNominal) vm[ch] = v - 0.30;
  ASSERT_TRUE(s.begin(kNominal, kNominal, vm, "t", 1, 1.0, 0.25));
  EXPECT_TRUE(s.offsetDriftOk(0, 66.00, 66.00 - 0.30 + 0.15));   // justo 0.15 V
  EXPECT_FALSE(s.offsetDriftOk(0, 66.00, 66.00 - 0.30 + 0.16));
  EXPECT_FALSE(s.offsetDriftOk(0, 66.00, 66.00 - 0.30 - 0.16));
  s.updateOffset(0, -0.20);  // el de referencia cambia; el del nominal no
  EXPECT_NEAR(s.nominalOffset(0), -0.30, 1e-9);
  BiasSession t;
  ASSERT_TRUE(t.restore(s.toJson(2)));
  EXPECT_NEAR(t.nominalOffset(0), -0.30, 1e-9);
  EXPECT_NEAR(t.offset(0), -0.20, 1e-9);
}

TEST(BiasScanPoller, OffsetDriftOverPointFifteenFromNominalIsAHazardAndRestoreStillWorks) {
  MemoryFileSystem fs;
  auto sc = scenario(nlohmann::json::array({{{"t_s", 100.0}, {"type", "vmon_offset"}, {"slot", 3}, {"ch", 0}, {"shift_V", 0.08}},
                                            {{"t_s", 200.0}, {"type", "vmon_offset"}, {"slot", 3}, {"ch", 0}, {"shift_V", 0.08}}}));
  Rig r(sc, fs);
  r.run(5s);
  ASSERT_TRUE(r.begin(kNominal, 1.0, 0.25).value("ok", false));
  ASSERT_EQ(r.v0set(0, 66.00).value("result", ""), "verified");
  r.run(110s);  // corrimiento de +0.08 V: 0.08 V del nominal, dentro de 0.10 V del escalón anterior
  const auto a = r.v0set(0, 66.25);
  ASSERT_EQ(a.value("result", ""), "verified") << a.dump();
  r.run(100s);  // otro +0.08 V: 0.16 V del nominal
  const auto b = r.v0set(0, 66.50);
  EXPECT_EQ(b.value("result", ""), "hazard") << b.dump();
  EXPECT_NE(b.value("detail", "").find("desvío al nominal"), std::string::npos) << b.dump();
  const auto res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_NEAR(r.v0(0), 65.75, 1e-6);
  EXPECT_FALSE(r.poller.bias().active());
}

// XPLUS-1: tras un corte de la fuente, el V0Set de un canal apagado vuelve al nominal que hvmon registró en la última sesión
// (estado persistido, sesión cerrada). Solo con el canal apagado, sin sesión activa y con allow_commands; el valor no viene en
// la petición.
namespace {
nlohmann::json nominalReq(Rig& r, int ch) {
  PollRequest q;
  q.kind = PollRequest::Kind::NominalV0Set;
  q.ch = ch;
  return r.send(std::move(q));
}
void recordNominal(MemoryFileSystem& fs) {
  Rig r(scenario(), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  ASSERT_TRUE(r.restore().value("ok", false));  // cierra la sesión: queda el nominal persistido con active = false
}
}  // namespace

TEST(BiasScanPoller, NominalV0SetRestoresALostValueOnlyOnAnOffChannel) {
  MemoryFileSystem fs;
  recordNominal(fs);
  auto sc = scenario();
  sc["slots"][0]["set"]["Pw"] = nlohmann::json{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};  // tras el corte: apagados, salvo CH07
  sc["slots"][0]["set"]["V0Set"][3] = 0.0;                                          // V0Set perdido en CH03
  sc["slots"][0]["set"]["V0Set"][7] = 10.0;                                         // y en CH07, que está encendido
  Rig r(sc, fs);
  EXPECT_FALSE(r.poller.bias().active());
  r.run(5s);
  const auto ok = nominalReq(r, 3);
  EXPECT_TRUE(ok.value("ok", false)) << ok.dump();
  EXPECT_EQ(ok.value("result", ""), "verified");
  EXPECT_NEAR(r.v0(3), kNominal.at(3), 1e-6);
  EXPECT_EQ(r.backend.v0setWrites(), 1u);
  const auto same = nominalReq(r, 0);  // ya en el nominal
  EXPECT_TRUE(same.value("ok", false));
  EXPECT_EQ(same.value("result", ""), "noop");
  const auto on = nominalReq(r, 7);  // encendido: rechazado, sin escritura
  EXPECT_FALSE(on.value("ok", true));
  EXPECT_EQ(on.value("result", ""), "rejected");
  const auto outside = nominalReq(r, 9);  // canal fuera de la sesión
  EXPECT_EQ(outside.value("result", ""), "rejected");
  EXPECT_EQ(r.backend.v0setWrites(), 1u);
  EXPECT_NEAR(r.v0(7), 10.0, 1e-6);
}

TEST(BiasScanPoller, NominalV0SetRefusedWithoutRecordWithoutAllowOrDuringSession) {
  {
    MemoryFileSystem fs;  // sin nominal registrado
    auto sc = scenario();
    sc["slots"][0]["set"]["Pw"] = nlohmann::json{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Rig r(sc, fs);
    r.run(5s);
    EXPECT_EQ(nominalReq(r, 3).value("result", ""), "rejected");
  }
  {
    MemoryFileSystem fs;
    recordNominal(fs);
    auto sc = scenario();
    sc["slots"][0]["set"]["Pw"] = nlohmann::json{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    sc["slots"][0]["set"]["V0Set"][3] = 0.0;
    Rig r(sc, fs, /*allow=*/false);
    r.run(5s);
    EXPECT_EQ(nominalReq(r, 3).value("result", ""), "rejected");
    EXPECT_EQ(r.backend.v0setWrites(), 0u);
  }
  {
    MemoryFileSystem fs;
    Rig r(scenario(), fs);
    r.run(5s);
    ASSERT_TRUE(r.begin().value("ok", false));  // sesión activa: la restauración del nominal es de la sesión
    EXPECT_EQ(nominalReq(r, 3).value("result", ""), "rejected");
  }
}

// ---- RUp dentro de la sesión de barrido (seguimiento 2026-10-31 (7)) -------------------------------------------------

TEST(BiasSession, RUpLimitsAndPreviousValueArePersisted) {
  BiasSession s;
  std::map<int, double> vm;
  for (const auto& [ch, v] : kNominal) vm[ch] = v - 0.3;
  EXPECT_FALSE(s.checkRUp(0, 1.0));  // sin sesión
  ASSERT_TRUE(s.begin(kNominal, kNominal, vm, "t", 1));
  EXPECT_TRUE(s.checkRUp(0, 0.05));
  EXPECT_TRUE(s.checkRUp(0, 5.0));
  EXPECT_FALSE(s.checkRUp(0, 0.04));
  EXPECT_FALSE(s.checkRUp(0, 5.01));
  EXPECT_FALSE(s.checkRUp(0, NAN));
  EXPECT_FALSE(s.checkRUp(8, 1.0));  // canal fuera de la sesión
  EXPECT_TRUE(s.notePreviousRUp(3, 1.0));
  EXPECT_FALSE(s.notePreviousRUp(3, 0.1));  // solo el primero de la sesión
  EXPECT_DOUBLE_EQ(s.rupPrevious().at(3), 1.0);
  BiasSession t;
  ASSERT_TRUE(t.restore(s.toJson(2)));
  EXPECT_DOUBLE_EQ(t.rupPrevious().at(3), 1.0);
  t.end();
  EXPECT_TRUE(t.rupPrevious().empty());
}

TEST(BiasScanPoller, RUpIsSetReadBackLoggedAndRestoredByBiasRestore) {
  MemoryFileSystem fs;
  Rig r(scenario(), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  EXPECT_NEAR(r.rupRead(0), 1.0, 1e-9);
  auto res = r.rup(0, 2.0);
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_EQ(res.value("result", ""), "verified");
  EXPECT_NEAR(r.rupRead(0), 2.0, 1e-9);
  EXPECT_EQ(r.rup(0, 2.0).value("result", ""), "noop");
  // RUp previo persistido con la sesión
  const auto st = nlohmann::json::parse(fs.content("/state/hvmon_bias_state.json"));
  EXPECT_DOUBLE_EQ(st["rup_previous_V_per_s"]["0"].get<double>(), 1.0);
  // un escalón de V0Set con el RUp de la sesión se verifica igual
  ASSERT_TRUE(r.v0set(0, 65.25).value("ok", false));
  // registro en hv_changes: request, written, verified con param RUp y el valor en V/s
  int verified = 0;
  for (const auto& c : r.cmds)
    if (c.param == "RUp" && c.phase == "verified" && c.ch == 0) {
      ++verified;
      EXPECT_DOUBLE_EQ(c.toJson()["value"].get<double>(), 2.0);
      EXPECT_FALSE(c.reason.empty());
    }
  EXPECT_EQ(verified, 1);
  EXPECT_EQ(r.backend.rupWrites(), 1u);
  res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  EXPECT_NEAR(r.rupRead(0), 1.0, 1e-9);
  EXPECT_NEAR(r.v0(0), 65.75, 1e-6);
  EXPECT_EQ(r.backend.rupWrites(), 2u);
  EXPECT_FALSE(r.poller.bias().active());
  EXPECT_TRUE(r.poller.bias().rupPrevious().empty());
}

TEST(BiasScanPoller, RUpOutOfRangeOrWithoutSessionIsRejectedWithoutWriting) {
  MemoryFileSystem fs;
  Rig q(scenario(), fs, false);
  q.run(5s);
  EXPECT_EQ(q.rup(0, 2.0).value("result", ""), "rejected");  // allow_commands = false
  Rig r(scenario(), fs);
  r.run(5s);
  EXPECT_EQ(r.rup(0, 2.0).value("result", ""), "rejected");  // sin sesión
  ASSERT_TRUE(r.begin().value("ok", false));
  for (double v : {0.04, 5.5, -1.0}) EXPECT_EQ(r.rup(0, v).value("result", ""), "rejected") << v;
  EXPECT_EQ(r.rup(8, 2.0).value("result", ""), "rejected");  // canal fuera de la sesión
  // 0.1 V/s está dentro de 0.05-5 V/s pero bajo el mínimo que informa la placa (A7030DP: 1 V/s)
  const auto low = r.rup(0, 0.1);
  EXPECT_EQ(low.value("result", ""), "rejected");
  EXPECT_NE(low.value("detail", "").find("placa"), std::string::npos) << low.dump();
  EXPECT_EQ(r.backend.rupWrites(), 0u);
  EXPECT_TRUE(r.poller.bias().rupPrevious().empty());
  EXPECT_TRUE(r.restore().value("ok", false));
}

TEST(BiasScanPoller, AbortMidScanRestoresRUpNominalAndPower) {
  MemoryFileSystem fs;
  // de 30 a 50 s CH01 recibe 30 µA extra: trip durante la rampa del escalón
  Rig r(scenario(nlohmann::json::array({{{"t_s", 30.0}, {"type", "current"}, {"slot", 3}, {"ch", 1}, {"extra_uA", 30.0}, {"duration_s", 20.0}}})), fs);
  r.run(5s);
  ASSERT_TRUE(r.begin().value("ok", false));
  for (int ch = 0; ch < 8; ++ch) ASSERT_TRUE(r.rup(ch, 2.5).value("ok", false)) << ch;
  ASSERT_TRUE(r.v0set(0, 65.25).value("ok", false));
  const auto res1 = r.v0set(1, 65.25);  // se dispara a mitad: el barrido aborta (como biasscan.py)
  r.run(30s);
  EXPECT_EQ(r.backend.readChannels(3, {1}, {"Pw"}).value().get("Pw", 1)->raw, 0u) << res1.dump();
  r.run(20s);  // vuelve la corriente normal
  const auto res = r.restore();
  ASSERT_TRUE(res.value("ok", false)) << res.dump();
  for (int ch = 0; ch < 8; ++ch) {
    EXPECT_NEAR(r.rupRead(ch), 1.0, 1e-9) << ch;
    EXPECT_NEAR(r.v0(ch), kNominal.at(ch), 1e-6) << ch;
    EXPECT_EQ(r.backend.readChannels(3, {ch}, {"Pw"}).value().get("Pw", ch)->raw, 1u) << ch;
  }
  EXPECT_EQ(r.backend.rupWrites(), 16u);  // 8 escrituras + 8 restauraciones
  EXPECT_FALSE(r.poller.bias().active());
}

TEST(BiasScanPoller, HvmonStartingWithActiveSessionRestoresPreviousRUp) {
  MemoryFileSystem fs;
  {
    Rig r(scenario(), fs);
    r.run(5s);
    ASSERT_TRUE(r.begin().value("ok", false));
    ASSERT_TRUE(r.rup(2, 3.0).value("ok", false));
    ASSERT_TRUE(r.v0set(2, 65.40).value("ok", false));
  }
  // "Reinicio" de hvmon (p. ej. el barrido murió y nadie llamó a bias restore): la placa quedó con RUp 3 V/s y V0Set bajo el
  // nominal en CH02; el estado persistido tiene el RUp previo.
  auto sc = scenario();
  sc["slots"][0]["set"]["V0Set"][2] = 65.40;
  sc["slots"][0]["set"]["RUp"][2] = 3.0;
  Rig r2(sc, fs);
  EXPECT_TRUE(r2.poller.bias().active());
  r2.run(300s);
  EXPECT_NEAR(r2.rupRead(2), 1.0, 1e-9);
  EXPECT_NEAR(r2.v0(2), 65.90, 1e-6);
  EXPECT_EQ(r2.backend.rupWrites(), 1u);
  EXPECT_FALSE(r2.poller.bias().active());
}
