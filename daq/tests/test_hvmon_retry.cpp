// Tests del reintento automático (fase 2) a través del hilo de sondeo completo, con
// FakeHVBackend y ManualClock (tiempo acelerado). Canales de prueba CH08/CH09 (libres) en
// test_mode salvo donde se indica.
#include <gtest/gtest.h>

#include <chrono>

#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/poller.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::hvmon;

namespace {

const std::string kStateFile = "/state/hvmon_retry_state.json";

// CH08 y CH09 como cargas de prueba: 10 V, RUp 5 V/s (rampa de 2 s), I0Set 5 µA, Trip 10 s, carga 1 µA.
nlohmann::json scenario(nlohmann::json events, nlohmann::json pw = nlohmann::json::array({0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0}),
                        double rup = 5.0) {
  auto sc = hv::FakeHVBackend::defaultScenario();
  auto& s = sc["slots"][0];
  s["imon_noise_uA"] = 0.0;
  s["vmon_noise_V"] = 0.0;
  s["set"]["Pw"] = std::move(pw);
  s["set"]["V0Set"] = nlohmann::json::array({65.75, 65.75, 65.90, 65.85, 65.75, 65.70, 66.15, 65.90, 10, 10, 0, 0});
  s["set"]["RUp"] = nlohmann::json::array({1, 1, 1, 1, 1, 1, 1, 1, rup, rup, 5, 5});
  s["load_uA"] = nlohmann::json::array({2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0});
  sc["events"] = std::move(events);
  return sc;
}

nlohmann::json retryCfg(nlohmann::json extra = nlohmann::json::object()) {
  nlohmann::json r{{"enabled_channels", {8, 9}},
                   {"test_mode", {{"enabled", true}, {"channels", {8, 9}}, {"wait_min", {1, 2, 3}}, {"stable_min", 1}}},
                   {"state_file", kStateFile}};
  for (auto it = extra.begin(); it != extra.end(); ++it) r[it.key()] = it.value();
  return r;
}

struct RRig {
  ManualClock clock;
  HvmonConfig cfg;
  hv::FakeHVBackend backend;
  BoundedQueue<WriterItem> writer{200000};
  BoundedQueue<PollRequest> requests{16};
  BoundedQueue<NotifyEvent> notify{1000};
  StateHub hub;
  MemoryFileSystem& fs;
  Poller poller;
  std::vector<ChangeRecord> retryChanges;
  std::vector<NotifyEvent> events;

  RRig(const nlohmann::json& sc, const nlohmann::json& rcfg, MemoryFileSystem& memfs, UtcNs startUtc = 1792281600LL * 1000000000LL)
      : clock(startUtc), cfg(makeCfg(rcfg)), backend(sc, clock), fs(memfs), poller(cfg, backend, clock, writer, hub, requests) {
    poller.enableRetryPersistence(fs);
    poller.setNotifyQueue(&notify);
  }
  static HvmonConfig makeCfg(const nlohmann::json& rcfg) {
    auto c = HvmonConfig::fromJson({{"allow_commands", true}, {"auto_retry", rcfg}});
    EXPECT_TRUE(c.ok()) << (c.ok() ? "" : c.error().message);
    return c.value();
  }
  void runUntil(std::chrono::milliseconds ms) {
    while (clock.now_steady() < SteadyTime{} + ms) poller.step();
    drain();
  }
  void drain() {
    while (auto item = writer.try_pop())
      if (auto* c = std::get_if<ChangeItem>(&*item); c && c->rec.origin == "hvmon_retry") retryChanges.push_back(c->rec);
    while (auto e = notify.try_pop()) events.push_back(std::move(*e));
  }
  std::vector<std::string> path(int ch) const {
    std::vector<std::string> out;
    for (const auto& c : retryChanges)
      if (c.ch == ch) out.push_back(c.oldValue.get<std::string>() + "->" + c.newValue.get<std::string>());
    return out;
  }
  int countEvents(const std::string& kind, int ch) const {
    int n = 0;
    for (const auto& e : events) n += (e.kind == kind && e.ch == ch) ? 1 : 0;
    return n;
  }
  RetryState st(int ch) const { return poller.retry().state(ch); }
};

using V = std::vector<std::string>;

}  // namespace

TEST(HvmonRetry, TripAndSuccessAtFirstAttempt) {
  MemoryFileSystem fs;
  // +20 µA en CH08 de t = 20 a 35 s: trip a t = 30 s (Trip = 10 s).
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}), retryCfg(), fs);
  r.runUntil(29s);
  EXPECT_EQ(r.st(8), RetryState::Normal);
  r.runUntil(32s);
  EXPECT_EQ(r.st(8), RetryState::Wait);
  r.runUntil(89s);
  EXPECT_EQ(r.backend.pwWrites(), 0u);  // la espera de test_mode es 1 min
  r.runUntil(95s);
  EXPECT_EQ(r.backend.pwWrites(), 1u);
  EXPECT_EQ(r.st(8), RetryState::Stabilizing);  // rampa de 2 s ya completa
  r.runUntil(200s);
  EXPECT_EQ(r.st(8), RetryState::Normal);
  EXPECT_EQ(r.path(8), (V{"UNINIT->NORMAL", "NORMAL->WAIT(1)", "WAIT(1)->RETRY(1)", "RETRY(1)->STABILIZING(1)",
                          "STABILIZING(1)->NORMAL"}));
  EXPECT_EQ(r.poller.retry().attempt(8), 1);  // el contador se conserva hasta reset_after_h
  EXPECT_EQ(r.countEvents("trip", 8), 1);
  EXPECT_EQ(r.countEvents("recovery", 8), 1);
  EXPECT_EQ(r.countEvents("gave_up", 8), 0);
  EXPECT_EQ(r.st(9), RetryState::Normal);
  EXPECT_EQ(r.st(3), RetryState::Disabled);
  // Persistencia en cada transición.
  const auto saved = nlohmann::json::parse(fs.content(kStateFile));
  EXPECT_EQ(saved["channels"]["8"]["state"], "NORMAL");
  // Cada transición es un registro con origen hvmon_retry, k y motivo.
  const ChangeRecord* toWait = nullptr;
  for (const auto& c : r.retryChanges)
    if (c.ch == 8 && c.newValue == "WAIT(1)") toWait = &c;
  ASSERT_NE(toWait, nullptr);
  EXPECT_EQ(toWait->detail["k"], 1);
  EXPECT_EQ(toWait->detail["reason"], "trip");
  // Estado publicado.
  const auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  EXPECT_EQ(st["channels"][8]["auto_retry"]["state"], "NORMAL");
  EXPECT_EQ(st["channels"][8]["auto_retry"]["last_outcome"], "recuperado en el intento 1");
}

TEST(HvmonRetry, TripsUntilGaveUpWithOneEmailEvent) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 100000.0}}}),
         retryCfg(), fs);
  r.runUntil(15min);
  EXPECT_EQ(r.st(8), RetryState::GaveUp);
  EXPECT_EQ(r.backend.pwWrites(), 3u);  // un intento por espera (1, 2, 3 min)
  EXPECT_EQ(r.path(8), (V{"UNINIT->NORMAL", "NORMAL->WAIT(1)", "WAIT(1)->RETRY(1)", "RETRY(1)->STABILIZING(1)",
                          "STABILIZING(1)->WAIT(2)", "WAIT(2)->RETRY(2)", "RETRY(2)->STABILIZING(2)",
                          "STABILIZING(2)->WAIT(3)", "WAIT(3)->RETRY(3)", "RETRY(3)->STABILIZING(3)",
                          "STABILIZING(3)->GAVE_UP"}));
  EXPECT_EQ(r.countEvents("gave_up", 8), 1);
  r.runUntil(40min);  // sigue apagado y sin más intentos ni correos
  EXPECT_EQ(r.st(8), RetryState::GaveUp);
  EXPECT_EQ(r.backend.pwWrites(), 3u);
  EXPECT_EQ(r.countEvents("gave_up", 8), 1);
  const NotifyEvent* g = nullptr;
  for (const auto& e : r.events)
    if (e.kind == "gave_up") g = &e;
  ASSERT_NE(g, nullptr);
  EXPECT_GE(g->timeline.size(), 7u);  // trips y reintentos del episodio
  EXPECT_TRUE(g->context.contains("VMon_V"));
  EXPECT_TRUE(g->context.contains("board_temperature_degC"));
}

TEST(HvmonRetry, ManualActionDuringWaitGoesToManualWithoutRetry) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}},
                   {{"t_s", 50.0}, {"type", "set"}, {"ch", 8}, {"param", "V0Set"}, {"value", 9.5}}}),
         retryCfg(), fs);
  r.runUntil(5min);
  EXPECT_EQ(r.st(8), RetryState::Manual);
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  EXPECT_EQ(r.path(8).back(), "WAIT(1)->MANUAL");
  // GAVE_UP o MANUAL solo se dejan por acción manual: encender a mano y quedar estable -> NORMAL.
}

TEST(HvmonRetry, ManualPwOnFromWaitThenStableReturnsToNormal) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}},
                   {{"t_s", 50.0}, {"type", "set"}, {"ch", 8}, {"param", "Pw"}, {"value", 1}}}),
         retryCfg(), fs);
  r.runUntil(55s);
  EXPECT_EQ(r.st(8), RetryState::Manual);
  r.runUntil(3min);
  EXPECT_EQ(r.st(8), RetryState::Normal);  // estable stable_min (1 min de test_mode)
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  EXPECT_EQ(r.poller.retry().attempt(8), 0);
}

TEST(HvmonRetry, ManualActionDuringStabilizing) {
  MemoryFileSystem fs;
  // Trip a 30 s, reintento a 90 s, estabilizando desde ~92 s; cambio externo de V0Set a 110 s.
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}},
                   {{"t_s", 110.0}, {"type", "set"}, {"ch", 8}, {"param", "V0Set"}, {"value", 9.0}}}),
         retryCfg(), fs);
  r.runUntil(105s);
  EXPECT_EQ(r.st(8), RetryState::Stabilizing);
  r.runUntil(115s);
  EXPECT_EQ(r.st(8), RetryState::Manual);
  EXPECT_EQ(r.path(8).back(), "STABILIZING(1)->MANUAL");
}

TEST(HvmonRetry, ConnectionLossDuringWaitDelaysRetryUntilFreshRead) {
  MemoryFileSystem fs;
  // Trip a 30 s; la espera vence a 90 s, pero la conexión se pierde de 60 a 200 s.
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}},
                   {{"t_s", 60.0}, {"type", "connection_loss"}, {"duration_s", 140.0}}}),
         retryCfg(), fs);
  r.runUntil(199s);
  EXPECT_EQ(r.st(8), RetryState::Wait);
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  EXPECT_FALSE(r.poller.connected());
  r.runUntil(6min);
  EXPECT_TRUE(r.poller.connected());
  EXPECT_EQ(r.backend.pwWrites(), 1u);  // tras la primera lectura fresca
  // La reconexión (último intento de conexión) fue después del fin del corte (200 s).
  const auto reconnectAt = r.poller.connectAttemptTimes().back();
  EXPECT_GE(reconnectAt, SteadyTime{} + 200s);
  EXPECT_EQ(r.st(8), RetryState::Normal);
}

TEST(HvmonRetry, TwoChannelsTrippingTogetherRampOneAtATime) {
  MemoryFileSystem fs;
  // Ambos disparan a 30 s; rampa lenta (0.5 V/s -> 20 s).
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}},
                   {{"t_s", 20.0}, {"type", "current"}, {"ch", 9}, {"extra_uA", 20.0}, {"duration_s", 15.0}}},
                  nlohmann::json::array({0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0}), 0.5),
         retryCfg({{"max_ramping_channels", 1}}), fs);
  r.runUntil(95s);
  EXPECT_EQ(r.backend.pwWrites(), 1u);
  EXPECT_EQ(r.st(8), RetryState::Retry);
  EXPECT_EQ(r.st(9), RetryState::Wait);
  const auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  EXPECT_NE(st["channels"][9]["auto_retry"]["block_reason"].get<std::string>().find("rampa"), std::string::npos);
  r.runUntil(10min);
  EXPECT_EQ(r.backend.pwWrites(), 2u);
  EXPECT_EQ(r.backend.writeLog(), (V{"Pw:8=1", "Pw:9=1"}));
  EXPECT_EQ(r.st(8), RetryState::Normal);
  EXPECT_EQ(r.st(9), RetryState::Normal);
}

TEST(HvmonRetry, RestartDuringWait2ResumesAndAfterWindowGoesManual) {
  MemoryFileSystem fs;
  UtcNs savedUtc = 0;
  {
    // Trip a 30 s, reintento a 90 s que vuelve a disparar (excursión larga) -> WAIT(2).
    RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 100000.0}}}),
           retryCfg(), fs);
    r.runUntil(2min);
    ASSERT_EQ(r.st(8), RetryState::Wait);
    ASSERT_EQ(r.poller.retry().attempt(8), 2);
    savedUtc = r.clock.now_utc_ns();
  }
  const auto saved = nlohmann::json::parse(fs.content(kStateFile));
  EXPECT_EQ(saved["channels"]["8"]["state"], "WAIT");
  EXPECT_EQ(saved["channels"]["8"]["k"], 2);
  // "Reinicio" 30 s después: CH08 sigue apagado (Pw = 0 en el guion nuevo).
  {
    RRig r(scenario(nlohmann::json::array(), nlohmann::json::array({0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0})), retryCfg(), fs,
           savedUtc + 30LL * 1000000000LL);
    r.runUntil(2s);
    EXPECT_EQ(r.st(8), RetryState::Wait);
    EXPECT_EQ(r.poller.retry().attempt(8), 2);
    r.runUntil(3min);  // la espera de 2 min desde el trip vence: reintento 2
    EXPECT_EQ(r.backend.pwWrites(), 1u);
    EXPECT_EQ(r.st(8), RetryState::Normal);
  }
  // Otro "reinicio" con un WAIT persistido viejo (más de wait + 10 min): MANUAL, sin acción.
  {
    MemoryFileSystem fs2;
    auto old = saved;
    old["channels"]["8"]["state"] = "WAIT";
    old["channels"]["8"]["k"] = 2;
    ASSERT_TRUE(fs2.writeAtomic(kStateFile, old.dump(), 0600).ok());
    RRig r(scenario(nlohmann::json::array(), nlohmann::json::array({0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0})), retryCfg(), fs2,
           savedUtc + 13LL * 60LL * 1000000000LL);
    r.runUntil(5min);
    EXPECT_EQ(r.st(8), RetryState::Manual);
    EXPECT_EQ(r.backend.pwWrites(), 0u);
  }
}

TEST(HvmonRetry, ChannelsOffAtStartAreManual) {
  MemoryFileSystem fs;
  RRig r(scenario(nlohmann::json::array(), nlohmann::json::array({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})), retryCfg(), fs);
  r.runUntil(10min);
  EXPECT_EQ(r.st(8), RetryState::Manual);
  EXPECT_EQ(r.st(9), RetryState::Manual);
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  EXPECT_EQ(r.path(8), (V{"UNINIT->MANUAL"}));
}

TEST(HvmonRetry, TestModeRefusedOnDetectorChannelsUnlessAllowed) {
  auto bad = HvmonConfig::fromJson(
      {{"allow_commands", true},
       {"auto_retry", {{"enabled_channels", {3}}, {"test_mode", {{"enabled", true}, {"channels", {3}}}}}}});
  ASSERT_FALSE(bad.ok());
  EXPECT_NE(bad.error().message.find("test_mode_allow_detector_channels"), std::string::npos);
  auto ok = HvmonConfig::fromJson({{"allow_commands", true},
                                   {"auto_retry",
                                    {{"enabled_channels", {3}},
                                     {"test_mode_allow_detector_channels", true},
                                     {"test_mode", {{"enabled", true}, {"channels", {3}}}}}}});
  EXPECT_TRUE(ok.ok());
  // auto_retry con canales habilitados exige allow_commands.
  EXPECT_FALSE(HvmonConfig::fromJson({{"auto_retry", {{"enabled_channels", {8}}}}}).ok());
  // Por defecto ningún canal tiene reintento.
  EXPECT_TRUE(HvmonConfig::fromJson(nlohmann::json::object()).value().retry.enabledChannels.empty());
}

TEST(HvmonRetry, GecoRunningWarnsButDoesNotInhibitByDefault) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}), retryCfg(), fs);
  r.hub.gecoRunning.store(true);
  r.runUntil(95s);
  EXPECT_EQ(r.backend.pwWrites(), 1u);
}

TEST(HvmonRetry, GecoRunningInhibitsWhenConfigured) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}),
         retryCfg({{"inhibit_when_geco_running", true}}), fs);
  r.hub.gecoRunning.store(true);
  r.runUntil(5min);
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  EXPECT_EQ(r.st(8), RetryState::Wait);
  const auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  EXPECT_NE(st["channels"][8]["auto_retry"]["block_reason"].get<std::string>().find("GECO"), std::string::npos);
  r.hub.gecoRunning.store(false);
  r.runUntil(6min);
  EXPECT_EQ(r.backend.pwWrites(), 1u);
}

TEST(HvmonRetry, ClearAlarmBeforePwWhenEnabled) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}),
         retryCfg({{"clear_alarm_before_retry", true}}), fs);
  r.runUntil(95s);
  EXPECT_EQ(r.backend.writeLog(), (V{"ClrAlarm", "Pw:8=1"}));
  EXPECT_EQ(r.backend.clearAlarmWrites(), 1u);
}

TEST(HvmonRetry, ClearAlarmNotSentByDefaultAndRefusedWithoutAllowCommands) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}), retryCfg(), fs);
  r.runUntil(95s);
  EXPECT_EQ(r.backend.clearAlarmWrites(), 0u);
  // Con allow_commands = false no hay auto-retry posible (la configuración se rechaza).
  EXPECT_FALSE(HvmonConfig::fromJson({{"auto_retry", {{"enabled_channels", {8}}, {"clear_alarm_before_retry", true}}}}).ok());
}

TEST(HvmonRetry, DisabledChannelTripIsOnlyLogged) {
  MemoryFileSystem fs;
  // CH09 dispara pero solo CH08 está habilitado.
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 9}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}),
         retryCfg({{"enabled_channels", {8}}, {"test_mode", {{"enabled", true}, {"channels", {8}}, {"wait_min", {1, 2, 3}}, {"stable_min", 1}}}}),
         fs);
  r.runUntil(10min);
  EXPECT_EQ(r.st(9), RetryState::Disabled);
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  EXPECT_TRUE(r.path(9).empty());
  // XPLUS-1: sin reintento, pero el trip sí genera un evento de aviso (uno solo por trip)
  EXPECT_EQ(r.countEvents("trip", 9), 1);
  EXPECT_EQ(r.countEvents("gave_up", 9), 0);
}

// XPLUS-1: con enabled_channels vacío (la configuración real) el trip de un canal se avisa igual, una vez por trip; si el
// canal se vuelve a encender y dispara otra vez, hay un segundo aviso. Ninguna escritura de hvmon.
TEST(HvmonRetry, TripEventWithoutAnyAutoRetryChannel) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}},
                   {{"t_s", 150.0}, {"type", "set"}, {"ch", 8}, {"param", "Pw"}, {"value", 1}},  // encendido manual (escritura externa)
                   {{"t_s", 200.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}),
         nlohmann::json{{"enabled_channels", nlohmann::json::array()}, {"state_file", kStateFile}}, fs);
  r.runUntil(100s);
  EXPECT_EQ(r.st(8), RetryState::Disabled);
  EXPECT_EQ(r.countEvents("trip", 8), 1);
  r.runUntil(150s);
  EXPECT_EQ(r.countEvents("trip", 8), 1);  // sigue apagado: sin aviso repetido
  // encendido manual a t = 150 s (evento "set" del guion) y segundo trip a t = 210 s
  r.runUntil(300s);
  EXPECT_EQ(r.countEvents("trip", 8), 2);
  EXPECT_EQ(r.backend.pwWrites(), 0u);  // hvmon no escribió
  EXPECT_EQ(r.countEvents("trip", 9), 0);
}

TEST(HvmonRetry, ResetAfterHoursOfNormal) {
  MemoryFileSystem fs;
  RRig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 8}, {"extra_uA", 20.0}, {"duration_s", 15.0}}}),
         retryCfg({{"reset_after_h", 0.05}}), fs);  // 3 min
  r.runUntil(200s);
  EXPECT_EQ(r.poller.retry().attempt(8), 1);
  r.runUntil(8min);
  EXPECT_EQ(r.poller.retry().attempt(8), 0);
  EXPECT_EQ(r.st(8), RetryState::Normal);
}
