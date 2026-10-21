// Tests del hilo de sondeo de hvmon con FakeHVBackend y ManualClock (tiempo acelerado): cambios y
// su origen, pérdida de conexión y espera exponencial, ruta de escritura de Pw, offsets.
#include <gtest/gtest.h>

#include <chrono>

#include "daq/core/clock.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/poller.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::hvmon;

namespace {

nlohmann::json scenario(nlohmann::json events, nlohmann::json pw = 0) {
  auto sc = hv::FakeHVBackend::defaultScenario();
  auto& s = sc["slots"][0];
  s["imon_noise_uA"] = 0.0;
  s["vmon_noise_V"] = 0.0;
  s["set"]["Pw"] = std::move(pw);
  sc["events"] = std::move(events);
  return sc;
}

struct Rig {
  ManualClock clock;
  HvmonConfig cfg;
  hv::FakeHVBackend backend;
  BoundedQueue<WriterItem> writer{100000};
  BoundedQueue<PollRequest> requests{16};
  StateHub hub;
  Poller poller;

  explicit Rig(const nlohmann::json& sc, bool allowCommands = false)
      : cfg(makeCfg(allowCommands)), backend(sc, clock), poller(cfg, backend, clock, writer, hub, requests) {}
  static HvmonConfig makeCfg(bool allow) {
    auto c = HvmonConfig::fromJson(nlohmann::json::object()).value();
    c.allowCommands = allow;
    return c;
  }
  SteadyTime t(std::chrono::milliseconds ms) const { return SteadyTime{} + ms; }
  void runUntil(std::chrono::milliseconds ms) {
    while (clock.now_steady() < t(ms)) poller.step();
  }
  std::vector<ChangeRecord> drainChanges(std::vector<CommandRecord>* commands = nullptr) {
    std::vector<ChangeRecord> out;
    while (auto item = writer.try_pop()) {
      if (auto* c = std::get_if<ChangeItem>(&*item)) out.push_back(c->rec);
      if (auto* k = std::get_if<CommandItem>(&*item); k && commands) commands->push_back(k->rec);
    }
    return out;
  }
  std::future<nlohmann::json> pw(int ch, bool on, bool internal = false) {
    PollRequest r;
    r.kind = PollRequest::Kind::Pw;
    r.ch = ch;
    r.on = on;
    r.reason = "test";
    r.requester = "gtest";
    r.internal = internal;
    auto f = r.reply->get_future();
    EXPECT_TRUE(requests.try_push(std::move(r), DropPolicy::DropNewest));
    return f;
  }
};

const ChangeRecord* find(const std::vector<ChangeRecord>& v, int ch, const std::string& param) {
  for (const auto& c : v)
    if (c.ch == ch && c.param == param) return &c;
  return nullptr;
}

}  // namespace

TEST(HvmonPoller, FastEverySecondFullEveryMinute) {
  Rig r(scenario(nlohmann::json::array()));
  r.runUntil(125s);  // termina justo antes de la lectura de t = 125 s
  EXPECT_EQ(r.poller.stats().fastReads, 125u);  // t = 0, 1, ..., 124
  EXPECT_EQ(r.poller.stats().fullReads, 3u);    // t = 0, 60, 120
  EXPECT_EQ(r.poller.stats().readErrors, 0u);
  EXPECT_TRUE(r.drainChanges().empty());
  const auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  EXPECT_TRUE(st["connection"]["connected"].get<bool>());
  EXPECT_FALSE(st["connection"]["stale"].get<bool>());
  EXPECT_EQ(st["channels"].size(), 12u);
  EXPECT_DOUBLE_EQ(st["channels"][6]["values"]["V0Set_V"].get<double>(), 66.15);
  EXPECT_DOUBLE_EQ(st["channels"][0]["values"]["I0Set_A"].get<double>(), 1.5e-05);
  EXPECT_EQ(st["channels"][0]["auto_retry"]["state"], "DISABLED");  // sin auto_retry configurado
  EXPECT_TRUE(st["temperature"].is_null());
}

TEST(HvmonPoller, ScriptedTripHasOriginTrip) {
  // CH2 encendido; +20 µA desde t = 20 s supera I0Set = 15 µA; Trip = 5 s.
  Rig r(scenario({{{"t_s", 20.0}, {"type", "current"}, {"ch", 2}, {"extra_uA", 20.0}, {"duration_s", 30.0}}},
                 nlohmann::json::array({0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0})));
  r.runUntil(30s);
  const auto ch = r.drainChanges();
  const ChangeRecord* pw = find(ch, 2, "Pw");
  ASSERT_NE(pw, nullptr);
  EXPECT_EQ(pw->origin, "trip");
  EXPECT_EQ(pw->oldValue, 1);
  EXPECT_EQ(pw->newValue, 0);
  bool statusTrip = false;
  for (const auto& c : ch)
    if (c.ch == 2 && c.param == "Status" && c.origin == "trip" && (c.newValue.get<unsigned>() & 0x200u)) statusTrip = true;
  EXPECT_TRUE(statusTrip);
  EXPECT_GT(pw->tUtcNs, pw->tPrevReadUtcNs);
  EXPECT_FALSE(pw->commandId.has_value());
}

TEST(HvmonPoller, ExternalChangesOfV0SetAndPw) {
  Rig r(scenario({{{"t_s", 10.5}, {"type", "set"}, {"ch", 4}, {"param", "V0Set"}, {"value", 65.5}},
                  {{"t_s", 12.5}, {"type", "set"}, {"ch", 5}, {"param", "Pw"}, {"value", 1}},
                  {{"t_s", 70.5}, {"type", "set"}, {"ch", 1}, {"param", "RUp"}, {"value", 2}}}));
  r.runUntil(125s);
  const auto ch = r.drainChanges();
  const ChangeRecord* v = find(ch, 4, "V0Set");
  ASSERT_NE(v, nullptr);
  EXPECT_EQ(v->origin, "external");
  EXPECT_DOUBLE_EQ(v->oldValue.get<double>(), 65.75);
  EXPECT_DOUBLE_EQ(v->newValue.get<double>(), 65.5);
  EXPECT_EQ(v->unit, "V");
  EXPECT_EQ(v->readKind, "fast");
  const ChangeRecord* p = find(ch, 5, "Pw");
  ASSERT_NE(p, nullptr);
  EXPECT_EQ(p->origin, "external");
  const ChangeRecord* rup = find(ch, 1, "RUp");  // solo en lecturas completas
  ASSERT_NE(rup, nullptr);
  EXPECT_EQ(rup->readKind, "full");
  EXPECT_EQ(rup->unit, "V/s");
  // El anillo del estado tiene los mismos cambios.
  const auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  EXPECT_EQ(st["changes"].size(), ch.size());
}

TEST(HvmonPoller, ConnectionLossBackoffStalenessAndNoSpuriousChanges) {
  // Pérdida de t = 10.5 a 40.5 s. La lectura de t = 11 falla; reintentos a +1, +2, +4, +8, +16 s.
  Rig r(scenario({{{"t_s", 10.5}, {"type", "connection_loss"}, {"duration_s", 30.0}}}));
  r.runUntil(20s);
  auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  EXPECT_FALSE(st["connection"]["connected"].get<bool>());
  EXPECT_TRUE(st["connection"]["stale"].get<bool>());
  EXPECT_FALSE(st["connection"]["last_error"].get<std::string>().empty());
  r.runUntil(60s);
  const auto& at = r.poller.connectAttemptTimes();
  std::vector<double> secs;
  for (auto t : at) secs.push_back(std::chrono::duration<double>(t - SteadyTime{}).count());
  ASSERT_EQ(secs.size(), 6u);
  EXPECT_NEAR(secs[0], 0.0, 1e-6);
  EXPECT_NEAR(secs[1], 12.0, 1e-6);
  EXPECT_NEAR(secs[2], 14.0, 1e-6);
  EXPECT_NEAR(secs[3], 18.0, 1e-6);
  EXPECT_NEAR(secs[4], 26.0, 1e-6);
  EXPECT_NEAR(secs[5], 42.0, 1e-6);
  EXPECT_EQ(r.poller.stats().readErrors, 1u);
  EXPECT_EQ(r.poller.stats().reconnections, 1u);
  EXPECT_EQ(r.poller.stats().connectFailures, 4u);
  st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 2, false);
  EXPECT_TRUE(st["connection"]["connected"].get<bool>());
  EXPECT_FALSE(st["connection"]["stale"].get<bool>());
  EXPECT_TRUE(r.drainChanges().empty());  // nada cambió durante el corte: sin cambios espurios
}

TEST(HvmonPoller, ChangeDuringDisconnectIsFlagged) {
  Rig r(scenario({{{"t_s", 10.5}, {"type", "connection_loss"}, {"duration_s", 5.0}},
                  {{"t_s", 12.0}, {"type", "set"}, {"ch", 7}, {"param", "V0Set"}, {"value", 60.0}}}));
  r.runUntil(30s);
  const auto ch = r.drainChanges();
  ASSERT_EQ(ch.size(), 1u);
  EXPECT_EQ(ch[0].param, "V0Set");
  EXPECT_TRUE(ch[0].duringDisconnect);
  EXPECT_EQ(ch[0].origin, "external");
}

TEST(HvmonPoller, PwRefusedWhenCommandsDisabled) {
  Rig r(scenario(nlohmann::json::array()), /*allowCommands=*/false);
  r.runUntil(3s);
  auto f = r.pw(3, true);
  r.runUntil(5s);
  ASSERT_EQ(f.wait_for(0s), std::future_status::ready);
  const auto j = f.get();
  EXPECT_FALSE(j["ok"].get<bool>());
  EXPECT_EQ(j["result"], "rejected");
  EXPECT_EQ(r.backend.pwWrites(), 0u);
  std::vector<CommandRecord> cmds;
  r.drainChanges(&cmds);
  ASSERT_EQ(cmds.size(), 2u);
  EXPECT_EQ(cmds[0].phase, "request");
  EXPECT_EQ(cmds[1].phase, "rejected");
  EXPECT_EQ(cmds[1].reason, "test");
}

TEST(HvmonPoller, PwOnRampsVerifiesAndChangesCarryCommandId) {
  Rig r(scenario(nlohmann::json::array()), /*allowCommands=*/true);
  r.runUntil(5s);
  auto f = r.pw(3, true);  // CH3: V0Set 65.85 V, RUp 1 V/s
  r.runUntil(60s);
  EXPECT_NE(f.wait_for(0s), std::future_status::ready);  // todavía en rampa
  r.runUntil(75s);
  ASSERT_EQ(f.wait_for(0s), std::future_status::ready);
  const auto j = f.get();
  EXPECT_TRUE(j["ok"].get<bool>()) << j.dump();
  EXPECT_EQ(j["result"], "verified");
  EXPECT_EQ(r.backend.pwWrites(), 1u);
  std::vector<CommandRecord> cmds;
  const auto ch = r.drainChanges(&cmds);
  const ChangeRecord* pw = find(ch, 3, "Pw");
  ASSERT_NE(pw, nullptr);
  EXPECT_EQ(pw->origin, "hvmon_command");
  ASSERT_TRUE(pw->commandId.has_value());
  EXPECT_EQ(*pw->commandId, j["command_id"].get<std::uint64_t>());
  for (const auto& c : ch) {
    if (c.ch == 3) {
      EXPECT_EQ(c.origin, "hvmon_command") << c.param;
    }
  }
  std::vector<std::string> phases;
  for (const auto& c : cmds) phases.push_back(c.phase);
  EXPECT_EQ(phases, (std::vector<std::string>{"request", "written", "verified"}));
  EXPECT_NEAR(r.backend.trueVoltage(3, 3), 65.85, 1e-6);
}

TEST(HvmonPoller, PwAbortedWhenChannelChangedBeforeWrite) {
  Rig r(scenario({{{"t_s", 10.3}, {"type", "set"}, {"ch", 3}, {"param", "V0Set"}, {"value", 50.0}}}), true);
  r.runUntil(10s);  // última lectura publicada a t = 10 s
  r.clock.advance(400ms);  // el cambio de t = 10.3 s todavía no fue leído
  auto f = r.pw(3, true);
  r.poller.step();  // la petición se atiende antes de la próxima lectura
  ASSERT_EQ(f.wait_for(0s), std::future_status::ready);
  const auto j = f.get();
  EXPECT_FALSE(j["ok"].get<bool>());
  EXPECT_EQ(j["result"], "aborted");
  EXPECT_NE(j["detail"].get<std::string>().find("V0Set"), std::string::npos);
  EXPECT_EQ(r.backend.pwWrites(), 0u);
}

TEST(HvmonPoller, VerificationFailureIsReported) {
  auto sc = scenario(nlohmann::json::array());
  sc["pw_write_ignored"] = true;  // la escritura "funciona" pero el canal no prende
  Rig r(sc, true);
  r.runUntil(3s);
  auto f = r.pw(0, true);  // timeout = 65.75 / 1 + 15 = 80.75 s
  r.runUntil(80s);
  EXPECT_NE(f.wait_for(0s), std::future_status::ready);
  r.runUntil(90s);
  ASSERT_EQ(f.wait_for(0s), std::future_status::ready);
  const auto j = f.get();
  EXPECT_FALSE(j["ok"].get<bool>());
  EXPECT_EQ(j["result"], "verify_failed");
  EXPECT_EQ(r.backend.pwWrites(), 1u);
}

TEST(HvmonPoller, InternalRequestRefusedWhileAnotherChannelRamps) {
  Rig r(scenario({{{"t_s", 2.0}, {"type", "set"}, {"ch", 1}, {"param", "Pw"}, {"value", 1}}}), true);
  r.runUntil(5s);  // CH1 en rampa (externa)
  auto f = r.pw(2, true, /*internal=*/true);
  r.runUntil(6s);
  const auto j = f.get();
  EXPECT_EQ(j["result"], "busy");  // fase 2: una petición interna espera en vez de fallar
  EXPECT_NE(j["detail"].get<std::string>().find("rampa"), std::string::npos);
  auto g = r.pw(2, true, /*internal=*/false);  // una petición del operador sí procede
  r.runUntil(7s);
  EXPECT_NE(g.wait_for(0s), std::future_status::ready);
  EXPECT_EQ(r.backend.pwWrites(), 1u);
}

TEST(HvmonPoller, MeasureOffsetsRefusesPwOnAndAverages) {
  Rig r(scenario(nlohmann::json::array(), nlohmann::json::array({1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})));
  r.runUntil(5s);
  PollRequest q;
  q.kind = PollRequest::Kind::Offsets;
  q.seconds = 60;
  auto f = q.reply->get_future();
  ASSERT_TRUE(r.requests.try_push(std::move(q), DropPolicy::DropNewest));
  r.runUntil(70s);
  ASSERT_EQ(f.wait_for(0s), std::future_status::ready);
  const auto j = f.get();
  ASSERT_TRUE(j["ok"].get<bool>()) << j.dump();
  const auto& off = j["offsets"];
  EXPECT_EQ(off["refused_pw_on"], nlohmann::json::array({0, 1}));
  EXPECT_EQ(off["channels"].size(), 10u);
  EXPECT_EQ(off["n_reads"].get<int>(), 61);  // lecturas de t0 a t0 + 60 s
  // CH5 del guion: IMon offset -0.046 µA, VMon offset 0.0 V, sin ruido.
  for (const auto& c : off["channels"]) {
    if (c["ch"] == 5) {
      EXPECT_NEAR(c["imon_offset_A"].get<double>(), -4.6e-8, 1e-15);
      EXPECT_NEAR(c["imon_std_A"].get<double>(), 0.0, 1e-15);
      EXPECT_NEAR(c["vmon_offset_V"].get<double>(), 0.0, 1e-12);
    }
    if (c["ch"] == 11) {
      EXPECT_NEAR(c["imon_offset_A"].get<double>(), -8.4e-8, 1e-15);
    }
  }
  bool written = false;
  while (auto item = r.writer.try_pop())
    if (std::get_if<OffsetItem>(&*item)) written = true;
  EXPECT_TRUE(written);
}

TEST(HvmonPoller, MeasureOffsetsFailsWhenAllRequestedChannelsOn) {
  Rig r(scenario(nlohmann::json::array(), 1));
  r.runUntil(3s);
  PollRequest q;
  q.kind = PollRequest::Kind::Offsets;
  q.seconds = 10;
  q.channels = {0, 1};
  auto f = q.reply->get_future();
  ASSERT_TRUE(r.requests.try_push(std::move(q), DropPolicy::DropNewest));
  r.runUntil(4s);
  const auto j = f.get();
  EXPECT_FALSE(j["ok"].get<bool>());
  EXPECT_EQ(j["refused"], nlohmann::json::array({0, 1}));
}

TEST(HvmonPoller, ShutdownAnswersPendingRequests) {
  Rig r(scenario(nlohmann::json::array()), true);
  r.runUntil(3s);
  auto f = r.pw(0, true);
  r.runUntil(4s);  // comando en verificación
  auto g = r.pw(1, true);  // queda en la cola
  r.poller.shutdown();
  EXPECT_EQ(f.get()["error"], "hvmon se está apagando");
  EXPECT_EQ(g.get()["error"], "hvmon se está apagando");
}

TEST(HvmonPoller, FullReadIncludesBoardParametersAndBoardTemperature) {
  auto sc = scenario({{{"t_s", 30.0}, {"type", "board"}, {"param", "Temp"}, {"value", 27.25}}});
  Rig r(sc);
  r.runUntil(61s);
  auto st = r.hub.compose(r.cfg, r.clock.now_utc_ns(), 1, false);
  ASSERT_TRUE(st["board"].is_object());
  EXPECT_DOUBLE_EQ(st["board_temperature_degC"].get<double>(), 27.25);  // lectura completa de t = 60 s
  EXPECT_DOUBLE_EQ(st["board"]["values"]["HVMax_V"].get<double>(), 3127.0);
  EXPECT_EQ(st["board"]["values"]["BdStatus"], 0);
  int fulls = 0;
  bool hasBoard = true;
  while (auto item = r.writer.try_pop())
    if (auto* f = std::get_if<FullRows>(&*item)) {
      ++fulls;
      hasBoard = hasBoard && f->board.has_value() && f->board->params.size() == 4 &&
                 !f->board->get("ClrAlarm");  // excluido por board_params_exclude
    }
  EXPECT_EQ(fulls, 2);
  EXPECT_TRUE(hasBoard);
}

// DAQ2-10: el apagado se verifica con el bit on de Status en 0 (y Pw = 0), no con VMon. En el
// mainframe real CH08 (V0Set 0) tiene 1.5 V de offset de VMon y los canales del detector descargan
// la capacidad del PCB con τ ≈ 115 s: con la regla anterior (VMon ≤ 2 V) CH07 dio verify_failed.
TEST(HvmonPoller, PwOffVerifiedByStatusOnBitNotByVMon) {
  auto sc = scenario(nlohmann::json::array(), 1);
  sc["slots"][0]["vmon_offset_V"] = 5.0;  // VMon nunca baja de 5 V aunque la tensión real sea 0
  Rig r(sc, /*allowCommands=*/true);
  r.runUntil(5s);
  auto f = r.pw(3, false);
  r.runUntil(20s);
  ASSERT_EQ(f.wait_for(0s), std::future_status::ready);
  const auto j = f.get();
  EXPECT_TRUE(j["ok"].get<bool>()) << j.dump();
  EXPECT_EQ(j["result"], "verified");
}
