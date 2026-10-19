// Tests varios de hvmon: configuración, serialización de lecturas, detección de procesos,
// enlace HTTP permitido, backend de temperatura fake.
#include <gtest/gtest.h>

#include <sys/stat.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>

#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hv/reading_json.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/system_probe.h"
#include "daq/temp/temp_backend.h"

using namespace std::chrono_literals;
using namespace std::string_literals;
using namespace daq;

TEST(HvmonConfig, DefaultsAreReadOnlyAndValidated) {
  auto c = hvmon::HvmonConfig::fromJson(nlohmann::json::object());
  ASSERT_TRUE(c.ok());
  EXPECT_FALSE(c.value().allowCommands);
  EXPECT_EQ(c.value().httpBind, "127.0.0.1");
  EXPECT_EQ(c.value().channels.size(), 12u);
  EXPECT_DOUBLE_EQ(c.value().fsyncPeriodS, 10.0);
  EXPECT_FALSE(hvmon::HvmonConfig::fromJson({{"http", {{"bind", "0.0.0.0"}}}}).ok());
  EXPECT_FALSE(hvmon::HvmonConfig::fromJson({{"http", {{"bind", "10.3.84.38"}}}}).ok());
  EXPECT_FALSE(hvmon::HvmonConfig::fromJson({{"fast_params", {"VMon"}}}).ok());
  EXPECT_FALSE(hvmon::HvmonConfig::fromJson({{"trip_flags", {"no_existe"}}}).ok());
  EXPECT_FALSE(hvmon::HvmonConfig::fromJson({{"fsync_period_s", 0}}).ok());
}

TEST(HvmonConfig, HttpBindResolution) {
  EXPECT_EQ(hvmon::resolveHttpBind("127.0.0.1").value(), "127.0.0.1");
  EXPECT_EQ(hvmon::resolveHttpBind("localhost").value(), "127.0.0.1");
  EXPECT_FALSE(hvmon::resolveHttpBind("0.0.0.0").ok());
  auto t = hvmon::resolveHttpBind("tailscale");
  if (t.ok()) {
    EXPECT_EQ(t.value().rfind("100.", 0), 0u);  // red CGNAT de Tailscale
  }
}

TEST(ReadingJson, RoundTripAndColumnNames) {
  ManualClock clock;
  hv::FakeHVBackend be(hv::FakeHVBackend::defaultScenario(), clock);
  ASSERT_TRUE(be.connect().ok());
  auto r = be.readChannels(3, {0, 8}, {"V0Set", "IMon", "RUp", "Trip", "Status", "Pw"}).value();
  const auto j = hv::readingToJson(r, hv::StatusBitMap::wrapperRev13());
  ASSERT_EQ(j["params"].size(), 6u);
  EXPECT_EQ(j["params"][0]["name"], "V0Set");  // orden preservado
  auto back = hv::readingFromJson(j);
  ASSERT_TRUE(back.ok()) << back.error().message;
  EXPECT_EQ(back.value().params.size(), 6u);
  EXPECT_DOUBLE_EQ(back.value().get("V0Set", 0)->boardValue, 65.75);
  EXPECT_DOUBLE_EQ(back.value().get("IMon", 8)->si, r.get("IMon", 8)->si);
  EXPECT_EQ(back.value().get("Pw", 8)->raw, 0u);
  EXPECT_EQ(hv::columnName(r.params[1]), "IMon_A");
  EXPECT_EQ(hv::columnName(r.params[2]), "RUp_V_per_s");
  EXPECT_EQ(hv::columnName(r.params[3]), "Trip_s");
  EXPECT_EQ(hv::columnName(r.params[4]), "Status");
  EXPECT_FALSE(hv::readingFromJson(nlohmann::json::object()).ok());
}

TEST(SystemProbe, MatchesExecutableOrArgv0Only) {
  const std::string root = ::testing::TempDir() + "daq_fake_proc";
  auto mkproc = [&](const std::string& pid, const std::string& cmd, const std::string& exe) {
    ::mkdir((root + "/" + pid).c_str(), 0700);
    std::ofstream o(root + "/" + pid + "/cmdline", std::ios::binary);
    o.write(cmd.data(), static_cast<std::streamsize>(cmd.size()));
    if (!exe.empty()) {
      ASSERT_EQ(::symlink(exe.c_str(), (root + "/" + pid + "/exe").c_str()), 0);
    }
  };
  auto clean = [&] {
    for (const char* pid : {"4242", "4243", "4244"}) {
      std::remove((root + "/" + pid + "/cmdline").c_str());
      std::remove((root + "/" + pid + "/exe").c_str());
      ::rmdir((root + "/" + pid).c_str());
    }
    ::rmdir(root.c_str());
  };
  clean();
  ::mkdir(root.c_str(), 0700);
  // Un proceso que solo menciona el patrón en sus argumentos (grep, script de muestreo) no cuenta.
  mkproc("4242", "/bin/bash\0-c\0pgrep bin/CAENGECO2020; geco_running\0"s, "/usr/bin/bash");
  EXPECT_FALSE(hvmon::processRunning("bin/CAENGECO2020", root));
  // argv[0] con el ejecutable del snap.
  mkproc("4243", "/snap/caengeco2020/x1/bin/CAENGECO2020\0--algo\0"s, "");
  EXPECT_TRUE(hvmon::processRunning("bin/CAENGECO2020", root));
  std::remove((root + "/4243/cmdline").c_str());
  // Solo el enlace exe (argv[0] cambiado por el programa).
  mkproc("4244", "geco\0"s, "/snap/caengeco2020/x1/bin/CAENGECO2020");
  EXPECT_TRUE(hvmon::processRunning("bin/CAENGECO2020", root));
  EXPECT_FALSE(hvmon::processRunning("bin/NOEXISTE", root));
  clean();
  EXPECT_FALSE(hvmon::processRunning("no-hay-tal-proceso-en-este-host"));
}

TEST(FakeTemp, SineAndFailure) {
  ManualClock clock;
  temp::FakeTempBackend t(clock, {{"sensors", {{{"id", "s1"}, {"base_degC", 20.0}, {"amplitude_degC", 2.0},
                                                {"period_s", 40.0}}}}});
  auto r0 = t.read();
  ASSERT_TRUE(r0.ok());
  EXPECT_NEAR(r0.value()[0].tDegC, 20.0, 1e-9);
  clock.advance(10s);
  EXPECT_NEAR(t.read().value()[0].tDegC, 22.0, 1e-9);
  EXPECT_EQ(t.read().value()[0].sensorId, "s1");
  temp::FakeTempBackend bad(clock, {{"fail", true}});
  EXPECT_FALSE(bad.read().ok());
  EXPECT_FALSE(temp::makeTempBackend("ds18b20", clock, {}).ok());
  EXPECT_EQ(temp::makeTempBackend("none", clock, {}).value(), nullptr);
}

TEST(HvmonConfig, SimulatedOutputNeverGoesToRealSlowcontrol) {
  using hvmon::HvmonConfig;
  // Fuente real sin test_mode: valores reales por omisión.
  auto real = HvmonConfig::fromJson(nlohmann::json::object());
  ASSERT_TRUE(real.ok());
  EXPECT_EQ(real.value().slowcontrolDir, "/data/daq/slowcontrol");
  EXPECT_FALSE(real.value().simulated);
  // Fake sin claves de salida: directorios de simulación.
  auto fake = HvmonConfig::fromJson(nlohmann::json::object(), true);
  ASSERT_TRUE(fake.ok()) << fake.error().message;
  EXPECT_EQ(fake.value().slowcontrolDir, "/data/daq/sim/slowcontrol");
  EXPECT_EQ(fake.value().stateDir, "/dev/shm/daqv2-sim");
  EXPECT_EQ(fake.value().runtimeDir, hvmon::defaultRuntimeDir() + "-sim");
  EXPECT_TRUE(fake.value().simulated);
  // Fake con la salida real (o debajo de ella, o por un camino con "..") se rechaza.
  for (const char* d : {"/data/daq/slowcontrol", "/data/daq/slowcontrol/", "/data/daq/slowcontrol/2026/x",
                        "/data/daq/sim/../slowcontrol"}) {
    auto r = HvmonConfig::fromJson({{"slowcontrol_dir", d}}, true);
    ASSERT_FALSE(r.ok()) << d;
    EXPECT_NE(r.error().message.find("FakeHVBackend"), std::string::npos);
  }
  // Fake con el estado real: rechazado (lo leen hvctl y daqrun).
  EXPECT_FALSE(HvmonConfig::fromJson({{"state_dir", "/dev/shm/daqv2"}}, true).ok());
  // test_mode con la fuente real: slowcontrol_dir ausente va a simulación; explícito real, rechazado.
  const nlohmann::json retry{{"enabled_channels", {8}}, {"test_mode", {{"enabled", true}, {"channels", {8}}}}};
  auto tm = HvmonConfig::fromJson({{"allow_commands", true}, {"auto_retry", retry}});
  ASSERT_TRUE(tm.ok()) << tm.error().message;
  EXPECT_EQ(tm.value().slowcontrolDir, "/data/daq/sim/slowcontrol");
  EXPECT_TRUE(tm.value().simulated);
  auto tmBad = HvmonConfig::fromJson(
      {{"allow_commands", true}, {"auto_retry", retry}, {"slowcontrol_dir", "/data/daq/slowcontrol"}});
  ASSERT_FALSE(tmBad.ok());
  EXPECT_NE(tmBad.error().message.find("test_mode"), std::string::npos);
  // Un directorio enlazado al real tampoco pasa.
  const std::string base = ::testing::TempDir() + "daq_simlink";
  ::mkdir(base.c_str(), 0700);
  ::mkdir((base + "/real").c_str(), 0700);
  ::unlink((base + "/link").c_str());
  ASSERT_EQ(::symlink((base + "/real").c_str(), (base + "/link").c_str()), 0);
  EXPECT_TRUE(pathIsWithin(base + "/link/2026/hv.csv", base + "/real"));
  EXPECT_TRUE(pathIsWithin(base + "/real", base + "/link"));
  EXPECT_FALSE(pathIsWithin(base + "/realx", base + "/real"));
  EXPECT_FALSE(pathIsWithin(base + "/other", base + "/real"));
  EXPECT_EQ(resolvePath("/a/./b/../c//d/"), "/a/c/d");
}
