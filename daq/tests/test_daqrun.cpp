// Tests de daqrun sin ROOT: configuración y reglas de simulación, campañas, cliente HV y flags,
// flags de evento, desenrollado del TTT y notas. Todo en directorios temporales.
#include <gtest/gtest.h>

#include <sys/stat.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <cmath>
#include <map>
#include <random>
#include <set>

#include "daq/core/file_system.h"
#include "daq/core/gzip.h"
#include "daq/daqrun/campaign.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/dash.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/notes.h"
#include "daq/dig/fake_digitizer.h"
#include "daq/core/spawn.h"
#include "daq/hvmon/mapping.h"
#include "daq/hvmon/notifier.h"
#include "daq/hvmon/state.h"
#include <iterator>

using namespace daq;
using namespace daq::daqrun;

namespace {

constexpr UtcNs kS = 1000000000LL;
constexpr UtcNs kT0 = 1792281600LL * kS;  // 2026-10-18T00:00:00Z

std::string tmpDir(const std::string& tag) {
  const std::string d = ::testing::TempDir() + "daqrun_" + tag + "_" + std::to_string(::getpid());
  PosixFileSystem().makeDirectories(d, 0755);
  return d;
}

DaqrunConfig cfgFor(nlohmann::json extra = nlohmann::json::object()) {
  nlohmann::json j{{"campaign", {{"name", "t1"}}}};
  for (auto it = extra.begin(); it != extra.end(); ++it) j[it.key()] = it.value();
  auto c = DaqrunConfig::fromJson(j);
  EXPECT_TRUE(c.ok()) << (c.ok() ? "" : c.error().message);
  return c.value();
}

hvmon::Mapping mapping() {
  auto m = hvmon::Mapping::fromFile(DAQ_SOURCE_DIR "/config/mapping/hodo_mapping_2026_10.json");
  EXPECT_TRUE(m.ok());
  return m.value();
}

// Estado de hvmon guionado: lectura k a los k s (t_req = k s − 10 ms), CH00-07 encendidos a
// 65.75 V; offCh apagado con internal_trip desde la lectura tripAt.
nlohmann::json state(int k, int offCh = -1, int tripAt = 1 << 30, bool final = false, bool stale = false) {
  nlohmann::json chs = nlohmann::json::array();
  for (int h = 0; h < 12; ++h) {
    const bool on = h < 8 && !(h == offCh && k >= tripAt);
    chs.push_back({{"ch", h},
                   {"values",
                    {{"Pw", on ? 1 : 0},
                     {"Status", on ? 1 : (h == offCh ? 512 : 0)},
                     {"V0Set_V", h < 8 ? 65.75 : 0.0},
                     {"VMon_V", on ? 65.75 : 0.2},
                     {"IMon_A", on ? 2e-6 : -4e-8},
                     {"I0Set_A", 1.5e-5}}},
                   {"auto_retry", {{"label", "DISABLED"}}}});
  }
  return {{"schema", "daqv2.hvmon_state/1"},
          {"seq", k},
          {"t_utc_ns", kT0 + k * kS + 5000000},
          {"pid", 1},
          {"final", final},
          {"connection", {{"connected", true}, {"stale", stale}}},
          {"last_read", {{"read_seq", k}, {"t_req_utc_ns", kT0 + k * kS - 10000000}, {"t_resp_utc_ns", kT0 + k * kS}}},
          {"channels", chs},
          {"changes", nlohmann::json::array()}};
}

HvSample sample(int k, int offCh = -1, int tripAt = 1 << 30, bool alive = true) {
  return parseHvState(state(k, offCh, tripAt), kT0 + k * kS + 20000000, 3.0, alive);
}

}  // namespace

TEST(DaqrunConfig, SimulatedRunsNeverWriteToRealRawOrManifests) {
  // Fake sin rutas: /data/daq/sim/...
  auto c = cfgFor();
  EXPECT_EQ(c.rawRoot, "/data/daq/sim/raw");
  EXPECT_EQ(c.manifestsDir, "/data/daq/sim/manifests");
  EXPECT_TRUE(c.simulated());
  // Fake (o replay) con rutas reales: rechazado, también por "..".
  for (const char* r : {"/data/daq/raw", "/data/daq/raw/hodoscope_extnim_2026", "/data/daq/sim/../raw"}) {
    auto bad = DaqrunConfig::fromJson({{"campaign", {{"name", "t1"}, {"raw_root", r}}}});
    ASSERT_FALSE(bad.ok()) << r;
    EXPECT_NE(bad.error().message.find("simulado"), std::string::npos);
    auto bad2 = DaqrunConfig::fromJson(
        {{"campaign", {{"name", "t1"}, {"raw_root", r}}}, {"digitizer", {{"backend", "replay"}}}});
    EXPECT_FALSE(bad2.ok());
  }
  EXPECT_FALSE(DaqrunConfig::fromJson({{"campaign", {{"name", "t1"}, {"manifests_dir", "/data/daq/manifests"}}}}).ok());
  // Fase 4: backend real con las rutas reales por omisión; la longitud pedida no necesita ser
  // múltiplo de 3 (la placa la ajusta: 625 -> 624).
  auto real = DaqrunConfig::fromJson({{"digitizer", {{"backend", "caen"}, {"record_length_samples", 625}}}});
  ASSERT_TRUE(real.ok()) << real.error().message;
  EXPECT_FALSE(real.value().simulated());
  EXPECT_EQ(real.value().rawRoot, "/data/daq/raw");
  EXPECT_EQ(real.value().manifestsDir, "/data/daq/manifests");
  EXPECT_FALSE(DaqrunConfig::fromJson({{"digitizer", {{"backend", "fake"}, {"record_length_samples", 625}}}}).ok());
  EXPECT_FALSE(DaqrunConfig::fromJson({{"digitizer", {{"backend", "otro"}}}}).ok());
  // Validación de la ventana y del nombre.
  EXPECT_FALSE(DaqrunConfig::fromJson({{"digitizer", {{"window_layout", {{"signal", {30, 58}}}}}}}).ok());
  EXPECT_FALSE(DaqrunConfig::fromJson({{"digitizer", {{"record_length_samples", 145}}}}).ok());
  EXPECT_FALSE(DaqrunConfig::fromJson({{"campaign", {{"name", "../x"}}}}).ok());
  EXPECT_FALSE(DaqrunConfig::fromJson({{"monitoring", {{"http", {{"bind", "0.0.0.0"}}}}}}).ok());
  EXPECT_TRUE(validCampaignName("hodo_2026-10.sim"));
  EXPECT_FALSE(validCampaignName("a/b"));
  EXPECT_FALSE(validCampaignName("a..b"));
  // Ida y vuelta de la configuración.
  auto rt = DaqrunConfig::fromJson(c.toJson());
  ASSERT_TRUE(rt.ok()) << rt.error().message;
  EXPECT_EQ(rt.value().toJson(), c.toJson());
}

TEST(Campaign, OpenRunCloseRules) {
  const std::string base = tmpDir("camp");
  const CampaignLocation loc{base + "/raw", base + "/man", "simc"};
  const CampaignLocation other{base + "/raw", base + "/man", "simd"};
  // Sin campaña abierta no hay run.
  auto none = RunReservation::reserve(loc, kT0);
  ASSERT_FALSE(none.ok());
  EXPECT_NE(none.error().message.find("no hay campaña abierta"), std::string::npos);
  // Una campaña simulada no va en /data/daq/raw (se rechaza antes de tocar nada).
  EXPECT_FALSE(openCampaign(CampaignLocation{"/data/daq/raw", base + "/man", "x"}, "", true, kT0).ok());
  ASSERT_TRUE(openCampaign(loc, "tester", true, kT0).ok());
  ASSERT_TRUE(openCampaign(other, "tester", true, kT0).ok());
  EXPECT_EQ(campaignState(loc), CampaignState::Open);
  EXPECT_FALSE(openCampaign(loc, "tester", true, kT0).ok());  // ya abierta
  auto r1 = RunReservation::reserve(loc, kT0 + kS);
  ASSERT_TRUE(r1.ok()) << r1.error().message;
  EXPECT_EQ(r1.value().run(), 1);
  EXPECT_EQ(r1.value().baseName(), "simc_20261018T000001Z_run001");
  EXPECT_EQ(r1.value().segmentFile(2), "simc_20261018T000001Z_run001_seg002.root");
  // Un segundo run concurrente en la misma campaña: rechazado (lock).
  EXPECT_FALSE(RunReservation::reserve(loc, kT0 + 2 * kS).ok());
  // Cerrar con un run en curso: rechazado.
  EXPECT_FALSE(closeCampaign(loc, kT0 + 3 * kS).ok());
  { std::ofstream(loc.dir() + "/" + r1.value().indexFile()) << "{}"; }
  ASSERT_TRUE(appendSha256Sum(loc.dir(), r1.value().indexFile(), std::string(64, 'a')).ok());
  ASSERT_TRUE(updateManifestRun(loc, {{"run", 1}, {"n_events", 5}}).ok());
  ASSERT_TRUE(updateManifestRun(loc, {{"run", 1}, {"n_events", 7}}).ok());  // reemplaza
  r1 = Result<RunReservation>::failure(0, "liberado");                    // suelta el lock
  auto r2 = RunReservation::reserve(loc, kT0 + 4 * kS);
  ASSERT_TRUE(r2.ok());
  EXPECT_EQ(r2.value().run(), 2);
  // El nombre de otro run ya presente en otra campaña se detecta.
  { std::ofstream(other.dir() + "/dup_seg001.root") << "x"; }
  EXPECT_TRUE(fileNameUsedInAnyCampaign(loc.rawRoot, "dup_seg001.root"));
  EXPECT_FALSE(fileNameUsedInAnyCampaign(loc.rawRoot, "nada.root"));
  r2 = Result<RunReservation>::failure(0, "liberado");
  ASSERT_TRUE(closeCampaign(loc, kT0 + 5 * kS).ok());
  EXPECT_EQ(campaignState(loc), CampaignState::Closed);
  EXPECT_FALSE(PosixFileSystem().exists(loc.dir() + "/CAMPAIGN_OPEN"));
  struct stat st {};
  ASSERT_EQ(::stat(loc.dir().c_str(), &st), 0);
  EXPECT_EQ(st.st_mode & 0777, 0555u);
  // Cerrada: ni runs nuevos ni reapertura.
  auto r3 = RunReservation::reserve(loc, kT0 + 6 * kS);
  ASSERT_FALSE(r3.ok());
  EXPECT_NE(r3.error().message.find("cerrada"), std::string::npos);
  EXPECT_FALSE(openCampaign(loc, "tester", true, kT0).ok());
  const auto status = campaignStatus(loc);
  EXPECT_EQ(status["state"], "closed");
  EXPECT_EQ(status["manifest"]["status"], "closed");
  ASSERT_EQ(status["manifest"]["runs"].size(), 1u);
  EXPECT_EQ(status["manifest"]["runs"][0]["n_events"], 7);
  EXPECT_EQ(status["manifest"]["sha256sums_sha256"].get<std::string>().size(), 64u);
  ::chmod(loc.dir().c_str(), 0755);  // para que el temporal se pueda limpiar
}

TEST(HvClient, ParseStalenessFinalAndDeadPid) {
  auto s = parseHvState(state(10), kT0 + 10 * kS + 100000000, 3.0, true);
  EXPECT_TRUE(s.fileOk);
  EXPECT_TRUE(s.alive);
  EXPECT_FALSE(s.stale);
  EXPECT_EQ(s.readSeq, 10u);
  EXPECT_EQ(s.tReqUtcNs, kT0 + 10 * kS - 10000000);
  EXPECT_EQ(s.ch[2].pw, 1u);
  EXPECT_EQ(s.ch[2].retryState, "DISABLED");
  EXPECT_FALSE(parseHvState(state(10), kT0 + 20 * kS, 3.0, true).alive);         // estado viejo
  EXPECT_FALSE(parseHvState(state(10, -1, 0, true), kT0 + 10 * kS, 3.0, true).alive);  // final
  EXPECT_FALSE(parseHvState(state(10), kT0 + 10 * kS, 3.0, false).alive);         // pid muerto
  EXPECT_TRUE(parseHvState(state(10, -1, 0, false, true), kT0 + 10 * kS, 3.0, true).stale);
  EXPECT_FALSE(parseHvState({{"schema", "otra"}}, kT0, 3.0, true).fileOk);
}

TEST(FlagEngine, GroupOffAmbiguityStalenessAndStaticBits) {
  const auto m = mapping();
  const DaqrunConfig cfg = cfgFor();
  FlagEngine fe(cfg, m);
  // Canales del grupo 3 (CH02) según el mapeo, y canales sin grupo conocido.
  std::set<int> g3, unknown;
  for (const auto& e : m.channels) {
    if (e.hvGroup && *e.hvGroup == 3) g3.insert(e.digCh);
    if (!e.hvGroup) unknown.insert(e.digCh);
  }
  ASSERT_EQ(g3, (std::set<int>{33, 38, 49}));
  EXPECT_EQ(unknown.size(), 23u);  // hod2018: 9 canales con grupo supuesto, 23 sin grupo
  for (int k = 0; k <= 70; ++k) fe.addSample(sample(k, 2, 50));
  auto at = [&](double s) { return fe.flagsFor(kT0 + static_cast<UtcNs>(s * 1e9)); };
  // Antes del trip: sin bits 0-3 ni 5 en los canales con grupo; bit 4 en los sin grupo; bit 9 en todos.
  const auto f30 = at(30.5);
  for (int c = 0; c < 64; ++c) {
    const std::uint16_t f = f30.ch[static_cast<std::size_t>(c)];
    // Solo ch48 tiene grupo y cara confirmados por medición (HV-DIAG-5).
    EXPECT_EQ((f & chflag::MapUnconfirmed) != 0, c != 48) << c;
    EXPECT_EQ((f & chflag::HvUnknown) != 0, unknown.count(c) == 1) << c;
    EXPECT_EQ(f & (chflag::HvOff | chflag::HvOvc | chflag::HvRamp | chflag::HvVdev | chflag::HvAmbiguous), 0) << c;
  }
  EXPECT_EQ(f30.hvSeq, 30u);
  // Trip en la lectura 50: intervalo (t_req(49) − 5 s, t_resp(50)] = (43.99 s, 50 s]: bit 5 y el
  // estado pesimista (bit 0) en el grupo 3; nada en los demás.
  for (double s : {45.0, 49.5}) {
    const auto f = at(s);
    for (int c = 0; c < 64; ++c) {
      const std::uint16_t x = f.ch[static_cast<std::size_t>(c)];
      if (g3.count(c)) {
        EXPECT_TRUE(x & chflag::HvAmbiguous) << s << " ch" << c;
        EXPECT_TRUE(x & chflag::HvOff) << s << " ch" << c;
      } else if (!unknown.count(c)) {
        EXPECT_EQ(x & (chflag::HvAmbiguous | chflag::HvOff), 0) << s << " ch" << c;
      }
    }
    EXPECT_TRUE(f.anyHv);
  }
  EXPECT_EQ(at(43.0).ch[33] & chflag::HvAmbiguous, 0);
  // Después: bit 0 exactamente en el grupo 3, sin ambigüedad.
  const auto f55 = at(55.5);
  for (int c = 0; c < 64; ++c) {
    const std::uint16_t x = f55.ch[static_cast<std::size_t>(c)];
    EXPECT_EQ((x & chflag::HvOff) != 0, g3.count(c) == 1) << c;
    EXPECT_EQ(x & chflag::HvAmbiguous, 0) << c;
  }
  // Sin lectura posterior (evento después de la última lectura): bit 4 en todos.
  const auto f80 = at(80.5);
  for (int c = 0; c < 64; ++c) EXPECT_TRUE(f80.ch[static_cast<std::size_t>(c)] & chflag::HvUnknown);
  // hvmon no disponible: bit 4 en todos aunque haya lecturas alrededor.
  FlagEngine fe2(cfg, m);
  for (int k = 0; k <= 20; ++k) fe2.addSample(sample(k, -1, 1 << 30, k < 10));
  EXPECT_FALSE(fe2.flagsFor(kT0 + 5 * kS + kS / 2).ch[33] & chflag::HvUnknown);
  EXPECT_TRUE(fe2.flagsFor(kT0 + 15 * kS + kS / 2).ch[33] & chflag::HvUnknown);
  // Historial: prune conserva lo necesario.
  fe.prune(kT0 + 60 * kS);
  EXPECT_LT(fe.historySize(), 15u);
  EXPECT_EQ(fe.hvChannelOf(33), 2);
  EXPECT_EQ(fe.hvChannelOf(2), -1);
}

TEST(EventFlagger, CounterJumpBoardFullMismatchGapBurst) {
  EventFlagger ef(1e-6, 1e-6, 1.0);
  std::uint32_t f = 0;
  UtcNs t = kT0;
  for (int i = 0; i < 30; ++i) {
    t += kS;
    f = ef.flags(t, t + 1000000, static_cast<std::uint32_t>(i), false, kS);
    EXPECT_EQ(f, 0u) << i;
  }
  EXPECT_NEAR(ef.rateHz(), 1.0, 1e-6);
  EXPECT_TRUE(ef.flags(t + kS, t + kS, 31, false, kS) & evtflag::HwCounterJump);  // 29 -> 31
  EXPECT_TRUE(ef.flags(t + 2 * kS, t + 2 * kS, 32, true, kS) & evtflag::BoardFull);
  EXPECT_TRUE(ef.flags(t + 3 * kS, t + 5 * kS, 33, false, kS) & evtflag::TttPcMismatch);
  // Hueco: T = −ln(1e−6)/r ≈ 13.8 s a 1 Hz.
  EXPECT_TRUE(ef.flags(t + 63 * kS, t + 63 * kS, 34, false, 60 * kS) & evtflag::Gap);
  // Ráfaga: 25 triggers en 0.25 s.
  std::uint32_t any = 0;
  for (int i = 0; i < 25; ++i) any |= ef.flags(t + 70 * kS + i * 10000000, t + 70 * kS + i * 10000000, 35u + static_cast<std::uint32_t>(i), false, 10000000);
  EXPECT_TRUE(any & evtflag::Burst);
  // Vuelta del contador de 24 bit: 0xFFFFFF -> 0 no es salto.
  EventFlagger e2(1e-6, 1e-6, 1.0);
  e2.flags(kT0, kT0, 0xFFFFFF, false, 0);
  EXPECT_FALSE(e2.flags(kT0 + kS, kT0 + kS, 0, false, kS) & evtflag::HwCounterJump);
}

TEST(TttUnwrapper, UsesPcTimeAcrossWraps) {
  const double tick = 8.0;
  const std::uint64_t wrap = 1ull << 31;  // 17.18 s
  TttUnwrapper u(kT0, tick, 31);
  for (double s : {0.5, 3.0, 16.9, 17.5, 40.0, 41.0, 300.0, 300.00001}) {
    const auto trueNs = static_cast<std::int64_t>(s * 1e9);
    const auto raw = static_cast<std::uint32_t>((static_cast<std::uint64_t>(trueNs / 8)) % wrap);
    const auto got = u.tttNs(raw, kT0 + trueNs + 3000000);  // lectura 3 ms después
    EXPECT_NEAR(static_cast<double>(got), static_cast<double>(trueNs), 8.0) << s;
  }
}

TEST(Notes, IntervalFromPlainAndGzDays) {
  const std::string d = tmpDir("notes");
  PosixFileSystem fs;
  fs.makeDirectories(d + "/2026", 0755);
  auto line = [](UtcNs t, const char* text) {
    return nlohmann::json{{"t_utc_ns", t}, {"user", "rene"}, {"host", "h"}, {"text", text}}.dump() + "\n";
  };
  {
    std::ofstream o(d + "/2026/notes_20261018.jsonl");
    o << line(kT0 + 10 * kS, "a") << line(kT0 + 86000 * kS, "b");
  }
  {
    std::ofstream o(d + "/2026/notes_20261019.jsonl");
    o << line(kT0 + 86410 * kS, "c") << line(kT0 + 90000 * kS, "d");
  }
  ASSERT_TRUE(gzipFileVerified(d + "/2026/notes_20261019.jsonl", d + "/2026/notes_20261019.jsonl.gz", true).ok());
  const auto n = readNotes(d, kT0 + 20 * kS, kT0 + 90000 * kS);
  ASSERT_EQ(n.size(), 2u);
  EXPECT_EQ(n[0].text, "b");
  EXPECT_EQ(n[1].text, "c");
  EXPECT_EQ(n[1].user, "rene");
  EXPECT_TRUE(readNotes(d + "/nada", kT0, kT0 + kS).empty());
}

namespace {
struct NoSpawn : Spawner {
  std::vector<std::string> bodies;
  Result<SpawnResult> run(const std::vector<std::string>&, const std::string& stdinData, std::chrono::seconds) override {
    const auto p = stdinData.find("upload-file = \"");
    const std::string path = stdinData.substr(p + 15, stdinData.find('"', p + 15) - p - 15);
    std::ifstream in(path);
    bodies.emplace_back((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return SpawnResult{};
  }
};
}  // namespace

TEST(HvmonWatchdog, OneMailPerStaleEpisodeAfterFiveMinutesMarkedWithRun) {
  HvmonWatchdog wd(300.0, "daqrun campaña commissioning_2026_10 run 3");
  std::vector<hvmon::NotifyEvent> events;
  // Vivo hasta 100 s; muerto (pid inexistente) de 101 a 1000 s; vivo 1001-1099; viejo desde 1100 s.
  for (int k = 0; k <= 1600; ++k) {
    const bool alive = k <= 100 || (k > 1000 && k < 1100);
    HvSample s;
    if (k < 1100 || k > 1500) s = parseHvState(state(k), kT0 + k * kS + 20000000, 3.0, alive);
    else s = parseHvState(state(1099), kT0 + k * kS, 3.0, true);  // archivo que ya no cambia
    if (auto e = wd.feed(s)) events.push_back(*e);
  }
  ASSERT_EQ(events.size(), 2u);
  EXPECT_EQ(events[0].kind, "hvmon_stale");
  EXPECT_EQ(events[0].tUtcNs, kT0 + 101 * kS + 20000000);  // inicio del episodio
  EXPECT_NE(events[0].reason.find("bit 4"), std::string::npos);
  EXPECT_NE(events[0].reason.find("el proceso de hvmon no existe"), std::string::npos);
  EXPECT_EQ(events[0].origin, "daqrun campaña commissioning_2026_10 run 3");
  EXPECT_NE(events[1].episode, events[0].episode);
  // El correo lleva el run y la campaña en el asunto y el cuerpo; un segundo aviso del mismo tipo
  // dentro de problem_interval_min (30 min) se descarta.
  ManualClock clock(kT0 + 401 * kS);
  hvmon::StateHub hub;
  NoSpawn sp;
  hvmon::SmtpSettings smtp;
  smtp.url = "smtps://x";
  smtp.from = "a@x";
  smtp.to = {"b@x"};
  auto nc = hvmon::NotifierConfig::fromJson({{"enabled", true}, {"outbox_dir", tmpDir("wdbox")}});
  hvmon::Notifier n(nc.value(), smtp, std::nullopt, sp, clock, ::testing::TempDir(), hub, "h");
  EXPECT_TRUE(n.submit(events[0]));
  EXPECT_FALSE(n.submit(events[0]));
  n.processDue();
  ASSERT_EQ(sp.bodies.size(), 1u);
  // El asunto tiene caracteres no ASCII: viaja codificado (RFC 2047); se mira el asunto sin codificar.
  std::string subj;
  n.composeMessage(events[0], &subj);
  EXPECT_EQ(subj, "[DAQ hodoscopio] hvmon_stale — daqrun campaña commissioning_2026_10 run 3");
  EXPECT_NE(sp.bodies[0].find("Subject: =?UTF-8?"), std::string::npos);
  EXPECT_NE(sp.bodies[0].find("Origen: daqrun campaña commissioning_2026_10 run 3"), std::string::npos) << sp.bodies[0];
  EXPECT_NE(sp.bodies[0].find("Retraso de entrega: 300 s"), std::string::npos) << sp.bodies[0];
}

TEST(TriggerWatchdog, NoTriggersAndResumedWithFakeDigitizerGap) {
  // FakeDigitizer a 0.6 Hz con un hueco de triggers de 12 min desde t = 600 s (como el crate NIM
  // apagado por temperatura): un correo "no_triggers" a los 5 min sin triggers y uno
  // "triggers_resumed" con la duración del hueco cuando vuelven.
  ManualClock clock;
  auto sp = dig::SyntheticParams::fromJson(
      {{"seed", 9}, {"trigger_rate_Hz", 0.6}, {"events", {{{"t_s", 600}, {"type", "trigger_gap"}, {"duration_s", 720}}}}});
  ASSERT_TRUE(sp.ok());
  dig::FakeDigitizer fd(sp.value(), clock);
  dig::DigitizerConfig dc;
  ASSERT_TRUE(fd.configure(dc).ok());
  const UtcNs t0 = fd.start().value().utcBeforeStart;
  TriggerWatchdog wd(300.0, "daqrun campaña commissioning_2026_10 run 4", t0);
  std::vector<hvmon::NotifyEvent> ev;
  UtcNs lastBeforeGap = 0;
  for (int s = 0; s < 2400; ++s) {
    clock.advance(std::chrono::seconds(1));
    dig::EventBlock b;
    ASSERT_TRUE(fd.readBlock(b, 100).ok());
    if (!b.events.empty() && s < 600) lastBeforeGap = b.tReadUtcNs;
    if (auto e = wd.onTriggers(b.tReadUtcNs, b.events.size(), 25.0)) ev.push_back(*e);
    if (auto e = wd.tick(clock.now_utc_ns(), 25.0)) ev.push_back(*e);
  }
  ASSERT_EQ(ev.size(), 2u);
  EXPECT_EQ(ev[0].kind, "no_triggers");
  EXPECT_EQ(ev[0].tUtcNs, lastBeforeGap);  // hora del último trigger
  EXPECT_NE(ev[0].reason.find("último trigger a las"), std::string::npos);
  EXPECT_EQ(ev[0].context["board_temperature_degC"], 25.0);
  EXPECT_EQ(ev[0].context["last_trigger_utc_ns"], lastBeforeGap);
  EXPECT_EQ(ev[1].kind, "triggers_resumed");
  const double gap = ev[1].context["gap_s"].get<double>();
  EXPECT_GT(gap, 720.0);   // el hueco del guion más los intervalos de Poisson a ambos lados
  EXPECT_LT(gap, 760.0);
  EXPECT_NE(ev[1].reason.find("hueco más largo"), std::string::npos);
  EXPECT_EQ(ev[1].origin, "daqrun campaña commissioning_2026_10 run 4");
  // Sin triggers desde el inicio del run: el aviso lo dice; un solo aviso por episodio.
  TriggerWatchdog wd2(300.0, "x", t0);
  EXPECT_FALSE(wd2.tick(t0 + 299 * kS, std::nullopt));
  auto e = wd2.tick(t0 + 301 * kS, std::nullopt);
  ASSERT_TRUE(e);
  EXPECT_NE(e->reason.find("ningún trigger desde el inicio"), std::string::npos);
  EXPECT_TRUE(e->context["board_temperature_degC"].is_null());
  EXPECT_FALSE(wd2.tick(t0 + 900 * kS, std::nullopt));
  // El correo lleva el texto, la temperatura de la placa HV y la hora del último trigger.
  ManualClock nclock(lastBeforeGap + 300 * kS);
  hvmon::StateHub hub;
  NoSpawn spn;
  hvmon::SmtpSettings smtp;
  smtp.url = "smtps://x";
  smtp.from = "a@x";
  smtp.to = {"b@x"};
  hvmon::Notifier n(hvmon::NotifierConfig::fromJson({{"enabled", true}, {"outbox_dir", tmpDir("trigbox")}}).value(), smtp,
                    std::nullopt, spn, nclock, ::testing::TempDir(), hub, "h");
  ASSERT_TRUE(n.submit(ev[0]));
  ASSERT_TRUE(n.submit(ev[1]));
  n.processDue();
  ASSERT_EQ(spn.bodies.size(), 2u);
  EXPECT_NE(spn.bodies[0].find("SIN TRIGGERS externos"), std::string::npos);
  EXPECT_NE(spn.bodies[0].find("Temperatura de la placa A7030DP (no del detector): 25.0 degC"), std::string::npos);
  EXPECT_NE(spn.bodies[1].find("VOLVIERON"), std::string::npos);
}

TEST(TriggerWatchdog, SparseTriggersCountAsLowRate) {
  // Como en el run005 de DAQ2-6 (cable de trigger fallando): tras 5 min normales, un trigger suelto cada ~3 min (nunca 5 min de
  // silencio): el criterio de tasa (< 0.1 Hz en 10 min) avisa; vuelve cuando la tasa se recupera.
  const UtcNs t0 = kT0;
  TriggerWatchdog wd(300.0, "run 5", t0, 0.1, 600.0);
  std::vector<hvmon::NotifyEvent> ev;
  for (int s = 1; s <= 3600; ++s) {
    const UtcNs now = t0 + s * kS;
    std::uint64_t n = 0;
    if (s < 300 || s > 2400) n = (s % 2 == 0) ? 1 : 0;  // 0.5 Hz normal
    else if (s % 180 == 0) n = 1;                      // suelto
    if (auto e = wd.onTriggers(now, n, 30.0)) ev.push_back(*e);
    if (auto e = wd.tick(now, 30.0)) ev.push_back(*e);
  }
  ASSERT_EQ(ev.size(), 2u);
  EXPECT_EQ(ev[0].kind, "no_triggers");
  EXPECT_EQ(ev[0].context["condition"], "tasa baja");
  EXPECT_LT(ev[0].context["window_rate_Hz"].get<double>(), 0.1);
  EXPECT_NE(ev[0].reason.find("tasa de triggers externos"), std::string::npos);
  EXPECT_EQ(ev[1].kind, "triggers_resumed");
  EXPECT_GT(ev[1].tUtcNs, t0 + 2400 * kS);
  // Sin criterio de tasa (min 0) los triggers sueltos no avisan: el comportamiento anterior.
  TriggerWatchdog old(300.0, "run 5", t0, 0.0, 600.0);
  int n = 0;
  for (int s = 1; s <= 3600; ++s) {
    const UtcNs now = t0 + s * kS;
    std::uint64_t k = (s < 300 || s > 2400) ? ((s % 2 == 0) ? 1 : 0) : (s % 180 == 0 ? 1 : 0);
    n += old.onTriggers(now, k, std::nullopt).has_value();
    n += old.tick(now, std::nullopt).has_value();
  }
  EXPECT_EQ(n, 0);
}

// Objeto "dash" de daqrun (DAQ2-8): tasas de 1 y 10 min, canales vivos (> 2 LSB en 10 min),
// mediana entre vivos y avisos.
TEST(DaqrunDash, RatesLiveChannelsMedianAndAlerts) {
  using namespace daq::daqrun;
  constexpr daq::UtcNs kS = 1000000000LL, t0 = 1792500000LL * kS;
  MonitorView v;
  v.events = 1000;
  v.lastEventUtcNs = t0 + 590 * kS;
  for (int k = 0; k < 10; ++k) {
    MonitorSummary s;
    s.tStartUtcNs = t0 + k * 60 * kS;
    s.tEndUtcNs = s.tStartUtcNs + 60 * kS;
    s.nEvents = 42 + static_cast<std::uint64_t>(k);
    s.triggerRateHz = static_cast<double>(s.nEvents) / 60.0;
    for (int c = 0; c < 64; ++c) {
      if (c == 3 || c == 38) continue;  // sin pulsos: no vivos
      s.lowCounts[static_cast<std::size_t>(c)] = 5;
      s.ootCounts[static_cast<std::size_t>(c)] = static_cast<std::uint32_t>(c < 32 ? c : 0);
      s.ootTimeS[static_cast<std::size_t>(c)] = 0.001;
    }
    v.recent.push_back(s);
    v.last = s;
  }
  HvSample hv;
  hv.alive = true;
  hv.stale = false;
  hv.stateUtcNs = t0 + 599 * kS;
  DashInputs in;
  in.campaign = "c";
  in.runName = "r";
  in.runStartUtcNs = t0;
  in.nowUtcNs = t0 + 600 * kS;
  in.view = &v;
  in.hv = &hv;
  in.diskFreeBytes = 400e9;
  in.watchTriggers = true;
  const auto j = dashJson(in);
  EXPECT_NEAR(j["rate_1min_Hz"].get<double>(), 51.0 / 60.0, 1e-12);
  EXPECT_NEAR(j["rate_10min_Hz"].get<double>(), 465.0 / 600.0, 1e-12);
  EXPECT_NEAR(j["since_last_trigger_s"].get<double>(), 10.0, 1e-9);
  EXPECT_EQ(j["stations"][0]["live"], 31);
  EXPECT_EQ(j["stations"][1]["live"], 31);
  // hod2018: canales 0-31 sin el 3, tasa = c / 0.001 s por período -> mediana de 31 valores = canal 16.
  EXPECT_NEAR(j["stations"][0]["median_oot_Hz"].get<double>(), 16000.0, 1e-6);
  EXPECT_NEAR(j["stations"][1]["median_oot_Hz"].get<double>(), 0.0, 1e-12);
  EXPECT_TRUE(j["alerts"].empty()) << j["alerts"].dump();
  hv.stale = true;
  in.noTriggerAlert = true;
  in.readout = {{"read_errors", 2}};
  EXPECT_EQ(dashJson(in)["alerts"].size(), 3u);
}

// TttClock (DAQ2-9): tick real 8.000071 ns (run 8), latencias de lectura uniformes de 0.7-7.1 ms
// (percentil 5 ~1 ms y mediana ~3.9 ms en el run 8) con 1 % de cola de 20-60 ms, 0.73 Hz durante
// 16 h, un hueco de 55 s y re-anclaje por segmento cada hora. Tras el primer ajuste la deriva
// (mediana del error de t_utc contra la hora real del trigger, por bin de 10 min) queda bajo 1 ms;
// con el tick fijo de 8 ns crece a ~0.4 s.
TEST(TttClock, FittedTickKeepsUtcWithinOneMillisecondOverSixteenHours) {
  using namespace daq::daqrun;
  constexpr double kTrue = 8.000071;
  const daq::UtcNs origin = 1792444430LL * 1000000000LL;
  TttClock::Params p;
  p.initialTickNs = 8.0;
  TttClock fit(origin, p);
  p.fit = false;
  TttClock fixed(origin, p);
  std::mt19937_64 rng(7);
  std::exponential_distribution<double> gap(0.73);
  std::uniform_real_distribution<double> latCore(0.7e-3, 7.1e-3), latTail(20e-3, 60e-3), u01(0, 1);
  auto lat = [&](std::mt19937_64& r) { return u01(r) < 0.01 ? latTail(r) : latCore(r); };
  double t = 0.5, maxErrFixed = 0, nextSeg = 3600;
  std::size_t n = 0;
  std::map<long, std::vector<double>> binErr;  // bin de 10 min -> errores (s)
  while (t < 16 * 3600.0) {
    t += gap(rng);
    if (n == 20000) t += 55.0;  // hueco de 55 s (más de tres vueltas del TTT)
    if (t > nextSeg) {
      fit.requestReanchor("segmento");
      nextSeg += 3600;
    }
    const double ticks = t * 1e9 / kTrue;
    const auto raw = static_cast<std::uint32_t>(static_cast<std::uint64_t>(ticks) & 0x7FFFFFFFu) | (t > 17.2 ? 0x80000000u : 0u);
    const daq::UtcNs tRead = origin + static_cast<daq::UtcNs>((t + lat(rng)) * 1e9);
    const daq::UtcNs truth = origin + static_cast<daq::UtcNs>(t * 1e9);
    const auto a = fit.map(raw, tRead);
    const auto b = fixed.map(raw, tRead);
    if (t > 600) binErr[static_cast<long>(t / 600)].push_back(static_cast<double>(a.second - truth) * 1e-9);  // tras el primer ajuste
    maxErrFixed = std::max(maxErrFixed, std::fabs(static_cast<double>(b.second - truth)) * 1e-9);
    ++n;
  }
  double worst = 0;
  long worstBin = -1;
  for (auto& [b, v] : binErr) {
    std::sort(v.begin(), v.end());
    if (std::fabs(v[v.size() / 2]) > worst) {
      worst = std::fabs(v[v.size() / 2]);
      worstBin = b;
    }
  }
  EXPECT_LT(worst, 1e-3) << worst << " en el bin " << worstBin;
  EXPECT_GT(maxErrFixed, 0.4) << maxErrFixed;
  EXPECT_NEAR(fit.tickNs(), kTrue, 5e-5);  // ventana de 300 s
  // Anclas: inicio, primer ajuste (salto: corrige el tick provisional), luego saltos del orden del
  // ruido del ajuste (< 2 ms); una por segmento.
  const auto& an = fit.anchors();
  ASSERT_GE(an.size(), 900u);
  EXPECT_EQ(an[0].reason, "inicio");
  EXPECT_NE(an[1].jumpNs, 0);
  for (std::size_t i = 2; i < an.size(); ++i) EXPECT_LT(std::llabs(an[i].jumpNs), 2000000) << i;
  std::size_t seg = 0;
  for (const auto& x : an) seg += x.reason == "segmento";
  EXPECT_EQ(seg, 16u);  // t = 1 h, 2 h, ..., 16 h
  EXPECT_EQ(fixed.anchors().size(), 1u);
}

// La razón TTT/PC varía unos ppm en 10-20 min (run 8, DAQ2-9): con un tick que oscila ±2.5 ppm en
// 30 min la deriva de t_utc (mediana por 10 min) sigue bajo 1 ms con la ventana de 300 s y
// re-anclaje cada 60 s.
TEST(TttClock, FollowsTickWanderOfAFewPpm) {
  using namespace daq::daqrun;
  const daq::UtcNs origin = 1792444430LL * 1000000000LL;
  TttClock fit(origin, TttClock::Params{});
  std::mt19937_64 rng(5);
  std::exponential_distribution<double> gap(0.73);
  std::uniform_real_distribution<double> lat(1e-3, 6e-3);
  double t = 0.5, ticks = 0.5e9 / 8.000071, prev = 0;
  std::map<long, std::vector<double>> binErr;
  while (t < 4 * 3600.0) {
    const double dt = gap(rng);
    t += dt;
    const double tick = 8.000071 * (1 + 2.5e-6 * std::sin(2 * M_PI * t / 1800));
    ticks += dt * 1e9 / tick;
    const auto raw = static_cast<std::uint32_t>(static_cast<std::uint64_t>(ticks) & 0x7FFFFFFFu);
    const auto a = fit.map(raw, origin + static_cast<daq::UtcNs>((t + lat(rng)) * 1e9));
    if (t > 600) binErr[static_cast<long>(t / 600)].push_back(static_cast<double>(a.second - origin) * 1e-9 - t);
    prev = t;
  }
  (void)prev;
  double worst = 0;
  for (auto& [b, v] : binErr) {
    std::sort(v.begin(), v.end());
    worst = std::max(worst, std::fabs(v[v.size() / 2]));
  }
  EXPECT_LT(worst, 1e-3) << worst;
}

// Episodio por tasa baja como el del 2026-10-20 (DAQ2-9): triggers a 0.7 Hz, silencio desde 600 s,
// aviso por tasa < 0.1 Hz, vuelven a 2340 s. El correo de vuelta mide desde el último trigger
// anterior al aviso (no desde now − ventana) y da el hueco más largo con su inicio y fin.
TEST(TriggerWatchdog, ResumedReportsEpisodeFromLastTriggerAndLongestGap) {
  using namespace daq::daqrun;
  constexpr daq::UtcNs kS = 1000000000LL, t0 = 1792514830LL * kS;
  TriggerWatchdog wd(300.0, "x", t0, 0.1, 600.0);
  std::vector<hvmon::NotifyEvent> ev;
  daq::UtcNs lastBefore = 0;
  for (double t = 1.0; t < 3000; t += 1.0) {
    const daq::UtcNs now = t0 + static_cast<daq::UtcNs>(t * kS);
    const bool trig = (t < 600 && std::fmod(t, 1.4) < 1.0) || (t >= 2340 && std::fmod(t, 1.4) < 1.0);
    if (trig) {
      if (t < 600) lastBefore = now;
      if (auto e = wd.onTriggers(now, 1, std::nullopt)) ev.push_back(*e);
    }
    if (auto e = wd.tick(now, std::nullopt)) ev.push_back(*e);
  }
  ASSERT_EQ(ev.size(), 2u);
  EXPECT_EQ(ev[0].kind, "no_triggers");
  EXPECT_EQ(ev[1].kind, "triggers_resumed");
  EXPECT_EQ(ev[1].context["episode_start_utc_ns"].get<daq::UtcNs>(), lastBefore);
  EXPECT_EQ(ev[1].context["longest_gap_from_utc_ns"].get<daq::UtcNs>(), lastBefore);
  EXPECT_NEAR(ev[1].context["gap_s"].get<double>(), 2340.0 - static_cast<double>(lastBefore - t0) / 1e9, 1.0);
  EXPECT_GE(ev[1].context["episode_s"].get<double>(), ev[1].context["gap_s"].get<double>());
}

// Corte de trigger poco después de empezar (run de prueba de la fase 5, DAQ2-9: 31 eventos antes de
// cortarse a los 570 s, vuelven a los 2340 s). La latencia de referencia sale de las lecturas cercanas
// al origen; tras la vuelta la deriva sigue bajo 1 ms (antes t_utc quedaba ~9 ms adelantado).
TEST(TttClock, TriggerCutBeforeTheFirstFitKeepsTheReference) {
  using namespace daq::daqrun;
  constexpr double kTrue = 8.000095;
  const daq::UtcNs origin = 1792514830LL * 1000000000LL;
  TttClock::Params p;
  p.initialTickNs = 8.000071;
  TttClock fit(origin, p);
  std::mt19937_64 rng(9);
  std::exponential_distribution<double> gap(0.73);
  std::uniform_real_distribution<double> lat(0.7e-3, 7.1e-3);
  std::map<long, std::vector<double>> binErr;
  double t = 0.5;
  while (t < 3600.0) {
    t += gap(rng);
    if (t > 570 && t < 2340) t = 2340 + gap(rng);
    const double ticks = t * 1e9 / kTrue;
    const auto raw = static_cast<std::uint32_t>(static_cast<std::uint64_t>(ticks) & 0x7FFFFFFFu);
    const auto a = fit.map(raw, origin + static_cast<daq::UtcNs>((t + lat(rng)) * 1e9));
    if (t > 2340 + 300) binErr[static_cast<long>(t / 600)].push_back(static_cast<double>(a.second - origin) * 1e-9 - t);
  }
  ASSERT_FALSE(binErr.empty());
  double worst = 0;
  for (auto& [b, v] : binErr) {
    std::sort(v.begin(), v.end());
    worst = std::max(worst, std::fabs(v[v.size() / 2]));
  }
  EXPECT_LT(worst, 1e-3) << worst;
}
