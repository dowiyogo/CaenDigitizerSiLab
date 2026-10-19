// Tests de daqrun sin ROOT: configuración y reglas de simulación, campañas, cliente HV y flags,
// flags de evento, desenrollado del TTT y notas. Todo en directorios temporales.
#include <gtest/gtest.h>

#include <sys/stat.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <set>

#include "daq/core/file_system.h"
#include "daq/core/gzip.h"
#include "daq/daqrun/campaign.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/notes.h"
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
