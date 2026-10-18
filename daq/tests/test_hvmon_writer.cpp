// Tests del escritor de hvmon: rotación a medianoche UTC, esquemas, período de fsync, errores.
#include <gtest/gtest.h>

#include <chrono>

#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hvmon/writer.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::hvmon;

namespace {

// 2026-10-18T23:59:58Z
constexpr UtcNs kNearMidnight = (1792281600LL + 86398LL) * 1000000000LL;

hv::ChannelReading reading(ManualClock& clock) {
  hv::FakeHVBackend be(hv::FakeHVBackend::defaultScenario(), clock);
  EXPECT_TRUE(be.connect().ok());
  return be.readChannels(3, {0, 1}, {"V0Set", "I0Set", "VMon", "IMon", "Status", "Pw"}).value();
}

struct WRig {
  ManualClock clock{kNearMidnight};
  MemoryFileSystem fs;
  StateHub hub;
  SlowControlWriter w{"/sc", 10.0, fs, clock, hub};
};

}  // namespace

TEST(HvmonWriter, DailyRotationAtUtcMidnightWithHeadersAndSchemas) {
  WRig r;
  auto rd = reading(r.clock);  // t = 23:59:58Z
  r.w.handle(FastRows{rd, 1});
  r.clock.advance(3s);  // 00:00:01Z del 19
  auto rd2 = reading(r.clock);
  r.w.handle(FastRows{rd2, 2});
  const auto files = r.fs.files();
  ASSERT_TRUE(files.count("/sc/2026/hv_20261018.csv"));
  ASSERT_TRUE(files.count("/sc/2026/hv_20261019.csv"));
  ASSERT_TRUE(files.count("/sc/2026/hv_20261018.csv.schema.json"));
  ASSERT_TRUE(files.count("/sc/2026/hv_20261019.csv.schema.json"));
  const std::string day1 = r.fs.content("/sc/2026/hv_20261018.csv");
  const std::string header =
      "t_req_utc_ns,t_resp_utc_ns,read_seq,slot,ch,V0Set_V,I0Set_A,VMon_V,IMon_A,Status,Pw,read_duration_ns\n";
  EXPECT_EQ(day1.substr(0, header.size()), header);
  EXPECT_EQ(std::count(day1.begin(), day1.end(), '\n'), 3);  // header + 2 canales
  EXPECT_EQ(r.fs.content("/sc/2026/hv_20261019.csv").substr(0, header.size()), header);
  // Al rotar se hace fsync del archivo del día anterior.
  EXPECT_EQ(r.fs.syncCount("/sc/2026/hv_20261018.csv"), 1u);
  const auto schema = nlohmann::json::parse(r.fs.content("/sc/2026/hv_20261018.csv.schema.json"));
  EXPECT_EQ(schema["schema"], "daqv2.hv_fast/1");
  EXPECT_EQ(schema["file"], "hv_20261018.csv");
  bool imon = false;
  for (const auto& c : schema["columns"])
    if (c["name"] == "IMon_A") imon = c["unit"] == "A" && c["type"] == "float64";
  EXPECT_TRUE(imon);
}

TEST(HvmonWriter, FsyncPeriodIsRespected) {
  WRig r;
  r.clock.advance(5s);  // lejos de la medianoche para no rotar: 00:00:03Z del 19
  auto rd = reading(r.clock);
  for (int i = 0; i < 9; ++i) {
    r.w.handle(FastRows{rd, static_cast<std::uint64_t>(i)});
    r.clock.advance(1s);
    r.w.maybeSync();
  }
  // El escritor se creó en t0 = kNearMidnight; datos de t0 + 5 a t0 + 13 s, maybeSync cada segundo
  // hasta t0 + 14 s: un único fsync, en t0 + 10 s.
  EXPECT_EQ(r.fs.totalSyncs(), 1u);
  for (int i = 0; i < 5; ++i) {
    r.clock.advance(1s);
    r.w.maybeSync();  // sin datos nuevos hasta t0 + 20 s
  }
  r.w.handle(FastRows{rd, 99});
  r.clock.advance(2s);
  r.w.maybeSync();  // t0 + 21 s: vence y hay datos nuevos
  EXPECT_EQ(r.fs.totalSyncs(), 2u);
  r.clock.advance(15s);
  r.w.maybeSync();  // venció pero sin datos nuevos: no hay fsync
  EXPECT_EQ(r.fs.totalSyncs(), 2u);
  EXPECT_EQ(r.w.health().fsyncs, 2u);
}

TEST(HvmonWriter, WriteErrorIsPublishedAndRecovers) {
  WRig r;
  r.clock.advance(5s);
  auto rd = reading(r.clock);
  r.w.handle(FastRows{rd, 1});
  r.fs.failAppendsTo("/sc/2026/hv_20261019.csv", true);
  r.w.handle(FastRows{rd, 2});
  auto h = r.hub.writer.read().first;
  EXPECT_FALSE(h.ok);
  EXPECT_EQ(h.errors, 2u);  // una escritura fallida por fila (2 canales)
  EXPECT_NE(h.lastError.find("No space left"), std::string::npos);
  r.fs.failAppendsTo("/sc/2026/hv_20261019.csv", false);
  r.w.handle(FastRows{rd, 3});
  h = r.hub.writer.read().first;
  EXPECT_TRUE(h.ok);
  EXPECT_EQ(h.errors, 2u);  // los errores quedan registrados
  const auto st = r.hub.compose(HvmonConfig::fromJson(nlohmann::json::object()).value(), r.clock.now_utc_ns(), 1, false);
  EXPECT_EQ(st["writer"]["errors"], 2);
  EXPECT_NE(st["writer"]["last_error"].get<std::string>().find("No space left"), std::string::npos);
}

TEST(HvmonWriter, ChangesNotesTempAndOffsetsFiles) {
  WRig r;
  r.clock.advance(5s);
  const UtcNs t = r.clock.now_utc_ns();
  ChangeRecord c;
  c.id = 7;
  c.tUtcNs = t;
  c.ch = 2;
  c.param = "Pw";
  c.oldValue = 1;
  c.newValue = 0;
  c.origin = "trip";
  r.w.handle(ChangeItem{c});
  CommandRecord k;
  k.id = 3;
  k.tUtcNs = t;
  k.phase = "rejected";
  k.reason = "prueba";
  r.w.handle(CommandItem{k});
  r.w.handle(NoteItem{t, "xplus", "host", "abrí GECO"});
  r.w.handle(TempItem{t, {{"fake-1", 21.5, t}}});
  r.w.handle(OffsetItem{{{"type", "offsets"}, {"n_reads", 61}}});
  const std::string ch = r.fs.content("/sc/2026/hv_changes_20261019.jsonl");
  const auto lines = std::count(ch.begin(), ch.end(), '\n');
  EXPECT_EQ(lines, 2);
  const auto first = nlohmann::json::parse(ch.substr(0, ch.find('\n')));
  EXPECT_EQ(first["type"], "change");
  EXPECT_EQ(first["origin"], "trip");
  EXPECT_NE(ch.find("\"phase\":\"rejected\""), std::string::npos);
  EXPECT_TRUE(r.fs.exists("/sc/2026/hv_changes_20261019.jsonl.schema.json"));
  EXPECT_NE(r.fs.content("/sc/2026/notes_20261019.jsonl").find("abrí GECO"), std::string::npos);
  EXPECT_NE(r.fs.content("/sc/2026/temp_20261019.csv").find("fake-1,21.5"), std::string::npos);
  EXPECT_NE(r.fs.content("/sc/offsets.jsonl").find("\"n_reads\":61"), std::string::npos);
  EXPECT_EQ(r.fs.syncCount("/sc/offsets.jsonl"), 1u);  // offsets: fsync inmediato
  r.w.closeAll();
  EXPECT_GE(r.fs.syncCount("/sc/2026/hv_changes_20261019.jsonl"), 1u);
}

TEST(HvmonWriter, DateKeyIsUtc) {
  EXPECT_EQ(SlowControlWriter::dateKey(kNearMidnight), "20261018");
  EXPECT_EQ(SlowControlWriter::dateKey(kNearMidnight + 2000000000LL), "20261019");
}
