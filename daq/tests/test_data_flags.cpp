// Tests de los bits de datos (DAQ2-9): fuga de luz por cara (bit 7), apagón inferido por grupo
// (bit 6), estadística insuficiente (bit 11), aviso de IMon, base explícita, correos por episodio y
// bit 8 (estabilización de un auto-retry). Períodos sintéticos de 60 s con cuentas de Poisson.
#include <gtest/gtest.h>

#include <random>
#include <set>

#include "daq/daqrun/data_flags.h"
#include "daq/core/clock.h"
#include "daq/core/spawn.h"
#include "daq/hvmon/notifier.h"
#include "daq/hvmon/state.h"

using namespace daq;
using namespace daq::daqrun;

namespace {

constexpr UtcNs kS = 1000000000LL;
constexpr UtcNs kT0 = 1792444430LL * kS;
// s por canal y período en los tests de lógica (cuentas holgadas); el real es 44 eventos × 135
// muestras × 16 ns = 9.5e-5 s (Gen::exposure).
constexpr double kExposure = 0.0948;

hvmon::Mapping mapping() {
  auto m = hvmon::Mapping::fromFile(DAQ_SOURCE_DIR "/config/mapping/hodo_mapping_2026_10.json");
  EXPECT_TRUE(m.ok());
  return m.value();
}

struct Gen {
  std::mt19937_64 rng{11};
  // Tasas (Hz por canal) por canal para las cuentas del bit 6 y del bit 7.
  std::array<double, 64> blk{}, leak{};
  double exposure = kExposure;
  MonitorSummary period(int k) {
    MonitorSummary s;
    s.tStartUtcNs = kT0 + k * 60 * kS;
    s.tEndUtcNs = s.tStartUtcNs + 60 * kS;
    s.nEvents = 44;
    for (std::size_t c = 0; c < 64; ++c) {
      s.ootTimeS[c] = exposure;
      s.blackoutCounts[c] = static_cast<std::uint32_t>(std::poisson_distribution<int>(std::max(1e-9, blk[c] * exposure))(rng));
      s.leakCounts[c] = static_cast<std::uint32_t>(std::poisson_distribution<int>(std::max(1e-9, leak[c] * exposure))(rng));
      s.lowCounts[c] = s.blackoutCounts[c] + 1;
    }
    return s;
  }
};

HvSample hvAll(double imonA, double i0setA = 15e-6, int h0 = -1, double imon0 = 0) {
  HvSample s;
  s.fileOk = true;
  s.connected = true;
  s.alive = true;
  s.stale = false;
  for (int h = 0; h < 12; ++h) {
    auto& c = s.ch[static_cast<std::size_t>(h)];
    c.present = true;
    c.pw = h < 8 ? 1 : 0;
    c.imonA = h == h0 ? imon0 : imonA;
    c.i0setA = i0setA;
  }
  return s;
}

const UnitStatus* unit(const DataFlagEngine& e, const std::string& key) {
  for (const auto& u : e.view().units)
    if (u.key == key) return &u;
  return nullptr;
}

std::vector<int> channelsOfFace(const hvmon::Mapping& m, const std::string& face) {
  std::vector<int> v;
  for (const auto& e : m.channels)
    if (e.station + " " + e.face == face) v.push_back(e.digCh);
  return v;
}

}  // namespace

// Cara Trasera de hod2019 a 6 kHz por canal; a la hora sube a 3× durante 20 min: bit 7 en sus 8
// canales, un correo light_leak al empezar, fin tras 3 evaluaciones sin condición. Antes de tener
// base (40 min), bit 11 y bit 6 de evento (provisional). Ningún falso positivo en las otras caras.
TEST(DataFlags, LightLeakOnFaceRaisesBit7OnceAndEnds) {
  const auto m = mapping();
  DataFlagConfig cfg;
  DataFlagEngine eng(cfg, m, "test");
  Gen g;
  for (int c = 32; c < 64; ++c) g.leak[static_cast<std::size_t>(c)] = 6000;
  const auto tras = channelsOfFace(m, "hod2019 Trasera");
  ASSERT_EQ(tras.size(), 8u);
  int mails = 0, firstActive = -1, lastActive = -1, otherFaceActive = 0;
  const HvSample hv = hvAll(4.5e-6);
  for (int k = 0; k < 150; ++k) {
    const bool leak = k >= 60 && k < 80;
    for (int c : tras) g.leak[static_cast<std::size_t>(c)] = leak ? 18000 : 6000;
    const auto ev = eng.addPeriod(g.period(k), &hv);
    for (const auto& e : ev) {
      EXPECT_EQ(e.kind, "light_leak");
      EXPECT_EQ(e.context["unit"], "face:hod2019 Trasera");
      ++mails;
    }
    if (k == 5) {
      EXPECT_TRUE(eng.bits().provisional);
      EXPECT_NE(eng.bits().ch[static_cast<std::size_t>(tras[0])] & chflag::StatInsufficient, 0);
    }
    const bool active = eng.bits().ch[static_cast<std::size_t>(tras[0])] & chflag::LightLeak;
    for (int c : tras) EXPECT_EQ((eng.bits().ch[static_cast<std::size_t>(c)] & chflag::LightLeak) != 0, active);
    if (active && firstActive < 0) firstActive = k;
    if (active) lastActive = k;
    for (const char* f : {"face:hod2019 Delantera", "face:hod2019 Frontal", "face:hod2019 Opuesta"})
      otherFaceActive += unit(eng, f)->active;
  }
  EXPECT_EQ(mails, 1);
  EXPECT_GE(firstActive, 60);
  EXPECT_LE(firstActive, 66);  // con 3× la ventana de 10 min pasa la razón 2 tras ~6 min
  EXPECT_GE(lastActive, 79);
  EXPECT_LE(lastActive, 92);  // la ventana se vacía de luz 10 min después y luego 3 evaluaciones
  EXPECT_EQ(otherFaceActive, 0);
  EXPECT_FALSE(eng.bits().provisional);
}

// Grupo 3 (canales 33, 38, 49) a 14 kHz sobre 1.5 mV; a los 60 min callan: bit 6 en esos canales,
// un correo data_blackout con el canal HV CH02. Un grupo con μ bajo el mínimo no decide (bit 11).
TEST(DataFlags, BlackoutOfGroupRaisesBit6AndLowStatisticsGivesBit11) {
  const auto m = mapping();
  DataFlagConfig cfg;
  DataFlagEngine eng(cfg, m, "test");
  Gen g;
  for (int c = 0; c < 32; ++c) g.blk[static_cast<std::size_t>(c)] = 1.6e6;
  for (int c = 32; c < 64; ++c) g.blk[static_cast<std::size_t>(c)] = 14000;
  for (int c : {50, 61}) g.blk[static_cast<std::size_t>(c)] = 0.5;  // grupo 4: casi sin cuentas
  const HvSample hv = hvAll(4.5e-6);
  int mails = 0, first = -1, falsePos = 0;
  for (int k = 0; k < 90; ++k) {
    for (int c : {33, 38, 49}) g.blk[static_cast<std::size_t>(c)] = k >= 60 ? 0.0 : 14000;
    for (const auto& e : eng.addPeriod(g.period(k), &hv)) {
      EXPECT_EQ(e.kind, "data_blackout");
      EXPECT_EQ(e.ch, 2);
      ++mails;
    }
    const bool on = eng.bits().ch[33] & chflag::DataBlackout;
    if (on && first < 0) first = k;
    for (int c = 0; c < 64; ++c)
      if (c != 33 && c != 38 && c != 49 && (eng.bits().ch[static_cast<std::size_t>(c)] & chflag::DataBlackout)) ++falsePos;
  }
  EXPECT_EQ(mails, 1);
  EXPECT_GE(first, 60);
  EXPECT_LE(first, 65);  // ventana de 5 min: ~4-5 min de silencio para bajar de 0.1
  EXPECT_EQ(falsePos, 0);
  // Grupo 4 (50, 61 de hod2019): μ < 20 en 5 min, no se decide.
  EXPECT_FALSE(unit(eng, "group:4")->decided);
  EXPECT_NE(eng.bits().ch[50] & chflag::StatInsufficient, 0);
}

// Base explícita: decide desde la primera ventana completa, sin esperar historia.
TEST(DataFlags, ExplicitBaselineDecidesWithoutHistory) {
  const auto m = mapping();
  DataFlagConfig cfg;
  cfg.baselines["face:hod2019 Trasera"] = 6000;
  DataFlagEngine eng(cfg, m, "test");
  Gen g;
  for (int c = 32; c < 64; ++c) g.leak[static_cast<std::size_t>(c)] = 6000;
  for (int c : channelsOfFace(m, "hod2019 Trasera")) g.leak[static_cast<std::size_t>(c)] = 15000;
  const HvSample hv = hvAll(4.5e-6);
  int mails = 0;
  for (int k = 0; k < 12; ++k) mails += static_cast<int>(eng.addPeriod(g.period(k), &hv).size());
  EXPECT_EQ(mails, 1);
  const auto* u = unit(eng, "face:hod2019 Trasera");
  EXPECT_TRUE(u->decided && u->active);
  EXPECT_DOUBLE_EQ(u->baselineHz, 6000);
  EXPECT_FALSE(unit(eng, "face:hod2019 Delantera")->decided);  // sin base todavía
}

// IMon: CH01 sube 1.6× su base (aviso por k = 1.5); CH05 pasa 60 % de I0Set sin base.
TEST(DataFlags, ImonWarningAboveKTimesBaselineOrSixtyPercentOfI0Set) {
  const auto m = mapping();
  DataFlagConfig cfg;
  DataFlagEngine eng(cfg, m, "test");
  Gen g;
  std::vector<std::string> warned;
  for (int k = 0; k < 80; ++k) {
    HvSample hv = hvAll(4.0e-6);
    hv.ch[1].imonA = k >= 60 ? 6.4e-6 : 4.0e-6;
    hv.ch[5].imonA = k >= 3 ? 9.5e-6 : 4.0e-6;  // > 9 µA = 60 % de 15 µA
    for (const auto& e : eng.addPeriod(g.period(k), &hv)) {
      EXPECT_EQ(e.kind, "imon_warning");
      warned.push_back(e.context["unit"].get<std::string>() + "@" + std::to_string(k));
    }
  }
  ASSERT_EQ(warned.size(), 2u);
  EXPECT_EQ(warned[0].substr(0, 10), "imon:CH05@");
  EXPECT_EQ(warned[1].substr(0, 10), "imon:CH01@");
  const int kCh01 = std::stoi(warned[1].substr(10));
  EXPECT_GE(kCh01, 60);
  EXPECT_LE(kCh01, 62);  // ventana de 2 min
}

// Notificador: un correo por unidad y episodio para los tipos de daqrun.
TEST(DataFlags, NotifierSendsOneMailPerUnitAndEpisode) {
  struct NullSpawner final : Spawner {
    int calls = 0;
    Result<SpawnResult> run(const std::vector<std::string>&, const std::string&, std::chrono::seconds) override {
      ++calls;
      return SpawnResult{};
    }
  } sp;
  ManualClock clock(kT0);
  hvmon::StateHub hub;
  auto cfg = hvmon::NotifierConfig::fromJson({{"enabled", true}, {"outbox_dir", ::testing::TempDir() + "dataflags_outbox"}}).value();
  hvmon::SmtpSettings smtp;
  smtp.url = "smtps://x:465";
  smtp.user = "u";
  smtp.password = "p";
  smtp.from = "a@b";
  smtp.to = {"c@d"};
  hvmon::Notifier n(cfg, smtp, std::nullopt, sp, clock, ::testing::TempDir(), hub, "h");
  auto ev = [](const std::string& unit, std::uint64_t ep) {
    hvmon::NotifyEvent e;
    e.kind = "light_leak";
    e.ch = -1;
    e.episode = ep;
    e.context = {{"unit", unit}, {"unit_name", unit}};
    return e;
  };
  EXPECT_TRUE(n.submit(ev("face:hod2019 Trasera", 1)));
  EXPECT_FALSE(n.submit(ev("face:hod2019 Trasera", 1)));  // mismo episodio
  EXPECT_TRUE(n.submit(ev("face:hod2019 Frontal", 1)));   // otra cara
  EXPECT_TRUE(n.submit(ev("face:hod2019 Trasera", 2)));   // episodio nuevo
  n.processDue();
  EXPECT_EQ(sp.calls, 3);
}

// Bit 8: canal HV en STABILIZING(k) de un auto-retry.
TEST(DataFlags, RetryStabilizationSetsBit8FromHvmonState) {
  const auto m = mapping();
  auto cj = DaqrunConfig::fromJson({{"campaign", {{"name", "t"}}}});
  ASSERT_TRUE(cj.ok());
  FlagEngine fe(cj.value(), m);
  for (int k = 0; k < 10; ++k) {
    HvSample s = hvAll(4e-6);
    s.readSeq = static_cast<std::uint64_t>(k + 1);
    s.tReqUtcNs = kT0 + k * kS;
    s.tRespUtcNs = s.tReqUtcNs + 3000000;
    s.receivedUtcNs = s.tRespUtcNs + 1000000;
    s.stateUtcNs = s.receivedUtcNs;
    for (int h = 0; h < 8; ++h) {
      auto& c = s.ch[static_cast<std::size_t>(h)];
      c.v0setV = c.vmonV = 65.7;
      c.status = 1;  // on
      c.retryState = h == 2 && k >= 5 ? "STABILIZING(1)" : "NORMAL";
    }
    fe.addSample(s);
  }
  const auto f = fe.flagsFor(kT0 + 7 * kS + 500000000);
  for (int c : {33, 38, 49}) EXPECT_NE(f.ch[static_cast<std::size_t>(c)] & chflag::RetryStab, 0) << c;
  EXPECT_EQ(f.ch[32] & chflag::RetryStab, 0);
  const auto f0 = fe.flagsFor(kT0 + 2 * kS + 500000000);
  EXPECT_EQ(f0.ch[33] & chflag::RetryStab, 0);
}

// Períodos reales de ~61 s (cierran con la lectura siguiente): la ventana de 10 min debe quedar
// completa y decidirse (DAQ2-9: con la poda a 601 s nunca se decidía en los runs 7-8).
TEST(DataFlags, SlightlyLongerRealPeriodsStillFillTheWindow) {
  const auto m = mapping();
  DataFlagConfig cfg;
  DataFlagEngine eng(cfg, m, "test");
  Gen g;
  for (int c = 32; c < 64; ++c) g.leak[static_cast<std::size_t>(c)] = 6000;
  const HvSample hv = hvAll(4.5e-6);
  for (int k = 0; k < 60; ++k) {
    MonitorSummary s = g.period(0);
    s.tStartUtcNs = kT0 + k * 61 * kS;
    s.tEndUtcNs = s.tStartUtcNs + 61 * kS;
    eng.addPeriod(s, &hv);
  }
  EXPECT_TRUE(unit(eng, "face:hod2019 Trasera")->decided);
  EXPECT_GT(unit(eng, "face:hod2019 Trasera")->baselineHz, 0);
}

// Base baja como la de la tarde (1.4 kHz por canal en la cara Trasera): μ ≈ 7 en 10 min, menor que
// el mínimo de 20, pero la luz da n ≥ 20 y la base tiene cientos de cuentas: el bit 7 se decide y
// se marca (runs 7-8, DAQ2-9). Con la base de la noche (6 kHz) nada cambia.
TEST(DataFlags, LowBaselineExcessDecidesWhenObservedCountsAreLarge) {
  const auto m = mapping();
  DataFlagConfig cfg;
  DataFlagEngine eng(cfg, m, "test");
  Gen g;
  g.exposure = 9.5e-5;  // real: μ ≈ 10 en 10 min con 1.4 kHz por canal
  const auto tras = channelsOfFace(m, "hod2019 Trasera");
  const HvSample hv = hvAll(4.5e-6);
  int first = -1;
  for (int k = 0; k < 120; ++k) {
    for (int c : tras) g.leak[static_cast<std::size_t>(c)] = k >= 90 ? 7000 : 1400;
    eng.addPeriod(g.period(k), &hv);
    const auto* u = unit(eng, "face:hod2019 Trasera");
    if (k == 80) {
      EXPECT_LT(u->mu, 20);
      EXPECT_FALSE(u->decided);  // n < 20 sin luz: no se decide
    }
    if (u->active && first < 0) first = k;
  }
  EXPECT_GE(first, 90);
  EXPECT_LE(first, 94);
}

// Mapeo r3 (DAQ2-11): coordenada medida de hod2018 X (orden de mayo) y la de la hoja; la unidad de
// luz empírica reemplaza a las caras de hod2019 en el bit 7.
TEST(MappingR3, MeasuredOrderSheetOrderAndLightUnit) {
  auto m = hvmon::Mapping::fromFile(DAQ_SOURCE_DIR "/config/mapping/hodo_mapping_2026_10_r3.json");
  ASSERT_TRUE(m.ok()) << m.error().message;
  const auto& r3 = m.value();
  EXPECT_EQ(r3.version, "2026_10-r3");
  ASSERT_EQ(r3.lightUnits.size(), 1u);
  EXPECT_EQ(r3.lightUnits.begin()->second, (std::vector<int>{33, 35, 36, 37, 39, 43, 47, 49, 53, 55, 57, 61, 63}));
  const std::vector<int> orderX{15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14};
  for (const auto& e : r3.channels) {
    if (e.station == "hod2018" && e.plane == "X") {
      EXPECT_EQ(orderX[static_cast<std::size_t>(e.strip)], e.digCh);
      EXPECT_EQ(e.stripSheet, e.digCh);
      EXPECT_EQ(e.label, "Canal " + std::to_string(e.digCh));
    } else {
      EXPECT_EQ(e.strip, e.stripSheet) << e.digCh;
    }
  }
  DataFlagConfig cfg;
  DataFlagEngine eng(cfg, r3, "t");
  Gen g;
  const auto hv = hvAll(4.5e-6);
  eng.addPeriod(g.period(0), &hv);
  int light = 0, faces = 0;
  for (const auto& u : eng.view().units) {
    light += u.key.rfind("light:", 0) == 0;
    faces += u.key.rfind("face:", 0) == 0;
  }
  EXPECT_EQ(light, 1);
  EXPECT_EQ(faces, 0);
}
