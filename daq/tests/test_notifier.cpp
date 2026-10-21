// Tests del notificador: argv sin secretos, contenido, límites por episodio, reintentos, logs.
// El Spawner es un fake que registra argv y stdin y lee el archivo del cuerpo antes de que se borre.
#include <gtest/gtest.h>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <fstream>
#include <iterator>

#include "daq/core/clock.h"
#include "daq/core/logger.h"
#include "daq/core/time.h"
#include "daq/hvmon/notifier.h"
#include "daq/hvmon/state.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::hvmon;

namespace {

const std::string kPassword = "S3cr3t-P4ss\"w\\rd";  // ficticia, con comilla y barra para el escape

struct FakeSpawner final : Spawner {
  struct Call {
    std::vector<std::string> argv;
    std::string stdinData;
    std::string body;
  };
  std::vector<Call> calls;
  std::vector<int> exitCodes;  // uno por llamada (si se acaban, 0)
  std::string stderrText;
  Result<SpawnResult> run(const std::vector<std::string>& argv, const std::string& stdinData,
                          std::chrono::seconds) override {
    Call c{argv, stdinData, ""};
    // upload-file = "<ruta>"
    const auto p = stdinData.find("upload-file = \"");
    if (p != std::string::npos) {
      const auto s = p + 15;
      const std::string path = stdinData.substr(s, stdinData.find('"', s) - s);
      std::ifstream in(path);
      c.body.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    calls.push_back(c);
    SpawnResult r;
    r.exitCode = calls.size() <= exitCodes.size() ? exitCodes[calls.size() - 1] : 0;
    r.stderrText = stderrText;
    return r;
  }
};

SmtpSettings smtp() {
  SmtpSettings s;
  s.url = "smtps://smtp.example.org:465";
  s.user = "daq@example.org";
  s.password = kPassword;
  s.from = "DAQ hodoscopio <daq@example.org>";
  s.to = {"a@example.org", "b@example.org"};
  return s;
}

std::optional<Mapping> mapping() {
  auto m = Mapping::fromFile(DAQ_SOURCE_DIR "/config/mapping/hodo_mapping_2026_10.json");
  EXPECT_TRUE(m.ok()) << (m.ok() ? "" : m.error().message);
  return m.ok() ? std::optional<Mapping>(m.value()) : std::nullopt;
}

std::string freshOutbox() {
  static int n = 0;
  const std::string d = ::testing::TempDir() + "daq_outbox_" + std::to_string(::getpid()) + "_" +
                        ::testing::UnitTest::GetInstance()->current_test_info()->name() + "_" + std::to_string(n++);
  return d;
}

int filesIn(const std::string& dir) {
  int n = 0;
  if (DIR* d = ::opendir(dir.c_str())) {
    while (dirent* e = ::readdir(d)) n += std::string(e->d_name).size() > 5 && std::string(e->d_name).find(".json") != std::string::npos;
    ::closedir(d);
  }
  return n;
}

// Primera línea del cuerpo (tras los encabezados).
std::string firstBodyLine(const std::string& msg) {
  const auto p = msg.find("\r\n\r\n");
  if (p == std::string::npos) return "";
  return msg.substr(p + 4, msg.find("\r\n", p + 4) - p - 4);
}

NotifyEvent gaveUp(int ch, std::uint64_t episode) {
  NotifyEvent e;
  e.kind = "gave_up";
  e.ch = ch;
  e.episode = episode;
  e.k = 3;
  e.reason = "trip durante la estabilización del intento 3";
  e.tUtcNs = 1792281600LL * 1000000000LL;
  e.timeline = {{{"t_utc", "2026-10-18T00:00:00Z"}, {"event", "trip"}, {"k", 0}}};
  e.context = {{"VMon_V", 0.2}, {"IMon_A", -4.8e-8}, {"Status", 512}, {"status_flags", {"internal_trip"}},
               {"board_temperature_degC", 25.0}, {"geco_running", false}};
  return e;
}

struct NRig {
  ManualClock clock;
  StateHub hub;
  FakeSpawner spawner;
  NotifierConfig cfg;
  std::unique_ptr<Notifier> n;
  explicit NRig(nlohmann::json c = {{"enabled", true}}) {
    // Outbox propio de cada test (nunca el real de ~/.local/state).
    if (!c.contains("outbox_dir")) c["outbox_dir"] = freshOutbox();
    cfg = NotifierConfig::fromJson(c).value();
    n = std::make_unique<Notifier>(cfg, smtp(), mapping(), spawner, clock, ::testing::TempDir(), hub, "host-test");
  }
};

}  // namespace

TEST(Notifier, GaveUpMailHasNoSecretInArgvOrBodyAndCorrectHeaders) {
  NRig r;
  ASSERT_TRUE(r.n->submit(gaveUp(3, 1)));
  r.n->processDue();
  ASSERT_EQ(r.spawner.calls.size(), 1u);
  const auto& c = r.spawner.calls[0];
  EXPECT_EQ(c.argv, (std::vector<std::string>{"curl", "--config", "-"}));
  for (const auto& a : c.argv) EXPECT_EQ(a.find("S3cr3t"), std::string::npos);
  // La contraseña va solo por stdin, escapada para el formato de configuración de curl.
  EXPECT_NE(c.stdinData.find("user = \"daq@example.org:S3cr3t-P4ss\\\"w\\\\rd\""), std::string::npos) << c.stdinData;
  EXPECT_NE(c.stdinData.find("mail-rcpt = \"a@example.org\""), std::string::npos);
  EXPECT_NE(c.stdinData.find("mail-rcpt = \"b@example.org\""), std::string::npos);
  EXPECT_NE(c.stdinData.find("mail-from = \"daq@example.org\""), std::string::npos);
  EXPECT_NE(c.stdinData.find("max-time = 30"), std::string::npos);
  // Cuerpo: sin secretos, con encabezados y contenido pedidos.
  ASSERT_FALSE(c.body.empty());
  EXPECT_EQ(c.body.find("S3cr3t"), std::string::npos);
  EXPECT_NE(c.body.find("Subject: [DAQ hodoscopio] gave_up CH03 (grupo 4)\r\n"), std::string::npos) << c.body;
  EXPECT_NE(c.body.find("To: a@example.org, b@example.org\r\n"), std::string::npos);
  // Date: "Sun, 18 Oct 2026 00:00:00 +0000" (RFC 5322, UTC).
  const auto dpos = c.body.find("Date: ");
  ASSERT_NE(dpos, std::string::npos);
  const std::string date = c.body.substr(dpos + 6, c.body.find("\r\n", dpos) - dpos - 6);
  EXPECT_EQ(date, "Sun, 18 Oct 2026 00:00:00 +0000");
  // Message-ID: "<ns.16hex@host-test>".
  const auto mpos = c.body.find("Message-ID: <");
  ASSERT_NE(mpos, std::string::npos);
  const std::string mid = c.body.substr(mpos + 13, c.body.find(">", mpos) - mpos - 13);
  const auto dot = mid.find('.'), at = mid.find('@');
  ASSERT_NE(dot, std::string::npos);
  ASSERT_NE(at, std::string::npos);
  EXPECT_EQ(mid.substr(0, dot).find_first_not_of("0123456789"), std::string::npos);
  EXPECT_EQ(at - dot - 1, 16u);
  EXPECT_EQ(mid.substr(dot + 1, 16).find_first_not_of("0123456789abcdef"), std::string::npos);
  EXPECT_EQ(mid.substr(at + 1), "host-test");
  EXPECT_NE(c.body.find("hod2019: ch50*"), std::string::npos) << c.body;  // grupo 4 = ch50, ch61 (derivados)
  EXPECT_NE(c.body.find("ch61*"), std::string::npos);
  EXPECT_NE(c.body.find("no confirmada por medición"), std::string::npos);
  EXPECT_NE(c.body.find("Temperatura de la placa A7030DP (no del detector): 25.0 degC"), std::string::npos);
  EXPECT_NE(c.body.find("GECO2020 en ejecución: no"), std::string::npos);
  EXPECT_NE(c.body.find("Host: host-test"), std::string::npos);
  // El archivo temporal del cuerpo se borró.
  const auto p = c.stdinData.find("upload-file = \"");
  const std::string path = c.stdinData.substr(p + 15, c.stdinData.find('"', p + 15) - p - 15);
  std::ifstream gone(path);
  EXPECT_FALSE(gone.good());
  EXPECT_EQ(r.n->status()["sent"], 1);
}

TEST(Notifier, OneMailPerChannelAndEpisode) {
  NRig r;
  EXPECT_TRUE(r.n->submit(gaveUp(3, 1)));
  EXPECT_FALSE(r.n->submit(gaveUp(3, 1)));  // mismo episodio
  EXPECT_TRUE(r.n->submit(gaveUp(3, 2)));   // episodio nuevo
  EXPECT_TRUE(r.n->submit(gaveUp(4, 1)));   // otro canal
  auto trip = gaveUp(3, 1);
  trip.kind = "trip";
  EXPECT_FALSE(r.n->submit(trip));  // notify_each_trip = false por defecto
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 3u);
}

TEST(Notifier, MaxMessagesCapsAcceptedMailsButNotTests) {
  NRig r({{"enabled", true}, {"notify_each_trip", true}, {"max_messages", 2}});
  auto trip = gaveUp(0, 1);
  trip.kind = "trip";
  EXPECT_TRUE(r.n->submit(trip));
  trip.ch = 1;
  EXPECT_TRUE(r.n->submit(trip));
  trip.ch = 2;
  EXPECT_FALSE(r.n->submit(trip));          // tercer correo: tope
  EXPECT_FALSE(r.n->submit(gaveUp(3, 1)));  // cualquier tipo
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 2u);
  EXPECT_EQ(r.n->status()["accepted"], 2);
  EXPECT_FALSE(NotifierConfig::fromJson({{"max_messages", -1}}).ok());
}

TEST(Notifier, DaqrunMailsAreSignedByDaqrunAndPointToItsDashboard) {
  NRig r({{"enabled", true}, {"page_ssh_host", "100.114.138.29"}});
  r.cfg.senderName = "daqrun";
  r.cfg.pagePort = 8080;
  r.cfg.pagePath = "/dash/";
  r.n = std::make_unique<Notifier>(r.cfg, smtp(), mapping(), r.spawner, r.clock, ::testing::TempDir(), r.hub, "host-test");
  auto e = gaveUp(-1, 1);
  e.kind = "no_triggers";
  e.ch = -1;
  ASSERT_TRUE(r.n->submit(e));
  r.n->processDue();
  ASSERT_EQ(r.spawner.calls.size(), 1u);
  const std::string& b = r.spawner.calls[0].body;
  EXPECT_NE(b.find("ssh -N -L 8080:127.0.0.1:8080"), std::string::npos);
  EXPECT_NE(b.find("http://localhost:8080/dash/"), std::string::npos);
  EXPECT_NE(b.find("daqrun (DAQ v2, hodoscopio)"), std::string::npos);
  EXPECT_EQ(b.find("hvmon (DAQ v2"), std::string::npos);
  EXPECT_EQ(b.find("8081"), std::string::npos);
}

TEST(Notifier, MaxMessagesPerDayResetsAtUtcMidnight) {
  NRig r({{"enabled", true}, {"max_messages_per_day", 2}});
  EXPECT_TRUE(r.n->submit(gaveUp(0, 1)));
  EXPECT_TRUE(r.n->submit(gaveUp(1, 1)));
  EXPECT_FALSE(r.n->submit(gaveUp(2, 1)));  // tercero del día
  r.clock.advance(std::chrono::hours(24));
  EXPECT_TRUE(r.n->submit(gaveUp(3, 1)));  // día nuevo
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 3u);
}

TEST(Notifier, OnlyKindsFiltersByKind) {
  NRig r({{"enabled", true}, {"notify_each_trip", true}, {"only_kinds", {"light_leak", "no_triggers"}}});
  auto trip = gaveUp(0, 1);
  trip.kind = "trip";
  EXPECT_FALSE(r.n->submit(trip));
  EXPECT_FALSE(r.n->submit(gaveUp(1, 1)));
  auto leak = gaveUp(-1, 1);
  leak.kind = "light_leak";
  leak.ch = -1;
  leak.context = {{"unit", "face:hod2019 Trasera"}, {"unit_name", "hod2019 Trasera"}};
  EXPECT_TRUE(r.n->submit(leak));
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 1u);
}

TEST(Notifier, ProblemMailsAtMostOnePer30Min) {
  NRig r;
  WriterHealth w;
  w.errors = 1;
  w.lastError = "No space left on device";
  r.hub.writer.publish(w);
  r.n->checkProblems();
  w.errors = 2;
  r.hub.writer.publish(w);
  r.clock.advance(10min);
  r.n->checkProblems();
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 1u);
  w.errors = 3;
  r.hub.writer.publish(w);
  r.clock.advance(25min);
  r.n->checkProblems();
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 2u);
}

TEST(Notifier, ConnectionLostAfterFiveMinutes) {
  NRig r;
  PollerState p;
  p.connected = false;
  p.everConnected = true;
  p.lastGoodReadUtcNs = r.clock.now_utc_ns();
  r.hub.poller.publish(p);
  r.clock.advance(4min);
  r.n->checkProblems();
  r.n->processDue();
  EXPECT_EQ(r.spawner.calls.size(), 0u);
  r.clock.advance(2min);
  r.n->checkProblems();
  r.n->checkProblems();  // una sola alerta por corte
  r.n->processDue();
  ASSERT_EQ(r.spawner.calls.size(), 1u);
  EXPECT_NE(r.spawner.calls[0].body.find("CONEXIÓN PERDIDA"), std::string::npos);
}

TEST(Notifier, OutageBackoffOutboxAndDeliveryDelay) {
  // La red cae 3 h: reintentos a 1, 5, 15, 60 min y luego cada hora; el correo sale al volver la red,
  // con la hora del evento y el retraso; el archivo del outbox desaparece al entregarse.
  NRig r;
  const std::string box = r.cfg.outboxDir;
  r.spawner.exitCodes = std::vector<int>(6, 6);  // curl: (6) no se pudo resolver el host
  r.spawner.stderrText = "curl: (6) Could not resolve host; pass " + kPassword;
  const std::string log = ::testing::TempDir() + "daq_notifier_outage.log";
  std::remove(log.c_str());
  Logger::instance().setStderr(false);
  ASSERT_TRUE(Logger::instance().openFile(log).ok());
  auto e = gaveUp(5, 1);
  e.tUtcNs = r.clock.now_utc_ns();
  ASSERT_TRUE(r.n->submit(std::move(e)));
  EXPECT_EQ(filesIn(box), 1);
  // Contenido del outbox: modo 600 y sin secretos.
  DIR* d = ::opendir(box.c_str());
  std::string f;
  while (dirent* x = ::readdir(d))
    if (std::string(x->d_name).find(".json") != std::string::npos) f = box + "/" + x->d_name;
  ::closedir(d);
  struct stat st {};
  ASSERT_EQ(::stat(f.c_str(), &st), 0);
  EXPECT_EQ(st.st_mode & 0777, 0600u);
  {
    std::ifstream in(f);
    const std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_EQ(all.find("S3cr3t"), std::string::npos);
    EXPECT_EQ(all.find("smtp"), std::string::npos);
  }
  // Minuto de cada intento, avanzando de a 1 min durante 4 h.
  std::vector<int> callMinutes;
  r.n->processDue();
  if (r.spawner.calls.size() == 1) callMinutes.push_back(0);
  for (int m = 1; m <= 240; ++m) {
    r.clock.advance(1min);
    const auto before = r.spawner.calls.size();
    r.n->processDue();
    if (r.spawner.calls.size() > before) callMinutes.push_back(m);
  }
  EXPECT_EQ(callMinutes, (std::vector<int>{0, 1, 6, 21, 81, 141, 201}));
  // Séptimo intento (t = 0+1+5+15+60+60+60 = 201 min) entregado: curl devuelve 0.
  ASSERT_EQ(r.spawner.calls.size(), 7u);
  EXPECT_EQ(r.n->pending(), 0u);
  EXPECT_EQ(filesIn(box), 0);
  const std::string& body = r.spawner.calls.back().body;
  EXPECT_NE(body.find("Hora del evento: 2026-10-18T00:00:00.000Z"), std::string::npos) << body;
  EXPECT_NE(body.find("Retraso de entrega: 12060 s (201.0 min)"), std::string::npos) << body;
  EXPECT_EQ(r.n->status()["sent"], 1);
  Logger::instance().closeFile();
  Logger::instance().setStderr(true);
  std::ifstream in(log);
  const std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(all.find("reintento en 60 min"), std::string::npos);
  EXPECT_EQ(all.find("S3cr3t"), std::string::npos);
  std::remove(log.c_str());
}

TEST(Notifier, OutboxSurvivesRestartAndAbandonsAfter48h) {
  const std::string box = freshOutbox();
  {
    NRig a({{"enabled", true}, {"outbox_dir", box}});
    a.spawner.exitCodes.assign(1, 7);
    ASSERT_TRUE(a.n->submit(gaveUp(3, 1)));
    a.n->processDue();  // falla
    EXPECT_EQ(filesIn(box), 1);
  }  // "hvmon" se reinicia: el correo pendiente sigue en el outbox
  NRig b({{"enabled", true}, {"outbox_dir", box}});
  EXPECT_EQ(b.n->restoredFromOutbox(), 1u);
  EXPECT_EQ(b.n->pending(), 1u);
  b.clock.advance(2min);
  b.n->processDue();
  ASSERT_EQ(b.spawner.calls.size(), 1u);
  EXPECT_NE(b.spawner.calls[0].body.find("gave_up CH03"), std::string::npos);
  EXPECT_EQ(filesIn(box), 0);
  // Sin red nunca: a las 48 h del evento se abandona y se aparta en failed/.
  const std::string box2 = freshOutbox();
  NRig c({{"enabled", true}, {"outbox_dir", box2}});
  c.spawner.exitCodes = std::vector<int>(100, 6);
  auto e = gaveUp(4, 2);
  e.tUtcNs = c.clock.now_utc_ns();
  ASSERT_TRUE(c.n->submit(std::move(e)));
  for (int h = 0; h < 50 * 4; ++h) {
    c.clock.advance(15min);
    c.n->processDue();
  }
  EXPECT_EQ(c.n->pending(), 0u);
  EXPECT_EQ(filesIn(box2), 0);
  EXPECT_EQ(filesIn(box2 + "/failed"), 1);
  EXPECT_EQ(c.n->status()["abandoned"], 1);
  // 1, 5, 15, 60 y luego cada hora hasta 48 h: 4 + ~47 intentos más el primero.
  EXPECT_GE(c.spawner.calls.size(), 48u);
  EXPECT_LE(c.spawner.calls.size(), 53u);
}

TEST(Notifier, TestEventIsMarkedAndRepliesWithExitCode) {
  NRig r({{"enabled", false}});  // notify-test funciona aunque el notificador esté deshabilitado
  auto e = gaveUp(11, 0);
  e.test = true;
  e.reply = std::make_shared<std::promise<nlohmann::json>>();
  auto fut = e.reply->get_future();
  ASSERT_TRUE(r.n->submit(std::move(e)));
  r.n->processDue();
  ASSERT_EQ(fut.wait_for(0s), std::future_status::ready);
  const auto j = fut.get();
  EXPECT_TRUE(j["ok"].get<bool>());
  EXPECT_EQ(j["curl_exit"], 0);
  ASSERT_EQ(r.spawner.calls.size(), 1u);
  EXPECT_NE(r.spawner.calls[0].body.find("Subject: [PRUEBA] [DAQ hodoscopio] gave_up CH11\r\n"), std::string::npos)
      << r.spawner.calls[0].body;
  EXPECT_EQ(firstBodyLine(r.spawner.calls[0].body).rfind("*** PRUEBA ***", 0), 0u);
  EXPECT_NE(r.spawner.calls[0].body.find("canal libre"), std::string::npos);  // CH11 no tiene MPPC
  // Un evento normal con el notificador deshabilitado se descarta.
  EXPECT_FALSE(r.n->submit(gaveUp(3, 9)));
}

TEST(Notifier, SmtpSettingsFromSecretsSplitsRecipients) {
  const std::string path = ::testing::TempDir() + "daq_test_smtp_secrets.env";
  {
    std::ofstream o(path);
    o << "DAQ_SMTP_URL=smtps://h:465\nDAQ_SMTP_USER=u\nDAQ_SMTP_PASSWORD=p\nDAQ_NOTIFY_FROM=f@x\n"
         "DAQ_NOTIFY_TO=a@x, b@x ,c@x\n";
  }
  ::chmod(path.c_str(), 0600);
  for (const char* v : {"DAQ_SMTP_URL", "DAQ_SMTP_USER", "DAQ_SMTP_PASSWORD", "DAQ_NOTIFY_FROM", "DAQ_NOTIFY_TO"}) ::unsetenv(v);
  Secrets s(path);
  auto r = SmtpSettings::fromSecrets(s);
  ASSERT_TRUE(r.ok()) << r.error().message;
  EXPECT_EQ(r.value().to, (std::vector<std::string>{"a@x", "b@x", "c@x"}));
  std::remove(path.c_str());
}

TEST(Notifier, SimulatedSourcesAreMarkedInSubjectAndFirstLine) {
  // Real: el asunto empieza con [DAQ hodoscopio] y no hay marca de simulación.
  {
    NRig r;
    ASSERT_TRUE(r.n->submit(gaveUp(3, 1)));
    r.n->processDue();
    ASSERT_EQ(r.spawner.calls.size(), 1u);
    EXPECT_NE(r.spawner.calls[0].body.find("Subject: [DAQ hodoscopio] gave_up CH03"), std::string::npos);
    EXPECT_EQ(r.spawner.calls[0].body.find("SIMULADO"), std::string::npos);
  }
  // FakeHVBackend: [SIMULADO] al inicio del asunto y primera línea del cuerpo.
  {
    NRig r;
    r.n->setSimulatedBackend(true);
    ASSERT_TRUE(r.n->submit(gaveUp(8, 1)));
    r.n->processDue();
    ASSERT_EQ(r.spawner.calls.size(), 1u);
    const auto& b = r.spawner.calls[0].body;
    EXPECT_NE(b.find("Subject: [SIMULADO] [DAQ hodoscopio] gave_up CH08\r\n"), std::string::npos) << b;
    EXPECT_EQ(firstBodyLine(b),
              "*** SIMULADO *** Evento simulado (FakeHVBackend): ningún canal real estuvo involucrado.");
  }
  // Canal en test_mode con la fuente real: también [SIMULADO].
  {
    NRig r;
    auto e = gaveUp(9, 1);
    e.simulated = true;
    ASSERT_TRUE(r.n->submit(std::move(e)));
    r.n->processDue();
    ASSERT_EQ(r.spawner.calls.size(), 1u);
    const auto& b = r.spawner.calls[0].body;
    EXPECT_NE(b.find("Subject: [SIMULADO] [DAQ hodoscopio] gave_up CH09\r\n"), std::string::npos) << b;
    EXPECT_EQ(firstBodyLine(b), "*** SIMULADO *** Evento simulado (test_mode del canal CH09): ningún canal real estuvo "
                                "involucrado.");
  }
  // notify-test con fuente simulada: [PRUEBA] manda en el asunto; el cuerpo lo dice igual.
  {
    NRig r({{"enabled", false}});
    r.n->setSimulatedBackend(true);
    auto e = gaveUp(11, 0);
    e.test = true;
    ASSERT_TRUE(r.n->submit(std::move(e)));
    r.n->processDue();
    ASSERT_EQ(r.spawner.calls.size(), 1u);
    const auto& b = r.spawner.calls[0].body;
    EXPECT_NE(b.find("Subject: [PRUEBA] [DAQ hodoscopio] gave_up CH11\r\n"), std::string::npos);
    EXPECT_EQ(firstBodyLine(b).rfind("*** PRUEBA ***", 0), 0u);
    EXPECT_NE(b.find("Fuente HV simulada (FakeHVBackend)."), std::string::npos);
  }
}

TEST(Notifier, BodyUnitsLocalTimeWaitsAndPageLine) {
  ASSERT_TRUE(setProcessTimeZone("America/Santiago").ok());
  EXPECT_FALSE(setProcessTimeZone("Marte/Olympus").ok());
  NRig r({{"enabled", true}, {"page_ssh_host", "tailnet-host"}});
  auto e = gaveUp(3, 1);
  const UtcNs t0 = 1792281600LL * 1000000000LL;  // 2026-10-18T00:00:00Z = 21:00 -03 del 17
  e.timeline = nlohmann::json::array(
      {{{"t_utc_ns", t0}, {"event", "trip"}, {"k", 0}},
       {{"t_utc_ns", t0 + 300LL * 1000000000LL}, {"event", "retry"}, {"k", 1}, {"wait_min", 5.0}},
       {{"t_utc_ns", t0 + 302LL * 1000000000LL + 5000000LL}, {"event", "ramp_done"}, {"k", 1}}});
  e.context["IMon_A"] = -1.22e-7;
  e.context["VMon_V"] = 65.7468;
  ASSERT_TRUE(r.n->submit(std::move(e)));
  r.n->processDue();
  ASSERT_EQ(r.spawner.calls.size(), 1u);
  const auto& b = r.spawner.calls[0].body;
  EXPECT_NE(b.find("IMon = -0.122 µA\r\n"), std::string::npos) << b;
  EXPECT_NE(b.find("VMon = 65.75 V\r\n"), std::string::npos);
  EXPECT_NE(b.find("Hora del evento: 2026-10-18T00:00:00.000Z (2026-10-17 21:00:00 -03)"), std::string::npos);
  EXPECT_NE(b.find("2026-10-18T00:05:00.000Z (2026-10-17 21:05:00 -03)  retry (k = 1), espera aplicada 5 min"),
            std::string::npos);
  EXPECT_NE(b.find("2026-10-18T00:05:02.005Z (2026-10-17 21:05:02 -03)  ramp_done (k = 1)\r\n"), std::string::npos);
  EXPECT_NE(b.find("ssh -N -L 8081:127.0.0.1:8081 xplus@tailnet-host y abrir http://localhost:8081"), std::string::npos);
  // Sin host configurado la línea lo dice; el puerto viene de http.port de hvmon.
  NotifierConfig c2 = NotifierConfig::fromJson({{"enabled", true}, {"outbox_dir", freshOutbox()}}).value();
  c2.pagePort = 8091;
  c2.pageSshHost = "otro";
  FakeSpawner sp;
  Notifier n2(c2, smtp(), mapping(), sp, r.clock, ::testing::TempDir(), r.hub, "h");
  std::string subj;
  EXPECT_NE(n2.composeMessage(gaveUp(3, 1), &subj).find("ssh -N -L 8091:127.0.0.1:8091 xplus@otro y abrir http://localhost:8091"),
            std::string::npos);
  NotifierConfig c3 = NotifierConfig::fromJson({{"enabled", true}, {"outbox_dir", freshOutbox()}}).value();
  Notifier n3(c3, smtp(), mapping(), sp, r.clock, ::testing::TempDir(), r.hub, "h");
  EXPECT_NE(n3.composeMessage(gaveUp(3, 1)).find("notifier.page_ssh_host sin configurar"), std::string::npos);
}
