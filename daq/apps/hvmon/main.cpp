// hvmon: daemon de control lento de HV (fase 1). Único dueño de la sesión de la Wrapper.
//   hvmon [--config <archivo.json>] [--fake <guion.json|default>] [--print-config]
// Hilos: sondeo (backend HV), escritor (archivos diarios), comandos (socket Unix), temperatura
// (si hay backend) y principal (estado en /dev/shm, THttpServer, detección de GECO). Todo objeto
// ROOT pertenece al hilo principal. SIGINT/SIGTERM: parada ordenada (sondeo, comandos,
// temperatura, escritor con fsync, principal), estado final marcado "final": true.
// Secretos: DAQ_HV_HOST/DAQ_HV_USER/DAQ_HV_PASSWORD del entorno o de secrets.env; nunca en argv,
// logs, estado ni HTTP.
#include <signal.h>
#include <sys/resource.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>

#include <ROOT/RConfig.hxx>
#include <TH1.h>
#include <TROOT.h>

#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/core/instance_lock.h"
#include "daq/core/jthread.h"
#include "daq/core/logger.h"
#include "daq/core/secrets.h"
#include "daq/core/time.h"
#include "daq/core/unix_socket.h"
#include "daq/hv/caen_hv_backend.h"
#include "daq/hv/credentials.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hvmon/commands.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/poller.h"
#include "daq/hvmon/state.h"
#include "daq/hvmon/system_probe.h"
#include "daq/hvmon/writer.h"
#include "daq/hvmon/compressor.h"
#include "daq/hvmon/mapping.h"
#include "daq/hvmon/notifier.h"
#include "daq/core/spawn.h"
#include "daq/temp/temp_backend.h"
#include "http_page.h"

using namespace daq;
using namespace daq::hvmon;

namespace {

std::atomic<bool> gStop{false};
static_assert(std::atomic<bool>::is_always_lock_free, "el flag de señal debe ser lock-free");

extern "C" void onSignal(int) { gStop.store(true); }

void installSignals() {
  struct sigaction sa {};
  sa.sa_handler = onSignal;
  sigemptyset(&sa.sa_mask);
  sigaction(SIGINT, &sa, nullptr);
  sigaction(SIGTERM, &sa, nullptr);
  signal(SIGPIPE, SIG_IGN);
}

void disableCoreDumps() {
  rlimit rl{0, 0};
  setrlimit(RLIMIT_CORE, &rl);  // un volcado de memoria podría contener la contraseña
}

const char* kUsage =
    "Uso: hvmon [--config archivo.json] [--fake guion.json|default] [--print-config]\n"
    "Sin --config usa ~/.config/daqv2/hvmon.json si existe, si no los valores por defecto.\n";

}  // namespace

int main(int argc, char** argv) {
  disableCoreDumps();
  Logger::setThreadName("main");
  std::string configPath, fake;
  bool printConfig = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--config" && i + 1 < argc) configPath = argv[++i];
    else if (a == "--fake" && i + 1 < argc) fake = argv[++i];
    else if (a == "--print-config") printConfig = true;
    else if (a == "--help" || a == "-h") {
      std::fputs(kUsage, stdout);
      return 0;
    } else {
      std::fprintf(stderr, "hvmon: argumento desconocido: %s\n%s", a.c_str(), kUsage);
      return 2;
    }
  }
  if (configPath.empty()) {
    const std::string def = hv::expandHome("~/.config/daqv2/hvmon.json");
    if (PosixFileSystem().exists(def)) configPath = def;
  }
  const bool fakeBackend = !fake.empty();
  auto cfgR = configPath.empty() ? HvmonConfig::fromJson(nlohmann::json::object(), fakeBackend)
                                 : HvmonConfig::fromFile(configPath, fakeBackend);
  if (!cfgR) {
    std::fprintf(stderr, "hvmon: %s\n", cfgR.error().message.c_str());
    return 2;
  }
  const HvmonConfig cfg = cfgR.value();
  if (printConfig) {
    std::printf("%s\n", redactSecrets(cfg.toJson()).dump(2).c_str());
    return 0;
  }
  if (!cfg.logFile.empty()) {
    if (auto st = Logger::instance().openFile(cfg.logFile); !st) std::fprintf(stderr, "hvmon: %s\n", st.error().message.c_str());
  }
  // Hora local de los correos: TZ del proceso, fijada antes de crear hilos.
  if (auto st = setProcessTimeZone(cfg.notifier.localTimezone); !st) {
    std::fprintf(stderr, "hvmon: notifier.local_timezone: %s\n", st.error().message.c_str());
    return 2;
  }
  logInfo("hvmon arranca (pid ", ::getpid(), "), configuración ", configPath.empty() ? "por defecto" : configPath, ": ",
          redactSecrets(cfg.toJson()).dump());

  for (const auto& d : {cfg.runtimeDir, cfg.stateDir}) {
    if (auto st = ensurePrivateDirectory(d, 0700); !st) {
      logError(st.error().message);
      return 1;
    }
  }
  auto lock = InstanceLock::acquire(cfg.lockPath());
  if (!lock) {
    std::fprintf(stderr, "hvmon: ya hay un hvmon (u otro dueño de la sesión HV) corriendo: %s\n",
                 lock.error().message.c_str());
    logError("arranque rechazado: ", lock.error().message);
    return 3;
  }
  PosixFileSystem fs;
  if (auto st = fs.makeDirectories(cfg.slowcontrolDir, 0755); !st) {
    logError(st.error().message);
    return 1;
  }

  SystemClock clock;
  std::unique_ptr<hv::HVBackend> backend;
  if (!fake.empty()) {
    auto f = hv::FakeHVBackend::fromFile(fake, clock);
    if (!f) {
      logError(f.error().message);
      return 1;
    }
    backend = std::move(f).value();
    logWarn("usando FakeHVBackend (", fake, "): no hay conexión con la fuente real");
  } else {
    Secrets secrets;
    auto login = hv::resolveHvLogin(secrets, cfg.host, cfg.legacyCredentials);
    if (!login) {
      logError(login.error().message);
      return 1;
    }
    logInfo("credenciales HV: fuente ", login.value().source, ", usuario ", login.value().credentials.user, ", host ",
            login.value().host);
    backend = std::make_unique<hv::CaenHVBackend>(
        hv::ConnectionConfig{cfg.system, login.value().host, login.value().credentials}, &clock);
  }

  BoundedQueue<WriterItem> writerQ(cfg.writerQueueCapacity);
  BoundedQueue<PollRequest> pollQ(64);
  StateHub hub;
  hub.startedUtcNs = clock.now_utc_ns();
  SlowControlWriter writer(cfg.slowcontrolDir, cfg.fsyncPeriodS, fs, clock, hub);
  Poller poller(cfg, *backend, clock, writerQ, hub, pollQ);
  poller.enableRetryPersistence(fs);
  if (!cfg.retry.enabledChannels.empty())
    logWarn("auto-retry HABILITADO en ", cfg.retry.enabledChannels.size(), " canales (allow_commands = true)",
            cfg.retry.testModeEnabled ? "; test_mode activo" : "");
  CommandHandler handler(cfg, hub, pollQ, writerQ, clock);

  // Fase 2: compresión de días cerrados.
  BoundedQueue<std::string> compressQ(1000);
  Compressor compressor(hub);
  if (cfg.compressClosedDays) {
    writer.setCompressionQueue(&compressQ);
    for (const auto& f : Compressor::closedDailyFiles(cfg.slowcontrolDir, SlowControlWriter::dateKey(clock.now_utc_ns())))
      compressQ.try_push(f, DropPolicy::DropNewest);
  }

  // Fase 2: notificador por correo (secretos DAQ_SMTP_*, DAQ_NOTIFY_* de secrets.env o del entorno).
  BoundedQueue<NotifyEvent> notifyQ(256);
  PosixSpawner spawner;
  std::unique_ptr<Notifier> notifier;
  {
    Secrets smtpSecrets;
    auto smtp = SmtpSettings::fromSecrets(smtpSecrets);
    if (!smtp) {
      if (cfg.notifier.enabled) {
        logError("notificador habilitado pero sin configuración SMTP: ", smtp.error().message);
        return 1;
      }
      logWarn("notificador sin SMTP (", smtp.error().message, "): notify-test no disponible");
    } else {
      std::optional<Mapping> mapping;
      auto m = Mapping::fromFile(cfg.notifier.mappingFile);
      if (m) mapping = m.value();
      else logWarn("notificador sin mapeo: ", m.error().message);
      notifier = std::make_unique<Notifier>(cfg.notifier, smtp.value(), mapping, spawner, clock, cfg.runtimeDir, hub);
      notifier->setSimulatedBackend(fakeBackend);
      poller.setNotifyQueue(&notifyQ);
      handler.setNotifyQueue(&notifyQ);
      logInfo("notificador: ", cfg.notifier.enabled ? "habilitado" : "deshabilitado (solo notify-test)", ", ",
              smtp.value().to.size(), " destinatarios");
    }
  }
  auto listener = UnixListener::listen(cfg.socketPath());
  if (!listener) {
    logError(listener.error().message);
    return 1;
  }
  CommandServer server(std::move(listener).value(), handler);

  auto tempR = temp::makeTempBackend(cfg.tempBackend, clock, cfg.tempConfig);
  if (!tempR) {
    logError(tempR.error().message);
    return 1;
  }
  std::unique_ptr<temp::TempBackend> tempBackend = std::move(tempR).value();

  ROOT::EnableThreadSafety();
  TH1::AddDirectory(kFALSE);
  std::unique_ptr<HttpPage> page;
  if (cfg.httpEnabled) {
    auto bind = resolveHttpBind(cfg.httpBind);
    if (!bind) {
      logError(bind.error().message);
      return 1;
    }
    auto p = HttpPage::create(cfg, bind.value());
    if (!p) {
      logError(p.error().message);
      return 1;
    }
    page = std::move(p).value();
    logInfo("página web de solo lectura en ", page->url());
  }

  installSignals();
  JThread writerThread([&] {
    Logger::setThreadName("writer");
    writer.run(writerQ);
  });
  JThread pollThread([&](StopToken st) { poller.run(st); });
  JThread cmdThread([&](StopToken st) { server.run(st); });
  JThread notifyThread;
  if (notifier) notifyThread = JThread([&](StopToken st) { notifier->run(notifyQ, st); });
  JThread compressThread([&](StopToken st) { compressor.run(compressQ, st); });
  JThread tempThread;
  if (tempBackend) {
    tempThread = JThread([&](StopToken st) {
      Logger::setThreadName("temp");
      while (!st.stop_requested()) {
        auto r = tempBackend->read();
        if (r) {
          nlohmann::json sensors = nlohmann::json::array();
          for (const auto& x : r.value())
            sensors.push_back({{"id", x.sensorId}, {"T_degC", x.tDegC}, {"t_utc_ns", x.tUtcNs}});
          hub.temperature.publish({{"backend", tempBackend->describe()}, {"sensors", sensors}});
          if (!writerQ.try_push(TempItem{clock.now_utc_ns(), r.value()}, DropPolicy::DropNewest)) ++hub.writerDropped;
        } else {
          logWarn("temperatura: ", r.error().message);
          hub.temperature.publish({{"backend", tempBackend->describe()}, {"error", r.error().message}});
        }
        clock.sleep_for(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(cfg.tempPeriodS)), st);
      }
    });
  }

  std::uint64_t stateSeq = 0;
  SteadyTime nextState = clock.now_steady();
  SteadyTime nextGeco = clock.now_steady();
  std::string lastStateError;
  bool geco = false;
  while (!gStop.load()) {
    const SteadyTime now = clock.now_steady();
    if (now >= nextGeco) {
      const bool running = processRunning(cfg.gecoPattern);
      if (running != geco) logInfo("GECO2020 ", running ? "detectado en ejecución" : "ya no está en ejecución");
      geco = running;
      hub.gecoRunning.store(running);
      nextGeco = now + std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(cfg.gecoPeriodS));
    }
    if (now >= nextState) {
      const UtcNs utc = clock.now_utc_ns();
      const auto state = hub.compose(cfg, utc, ++stateSeq, false);
      if (auto st = fs.writeAtomic(cfg.statePath(), state.dump() + "\n", 0600); !st) {
        if (st.error().message != lastStateError) logError("estado: ", st.error().message);
        lastStateError = st.error().message;
      }
      if (page) page->update(hub.poller.read().first, state, utc);
      nextState += std::chrono::seconds(1);
      if (nextState <= now) nextState = now + std::chrono::seconds(1);
    }
    if (page) page->processRequests();
    clock.sleep_for(std::chrono::milliseconds(66));  // ~15 Hz
  }

  logInfo("señal recibida: parada ordenada");
  pollThread.request_stop();
  pollThread.join();
  cmdThread.request_stop();
  cmdThread.join();
  if (tempThread.joinable()) {
    tempThread.request_stop();
    tempThread.join();
  }
  writerQ.close();
  writerThread.join();
  compressQ.close();
  compressThread.request_stop();
  compressThread.join();
  if (notifyThread.joinable()) {
    notifyQ.close();
    notifyThread.request_stop();
    notifyThread.join();
  }
  const auto finalState = hub.compose(cfg, clock.now_utc_ns(), ++stateSeq, true);
  if (auto st = fs.writeAtomic(cfg.statePath(), finalState.dump() + "\n", 0600); !st)
    logError("estado final: ", st.error().message);
  page.reset();
  const auto w = hub.writer.read().first;
  logInfo("hvmon termina: fsync ", w.fsyncs, ", bytes ", w.bytes, ", errores de escritura ", w.errors);
  return 0;
}
