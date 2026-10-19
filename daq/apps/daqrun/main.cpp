// daqrun: un run de adquisición (fase 3: solo digitizador fake o replay).
//   daqrun --config <archivo.json> [--campaign NOMBRE] [--max-events N] [--max-duration S]
//          [--no-http] [--print-config]
// Hilos (diseño 2.1): lectura (único dueño del backend), escritor (único dueño del TFile), monitor
// (datos planos), cliente HV (lee el estado de hvmon) y principal (THttpServer, señales). Colas:
// escritura con push bloqueante (nunca se descarta: si se llena, la lectura espera y se cuenta),
// monitor con descarte del más antiguo, HV y resúmenes hacia el escritor.
// Apagado (diseño 3.2): SIGINT/SIGTERM solo marcan un flag; el principal detiene la lectura (que
// vacía la placa y cierra las colas), el escritor vacía la cola, espera flag_delay_s para fijar los
// bits HV, cierra el segmento (SHA256, a-w), escribe el índice completo y el manifiesto; luego se
// detiene el cliente HV y se libera el lock del run. Si el escritor no termina en
// writer_shutdown_timeout_s (120 s), sale con código 4 y el índice queda con complete = false.
#include <signal.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>

#include <TH1.h>
#include <TROOT.h>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/core/instance_lock.h"
#include "daq/core/jthread.h"
#include "daq/core/logger.h"
#include "daq/core/sha256.h"
#include "daq/daqrun/campaign.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/monitor.h"
#include "daq/dig/fake_digitizer.h"
#include "daq/dig/readout.h"
#include "daq/dig/replay_digitizer.h"
#include "daq/core/secrets.h"
#include "daq/core/spawn.h"
#include "daq/core/time.h"
#include "daq/hvmon/mapping.h"
#include "daq/hvmon/notifier.h"
#include "daq/hvmon/state.h"
#include "daq/hvmon/system_probe.h"
#include "daq/io/run_writer.h"
#include "daq_build_info.h"
#include "http_page.h"

using namespace daq;
using namespace daq::daqrun;

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

const char* kUsage =
    "Uso: daqrun --config archivo.json [--campaign NOMBRE] [--max-events N] [--max-duration S] [--no-http]\n"
    "            [--print-config]\n"
    "Fase 3: digitizer.backend = fake | replay (sin acceso al digitizador real).\n";

std::string hostName() {
  char h[256] = {0};
  return ::gethostname(h, sizeof h - 1) == 0 ? h : "unknown";
}

Result<nlohmann::json> loadScenario(const dig::DigitizerConfig& d) {
  if (!d.scenario.is_null()) return d.scenario;
  if (d.scenarioFile.empty()) return nlohmann::json::object();
  std::ifstream in(d.scenarioFile);
  if (!in) return Result<nlohmann::json>::failure(-1, "no se pudo abrir el escenario " + d.scenarioFile);
  try {
    nlohmann::json j;
    in >> j;
    return j;
  } catch (const std::exception& e) {
    return Result<nlohmann::json>::failure(-1, "escenario inválido " + d.scenarioFile + ": " + e.what());
  }
}

std::vector<dig::ChannelInfo> channelInfo(const hvmon::Mapping& m) {
  auto v = dig::defaultChannelInfo();
  for (const auto& e : m.channels) {
    if (e.digCh < 0 || e.digCh >= dig::kChannels) continue;
    auto& c = v[static_cast<std::size_t>(e.digCh)];
    c.station = e.station == "hod2019" ? 1 : 0;
    c.plane = e.plane == "Y" ? 'Y' : 'X';
    if (e.strip >= 0) c.strip = e.strip;
    c.hvGroup = e.hvGroup ? *e.hvGroup : 0;
    c.face = e.face;
  }
  return v;
}

}  // namespace

int main(int argc, char** argv) {
  Logger::setThreadName("main");
  std::string configPath, campaign;
  long long maxEvents = -1;
  double maxDuration = -1;
  bool printConfig = false, noHttp = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--config" && i + 1 < argc) configPath = argv[++i];
    else if (a == "--campaign" && i + 1 < argc) campaign = argv[++i];
    else if (a == "--max-events" && i + 1 < argc) maxEvents = std::atoll(argv[++i]);
    else if (a == "--max-duration" && i + 1 < argc) maxDuration = std::atof(argv[++i]);
    else if (a == "--print-config") printConfig = true;
    else if (a == "--no-http") noHttp = true;
    else if (a == "--help" || a == "-h") {
      std::fputs(kUsage, stdout);
      return 0;
    } else {
      std::fprintf(stderr, "daqrun: argumento desconocido: %s\n%s", a.c_str(), kUsage);
      return 2;
    }
  }
  if (configPath.empty()) {
    std::fputs(kUsage, stderr);
    return 2;
  }
  nlohmann::json cj;
  {
    std::ifstream in(configPath);
    if (!in) {
      std::fprintf(stderr, "daqrun: no se pudo abrir %s\n", configPath.c_str());
      return 2;
    }
    try {
      in >> cj;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "daqrun: JSON inválido en %s: %s\n", configPath.c_str(), e.what());
      return 2;
    }
  }
  if (!campaign.empty()) cj["campaign"]["name"] = campaign;
  if (maxEvents >= 0) cj["run"]["max_events"] = maxEvents;
  if (maxDuration >= 0) cj["run"]["max_duration_s"] = maxDuration;
  if (noHttp) cj["monitoring"]["http"]["enabled"] = false;
  auto cfgR = DaqrunConfig::fromJson(cj);
  if (!cfgR) {
    std::fprintf(stderr, "daqrun: %s\n", cfgR.error().message.c_str());
    return 2;
  }
  DaqrunConfig cfg = cfgR.value();
  auto scen = loadScenario(cfg.digitizer);
  if (!scen) {
    std::fprintf(stderr, "daqrun: %s\n", scen.error().message.c_str());
    return 2;
  }
  cfg.digitizer.scenario = scen.value();
  if (printConfig) {
    std::printf("%s\n", cfg.toJson().dump(2).c_str());
    return 0;
  }
  if (cfg.campaign.empty()) {
    std::fprintf(stderr, "daqrun: falta la campaña (campaign.name o --campaign)\n");
    return 2;
  }
  auto mapping = hvmon::Mapping::fromFile(cfg.mappingFile);
  if (!mapping) {
    std::fprintf(stderr, "daqrun: %s\n", mapping.error().message.c_str());
    return 2;
  }
  std::string mappingText;
  {
    std::ifstream in(cfg.mappingFile, std::ios::binary);
    mappingText.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }

  if (auto st = setProcessTimeZone(cfg.notifier.value("local_timezone", "America/Santiago")); !st)
    logWarn("zona horaria de los correos: ", st.error().message);
  SystemClock clock;
  const CampaignLocation loc{cfg.rawRoot, cfg.manifestsDir, cfg.campaign};
  auto res = RunReservation::reserve(loc, clock.now_utc_ns());
  if (!res) {
    std::fprintf(stderr, "daqrun: %s\n", res.error().message.c_str());
    return 3;
  }
  RunReservation run = std::move(res).value();
  logInfo("daqrun arranca (pid ", ::getpid(), "): campaña ", cfg.campaign, ", run ", run.run(), " (", run.baseName(),
          "), digitizador ", cfg.digitizer.backend, cfg.simulated() ? " [SIMULADO]" : "", ", configuración ", configPath);

  // Backend del digitizador.
  std::unique_ptr<dig::DigitizerBackend> backend;
  if (cfg.digitizer.backend == "fake") {
    auto sp = dig::SyntheticParams::fromJson(cfg.digitizer.scenario);
    if (!sp) {
      logError(sp.error().message);
      return 2;
    }
    backend = std::make_unique<dig::FakeDigitizer>(sp.value(), clock, channelInfo(mapping.value()));
  } else {
    auto rp = dig::ReplayParams::fromJson(cfg.digitizer.scenario);
    if (!rp) {
      logError(rp.error().message);
      return 2;
    }
    backend = std::make_unique<dig::ReplayDigitizer>(rp.value(), clock);
  }
  if (auto s = backend->configure(cfg.digitizer); !s) {
    logError("configuración del digitizador: ", s.error().message);
    return 1;
  }

  // Colas e hilos.
  BoundedQueue<dig::EventBlockPtr> writerQ(cfg.writerQueueCapacity);
  BoundedQueue<dig::EventBlockPtr> monitorQ(cfg.monitorQueueCapacity);
  BoundedQueue<HvSample> hvQ(cfg.hvQueueCapacity);
  BoundedQueue<MonitorSummary> sumQ(64);
  HvClient hvClient(cfg, clock, hvQ);
  HvSample hvStart = hvClient.pollOnce();
  if (!hvStart.alive) logWarn("hvmon no disponible al arrancar (", hvStart.error, "): bit 4 en todos los canales");

  // Configuración resuelta junto a los datos.
  const std::string configText = cfg.toJson().dump(2) + "\n";
  const std::string configSha = sha256Hex(configText);
  {
    PosixFileSystem fs;
    const std::string cfgFile = run.baseName() + ".config.json";
    if (auto s = fs.writeAtomic(loc.dir() + "/" + cfgFile, configText, 0444); !s) {
      logError(s.error().message);
      return 1;
    }
    appendSha256Sum(loc.dir(), cfgFile, configSha);
  }

  auto origin = backend->start();
  if (!origin) {
    logError("start del digitizador: ", origin.error().message);
    return 1;
  }
  ROOT::EnableThreadSafety();
  TH1::AddDirectory(kFALSE);
  io::RunStatic rs;
  rs.configJson = cfg.toJson();
  rs.configSha256 = configSha;
  rs.mappingText = mappingText;
  rs.digitizerInfo = backend->info();
  rs.origin = origin.value();
  rs.hvAtStart = hvStart.raw.is_null() ? nlohmann::json{{"error", hvStart.error}} : hvStart.raw;
  rs.gitCommit = DAQ_GIT_COMMIT;
  rs.gitDirty = DAQ_GIT_DIRTY;
  rs.host = hostName();
  rs.rootVersion = gROOT->GetVersion();
  io::WriterStats wstats;
  auto writer = std::make_unique<io::RunWriter>(cfg, run, mapping.value(), rs, clock, wstats);
  if (auto s = writer->open(); !s) {
    logError(s.error().message);
    return 1;
  }
  hvStart.raw = nlohmann::json();
  hvQ.push(hvStart);

  dig::ReadoutStats rstats;
  dig::ReadoutOptions ro;
  ro.pollSleepS = cfg.digitizer.pollSleepS;
  ro.maxEventsBlt = static_cast<std::size_t>(cfg.digitizer.maxEventsBlt);
  ro.maxEvents = cfg.maxEvents;
  ro.maxDurationS = cfg.maxDurationS;
  ro.softwareTrigger = cfg.digitizer.softwareTrigger;
  ro.swTriggerRateHz = cfg.digitizer.swTriggerRateHz;
  dig::ReadoutLoop readout(*backend, clock, writerQ, monitorQ, ro, rstats);
  auto queueLevels = [&] {
    return nlohmann::json{{"writer", writerQ.size()}, {"monitor", monitorQ.size()}, {"hv", hvQ.size()}};
  };
  Monitor monitor(cfg, rs.origin.utcBeforeStart, queueLevels);

  std::unique_ptr<DaqrunPage> page;
  if (cfg.httpEnabled) {
    auto bind = hvmon::resolveHttpBind(cfg.httpBind);
    if (!bind) {
      logError(bind.error().message);
      return 1;
    }
    auto p = DaqrunPage::create(cfg, bind.value());
    if (!p) {
      logError(p.error().message);
      return 1;
    }
    page = std::move(p).value();
    logInfo("página de solo lectura en ", page->url());
  }

  // Aviso por correo si hvmon deja de responder durante el run (notificador de hvmon, outbox propio).
  BoundedQueue<hvmon::NotifyEvent> notifyQ(64);
  hvmon::StateHub notifyHub;  // el notificador de daqrun no mira escritor ni conexión de hvmon
  PosixSpawner spawner;
  std::unique_ptr<hvmon::Notifier> notifier;
  nlohmann::json ncfgJson = cfg.notifier;
  if (!ncfgJson.contains("outbox_dir")) ncfgJson["outbox_dir"] = hvmon::defaultOutboxDir() + "-daqrun";
  ncfgJson.erase("hvmon_stale_alert_min");
  auto ncfg = hvmon::NotifierConfig::fromJson(ncfgJson);
  if (!ncfg) {
    logError("notifier: ", ncfg.error().message);
    return 2;
  }
  if (ncfg.value().enabled) {
    Secrets secrets;
    auto smtp = hvmon::SmtpSettings::fromSecrets(secrets);
    if (!smtp) {
      logWarn("notificador de daqrun habilitado pero sin SMTP (", smtp.error().message, "): sin avisos por correo");
    } else {
      const std::string tmp = std::string(std::getenv("XDG_RUNTIME_DIR") ? std::getenv("XDG_RUNTIME_DIR") : "/tmp") + "/daqv2-daqrun";
      ensurePrivateDirectory(tmp, 0700);
      notifier = std::make_unique<hvmon::Notifier>(ncfg.value(), smtp.value(), mapping.value(), spawner, clock, tmp, notifyHub);
      notifier->setSimulatedBackend(cfg.simulated());
      logInfo("notificador de daqrun: aviso si hvmon no responde más de ", cfg.hvmonStaleAlertMin, " min");
    }
  }
  HvmonWatchdog watchdog(cfg.hvmonStaleAlertMin * 60.0,
                         "daqrun campaña " + cfg.campaign + " run " + std::to_string(run.run()) + " (" + run.baseName() + ")");
  if (notifier)
    hvClient.setOnSample([&](const HvSample& s) {
      if (auto e = watchdog.feed(s)) {
        logWarn("hvmon no responde desde hace más de ", cfg.hvmonStaleAlertMin, " min: aviso por correo");
        notifyQ.try_push(std::move(*e), DropPolicy::DropNewest);
      }
    });

  installSignals();
  std::atomic<bool> writerDone{false}, writerOk{true};
  JThread notifyThread;
  if (notifier) notifyThread = JThread([&](StopToken st) { notifier->run(notifyQ, st); });
  JThread hvThread([&](StopToken st) { hvClient.run(st); });
  JThread monitorThread([&] {
    Logger::setThreadName("monitor");
    SteadyTime nextView = clock.now_steady();
    for (;;) {
      auto b = monitorQ.pop_until(std::chrono::steady_clock::now() + std::chrono::seconds(1));
      if (b) monitor.process(**b);
      MonitorSummary s;
      if (monitor.takeSummary(clock.now_utc_ns(), s, rstats.monitorDropped.load())) sumQ.try_push(std::move(s), DropPolicy::DropOldest);
      if (clock.now_steady() >= nextView) {
        monitor.publishView();
        nextView = clock.now_steady() + std::chrono::milliseconds(500);
      }
      if (!b && monitorQ.closed() && monitorQ.size() == 0) break;
    }
    monitor.publishView();
  });
  JThread writerThread([&] {
    Logger::setThreadName("writer");
    try {
      auto drainSide = [&] {
        while (auto h = hvQ.try_pop()) writer->addHv(*h);
        while (auto m = sumQ.try_pop()) writer->addSummary(*m);
      };
      for (;;) {
        drainSide();
        auto b = writerQ.pop_until(std::chrono::steady_clock::now() + std::chrono::milliseconds(100));
        if (b) {
          writer->addBlock(*b);
          for (int k = 0; k < 64; ++k) {
            auto more = writerQ.try_pop();
            if (!more) break;
            writer->addBlock(*more);
          }
        } else if (writerQ.closed() && writerQ.size() == 0) {
          break;
        }
        if (auto s = writer->tick(); !s) {
          logError("escritor: ", s.error().message);
          writerOk = false;
          break;
        }
      }
      // La lectura cerró la cola después de su último uso del backend: ya se puede consultar.
      if (writerQ.closed()) writer->setFinalDigitizerInfo(backend->info());
      // Cierre: se espera flag_delay_s (con lecturas HV nuevas) para fijar los bits del final.
      const SteadyTime until = clock.now_steady() + std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                          std::chrono::duration<double>(cfg.flagDelayS + 1.0));
      while (clock.now_steady() < until && writerOk) {
        drainSide();
        writer->tick(false);
        clock.sleep_for(std::chrono::milliseconds(100));
      }
      drainSide();
      if (auto s = writer->finish(writerOk.load()); !s) {
        logError("cierre del run: ", s.error().message);
        writerOk = false;
      }
    } catch (const std::exception& e) {
      logError("escritor: excepción ", e.what());
      writerOk = false;
    }
    writerQ.close();  // si el escritor terminó por un error, la lectura no queda bloqueada en push
    writerDone = true;
  });
  Status readoutResult = okStatus();
  JThread readoutThread([&](StopToken st) { readoutResult = readout.run(st); });

  SteadyTime nextPage = clock.now_steady();
  SteadyTime deadline{};
  bool stopping = false;
  while (!writerDone.load()) {
    if (!stopping && (gStop.load() || rstats.finished.load())) {
      logInfo(gStop.load() ? "señal recibida: parada ordenada" : "lectura terminada: cierre del run");
      readoutThread.request_stop();
      stopping = true;
      deadline = clock.now_steady() + std::chrono::duration_cast<std::chrono::nanoseconds>(
                                          std::chrono::duration<double>(cfg.writerShutdownTimeoutS));
    }
    if (page && clock.now_steady() >= nextPage) {
      const nlohmann::json status{{"run_name", run.baseName()},
                                  {"simulated", cfg.simulated()},
                                  {"readout", rstats.toJson()},
                                  {"writer", wstats.toJson()},
                                  {"writer_queue", writerQ.size()},
                                  {"writer_capacity", writerQ.capacity()},
                                  {"monitor_queue", monitorQ.size()},
                                  {"monitor_capacity", monitorQ.capacity()}};
      page->update(monitor.view().first, hvClient.latest().first, status);
      nextPage = clock.now_steady() + std::chrono::seconds(1);
    }
    if (page) page->processRequests();
    clock.sleep_for(std::chrono::milliseconds(66));  // ~15 Hz
    if (stopping) {
      if (clock.now_steady() > deadline) {
        logError("el escritor no terminó en ", cfg.writerShutdownTimeoutS, " s: salida con el índice incompleto");
        std::_Exit(4);
      }
    }
  }
  readoutThread.request_stop();
  readoutThread.join();
  monitorThread.join();
  writerThread.join();
  hvThread.request_stop();
  hvQ.close();
  hvThread.join();
  if (notifyThread.joinable()) {
    notifyQ.close();
    notifyThread.request_stop();
    notifyThread.join();
  }
  page.reset();
  logInfo("daqrun termina: ", wstats.eventsWritten.load(), " eventos, ", wstats.segmentsClosed.load(), " segmentos, ",
          wstats.bytesClosed.load(), " bytes; lectura: ", rstats.toJson().dump());
  if (!readoutResult) return 1;
  return writerOk.load() ? 0 : 1;
}
