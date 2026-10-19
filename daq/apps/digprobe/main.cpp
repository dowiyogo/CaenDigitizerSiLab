// digprobe: pruebas del VX1740 fuera de un run (fase 4, DAQ2-6; requiere autorización del operador).
//   digprobe info [--config daqrun.json]       abre, GetInfo, versiones, volcado de registros, cierra
//                                              (sin Reset ni configuración: lo que la placa tiene)
//   digprobe memfull --config daqrun.json [--triggers N] [--rate HZ]
//                                              configura (SW trigger, ventana de la config), arranca,
//                                              envía N SW triggers SIN leer, vigila 0x8104 y 0x812C,
//                                              luego lee todo: eventos, contadores, saltos; para y cierra
// Salida JSON por stdout; log por stderr.
#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <thread>

#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/core/logger.h"
#include "daq/dig/caen_digitizer.h"

using namespace daq;
using namespace daq::dig;

namespace {
Result<DigitizerConfig> loadCfg(const std::string& path) {
  if (path.empty()) return DigitizerConfig::fromJson({{"backend", "caen"}});
  std::ifstream in(path);
  if (!in) return Result<DigitizerConfig>::failure(-1, "no se pudo abrir " + path);
  nlohmann::json j;
  in >> j;
  nlohmann::json d = j.value("digitizer", nlohmann::json::object());
  d["backend"] = "caen";
  return DigitizerConfig::fromJson(d);
}
}  // namespace

int main(int argc, char** argv) {
  Logger::setThreadName("digprobe");
  std::string cmd = argc > 1 ? argv[1] : "", cfgPath;
  int triggers = 1500;
  double rate = 200;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--config" && i + 1 < argc) cfgPath = argv[++i];
    else if (a == "--triggers" && i + 1 < argc) triggers = std::atoi(argv[++i]);
    else if (a == "--rate" && i + 1 < argc) rate = std::atof(argv[++i]);
  }
  if (cmd != "info" && cmd != "memfull") {
    std::fprintf(stderr, "Uso: digprobe info|memfull [--config archivo.json] [--triggers N] [--rate HZ]\n");
    return 2;
  }
  auto cfg = loadCfg(cfgPath);
  if (!cfg) {
    std::fprintf(stderr, "digprobe: %s\n", cfg.error().message.c_str());
    return 2;
  }
  SystemClock clock;
  CaenDigitizerBackend dig(clock);
  if (auto s = dig.open(cfg.value()); !s) {
    std::fprintf(stderr, "digprobe: %s\n", s.error().message.c_str());
    return 1;
  }
  nlohmann::json out;
  if (cmd == "info") {
    out["info"] = dig.info();
    out["registers"] = dig.dumpRegisters(CaenDigitizerBackend::standardRegisters());
    dig.close();
    std::printf("%s\n", out.dump(2).c_str());
    return 0;
  }
  // memfull: SW trigger, ext deshabilitado.
  DigitizerConfig c = cfg.value();
  c.softwareTrigger = true;
  c.externalTrigger = false;
  if (auto s = dig.configure(c); !s) {
    std::fprintf(stderr, "digprobe: %s\n", s.error().message.c_str());
    return 1;
  }
  out["config"] = dig.info();
  auto bo = dig.readRegister(0x800C);
  out["buffer_organization_0x800C"] = bo ? bo.value() : 0u;
  if (!dig.start()) return 1;
  nlohmann::json samples = nlohmann::json::array();
  const auto period = std::chrono::duration<double>(1.0 / rate);
  int sent = 0, sendErrors = 0;
  for (int i = 0; i < triggers; ++i) {
    if (!dig.sendSoftwareTrigger()) ++sendErrors;
    else ++sent;
    if (i % 50 == 49 || i == triggers - 1) {
      auto st = dig.readRegister(0x8104);
      auto ev = dig.readRegister(0x812C);
      samples.push_back({{"triggers_sent", sent}, {"status_0x8104", st ? st.value() : 0xFFFFFFFFu},
                         {"event_stored_0x812C", ev ? ev.value() : 0xFFFFFFFFu}});
    }
    std::this_thread::sleep_for(period);
  }
  out["while_not_reading"] = samples;
  out["triggers_sent"] = sent;
  out["send_errors"] = sendErrors;
  // Vaciado: se leen todos los eventos.
  dig.stop();
  std::uint64_t events = 0, jumps = 0;
  std::uint32_t first = 0, last = 0;
  bool have = false, full = false;
  for (int k = 0; k < 100000; ++k) {
    EventBlock b;
    auto r = dig.readBlock(b, 1000);
    if (!r || r.value() == 0) break;
    full = full || b.boardFull;
    for (const auto& e : b.events) {
      if (!have) first = e.counterHw;
      else if (((e.counterHw - last) & 0xFFFFFF) != 1) ++jumps;
      last = e.counterHw;
      have = true;
      ++events;
    }
  }
  auto st = dig.readRegister(0x8104);
  out["after_drain_status_0x8104"] = st ? st.value() : 0u;
  out["events_read"] = events;
  out["first_counter"] = first;
  out["last_counter"] = last;
  out["counter_jumps"] = jumps;
  out["board_full_seen_on_first_read"] = full;
  out["final_info"] = dig.info();
  dig.close();
  std::printf("%s\n", out.dump(2).c_str());
  return 0;
}
