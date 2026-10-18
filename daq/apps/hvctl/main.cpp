// hvctl: herramienta de línea de comandos para la fuente HV CAEN.
// Con hvmon corriendo (tiene el lock $XDG_RUNTIME_DIR/daqv2/hvmon.lock) todo pasa por su socket y
// hvctl nunca abre una sesión de la Wrapper:
//   hvctl read [--json] [--decode]          última lectura de hvmon (misma tabla que hv_readall)
//   hvctl status [--json]                   estado de hvmon
//   hvctl changes [--n 20] [--json]         últimos cambios
//   hvctl note "<texto>"                    nota del operador
//   hvctl measure-offsets --seconds N [--ch 0-11] [--json]   solo canales con Pw = 0
//   hvctl pw --ch N on|off --authorize "<motivo>"   única escritura; hvmon la rechaza si
//                                                   allow_commands = false
// Sin hvmon, modo directo de solo lectura (toma el mismo lock mientras dura):
//   hvctl map | params --slot 3 [--ch 0] | read --slot 3 [--ch 0-11] [--params A,B] | bench --slot 3 --n 100
// Opciones globales: --config <archivo>, --host <ip>, --credentials <json heredado>,
//                    --fake <guion.json | default> (modo directo contra el fake), --socket <ruta>
// Credenciales del modo directo: DAQ_HV_USER/DAQ_HV_PASSWORD (entorno o secrets.env); el JSON
// heredado ~/.config/daqv2/hv_credentials.json solo como respaldo obsoleto.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <sys/resource.h>
#include <sys/stat.h>

#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/core/instance_lock.h"
#include "daq/core/logger.h"
#include "daq/core/secrets.h"
#include "daq/core/unix_socket.h"
#include "daq/core/time.h"
#include "daq/hv/caen_hv_backend.h"
#include "daq/hv/credentials.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hv/reading_json.h"
#include "daq/hv/status_bits.h"
#include "daq/hvmon/config.h"

using namespace daq;
using namespace daq::hv;

namespace {

constexpr int kExitOk = 0;
constexpr int kExitError = 1;
constexpr int kExitUsage = 2;

const char* kUsage =
    "Uso: hvctl [opciones] <subcomando> [opciones del subcomando]\n"
    "Con hvmon corriendo (vía socket):\n"
    "  read [--json] [--decode]              última lectura de hvmon\n"
    "  status [--json]                       estado de hvmon\n"
    "  changes [--n 20] [--json]             últimos cambios\n"
    "  note \"<texto>\"                       nota del operador\n"
    "  measure-offsets --seconds N [--ch 0-11] [--json]   offsets con Pw = 0\n"
    "  pw --ch N on|off --authorize \"<motivo>\"          escritura de Pw (requiere allow_commands)\n"
    "Sin hvmon (modo directo, solo lectura):\n"
    "  map                                   mapa del crate\n"
    "  params --slot S [--ch C]              parámetros de un canal: tipo, modo, unidades\n"
    "  read   --slot S [--ch 0-11] [--params A,B] [--json] [--decode]\n"
    "  bench  --slot S --n N [--ch 0-11] [--params V0Set,I0Set,VMon,IMon,Status,Pw] [--json]\n"
    "Opciones: --config F  --host IP  --credentials F  --fake <guion.json|default>  --socket RUTA\n";

struct Options {
  std::string command;
  std::vector<std::string> args;  // posicionales tras el subcomando (note, pw)
  std::string configPath;
  std::optional<std::string> host, credentials, fake, chSpec, paramSpec, authorize, socket;
  std::optional<int> slot, ch, n;
  std::optional<double> seconds;
  bool json = false;
  bool decode = false;
};

struct Config {
  std::string system = "SY5527";
  std::string host = "192.168.1.118";
  std::string credentials = defaultCredentialsPath();
  int slot = 3;
  StatusBitMap statusBits = StatusBitMap::wrapperRev13();
};

bool fileExists(const std::string& p) {
  struct stat st {};
  return ::stat(p.c_str(), &st) == 0;
}

std::optional<int> toInt(const std::string& s) {
  try {
    std::size_t pos = 0;
    const int v = std::stoi(s, &pos);
    if (pos != s.size()) return std::nullopt;
    return v;
  } catch (...) {
    return std::nullopt;
  }
}

// "0-11", "0,2,5", "3" -> lista de canales.
std::optional<std::vector<int>> parseChannels(const std::string& spec) {
  std::vector<int> out;
  std::stringstream ss(spec);
  std::string part;
  while (std::getline(ss, part, ',')) {
    const auto dash = part.find('-');
    if (dash == std::string::npos) {
      auto v = toInt(part);
      if (!v || *v < 0) return std::nullopt;
      out.push_back(*v);
    } else {
      auto a = toInt(part.substr(0, dash)), b = toInt(part.substr(dash + 1));
      if (!a || !b || *a < 0 || *b < *a) return std::nullopt;
      for (int c = *a; c <= *b; ++c) out.push_back(c);
    }
  }
  if (out.empty()) return std::nullopt;
  return out;
}

std::vector<std::string> splitList(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string part;
  while (std::getline(ss, part, ','))
    if (!part.empty()) out.push_back(part);
  return out;
}

std::optional<Options> parseArgs(int argc, char** argv, std::string& err) {
  Options o;
  auto need = [&](int& i) -> std::optional<std::string> {
    if (i + 1 >= argc) {
      err = std::string("falta el valor de ") + argv[i];
      return std::nullopt;
    }
    return std::string(argv[++i]);
  };
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto intOpt = [&](std::optional<int>& dst) {
      auto v = need(i);
      if (!v) return false;
      auto iv = toInt(*v);
      if (!iv) {
        err = "valor no entero para " + a + ": " + *v;
        return false;
      }
      dst = *iv;
      return true;
    };
    auto strOpt = [&](std::optional<std::string>& dst) {
      auto v = need(i);
      if (!v) return false;
      dst = *v;
      return true;
    };
    if (a == "--help" || a == "-h") {
      o.command = "help";
    } else if (a == "--config") {
      auto v = need(i);
      if (!v) return std::nullopt;
      o.configPath = *v;
    } else if (a == "--host") {
      if (!strOpt(o.host)) return std::nullopt;
    } else if (a == "--credentials") {
      if (!strOpt(o.credentials)) return std::nullopt;
    } else if (a == "--fake") {
      if (!strOpt(o.fake)) return std::nullopt;
    } else if (a == "--slot") {
      if (!intOpt(o.slot)) return std::nullopt;
    } else if (a == "--ch") {
      if (!strOpt(o.chSpec)) return std::nullopt;
    } else if (a == "--params") {
      if (!strOpt(o.paramSpec)) return std::nullopt;
    } else if (a == "--n") {
      if (!intOpt(o.n)) return std::nullopt;
    } else if (a == "--authorize") {
      if (!strOpt(o.authorize)) return std::nullopt;
    } else if (a == "--socket") {
      if (!strOpt(o.socket)) return std::nullopt;
    } else if (a == "--seconds") {
      auto v = need(i);
      if (!v) return std::nullopt;
      try {
        o.seconds = std::stod(*v);
      } catch (...) {
        err = "valor no numérico para --seconds: " + *v;
        return std::nullopt;
      }
    } else if (a == "--json") {
      o.json = true;
    } else if (a == "--decode") {
      o.decode = true;
    } else if (!a.empty() && a[0] == '-') {
      err = "opción desconocida: " + a;
      return std::nullopt;
    } else if (o.command.empty()) {
      o.command = a;
    } else {
      o.args.push_back(a);
    }
  }
  return o;
}

Result<Config> loadConfig(const std::string& explicitPath) {
  Config c;
  std::string path = explicitPath.empty() ? expandHome("~/.config/daqv2/hvctl.json") : expandHome(explicitPath);
  if (explicitPath.empty() && !fileExists(path)) return c;
  std::ifstream in(path);
  if (!in) return Result<Config>::failure(-1, "no se pudo abrir la configuración " + path);
  nlohmann::json j;
  try {
    in >> j;
  } catch (const std::exception& e) {
    return Result<Config>::failure(-1, "configuración inválida " + path + ": " + e.what());
  }
  c.system = j.value("system", c.system);
  c.host = j.value("host", c.host);
  c.credentials = j.value("credentials", c.credentials);
  c.slot = j.value("slot", c.slot);
  if (j.contains("status_bits")) {
    auto sb = StatusBitMap::fromJson(j["status_bits"]);
    if (!sb) return Result<Config>(sb.error());
    c.statusBits = sb.value();
  }
  return c;
}

std::string fmt(const char* f, double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, f, v);
  return buf;
}

int cmdMap(HVBackend& be, bool json) {
  auto m = be.crateMap();
  if (!m) {
    logError("crateMap: ", m.error().message);
    return kExitError;
  }
  if (json) {
    nlohmann::json j = nlohmann::json::array();
    for (const auto& s : m.value().slots)
      j.push_back({{"slot", s.slot}, {"model", s.model}, {"description", s.description}, {"channels", s.channels},
                   {"serial", s.serial}, {"fw", std::to_string(s.fwMajor) + "." + std::to_string(s.fwMinor)}});
    std::cout << nlohmann::json{{"backend", be.describe()}, {"slots", j}}.dump(2) << "\n";
    return kExitOk;
  }
  for (const auto& s : m.value().slots) {
    if (s.empty()) {
      std::printf("Slot %2d: vacio\n", s.slot);
    } else {
      std::printf("Slot %2d: %-8s %-30s ch=%d SN=%u FW=%d.%02d\n", s.slot, s.model.c_str(), s.description.c_str(),
                  s.channels, s.serial, s.fwMajor, s.fwMinor);
    }
  }
  return kExitOk;
}

int cmdParams(HVBackend& be, int slot, int ch, bool json) {
  auto r = be.channelParamNames(slot, ch);
  if (!r) {
    logError("channelParamNames: ", r.error().message);
    return kExitError;
  }
  if (json) {
    nlohmann::json j = nlohmann::json::array();
    for (const auto& p : r.value()) {
      nlohmann::json e{{"name", p.name}, {"type", toString(p.type)}, {"type_code", static_cast<int>(p.type)},
                       {"mode", toString(p.mode)}};
      if (p.type == ParamType::Numeric) {
        e["unit_code"] = p.unitCode;
        e["exp"] = p.exponent;
        e["unit_board"] = p.boardUnit();
        e["unit_si"] = p.siUnit();
        e["unit_from_board"] = p.unitFromBoard;
        if (p.minBoard) e["min_board"] = *p.minBoard;
        if (p.maxBoard) e["max_board"] = *p.maxBoard;
      }
      if (p.type == ParamType::OnOff) {
        e["on_state"] = p.onState;
        e["off_state"] = p.offState;
      }
      j.push_back(e);
    }
    std::cout << nlohmann::json{{"backend", be.describe()}, {"slot", slot}, {"channel", ch}, {"params", j}}.dump(2)
              << "\n";
    return kExitOk;
  }
  std::printf("Slot %d, canal %d: %zu parametros\n\n", slot, ch, r.value().size());
  std::printf("%-10s %-9s %-7s %-6s %4s %-6s %10s %10s  %s\n", "nombre", "tipo", "modo", "unidad", "exp", "SI",
              "min", "max", "on/off o fuente de unidad");
  for (const auto& p : r.value()) {
    if (p.type == ParamType::Numeric) {
      std::printf("%-10s %-9s %-7s %-6s %4d %-6s %10s %10s  %s\n", p.name.c_str(), toString(p.type), toString(p.mode),
                  p.boardUnit().c_str(), p.exponent, p.siUnit().c_str(),
                  p.minBoard ? fmt("%.3f", *p.minBoard).c_str() : "-", p.maxBoard ? fmt("%.3f", *p.maxBoard).c_str() : "-",
                  p.unitFromBoard ? "placa" : "tabla de respaldo");
    } else {
      std::string st = p.type == ParamType::OnOff ? p.onState + "/" + p.offState : "";
      std::printf("%-10s %-9s %-7s %-6s %4s %-6s %10s %10s  %s\n", p.name.c_str(), toString(p.type), toString(p.mode),
                  "-", "-", "-", "-", "-", st.c_str());
    }
  }
  return kExitOk;
}

Result<std::vector<int>> resolveChannels(HVBackend& be, int slot, const std::optional<std::string>& spec) {
  if (spec) {
    auto c = parseChannels(*spec);
    if (!c) return Result<std::vector<int>>::failure(-1, "lista de canales inválida: " + *spec);
    return c.value();
  }
  auto m = be.crateMap();
  if (!m) return Result<std::vector<int>>(m.error());
  for (const auto& s : m.value().slots) {
    if (s.slot == slot && !s.empty()) {
      std::vector<int> all(static_cast<std::size_t>(s.channels));
      for (int i = 0; i < s.channels; ++i) all[static_cast<std::size_t>(i)] = i;
      return all;
    }
  }
  return Result<std::vector<int>>::failure(-1, "slot " + std::to_string(slot) + " vacío o inexistente");
}

void printReadingTable(const ChannelReading& r, const StatusBitMap& bits, bool decode) {
  std::printf("Slot %d: %zu canales, %zu parametros\n\n", r.slot, r.channels.size(), r.params.size());
  for (std::size_t p = 0; p < r.params.size(); ++p) {
    const auto& info = r.params[p];
    std::printf("%-10s tipo=%d |", info.name.c_str(), static_cast<int>(info.type));
    for (const auto& v : r.values[p]) {
      if (info.type == ParamType::Numeric)
        std::printf(" %9.3f", v.boardValue);
      else if (info.type == ParamType::ChStatus)
        std::printf("    0x%04X", v.raw);
      else
        std::printf(" %9u", v.raw);
    }
    std::printf("\n");
  }
  if (decode) {
    for (std::size_t p = 0; p < r.params.size(); ++p) {
      if (r.params[p].type != ParamType::ChStatus) continue;
      std::printf("\nStatus decodificado:\n");
      for (std::size_t c = 0; c < r.channels.size(); ++c) {
        auto flags = bits.decode(r.values[p][c].raw);
        std::string s;
        for (const auto& f : flags) s += (s.empty() ? "" : ",") + f;
        std::printf("  CH%02d 0x%04X %s\n", r.channels[c], r.values[p][c].raw, s.empty() ? "(ninguna)" : s.c_str());
      }
    }
  }
  std::printf("\n# t_req_utc=%s duracion=%.3f ms\n", formatUtcIso(r.tReqUtcNs).c_str(),
              static_cast<double>(r.durationNs) / 1e6);
}

int cmdRead(HVBackend& be, int slot, const Options& o, const Config& cfg) {
  auto chans = resolveChannels(be, slot, o.chSpec);
  if (!chans) {
    logError(chans.error().message);
    return kExitError;
  }
  std::vector<std::string> names;
  if (o.paramSpec) {
    names = splitList(*o.paramSpec);
  } else {
    auto info = be.channelParamNames(slot, chans.value().front());
    if (!info) {
      logError("channelParamNames: ", info.error().message);
      return kExitError;
    }
    for (const auto& p : info.value()) names.push_back(p.name);
  }
  auto r = be.readChannels(slot, chans.value(), names);
  if (!r) {
    logError("readChannels: ", r.error().message);
    return kExitError;
  }
  if (o.json) {
    auto j = readingToJson(r.value(), cfg.statusBits);
    j["backend"] = be.describe();
    std::cout << j.dump(2) << "\n";
  } else {
    printReadingTable(r.value(), cfg.statusBits, o.decode);
  }
  return kExitOk;
}

int cmdBench(HVBackend& be, int slot, const Options& o) {
  const int n = o.n.value_or(100);
  if (n <= 0) {
    logError("--n debe ser positivo");
    return kExitUsage;
  }
  auto chans = resolveChannels(be, slot, o.chSpec ? o.chSpec : std::optional<std::string>("0-11"));
  if (!chans) {
    logError(chans.error().message);
    return kExitError;
  }
  const auto names = splitList(o.paramSpec.value_or("V0Set,I0Set,VMon,IMon,Status,Pw"));
  std::vector<double> ms;
  int errors = 0;
  std::string firstError;
  const UtcNs tStart = utcNowNs();
  for (int i = 0; i < n; ++i) {
    auto r = be.readChannels(slot, chans.value(), names);
    if (!r) {
      ++errors;
      if (firstError.empty()) firstError = r.error().message;
      continue;
    }
    ms.push_back(static_cast<double>(r.value().durationNs) / 1e6);
  }
  const UtcNs tEnd = utcNowNs();
  std::vector<double> sorted = ms;
  std::sort(sorted.begin(), sorted.end());
  auto rank = [&](double q) -> double {
    if (sorted.empty()) return NAN;
    auto idx = static_cast<std::size_t>(std::ceil(q * static_cast<double>(sorted.size())));
    idx = std::clamp<std::size_t>(idx, 1, sorted.size());
    return sorted[idx - 1];
  };
  double mean = NAN;
  if (!ms.empty()) {
    double s = 0;
    for (double x : ms) s += x;
    mean = s / static_cast<double>(ms.size());
  }
  const double median = sorted.empty() ? NAN
                        : (sorted.size() % 2 ? sorted[sorted.size() / 2]
                                             : 0.5 * (sorted[sorted.size() / 2 - 1] + sorted[sorted.size() / 2]));
  std::string plist;
  for (const auto& p : names) plist += (plist.empty() ? "" : ",") + p;
  if (o.json) {
    nlohmann::json j{{"backend", be.describe()}, {"slot", slot}, {"channels", chans.value()}, {"params", names},
                     {"n", n}, {"ok", ms.size()}, {"errors", errors}, {"t_start_utc", formatUtcIso(tStart)},
                     {"t_end_utc", formatUtcIso(tEnd)}, {"min_ms", sorted.empty() ? NAN : sorted.front()},
                     {"median_ms", median}, {"p95_ms", rank(0.95)}, {"max_ms", sorted.empty() ? NAN : sorted.back()},
                     {"mean_ms", mean}, {"durations_ms", ms}};
    if (!firstError.empty()) j["first_error"] = firstError;
    std::cout << j.dump(2) << "\n";
  } else {
    std::printf("bench: slot %d, %zu canales, %zu parametros (%s), n=%d\n", slot, chans.value().size(), names.size(),
                plist.c_str(), n);
    std::printf("  inicio %s  fin %s\n", formatUtcIso(tStart).c_str(), formatUtcIso(tEnd).c_str());
    std::printf("  lecturas ok %zu, errores %d\n", ms.size(), errors);
    if (!sorted.empty()) {
      std::printf("  min    %8.3f ms\n  median %8.3f ms\n  p95    %8.3f ms\n  max    %8.3f ms\n  media  %8.3f ms\n",
                  sorted.front(), median, rank(0.95), sorted.back(), mean);
    }
    if (!firstError.empty()) std::printf("  primer error: %s\n", firstError.c_str());
  }
  return errors == 0 ? kExitOk : kExitError;
}

// ---- Cliente de hvmon -------------------------------------------------------------------

Result<nlohmann::json> request(const std::string& sock, const nlohmann::json& req, double timeoutS) {
  auto r = unixRequest(sock, req.dump(), std::chrono::milliseconds(static_cast<long>(timeoutS * 1000)));
  if (!r) return Result<nlohmann::json>(r.error());
  try {
    auto j = nlohmann::json::parse(r.value());
    if (!j.value("ok", false)) {
      std::string msg = j.value("error", "");
      if (msg.empty() && j.contains("result"))
        msg = "comando " + std::to_string(j.value("command_id", 0ULL)) + ": " + j.value("result", "") + " (" +
              j.value("detail", "") + ")";
      return Result<nlohmann::json>::failure(-1, msg.empty() ? j.dump() : msg);
    }
    return j;
  } catch (const std::exception& e) {
    return Result<nlohmann::json>::failure(-1, std::string("respuesta inválida de hvmon: ") + e.what());
  }
}

int fail(const Error& e) {
  logError(e.message);
  return kExitError;
}

int cmdSocketRead(const std::string& sock, const Options& o, const Config& cfg) {
  auto r = request(sock, {{"cmd", "read"}}, 10);
  if (!r) return fail(r.error());
  const auto& j = r.value()["reading"];
  if (o.json) {
    std::cout << j.dump(2) << "\n";
    return kExitOk;
  }
  auto rd = readingFromJson(j);
  if (!rd) return fail(rd.error());
  printReadingTable(rd.value(), cfg.statusBits, o.decode);
  std::printf("# fuente: hvmon (sin sesión propia de la Wrapper), conectado=%s\n",
              j.value("connected", false) ? "si" : "NO");
  return kExitOk;
}

int cmdStatus(const std::string& sock, bool json) {
  auto r = request(sock, {{"cmd", "status"}}, 10);
  if (!r) return fail(r.error());
  const auto& s = r.value()["state"];
  if (json) {
    std::cout << s.dump(2) << "\n";
    return kExitOk;
  }
  const auto& c = s["connection"];
  const auto& st = s["stats"];
  const auto& w = s["writer"];
  std::printf("hvmon pid %ld  %s\n", s.value("pid", 0L), s.value("t_utc", "").c_str());
  std::printf("conexión: conectado=%s stale=%s edad=%.3f s reconexiones=%llu último error: %s\n",
              c.value("connected", false) ? "si" : "NO", c.value("stale", true) ? "SI" : "no", c.value("age_s", -1.0),
              static_cast<unsigned long long>(st.value("reconnections", 0ULL)), c.value("last_error", "").c_str());
  std::printf("lecturas: rápidas %llu, completas %llu, errores %llu, media %.3f ms, máx %.3f ms\n",
              static_cast<unsigned long long>(st.value("fast_reads", 0ULL)),
              static_cast<unsigned long long>(st.value("full_reads", 0ULL)),
              static_cast<unsigned long long>(st.value("read_errors", 0ULL)), st.value("fast_read_mean_ms", 0.0),
              st.value("fast_read_max_ms", 0.0));
  std::printf("escritor: ok=%s fsync=%llu bytes=%llu errores=%llu descartados=%llu %s\n", w.value("ok", false) ? "si" : "NO",
              static_cast<unsigned long long>(w.value("fsyncs", 0ULL)), static_cast<unsigned long long>(w.value("bytes", 0ULL)),
              static_cast<unsigned long long>(w.value("errors", 0ULL)), static_cast<unsigned long long>(w.value("dropped", 0ULL)),
              w.value("last_error", "").c_str());
  std::printf("GECO2020 en ejecución: %s   allow_commands: %s\n\n", s.value("geco_running", false) ? "si" : "no",
              s.value("allow_commands", false) ? "si" : "no");
  std::printf("CH   V0Set[V]   I0Set[A]    VMon[V]    IMon[A]  Pw  Status  flags\n");
  for (const auto& ch : s["channels"]) {
    const auto& v = ch["values"];
    std::string flags;
    for (const auto& f : ch["status_flags"]) flags += (flags.empty() ? "" : ",") + f.get<std::string>();
    std::printf("%02d %9.3f %10.3e %9.3f %10.3e  %2u  0x%04X  %s\n", ch.value("ch", -1), v.value("V0Set_V", 0.0),
                v.value("I0Set_A", 0.0), v.value("VMon_V", 0.0), v.value("IMon_A", 0.0), v.value("Pw", 0u),
                v.value("Status", 0u), flags.empty() ? "-" : flags.c_str());
  }
  return kExitOk;
}

int cmdChanges(const std::string& sock, int n, bool json) {
  auto r = request(sock, {{"cmd", "changes"}, {"n", n}}, 10);
  if (!r) return fail(r.error());
  const auto& ch = r.value()["changes"];
  if (json) {
    std::cout << ch.dump(2) << "\n";
    return kExitOk;
  }
  if (ch.empty()) std::printf("(sin cambios desde que arrancó hvmon)\n");
  for (const auto& c : ch)
    std::printf("%s  CH%02d %-8s %s -> %s %s  origen=%s%s\n", c.value("t_utc", "").c_str(), c.value("ch", -1),
                c.value("param", "").c_str(), c["old"].dump().c_str(), c["new"].dump().c_str(), c.value("unit", "").c_str(),
                c.value("origin", "").c_str(), c.value("during_disconnect", false) ? " (tras desconexión)" : "");
  return kExitOk;
}

int cmdNote(const std::string& sock, const Options& o) {
  if (o.args.size() != 1) {
    std::fprintf(stderr, "hvctl note: se espera un texto entre comillas\n");
    return kExitUsage;
  }
  auto r = request(sock, {{"cmd", "note"}, {"text", o.args[0]}}, 10);
  if (!r) return fail(r.error());
  std::printf("nota registrada %s\n", r.value().value("t_utc", "").c_str());
  return kExitOk;
}

int cmdOffsets(const std::string& sock, const Options& o) {
  if (!o.seconds) {
    std::fprintf(stderr, "hvctl measure-offsets: falta --seconds N\n");
    return kExitUsage;
  }
  nlohmann::json req{{"cmd", "measure-offsets"}, {"seconds", *o.seconds}};
  if (o.chSpec) {
    auto c = parseChannels(*o.chSpec);
    if (!c) {
      std::fprintf(stderr, "hvctl: lista de canales inválida\n");
      return kExitUsage;
    }
    req["channels"] = *c;
  }
  auto r = request(sock, req, *o.seconds + 90);
  if (!r) return fail(r.error());
  const auto& off = r.value()["offsets"];
  if (o.json) {
    std::cout << off.dump(2) << "\n";
    return kExitOk;
  }
  std::printf("offsets %s, %.0f s, %llu lecturas; rechazados por Pw = 1: %s\n", off.value("t_start_utc", "").c_str(),
              off.value("seconds", 0.0), static_cast<unsigned long long>(off.value("n_reads", 0ULL)),
              off["refused_pw_on"].dump().c_str());
  std::printf("CH   IMon offset [uA]   std [uA]   VMon offset [V]   std [V]\n");
  for (const auto& c : off["channels"])
    std::printf("%02d   %14.4f   %8.4f   %15.3f   %7.3f\n", c.value("ch", -1), c.value("imon_offset_A", 0.0) * 1e6,
                c.value("imon_std_A", 0.0) * 1e6, c.value("vmon_offset_V", 0.0), c.value("vmon_std_V", 0.0));
  return kExitOk;
}

int cmdPw(const std::string& sock, const Options& o) {
  if (!o.chSpec || o.args.size() != 1 || (o.args[0] != "on" && o.args[0] != "off") || !o.authorize) {
    std::fprintf(stderr, "Uso: hvctl pw --ch N on|off --authorize \"<motivo>\"\n");
    return kExitUsage;
  }
  auto ch = toInt(*o.chSpec);
  if (!ch) {
    std::fprintf(stderr, "hvctl pw: --ch debe ser un único canal\n");
    return kExitUsage;
  }
  auto r = request(sock, {{"cmd", "pw"}, {"ch", *ch}, {"on", o.args[0] == "on"}, {"authorize", *o.authorize}}, 3700);
  if (!r) return fail(r.error());
  std::printf("comando %llu: %s (%s)\n", static_cast<unsigned long long>(r.value().value("command_id", 0ULL)),
              r.value().value("result", "").c_str(), r.value().value("detail", "").c_str());
  return kExitOk;
}

}  // namespace

int main(int argc, char** argv) {
  rlimit noCore{0, 0};
  setrlimit(RLIMIT_CORE, &noCore);  // sin volcados de memoria que puedan contener la contraseña
  Logger::setThreadName("hvctl");
  std::string err;
  auto opt = parseArgs(argc, argv, err);
  if (!opt) {
    std::fprintf(stderr, "hvctl: %s\n\n%s", err.c_str(), kUsage);
    return kExitUsage;
  }
  const Options& o = *opt;
  if (o.command.empty() || o.command == "help") {
    std::fputs(kUsage, o.command == "help" ? stdout : stderr);
    return o.command == "help" ? kExitOk : kExitUsage;
  }
  const std::vector<std::string> socketOnly{"status", "changes", "note", "measure-offsets", "pw"};
  const std::vector<std::string> directOnly{"map", "params", "bench"};
  const bool isSocketOnly = std::find(socketOnly.begin(), socketOnly.end(), o.command) != socketOnly.end();
  const bool isDirectOnly = std::find(directOnly.begin(), directOnly.end(), o.command) != directOnly.end();
  if (!isSocketOnly && !isDirectOnly && o.command != "read") {
    std::fprintf(stderr, "hvctl: subcomando desconocido: %s\n\n%s", o.command.c_str(), kUsage);
    return kExitUsage;
  }
  auto cfg = loadConfig(o.configPath);
  if (!cfg) {
    logError(cfg.error().message);
    return kExitError;
  }
  Config c = cfg.value();
  if (o.host) c.host = *o.host;
  if (o.credentials) c.credentials = *o.credentials;
  const int slot = o.slot.value_or(c.slot);

  const std::string runtimeDir = daq::hvmon::defaultRuntimeDir();
  const std::string sock = o.socket.value_or(runtimeDir + "/hvmon.sock");
  const std::string lockPath = runtimeDir + "/hvmon.lock";
  const bool hvmonRunning = !o.fake && InstanceLock::isHeld(lockPath);
  if (isSocketOnly) {
    if (!hvmonRunning) {
      logError("hvmon no está corriendo (lock libre: ", lockPath, "); '", o.command, "' necesita hvmon");
      return kExitError;
    }
    if (o.command == "status") return cmdStatus(sock, o.json);
    if (o.command == "changes") return cmdChanges(sock, o.n.value_or(20), o.json);
    if (o.command == "note") return cmdNote(sock, o);
    if (o.command == "measure-offsets") return cmdOffsets(sock, o);
    return cmdPw(sock, o);
  }
  if (hvmonRunning) {
    if (o.command == "read") return cmdSocketRead(sock, o, c);
    logError("hvmon tiene la sesión de la Wrapper; '", o.command, "' no está disponible mientras corre (use read/status)");
    return kExitError;
  }

  SystemClock clock;
  std::unique_ptr<HVBackend> be;
  InstanceLock directLock;
  if (o.fake) {
    auto f = FakeHVBackend::fromFile(*o.fake, clock);
    if (!f) {
      logError(f.error().message);
      return kExitError;
    }
    be = std::move(f).value();
  } else {
    // Modo directo: tomar el lock de dueño de la Wrapper mientras dure (hvmon no puede arrancar).
    if (auto st = ensurePrivateDirectory(runtimeDir, 0700); !st) {
      logError(st.error().message);
      return kExitError;
    }
    auto l = InstanceLock::acquire(lockPath);
    if (!l) {
      logError("no se puede usar la Wrapper: ", l.error().message);
      return kExitError;
    }
    directLock = std::move(l).value();
    Secrets secrets;
    auto login = resolveHvLogin(secrets, c.host, c.credentials);
    if (!login) {
      logError(login.error().message);
      return kExitError;
    }
    const std::string host = o.host.value_or(login.value().host);
    be = std::make_unique<CaenHVBackend>(ConnectionConfig{c.system, host, login.value().credentials}, &clock);
  }
  if (auto st = be->connect(); !st) {
    logError("connect: ", st.error().message);
    return kExitError;
  }

  int rc = kExitOk;
  if (o.command == "map") rc = cmdMap(*be, o.json);
  else if (o.command == "params") rc = cmdParams(*be, slot, o.chSpec ? toInt(*o.chSpec).value_or(0) : 0, o.json);
  else if (o.command == "read") rc = cmdRead(*be, slot, o, c);
  else if (o.command == "bench") rc = cmdBench(*be, slot, o);
  be->disconnect();
  return rc;
}
