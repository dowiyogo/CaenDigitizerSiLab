// hvctl: herramienta de línea de comandos de SOLO LECTURA para la fuente HV CAEN (fase 0).
//   hvctl map
//   hvctl params --slot 3 [--ch 0]
//   hvctl read   --slot 3 [--ch 0-11] [--params A,B] [--json] [--decode]
//   hvctl bench  --slot 3 --n 100 [--ch 0-11] [--params V0Set,I0Set,VMon,IMon,Status,Pw] [--json]
// Opciones globales: --config <archivo>, --host <ip>, --credentials <archivo>,
//                    --fake <guion.json | default>
// Configuración por defecto: ~/.config/daqv2/hvctl.json si existe (ver daq/config/hvctl.example.json).
// Ningún subcomando escribe en la fuente.
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

#include <sys/stat.h>

#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/core/logger.h"
#include "daq/core/time.h"
#include "daq/hv/caen_hv_backend.h"
#include "daq/hv/credentials.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hv/status_bits.h"

using namespace daq;
using namespace daq::hv;

namespace {

constexpr int kExitOk = 0;
constexpr int kExitError = 1;
constexpr int kExitUsage = 2;

const char* kUsage =
    "Uso: hvctl [opciones] <map|params|read|bench> [opciones del subcomando]\n"
    "  map                                   mapa del crate\n"
    "  params --slot S [--ch C]              parámetros de un canal: tipo, modo, unidades\n"
    "  read   --slot S [--ch 0-11] [--params A,B] [--json] [--decode]\n"
    "                                        una lectura de todos los parámetros\n"
    "  bench  --slot S --n N [--ch 0-11] [--params V0Set,I0Set,VMon,IMon,Status,Pw] [--json]\n"
    "                                        latencia de N lecturas (ms)\n"
    "Opciones: --config F  --host IP  --credentials F  --fake <guion.json|default>\n"
    "Solo lectura: ningún subcomando escribe en la fuente.\n";

struct Options {
  std::string command;
  std::string configPath;
  std::optional<std::string> host, credentials, fake, chSpec, paramSpec;
  std::optional<int> slot, ch, n;
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
      err = "argumento inesperado: " + a;
      return std::nullopt;
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

nlohmann::json readingJson(const ChannelReading& r, const StatusBitMap& bits) {
  nlohmann::json params = nlohmann::json::object();
  for (std::size_t p = 0; p < r.params.size(); ++p) {
    const auto& info = r.params[p];
    nlohmann::json e{{"type", toString(info.type)}};
    if (info.type == ParamType::Numeric) {
      nlohmann::json si = nlohmann::json::array(), board = nlohmann::json::array();
      for (const auto& v : r.values[p]) {
        si.push_back(v.si);
        board.push_back(v.boardValue);
      }
      e["unit_si"] = info.siUnit();
      e["unit_board"] = info.boardUnit();
      e["exp"] = info.exponent;
      e["unit_from_board"] = info.unitFromBoard;
      e["values_si"] = si;
      e["values_board"] = board;
    } else {
      nlohmann::json vals = nlohmann::json::array();
      for (const auto& v : r.values[p]) vals.push_back(v.raw);
      e["values"] = vals;
      if (info.type == ParamType::ChStatus) {
        nlohmann::json dec = nlohmann::json::array();
        for (const auto& v : r.values[p]) dec.push_back(bits.decode(v.raw));
        e["decoded"] = dec;
      }
    }
    params[info.name] = e;
  }
  return nlohmann::json{{"slot", r.slot},
                        {"channels", r.channels},
                        {"t_req_utc_ns", r.tReqUtcNs},
                        {"t_resp_utc_ns", r.tRespUtcNs},
                        {"t_req_utc", formatUtcIso(r.tReqUtcNs)},
                        {"duration_ns", r.durationNs},
                        {"params", params}};
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
    auto j = readingJson(r.value(), cfg.statusBits);
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

}  // namespace

int main(int argc, char** argv) {
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
  if (o.command != "map" && o.command != "params" && o.command != "read" && o.command != "bench") {
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

  SystemClock clock;
  std::unique_ptr<HVBackend> be;
  if (o.fake) {
    auto f = FakeHVBackend::fromFile(*o.fake, clock);
    if (!f) {
      logError(f.error().message);
      return kExitError;
    }
    be = std::move(f).value();
  } else {
    auto cred = loadCredentials(c.credentials);
    if (!cred) {
      logError(cred.error().message);
      return kExitError;
    }
    be = std::make_unique<CaenHVBackend>(ConnectionConfig{c.system, c.host, cred.value()}, &clock);
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
