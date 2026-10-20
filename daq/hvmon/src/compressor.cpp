// compressor.cpp
#include "daq/hvmon/compressor.h"

#include <sys/stat.h>

#include <dirent.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <algorithm>
#include <cctype>

#include "daq/core/logger.h"

namespace daq::hvmon {

namespace {

bool endsWith(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

// "hv_20261018.csv" o "hv_full_20261018.2.csv" (parte con otras columnas) -> "20261018";
// "" si no tiene la forma <tipo>_YYYYMMDD[.N].<csv|jsonl>.
std::string dateOf(const std::string& name) {
  auto dot = name.rfind('.');
  if (dot == std::string::npos || dot < 9) return "";
  const std::string ext = name.substr(dot);
  if (ext != ".csv" && ext != ".jsonl") return "";
  const auto prev = name.rfind('.', dot - 1);
  if (prev != std::string::npos && prev + 1 < dot &&
      std::all_of(name.begin() + static_cast<std::ptrdiff_t>(prev) + 1, name.begin() + static_cast<std::ptrdiff_t>(dot),
                  [](unsigned char c) { return std::isdigit(c); }))
    dot = prev;  // parte .N
  if (dot < 9) return "";
  const std::string d = name.substr(dot - 8, 8);
  if (name[dot - 9] != '_' || !std::all_of(d.begin(), d.end(), [](unsigned char c) { return std::isdigit(c); }))
    return "";
  return d;
}

}  // namespace

Compressor::Compressor(StateHub& hub) : hub_(hub) { hub_.compressor.publish(status()); }

std::vector<std::string> Compressor::closedDailyFiles(const std::string& baseDir, const std::string& todayKey) {
  std::vector<std::string> out;
  DIR* top = ::opendir(baseDir.c_str());
  if (!top) return out;
  while (dirent* y = ::readdir(top)) {
    const std::string year = y->d_name;
    if (year.size() != 4 || !std::all_of(year.begin(), year.end(), [](unsigned char c) { return std::isdigit(c); }))
      continue;
    const std::string dir = baseDir + "/" + year;
    DIR* d = ::opendir(dir.c_str());
    if (!d) continue;
    while (dirent* f = ::readdir(d)) {
      const std::string name = f->d_name;
      const std::string date = dateOf(name);
      if (!date.empty() && date < todayKey && !endsWith(name, ".schema.json")) out.push_back(dir + "/" + name);
    }
    ::closedir(d);
  }
  ::closedir(top);
  std::sort(out.begin(), out.end());
  return out;
}

Result<GzipResult> Compressor::compressOne(const std::string& path) {
  // Encolado dos veces (rotación del escritor y búsqueda de días cerrados): ya está comprimido.
  struct stat sb {};
  if (::stat(path.c_str(), &sb) != 0 && ::stat((path + ".gz").c_str(), &sb) == 0) {
    ++skipped_;
    logDebug("compresión: ", path, " ya comprimido, se omite");
    hub_.compressor.publish(status());
    return GzipResult{};
  }
  auto r = gzipFileVerified(path, path + ".gz", /*removeSource=*/true);
  if (r) {
    ++done_;
    inBytes_ += r.value().source.bytes;
    outBytes_ += r.value().compressedBytes;
    last_ = {{"file", path}, {"ok", true}, {"bytes", r.value().source.bytes}, {"gz_bytes", r.value().compressedBytes},
             {"lines", r.value().source.lines}, {"ratio", r.value().ratio()}};
    logInfo("comprimido y verificado: ", path, " (", r.value().source.bytes, " -> ", r.value().compressedBytes, " bytes, ",
            r.value().source.lines, " líneas)");
  } else {
    ++failed_;
    last_ = {{"file", path}, {"ok", false}, {"error", r.error().message}};
    logError("compresión de ", path, ": ", r.error().message);
  }
  hub_.compressor.publish(status());
  return r;
}

nlohmann::json Compressor::status() const {
  return nlohmann::json{{"compressed", done_}, {"failed", failed_}, {"skipped", skipped_}, {"in_bytes", inBytes_}, {"gz_bytes", outBytes_},
                        {"last", last_}};
}

void Compressor::run(BoundedQueue<std::string>& queue, const StopToken& st) {
  Logger::setThreadName("gzip");
  // Prioridad baja solo para este hilo (en Linux setpriority acepta el id del hilo).
  ::setpriority(PRIO_PROCESS, static_cast<id_t>(::syscall(SYS_gettid)), 19);
  while (!st.stop_requested()) {
    auto p = queue.pop_until(std::chrono::steady_clock::now() + std::chrono::milliseconds(200));
    if (p) compressOne(*p);
    else if (queue.closed()) break;
  }
}

}  // namespace daq::hvmon
