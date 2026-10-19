// notes.cpp
#include "daq/daqrun/notes.h"

#include <sys/stat.h>

#include <algorithm>
#include <ctime>

#include <nlohmann/json.hpp>

#include "daq/core/gzip.h"

namespace daq::daqrun {

namespace {

std::string dayKey(UtcNs t) {
  const std::time_t s = static_cast<std::time_t>(t / 1000000000LL);
  std::tm tm{};
  gmtime_r(&s, &tm);
  char b[16];
  std::snprintf(b, sizeof b, "%04d%02d%02d", (tm.tm_year + 1900) % 10000, (tm.tm_mon + 1) % 100, tm.tm_mday % 100);
  return b;
}

bool exists(const std::string& p) {
  struct stat st {};
  return ::stat(p.c_str(), &st) == 0;
}

}  // namespace

std::vector<Note> readNotes(const std::string& dir, UtcNs fromUtc, UtcNs toUtc) {
  std::vector<Note> out;
  if (toUtc < fromUtc) return out;
  constexpr UtcNs kDay = 86400LL * 1000000000LL;
  std::vector<std::string> days;
  for (UtcNs t = fromUtc - (fromUtc % kDay); t <= toUtc; t += kDay) days.push_back(dayKey(t));
  for (const auto& d : days) {
    const std::string base = dir + "/" + d.substr(0, 4) + "/notes_" + d;
    for (int part = 1; part < 100; ++part) {
      const std::string stem = part == 1 ? base : base + "." + std::to_string(part);
      bool found = false;
      for (const std::string& p : {stem + ".jsonl", stem + ".jsonl.gz"}) {
        if (!exists(p)) continue;
        found = true;
        auto lr = LineReader::open(p);
        if (!lr) continue;
        while (auto line = lr.value()->next()) {
          try {
            const auto j = nlohmann::json::parse(*line);
            const UtcNs t = j.value("t_utc_ns", UtcNs{0});
            if (t < fromUtc || t >= toUtc) continue;
            out.push_back(Note{t, j.value("user", ""), j.value("host", ""), j.value("text", "")});
          } catch (const std::exception&) {
            // línea ilegible: se ignora (la nota sigue en el archivo de slow control)
          }
        }
      }
      if (!found && part > 1) break;
    }
  }
  std::sort(out.begin(), out.end(), [](const Note& a, const Note& b) { return a.tUtcNs < b.tUtcNs; });
  return out;
}

}  // namespace daq::daqrun
