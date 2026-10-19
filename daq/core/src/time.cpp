// time.cpp
#include "daq/core/time.h"

#include <sys/stat.h>

#include <cstdio>
#include <cstdlib>
#include <ctime>

namespace daq {

UtcNs utcNowNs() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

SteadyTime steadyNow() { return std::chrono::steady_clock::now(); }

namespace {

// Separa en segundos (piso) y nanosegundos en [0, 1e9), válido también para t < 0.
void splitNs(UtcNs t, std::time_t& sec, long& nsec) {
  constexpr UtcNs kNs = 1000000000;
  UtcNs s = t / kNs;
  UtcNs r = t % kNs;
  if (r < 0) {
    r += kNs;
    --s;
  }
  sec = static_cast<std::time_t>(s);
  nsec = static_cast<long>(r);
}

}  // namespace

std::string formatUtcCompact(UtcNs t) {
  std::time_t sec;
  long nsec;
  splitNs(t, sec, nsec);
  std::tm tm{};
  gmtime_r(&sec, &tm);
  char buf[80];
  std::snprintf(buf, sizeof buf, "%04d%02d%02dT%02d%02d%02dZ", tm.tm_year + 1900, tm.tm_mon + 1,
                tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  return buf;
}

std::string formatUtcIso(UtcNs t) {
  std::time_t sec;
  long nsec;
  splitNs(t, sec, nsec);
  std::tm tm{};
  gmtime_r(&sec, &tm);
  char buf[96];
  std::snprintf(buf, sizeof buf, "%04d-%02d-%02dT%02d:%02d:%02d.%09ldZ", tm.tm_year + 1900,
                tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, nsec);
  return buf;
}

Status setProcessTimeZone(const std::string& tz) {
  if (tz.empty() || tz.find("..") != std::string::npos || tz[0] == '/')
    return Status::failure(-1, "zona horaria inválida: '" + tz + "'");
  struct stat st {};
  if (::stat(("/usr/share/zoneinfo/" + tz).c_str(), &st) != 0 || !S_ISREG(st.st_mode))
    return Status::failure(-1, "zona horaria desconocida: " + tz + " (no está en /usr/share/zoneinfo)");
  if (::setenv("TZ", tz.c_str(), 1) != 0) return Status::failure(-1, "setenv TZ falló");
  ::tzset();
  return okStatus();
}

std::string formatLocalShort(UtcNs t) {
  std::time_t sec;
  long nsec;
  splitNs(t, sec, nsec);
  std::tm tm{};
  localtime_r(&sec, &tm);
  char buf[64];
  std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S %Z", &tm);
  return buf;
}

}  // namespace daq
