// campaign.cpp
#include "daq/daqrun/campaign.h"

#include <dirent.h>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "daq/core/file_system.h"
#include "daq/core/sha256.h"
#include "daq/daqrun/config.h"

namespace daq::daqrun {

namespace {

constexpr const char* kOpenMarker = "CAMPAIGN_OPEN";
constexpr const char* kClosedMarker = "CAMPAIGN_CLOSED";
constexpr const char* kRunLock = ".run.lock";
constexpr const char* kSums = "SHA256SUMS";

bool exists(const std::string& p) {
  struct stat st {};
  return ::stat(p.c_str(), &st) == 0;
}

std::string hostName() {
  char h[256] = {0};
  if (::gethostname(h, sizeof h - 1) != 0) return "unknown";
  return h;
}

Result<nlohmann::json> readJson(const std::string& path) {
  PosixFileSystem fs;
  auto t = fs.readFile(path);
  if (!t) return Result<nlohmann::json>(t.error());
  try {
    return nlohmann::json::parse(t.value());
  } catch (const std::exception& e) {
    return Result<nlohmann::json>::failure(-1, "JSON inválido en " + path + ": " + e.what());
  }
}

// Lock exclusivo bloqueante de un archivo auxiliar (manifiesto, SHA256SUMS); se libera al cerrar.
class FileLock {
 public:
  explicit FileLock(const std::string& path) {
    fd_ = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
    if (fd_ >= 0 && ::flock(fd_, LOCK_EX) != 0) {
      ::close(fd_);
      fd_ = -1;
    }
  }
  ~FileLock() {
    if (fd_ >= 0) ::close(fd_);
  }
  bool ok() const { return fd_ >= 0; }

 private:
  int fd_ = -1;
};

// Número de run de un nombre "<campaña>_<ts>_runNNN..." (0 si no coincide).
int runNumberOf(const std::string& campaign, const std::string& name) {
  if (name.rfind(campaign + "_", 0) != 0) return 0;
  const auto p = name.find("_run", campaign.size());
  if (p == std::string::npos) return 0;
  int n = 0, digits = 0;
  for (std::size_t i = p + 4; i < name.size() && std::isdigit(static_cast<unsigned char>(name[i])); ++i, ++digits)
    n = n * 10 + (name[i] - '0');
  return digits >= 3 ? n : 0;
}

}  // namespace

CampaignState campaignState(const CampaignLocation& loc) {
  const std::string d = loc.dir();
  struct stat st {};
  if (::stat(d.c_str(), &st) != 0) return CampaignState::Missing;
  if (!S_ISDIR(st.st_mode)) return CampaignState::Invalid;
  const bool open = exists(d + "/" + kOpenMarker), closed = exists(d + "/" + kClosedMarker);
  if (closed) return CampaignState::Closed;
  if (open) return CampaignState::Open;
  return CampaignState::Invalid;
}

Status openCampaign(const CampaignLocation& loc, const std::string& operatorName, bool simulated, UtcNs now) {
  if (!validCampaignName(loc.name)) return Status::failure(-1, "nombre de campaña inválido: " + loc.name);
  if (!simulated && (pathIsWithin(loc.rawRoot, "/data/daq/sim") || pathIsWithin(loc.manifestsDir, "/data/daq/sim")))
    return Status::failure(-1, "una campaña real no va en /data/daq/sim");
  if (simulated && (pathIsWithin(loc.rawRoot, kRealRawRoot) || pathIsWithin(loc.manifestsDir, kRealManifests)))
    return Status::failure(-1, std::string("una campaña simulada no puede ir en ") + kRealRawRoot + " ni en " + kRealManifests);
  switch (campaignState(loc)) {
    case CampaignState::Open: return Status::failure(-1, "la campaña " + loc.name + " ya está abierta");
    case CampaignState::Closed: return Status::failure(-1, "la campaña " + loc.name + " está cerrada; no se reabre");
    case CampaignState::Invalid:
      return Status::failure(-1, loc.dir() + " existe sin marcador de campaña: no se usa");
    case CampaignState::Missing: break;
  }
  // Ningún otro directorio de raw_root puede tener el mismo nombre de campaña (único globalmente).
  PosixFileSystem fs;
  if (auto st = fs.makeDirectories(loc.rawRoot, 0775); !st) return st;
  if (auto st = fs.makeDirectories(loc.manifestsDir, 0775); !st) return st;
  if (exists(loc.manifestPath())) return Status::failure(-1, "ya existe el manifiesto " + loc.manifestPath());
  if (::mkdir(loc.dir().c_str(), 0775) != 0)
    return Status::failure(errno, "mkdir " + loc.dir() + ": " + std::strerror(errno));
  const char* user = std::getenv("USER");
  const nlohmann::json marker{{"schema", "daqv2.campaign_marker/1"},
                              {"campaign", loc.name},
                              {"opened_utc_ns", now},
                              {"opened_utc", formatUtcIso(now)},
                              {"operator", operatorName.empty() ? std::string(user ? user : "") : operatorName},
                              {"host", hostName()},
                              {"simulated", simulated},
                              {"raw_root", resolvePath(loc.rawRoot)}};
  if (auto st = fs.writeAtomic(loc.dir() + "/" + kOpenMarker, marker.dump(2) + "\n", 0444); !st) return st;
  nlohmann::json manifest = marker;
  manifest["schema"] = "daqv2.campaign_manifest/1";
  manifest["status"] = "open";
  manifest["campaign_dir"] = resolvePath(loc.dir());
  manifest["runs"] = nlohmann::json::array();
  return fs.writeAtomic(loc.manifestPath(), manifest.dump(2) + "\n", 0644);
}

Status closeCampaign(const CampaignLocation& loc, UtcNs now) {
  if (campaignState(loc) != CampaignState::Open) return Status::failure(-1, "la campaña " + loc.name + " no está abierta");
  const std::string d = loc.dir();
  if (InstanceLock::isHeld(d + "/" + kRunLock)) return Status::failure(-1, "hay un run en curso en " + loc.name);
  auto marker = readJson(d + "/" + kOpenMarker);
  if (!marker) return Status(marker.error());
  std::string sumsSha;
  if (exists(d + "/" + kSums)) {
    auto s = sha256File(d + "/" + kSums);
    if (!s) return Status(s.error());
    sumsSha = s.value();
  }
  nlohmann::json closed = marker.value();
  closed["closed_utc_ns"] = now;
  closed["closed_utc"] = formatUtcIso(now);
  closed["sha256sums_sha256"] = sumsSha;
  {
    FileLock lk(loc.manifestPath() + ".lock");
    auto man = readJson(loc.manifestPath());
    nlohmann::json m = man ? man.value() : nlohmann::json::object();
    closed["runs"] = m.value("runs", nlohmann::json::array());
    m["status"] = "closed";
    m["closed_utc_ns"] = now;
    m["closed_utc"] = formatUtcIso(now);
    m["sha256sums_sha256"] = sumsSha;
    PosixFileSystem fs;
    if (auto st = fs.writeAtomic(loc.manifestPath(), m.dump(2) + "\n", 0644); !st) return st;
  }
  // El marcador se renombra (no se borra) y luego se reescribe con el cierre.
  if (::chmod(d.c_str(), 0775) != 0) return Status::failure(errno, "chmod " + d);
  if (::rename((d + "/" + kOpenMarker).c_str(), (d + "/" + kClosedMarker).c_str()) != 0)
    return Status::failure(errno, "rename del marcador: " + std::string(std::strerror(errno)));
  PosixFileSystem fs;
  if (auto st = fs.writeAtomic(d + "/" + kClosedMarker, closed.dump(2) + "\n", 0444); !st) return st;
  if (::chmod(d.c_str(), 0555) != 0) return Status::failure(errno, "chmod a-w " + d);
  return okStatus();
}

nlohmann::json campaignStatus(const CampaignLocation& loc) {
  const char* names[] = {"missing", "open", "closed", "invalid"};
  nlohmann::json j{{"campaign", loc.name}, {"dir", loc.dir()}, {"state", names[static_cast<int>(campaignState(loc))]}};
  if (auto m = readJson(loc.manifestPath())) j["manifest"] = m.value();
  j["run_in_progress"] = InstanceLock::isHeld(loc.dir() + "/" + kRunLock);
  return j;
}

bool fileNameUsedInAnyCampaign(const std::string& rawRoot, const std::string& fileName) {
  if (exists(rawRoot + "/" + fileName)) return true;
  DIR* d = ::opendir(rawRoot.c_str());
  if (!d) return false;
  bool used = false;
  while (dirent* e = ::readdir(d)) {
    const std::string n = e->d_name;
    if (n == "." || n == "..") continue;
    if (exists(rawRoot + "/" + n + "/" + fileName)) {
      used = true;
      break;
    }
  }
  ::closedir(d);
  return used;
}

Result<RunReservation> RunReservation::reserve(const CampaignLocation& loc, UtcNs startUtc) {
  using R = Result<RunReservation>;
  switch (campaignState(loc)) {
    case CampaignState::Open: break;
    case CampaignState::Closed: return R::failure(-1, "la campaña " + loc.name + " está cerrada: no se escriben runs nuevos");
    case CampaignState::Missing:
      return R::failure(-1, "no hay campaña abierta " + loc.name + " en " + loc.rawRoot + " (daqcampaign open " + loc.name + ")");
    case CampaignState::Invalid: return R::failure(-1, loc.dir() + " no es una campaña");
  }
  RunReservation r;
  r.loc_ = loc;
  auto lk = InstanceLock::acquire(loc.dir() + "/" + kRunLock);
  if (!lk) return R::failure(-1, "ya hay un run en curso en la campaña " + loc.name + ": " + lk.error().message);
  r.lock_ = std::move(lk).value();
  int maxRun = 0;
  if (DIR* d = ::opendir(loc.dir().c_str())) {
    while (dirent* e = ::readdir(d)) maxRun = std::max(maxRun, runNumberOf(loc.name, e->d_name));
    ::closedir(d);
  }
  r.run_ = maxRun + 1;
  char num[16];
  std::snprintf(num, sizeof num, "%03d", r.run_);
  r.base_ = loc.name + "_" + formatUtcCompact(startUtc) + "_run" + num;
  for (const auto& f : {r.segmentFile(1), r.indexFile()})
    if (fileNameUsedInAnyCampaign(loc.rawRoot, f)) return R::failure(-1, "el nombre " + f + " ya existe bajo " + loc.rawRoot);
  return r;
}

std::string RunReservation::segmentFile(int seg) const {
  char num[16];
  std::snprintf(num, sizeof num, "%03d", seg);
  return base_ + "_seg" + num + ".root";
}

Status appendSha256Sum(const std::string& campaignDir, const std::string& fileName, const std::string& sha256) {
  const std::string path = campaignDir + "/" + kSums;
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC, 0644);
  if (fd < 0) return Status::failure(errno, "open " + path + ": " + std::strerror(errno));
  if (::flock(fd, LOCK_EX) != 0) {
    ::close(fd);
    return Status::failure(errno, "flock " + path);
  }
  const std::string line = sha256 + "  " + fileName + "\n";
  const ssize_t n = ::write(fd, line.data(), line.size());
  const int e = errno;
  ::fsync(fd);
  ::close(fd);
  if (n != static_cast<ssize_t>(line.size())) return Status::failure(e, "write " + path);
  return okStatus();
}

Status updateManifestRun(const CampaignLocation& loc, const nlohmann::json& runEntry) {
  FileLock lk(loc.manifestPath() + ".lock");
  if (!lk.ok()) return Status::failure(-1, "no se pudo bloquear " + loc.manifestPath());
  auto man = readJson(loc.manifestPath());
  if (!man) return Status(man.error());
  nlohmann::json m = man.value();
  nlohmann::json runs = m.value("runs", nlohmann::json::array());
  bool replaced = false;
  for (auto& r : runs)
    if (r.value("run", -1) == runEntry.value("run", -2)) {
      r = runEntry;
      replaced = true;
    }
  if (!replaced) runs.push_back(runEntry);
  m["runs"] = runs;
  PosixFileSystem fs;
  return fs.writeAtomic(loc.manifestPath(), m.dump(2) + "\n", 0644);
}

}  // namespace daq::daqrun
