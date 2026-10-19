// file_system.cpp
#include "daq/core/file_system.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace daq {

namespace {

Error sysError(const std::string& what, const std::string& path) {
  const int e = errno;
  return Error{e, what + " " + path + ": " + std::strerror(e)};
}

class PosixFile final : public WritableFile {
 public:
  PosixFile(int fd, std::string path) : fd_(fd), path_(std::move(path)) {}
  ~PosixFile() override {
    if (fd_ >= 0) ::close(fd_);
  }
  Status append(std::string_view data) override {
    const char* p = data.data();
    std::size_t left = data.size();
    while (left > 0) {
      const ssize_t n = ::write(fd_, p, left);
      if (n < 0) {
        if (errno == EINTR) continue;
        return Status(sysError("write", path_));
      }
      p += n;
      left -= static_cast<std::size_t>(n);
    }
    return okStatus();
  }
  Status sync() override {
    if (::fsync(fd_) != 0) return Status(sysError("fsync", path_));
    return okStatus();
  }
  std::uint64_t size() const override {
    struct stat st {};
    if (::fstat(fd_, &st) != 0) return 0;
    return static_cast<std::uint64_t>(st.st_size);
  }

 private:
  int fd_;
  std::string path_;
};

}  // namespace

namespace {
std::string lexicalNormalize(const std::string& path) {
  std::vector<std::string> parts;
  std::size_t a = 0;
  while (a <= path.size()) {
    const auto b = path.find('/', a);
    const std::string s = path.substr(a, b == std::string::npos ? std::string::npos : b - a);
    if (s == "..") {
      if (!parts.empty()) parts.pop_back();
    } else if (!s.empty() && s != ".") {
      parts.push_back(s);
    }
    if (b == std::string::npos) break;
    a = b + 1;
  }
  std::string out;
  for (const auto& s : parts) out += "/" + s;
  return out.empty() ? "/" : out;
}
}  // namespace

std::string resolvePath(const std::string& path) {
  std::string abs = path;
  if (abs.empty() || abs[0] != '/') {
    char cwd[4096];
    abs = std::string(::getcwd(cwd, sizeof cwd) ? cwd : "/") + "/" + abs;
  }
  abs = lexicalNormalize(abs);
  // Prefijo existente más largo con realpath; el resto se agrega tal cual.
  std::string head = abs, tail;
  for (;;) {
    char buf[PATH_MAX];
    if (::realpath(head.c_str(), buf)) return lexicalNormalize(std::string(buf) + tail);
    if (head == "/" || head.empty()) return abs;
    const auto pos = head.find_last_of('/');
    tail = head.substr(pos) + tail;
    head = pos == 0 ? "/" : head.substr(0, pos);
  }
}

bool pathIsWithin(const std::string& path, const std::string& dir) {
  const std::string p = resolvePath(path), d = resolvePath(dir);
  if (d == "/") return true;
  return p == d || (p.size() > d.size() && p.compare(0, d.size(), d) == 0 && p[d.size()] == '/');
}

std::string parentDirectory(const std::string& path) {
  const auto pos = path.find_last_of('/');
  if (pos == std::string::npos) return "";
  if (pos == 0) return "/";
  return path.substr(0, pos);
}

Result<std::unique_ptr<WritableFile>> PosixFileSystem::openAppend(const std::string& path) {
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC, 0644);
  if (fd < 0) return Result<std::unique_ptr<WritableFile>>(sysError("open", path));
  return std::unique_ptr<WritableFile>(std::make_unique<PosixFile>(fd, path));
}

Status PosixFileSystem::makeDirectories(const std::string& path, unsigned mode) {
  if (path.empty()) return okStatus();
  std::string cur;
  std::size_t pos = 0;
  while (pos != std::string::npos) {
    pos = path.find('/', pos + 1);
    cur = path.substr(0, pos);
    if (cur.empty()) continue;
    if (::mkdir(cur.c_str(), static_cast<mode_t>(mode)) != 0 && errno != EEXIST)
      return Status(sysError("mkdir", cur));
  }
  struct stat st {};
  if (::stat(path.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) return Status::failure(-1, "no es un directorio: " + path);
  return okStatus();
}

Status PosixFileSystem::writeAtomic(const std::string& path, std::string_view content, unsigned mode) {
  const std::string tmp = path + ".tmp." + std::to_string(::getpid());
  const int fd = ::open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, static_cast<mode_t>(mode));
  if (fd < 0) return Status(sysError("open", tmp));
  PosixFile f(fd, tmp);
  if (auto st = f.append(content); !st) {
    ::unlink(tmp.c_str());
    return st;
  }
  if (::fchmod(fd, static_cast<mode_t>(mode)) != 0) {
    ::unlink(tmp.c_str());
    return Status(sysError("fchmod", tmp));
  }
  if (::rename(tmp.c_str(), path.c_str()) != 0) {
    const Error e = sysError("rename", tmp);
    ::unlink(tmp.c_str());
    return Status(e);
  }
  return okStatus();
}

Result<std::string> PosixFileSystem::readFile(const std::string& path) {
  const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
  if (fd < 0) return Result<std::string>(sysError("open", path));
  std::string out;
  char buf[65536];
  for (;;) {
    const ssize_t n = ::read(fd, buf, sizeof buf);
    if (n < 0) {
      if (errno == EINTR) continue;
      const Error e = sysError("read", path);
      ::close(fd);
      return Result<std::string>(e);
    }
    if (n == 0) break;
    out.append(buf, static_cast<std::size_t>(n));
  }
  ::close(fd);
  return out;
}

Result<std::string> PosixFileSystem::readHead(const std::string& path, std::size_t maxBytes) {
  const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
  if (fd < 0) return Result<std::string>(sysError("open", path));
  std::string out(maxBytes, '\0');
  std::size_t got = 0;
  while (got < maxBytes) {
    const ssize_t n = ::read(fd, &out[got], maxBytes - got);
    if (n < 0) {
      if (errno == EINTR) continue;
      const Error e = sysError("read", path);
      ::close(fd);
      return Result<std::string>(e);
    }
    if (n == 0) break;
    got += static_cast<std::size_t>(n);
  }
  ::close(fd);
  out.resize(got);
  return out;
}

bool PosixFileSystem::exists(const std::string& path) {
  struct stat st {};
  return ::stat(path.c_str(), &st) == 0;
}

// ---- MemoryFileSystem -------------------------------------------------------------------------

class MemoryFile final : public WritableFile {
 public:
  MemoryFile(MemoryFileSystem& fs, std::string path) : fs_(fs), path_(std::move(path)) {}
  Status append(std::string_view data) override {
    std::lock_guard<std::mutex> lk(fs_.m_);
    if (fs_.failing_.count(path_)) return Status::failure(28, "write " + path_ + ": No space left on device (simulado)");
    fs_.data_[path_].append(data.data(), data.size());
    return okStatus();
  }
  Status sync() override {
    std::lock_guard<std::mutex> lk(fs_.m_);
    ++fs_.syncs_[path_];
    return okStatus();
  }
  std::uint64_t size() const override {
    std::lock_guard<std::mutex> lk(fs_.m_);
    auto it = fs_.data_.find(path_);
    return it == fs_.data_.end() ? 0 : it->second.size();
  }

 private:
  MemoryFileSystem& fs_;
  std::string path_;
};

Result<std::unique_ptr<WritableFile>> MemoryFileSystem::openAppend(const std::string& path) {
  {
    std::lock_guard<std::mutex> lk(m_);
    data_.emplace(path, std::string());
  }
  return std::unique_ptr<WritableFile>(std::make_unique<MemoryFile>(*this, path));
}

Status MemoryFileSystem::makeDirectories(const std::string&, unsigned) { return okStatus(); }

Status MemoryFileSystem::writeAtomic(const std::string& path, std::string_view content, unsigned) {
  std::lock_guard<std::mutex> lk(m_);
  data_[path] = std::string(content);
  return okStatus();
}

bool MemoryFileSystem::exists(const std::string& path) {
  std::lock_guard<std::mutex> lk(m_);
  return data_.count(path) > 0;
}

Result<std::string> MemoryFileSystem::readFile(const std::string& path) {
  std::lock_guard<std::mutex> lk(m_);
  auto it = data_.find(path);
  if (it == data_.end()) return Result<std::string>::failure(2, "no existe " + path);
  return it->second;
}

Result<std::string> MemoryFileSystem::readHead(const std::string& path, std::size_t maxBytes) {
  std::lock_guard<std::mutex> lk(m_);
  auto it = data_.find(path);
  if (it == data_.end()) return Result<std::string>::failure(2, "no existe " + path);
  return it->second.substr(0, maxBytes);
}

std::string MemoryFileSystem::content(const std::string& path) const {
  std::lock_guard<std::mutex> lk(m_);
  auto it = data_.find(path);
  return it == data_.end() ? std::string() : it->second;
}

std::uint64_t MemoryFileSystem::syncCount(const std::string& path) const {
  std::lock_guard<std::mutex> lk(m_);
  auto it = syncs_.find(path);
  return it == syncs_.end() ? 0 : it->second;
}

std::uint64_t MemoryFileSystem::totalSyncs() const {
  std::lock_guard<std::mutex> lk(m_);
  std::uint64_t n = 0;
  for (const auto& kv : syncs_) n += kv.second;
  return n;
}

std::set<std::string> MemoryFileSystem::files() const {
  std::lock_guard<std::mutex> lk(m_);
  std::set<std::string> out;
  for (const auto& kv : data_) out.insert(kv.first);
  return out;
}

void MemoryFileSystem::failAppendsTo(const std::string& path, bool fail) {
  std::lock_guard<std::mutex> lk(m_);
  if (fail)
    failing_.insert(path);
  else
    failing_.erase(path);
}

}  // namespace daq
