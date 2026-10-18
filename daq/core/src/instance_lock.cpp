// instance_lock.cpp
#include "daq/core/instance_lock.h"

#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>

namespace daq {

Result<InstanceLock> InstanceLock::acquire(const std::string& path) {
  const int fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0600);
  if (fd < 0) return Result<InstanceLock>::failure(errno, "open " + path + ": " + std::strerror(errno));
  if (::flock(fd, LOCK_EX | LOCK_NB) != 0) {
    const int e = errno;
    ::close(fd);
    if (e == EWOULDBLOCK) return Result<InstanceLock>::failure(e, "lock tomado por otro proceso: " + path);
    return Result<InstanceLock>::failure(e, "flock " + path + ": " + std::strerror(e));
  }
  return InstanceLock(fd);
}

bool InstanceLock::isHeld(const std::string& path) {
  const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
  if (fd < 0) return false;  // sin archivo no hay dueño
  const bool held = ::flock(fd, LOCK_EX | LOCK_NB) != 0 && errno == EWOULDBLOCK;
  ::close(fd);  // si lo tomamos, cerrar lo libera
  return held;
}

InstanceLock::~InstanceLock() { release(); }

InstanceLock& InstanceLock::operator=(InstanceLock&& o) noexcept {
  if (this != &o) {
    release();
    fd_ = o.fd_;
    o.fd_ = -1;
  }
  return *this;
}

void InstanceLock::release() {
  if (fd_ >= 0) {
    ::flock(fd_, LOCK_UN);
    ::close(fd_);
    fd_ = -1;
  }
}

Status ensurePrivateDirectory(const std::string& path, unsigned mode) {
  if (::mkdir(path.c_str(), static_cast<mode_t>(mode)) != 0 && errno != EEXIST)
    return Status::failure(errno, "mkdir " + path + ": " + std::strerror(errno));
  struct stat st {};
  if (::stat(path.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) return Status::failure(-1, "no es un directorio: " + path);
  if (st.st_uid != ::getuid()) return Status::failure(-1, "el directorio no pertenece al usuario actual: " + path);
  if ((static_cast<unsigned>(st.st_mode) & 0777u) != mode) ::chmod(path.c_str(), static_cast<mode_t>(mode));
  return okStatus();
}

}  // namespace daq
