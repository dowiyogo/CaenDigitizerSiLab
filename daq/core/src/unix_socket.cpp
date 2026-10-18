// unix_socket.cpp
#include "daq/core/unix_socket.h"

#include <poll.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/un.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>

namespace daq {

namespace {

Error sysError(const std::string& what) {
  const int e = errno;
  return Error{e, what + ": " + std::strerror(e)};
}

Result<sockaddr_un> makeAddress(const std::string& path) {
  sockaddr_un addr{};
  addr.sun_family = AF_UNIX;
  if (path.size() >= sizeof(addr.sun_path)) return Result<sockaddr_un>::failure(-1, "ruta de socket demasiado larga: " + path);
  std::memcpy(addr.sun_path, path.c_str(), path.size() + 1);
  return addr;
}

int remainingMs(std::chrono::steady_clock::time_point deadline) {
  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::steady_clock::now()).count();
  return ms < 0 ? 0 : static_cast<int>(ms);
}

}  // namespace

FdGuard::~FdGuard() {
  if (fd_ >= 0) ::close(fd_);
}

FdGuard& FdGuard::operator=(FdGuard&& o) noexcept {
  if (this != &o) {
    if (fd_ >= 0) ::close(fd_);
    fd_ = o.fd_;
    o.fd_ = -1;
  }
  return *this;
}

Result<UnixListener> UnixListener::listen(const std::string& path) {
  auto addr = makeAddress(path);
  if (!addr) return Result<UnixListener>(addr.error());
  ::unlink(path.c_str());
  const int fd = ::socket(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0);
  if (fd < 0) return Result<UnixListener>(sysError("socket"));
  UnixListener l;
  l.fd_ = fd;
  l.path_ = path;
  if (::bind(fd, reinterpret_cast<const sockaddr*>(&addr.value()), sizeof(sockaddr_un)) != 0)
    return Result<UnixListener>(sysError("bind " + path));
  if (::chmod(path.c_str(), 0600) != 0) return Result<UnixListener>(sysError("chmod " + path));
  if (::listen(fd, 16) != 0) return Result<UnixListener>(sysError("listen " + path));
  return l;
}

UnixListener::~UnixListener() {
  if (fd_ >= 0) {
    ::close(fd_);
    ::unlink(path_.c_str());
  }
}

UnixListener::UnixListener(UnixListener&& o) noexcept : fd_(o.fd_), path_(std::move(o.path_)) { o.fd_ = -1; }

UnixListener& UnixListener::operator=(UnixListener&& o) noexcept {
  if (this != &o) {
    if (fd_ >= 0) {
      ::close(fd_);
      ::unlink(path_.c_str());
    }
    fd_ = o.fd_;
    path_ = std::move(o.path_);
    o.fd_ = -1;
  }
  return *this;
}

Result<FdGuard> UnixListener::acceptFor(std::chrono::milliseconds timeout) {
  pollfd p{fd_, POLLIN, 0};
  const int r = ::poll(&p, 1, static_cast<int>(timeout.count()));
  if (r < 0) {
    if (errno == EINTR) return FdGuard(-1);
    return Result<FdGuard>(sysError("poll"));
  }
  if (r == 0) return FdGuard(-1);
  const int c = ::accept4(fd_, nullptr, nullptr, SOCK_CLOEXEC);
  if (c < 0) {
    if (errno == EINTR || errno == EAGAIN || errno == ECONNABORTED) return FdGuard(-1);
    return Result<FdGuard>(sysError("accept"));
  }
  return FdGuard(c);
}

PeerCredentials peerCredentials(int fd) {
  ucred cred{};
  socklen_t len = sizeof cred;
  PeerCredentials out;
  if (::getsockopt(fd, SOL_SOCKET, SO_PEERCRED, &cred, &len) == 0) {
    out.pid = cred.pid;
    out.uid = cred.uid;
  }
  return out;
}

Result<std::string> readLine(int fd, std::size_t maxBytes, std::chrono::milliseconds timeout) {
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  std::string out;
  char buf[4096];
  for (;;) {
    pollfd p{fd, POLLIN, 0};
    const int r = ::poll(&p, 1, remainingMs(deadline));
    if (r < 0) {
      if (errno == EINTR) continue;
      return Result<std::string>(sysError("poll"));
    }
    if (r == 0) return Result<std::string>::failure(-1, "timeout leyendo la petición");
    const ssize_t n = ::recv(fd, buf, sizeof buf, 0);
    if (n < 0) {
      if (errno == EINTR) continue;
      return Result<std::string>(sysError("recv"));
    }
    if (n == 0) return Result<std::string>::failure(-1, "el otro extremo cerró antes del fin de línea");
    out.append(buf, static_cast<std::size_t>(n));
    const auto nl = out.find('\n');
    if ((nl != std::string::npos && nl > maxBytes) || (nl == std::string::npos && out.size() > maxBytes))
      return Result<std::string>::failure(-1, "línea demasiado larga");
    if (nl != std::string::npos) {
      out.resize(nl);
      return out;
    }
  }
}

Status writeAll(int fd, std::string_view data) {
  const char* p = data.data();
  std::size_t left = data.size();
  while (left > 0) {
    const ssize_t n = ::send(fd, p, left, MSG_NOSIGNAL);
    if (n < 0) {
      if (errno == EINTR) continue;
      return Status(sysError("send"));
    }
    p += n;
    left -= static_cast<std::size_t>(n);
  }
  return okStatus();
}

Result<std::string> unixRequest(const std::string& path, const std::string& line, std::chrono::milliseconds timeout) {
  auto addr = makeAddress(path);
  if (!addr) return Result<std::string>(addr.error());
  FdGuard fd(::socket(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0));
  if (fd.get() < 0) return Result<std::string>(sysError("socket"));
  if (::connect(fd.get(), reinterpret_cast<const sockaddr*>(&addr.value()), sizeof(sockaddr_un)) != 0)
    return Result<std::string>(sysError("connect " + path));
  if (auto st = writeAll(fd.get(), line + "\n"); !st) return Result<std::string>(st.error());
  return readLine(fd.get(), 16u << 20, timeout);
}

}  // namespace daq
