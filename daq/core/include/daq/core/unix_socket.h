// unix_socket.h
// Socket Unix de flujo para el protocolo de comandos de hvmon: una petición JSON por línea y una
// respuesta JSON por línea, una petición por conexión.
// - UnixListener: bind + listen; el socket queda con modo 600. acceptFor() espera con timeout
//   para poder revisar la parada. El destructor cierra y borra el socket.
// - readLine(): lee hasta '\n' con timeout total y tamaño máximo; un cliente que corta a mitad de
//   petición produce un error, no un cuelgue.
// - writeAll(): send con MSG_NOSIGNAL (sin SIGPIPE si el cliente ya cerró).
// - unixRequest(): cliente completo (conecta, envía una línea, lee una línea).
// Seguridad entre hilos: cada descriptor pertenece a un hilo.
#pragma once

#include <chrono>
#include <string>
#include <string_view>

#include "daq/core/result.h"

namespace daq {

class FdGuard {
 public:
  explicit FdGuard(int fd = -1) : fd_(fd) {}
  ~FdGuard();
  FdGuard(FdGuard&& o) noexcept : fd_(o.fd_) { o.fd_ = -1; }
  FdGuard& operator=(FdGuard&& o) noexcept;
  FdGuard(const FdGuard&) = delete;
  FdGuard& operator=(const FdGuard&) = delete;
  int get() const { return fd_; }

 private:
  int fd_;
};

class UnixListener {
 public:
  UnixListener() = default;
  // Borra un socket viejo en path (el llamador debe tener el lock de instancia).
  static Result<UnixListener> listen(const std::string& path);
  ~UnixListener();
  UnixListener(UnixListener&& o) noexcept;
  UnixListener& operator=(UnixListener&& o) noexcept;
  UnixListener(const UnixListener&) = delete;
  UnixListener& operator=(const UnixListener&) = delete;

  // Conexión aceptada, o FdGuard(-1) si venció el timeout.
  Result<FdGuard> acceptFor(std::chrono::milliseconds timeout);
  const std::string& path() const { return path_; }

 private:
  int fd_ = -1;
  std::string path_;
};

struct PeerCredentials {
  int pid = -1;
  unsigned uid = 0;
};
PeerCredentials peerCredentials(int fd);

Result<std::string> readLine(int fd, std::size_t maxBytes, std::chrono::milliseconds timeout);
Status writeAll(int fd, std::string_view data);
Result<std::string> unixRequest(const std::string& path, const std::string& line, std::chrono::milliseconds timeout);

}  // namespace daq
