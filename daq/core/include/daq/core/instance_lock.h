// instance_lock.h
// Lock de instancia única con flock(LOCK_EX | LOCK_NB) sobre un archivo. Lo usa hvmon (dueño de
// la sesión de la Wrapper) y hvctl en modo directo. El lock lo libera el kernel si el proceso
// muere, así que no quedan locks huérfanos. flock es por descripción de archivo abierta: dos
// aperturas en el mismo proceso también se excluyen (útil en tests).
// Seguridad entre hilos: un InstanceLock pertenece a un hilo; isHeld() es segura.
#pragma once

#include <string>

#include "daq/core/result.h"

namespace daq {

class InstanceLock {
 public:
  InstanceLock() = default;
  // Crea el archivo (modo 600) si no existe y toma el lock; falla si otro lo tiene.
  static Result<InstanceLock> acquire(const std::string& path);
  // true si otro proceso (u otra descripción) tiene el lock ahora.
  static bool isHeld(const std::string& path);

  ~InstanceLock();
  InstanceLock(InstanceLock&& o) noexcept : fd_(o.fd_) { o.fd_ = -1; }
  InstanceLock& operator=(InstanceLock&& o) noexcept;
  InstanceLock(const InstanceLock&) = delete;
  InstanceLock& operator=(const InstanceLock&) = delete;

  bool held() const { return fd_ >= 0; }
  void release();

 private:
  explicit InstanceLock(int fd) : fd_(fd) {}
  int fd_ = -1;
};

// Crea el directorio con el modo dado si no existe y comprueba que sea del usuario actual.
Status ensurePrivateDirectory(const std::string& path, unsigned mode);

}  // namespace daq
