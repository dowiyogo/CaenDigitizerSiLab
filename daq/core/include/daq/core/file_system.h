// file_system.h
// Capa de archivos inyectable. El código que escribe archivos de control lento la usa para que
// los tests puedan contar fsync y simular errores (MemoryFileSystem) sin tocar el disco.
// - WritableFile::append() hace una sola llamada write() con O_APPEND (sin búfer de usuario), así
//   que cada línea llega al kernel de inmediato; sync() hace fsync().
// - writeAtomic(): escribe en un temporal del mismo directorio y hace rename(), de modo que un
//   lector nunca ve un archivo a medio escribir.
// Seguridad entre hilos: PosixFileSystem no tiene estado y se puede usar desde varios hilos; cada
// WritableFile pertenece a un solo hilo. MemoryFileSystem protege su estado con un mutex.
#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <string_view>

#include "daq/core/result.h"

namespace daq {

class WritableFile {
 public:
  virtual ~WritableFile() = default;
  virtual Status append(std::string_view data) = 0;
  virtual Status sync() = 0;
  virtual std::uint64_t size() const = 0;  // tamaño actual del archivo en bytes
};

class FileSystem {
 public:
  virtual ~FileSystem() = default;
  // Abre (o crea con modo 0644) para agregar al final.
  virtual Result<std::unique_ptr<WritableFile>> openAppend(const std::string& path) = 0;
  virtual Status makeDirectories(const std::string& path, unsigned mode) = 0;
  virtual Status writeAtomic(const std::string& path, std::string_view content, unsigned mode) = 0;
  virtual bool exists(const std::string& path) = 0;
  virtual Result<std::string> readFile(const std::string& path) = 0;
  // Primeros maxBytes bytes (menos si el archivo es más corto); p. ej. para leer un encabezado.
  virtual Result<std::string> readHead(const std::string& path, std::size_t maxBytes) = 0;
};

class PosixFileSystem final : public FileSystem {
 public:
  Result<std::unique_ptr<WritableFile>> openAppend(const std::string& path) override;
  Status makeDirectories(const std::string& path, unsigned mode) override;
  Status writeAtomic(const std::string& path, std::string_view content, unsigned mode) override;
  bool exists(const std::string& path) override;
  Result<std::string> readFile(const std::string& path) override;
  Result<std::string> readHead(const std::string& path, std::size_t maxBytes) override;
};

// Sistema de archivos en memoria para tests: cuenta fsync por archivo y permite inyectar fallos.
class MemoryFileSystem final : public FileSystem {
 public:
  Result<std::unique_ptr<WritableFile>> openAppend(const std::string& path) override;
  Status makeDirectories(const std::string& path, unsigned mode) override;
  Status writeAtomic(const std::string& path, std::string_view content, unsigned mode) override;
  bool exists(const std::string& path) override;
  Result<std::string> readFile(const std::string& path) override;
  Result<std::string> readHead(const std::string& path, std::size_t maxBytes) override;

  std::string content(const std::string& path) const;
  std::uint64_t syncCount(const std::string& path) const;
  std::uint64_t totalSyncs() const;
  std::set<std::string> files() const;
  void failAppendsTo(const std::string& path, bool fail);  // las próximas append fallan

 private:
  friend class MemoryFile;
  mutable std::mutex m_;
  std::map<std::string, std::string> data_;
  std::map<std::string, std::uint64_t> syncs_;
  std::set<std::string> failing_;
};

// Directorio que contiene path ("" si no hay '/').
std::string parentDirectory(const std::string& path);

// Ruta absoluta sin "." ni "..", con los enlaces simbólicos de la parte que ya existe resueltos
// (realpath del prefijo existente más largo; el resto se normaliza léxicamente).
std::string resolvePath(const std::string& path);
// true si path es dir o está debajo de dir (ambas resueltas con resolvePath).
bool pathIsWithin(const std::string& path, const std::string& dir);

}  // namespace daq
