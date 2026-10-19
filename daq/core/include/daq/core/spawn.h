// spawn.h
// Ejecución de un programa externo (curl para el correo) con datos por stdin y timeout.
// - Spawner: interfaz inyectable; los tests usan una implementación que registra argv y stdin para
//   comprobar que ningún secreto va en argv.
// - PosixSpawner: posix_spawnp con stdin y stderr por pipes (stdout descartado). Escribe stdin y lo
//   cierra, lee stderr y espera la salida; tras timeout manda SIGTERM y luego SIGKILL.
// Seguridad entre hilos: PosixSpawner no tiene estado; se puede usar desde varios hilos.
#pragma once

#include <chrono>
#include <string>
#include <vector>

#include "daq/core/result.h"

namespace daq {

struct SpawnResult {
  int exitCode = -1;      // código de salida; -1 si terminó por señal
  bool timedOut = false;
  std::string stderrText;  // primeros 64 KiB
};

class Spawner {
 public:
  virtual ~Spawner() = default;
  virtual Result<SpawnResult> run(const std::vector<std::string>& argv, const std::string& stdinData,
                                  std::chrono::seconds timeout) = 0;
};

class PosixSpawner final : public Spawner {
 public:
  Result<SpawnResult> run(const std::vector<std::string>& argv, const std::string& stdinData,
                          std::chrono::seconds timeout) override;
};

}  // namespace daq
