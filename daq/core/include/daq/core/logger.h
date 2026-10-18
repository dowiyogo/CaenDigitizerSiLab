// logger.h
// Registro de líneas seguro entre hilos, a stderr y opcionalmente a un archivo.
// Formato de cada línea: "<UTC ISO-8601 con ns> <NIVEL> [<hilo>] <mensaje>".
// Cada llamada a log() escribe una línea completa bajo un mutex (C++17: sin osyncstream), así
// que las líneas de distintos hilos nunca se mezclan. El nombre del hilo se fija con
// setThreadName() en el propio hilo (por defecto "main" o el id numérico).
#pragma once

#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <string_view>

#include "daq/core/clock.h"
#include "daq/core/result.h"

namespace daq {

enum class LogLevel { Debug = 0, Info = 1, Warn = 2, Error = 3 };
const char* toString(LogLevel level);

class Logger {
 public:
  static Logger& instance();

  void setMinLevel(LogLevel level);
  void setStderr(bool enabled);
  Status openFile(const std::string& path);  // agrega al final; crea si no existe
  void closeFile();
  // Reloj para la marca de tiempo (por defecto SystemClock). El reloj debe vivir más que el uso.
  void setClock(const Clock* clock);

  void log(LogLevel level, std::string_view message);

  static void setThreadName(std::string name);
  static const std::string& threadName();

 private:
  Logger() = default;
  std::mutex m_;
  LogLevel minLevel_ = LogLevel::Info;
  bool toStderr_ = true;
  std::ofstream file_;
  const Clock* clock_ = nullptr;
};

namespace detail {
template <class... Args>
std::string concat(Args&&... args) {
  std::ostringstream os;
  (os << ... << args);
  return os.str();
}
}  // namespace detail

template <class... Args>
void logDebug(Args&&... a) { Logger::instance().log(LogLevel::Debug, detail::concat(a...)); }
template <class... Args>
void logInfo(Args&&... a) { Logger::instance().log(LogLevel::Info, detail::concat(a...)); }
template <class... Args>
void logWarn(Args&&... a) { Logger::instance().log(LogLevel::Warn, detail::concat(a...)); }
template <class... Args>
void logError(Args&&... a) { Logger::instance().log(LogLevel::Error, detail::concat(a...)); }

}  // namespace daq
