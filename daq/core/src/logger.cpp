// logger.cpp
#include "daq/core/logger.h"

#include <functional>
#include <iostream>
#include <thread>

namespace daq {

namespace {
thread_local std::string tlsThreadName;
}

const char* toString(LogLevel level) {
  switch (level) {
    case LogLevel::Debug: return "DEBUG";
    case LogLevel::Info: return "INFO";
    case LogLevel::Warn: return "WARN";
    case LogLevel::Error: return "ERROR";
  }
  return "?";
}

Logger& Logger::instance() {
  static Logger logger;
  return logger;
}

void Logger::setMinLevel(LogLevel level) {
  std::lock_guard<std::mutex> lk(m_);
  minLevel_ = level;
}

void Logger::setStderr(bool enabled) {
  std::lock_guard<std::mutex> lk(m_);
  toStderr_ = enabled;
}

Status Logger::openFile(const std::string& path) {
  std::lock_guard<std::mutex> lk(m_);
  if (file_.is_open()) file_.close();
  file_.open(path, std::ios::out | std::ios::app);
  if (!file_) return Status::failure(-1, "no se pudo abrir el log " + path);
  return okStatus();
}

void Logger::closeFile() {
  std::lock_guard<std::mutex> lk(m_);
  if (file_.is_open()) file_.close();
}

void Logger::setClock(const Clock* clock) {
  std::lock_guard<std::mutex> lk(m_);
  clock_ = clock;
}

void Logger::setThreadName(std::string name) { tlsThreadName = std::move(name); }

const std::string& Logger::threadName() {
  if (tlsThreadName.empty()) {
    tlsThreadName = "t" + std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()) % 100000);
  }
  return tlsThreadName;
}

void Logger::log(LogLevel level, std::string_view message) {
  const std::string& tname = threadName();
  std::lock_guard<std::mutex> lk(m_);
  if (level < minLevel_) return;
  const UtcNs t = clock_ ? clock_->now_utc_ns() : utcNowNs();
  std::string line = formatUtcIso(t);
  line += ' ';
  line += toString(level);
  line += " [";
  line += tname;
  line += "] ";
  line += message;
  line += '\n';
  if (toStderr_) {
    std::cerr << line;
    std::cerr.flush();
  }
  if (file_.is_open()) {
    file_ << line;
    file_.flush();
  }
}

}  // namespace daq
