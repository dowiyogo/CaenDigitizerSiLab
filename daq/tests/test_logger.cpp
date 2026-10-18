// Tests del Logger: una línea completa por llamada, también con varios hilos.
#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "daq/core/logger.h"

TEST(Logger, ConcurrentLinesAreNeverInterleaved) {
  const std::string path = ::testing::TempDir() + "daq_logger_test.log";
  std::remove(path.c_str());
  auto& log = daq::Logger::instance();
  daq::ManualClock clock;
  log.setClock(&clock);
  log.setStderr(false);
  ASSERT_TRUE(log.openFile(path).ok());
  constexpr int kThreads = 4, kLines = 2000;
  std::vector<std::thread> ts;
  for (int t = 0; t < kThreads; ++t)
    ts.emplace_back([t] {
      daq::Logger::setThreadName("w" + std::to_string(t));
      for (int i = 0; i < kLines; ++i) daq::logInfo("linea ", i, " del hilo ", t);
    });
  for (auto& th : ts) th.join();
  log.closeFile();
  log.setStderr(true);
  log.setClock(nullptr);

  std::ifstream in(path);
  std::string line;
  int n = 0;
  // Cada línea debe ser exactamente "<marca> INFO [wT] linea I del hilo T" para algún T, I.
  const std::string prefix = "2026-10-18T00:00:00.000000000Z INFO [w";
  while (std::getline(in, line)) {
    bool ok = line.compare(0, prefix.size(), prefix) == 0 && line.size() > prefix.size() + 2;
    if (ok) {
      const char t = line[prefix.size()];
      const std::string rest = line.substr(prefix.size() + 1);
      const std::string head = "] linea ";
      const std::string tail = std::string(" del hilo ") + t;
      ok = t >= '0' && t <= '3' && rest.compare(0, head.size(), head) == 0 && rest.size() > head.size() + tail.size() &&
           rest.compare(rest.size() - tail.size(), tail.size(), tail) == 0;
      if (ok) {
        const std::string num = rest.substr(head.size(), rest.size() - head.size() - tail.size());
        ok = !num.empty() && num.find_first_not_of("0123456789") == std::string::npos;
      }
    }
    EXPECT_TRUE(ok) << line;
    ++n;
  }
  EXPECT_EQ(n, kThreads * kLines);
  std::remove(path.c_str());
}

TEST(Logger, MinLevelFilters) {
  const std::string path = ::testing::TempDir() + "daq_logger_level.log";
  std::remove(path.c_str());
  auto& log = daq::Logger::instance();
  log.setStderr(false);
  ASSERT_TRUE(log.openFile(path).ok());
  log.setMinLevel(daq::LogLevel::Warn);
  daq::logInfo("no");
  daq::logWarn("si");
  log.setMinLevel(daq::LogLevel::Info);
  log.closeFile();
  log.setStderr(true);
  std::ifstream in(path);
  std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_EQ(all.find(" INFO "), std::string::npos);
  EXPECT_NE(all.find(" WARN "), std::string::npos);
  std::remove(path.c_str());
}
