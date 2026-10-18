// Tests de infraestructura: lock de instancia, socket Unix, escritura atómica y archivos en memoria.
#include <gtest/gtest.h>

#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/un.h>
#include <unistd.h>

#include <atomic>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <thread>

#include <nlohmann/json.hpp>

#include "daq/core/file_system.h"
#include "daq/core/instance_lock.h"
#include "daq/core/jthread.h"
#include "daq/core/unix_socket.h"

using namespace std::chrono_literals;

TEST(InstanceLock, SecondAcquireFailsUntilRelease) {
  const std::string path = ::testing::TempDir() + "daq_test.lock";
  std::remove(path.c_str());
  EXPECT_FALSE(daq::InstanceLock::isHeld(path));
  auto a = daq::InstanceLock::acquire(path);
  ASSERT_TRUE(a.ok());
  EXPECT_TRUE(daq::InstanceLock::isHeld(path));
  auto b = daq::InstanceLock::acquire(path);
  ASSERT_FALSE(b.ok());
  EXPECT_NE(b.error().message.find("lock tomado"), std::string::npos);
  a.value().release();
  EXPECT_FALSE(daq::InstanceLock::isHeld(path));
  auto c = daq::InstanceLock::acquire(path);
  EXPECT_TRUE(c.ok());
  struct stat st {};
  ASSERT_EQ(::stat(path.c_str(), &st), 0);
  EXPECT_EQ(st.st_mode & 0777, 0600u);
}

TEST(InstanceLock, PrivateDirectoryGetsMode700) {
  const std::string dir = ::testing::TempDir() + "daq_test_private_dir";
  ::rmdir(dir.c_str());
  ASSERT_TRUE(daq::ensurePrivateDirectory(dir, 0700).ok());
  struct stat st {};
  ASSERT_EQ(::stat(dir.c_str(), &st), 0);
  EXPECT_EQ(st.st_mode & 0777, 0700u);
  ::rmdir(dir.c_str());
}

namespace {
// Servidor de eco con una petición por conexión, en su propio hilo.
struct EchoServer {
  std::string path;
  daq::UnixListener listener;
  std::atomic<int> served{0}, failed{0};
  daq::JThread thread;
  explicit EchoServer(const std::string& p) : path(p) {
    auto l = daq::UnixListener::listen(p);
    EXPECT_TRUE(l.ok());
    listener = std::move(l).value();
    thread = daq::JThread([this](daq::StopToken st) {
      while (!st.stop_requested()) {
        auto c = listener.acceptFor(50ms);
        if (!c || c.value().get() < 0) continue;
        auto line = daq::readLine(c.value().get(), 1024, 500ms);
        if (!line) {
          ++failed;
          continue;
        }
        daq::writeAll(c.value().get(), "eco:" + line.value() + "\n");
        ++served;
      }
    });
  }
};

int rawConnect(const std::string& path) {
  const int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
  sockaddr_un a{};
  a.sun_family = AF_UNIX;
  std::strncpy(a.sun_path, path.c_str(), sizeof(a.sun_path) - 1);
  EXPECT_EQ(::connect(fd, reinterpret_cast<sockaddr*>(&a), sizeof a), 0);
  return fd;
}
}  // namespace

TEST(UnixSocket, RequestResponseAndMode600) {
  EchoServer s(::testing::TempDir() + "daq_test.sock");
  struct stat st {};
  ASSERT_EQ(::stat(s.path.c_str(), &st), 0);
  EXPECT_EQ(st.st_mode & 0777, 0600u);
  auto r = daq::unixRequest(s.path, "hola", 2s);
  ASSERT_TRUE(r.ok()) << r.error().message;
  EXPECT_EQ(r.value(), "eco:hola");
}

TEST(UnixSocket, ClientDisconnectingMidRequestDoesNotStopServer) {
  EchoServer s(::testing::TempDir() + "daq_test2.sock");
  const int fd = rawConnect(s.path);
  ASSERT_GE(fd, 0);
  ::send(fd, "{\"cmd\":\"sta", 11, MSG_NOSIGNAL);  // sin fin de línea
  ::close(fd);
  for (int i = 0; i < 100 && s.failed.load() == 0; ++i) std::this_thread::sleep_for(10ms);
  EXPECT_EQ(s.failed.load(), 1);
  auto r = daq::unixRequest(s.path, "sigue", 2s);
  ASSERT_TRUE(r.ok());
  EXPECT_EQ(r.value(), "eco:sigue");
}

TEST(UnixSocket, OverlongLineIsRejected) {
  EchoServer s(::testing::TempDir() + "daq_test3.sock");
  auto r = daq::unixRequest(s.path, std::string(5000, 'x'), 2s);
  EXPECT_FALSE(r.ok());  // el servidor no responde a una línea de más de 1024 bytes
  for (int i = 0; i < 100 && s.failed.load() == 0; ++i) std::this_thread::sleep_for(10ms);
  EXPECT_EQ(s.failed.load(), 1);
}

TEST(UnixSocket, ConnectToMissingSocketFails) {
  EXPECT_FALSE(daq::unixRequest(::testing::TempDir() + "no_existe.sock", "x", 1s).ok());
}

TEST(StateFile, ConcurrentReaderNeverSeesPartialFile) {
  const std::string path = ::testing::TempDir() + "daq_test_state.json";
  daq::PosixFileSystem fs;
  ASSERT_TRUE(fs.writeAtomic(path, nlohmann::json{{"seq", 0}}.dump(), 0600).ok());
  std::atomic<bool> done{false};
  std::atomic<long> reads{0}, bad{0};
  std::thread reader([&] {
    while (!done.load()) {
      std::ifstream in(path);
      const std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      try {
        auto j = nlohmann::json::parse(s);
        if (!j.contains("seq")) ++bad;
      } catch (...) {
        ++bad;
      }
      ++reads;
    }
  });
  for (int i = 1; i <= 3000; ++i) {
    // Tamaño variable: un archivo truncado o mezclado no sería JSON válido.
    nlohmann::json j{{"seq", i}, {"pad", std::string(static_cast<std::size_t>((i % 7) * 3000), 'a')}};
    ASSERT_TRUE(fs.writeAtomic(path, j.dump(), 0600).ok());
  }
  done = true;
  reader.join();
  EXPECT_GT(reads.load(), 0);
  EXPECT_EQ(bad.load(), 0);
  struct stat st {};
  ASSERT_EQ(::stat(path.c_str(), &st), 0);
  EXPECT_EQ(st.st_mode & 0777, 0600u);
  std::remove(path.c_str());
}

TEST(MemoryFileSystem, CountsSyncsAndInjectsFailures) {
  daq::MemoryFileSystem fs;
  auto f = fs.openAppend("/x/a.csv");
  ASSERT_TRUE(f.ok());
  EXPECT_TRUE(f.value()->append("l1\n").ok());
  EXPECT_TRUE(f.value()->sync().ok());
  EXPECT_TRUE(f.value()->sync().ok());
  EXPECT_EQ(fs.syncCount("/x/a.csv"), 2u);
  fs.failAppendsTo("/x/a.csv", true);
  EXPECT_FALSE(f.value()->append("l2\n").ok());
  fs.failAppendsTo("/x/a.csv", false);
  EXPECT_TRUE(f.value()->append("l3\n").ok());
  EXPECT_EQ(fs.content("/x/a.csv"), "l1\nl3\n");
  EXPECT_EQ(f.value()->size(), 6u);
}
