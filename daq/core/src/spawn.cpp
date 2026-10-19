// spawn.cpp
#include "daq/core/spawn.h"

#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>

extern char** environ;

namespace daq {

Result<SpawnResult> PosixSpawner::run(const std::vector<std::string>& argv, const std::string& stdinData,
                                      std::chrono::seconds timeout) {
  using R = Result<SpawnResult>;
  if (argv.empty()) return R::failure(-1, "spawn: argv vacío");
  int inPipe[2], errPipe[2];
  if (::pipe2(inPipe, O_CLOEXEC) != 0) return R::failure(errno, "pipe");
  if (::pipe2(errPipe, O_CLOEXEC) != 0) {
    ::close(inPipe[0]);
    ::close(inPipe[1]);
    return R::failure(errno, "pipe");
  }
  posix_spawn_file_actions_t fa;
  posix_spawn_file_actions_init(&fa);
  posix_spawn_file_actions_adddup2(&fa, inPipe[0], 0);
  posix_spawn_file_actions_addopen(&fa, 1, "/dev/null", O_WRONLY, 0);
  posix_spawn_file_actions_adddup2(&fa, errPipe[1], 2);
  std::vector<char*> args;
  for (const auto& a : argv) args.push_back(const_cast<char*>(a.c_str()));
  args.push_back(nullptr);
  pid_t pid = -1;
  const int rc = posix_spawnp(&pid, argv[0].c_str(), &fa, nullptr, args.data(), environ);
  posix_spawn_file_actions_destroy(&fa);
  ::close(inPipe[0]);
  ::close(errPipe[1]);
  if (rc != 0) {
    ::close(inPipe[1]);
    ::close(errPipe[0]);
    return R::failure(rc, std::string("posix_spawnp ") + argv[0] + ": " + std::strerror(rc));
  }
  // stdin: escribir todo y cerrar (los datos son pocos KiB; el hijo los lee enseguida).
  {
    const char* p = stdinData.data();
    std::size_t left = stdinData.size();
    signal(SIGPIPE, SIG_IGN);
    while (left > 0) {
      const ssize_t n = ::write(inPipe[1], p, left);
      if (n < 0) {
        if (errno == EINTR) continue;
        break;  // el hijo cerró stdin; se informa por su código de salida
      }
      p += n;
      left -= static_cast<std::size_t>(n);
    }
    ::close(inPipe[1]);
  }
  SpawnResult res;
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  bool termSent = false;
  int status = 0;
  char buf[4096];
  bool errOpen = true;
  for (;;) {
    const pid_t w = ::waitpid(pid, &status, WNOHANG);
    if (w == pid) break;
    if (std::chrono::steady_clock::now() >= deadline) {
      res.timedOut = true;
      ::kill(pid, termSent ? SIGKILL : SIGTERM);
      termSent = true;
    }
    if (errOpen) {
      pollfd pf{errPipe[0], POLLIN, 0};
      if (::poll(&pf, 1, 100) > 0) {
        const ssize_t n = ::read(errPipe[0], buf, sizeof buf);
        if (n > 0 && res.stderrText.size() < 65536) res.stderrText.append(buf, static_cast<std::size_t>(n));
        if (n == 0) errOpen = false;
      }
    } else {
      ::usleep(50000);
    }
  }
  while (errOpen) {
    const ssize_t n = ::read(errPipe[0], buf, sizeof buf);
    if (n <= 0) break;
    if (res.stderrText.size() < 65536) res.stderrText.append(buf, static_cast<std::size_t>(n));
  }
  ::close(errPipe[0]);
  res.exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
  return res;
}

}  // namespace daq
