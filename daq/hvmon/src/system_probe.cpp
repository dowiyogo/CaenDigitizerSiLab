// system_probe.cpp
#include "daq/hvmon/system_probe.h"

#include <arpa/inet.h>
#include <dirent.h>
#include <ifaddrs.h>
#include <netinet/in.h>
#include <unistd.h>

#include <cctype>
#include <fstream>
#include <iterator>

namespace daq::hvmon {

bool processRunning(const std::string& pattern, const std::string& procRoot) {
  DIR* d = ::opendir(procRoot.c_str());
  if (!d) return false;
  const std::string self = std::to_string(::getpid());
  bool found = false;
  while (dirent* e = ::readdir(d)) {
    const std::string name = e->d_name;
    if (name.empty() || !std::isdigit(static_cast<unsigned char>(name[0])) || name == self) continue;
    std::ifstream in(procRoot + "/" + name + "/cmdline", std::ios::binary);
    if (!in) continue;
    std::string cmd((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (char& c : cmd)
      if (c == '\0') c = ' ';
    if (cmd.find(pattern) != std::string::npos) {
      found = true;
      break;
    }
  }
  ::closedir(d);
  return found;
}

Result<std::string> resolveHttpBind(const std::string& bind) {
  if (bind == "127.0.0.1" || bind == "localhost") return std::string("127.0.0.1");
  if (bind != "tailscale") return Result<std::string>::failure(-1, "http.bind no permitido: " + bind);
  ifaddrs* ifs = nullptr;
  if (::getifaddrs(&ifs) != 0) return Result<std::string>::failure(-1, "getifaddrs falló");
  std::string out;
  for (ifaddrs* i = ifs; i; i = i->ifa_next) {
    if (!i->ifa_addr || i->ifa_addr->sa_family != AF_INET || std::string(i->ifa_name) != "tailscale0") continue;
    char buf[INET_ADDRSTRLEN];
    const auto* sin = reinterpret_cast<const sockaddr_in*>(i->ifa_addr);
    if (::inet_ntop(AF_INET, &sin->sin_addr, buf, sizeof buf)) out = buf;
  }
  ::freeifaddrs(ifs);
  if (out.empty()) return Result<std::string>::failure(-1, "tailscale0 sin dirección IPv4");
  return out;
}

}  // namespace daq::hvmon
