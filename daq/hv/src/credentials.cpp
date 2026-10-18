// credentials.cpp
#include "daq/hv/credentials.h"

#include <sys/stat.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>

#include <nlohmann/json.hpp>

#include "daq/core/logger.h"

namespace daq::hv {

std::string expandHome(const std::string& path) {
  if (path.empty() || path[0] != '~') return path;
  const char* home = std::getenv("HOME");
  return std::string(home ? home : "") + path.substr(1);
}

std::string defaultCredentialsPath() { return expandHome("~/.config/daqv2/hv_credentials.json"); }

Result<Credentials> loadCredentials(const std::string& rawPath) {
  const std::string path = expandHome(rawPath);
  struct stat st {};
  if (::stat(path.c_str(), &st) != 0)
    return Result<Credentials>::failure(-1, "credenciales: no existe " + path);
  if (!S_ISREG(st.st_mode)) return Result<Credentials>::failure(-1, "credenciales: no es un archivo regular: " + path);
  if (st.st_uid != ::getuid())
    return Result<Credentials>::failure(-1, "credenciales: el archivo no pertenece al usuario actual: " + path);
  const unsigned mode = static_cast<unsigned>(st.st_mode) & 0777u;
  if (mode != 0600u) {
    char buf[8];
    std::snprintf(buf, sizeof buf, "%03o", mode);
    return Result<Credentials>::failure(-1, "credenciales: modo " + std::string(buf) + " (se exige 600): " + path);
  }
  std::ifstream in(path);
  nlohmann::json j;
  try {
    in >> j;
  } catch (const std::exception& e) {
    return Result<Credentials>::failure(-1, "credenciales: JSON inválido en " + path + ": " + e.what());
  }
  if (!j.is_object() || !j.contains("user") || !j.contains("password") || !j["user"].is_string() ||
      !j["password"].is_string())
    return Result<Credentials>::failure(-1, "credenciales: faltan 'user' o 'password' en " + path);
  return Credentials{j["user"].get<std::string>(), j["password"].get<std::string>()};
}

Result<HvLogin> resolveHvLogin(Secrets& secrets, const std::string& defaultHost, const std::string& legacyPath) {
  using R = Result<HvLogin>;
  HvLogin login;
  auto host = secrets.get("DAQ_HV_HOST");
  if (!host) return R(host.error());
  login.host = host.value().value_or(defaultHost);
  auto user = secrets.get("DAQ_HV_USER");
  if (!user) return R(user.error());
  auto pass = secrets.get("DAQ_HV_PASSWORD");
  if (!pass) return R(pass.error());
  if (user.value() && pass.value()) {
    login.credentials = Credentials{*user.value(), *pass.value()};
    login.source = "secrets";
    return login;
  }
  auto legacy = loadCredentials(legacyPath);
  if (!legacy)
    return R::failure(-1, "faltan DAQ_HV_USER/DAQ_HV_PASSWORD (entorno o " + secrets.path() +
                              ") y el respaldo falló: " + legacy.error().message);
  logWarn("credenciales HV desde ", expandHome(legacyPath),
          ": mecanismo OBSOLETO, usar DAQ_HV_USER/DAQ_HV_PASSWORD en ", secrets.path());
  login.credentials = legacy.value();
  login.source = "hv_credentials.json (obsoleto)";
  return login;
}

}  // namespace daq::hv
