// secrets.cpp
#include "daq/core/secrets.h"

#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

namespace daq {

namespace {

std::string homeExpand(const std::string& p) {
  if (p.empty() || p[0] != '~') return p;
  const char* home = std::getenv("HOME");
  return std::string(home ? home : "") + p.substr(1);
}

std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t");
  return s.substr(b, e - b + 1);
}

}  // namespace

Status checkPrivateFile(const std::string& path) {
  struct stat st {};
  if (::stat(path.c_str(), &st) != 0)
    return Status::failure(errno, "no existe o no se puede leer: " + path);
  if (!S_ISREG(st.st_mode)) return Status::failure(-1, "no es un archivo regular: " + path);
  if (st.st_uid != ::getuid()) return Status::failure(-1, "el archivo no pertenece al usuario actual: " + path);
  const unsigned mode = static_cast<unsigned>(st.st_mode) & 0777u;
  if (mode != 0600u) {
    char buf[8];
    std::snprintf(buf, sizeof buf, "%03o", mode);
    return Status::failure(-1, "modo " + std::string(buf) + " (se exige 600): " + path);
  }
  return okStatus();
}

Result<std::map<std::string, std::string>> parseEnvFile(const std::string& path) {
  using R = Result<std::map<std::string, std::string>>;
  if (auto st = checkPrivateFile(path); !st) return R(st.error());
  std::ifstream in(path);
  if (!in) return R::failure(-1, "no se pudo abrir " + path);
  std::map<std::string, std::string> out;
  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    ++n;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::string t = trim(line);
    if (t.empty() || t[0] == '#') continue;
    const auto eq = line.find('=');
    if (eq == std::string::npos) return R::failure(-1, path + ":" + std::to_string(n) + ": falta '='");
    const std::string key = trim(line.substr(0, eq));
    if (key.empty() || !std::all_of(key.begin(), key.end(), [](unsigned char c) { return std::isalnum(c) || c == '_'; }))
      return R::failure(-1, path + ":" + std::to_string(n) + ": clave inválida");
    out[key] = line.substr(eq + 1);  // valor literal, sin expansión
  }
  return out;
}

std::string Secrets::defaultPath() {
  const char* env = std::getenv("DAQV2_SECRETS");
  if (env && *env) return homeExpand(env);
  return homeExpand("~/.config/daqv2/secrets.env");
}

Secrets::Secrets(std::string filePath) : path_(homeExpand(std::move(filePath))) {}

Status Secrets::loadFile() {
  if (loaded_) return fileError_ ? Status(*fileError_) : okStatus();
  loaded_ = true;
  struct stat st {};
  if (::stat(path_.c_str(), &st) != 0) return okStatus();  // sin archivo: solo entorno
  auto r = parseEnvFile(path_);
  if (!r) {
    fileError_ = r.error();
    return Status(r.error());
  }
  values_ = std::move(r).value();
  return okStatus();
}

Result<std::optional<std::string>> Secrets::get(const std::string& key) {
  using R = Result<std::optional<std::string>>;
  if (const char* v = std::getenv(key.c_str()); v != nullptr && *v != '\0') return std::optional<std::string>(v);
  if (auto st = loadFile(); !st) return R(st.error());
  auto it = values_.find(key);
  if (it == values_.end() || it->second.empty()) return std::optional<std::string>();
  return std::optional<std::string>(it->second);
}

bool isSecretKey(const std::string& key) {
  std::string u = key;
  std::transform(u.begin(), u.end(), u.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return u.find("PASSWORD") != std::string::npos;
}

nlohmann::json redactSecrets(const nlohmann::json& j) {
  if (j.is_object()) {
    nlohmann::json out = nlohmann::json::object();
    for (auto it = j.begin(); it != j.end(); ++it)
      out[it.key()] = isSecretKey(it.key()) ? nlohmann::json("***") : redactSecrets(it.value());
    return out;
  }
  if (j.is_array()) {
    nlohmann::json out = nlohmann::json::array();
    for (const auto& e : j) out.push_back(redactSecrets(e));
    return out;
  }
  return j;
}

}  // namespace daq
