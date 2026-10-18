// credentials.h
// Credenciales del SY5527 desde un archivo JSON fuera del repositorio:
//   {"user": "...", "password": "..."}
// Ruta por defecto: ~/.config/daqv2/hv_credentials.json. Se rechaza el archivo si no es un
// archivo regular del usuario actual o si su modo no es exactamente 600. La contraseña nunca
// se escribe en logs ni en la salida.
#pragma once

#include <string>

#include "daq/core/result.h"

namespace daq::hv {

struct Credentials {
  std::string user;
  std::string password;
};

std::string defaultCredentialsPath();
std::string expandHome(const std::string& path);  // "~/x" -> "$HOME/x"
Result<Credentials> loadCredentials(const std::string& path);

}  // namespace daq::hv
