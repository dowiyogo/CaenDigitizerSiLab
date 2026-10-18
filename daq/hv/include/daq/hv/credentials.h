// credentials.h
// Credenciales del SY5527. Mecanismo vigente (DAQ2-3): variables DAQ_HV_HOST, DAQ_HV_USER y
// DAQ_HV_PASSWORD del entorno o de ~/.config/daqv2/secrets.env (ver daq/core/secrets.h), que crea
// René. Respaldo OBSOLETO: el JSON {"user": "...", "password": "..."} de
// ~/.config/daqv2/hv_credentials.json, con un aviso; se quitará cuando secrets.env funcione.
// Ambos archivos se rechazan si no son del usuario actual con modo 600. La contraseña nunca se
// escribe en logs, argv ni salidas.
#pragma once

#include <string>

#include "daq/core/result.h"
#include "daq/core/secrets.h"

namespace daq::hv {

struct Credentials {
  std::string user;
  std::string password;
};

std::string defaultCredentialsPath();
std::string expandHome(const std::string& path);  // "~/x" -> "$HOME/x"
Result<Credentials> loadCredentials(const std::string& path);

struct HvLogin {
  std::string host;
  Credentials credentials;
  std::string source;  // "secrets" o "hv_credentials.json (obsoleto)"
};
// host: DAQ_HV_HOST si existe, si no defaultHost. Usuario y contraseña: DAQ_HV_USER y
// DAQ_HV_PASSWORD; si faltan, el archivo heredado con aviso de obsolescencia.
Result<HvLogin> resolveHvLogin(Secrets& secrets, const std::string& defaultHost, const std::string& legacyPath);

}  // namespace daq::hv
