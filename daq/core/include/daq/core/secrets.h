// secrets.h
// Secretos del sistema (usuario y contraseña del SY5527, cuenta SMTP). Los crea René en
// ~/.config/daqv2/secrets.env (modo 600, directorio 700, nunca en el repositorio).
// - get(): primero el entorno del proceso (getenv); si falta la variable, el archivo indicado por
//   DAQV2_SECRETS (por defecto ~/.config/daqv2/secrets.env), que se lee una sola vez.
// - El archivo se rechaza si no es un archivo regular del usuario actual con modo exactamente 600.
// - Formato: líneas KEY=VALUE, comentarios con '#', sin expansión de shell. El valor es el texto
//   exacto tras el primer '=' (solo se quita un '\r' final); la clave se recorta.
// - Los secretos nunca se escriben en logs, argv, archivos de estado ni respuestas: para volcar
//   configuraciones usar redactSecrets(), que reemplaza todo valor cuya clave contenga
//   "PASSWORD" (sin distinguir mayúsculas) por "***".
// Seguridad entre hilos: una instancia no es segura entre hilos; usarla durante el arranque.
#pragma once

#include <map>
#include <optional>
#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"

namespace daq {

// Comprueba que path sea un archivo regular del usuario actual con modo 600.
Status checkPrivateFile(const std::string& path);
// Lee un archivo KEY=VALUE (tras checkPrivateFile).
Result<std::map<std::string, std::string>> parseEnvFile(const std::string& path);

class Secrets {
 public:
  static std::string defaultPath();  // $DAQV2_SECRETS o ~/.config/daqv2/secrets.env
  explicit Secrets(std::string filePath = defaultPath());

  // nullopt si la variable no está ni en el entorno ni en el archivo (o el archivo no existe).
  // Error si el archivo existe pero es inseguro o ilegible.
  Result<std::optional<std::string>> get(const std::string& key);
  const std::string& path() const { return path_; }

 private:
  Status loadFile();
  std::string path_;
  bool loaded_ = false;
  std::optional<Error> fileError_;
  std::map<std::string, std::string> values_;
};

bool isSecretKey(const std::string& key);  // contiene "PASSWORD"
nlohmann::json redactSecrets(const nlohmann::json& j);

}  // namespace daq
