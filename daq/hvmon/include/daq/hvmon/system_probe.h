// system_probe.h
// Consultas al sistema operativo que hvmon necesita:
// - GECO2020 en ejecución: busca el patrón (por defecto "bin/CAENGECO2020", el ejecutable ELF del
//   snap caengeco2020) en el ejecutable de cada proceso (/proc/<pid>/exe) o en su argv[0]; nunca en
//   el resto de la línea de comandos, para no confundir con procesos que solo mencionan el texto
//   (p. ej. un grep o un script de muestreo). Excluye el proceso propio. Raíz de /proc inyectable.
// - Dirección de enlace del servidor HTTP: "127.0.0.1"/"localhost" o la IPv4 de tailscale0;
//   nunca 0.0.0.0 ni la interfaz de la universidad.
// Funciones sin estado, seguras entre hilos.
#pragma once

#include <string>

#include "daq/core/result.h"

namespace daq::hvmon {

bool processRunning(const std::string& pattern, const std::string& procRoot = "/proc");
Result<std::string> resolveHttpBind(const std::string& bind);

}  // namespace daq::hvmon
