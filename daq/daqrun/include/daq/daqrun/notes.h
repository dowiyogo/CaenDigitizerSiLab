// notes.h
// Notas del operador de un intervalo (diseño 4.3): lee <dir>/YYYY/notes_YYYYMMDD[.N].jsonl[.gz] de
// los días UTC que toca [from, to) y devuelve las notas con t_utc_ns en el intervalo, ordenadas.
// daqrun las copia al árbol `notes` de cada segmento al cerrarlo y al final del run.
// Función sin estado: segura entre hilos.
#pragma once

#include <string>
#include <vector>

#include "daq/core/time.h"

namespace daq::daqrun {

struct Note {
  UtcNs tUtcNs = 0;
  std::string user, host, text;
};

std::vector<Note> readNotes(const std::string& slowcontrolDir, UtcNs fromUtc, UtcNs toUtc);

}  // namespace daq::daqrun
