// sc_reader.h
// Lector de los archivos de control lento de hvmon (diseño 4.3) a un TTree de ROOT (y, a través de
// él, a RDataFrame: ROOT::RDataFrame df(*tree)). Biblioteca compartida libdaq_screader.so, para
// macros de análisis:
//   R__ADD_INCLUDE_PATH(daq/io/include)
//   R__LOAD_LIBRARY(daq/build/io/libdaq_screader.so)
//   TTree* t = daq::io::readSlowControl("/data/daq/slowcontrol", "hv", "20261019", "20261019");
// Acepta archivos planos y .gz, partes .N (otras columnas) y los .schema.json de cada parte:
//  - CSV (hv, hv_full, temp): ramas con el tipo del esquema (int64 → L, uint64 → l, int → I,
//    uint32 → i, float64 → D); la unión de las columnas de todas las partes; una columna ausente
//    en una parte vale NaN (D) o 0 (enteros). Rama extra `part` (1, 2, ...).
//  - JSONL hv_changes: t_utc_ns, t_prev_read_utc_ns, id, ch, type, param, old_value, new_value,
//    origin, command_id; notes: t_utc_ns, user, host, text.
// Las filas de las partes de un día se intercalan por t_resp_utc_ns (o t_utc_ns) leyendo en
// paralelo línea por línea: nunca se carga un archivo entero. Las fechas van de from a to (YYYYMMDD,
// inclusive). Devuelve un TTree sin directorio (en memoria) que pasa a ser del llamador, o nullptr y
// el motivo en *error.
// Seguridad entre hilos: sin estado global; el TTree devuelto pertenece al hilo que llama.
#pragma once

#include <string>

class TTree;

namespace daq::io {

TTree* readSlowControl(const std::string& dir, const std::string& kind, const std::string& fromDate,
                       const std::string& toDate, std::string* error = nullptr);

}  // namespace daq::io
