// compressor.h
// Compresión de los días cerrados de control lento (fase 2). Al rotar a medianoche UTC el escritor
// encola los CSV/JSONL del día anterior; este hilo, con prioridad baja (nice 19), los comprime con
// gzipFileVerified (zlib en proceso): escribe <archivo>.gz, lo descomprime y compara bytes, líneas
// y CRC32, y solo entonces borra el archivo plano. Los .schema.json quedan planos. Si la
// verificación falla, el archivo plano queda y el .gz se renombra a .gz.failed.
// Al arrancar y en cada cambio de fecha (gancho del escritor, DAQ2-9) también se encolan los días
// cerrados que hayan quedado sin comprimir; uno ya comprimido (plano ausente, .gz presente) se omite.
// Seguridad entre hilos: una instancia pertenece al hilo compresor; publica en StateHub::compressor.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/gzip.h"
#include "daq/core/jthread.h"
#include "daq/hvmon/state.h"

namespace daq::hvmon {

class Compressor {
 public:
  explicit Compressor(StateHub& hub);
  void run(BoundedQueue<std::string>& queue, const StopToken& st);
  Result<GzipResult> compressOne(const std::string& path);
  // Archivos diarios planos (*.csv, *.jsonl con fecha YYYYMMDD) con fecha anterior a todayKey.
  static std::vector<std::string> closedDailyFiles(const std::string& baseDir, const std::string& todayKey);
  nlohmann::json status() const;

 private:
  StateHub& hub_;
  std::uint64_t done_ = 0, failed_ = 0, skipped_ = 0, inBytes_ = 0, outBytes_ = 0;
  nlohmann::json last_;
};

}  // namespace daq::hvmon
