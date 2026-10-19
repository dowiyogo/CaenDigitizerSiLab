// gzip.h
// Compresión verificada de archivos de control lento cerrados y lectura transparente de archivos
// planos o .gz.
// - gzipFileVerified(): comprime src a dst (formato gzip, zlib en proceso) escribiendo primero
//   dst + ".tmp" y renombrando; luego descomprime dst y compara bytes, líneas y CRC32 con src. Solo
//   si todo coincide y removeSource es true borra src. Si la verificación falla, src queda intacto
//   y el .gz defectuoso se renombra a dst + ".failed" (nunca se borra nada sin verificar).
// - LineReader: lee líneas de un archivo plano o .gz (gzopen lee ambos de forma transparente).
// Funciones sin estado compartido; seguras entre hilos para archivos distintos.
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>

#include "daq/core/result.h"

namespace daq {

struct FileDigest {
  std::uint64_t bytes = 0;
  std::uint64_t lines = 0;  // cantidad de '\n'
  std::uint32_t crc32 = 0;
};

struct GzipResult {
  FileDigest source;
  std::uint64_t compressedBytes = 0;
  std::string output;
  double ratio() const {
    return source.bytes ? static_cast<double>(compressedBytes) / static_cast<double>(source.bytes) : 0.0;
  }
};

// Resumen de un archivo plano o .gz (contenido descomprimido).
Result<FileDigest> digestFile(const std::string& path);

// beforeVerify: solo para tests (permite corromper el .gz antes de la verificación).
Result<GzipResult> gzipFileVerified(const std::string& src, const std::string& dst, bool removeSource,
                                    int level = 6,
                                    const std::function<void(const std::string&)>& beforeVerify = {});

class LineReader {
 public:
  static Result<std::unique_ptr<LineReader>> open(const std::string& path);
  ~LineReader();
  LineReader(const LineReader&) = delete;
  LineReader& operator=(const LineReader&) = delete;
  // Siguiente línea sin '\n'; nullopt al final del archivo.
  std::optional<std::string> next();

 private:
  LineReader() = default;
  void* gz_ = nullptr;  // gzFile
};

}  // namespace daq
