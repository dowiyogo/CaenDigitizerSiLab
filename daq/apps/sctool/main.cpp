// sctool: utilidades para los archivos de control lento de hvmon.
//   sctool compress [--keep] ARCHIVO...  gzip verificado (bytes, líneas y CRC32); sin --keep borra
//                                        el plano solo tras verificar. Misma rutina que hvmon.
//   sctool cat ARCHIVO                   imprime un archivo plano o .gz
//   sctool digest ARCHIVO                bytes, líneas y CRC32 del contenido (plano o .gz)
#include <cstdio>
#include <string>
#include <vector>

#include "daq/core/gzip.h"

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "Uso: sctool compress [--keep] ARCHIVO... | cat ARCHIVO | digest ARCHIVO\n");
    return 2;
  }
  const std::string cmd = argv[1];
  if (cmd == "cat" || cmd == "digest") {
    if (cmd == "digest") {
      auto d = daq::digestFile(argv[2]);
      if (!d) {
        std::fprintf(stderr, "sctool: %s\n", d.error().message.c_str());
        return 1;
      }
      std::printf("%s bytes=%llu lineas=%llu crc32=%08x\n", argv[2], static_cast<unsigned long long>(d.value().bytes),
                  static_cast<unsigned long long>(d.value().lines), d.value().crc32);
      return 0;
    }
    auto r = daq::LineReader::open(argv[2]);
    if (!r) {
      std::fprintf(stderr, "sctool: %s\n", r.error().message.c_str());
      return 1;
    }
    while (auto l = r.value()->next()) std::printf("%s\n", l->c_str());
    return 0;
  }
  if (cmd != "compress") {
    std::fprintf(stderr, "sctool: subcomando desconocido %s\n", cmd.c_str());
    return 2;
  }
  bool keep = false;
  std::vector<std::string> files;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--keep") keep = true;
    else files.push_back(a);
  }
  int rc = 0;
  unsigned long long in = 0, out = 0;
  for (const auto& f : files) {
    auto r = daq::gzipFileVerified(f, f + ".gz", !keep);
    if (!r) {
      std::fprintf(stderr, "sctool: %s: %s\n", f.c_str(), r.error().message.c_str());
      rc = 1;
      continue;
    }
    in += r.value().source.bytes;
    out += r.value().compressedBytes;
    std::printf("%s: %llu -> %llu bytes (%.4f), %llu lineas, crc32 %08x verificado%s\n", f.c_str(),
                static_cast<unsigned long long>(r.value().source.bytes),
                static_cast<unsigned long long>(r.value().compressedBytes), r.value().ratio(),
                static_cast<unsigned long long>(r.value().source.lines), r.value().source.crc32,
                keep ? "" : ", plano borrado");
  }
  if (files.size() > 1 && in > 0)
    std::printf("total: %llu -> %llu bytes (%.4f)\n", in, out, static_cast<double>(out) / static_cast<double>(in));
  return rc;
}
