// gzip.cpp
#include "daq/core/gzip.h"

#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <vector>

#include "daq/core/logger.h"

namespace daq {

namespace {

constexpr std::size_t kChunk = 1 << 16;

std::string zerr(gzFile f) {
  int code = 0;
  const char* m = gzerror(f, &code);
  return m ? m : "error de zlib";
}

}  // namespace

Result<FileDigest> digestFile(const std::string& path) {
  gzFile f = gzopen(path.c_str(), "rb");
  if (!f) return Result<FileDigest>::failure(errno, "no se pudo abrir " + path + ": " + std::strerror(errno));
  FileDigest d;
  uLong crc = crc32(0L, Z_NULL, 0);
  std::vector<unsigned char> buf(kChunk);
  for (;;) {
    const int n = gzread(f, buf.data(), static_cast<unsigned>(buf.size()));
    if (n < 0) {
      const std::string m = zerr(f);
      gzclose(f);
      return Result<FileDigest>::failure(-1, "lectura de " + path + ": " + m);
    }
    if (n == 0) break;
    const auto un = static_cast<std::size_t>(n);
    d.bytes += un;
    for (std::size_t i = 0; i < un; ++i)
      if (buf[i] == '\n') ++d.lines;
    crc = crc32(crc, buf.data(), static_cast<uInt>(un));
  }
  gzclose(f);
  d.crc32 = static_cast<std::uint32_t>(crc);
  return d;
}

Result<GzipResult> gzipFileVerified(const std::string& src, const std::string& dst, bool removeSource, int level,
                                    const std::function<void(const std::string&)>& beforeVerify) {
  using R = Result<GzipResult>;
  auto srcDigest = digestFile(src);
  if (!srcDigest) return R(srcDigest.error());
  FILE* in = std::fopen(src.c_str(), "rb");
  if (!in) return R::failure(errno, "no se pudo abrir " + src + ": " + std::strerror(errno));
  const std::string tmp = dst + ".tmp";
  const std::string mode = "wb" + std::to_string(level < 1 || level > 9 ? 6 : level);
  gzFile out = gzopen(tmp.c_str(), mode.c_str());
  if (!out) {
    std::fclose(in);
    return R::failure(errno, "no se pudo crear " + tmp);
  }
  std::vector<char> buf(kChunk);
  bool ok = true;
  std::string why;
  for (;;) {
    const std::size_t n = std::fread(buf.data(), 1, buf.size(), in);
    if (n > 0 && gzwrite(out, buf.data(), static_cast<unsigned>(n)) != static_cast<int>(n)) {
      ok = false;
      why = "escritura de " + tmp + ": " + zerr(out);
      break;
    }
    if (n < buf.size()) {
      if (std::ferror(in)) {
        ok = false;
        why = "lectura de " + src;
      }
      break;
    }
  }
  std::fclose(in);
  if (gzclose(out) != Z_OK && ok) {
    ok = false;
    why = "cierre de " + tmp;
  }
  if (!ok) {
    std::rename(tmp.c_str(), (dst + ".failed").c_str());
    return R::failure(-1, why);
  }
  // El .gz conserva el modo (y, si se puede, el dueño) del original; gzopen crea con 0666 & umask.
  struct stat st {};
  if (::stat(src.c_str(), &st) == 0) {
    if (::chmod(tmp.c_str(), st.st_mode & 07777) != 0) {
      std::rename(tmp.c_str(), (dst + ".failed").c_str());
      return R::failure(errno, "chmod " + tmp + ": " + std::strerror(errno));
    }
    if (::chown(tmp.c_str(), st.st_uid, st.st_gid) != 0 && errno != EPERM)
      logWarn("gzip: no se pudo copiar el dueño de ", src, " a ", tmp, ": ", std::strerror(errno));
  }
  if (std::rename(tmp.c_str(), dst.c_str()) != 0) return R::failure(errno, "rename " + tmp + ": " + std::strerror(errno));
  if (beforeVerify) beforeVerify(dst);

  // Verificación: el contenido descomprimido debe ser idéntico (bytes, líneas, CRC32).
  auto gzDigest = digestFile(dst);
  const FileDigest& s = srcDigest.value();
  if (!gzDigest || gzDigest.value().bytes != s.bytes || gzDigest.value().lines != s.lines ||
      gzDigest.value().crc32 != s.crc32) {
    std::rename(dst.c_str(), (dst + ".failed").c_str());
    return R::failure(-1, "verificación de " + dst + " falló: el original " + src + " se conserva" +
                              (gzDigest ? "" : " (" + gzDigest.error().message + ")"));
  }
  GzipResult res;
  res.source = s;
  res.output = dst;
  if (FILE* f = std::fopen(dst.c_str(), "rb")) {
    std::fseek(f, 0, SEEK_END);
    res.compressedBytes = static_cast<std::uint64_t>(std::ftell(f));
    std::fclose(f);
  }
  if (removeSource && ::unlink(src.c_str()) != 0)
    return R::failure(errno, "comprimido y verificado, pero no se pudo borrar " + src + ": " + std::strerror(errno));
  return res;
}

Result<std::unique_ptr<LineReader>> LineReader::open(const std::string& path) {
  gzFile f = gzopen(path.c_str(), "rb");
  if (!f) return Result<std::unique_ptr<LineReader>>::failure(errno, "no se pudo abrir " + path);
  std::unique_ptr<LineReader> r(new LineReader());
  r->gz_ = f;
  return r;
}

LineReader::~LineReader() {
  if (gz_) gzclose(static_cast<gzFile>(gz_));
}

std::optional<std::string> LineReader::next() {
  std::string line;
  char buf[4096];
  for (;;) {
    char* p = gzgets(static_cast<gzFile>(gz_), buf, sizeof buf);
    if (!p) return line.empty() ? std::nullopt : std::optional<std::string>(line);
    line += buf;
    if (!line.empty() && line.back() == '\n') {
      line.pop_back();
      return line;
    }
  }
}

}  // namespace daq
