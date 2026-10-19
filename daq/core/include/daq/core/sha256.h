// sha256.h
// SHA-256 (FIPS 180-4) propio, sin dependencias: hash de cada segmento ROOT cerrado, del índice
// del run, de la configuración y del mapeo. sha256File() lee el archivo por bloques (nunca entero).
// Seguridad entre hilos: cada instancia de Sha256 pertenece a un hilo; las funciones libres no
// tienen estado.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include "daq/core/result.h"

namespace daq {

class Sha256 {
 public:
  Sha256();
  void update(const void* data, std::size_t n);
  void update(std::string_view s) { update(s.data(), s.size()); }
  std::array<std::uint8_t, 32> digest();  // termina el cálculo; la instancia no se reutiliza
  std::string hexDigest();                 // 64 caracteres hexadecimales en minúscula

 private:
  void block(const std::uint8_t* p);
  std::array<std::uint32_t, 8> h_;
  std::array<std::uint8_t, 64> buf_{};
  std::size_t bufLen_ = 0;
  std::uint64_t totalBytes_ = 0;
};

std::string sha256Hex(std::string_view data);
Result<std::string> sha256File(const std::string& path);

}  // namespace daq
