// result.h
// Resultado de una operación que puede fallar: valor o error con código y mensaje.
// Los errores del hardware se propagan con este tipo; nunca como ceros silenciosos.
// Seguridad entre hilos: la de un valor; no hay estado compartido.
#pragma once

#include <string>
#include <utility>
#include <variant>

namespace daq {

struct Error {
  int code = 0;          // código del origen (p. ej. CAENHVRESULT); -1 si es interno
  std::string message;   // descripción legible, incluye el texto del Wrapper si lo hay
};

template <class T>
class Result {
 public:
  Result(T value) : v_(std::move(value)) {}            // NOLINT: conversión implícita deseada
  Result(Error error) : v_(std::move(error)) {}        // NOLINT
  static Result failure(int code, std::string message) { return Result(Error{code, std::move(message)}); }

  bool ok() const noexcept { return std::holds_alternative<T>(v_); }
  explicit operator bool() const noexcept { return ok(); }

  const T& value() const& { return std::get<T>(v_); }
  T& value() & { return std::get<T>(v_); }
  T&& value() && { return std::get<T>(std::move(v_)); }
  const Error& error() const { return std::get<Error>(v_); }

 private:
  std::variant<T, Error> v_;
};

using Status = Result<std::monostate>;
inline Status okStatus() { return Status(std::monostate{}); }

}  // namespace daq
