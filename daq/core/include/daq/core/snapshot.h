// snapshot.h
// Último valor publicado, con número de secuencia creciente (0 = nunca publicado).
// publish() y read() copian el valor bajo un std::scoped_lock corto; un lector siempre ve un
// par (valor, secuencia) consistente. Seguro entre hilos. T debe ser copiable.
#pragma once

#include <cstdint>
#include <mutex>
#include <utility>

namespace daq {

template <class T>
class Snapshot {
 public:
  Snapshot() = default;
  explicit Snapshot(T initial) : value_(std::move(initial)) {}

  // Devuelve la secuencia asignada al valor publicado.
  std::uint64_t publish(T value) {
    std::scoped_lock lk(m_);
    value_ = std::move(value);
    return ++seq_;
  }

  std::pair<T, std::uint64_t> read() const {
    std::scoped_lock lk(m_);
    return {value_, seq_};
  }

  std::uint64_t sequence() const {
    std::scoped_lock lk(m_);
    return seq_;
  }

 private:
  mutable std::mutex m_;
  T value_{};
  std::uint64_t seq_ = 0;
};

}  // namespace daq
