// bounded_queue.h
// Cola FIFO de capacidad fija entre hilos (productores múltiples, consumidores múltiples).
// - push(): bloquea mientras la cola esté llena; devuelve false si la cola está cerrada.
// - try_push(): nunca bloquea; con la cola llena aplica la política de descarte explícita
//   (DropNewest descarta el elemento nuevo, DropOldest el más antiguo) y cuenta el descarte.
// - pop(): bloquea hasta que haya un elemento; tras close() entrega lo que quede y luego
//   devuelve std::nullopt.
// - close(): despierta a todos los que esperan; push posteriores fallan.
// Seguridad entre hilos: todos los métodos se pueden llamar desde cualquier hilo. Las esperas
// usan predicado y las notificaciones se hacen después de soltar el mutex.
#pragma once

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <utility>

namespace daq {

enum class DropPolicy { DropNewest, DropOldest };

template <class T>
class BoundedQueue {
 public:
  explicit BoundedQueue(std::size_t capacity) : capacity_(capacity) {
    if (capacity_ == 0) throw std::invalid_argument("BoundedQueue: capacidad 0");
  }
  BoundedQueue(const BoundedQueue&) = delete;
  BoundedQueue& operator=(const BoundedQueue&) = delete;

  bool push(T item) {
    {
      std::unique_lock<std::mutex> lk(m_);
      notFull_.wait(lk, [&] { return closed_ || items_.size() < capacity_; });
      if (closed_) return false;
      items_.push_back(std::move(item));
    }
    notEmpty_.notify_one();
    return true;
  }

  // Devuelve true si el elemento nuevo quedó en la cola.
  bool try_push(T item, DropPolicy policy) {
    {
      std::lock_guard<std::mutex> lk(m_);
      if (closed_) return false;
      if (items_.size() >= capacity_) {
        ++dropped_;
        if (policy == DropPolicy::DropNewest) return false;
        items_.pop_front();
      }
      items_.push_back(std::move(item));
    }
    notEmpty_.notify_one();
    return true;
  }

  std::optional<T> pop() {
    std::optional<T> out;
    {
      std::unique_lock<std::mutex> lk(m_);
      notEmpty_.wait(lk, [&] { return closed_ || !items_.empty(); });
      if (items_.empty()) return std::nullopt;  // cerrada y vacía
      out.emplace(std::move(items_.front()));
      items_.pop_front();
    }
    notFull_.notify_one();
    return out;
  }

  std::optional<T> try_pop() {
    std::optional<T> out;
    {
      std::lock_guard<std::mutex> lk(m_);
      if (items_.empty()) return std::nullopt;
      out.emplace(std::move(items_.front()));
      items_.pop_front();
    }
    notFull_.notify_one();
    return out;
  }

  void close() {
    {
      std::lock_guard<std::mutex> lk(m_);
      closed_ = true;
    }
    notEmpty_.notify_all();
    notFull_.notify_all();
  }

  bool closed() const {
    std::lock_guard<std::mutex> lk(m_);
    return closed_;
  }
  std::size_t size() const {
    std::lock_guard<std::mutex> lk(m_);
    return items_.size();
  }
  std::size_t capacity() const noexcept { return capacity_; }
  std::uint64_t dropped() const {
    std::lock_guard<std::mutex> lk(m_);
    return dropped_;
  }

 private:
  const std::size_t capacity_;
  mutable std::mutex m_;
  std::condition_variable notFull_;
  std::condition_variable notEmpty_;
  std::deque<T> items_;
  bool closed_ = false;
  std::uint64_t dropped_ = 0;
};

}  // namespace daq
