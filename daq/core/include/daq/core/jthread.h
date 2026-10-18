// jthread.h
// Reemplazo C++17 de std::jthread, std::stop_source/stop_token/stop_callback y de la espera
// de condition_variable_any con stop_token. ROOT 6.32.10 está compilado con C++17, por eso no
// se usa la versión estándar (docs/daq_v2_design.md, sección 15, punto 1).
//
// Garantías:
// - StopSource/StopToken: request_stop() es seguro desde cualquier hilo e idempotente; solo la
//   primera llamada devuelve true y ejecuta los callbacks registrados, en el hilo que la llama.
// - StopCallback: si la parada ya fue pedida al construirlo, el callback corre en el
//   constructor. El destructor espera a que termine un callback en curso en otro hilo.
// - JThread: el destructor pide parada y hace join. Si el callable acepta un StopToken como
//   primer argumento, se le pasa el token de la hebra.
// - wait_with_stop(): espera en un std::condition_variable hasta que el predicado se cumpla o
//   se pida parada. Regla: no llamar request_stop() mientras se tiene tomado el mutex usado en
//   esa espera (el callback lo toma para notificar sin perder el despertar).
#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <list>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>

namespace daq {

namespace detail {

class StopCallbackBase {
 public:
  virtual void invoke() noexcept = 0;

 protected:
  ~StopCallbackBase() = default;
};

class StopState;

}  // namespace detail

class StopToken {
 public:
  StopToken() noexcept = default;
  bool stop_requested() const noexcept;
  bool stop_possible() const noexcept { return static_cast<bool>(state_); }

 private:
  friend class StopSource;
  template <class F>
  friend class StopCallback;
  explicit StopToken(std::shared_ptr<detail::StopState> s) noexcept : state_(std::move(s)) {}
  std::shared_ptr<detail::StopState> state_;
};

class StopSource {
 public:
  StopSource();
  bool request_stop() noexcept;
  bool stop_requested() const noexcept;
  StopToken get_token() const noexcept { return StopToken(state_); }

 private:
  std::shared_ptr<detail::StopState> state_;
};

namespace detail {

class StopState {
 public:
  bool stopRequested() const noexcept { return stopped_.load(std::memory_order_acquire); }
  bool requestStop() noexcept;
  // Devuelve false si la parada ya fue pedida (el llamador ejecuta el callback).
  bool add(StopCallbackBase* cb);
  void remove(StopCallbackBase* cb) noexcept;

 private:
  std::atomic<bool> stopped_{false};
  std::mutex m_;
  std::condition_variable done_;
  std::list<StopCallbackBase*> callbacks_;
  StopCallbackBase* running_ = nullptr;
  std::thread::id runningThread_;
};

}  // namespace detail

template <class F>
class StopCallback final : private detail::StopCallbackBase {
 public:
  StopCallback(const StopToken& token, F f) : f_(std::move(f)), state_(token.state_) {
    if (state_ && !state_->add(this)) {
      state_.reset();
      f_();
    }
  }
  ~StopCallback() {
    if (state_) state_->remove(this);
  }
  StopCallback(const StopCallback&) = delete;
  StopCallback& operator=(const StopCallback&) = delete;

 private:
  void invoke() noexcept override { f_(); }
  F f_;
  std::shared_ptr<detail::StopState> state_;
};

class JThread {
 public:
  JThread() noexcept = default;

  template <class F, class... Args,
            class = std::enable_if_t<!std::is_same_v<std::decay_t<F>, JThread>>>
  explicit JThread(F&& f, Args&&... args) : source_() {
    if constexpr (std::is_invocable_v<std::decay_t<F>, StopToken, std::decay_t<Args>...>) {
      thread_ = std::thread(std::forward<F>(f), source_.get_token(), std::forward<Args>(args)...);
    } else {
      thread_ = std::thread(std::forward<F>(f), std::forward<Args>(args)...);
    }
  }

  ~JThread() { stopAndJoin(); }

  JThread(JThread&&) noexcept = default;
  JThread& operator=(JThread&& other) noexcept {
    if (this != &other) {
      stopAndJoin();
      source_ = std::move(other.source_);
      thread_ = std::move(other.thread_);
    }
    return *this;
  }
  JThread(const JThread&) = delete;
  JThread& operator=(const JThread&) = delete;

  bool joinable() const noexcept { return thread_.joinable(); }
  void join() { thread_.join(); }
  bool request_stop() noexcept { return source_.request_stop(); }
  StopSource get_stop_source() const noexcept { return source_; }
  StopToken get_stop_token() const noexcept { return source_.get_token(); }

 private:
  void stopAndJoin() {
    if (thread_.joinable()) {
      source_.request_stop();
      thread_.join();
    }
  }
  StopSource source_;
  std::thread thread_;
};

// Espera hasta que pred() sea verdadero o se pida parada. Devuelve pred() con el lock tomado.
template <class Pred>
bool wait_with_stop(std::condition_variable& cv, std::unique_lock<std::mutex>& lock,
                    const StopToken& token, Pred pred) {
  for (;;) {
    if (pred()) return true;
    if (token.stop_requested()) return false;
    {
      // El callback toma el mutex antes de notificar: así el aviso no se pierde entre la
      // comprobación del predicado y el bloqueo en wait().
      std::mutex* m = lock.mutex();
      StopCallback cb(token, [m, &cv] {
        { std::lock_guard<std::mutex> g(*m); }
        cv.notify_all();
      });
      cv.wait(lock, [&] { return pred() || token.stop_requested(); });
      lock.unlock();  // el callback se destruye sin el mutex tomado (evita deadlock)
    }
    lock.lock();
  }
}

// Como wait_with_stop, con plazo. Devuelve pred() al despertar (por condición, parada o plazo).
template <class Clock, class Duration, class Pred>
bool wait_until_with_stop(std::condition_variable& cv, std::unique_lock<std::mutex>& lock,
                          const StopToken& token,
                          const std::chrono::time_point<Clock, Duration>& deadline, Pred pred) {
  for (;;) {
    if (pred()) return true;
    if (token.stop_requested() || Clock::now() >= deadline) return pred();
    {
      std::mutex* m = lock.mutex();
      StopCallback cb(token, [m, &cv] {
        { std::lock_guard<std::mutex> g(*m); }
        cv.notify_all();
      });
      cv.wait_until(lock, deadline, [&] { return pred() || token.stop_requested(); });
      lock.unlock();
    }
    lock.lock();
  }
}

// Duerme hasta el instante dado o hasta que se pida parada. Devuelve false si hubo parada.
bool sleep_until_with_stop(std::chrono::steady_clock::time_point deadline, const StopToken& token);

}  // namespace daq
