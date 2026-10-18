// stop_token.cpp
// Implementación de StopState (ver jthread.h).
#include "daq/core/jthread.h"

#include <algorithm>

namespace daq {

bool StopToken::stop_requested() const noexcept { return state_ && state_->stopRequested(); }

StopSource::StopSource() : state_(std::make_shared<detail::StopState>()) {}

bool StopSource::request_stop() noexcept { return state_ && state_->requestStop(); }

bool StopSource::stop_requested() const noexcept { return state_ && state_->stopRequested(); }

namespace detail {

bool StopState::requestStop() noexcept {
  std::unique_lock<std::mutex> lk(m_);
  if (stopped_.load(std::memory_order_relaxed)) return false;
  stopped_.store(true, std::memory_order_release);
  runningThread_ = std::this_thread::get_id();
  while (!callbacks_.empty()) {
    StopCallbackBase* cb = callbacks_.front();
    callbacks_.pop_front();
    running_ = cb;
    lk.unlock();
    cb->invoke();
    lk.lock();
    running_ = nullptr;
    done_.notify_all();
  }
  return true;
}

bool StopState::add(StopCallbackBase* cb) {
  std::lock_guard<std::mutex> lk(m_);
  if (stopped_.load(std::memory_order_relaxed)) return false;
  callbacks_.push_back(cb);
  return true;
}

void StopState::remove(StopCallbackBase* cb) noexcept {
  std::unique_lock<std::mutex> lk(m_);
  auto it = std::find(callbacks_.begin(), callbacks_.end(), cb);
  if (it != callbacks_.end()) {
    callbacks_.erase(it);
    return;
  }
  // El callback está corriendo en otro hilo: esperar a que termine antes de destruirlo.
  if (running_ == cb && runningThread_ != std::this_thread::get_id()) {
    done_.wait(lk, [&] { return running_ != cb; });
  }
}

}  // namespace detail

bool sleep_until_with_stop(std::chrono::steady_clock::time_point deadline, const StopToken& token) {
  std::mutex m;
  std::condition_variable cv;
  std::unique_lock<std::mutex> lk(m);
  wait_until_with_stop(cv, lk, token, deadline, [] { return false; });
  return !token.stop_requested();
}

}  // namespace daq
