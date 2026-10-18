// clock.cpp
#include "daq/core/clock.h"

namespace daq {

bool SystemClock::sleep_until(SteadyTime t, const StopToken& token) {
  return sleep_until_with_stop(t, token);
}

ManualClock::ManualClock(UtcNs startUtc) : utc_(startUtc), steady_(SteadyTime{}) {}

UtcNs ManualClock::now_utc_ns() const {
  std::lock_guard<std::mutex> lk(m_);
  return utc_;
}

SteadyTime ManualClock::now_steady() const {
  std::lock_guard<std::mutex> lk(m_);
  return steady_;
}

bool ManualClock::sleep_until(SteadyTime t, const StopToken& token) {
  if (token.stop_requested()) return false;
  std::lock_guard<std::mutex> lk(m_);
  if (t > steady_) {
    const auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(t - steady_);
    steady_ = t;
    utc_ += d.count();
  }
  return true;
}

void ManualClock::advance(std::chrono::nanoseconds d) {
  std::lock_guard<std::mutex> lk(m_);
  if (d.count() <= 0) return;
  steady_ += d;
  utc_ += d.count();
}

}  // namespace daq
