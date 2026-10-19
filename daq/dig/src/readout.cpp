// readout.cpp
#include "daq/dig/readout.h"

#include <algorithm>
#include <chrono>

#include "daq/core/logger.h"

namespace daq::dig {

namespace {
void atomicMax(std::atomic<std::uint64_t>& a, std::uint64_t v) {
  std::uint64_t cur = a.load();
  while (v > cur && !a.compare_exchange_weak(cur, v)) {
  }
}
}  // namespace

std::string ReadoutStats::finishReason() const {
  std::lock_guard<std::mutex> lk(m_);
  return finishReason_;
}
void ReadoutStats::setFinishReason(const std::string& s) {
  std::lock_guard<std::mutex> lk(m_);
  finishReason_ = s;
}
std::string ReadoutStats::lastError() const {
  std::lock_guard<std::mutex> lk(m_);
  return lastError_;
}
void ReadoutStats::setLastError(const std::string& s) {
  std::lock_guard<std::mutex> lk(m_);
  lastError_ = s;
}

nlohmann::json ReadoutStats::toJson() const {
  return {{"blocks", blocks.load()},
          {"events", events.load()},
          {"empty_polls", emptyPolls.load()},
          {"read_errors", readErrors.load()},
          {"writer_blocked", writerBlocked.load()},
          {"writer_blocked_s", static_cast<double>(writerBlockedNs.load()) / 1e9},
          {"writer_max_occupancy", writerMaxOccupancy.load()},
          {"monitor_dropped", monitorDropped.load()},
          {"monitor_max_occupancy", monitorMaxOccupancy.load()},
          {"board_full_blocks", boardFullBlocks.load()},
          {"max_block_events", maxBlockEvents.load()},
          {"sw_triggers_sent", swTriggersSent.load()},
          {"finished", finished.load()},
          {"finish_reason", finishReason()},
          {"last_error", lastError()}};
}

ReadoutLoop::ReadoutLoop(DigitizerBackend& backend, Clock& clock, BoundedQueue<EventBlockPtr>& writerQ,
                         BoundedQueue<EventBlockPtr>& monitorQ, ReadoutOptions opt, ReadoutStats& stats)
    : be_(backend), clock_(clock), writerQ_(writerQ), monitorQ_(monitorQ), opt_(opt), st_(stats) {}

void ReadoutLoop::publish(EventBlockPtr b) {
  ++st_.blocks;
  st_.events += b->events.size();
  atomicMax(st_.maxBlockEvents, b->events.size());
  if (b->boardFull) ++st_.boardFullBlocks;
  // Monitor primero (sin bloqueo); luego el escritor, que puede bloquear.
  if (monitorQ_.size() >= monitorQ_.capacity()) ++st_.monitorDropped;
  monitorQ_.try_push(b, DropPolicy::DropOldest);
  atomicMax(st_.monitorMaxOccupancy, monitorQ_.size());
  const bool full = writerQ_.size() >= writerQ_.capacity();
  const auto t0 = std::chrono::steady_clock::now();
  writerQ_.push(std::move(b));
  if (full) {
    ++st_.writerBlocked;
    st_.writerBlockedNs += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count());
  }
  atomicMax(st_.writerMaxOccupancy, writerQ_.size());
}

Status ReadoutLoop::run(const StopToken& st) {
  Logger::setThreadName("readout");
  const auto pollSleep = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(opt_.pollSleepS));
  const SteadyTime t0 = clock_.now_steady();
  SteadyTime nextSw = t0;
  const auto swPeriod = std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::duration<double>(opt_.swTriggerRateHz > 0 ? 1.0 / opt_.swTriggerRateHz : 1.0));
  Status result = okStatus();
  int consecutiveErrors = 0;
  std::string reason = "parada pedida";
  while (!st.stop_requested()) {
    if (opt_.maxEvents > 0 && st_.events.load() >= opt_.maxEvents) {
      reason = "max_events";
      break;
    }
    if (opt_.maxDurationS > 0 && toSeconds(clock_.now_steady() - t0) >= opt_.maxDurationS) {
      reason = "max_duration_s";
      break;
    }
    if (opt_.softwareTrigger) {
      const SteadyTime now = clock_.now_steady();
      while (nextSw <= now) {
        if (be_.sendSoftwareTrigger()) ++st_.swTriggersSent;
        nextSw += swPeriod;
      }
    }
    auto blk = std::make_shared<EventBlock>();
    std::size_t want = opt_.maxEventsBlt;
    if (opt_.maxEvents > 0) want = std::min<std::uint64_t>(want, opt_.maxEvents - st_.events.load());
    auto r = be_.readBlock(*blk, want);
    if (!r) {
      ++st_.readErrors;
      st_.setLastError(r.error().message);
      logError("lectura del digitizador: ", r.error().message);
      if (++consecutiveErrors >= 10) {
        result = Status(r.error());
        reason = "error de lectura";
        break;
      }
      clock_.sleep_for(std::chrono::milliseconds(100), st);
      continue;
    }
    consecutiveErrors = 0;
    if (r.value() == 0) {
      if (be_.exhausted()) {
        reason = "fuente agotada";
        break;
      }
      ++st_.emptyPolls;
      SteadyTime wake = clock_.now_steady() + pollSleep;
      if (opt_.softwareTrigger) wake = std::min(wake, nextSw);
      clock_.sleep_until(wake, st);
      continue;
    }
    publish(std::move(blk));
  }
  // Parada: stop() y vaciado de la memoria de la placa hasta que no quede nada.
  if (auto s = be_.stop(); !s) logError("stop del digitizador: ", s.error().message);
  for (int i = 0; i < 100000; ++i) {
    auto blk = std::make_shared<EventBlock>();
    auto r = be_.readBlock(*blk, opt_.maxEventsBlt);
    if (!r || r.value() == 0) break;
    publish(std::move(blk));
  }
  st_.setFinishReason(reason);
  st_.finished = true;
  writerQ_.close();
  monitorQ_.close();
  logInfo("lectura terminada (", reason, "): ", st_.events.load(), " eventos en ", st_.blocks.load(), " bloques");
  return result;
}

}  // namespace daq::dig
