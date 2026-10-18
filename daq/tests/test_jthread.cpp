// Tests de JThread, StopSource/StopToken/StopCallback y esperas con parada.
#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>

#include "daq/core/jthread.h"

using namespace std::chrono_literals;
using daq::JThread;
using daq::StopCallback;
using daq::StopSource;
using daq::StopToken;

TEST(StopToken, DefaultTokenCannotStop) {
  StopToken t;
  EXPECT_FALSE(t.stop_possible());
  EXPECT_FALSE(t.stop_requested());
}

TEST(StopSource, RequestStopIsIdempotent) {
  StopSource s;
  auto t = s.get_token();
  EXPECT_TRUE(t.stop_possible());
  EXPECT_FALSE(t.stop_requested());
  EXPECT_TRUE(s.request_stop());
  EXPECT_FALSE(s.request_stop());
  EXPECT_TRUE(t.stop_requested());
}

TEST(StopCallback, RunsOnStopAndImmediatelyIfAlreadyStopped) {
  StopSource s;
  int calls = 0;
  {
    StopCallback cb(s.get_token(), [&] { ++calls; });
    EXPECT_EQ(calls, 0);
    s.request_stop();
    EXPECT_EQ(calls, 1);
  }
  StopCallback late(s.get_token(), [&] { ++calls; });
  EXPECT_EQ(calls, 2);
}

TEST(StopCallback, DeregisteredCallbackDoesNotRun) {
  StopSource s;
  int calls = 0;
  { StopCallback cb(s.get_token(), [&] { ++calls; }); }
  s.request_stop();
  EXPECT_EQ(calls, 0);
}

TEST(JThread, DestructorRequestsStopAndJoins) {
  std::atomic<bool> finished{false};
  {
    JThread t([&](StopToken st) {
      while (!st.stop_requested()) std::this_thread::sleep_for(1ms);
      finished = true;
    });
    EXPECT_TRUE(t.joinable());
  }
  EXPECT_TRUE(finished.load());
}

TEST(JThread, CallableWithoutTokenIsJoined) {
  std::atomic<int> v{0};
  {
    JThread t([&](int x) { v = x; }, 42);
  }
  EXPECT_EQ(v.load(), 42);
}

TEST(JThread, StopBeforeStartIsSeenByThread) {
  std::atomic<bool> sawStop{false};
  std::mutex m;
  std::condition_variable cv;
  bool go = false;
  JThread t([&](StopToken st) {
    std::unique_lock<std::mutex> lk(m);
    cv.wait(lk, [&] { return go; });  // el cuerpo empieza después de request_stop()
    sawStop = st.stop_requested();
  });
  t.request_stop();
  {
    std::lock_guard<std::mutex> lk(m);
    go = true;
  }
  cv.notify_all();
  t.join();
  EXPECT_TRUE(sawStop.load());
}

TEST(JThread, StopWakesThreadWaitingOnConditionVariable) {
  std::mutex m;
  std::condition_variable cv;
  std::atomic<bool> returned{false};
  std::atomic<bool> predResult{true};
  JThread t([&](StopToken st) {
    std::unique_lock<std::mutex> lk(m);
    predResult = daq::wait_with_stop(cv, lk, st, [] { return false; });
    returned = true;
  });
  std::this_thread::sleep_for(50ms);
  EXPECT_FALSE(returned.load());
  const auto t0 = std::chrono::steady_clock::now();
  t.request_stop();
  t.join();
  EXPECT_TRUE(returned.load());
  EXPECT_FALSE(predResult.load());
  EXPECT_LT(std::chrono::steady_clock::now() - t0, 1s);
}

TEST(JThread, WaitWithStopReturnsTrueWhenPredicateHolds) {
  std::mutex m;
  std::condition_variable cv;
  bool ready = false;
  std::atomic<bool> result{false};
  JThread t([&](StopToken st) {
    std::unique_lock<std::mutex> lk(m);
    result = daq::wait_with_stop(cv, lk, st, [&] { return ready; });
  });
  std::this_thread::sleep_for(20ms);
  {
    std::lock_guard<std::mutex> lk(m);
    ready = true;
  }
  cv.notify_all();
  t.join();
  EXPECT_TRUE(result.load());
}

TEST(JThread, SleepUntilWithStopWakesEarly) {
  std::atomic<bool> slept{true};
  const auto t0 = std::chrono::steady_clock::now();
  JThread t([&](StopToken st) { slept = daq::sleep_until_with_stop(std::chrono::steady_clock::now() + 10s, st); });
  std::this_thread::sleep_for(20ms);
  t.request_stop();
  t.join();
  EXPECT_FALSE(slept.load());
  EXPECT_LT(std::chrono::steady_clock::now() - t0, 2s);
}

TEST(JThread, SleepUntilWithStopReachesDeadline) {
  StopSource s;
  const auto deadline = std::chrono::steady_clock::now() + 30ms;
  EXPECT_TRUE(daq::sleep_until_with_stop(deadline, s.get_token()));
  EXPECT_GE(std::chrono::steady_clock::now(), deadline);
}

TEST(JThread, MoveAssignmentStopsPreviousThread) {
  std::atomic<bool> firstDone{false};
  JThread a([&](StopToken st) {
    while (!st.stop_requested()) std::this_thread::sleep_for(1ms);
    firstDone = true;
  });
  a = JThread([](StopToken) {});
  EXPECT_TRUE(firstDone.load());
}
