// Tests de BoundedQueue: bloqueo, cierre, vaciado, políticas de descarte y estrés MPSC.
#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

#include "daq/core/bounded_queue.h"

using namespace std::chrono_literals;
using daq::BoundedQueue;
using daq::DropPolicy;

TEST(BoundedQueue, RejectsZeroCapacity) { EXPECT_THROW(BoundedQueue<int>(0), std::invalid_argument); }

TEST(BoundedQueue, FifoOrderAndSize) {
  BoundedQueue<int> q(3);
  EXPECT_EQ(q.capacity(), 3u);
  EXPECT_TRUE(q.push(1));
  EXPECT_TRUE(q.push(2));
  EXPECT_EQ(q.size(), 2u);
  EXPECT_EQ(q.pop().value(), 1);
  EXPECT_EQ(q.pop().value(), 2);
  EXPECT_FALSE(q.try_pop().has_value());
}

TEST(BoundedQueue, PushBlocksWhileFull) {
  BoundedQueue<int> q(1);
  ASSERT_TRUE(q.push(1));
  std::atomic<bool> pushed{false};
  std::thread producer([&] {
    q.push(2);
    pushed = true;
  });
  std::this_thread::sleep_for(50ms);
  EXPECT_FALSE(pushed.load());
  EXPECT_EQ(q.pop().value(), 1);
  producer.join();
  EXPECT_TRUE(pushed.load());
  EXPECT_EQ(q.pop().value(), 2);
}

TEST(BoundedQueue, PopBlocksUntilPush) {
  BoundedQueue<int> q(2);
  std::atomic<bool> got{false};
  std::thread consumer([&] {
    auto v = q.pop();
    got = v.has_value() && *v == 7;
  });
  std::this_thread::sleep_for(50ms);
  EXPECT_FALSE(got.load());
  q.push(7);
  consumer.join();
  EXPECT_TRUE(got.load());
}

TEST(BoundedQueue, CloseWakesAllWaitingConsumers) {
  BoundedQueue<int> q(2);
  std::atomic<int> nullopts{0};
  std::vector<std::thread> consumers;
  for (int i = 0; i < 3; ++i)
    consumers.emplace_back([&] {
      if (!q.pop().has_value()) ++nullopts;
    });
  std::this_thread::sleep_for(50ms);
  q.close();
  for (auto& t : consumers) t.join();
  EXPECT_EQ(nullopts.load(), 3);
}

TEST(BoundedQueue, CloseWakesBlockedProducer) {
  BoundedQueue<int> q(1);
  q.push(1);
  std::atomic<int> result{-1};
  std::thread producer([&] { result = q.push(2) ? 1 : 0; });
  std::this_thread::sleep_for(50ms);
  q.close();
  producer.join();
  EXPECT_EQ(result.load(), 0);
}

TEST(BoundedQueue, DrainAfterClose) {
  BoundedQueue<int> q(4);
  q.push(1);
  q.push(2);
  q.close();
  EXPECT_TRUE(q.closed());
  EXPECT_FALSE(q.push(3));
  EXPECT_FALSE(q.try_push(3, DropPolicy::DropOldest));
  EXPECT_EQ(q.pop().value(), 1);
  EXPECT_EQ(q.pop().value(), 2);
  EXPECT_FALSE(q.pop().has_value());
}

TEST(BoundedQueue, DropNewestKeepsOldItems) {
  BoundedQueue<int> q(2);
  EXPECT_TRUE(q.try_push(1, DropPolicy::DropNewest));
  EXPECT_TRUE(q.try_push(2, DropPolicy::DropNewest));
  EXPECT_FALSE(q.try_push(3, DropPolicy::DropNewest));
  EXPECT_FALSE(q.try_push(4, DropPolicy::DropNewest));
  EXPECT_EQ(q.dropped(), 2u);
  EXPECT_EQ(q.pop().value(), 1);
  EXPECT_EQ(q.pop().value(), 2);
}

TEST(BoundedQueue, DropOldestKeepsNewItems) {
  BoundedQueue<int> q(2);
  for (int i = 1; i <= 5; ++i) EXPECT_TRUE(q.try_push(i, DropPolicy::DropOldest));
  EXPECT_EQ(q.dropped(), 3u);
  EXPECT_EQ(q.pop().value(), 4);
  EXPECT_EQ(q.pop().value(), 5);
}

TEST(BoundedQueue, MultiProducerStressNoLossOrderPerProducer) {
  constexpr int kProducers = 4;
  constexpr int kItems = 100000;
  BoundedQueue<std::pair<int, int>> q(64);
  std::vector<std::thread> producers;
  for (int p = 0; p < kProducers; ++p)
    producers.emplace_back([&q, p] {
      for (int i = 0; i < kItems; ++i) ASSERT_TRUE(q.push({p, i}));
    });
  std::vector<int> next(kProducers, 0);
  long total = 0;
  bool ordered = true;
  std::thread consumer([&] {
    while (auto item = q.pop()) {
      auto [p, i] = *item;
      if (i != next[static_cast<std::size_t>(p)]) ordered = false;
      next[static_cast<std::size_t>(p)] = i + 1;
      ++total;
    }
  });
  for (auto& t : producers) t.join();
  q.close();
  consumer.join();
  EXPECT_TRUE(ordered);
  EXPECT_EQ(total, static_cast<long>(kProducers) * kItems);
  EXPECT_EQ(q.dropped(), 0u);
  for (int p = 0; p < kProducers; ++p) EXPECT_EQ(next[static_cast<std::size_t>(p)], kItems);
}
