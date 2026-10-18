// Tests de Snapshot: secuencia creciente y lecturas consistentes con lectores concurrentes.
#include <gtest/gtest.h>

#include <atomic>
#include <thread>
#include <vector>

#include "daq/core/snapshot.h"

namespace {
// Valor con redundancia: si un lector viera un valor a medio escribir, a != b o c != a + b.
struct Triple {
  long a = 0, b = 0, c = 0;
};
}  // namespace

TEST(Snapshot, SequenceStartsAtZeroAndIncreases) {
  daq::Snapshot<int> s;
  EXPECT_EQ(s.read().second, 0u);
  EXPECT_EQ(s.publish(5), 1u);
  EXPECT_EQ(s.publish(6), 2u);
  auto [v, seq] = s.read();
  EXPECT_EQ(v, 6);
  EXPECT_EQ(seq, 2u);
  EXPECT_EQ(s.sequence(), 2u);
}

TEST(Snapshot, ConcurrentReadersSeeConsistentValuesAndMonotonicSequence) {
  daq::Snapshot<Triple> s;
  constexpr long kWrites = 200000;
  std::atomic<bool> done{false};
  std::atomic<long> inconsistent{0}, backwards{0};
  std::vector<std::thread> readers;
  for (int r = 0; r < 3; ++r)
    readers.emplace_back([&] {
      std::uint64_t last = 0;
      while (!done.load()) {
        auto [v, seq] = s.read();
        if (v.a != v.b || v.c != v.a + v.b || static_cast<std::uint64_t>(v.a) != seq) ++inconsistent;
        if (seq < last) ++backwards;
        last = seq;
      }
    });
  for (long i = 1; i <= kWrites; ++i) s.publish(Triple{i, i, 2 * i});
  done = true;
  for (auto& t : readers) t.join();
  EXPECT_EQ(inconsistent.load(), 0);
  EXPECT_EQ(backwards.load(), 0);
  EXPECT_EQ(s.sequence(), static_cast<std::uint64_t>(kWrites));
}
