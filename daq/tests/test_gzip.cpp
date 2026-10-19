// Tests de compresión verificada, lectura de .gz y selección de días cerrados.
#include <gtest/gtest.h>

#include <sys/stat.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>

#include "daq/core/gzip.h"
#include "daq/hvmon/compressor.h"

namespace {
std::string dir() {
  const std::string d = ::testing::TempDir() + "daq_gzip_" + ::testing::UnitTest::GetInstance()->current_test_info()->name();
  ::mkdir(d.c_str(), 0700);
  return d;
}
bool exists(const std::string& p) {
  struct stat st {};
  return ::stat(p.c_str(), &st) == 0;
}
void writeLines(const std::string& p, int n) {
  std::ofstream o(p);
  for (int i = 0; i < n; ++i) o << "1792281600000000000," << i << ",3," << i % 12 << ",65.75,1.5e-05,0.18,-4.8e-08,0,0\n";
}
}  // namespace

TEST(Gzip, RoundTripVerifiesAndRemovesPlainFile) {
  const std::string d = dir();
  const std::string f = d + "/hv_20261018.csv";
  writeLines(f, 5000);
  auto before = daq::digestFile(f);
  ASSERT_TRUE(before.ok());
  auto r = daq::gzipFileVerified(f, f + ".gz", true);
  ASSERT_TRUE(r.ok()) << r.error().message;
  EXPECT_FALSE(exists(f));
  EXPECT_TRUE(exists(f + ".gz"));
  EXPECT_EQ(r.value().source.lines, 5000u);
  EXPECT_LT(r.value().ratio(), 0.3);
  auto after = daq::digestFile(f + ".gz");  // lee el contenido descomprimido
  ASSERT_TRUE(after.ok());
  EXPECT_EQ(after.value().bytes, before.value().bytes);
  EXPECT_EQ(after.value().crc32, before.value().crc32);
  // LineReader lee el .gz y también archivos planos.
  auto lr = daq::LineReader::open(f + ".gz");
  ASSERT_TRUE(lr.ok());
  int n = 0;
  std::string first;
  while (auto l = lr.value()->next()) {
    if (n == 0) first = *l;
    ++n;
  }
  EXPECT_EQ(n, 5000);
  EXPECT_EQ(first.rfind("1792281600000000000,0,", 0), 0u);
  writeLines(d + "/plano.csv", 3);
  auto lp = daq::LineReader::open(d + "/plano.csv");
  ASSERT_TRUE(lp.ok());
  int m = 0;
  while (lp.value()->next()) ++m;
  EXPECT_EQ(m, 3);
}

TEST(Gzip, VerificationFailureKeepsPlainFile) {
  const std::string d = dir();
  const std::string f = d + "/hv_changes_20261018.jsonl";
  writeLines(f, 200);
  auto r = daq::gzipFileVerified(f, f + ".gz", true, 6, [](const std::string& gz) {
    // Corrompe el .gz antes de la verificación (trunca).
    ASSERT_EQ(::truncate(gz.c_str(), 40), 0);
  });
  ASSERT_FALSE(r.ok());
  EXPECT_NE(r.error().message.find("se conserva"), std::string::npos);
  EXPECT_TRUE(exists(f));                // el plano NO se borró
  EXPECT_FALSE(exists(f + ".gz"));
  EXPECT_TRUE(exists(f + ".gz.failed"));  // el .gz defectuoso queda apartado
  auto d2 = daq::digestFile(f);
  ASSERT_TRUE(d2.ok());
  EXPECT_EQ(d2.value().lines, 200u);
}

TEST(Gzip, MissingSourceFails) { EXPECT_FALSE(daq::gzipFileVerified("/no/existe.csv", "/no/existe.csv.gz", true).ok()); }

TEST(Compressor, SelectsOnlyClosedDailyFilesAndKeepsSchemas) {
  const std::string base = dir();
  ::mkdir((base + "/2026").c_str(), 0700);
  for (const char* n : {"hv_20261018.csv", "hv_full_20261018.csv", "notes_20261018.jsonl", "hv_20261019.csv",
                        "hv_20261018.csv.schema.json", "hv_20261017.csv.gz", "offsets.jsonl"}) {
    std::ofstream o(base + "/2026/" + n);
    o << "x\n";
  }
  const auto files = daq::hvmon::Compressor::closedDailyFiles(base, "20261019");
  EXPECT_EQ(files, (std::vector<std::string>{base + "/2026/hv_20261018.csv", base + "/2026/hv_full_20261018.csv",
                                             base + "/2026/notes_20261018.jsonl"}));
  daq::hvmon::StateHub hub;
  daq::hvmon::Compressor c(hub);
  for (const auto& f : files) ASSERT_TRUE(c.compressOne(f).ok());
  EXPECT_TRUE(exists(base + "/2026/hv_20261018.csv.gz"));
  EXPECT_TRUE(exists(base + "/2026/hv_20261018.csv.schema.json"));  // el esquema queda plano
  EXPECT_TRUE(exists(base + "/2026/hv_20261019.csv"));              // el día de hoy no se toca
  EXPECT_EQ(hub.compressor.read().first["compressed"], 3);
}
