// Tests de carga de credenciales: modo 600 obligatorio.
#include <gtest/gtest.h>

#include <sys/stat.h>

#include <cstdio>
#include <fstream>

#include "daq/hv/credentials.h"

namespace {
std::string writeFile(const std::string& name, const std::string& content, mode_t mode) {
  const std::string p = ::testing::TempDir() + name;
  std::remove(p.c_str());
  {
    std::ofstream o(p);
    o << content;
  }
  ::chmod(p.c_str(), mode);
  return p;
}
}  // namespace

TEST(Credentials, LoadsWithMode600) {
  const auto p = writeFile("cred_ok.json", R"({"user":"u","password":"p"})", 0600);
  auto c = daq::hv::loadCredentials(p);
  ASSERT_TRUE(c.ok()) << c.error().message;
  EXPECT_EQ(c.value().user, "u");
  EXPECT_EQ(c.value().password, "p");
  std::remove(p.c_str());
}

TEST(Credentials, RefusesOtherModes) {
  for (mode_t m : {0644, 0640, 0400, 0700}) {
    const auto p = writeFile("cred_mode.json", R"({"user":"u","password":"p"})", m);
    auto c = daq::hv::loadCredentials(p);
    EXPECT_FALSE(c.ok()) << "modo " << std::oct << m;
    std::remove(p.c_str());
  }
}

TEST(Credentials, RefusesMissingFieldsAndMissingFile) {
  const auto p = writeFile("cred_bad.json", R"({"user":"u"})", 0600);
  EXPECT_FALSE(daq::hv::loadCredentials(p).ok());
  std::remove(p.c_str());
  EXPECT_FALSE(daq::hv::loadCredentials(::testing::TempDir() + "no_existe.json").ok());
}

TEST(Credentials, ExpandsHome) {
  EXPECT_EQ(daq::hv::expandHome("/abs"), "/abs");
  const char* home = std::getenv("HOME");
  ASSERT_NE(home, nullptr);
  EXPECT_EQ(daq::hv::expandHome("~/x"), std::string(home) + "/x");
}
