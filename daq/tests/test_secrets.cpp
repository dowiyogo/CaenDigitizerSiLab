// Tests de secretos: archivo KEY=VALUE con modo 600, prioridad del entorno, redacción.
// Nunca leen el secrets.env real: DAQV2_SECRETS apunta a un temporal y las variables DAQ_* se
// borran del entorno del test. Los valores son ficticios.
#include <gtest/gtest.h>

#include <sys/stat.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>

#include "daq/core/logger.h"
#include "daq/core/secrets.h"
#include "daq/hv/credentials.h"

namespace {

const char* kVars[] = {"DAQ_HV_HOST", "DAQ_HV_USER", "DAQ_HV_PASSWORD", "DAQ_SMTP_URL", "DAQ_SMTP_USER",
                       "DAQ_SMTP_PASSWORD", "DAQ_NOTIFY_FROM", "DAQ_NOTIFY_TO"};

class SecretsTest : public ::testing::Test {
 protected:
  void SetUp() override {
    for (const char* v : kVars) ::unsetenv(v);
    path_ = ::testing::TempDir() + "daq_test_secrets.env";
    ::setenv("DAQV2_SECRETS", path_.c_str(), 1);
  }
  void TearDown() override {
    std::remove(path_.c_str());
    ::unsetenv("DAQV2_SECRETS");
    for (const char* v : kVars) ::unsetenv(v);
  }
  void write(const std::string& content, mode_t mode = 0600) {
    std::remove(path_.c_str());
    {
      std::ofstream o(path_);
      o << content;
    }
    ::chmod(path_.c_str(), mode);
  }
  std::string path_;
};

}  // namespace

TEST_F(SecretsTest, DefaultPathFollowsEnvironment) { EXPECT_EQ(daq::Secrets::defaultPath(), path_); }

TEST_F(SecretsTest, ParsesKeyValueCommentsAndLiteralValues) {
  write("# comentario\n\nDAQ_HV_USER=operador\n  DAQ_HV_PASSWORD=p@ss=word con espacios\r\nDAQ_SMTP_URL=smtps://h:465\n"
        "DAQ_NOTIFY_TO=$HOME\n");
  daq::Secrets s;
  EXPECT_EQ(s.get("DAQ_HV_USER").value().value(), "operador");
  EXPECT_EQ(s.get("DAQ_HV_PASSWORD").value().value(), "p@ss=word con espacios");
  EXPECT_EQ(s.get("DAQ_SMTP_URL").value().value(), "smtps://h:465");
  EXPECT_EQ(s.get("DAQ_NOTIFY_TO").value().value(), "$HOME");  // sin expansión de shell
  EXPECT_FALSE(s.get("DAQ_SMTP_USER").value().has_value());
}

TEST_F(SecretsTest, EnvironmentWinsOverFile) {
  write("DAQ_HV_USER=desde_archivo\n");
  ::setenv("DAQ_HV_USER", "desde_entorno", 1);
  daq::Secrets s;
  EXPECT_EQ(s.get("DAQ_HV_USER").value().value(), "desde_entorno");
}

TEST_F(SecretsTest, MissingFileMeansNoValues) {
  std::remove(path_.c_str());
  daq::Secrets s;
  auto v = s.get("DAQ_HV_USER");
  ASSERT_TRUE(v.ok());
  EXPECT_FALSE(v.value().has_value());
}

TEST_F(SecretsTest, RefusesFileWithWrongMode) {
  for (mode_t m : {0644, 0640, 0660, 0400, 0700}) {
    write("DAQ_HV_USER=x\n", m);
    daq::Secrets s;
    EXPECT_FALSE(s.get("DAQ_HV_USER").ok()) << std::oct << m;
  }
}

TEST_F(SecretsTest, RefusesMalformedLines) {
  write("DAQ_HV_USER\n");
  EXPECT_FALSE(daq::parseEnvFile(path_).ok());
  write("MAL CLAVE=x\n");
  EXPECT_FALSE(daq::parseEnvFile(path_).ok());
}

TEST_F(SecretsTest, RedactsEveryPasswordKey) {
  const nlohmann::json cfg{{"DAQ_HV_PASSWORD", "secreto1"},
                           {"smtp", {{"DAQ_SMTP_PASSWORD", "secreto2"}, {"DAQ_SMTP_USER", "u"}}},
                           {"list", nlohmann::json::array({{{"password", "secreto3"}}})},
                           {"host", "192.168.1.118"}};
  const auto red = daq::redactSecrets(cfg);
  const std::string dump = red.dump();
  EXPECT_EQ(dump.find("secreto"), std::string::npos) << dump;
  EXPECT_EQ(red["DAQ_HV_PASSWORD"], "***");
  EXPECT_EQ(red["smtp"]["DAQ_SMTP_PASSWORD"], "***");
  EXPECT_EQ(red["smtp"]["DAQ_SMTP_USER"], "u");
  EXPECT_EQ(red["list"][0]["password"], "***");
  EXPECT_EQ(red["host"], "192.168.1.118");
}

TEST_F(SecretsTest, LoggedConfigDumpContainsNoPassword) {
  const std::string log = ::testing::TempDir() + "daq_test_redact.log";
  std::remove(log.c_str());
  auto& L = daq::Logger::instance();
  L.setStderr(false);
  ASSERT_TRUE(L.openFile(log).ok());
  const nlohmann::json cfg{{"DAQ_HV_USER", "admin"}, {"DAQ_HV_PASSWORD", "NoDebeAparecer1"},
                           {"DAQ_SMTP_PASSWORD", "NoDebeAparecer2"}};
  daq::logInfo("configuración: ", daq::redactSecrets(cfg).dump());
  L.closeFile();
  L.setStderr(true);
  std::ifstream in(log);
  const std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_EQ(all.find("NoDebeAparecer"), std::string::npos);
  EXPECT_NE(all.find("\"DAQ_HV_PASSWORD\":\"***\""), std::string::npos);
  EXPECT_NE(all.find("admin"), std::string::npos);
  std::remove(log.c_str());
}

TEST_F(SecretsTest, HvLoginPrefersSecretsAndFallsBackToLegacyJson) {
  write("DAQ_HV_HOST=10.0.0.5\nDAQ_HV_USER=u1\nDAQ_HV_PASSWORD=p1\n");
  {
    daq::Secrets s;
    auto l = daq::hv::resolveHvLogin(s, "192.168.1.118", "/no/existe.json");
    ASSERT_TRUE(l.ok()) << l.error().message;
    EXPECT_EQ(l.value().host, "10.0.0.5");
    EXPECT_EQ(l.value().credentials.user, "u1");
    EXPECT_EQ(l.value().source, "secrets");
  }
  write("# sin credenciales HV\n");
  const std::string legacy = ::testing::TempDir() + "daq_test_legacy.json";
  {
    std::ofstream o(legacy);
    o << R"({"user":"u2","password":"p2"})";
  }
  ::chmod(legacy.c_str(), 0600);
  {
    daq::Secrets s;
    auto l = daq::hv::resolveHvLogin(s, "192.168.1.118", legacy);
    ASSERT_TRUE(l.ok()) << l.error().message;
    EXPECT_EQ(l.value().host, "192.168.1.118");
    EXPECT_EQ(l.value().credentials.user, "u2");
    EXPECT_EQ(l.value().source, "hv_credentials.json (obsoleto)");
  }
  {
    daq::Secrets s;
    auto l = daq::hv::resolveHvLogin(s, "192.168.1.118", "/no/existe.json");
    ASSERT_FALSE(l.ok());
    EXPECT_EQ(l.error().message.find("p2"), std::string::npos);
  }
  std::remove(legacy.c_str());
}
