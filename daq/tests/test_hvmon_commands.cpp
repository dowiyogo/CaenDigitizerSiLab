// Tests del protocolo de comandos de hvmon: cada comando por un socket real, peticiones mal
// formadas y clientes que cortan a mitad de la petición. El hilo de sondeo se reemplaza por un
// "respondedor" que contesta las peticiones largas.
#include <gtest/gtest.h>

#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

#include <cstring>

#include "daq/core/clock.h"
#include "daq/core/jthread.h"
#include "daq/core/unix_socket.h"
#include "daq/hv/fake_hv_backend.h"
#include "daq/hv/reading_json.h"
#include "daq/hvmon/commands.h"

using namespace std::chrono_literals;
using namespace daq;
using namespace daq::hvmon;

namespace {

struct CRig {
  SystemClock clock;
  HvmonConfig cfg = HvmonConfig::fromJson(nlohmann::json::object()).value();
  StateHub hub;
  BoundedQueue<PollRequest> poll{16};
  BoundedQueue<WriterItem> writer{1000};
  CommandHandler handler{cfg, hub, poll, writer, clock};
  // Un socket por test: ctest corre los tests en paralelo.
  std::string sock = ::testing::TempDir() + "daq_hvmon_" +
                     ::testing::UnitTest::GetInstance()->current_test_info()->name() + ".sock";
  std::unique_ptr<CommandServer> server;
  JThread serverThread, responder;

  CRig() {
    // Estado con una lectura completa del fake.
    hv::FakeHVBackend be(hv::FakeHVBackend::defaultScenario(), clock);
    EXPECT_TRUE(be.connect().ok());
    auto names = be.channelParamNames(3, 0).value();
    std::vector<std::string> pn;
    for (const auto& p : names) pn.push_back(p.name);
    auto rd = be.readChannels(3, cfg.channels, pn).value();
    PollerState ps;
    ps.connected = true;
    ps.params = rd.params;
    ps.lastGoodReadUtcNs = rd.tRespUtcNs;
    for (std::size_t p = 0; p < rd.params.size(); ++p)
      for (std::size_t c = 0; c < rd.channels.size(); ++c) ps.values[rd.channels[c]][rd.params[p].name] = rd.values[p][c];
    for (int i = 0; i < 30; ++i) {
      ChangeRecord r;
      r.id = static_cast<std::uint64_t>(i + 1);
      r.ch = i % 12;
      r.param = "V0Set";
      r.oldValue = 1.0;
      r.newValue = 2.0;
      r.origin = "external";
      ps.recentChanges.push_back(r);
    }
    hub.poller.publish(ps);
    auto l = UnixListener::listen(sock);
    EXPECT_TRUE(l.ok());
    server = std::make_unique<CommandServer>(std::move(l).value(), handler);
    serverThread = JThread([this](StopToken st) { server->run(st); });
    responder = JThread([this](StopToken st) {
      while (!st.stop_requested()) {
        auto r = poll.pop_until(std::chrono::steady_clock::now() + 20ms);
        if (!r) continue;
        if (r->kind == PollRequest::Kind::Pw)
          r->reply->set_value({{"ok", false}, {"result", "rejected"}, {"detail", "allow_commands = false"},
                               {"ch", r->ch}, {"echo_reason", r->reason}, {"echo_requester", r->requester}});
        else
          r->reply->set_value({{"ok", true}, {"offsets", {{"seconds", r->seconds}}}});
      }
    });
  }
  ~CRig() {
    serverThread.request_stop();
    serverThread.join();
    responder.request_stop();
    responder.join();
  }
  nlohmann::json ask(const nlohmann::json& req) { return askRaw(req.dump()); }
  nlohmann::json askRaw(const std::string& line) {
    auto r = unixRequest(sock, line, 5s);
    EXPECT_TRUE(r.ok()) << (r.ok() ? "" : r.error().message);
    return r.ok() ? nlohmann::json::parse(r.value()) : nlohmann::json();
  }
};

}  // namespace

TEST(HvmonCommands, PingStatusRead) {
  CRig r;
  EXPECT_TRUE(r.ask({{"cmd", "ping"}})["ok"].get<bool>());
  const auto st = r.ask({{"cmd", "status"}});
  ASSERT_TRUE(st["ok"].get<bool>());
  EXPECT_EQ(st["state"]["schema"], "daqv2.hvmon_state/1");
  EXPECT_EQ(st["state"]["channels"].size(), 12u);
  EXPECT_FALSE(st["state"]["allow_commands"].get<bool>());
  const auto rd = r.ask({{"cmd", "read"}});
  ASSERT_TRUE(rd["ok"].get<bool>());
  auto parsed = hv::readingFromJson(rd["reading"]);
  ASSERT_TRUE(parsed.ok()) << parsed.error().message;
  EXPECT_EQ(parsed.value().params.size(), 18u);
  EXPECT_DOUBLE_EQ(parsed.value().get("V0Set", 6)->boardValue, 66.15);
}

TEST(HvmonCommands, ChangesReturnsLastN) {
  CRig r;
  const auto c = r.ask({{"cmd", "changes"}, {"n", 5}});
  ASSERT_TRUE(c["ok"].get<bool>());
  ASSERT_EQ(c["changes"].size(), 5u);
  EXPECT_EQ(c["changes"][4]["id"], 30);
  EXPECT_EQ(r.ask({{"cmd", "changes"}, {"n", 100}})["changes"].size(), 30u);
  EXPECT_FALSE(r.ask({{"cmd", "changes"}, {"n", 0}})["ok"].get<bool>());
}

TEST(HvmonCommands, NoteGoesToWriter) {
  CRig r;
  const auto n = r.ask({{"cmd", "note"}, {"text", "GECO abierto por René"}});
  ASSERT_TRUE(n["ok"].get<bool>());
  auto item = r.writer.try_pop();
  ASSERT_TRUE(item.has_value());
  const auto* note = std::get_if<NoteItem>(&*item);
  ASSERT_NE(note, nullptr);
  EXPECT_EQ(note->text, "GECO abierto por René");
  EXPECT_NE(note->user.find("uid"), std::string::npos);
  EXPECT_FALSE(r.ask({{"cmd", "note"}, {"text", ""}})["ok"].get<bool>());
  EXPECT_FALSE(r.ask({{"cmd", "note"}})["ok"].get<bool>());
}

TEST(HvmonCommands, LongCommandsAreForwarded) {
  CRig r;
  const auto pw = r.ask({{"cmd", "pw"}, {"ch", 3}, {"on", true}, {"authorize", "prueba de protocolo"}});
  EXPECT_FALSE(pw["ok"].get<bool>());
  EXPECT_EQ(pw["result"], "rejected");
  EXPECT_EQ(pw["echo_reason"], "prueba de protocolo");
  EXPECT_NE(pw["echo_requester"].get<std::string>().find("pid"), std::string::npos);
  EXPECT_FALSE(r.ask({{"cmd", "pw"}, {"ch", 3}, {"on", true}})["ok"].get<bool>());  // falta authorize
  EXPECT_FALSE(r.ask({{"cmd", "pw"}, {"ch", "3"}, {"on", true}, {"authorize", "x"}})["ok"].get<bool>());
  const auto off = r.ask({{"cmd", "measure-offsets"}, {"seconds", 2}});
  EXPECT_TRUE(off["ok"].get<bool>());
  EXPECT_EQ(off["offsets"]["seconds"], 2);
  EXPECT_FALSE(r.ask({{"cmd", "measure-offsets"}})["ok"].get<bool>());
}

TEST(HvmonCommands, MalformedRequests) {
  CRig r;
  EXPECT_FALSE(r.askRaw("esto no es json")["ok"].get<bool>());
  EXPECT_FALSE(r.askRaw("[1,2,3]")["ok"].get<bool>());
  EXPECT_FALSE(r.askRaw("{\"sin\":\"cmd\"}")["ok"].get<bool>());
  EXPECT_FALSE(r.askRaw("{\"cmd\":42}")["ok"].get<bool>());
  const auto u = r.ask({{"cmd", "set-v0"}});
  EXPECT_FALSE(u["ok"].get<bool>());
  EXPECT_NE(u["error"].get<std::string>().find("desconocido"), std::string::npos);
}

TEST(HvmonCommands, ClientDisconnectingMidRequest) {
  CRig r;
  for (int i = 0; i < 3; ++i) {
    const int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
    sockaddr_un a{};
    a.sun_family = AF_UNIX;
    std::strncpy(a.sun_path, r.sock.c_str(), sizeof(a.sun_path) - 1);
    ASSERT_EQ(::connect(fd, reinterpret_cast<sockaddr*>(&a), sizeof a), 0);
    ::send(fd, "{\"cmd\":\"stat", 12, MSG_NOSIGNAL);
    ::close(fd);
  }
  // El servidor sigue atendiendo.
  EXPECT_TRUE(r.ask({{"cmd", "ping"}})["ok"].get<bool>());
  // Un cliente que pide y se va sin leer la respuesta tampoco lo afecta.
  {
    const int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
    sockaddr_un a{};
    a.sun_family = AF_UNIX;
    std::strncpy(a.sun_path, r.sock.c_str(), sizeof(a.sun_path) - 1);
    ASSERT_EQ(::connect(fd, reinterpret_cast<sockaddr*>(&a), sizeof a), 0);
    ::send(fd, "{\"cmd\":\"status\"}\n", 17, MSG_NOSIGNAL);
    ::close(fd);
  }
  EXPECT_TRUE(r.ask({{"cmd", "status"}})["ok"].get<bool>());
}
