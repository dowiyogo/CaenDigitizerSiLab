// commands.cpp
#include "daq/hvmon/commands.h"

#include <pwd.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>

#include "daq/core/logger.h"
#include "daq/hv/reading_json.h"

namespace daq::hvmon {

namespace {
constexpr std::size_t kMaxRequestBytes = 64 * 1024;
constexpr std::size_t kMaxNoteChars = 4000;

nlohmann::json error(const std::string& msg) { return nlohmann::json{{"ok", false}, {"error", msg}}; }

std::string hostName() {
  char buf[256] = {0};
  if (::gethostname(buf, sizeof buf - 1) != 0) return "?";
  return buf;
}
}  // namespace

std::string describePeer(const PeerCredentials& peer) {
  std::string user = std::to_string(peer.uid);
  if (passwd* pw = ::getpwuid(peer.uid)) user = pw->pw_name;
  return user + " (uid " + std::to_string(peer.uid) + ", pid " + std::to_string(peer.pid) + ")";
}

CommandHandler::CommandHandler(const HvmonConfig& cfg, StateHub& hub, BoundedQueue<PollRequest>& poll,
                               BoundedQueue<WriterItem>& writer, const Clock& clock)
    : cfg_(cfg), hub_(hub), poll_(poll), writer_(writer), clock_(clock) {}

bool CommandHandler::isLong(const std::string& line) {
  try {
    const auto j = nlohmann::json::parse(line);
    const std::string cmd = j.value("cmd", "");
    return cmd == "pw" || cmd == "measure-offsets" || cmd == "notify-test" || cmd == "bias-begin" || cmd == "v0set" ||
           cmd == "bias-restore" || cmd == "notify-summary";
  } catch (...) {
    return false;
  }
}

nlohmann::json CommandHandler::handle(const std::string& line, const PeerCredentials& peer) {
  try {
    const auto req = nlohmann::json::parse(line);
    if (!req.is_object() || !req.contains("cmd") || !req["cmd"].is_string())
      return error("petición inválida: se espera un objeto JSON con \"cmd\"");
    return dispatch(req, peer);
  } catch (const nlohmann::json::exception& e) {
    return error(std::string("petición inválida: ") + e.what());
  } catch (const std::exception& e) {
    return error(std::string("error interno: ") + e.what());
  }
}

nlohmann::json CommandHandler::forward(PollRequest r, double timeoutS) {
  auto fut = r.reply->get_future();
  if (!poll_.try_push(std::move(r), DropPolicy::DropNewest)) return error("cola de peticiones llena o cerrada");
  if (fut.wait_for(std::chrono::duration<double>(timeoutS)) != std::future_status::ready)
    return error("sin respuesta del hilo de sondeo a tiempo");
  return fut.get();
}

nlohmann::json CommandHandler::readJson() const {
  const auto [p, seq] = hub_.poller.read();
  (void)seq;
  if (p.params.empty() || p.values.empty()) return error("todavía no hay lecturas");
  hv::ChannelReading r;
  r.slot = cfg_.slot;
  r.channels = cfg_.channels;
  r.tReqUtcNs = p.lastGoodReadUtcNs;
  r.tRespUtcNs = p.lastGoodReadUtcNs;
  for (const auto& info : p.params) {
    std::vector<hv::ParamValue> row;
    for (int ch : cfg_.channels) {
      auto c = p.values.find(ch);
      if (c == p.values.end()) return error("falta el canal " + std::to_string(ch) + " en el estado");
      auto v = c->second.find(info.name);
      if (v == c->second.end()) return error("falta " + info.name + " en el estado");
      row.push_back(v->second);
    }
    r.params.push_back(info);
    r.values.push_back(std::move(row));
  }
  auto j = hv::readingToJson(r, cfg_.statusBits);
  j["source"] = "hvmon";
  j["connected"] = p.connected;
  j["last_good_read_utc_ns"] = p.lastGoodReadUtcNs;
  return nlohmann::json{{"ok", true}, {"reading", j}};
}

nlohmann::json CommandHandler::dispatch(const nlohmann::json& req, const PeerCredentials& peer) {
  const std::string cmd = req["cmd"].get<std::string>();
  if (cmd == "ping") return {{"ok", true}, {"pid", static_cast<long>(::getpid())}};
  if (cmd == "status")
    return {{"ok", true}, {"state", hub_.compose(cfg_, clock_.now_utc_ns(), ++statusSeq_, false)}};
  if (cmd == "read") return readJson();
  if (cmd == "changes") {
    const int n = req.value("n", 20);
    if (n <= 0) return error("n debe ser positivo");
    const auto [p, seq] = hub_.poller.read();
    (void)seq;
    nlohmann::json out = nlohmann::json::array();
    const std::size_t total = p.recentChanges.size();
    const std::size_t start = total > static_cast<std::size_t>(n) ? total - static_cast<std::size_t>(n) : 0;
    for (std::size_t i = start; i < total; ++i) out.push_back(p.recentChanges[i].toJson());
    return {{"ok", true}, {"changes", out}, {"ring_size", cfg_.changeRingSize}};
  }
  if (cmd == "note") {
    if (!req.contains("text") || !req["text"].is_string()) return error("note: falta \"text\"");
    const std::string text = req["text"].get<std::string>();
    if (text.empty()) return error("note: texto vacío");
    if (text.size() > kMaxNoteChars) return error("note: texto demasiado largo");
    NoteItem n;
    n.tUtcNs = clock_.now_utc_ns();
    n.user = describePeer(peer);
    n.host = hostName();
    n.text = text;
    if (!writer_.try_push(WriterItem{n}, DropPolicy::DropNewest)) return error("note: cola del escritor llena");
    logInfo("nota de ", n.user, ": ", text);
    return {{"ok", true}, {"t_utc_ns", n.tUtcNs}, {"t_utc", formatUtcIso(n.tUtcNs)}};
  }
  if (cmd == "measure-offsets") {
    PollRequest r;
    r.kind = PollRequest::Kind::Offsets;
    if (!req.contains("seconds") || !req["seconds"].is_number()) return error("measure-offsets: falta \"seconds\"");
    r.seconds = req["seconds"].get<double>();
    if (req.contains("channels")) r.channels = req["channels"].get<std::vector<int>>();
    r.requester = describePeer(peer);
    const double timeoutS = r.seconds + 60.0;
    return forward(std::move(r), timeoutS);
  }
  if (cmd == "pw") {
    if (!req.contains("ch") || !req["ch"].is_number_integer()) return error("pw: falta \"ch\"");
    if (!req.contains("on") || !req["on"].is_boolean()) return error("pw: falta \"on\" (true/false)");
    const std::string reason = req.value("authorize", "");
    if (reason.empty()) return error("pw: falta \"authorize\" (motivo de la escritura)");
    PollRequest r;
    r.kind = PollRequest::Kind::Pw;
    r.ch = req["ch"].get<int>();
    r.on = req["on"].get<bool>();
    r.reason = reason;
    r.requester = describePeer(peer);
    // El rechazo con allow_commands = false lo hace (y lo registra) el hilo de sondeo.
    return forward(std::move(r), 3600.0);
  }
  if (cmd == "bias-begin" || cmd == "v0set" || cmd == "bias-restore" || cmd == "bias-status") {
    PollRequest r;
    r.requester = describePeer(peer);
    r.reason = req.value("authorize", "");
    if (cmd != "bias-status" && r.reason.empty()) return error(cmd + ": falta \"authorize\" (motivo de la escritura)");
    if (cmd == "bias-begin") {
      r.kind = PollRequest::Kind::BiasBegin;
      if (!req.contains("expected") || !req["expected"].is_object()) return error("bias-begin: falta \"expected\" {canal: V0Set}");
      for (const auto& [k, v] : req["expected"].items()) r.expected[std::stoi(k)] = v.get<double>();
      // DAQ2-17: límite sobre el nominal y escalón de la sesión (los topes los aplica BiasSession::begin)
      if (req.contains("max_above_V")) {
        if (!req["max_above_V"].is_number()) return error("bias-begin: \"max_above_V\" debe ser un número (V)");
        r.maxAboveV = req["max_above_V"].get<double>();
      }
      if (req.contains("max_step_V")) {
        if (!req["max_step_V"].is_number()) return error("bias-begin: \"max_step_V\" debe ser un número (V)");
        r.maxStepV = req["max_step_V"].get<double>();
      }
    } else if (cmd == "v0set") {
      r.kind = PollRequest::Kind::V0Set;
      if (!req.contains("ch") || !req["ch"].is_number_integer()) return error("v0set: falta \"ch\"");
      if (!req.contains("value") || !req["value"].is_number()) return error("v0set: falta \"value\" (V)");
      r.ch = req["ch"].get<int>();
      r.value = req["value"].get<double>();
    } else if (cmd == "bias-restore") {
      r.kind = PollRequest::Kind::BiasRestore;
    } else {
      r.kind = PollRequest::Kind::BiasStatus;
    }
    return forward(std::move(r), 3600.0);
  }
  if (cmd == "notify-summary") {
    // DAQ2-14: correo de resumen (p. ej. del barrido de tensión) con el texto pedido; cuenta en los límites del notificador.
    if (!notifyQ_) return error("notificador no disponible en este hvmon");
    if (!req.contains("text") || !req["text"].is_string()) return error("notify-summary: falta \"text\"");
    NotifyEvent e;
    e.kind = req.value("kind", std::string("summary"));
    e.ch = -1;
    e.slot = cfg_.slot;
    e.reason = req["text"].get<std::string>();
    e.origin = req.value("origin", describePeer(peer));
    e.tUtcNs = clock_.now_utc_ns();
    // Contexto como en notify-test (sin él, el cuerpo decía "Temperatura de la placa: sin dato").
    const auto state = hub_.compose(cfg_, e.tUtcNs, ++statusSeq_, false);
    e.context["board_temperature_degC"] = state["board_temperature_degC"];
    e.context["geco_running"] = state["geco_running"];
    e.reply = std::make_shared<std::promise<nlohmann::json>>();
    auto fut = e.reply->get_future();
    if (!notifyQ_->try_push(std::move(e), DropPolicy::DropNewest)) return error("cola del notificador llena");
    if (fut.wait_for(std::chrono::seconds(180)) != std::future_status::ready) return error("sin respuesta del notificador en 180 s");
    return fut.get();
  }
  if (cmd == "notify-test") {
    if (!notifyQ_) return error("notificador no disponible en este hvmon");
    NotifyEvent e;
    e.kind = "gave_up";
    e.test = true;
    e.ch = req.value("ch", 11);
    e.slot = cfg_.slot;
    e.k = 3;
    e.reason = "PRUEBA de hvctl notify-test pedida por " + describePeer(peer) + " (evento sintético)";
    const UtcNs now = clock_.now_utc_ns();
    e.tUtcNs = now;
    for (int i = 0; i < 4; ++i)
      e.timeline.push_back(
          {{"t_utc_ns", now}, {"t_utc", formatUtcIso(now)}, {"event", i < 3 ? "trip (PRUEBA)" : "GAVE_UP (PRUEBA)"}, {"k", i}});
    const auto state = hub_.compose(cfg_, now, ++statusSeq_, false);
    for (const auto& c : state["channels"])
      if (c.value("ch", -1) == e.ch) {
        const auto& v = c["values"];
        e.context = {{"VMon_V", v.value("VMon_V", 0.0)}, {"IMon_A", v.value("IMon_A", 0.0)}, {"Status", v.value("Status", 0u)},
                     {"status_flags", c["status_flags"]}};
      }
    e.context["board_temperature_degC"] = state["board_temperature_degC"];
    e.context["geco_running"] = state["geco_running"];
    e.reply = std::make_shared<std::promise<nlohmann::json>>();
    auto fut = e.reply->get_future();
    if (!notifyQ_->try_push(std::move(e), DropPolicy::DropNewest)) return error("cola del notificador llena");
    if (fut.wait_for(std::chrono::seconds(180)) != std::future_status::ready)
      return error("sin respuesta del notificador en 180 s");
    return fut.get();
  }
  return error("comando desconocido: " + cmd);
}

CommandServer::CommandServer(UnixListener listener, CommandHandler& handler)
    : listener_(std::move(listener)), handler_(handler) {}

void CommandServer::run(const StopToken& st) {
  Logger::setThreadName("cmd");
  while (!st.stop_requested()) {
    workers_.remove_if([](const Worker& w) { return w.done->load(); });  // JThread se une al destruirse
    auto c = listener_.acceptFor(std::chrono::milliseconds(200));
    if (!c) {
      logError("socket: ", c.error().message);
      continue;
    }
    FdGuard fd = std::move(c).value();
    if (fd.get() < 0) continue;
    // Las peticiones largas se atienden en su propio hilo; hay que leer la línea primero.
    auto line = readLine(fd.get(), kMaxRequestBytes, std::chrono::seconds(5));
    if (!line) {
      logDebug("socket: petición incompleta descartada: ", line.error().message);
      continue;
    }
    const auto peer = peerCredentials(fd.get());
    if (CommandHandler::isLong(line.value())) {
      auto done = std::make_shared<std::atomic<bool>>(false);
      auto shared = std::make_shared<FdGuard>(std::move(fd));
      std::string req = line.value();
      workers_.push_back(Worker{done, JThread([this, done, shared, req, peer] {
                                  Logger::setThreadName("cmd-long");
                                  const auto resp = handler_.handle(req, peer);
                                  writeAll(shared->get(), resp.dump() + "\n");
                                  done->store(true);
                                })});
    } else {
      const auto resp = handler_.handle(line.value(), peer);
      if (auto w = writeAll(fd.get(), resp.dump() + "\n"); !w) logDebug("socket: el cliente cerró: ", w.error().message);
    }
  }
  workers_.clear();  // une los hilos largos (sus futuros ya recibieron respuesta al apagar el sondeo)
}

}  // namespace daq::hvmon
