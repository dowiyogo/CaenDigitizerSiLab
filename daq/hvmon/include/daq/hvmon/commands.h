// commands.h
// Protocolo de comandos de hvmon por socket Unix ($XDG_RUNTIME_DIR/daqv2/hvmon.sock, modo 600):
// una petición JSON por línea, una respuesta JSON por línea, una petición por conexión.
//   {"cmd":"ping"}
//   {"cmd":"status"}                                  estado completo (como hvmon_state.json)
//   {"cmd":"read"}                                    última lectura de todos los parámetros
//   {"cmd":"changes","n":20}                          últimos n cambios (máximo el anillo, 50)
//   {"cmd":"note","text":"..."}                       nota del operador (notes_YYYYMMDD.jsonl)
//   {"cmd":"measure-offsets","seconds":60,"channels":[...]}   solo canales con Pw = 0
//   {"cmd":"pw","ch":3,"on":true,"authorize":"motivo"} única escritura; requiere allow_commands
// Respuesta: {"ok":true,...} o {"ok":false,"error":"..."}. Las peticiones largas (pw,
// measure-offsets) se atienden en un hilo propio para no bloquear status/read.
// CommandHandler::handle() nunca lanza excepciones; CommandServer es el transporte.
// Seguridad entre hilos: handle() se puede llamar desde varios hilos a la vez.
#pragma once

#include <atomic>
#include <list>
#include <memory>
#include <string>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/core/jthread.h"
#include "daq/core/unix_socket.h"
#include "daq/hvmon/config.h"
#include "daq/hvmon/poller.h"
#include "daq/hvmon/state.h"
#include "daq/hvmon/writer.h"

namespace daq::hvmon {

class CommandHandler {
 public:
  CommandHandler(const HvmonConfig& cfg, StateHub& hub, BoundedQueue<PollRequest>& poll, BoundedQueue<WriterItem>& writer,
                 const Clock& clock);
  nlohmann::json handle(const std::string& line, const PeerCredentials& peer);
  static bool isLong(const std::string& line);

 private:
  nlohmann::json dispatch(const nlohmann::json& req, const PeerCredentials& peer);
  nlohmann::json forward(PollRequest r, double timeoutS);
  nlohmann::json readJson() const;

  const HvmonConfig& cfg_;
  StateHub& hub_;
  BoundedQueue<PollRequest>& poll_;
  BoundedQueue<WriterItem>& writer_;
  const Clock& clock_;
  std::atomic<std::uint64_t> statusSeq_{0};
};

class CommandServer {
 public:
  CommandServer(UnixListener listener, CommandHandler& handler);
  void run(const StopToken& st);

 private:
  UnixListener listener_;
  CommandHandler& handler_;
  struct Worker {
    std::shared_ptr<std::atomic<bool>> done;
    JThread thread;
  };
  std::list<Worker> workers_;
};

std::string describePeer(const PeerCredentials& peer);

}  // namespace daq::hvmon
