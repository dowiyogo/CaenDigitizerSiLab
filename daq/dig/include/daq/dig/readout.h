// readout.h
// Bucle de lectura de la placa (diseño 3.1), escrito contra DigitizerBackend: la fase 4 solo agrega
// el backend real. Sondeo: readBlock(max_events_blt); si no hay datos, espera poll_sleep_s (5 ms).
// Cada bloque se publica como shared_ptr<const EventBlock>:
//  - cola de escritura: push bloqueante, nunca se descarta; si está llena la lectura espera (el
//    digitizador acumula en su memoria) y se cuenta (writer_blocked, writer_blocked_s);
//  - cola del monitor: try_push con descarte del más antiguo (monitor_dropped).
// Modo calibración: envía SW triggers a sw_trigger_rate_Hz entre lecturas.
// Fin: parada pedida, max_events, max_duration_s o fuente agotada (replay). Al terminar: stop(),
// vaciado de la placa hasta 0 eventos, y cierre de ambas colas.
// Seguridad entre hilos: run() corre en el hilo de lectura (único dueño del backend); las
// estadísticas son atómicas y se leen desde cualquier hilo.
#pragma once

#include <atomic>
#include <cstdint>
#include <mutex>
#include <string>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/core/jthread.h"
#include "daq/dig/digitizer.h"

namespace daq::dig {

struct ReadoutOptions {
  double pollSleepS = 0.005;
  std::size_t maxEventsBlt = 16;
  std::uint64_t maxEvents = 0;  // 0 = sin límite
  double maxDurationS = 0;      // 0 = sin límite
  bool softwareTrigger = false;
  double swTriggerRateHz = 100;
};

struct ReadoutStats {
  std::atomic<std::uint64_t> blocks{0}, events{0}, emptyPolls{0}, readErrors{0};
  std::atomic<std::uint64_t> writerBlocked{0};     // bloques que encontraron la cola de escritura llena
  std::atomic<std::uint64_t> writerBlockedNs{0};   // tiempo total esperando en push
  std::atomic<std::uint64_t> writerMaxOccupancy{0};
  std::atomic<std::uint64_t> monitorDropped{0};
  std::atomic<std::uint64_t> monitorMaxOccupancy{0};
  std::atomic<std::uint64_t> boardFullBlocks{0};
  std::atomic<std::uint64_t> maxBlockEvents{0};
  std::atomic<std::uint64_t> swTriggersSent{0};
  std::atomic<bool> finished{false};
  std::string finishReason() const;
  void setFinishReason(const std::string& s);
  std::string lastError() const;
  void setLastError(const std::string& s);
  nlohmann::json toJson() const;

 private:
  mutable std::mutex m_;
  std::string finishReason_, lastError_;
};

class ReadoutLoop {
 public:
  ReadoutLoop(DigitizerBackend& backend, Clock& clock, BoundedQueue<EventBlockPtr>& writerQ,
              BoundedQueue<EventBlockPtr>& monitorQ, ReadoutOptions opt, ReadoutStats& stats);
  // Corre hasta parada o fin; devuelve error si la placa falla de forma irrecuperable.
  Status run(const StopToken& st);

 private:
  void publish(EventBlockPtr b);
  DigitizerBackend& be_;
  Clock& clock_;
  BoundedQueue<EventBlockPtr>& writerQ_;
  BoundedQueue<EventBlockPtr>& monitorQ_;
  ReadoutOptions opt_;
  ReadoutStats& st_;
};

}  // namespace daq::dig
