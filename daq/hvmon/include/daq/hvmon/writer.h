// writer.h
// Hilo escritor de hvmon: archivos diarios de control lento bajo <slowcontrol_dir>/YYYY/, con
// rotación a medianoche UTC según la marca de tiempo de cada dato (no la hora de escritura):
//   hv_YYYYMMDD.csv          una fila por canal por lectura rápida
//   hv_full_YYYYMMDD.csv     una fila por canal por lectura completa (cada 60 s), con los
//                            parámetros de placa (bd_BdStatus, bd_HVMax_V, bd_HIMax_A, bd_Temp_degC)
//   hv_changes_YYYYMMDD.jsonl cambios y registros de comandos, uno por línea
//   temp_YYYYMMDD.csv        temperaturas
//   notes_YYYYMMDD.jsonl     notas del operador
// y <slowcontrol_dir>/offsets.jsonl (medidas de offsets, sin rotación).
// Cada archivo diario tiene un hermano <archivo>.schema.json (columnas, tipos, unidades SI).
// Cada línea se escribe con una llamada write() (O_APPEND) en cuanto llega; fsync cada
// fsync_period_s (por defecto 10 s) sobre los archivos con datos nuevos, al rotar y al cerrar.
// Los errores de escritura se registran en el log y en WriterHealth (publicado en el estado);
// nunca se descartan en silencio.
// Seguridad entre hilos: la instancia pertenece al hilo escritor; recibe datos por BoundedQueue.
#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/bounded_queue.h"
#include "daq/core/clock.h"
#include "daq/core/file_system.h"
#include "daq/hv/hv_types.h"
#include "daq/hvmon/state.h"
#include "daq/temp/temp_backend.h"

namespace daq::hvmon {

struct FastRows {
  hv::ChannelReading reading;
  std::uint64_t seq = 0;
};
struct FullRows {
  hv::ChannelReading reading;
  std::uint64_t seq = 0;
  std::optional<hv::BoardReading> board;  // columnas bd_* repetidas en cada fila de canal
};
struct ChangeItem {
  ChangeRecord rec;
};
struct CommandItem {
  CommandRecord rec;
};
struct OffsetItem {
  nlohmann::json j;
};
struct NoteItem {
  UtcNs tUtcNs = 0;
  std::string user;
  std::string host;
  std::string text;
};
struct TempItem {
  UtcNs tUtcNs = 0;
  std::vector<temp::TempReading> readings;
};
using WriterItem = std::variant<FastRows, FullRows, ChangeItem, CommandItem, OffsetItem, NoteItem, TempItem>;

class SlowControlWriter {
 public:
  SlowControlWriter(std::string baseDir, double fsyncPeriodS, FileSystem& fs, const Clock& clock, StateHub& hub);
  ~SlowControlWriter();

  void handle(const WriterItem& item);
  void maybeSync();  // fsync si venció el período
  void syncAll();    // fsync de todo lo que tenga datos nuevos
  void closeAll();   // syncAll + cerrar
  // Bucle del hilo: atiende la cola hasta que se cierre y se vacíe; luego syncAll y closeAll.
  void run(BoundedQueue<WriterItem>& queue);

  const WriterHealth& health() const { return health_; }
  // Fase 2: al rotar, encola los archivos del día cerrado para comprimirlos (nullptr = no).
  void setCompressionQueue(BoundedQueue<std::string>* q) { compressQ_ = q; }
  static std::string dateKey(UtcNs t);  // "YYYYMMDD" en UTC

 private:
  struct Stream {
    std::string path;
    std::string date;
    std::unique_ptr<WritableFile> file;
    bool dirty = false;
    std::uint64_t bytes = 0;
  };
  // Escribe una línea en el archivo diario `kind` para la fecha de t; header/schema al crearlo.
  void writeDaily(const std::string& kind, const std::string& ext, UtcNs t, const std::string& header,
                  const nlohmann::json& schema, const std::string& line);
  void writeTo(Stream& s, const std::string& line);
  bool openStream(Stream& s, const std::string& path, const std::string& header, const nlohmann::json& schema);
  void syncStream(Stream& s);
  void fail(const std::string& message);
  void publish();

  std::string baseDir_;
  double fsyncPeriodS_;
  FileSystem& fs_;
  const Clock& clock_;
  StateHub& hub_;
  SteadyTime nextSync_;
  std::map<std::string, Stream> streams_;  // kind -> archivo abierto
  Stream offsets_;
  WriterHealth health_;
  std::string lastLoggedError_;
  BoundedQueue<std::string>* compressQ_ = nullptr;
};

// Esquemas (también se escriben como .schema.json).
nlohmann::json readingSchema(const std::string& schemaName, const hv::ChannelReading& r,
                             const hv::BoardReading* board = nullptr);
std::string readingHeader(const hv::ChannelReading& r, const hv::BoardReading* board = nullptr);
std::vector<std::string> readingRows(const hv::ChannelReading& r, std::uint64_t seq,
                                     const hv::BoardReading* board = nullptr);

}  // namespace daq::hvmon
