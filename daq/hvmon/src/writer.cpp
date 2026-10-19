// writer.cpp
#include "daq/hvmon/writer.h"

#include <cstdio>
#include <ctime>

#include "daq/core/logger.h"
#include "daq/hv/reading_json.h"

namespace daq::hvmon {

namespace {

std::string fmtDouble(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.9g", v);
  return buf;
}

nlohmann::json col(const std::string& name, const std::string& type, const std::string& unit, const std::string& desc) {
  return nlohmann::json{{"name", name}, {"type", type}, {"unit", unit}, {"description", desc}};
}

nlohmann::json jsonlSchema(const std::string& name, const std::string& description, nlohmann::json fields) {
  return nlohmann::json{{"schema", name}, {"format", "jsonl"}, {"time", "UTC ns desde la época Unix"},
                        {"description", description}, {"fields", std::move(fields)}};
}

const nlohmann::json& changesSchema() {
  static const nlohmann::json s = jsonlSchema(
      "daqv2.hv_changes/1",
      "Cambios de parámetros HV (type=change) y registros de comandos de escritura (type=command)",
      nlohmann::json::array(
          {col("type", "string", "", "change | command"),
           col("id", "uint64", "", "id del cambio (type=change)"),
           col("command_id", "uint64|null", "", "id del comando (command; change con origin=hvmon_command)"),
           col("t_utc_ns", "int64", "ns", "UTC de la lectura que mostró el valor nuevo / del registro"),
           col("t_prev_read_utc_ns", "int64", "ns", "UTC de la lectura anterior (el cambio ocurrió entre ambas)"),
           col("slot", "int", "", "slot del crate"), col("ch", "int", "", "canal HV"),
           col("param", "string", "", "parámetro"), col("unit", "string", "", "unidad SI de old/new"),
           col("old", "number", "SI", "valor anterior"), col("new", "number", "SI", "valor nuevo"),
           col("origin", "string", "", "trip | hvmon_command | external"),
           col("during_disconnect", "bool", "", "la lectura anterior es previa a una pérdida de conexión"),
           col("read_kind", "string", "", "fast | full"),
           col("phase", "string", "", "command: request|rejected|aborted|noop|written|write_failed|verified|verify_failed"),
           col("value", "int", "", "command: valor pedido de Pw"),
           col("reason", "string", "", "command: texto de --authorize"),
           col("requester", "string", "", "command: uid/pid del cliente"),
           col("detail", "string", "", "command: detalle")}));
  return s;
}

const nlohmann::json& notesSchema() {
  static const nlohmann::json s = jsonlSchema(
      "daqv2.notes/1", "Notas del operador (hvctl note / daqnote)",
      nlohmann::json::array({col("t_utc_ns", "int64", "ns", "UTC de la nota"), col("user", "string", "", "usuario"),
                             col("host", "string", "", "host"), col("text", "string", "", "texto")}));
  return s;
}

const nlohmann::json& tempSchema() {
  static const nlohmann::json s{
      {"schema", "daqv2.temp/1"},
      {"format", "csv"},
      {"columns", nlohmann::json::array({col("t_utc_ns", "int64", "ns", "UTC de la lectura"),
                                         col("sensor_id", "string", "", "identificador del sensor"),
                                         col("T_degC", "float64", "degC", "temperatura")})}};
  return s;
}

}  // namespace

std::string readingHeader(const hv::ChannelReading& r, const hv::BoardReading* board) {
  std::string h = "t_req_utc_ns,t_resp_utc_ns,read_seq,slot,ch";
  for (const auto& p : r.params) h += "," + hv::columnName(p);
  h += ",read_duration_ns";
  if (board)
    for (const auto& p : board->params) h += ",bd_" + hv::columnName(p);
  return h;
}

nlohmann::json readingSchema(const std::string& schemaName, const hv::ChannelReading& r,
                             const hv::BoardReading* board) {
  nlohmann::json cols = nlohmann::json::array(
      {col("t_req_utc_ns", "int64", "ns", "UTC antes de la primera llamada de la lectura"),
       col("t_resp_utc_ns", "int64", "ns", "UTC después de la última llamada"),
       col("read_seq", "uint64", "", "número de lectura de hvmon"), col("slot", "int", "", "slot del crate"),
       col("ch", "int", "", "canal HV")});
  for (const auto& p : r.params) {
    if (p.type == hv::ParamType::Numeric)
      cols.push_back(col(hv::columnName(p), "float64", p.siUnit(), p.name + " (placa: " + p.boardUnit() + ")"));
    else
      cols.push_back(col(p.name, "uint32", "", std::string(p.name) + " (" + hv::toString(p.type) + ")"));
  }
  cols.push_back(col("read_duration_ns", "int64", "ns", "duración de la lectura (steady_clock)"));
  if (board) {
    for (const auto& p : board->params) {
      const std::string extra = p.name == "Temp" ? " (temperatura de la placa A7030DP, no del detector)" : "";
      if (p.type == hv::ParamType::Numeric)
        cols.push_back(col("bd_" + hv::columnName(p), "float64", p.siUnit(),
                           "placa: " + p.name + " (" + p.boardUnit() + ")" + extra));
      else
        cols.push_back(col("bd_" + p.name, "uint32", "", "placa: " + p.name + " (" + hv::toString(p.type) + ")"));
    }
  }
  return nlohmann::json{{"schema", schemaName}, {"format", "csv"}, {"columns", cols}};
}

std::vector<std::string> readingRows(const hv::ChannelReading& r, std::uint64_t seq, const hv::BoardReading* board) {
  std::vector<std::string> rows;
  for (std::size_t c = 0; c < r.channels.size(); ++c) {
    std::string row = std::to_string(r.tReqUtcNs) + "," + std::to_string(r.tRespUtcNs) + "," + std::to_string(seq) +
                      "," + std::to_string(r.slot) + "," + std::to_string(r.channels[c]);
    for (std::size_t p = 0; p < r.params.size(); ++p) {
      const auto& v = r.values[p][c];
      row += "," + (v.isNumeric() ? fmtDouble(v.si) : std::to_string(v.raw));
    }
    row += "," + std::to_string(r.durationNs);
    if (board)
      for (const auto& v : board->values) row += "," + (v.isNumeric() ? fmtDouble(v.si) : std::to_string(v.raw));
    rows.push_back(std::move(row));
  }
  return rows;
}

std::string SlowControlWriter::dateKey(UtcNs t) {
  std::time_t sec = static_cast<std::time_t>(t / 1000000000LL);
  if (t < 0 && t % 1000000000LL != 0) --sec;
  std::tm tm{};
  gmtime_r(&sec, &tm);
  char buf[16];
  std::snprintf(buf, sizeof buf, "%04d%02d%02d", (tm.tm_year + 1900) % 10000, (tm.tm_mon + 1) % 100, tm.tm_mday % 100);
  return buf;
}

SlowControlWriter::SlowControlWriter(std::string baseDir, double fsyncPeriodS, FileSystem& fs, const Clock& clock,
                                     StateHub& hub)
    : baseDir_(std::move(baseDir)), fsyncPeriodS_(fsyncPeriodS), fs_(fs), clock_(clock), hub_(hub) {
  nextSync_ = clock_.now_steady() + std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::duration<double>(fsyncPeriodS_));
}

SlowControlWriter::~SlowControlWriter() { closeAll(); }

void SlowControlWriter::fail(const std::string& message) {
  health_.ok = false;
  health_.lastError = message;
  health_.lastErrorUtcNs = clock_.now_utc_ns();
  ++health_.errors;
  if (message != lastLoggedError_) {
    logError("escritor: ", message);
    lastLoggedError_ = message;
  }
}

void SlowControlWriter::publish() { hub_.writer.publish(health_); }

bool SlowControlWriter::openStream(Stream& s, const std::string& path, const std::string& header,
                                   const nlohmann::json& schema) {
  if (auto st = fs_.makeDirectories(parentDirectory(path), 0755); !st) {
    fail(st.error().message);
    return false;
  }
  auto f = fs_.openAppend(path);
  if (!f) {
    fail(f.error().message);
    return false;
  }
  s.file = std::move(f).value();
  s.path = path;
  s.dirty = false;
  if (!header.empty() && s.file->size() == 0) writeTo(s, header);
  if (!schema.is_null()) {
    const std::string schemaPath = path + ".schema.json";
    if (!fs_.exists(schemaPath)) {
      auto j = schema;
      j["file"] = path.substr(path.find_last_of('/') + 1);
      if (auto st = fs_.writeAtomic(schemaPath, j.dump(2) + "\n", 0644); !st) fail(st.error().message);
    }
  }
  return s.file != nullptr;
}

void SlowControlWriter::writeTo(Stream& s, const std::string& line) {
  if (!s.file) return;
  const std::string data = line + "\n";
  if (auto st = s.file->append(data); !st) {
    fail(st.error().message);
    s.file.reset();  // se reabre con el próximo dato
    return;
  }
  s.dirty = true;
  s.bytes += data.size();
  health_.bytes += data.size();
  ++health_.lines;
  health_.openFiles[s.path] = s.bytes;
  health_.ok = true;  // ok = la última escritura funcionó; errors/last_error conservan la historia
}

void SlowControlWriter::syncStream(Stream& s) {
  if (!s.file || !s.dirty) return;
  if (auto st = s.file->sync(); !st) {
    fail(st.error().message);
    return;
  }
  s.dirty = false;
  ++health_.fsyncs;
}

void SlowControlWriter::writeDaily(const std::string& kind, const std::string& ext, UtcNs t, const std::string& header,
                                   const nlohmann::json& schema, const std::string& line) {
  const std::string date = dateKey(t);
  Stream& s = streams_[kind];
  if (s.file && s.date != date) {  // rotación a medianoche UTC
    syncStream(s);
    health_.openFiles.erase(s.path);
    s.file.reset();
    if (compressQ_ && !compressQ_->try_push(s.path, DropPolicy::DropNewest))
      logWarn("cola de compresión llena: ", s.path, " queda sin comprimir");
  }
  if (!s.file) {
    const std::string path = baseDir_ + "/" + date.substr(0, 4) + "/" + kind + "_" + date + "." + ext;
    if (s.path != path) s.bytes = 0;
    s.date = date;
    if (!openStream(s, path, header, schema)) return;
  }
  writeTo(s, line);
}

void SlowControlWriter::handle(const WriterItem& item) {
  if (const auto* f = std::get_if<FastRows>(&item)) {
    for (const auto& row : readingRows(f->reading, f->seq))
      writeDaily("hv", "csv", f->reading.tRespUtcNs, readingHeader(f->reading),
                 readingSchema("daqv2.hv_fast/1", f->reading), row);
  } else if (const auto* g = std::get_if<FullRows>(&item)) {
    const hv::BoardReading* bd = g->board ? &*g->board : nullptr;
    for (const auto& row : readingRows(g->reading, g->seq, bd))
      writeDaily("hv_full", "csv", g->reading.tRespUtcNs, readingHeader(g->reading, bd),
                 readingSchema("daqv2.hv_full/2", g->reading, bd), row);
  } else if (const auto* c = std::get_if<ChangeItem>(&item)) {
    writeDaily("hv_changes", "jsonl", c->rec.tUtcNs, "", changesSchema(), c->rec.toJson().dump());
  } else if (const auto* k = std::get_if<CommandItem>(&item)) {
    writeDaily("hv_changes", "jsonl", k->rec.tUtcNs, "", changesSchema(), k->rec.toJson().dump());
  } else if (const auto* n = std::get_if<NoteItem>(&item)) {
    nlohmann::json j{{"t_utc_ns", n->tUtcNs}, {"t_utc", formatUtcIso(n->tUtcNs)}, {"user", n->user},
                     {"host", n->host}, {"text", n->text}};
    writeDaily("notes", "jsonl", n->tUtcNs, "", notesSchema(), j.dump());
  } else if (const auto* t = std::get_if<TempItem>(&item)) {
    for (const auto& r : t->readings)
      writeDaily("temp", "csv", r.tUtcNs, "t_utc_ns,sensor_id,T_degC", tempSchema(),
                 std::to_string(r.tUtcNs) + "," + r.sensorId + "," + fmtDouble(r.tDegC));
  } else if (const auto* o = std::get_if<OffsetItem>(&item)) {
    if (!offsets_.file) openStream(offsets_, baseDir_ + "/offsets.jsonl", "", nlohmann::json());
    writeTo(offsets_, o->j.dump());
    syncStream(offsets_);  // las medidas de offsets son pocas y valiosas
  }
  publish();
}

void SlowControlWriter::maybeSync() {
  const SteadyTime now = clock_.now_steady();
  if (now < nextSync_) return;
  syncAll();
  nextSync_ = now + std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(fsyncPeriodS_));
}

void SlowControlWriter::syncAll() {
  for (auto& kv : streams_) syncStream(kv.second);
  syncStream(offsets_);
  publish();
}

void SlowControlWriter::closeAll() {
  syncAll();
  for (auto& kv : streams_) kv.second.file.reset();
  offsets_.file.reset();
  health_.openFiles.clear();
  publish();
}

void SlowControlWriter::run(BoundedQueue<WriterItem>& queue) {
  for (;;) {
    auto item = queue.pop_until(std::chrono::steady_clock::now() + std::chrono::milliseconds(200));
    if (item) handle(*item);
    else if (queue.closed()) break;
    maybeSync();
  }
  closeAll();
}

}  // namespace daq::hvmon
