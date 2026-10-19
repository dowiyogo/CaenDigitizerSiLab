// notifier.cpp
#include "daq/hvmon/notifier.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <dirent.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <sstream>

#include "daq/core/file_system.h"

#include "daq/core/logger.h"
#include "daq/hvmon/state.h"

namespace daq::hvmon {

namespace {

std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t");
  if (b == std::string::npos) return "";
  return s.substr(b, s.find_last_not_of(" \t") - b + 1);
}

// Valor entre comillas para el archivo de configuración de curl (escapa \ y ").
std::string curlQuote(const std::string& v) {
  std::string out = "\"";
  for (char c : v) {
    if (c == '\\' || c == '"') out += '\\';
    if (c == '\n' || c == '\r') continue;
    out += c;
  }
  return out + "\"";
}

// Dirección pura de "Nombre <dir>" o "dir".
std::string bareAddress(const std::string& a) {
  const auto l = a.find('<'), r = a.find('>');
  if (l != std::string::npos && r != std::string::npos && r > l) return trim(a.substr(l + 1, r - l - 1));
  return trim(a);
}

bool isAscii(const std::string& s) {
  for (unsigned char c : s)
    if (c >= 0x80) return false;
  return true;
}

std::string base64(const std::string& in) {
  static const char* t = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string out;
  std::size_t i = 0;
  while (i + 2 < in.size()) {
    const unsigned v = (static_cast<unsigned char>(in[i]) << 16) | (static_cast<unsigned char>(in[i + 1]) << 8) |
                       static_cast<unsigned char>(in[i + 2]);
    out += t[(v >> 18) & 63];
    out += t[(v >> 12) & 63];
    out += t[(v >> 6) & 63];
    out += t[v & 63];
    i += 3;
  }
  if (i < in.size()) {
    unsigned v = static_cast<unsigned char>(in[i]) << 16;
    if (i + 1 < in.size()) v |= static_cast<unsigned char>(in[i + 1]) << 8;
    out += t[(v >> 18) & 63];
    out += t[(v >> 12) & 63];
    out += (i + 1 < in.size()) ? t[(v >> 6) & 63] : '=';
    out += '=';
  }
  return out;
}

// Encabezado con UTF-8 codificado según RFC 2047 si hace falta.
std::string headerText(const std::string& s) { return isAscii(s) ? s : "=?UTF-8?B?" + base64(s) + "?="; }

std::string describeEvent(const std::string& kind) {
  if (kind == "gave_up") return "GAVE_UP: reintentos agotados, el canal queda apagado";
  if (kind == "trip") return "TRIP del canal HV";
  if (kind == "recovery") return "RECUPERADO tras un reintento automático";
  if (kind == "writer_error") return "ERROR de escritura de los archivos de control lento";
  if (kind == "connection_lost") return "CONEXIÓN PERDIDA con el SY5527";
  if (kind == "hvmon_stale") return "hvmon NO RESPONDE durante un run de daqrun (bit 4 en todos los canales)";
  if (kind == "no_triggers") return "SIN TRIGGERS externos durante un run de daqrun";
  if (kind == "triggers_resumed") return "Los triggers externos VOLVIERON";
  return kind;
}

// "2026-10-19T09:27:13.740Z (2026-10-19 06:27:13 -03)": UTC con ms y hora local del proceso.
std::string utcAndLocal(UtcNs t) {
  const std::string iso = formatUtcIso(t);  // YYYY-MM-DDThh:mm:ss.nnnnnnnnnZ
  return iso.substr(0, 23) + "Z (" + formatLocalShort(t) + ")";
}

}  // namespace

std::string defaultOutboxDir() {
  const char* x = std::getenv("XDG_STATE_HOME");
  if (x && *x) return std::string(x) + "/daqv2/outbox";
  const char* home = std::getenv("HOME");
  return std::string(home ? home : "") + "/.local/state/daqv2/outbox";
}

nlohmann::json NotifyEvent::toJson() const {
  return {{"kind", kind},       {"ch", ch},           {"slot", slot},       {"episode", episode},
          {"k", k},             {"reason", reason},   {"timeline", timeline}, {"context", context},
          {"t_utc_ns", tUtcNs}, {"simulated", simulated}, {"origin", origin}};
}

NotifyEvent NotifyEvent::fromJson(const nlohmann::json& j) {
  NotifyEvent e;
  e.kind = j.value("kind", "");
  e.ch = j.value("ch", -1);
  e.slot = j.value("slot", 3);
  e.episode = j.value("episode", std::uint64_t{0});
  e.k = j.value("k", 0);
  e.reason = j.value("reason", "");
  e.timeline = j.value("timeline", nlohmann::json::array());
  e.context = j.value("context", nlohmann::json::object());
  e.tUtcNs = j.value("t_utc_ns", UtcNs{0});
  e.simulated = j.value("simulated", false);
  e.origin = j.value("origin", "");
  return e;
}

Result<NotifierConfig> NotifierConfig::fromJson(const nlohmann::json& j) {
  using R = Result<NotifierConfig>;
  NotifierConfig c;
  c.outboxDir = defaultOutboxDir();
  try {
    c.enabled = j.value("enabled", c.enabled);
    c.notifyEachTrip = j.value("notify_each_trip", c.notifyEachTrip);
    c.notifyRecovery = j.value("notify_recovery", c.notifyRecovery);
    c.notifyWriterErrors = j.value("notify_writer_errors", c.notifyWriterErrors);
    c.notifyConnectionLost = j.value("notify_connection_lost", c.notifyConnectionLost);
    c.connectionLostMin = j.value("connection_lost_min", c.connectionLostMin);
    c.problemIntervalMin = j.value("problem_interval_min", c.problemIntervalMin);
    c.maxMessages = j.value("max_messages", c.maxMessages);
    c.retryDelaysMin = j.value("retry_delays_min", c.retryDelaysMin);
    c.curlPath = j.value("curl", c.curlPath);
    c.maxTimeS = j.value("max_time_s", c.maxTimeS);
    c.mappingFile = j.value("mapping_file", c.mappingFile);
    c.retryHourlyMin = j.value("retry_hourly_min", c.retryHourlyMin);
    c.maxAgeH = j.value("max_age_h", c.maxAgeH);
    c.outboxDir = j.value("outbox_dir", c.outboxDir);
    c.localTimezone = j.value("local_timezone", c.localTimezone);
    c.pageSshUser = j.value("page_ssh_user", c.pageSshUser);
    c.pageSshHost = j.value("page_ssh_host", c.pageSshHost);
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("notifier inválido: ") + e.what());
  }
  if (c.maxMessages < 0) return R::failure(-1, "notifier: max_messages debe ser >= 0");
  if (c.maxTimeS <= 0 || c.connectionLostMin <= 0 || c.problemIntervalMin <= 0 || c.retryHourlyMin <= 0 || c.maxAgeH <= 0)
    return R::failure(-1, "notifier: max_time_s, connection_lost_min y problem_interval_min deben ser positivos");
  return c;
}

nlohmann::json NotifierConfig::toJson() const {
  return nlohmann::json{{"enabled", enabled},           {"notify_each_trip", notifyEachTrip},
                        {"notify_recovery", notifyRecovery}, {"notify_writer_errors", notifyWriterErrors},
                        {"notify_connection_lost", notifyConnectionLost}, {"connection_lost_min", connectionLostMin},
                        {"problem_interval_min", problemIntervalMin}, {"retry_delays_min", retryDelaysMin},
                        {"curl", curlPath},             {"max_time_s", maxTimeS},
                        {"mapping_file", mappingFile}, {"local_timezone", localTimezone},
                        {"page_ssh_user", pageSshUser}, {"page_ssh_host", pageSshHost},
                        {"retry_hourly_min", retryHourlyMin}, {"max_age_h", maxAgeH}, {"outbox_dir", outboxDir},
                        {"max_messages", maxMessages}};
}

Result<SmtpSettings> SmtpSettings::fromSecrets(Secrets& s) {
  using R = Result<SmtpSettings>;
  SmtpSettings out;
  struct Item {
    const char* key;
    std::string* dst;
  };
  std::string to;
  for (Item it : {Item{"DAQ_SMTP_URL", &out.url}, Item{"DAQ_SMTP_USER", &out.user},
                  Item{"DAQ_SMTP_PASSWORD", &out.password}, Item{"DAQ_NOTIFY_FROM", &out.from},
                  Item{"DAQ_NOTIFY_TO", &to}}) {
    auto v = s.get(it.key);
    if (!v) return R(v.error());
    if (!v.value()) return R::failure(-1, std::string("falta ") + it.key + " (entorno o " + s.path() + ")");
    *it.dst = *v.value();
  }
  std::stringstream ss(to);
  std::string a;
  while (std::getline(ss, a, ','))
    if (!trim(a).empty()) out.to.push_back(trim(a));
  if (out.to.empty()) return R::failure(-1, "DAQ_NOTIFY_TO sin direcciones");
  return out;
}

Notifier::Notifier(NotifierConfig cfg, SmtpSettings smtp, std::optional<Mapping> mapping, Spawner& spawner,
                   Clock& clock, std::string tmpDir, StateHub& hub, std::string hostName)
    : cfg_(std::move(cfg)),
      smtp_(std::move(smtp)),
      mapping_(std::move(mapping)),
      spawner_(spawner),
      clock_(clock),
      tmpDir_(std::move(tmpDir)),
      hub_(hub),
      host_(std::move(hostName)) {
  if (host_.empty()) {
    char buf[256] = {0};
    host_ = ::gethostname(buf, sizeof buf - 1) == 0 ? buf : "localhost";
  }
  lastWriterErrors_ = hub_.writer.read().first.errors;
  loadOutbox();
  publish();
}

// Un archivo JSON por correo pendiente, modo 600, sin secretos (el evento; el mensaje y la
// configuración de curl se arman al enviar). Sobrevive a reinicios de hvmon.
void Notifier::persist(const Pending& p) {
  if (p.path.empty()) return;
  PosixFileSystem fs;
  const nlohmann::json j{{"schema", "daqv2.notify_outbox/1"}, {"event", p.ev.toJson()},  {"attempts", p.attempts},
                         {"queued_utc_ns", p.queuedUtcNs},     {"due_utc_ns", p.dueUtcNs}};
  if (auto st = fs.writeAtomic(p.path, j.dump() + "\n", 0600); !st) logError("outbox: ", st.error().message);
}

void Notifier::loadOutbox() {
  if (cfg_.outboxDir.empty()) return;
  PosixFileSystem fs;
  if (auto st = fs.makeDirectories(cfg_.outboxDir, 0700); !st) {
    logError("outbox: ", st.error().message);
    return;
  }
  DIR* d = ::opendir(cfg_.outboxDir.c_str());
  if (!d) return;
  std::vector<std::string> files;
  while (dirent* e = ::readdir(d)) {
    const std::string n = e->d_name;
    if (n.size() > 5 && n.compare(n.size() - 5, 5, ".json") == 0) files.push_back(cfg_.outboxDir + "/" + n);
  }
  ::closedir(d);
  std::sort(files.begin(), files.end());
  const UtcNs now = clock_.now_utc_ns();
  for (const auto& f : files) {
    auto t = fs.readFile(f);
    if (!t) continue;
    try {
      const auto j = nlohmann::json::parse(t.value());
      Pending p;
      p.ev = NotifyEvent::fromJson(j.at("event"));
      p.attempts = j.value("attempts", 0);
      p.queuedUtcNs = j.value("queued_utc_ns", now);
      p.dueUtcNs = std::max<UtcNs>(now, j.value("due_utc_ns", now));
      p.path = f;
      pending_.push_back(std::move(p));
      ++restored_;
      ++accepted_;
    } catch (const std::exception& ex) {
      logError("outbox: archivo ilegible ", f, ": ", ex.what());
    }
  }
  if (restored_) logInfo("notificador: ", restored_, " correos pendientes recuperados del outbox ", cfg_.outboxDir);
}

std::string Notifier::redact(std::string s) const {
  if (smtp_.password.empty()) return s;
  for (auto p = s.find(smtp_.password); p != std::string::npos; p = s.find(smtp_.password))
    s.replace(p, smtp_.password.size(), "***");
  return s;
}

std::string Notifier::rfc5322Date(UtcNs t) const {
  std::time_t sec = static_cast<std::time_t>(t / 1000000000LL);
  std::tm tm{};
  gmtime_r(&sec, &tm);
  static const char* days[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};
  static const char* months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  char buf[64];
  std::snprintf(buf, sizeof buf, "%s, %02d %s %04d %02d:%02d:%02d +0000", days[tm.tm_wday % 7], tm.tm_mday,
                months[tm.tm_mon % 12], (tm.tm_year + 1900) % 10000, tm.tm_hour, tm.tm_min, tm.tm_sec);
  return buf;
}

std::string Notifier::groupText(int hvCh) const {
  std::ostringstream os;
  if (!mapping_) {
    os << "  (mapeo no disponible: " << cfg_.mappingFile << ")\n";
    return os.str();
  }
  const auto g = mapping_->groupOfHvChannel(hvCh);
  if (!g) {
    os << "  (CH" << hvCh << " no corresponde a un grupo de MPPC: canal libre)\n";
    return os.str();
  }
  for (const char* st : {"hod2018", "hod2019"}) {
    os << "  " << st << ":";
    int n = 0;
    for (const auto& e : mapping_->forGroup(*g)) {
      if (e.station != st) continue;
      os << " ch" << e.digCh << (e.groupConfirmed() ? "" : "*") << "(" << e.plane << e.strip
         << (e.face.empty() ? "" : " " + e.face) << ")";
      ++n;
    }
    if (n == 0) os << " (ninguno conocido)";
    const int unknown = mapping_->unknownGroupCount(st);
    if (unknown > 0) os << "  [" << unknown << " canales de " << st << " sin grupo conocido]";
    os << "\n";
  }
  os << "  * = asignación no confirmada por medición (mapeo " << mapping_->version << ")\n";
  return os.str();
}

std::string Notifier::composeMessage(const NotifyEvent& e, std::string* subjectOut) const {
  // Marca de simulación: [PRUEBA] (hvctl notify-test) o [SIMULADO] (FakeHVBackend o canal en
  // test_mode) al inicio del asunto y en la primera línea del cuerpo.
  const bool simulated = e.simulated || simulatedBackend_;
  std::ostringstream subj;
  if (e.test)
    subj << "[PRUEBA] ";
  else if (simulated)
    subj << "[SIMULADO] ";
  subj << "[DAQ hodoscopio] ";
  int group = -1;  // -1: canal sin grupo de MPPC o sin mapeo
  if (mapping_ && e.ch >= 0)
    if (auto g = mapping_->groupOfHvChannel(e.ch)) group = *g;
  if (e.ch >= 0) {
    subj << e.kind << " CH" << (e.ch < 10 ? "0" : "") << e.ch;
    if (group > 0) subj << " (grupo " << group << ")";
  } else {
    subj << e.kind;
  }
  if (!e.origin.empty()) subj << " — " << e.origin;
  const std::string subject = subj.str();
  if (subjectOut) *subjectOut = subject;

  std::ostringstream body;
  if (e.test) {
    body << "*** PRUEBA *** Correo de prueba del notificador (hvctl notify-test): datos sintéticos, "
            "ningún canal real estuvo involucrado.\n";
    if (simulatedBackend_) body << "Fuente HV simulada (FakeHVBackend).\n";
    body << "\n";
  } else if (simulated) {
    body << "*** SIMULADO *** Evento simulado ("
         << (simulatedBackend_ ? std::string("FakeHVBackend")
                               : "test_mode del canal CH" + std::string(e.ch < 10 ? "0" : "") + std::to_string(e.ch))
         << "): ningún canal real estuvo involucrado.\n\n";
  }
  body << "Evento: " << describeEvent(e.kind) << "\n";
  body << "Hora del evento: " << utcAndLocal(e.tUtcNs) << "\n";
  body << "Hora de envío: " << utcAndLocal(clock_.now_utc_ns()) << "\n";
  {
    const double delayS = static_cast<double>(clock_.now_utc_ns() - e.tUtcNs) / 1e9;
    char dl[64];
    std::snprintf(dl, sizeof dl, "%.0f s (%.1f min)", delayS, delayS / 60.0);
    body << "Retraso de entrega: " << dl << "\n";
  }
  if (!e.origin.empty()) body << "Origen: " << e.origin << "\n";
  body << "Host: " << host_ << "\n";
  if (e.ch >= 0) {
    body << "Canal HV: CH" << (e.ch < 10 ? "0" : "") << e.ch << " del slot " << e.slot;
    if (group > 0) body << " (grupo HV " << group << ")";
    body << "\n";
    if (e.k > 0) body << "Intento: " << e.k << "\n";
    body << "\nCanales del digitizador alimentados por este canal HV:\n" << groupText(e.ch);
  }
  if (!e.reason.empty()) body << "\nMotivo: " << e.reason << "\n";
  if (e.timeline.is_array() && !e.timeline.empty()) {
    body << "\nSecuencia del episodio:\n";
    for (const auto& t : e.timeline) {
      const UtcNs tn = t.value("t_utc_ns", static_cast<UtcNs>(0));
      body << "  " << (tn != 0 ? utcAndLocal(tn) : t.value("t_utc", "")) << "  " << t.value("event", "")
           << " (k = " << t.value("k", 0) << ")";
      if (t.contains("wait_min") && t["wait_min"].is_number()) {
        char w[32];
        std::snprintf(w, sizeof w, "%g", t["wait_min"].get<double>());
        body << ", espera aplicada " << w << " min";
      }
      body << "\n";
    }
  }
  const auto& c = e.context;
  body << "\nÚltimo estado leído:\n";
  char num[64];
  if (c.contains("VMon_V") && c["VMon_V"].is_number()) {
    std::snprintf(num, sizeof num, "%.2f", c["VMon_V"].get<double>());
    body << "  VMon = " << num << " V\n";
  }
  if (c.contains("IMon_A") && c["IMon_A"].is_number()) {
    std::snprintf(num, sizeof num, "%.3f", c["IMon_A"].get<double>() * 1e6);
    body << "  IMon = " << num << " µA\n";
  }
  if (c.contains("Status")) {
    char hex[16];
    std::snprintf(hex, sizeof hex, "0x%04X", c["Status"].get<unsigned>());
    body << "  Status = " << hex << " " << (c.contains("status_flags") ? c["status_flags"].dump() : "") << "\n";
  }
  body << "  Temperatura de la placa A7030DP (no del detector): "
       << (c.contains("board_temperature_degC") && !c["board_temperature_degC"].is_null()
               ? c["board_temperature_degC"].dump() + " degC"
               : "sin dato")
       << "\n";
  body << "  GECO2020 en ejecución: " << (c.value("geco_running", false) ? "sí" : "no") << "\n";
  const std::string port = std::to_string(cfg_.pagePort);
  body << "\nPágina de estado (solo lectura): ";
  if (cfg_.pageSshHost.empty())
    body << "notifier.page_ssh_host sin configurar en hvmon.json.\n";
  else
    body << "ssh -N -L " << port << ":127.0.0.1:" << port << " " << cfg_.pageSshUser << "@" << cfg_.pageSshHost
         << " y abrir http://localhost:" << port << "\n";
  body << "\n-- \nhvmon (DAQ v2, hodoscopio). Mensaje automático; no responder.\n";

  // Encabezados y cuerpo con CRLF (SMTP).
  std::string to;
  for (const auto& a : smtp_.to) to += (to.empty() ? "" : ", ") + a;
  std::mt19937_64 rng(static_cast<std::uint64_t>(clock_.now_utc_ns()) ^ (seq_ * 0x9E3779B97F4A7C15ull));
  char msgid[128];
  std::snprintf(msgid, sizeof msgid, "<%lld.%016llx@%s>", static_cast<long long>(clock_.now_utc_ns()),
                static_cast<unsigned long long>(rng()), host_.c_str());
  std::ostringstream m;
  m << "Date: " << rfc5322Date(clock_.now_utc_ns()) << "\r\n"
    << "From: " << smtp_.from << "\r\n"
    << "To: " << to << "\r\n"
    << "Subject: " << headerText(subject) << "\r\n"
    << "Message-ID: " << msgid << "\r\n"
    << "MIME-Version: 1.0\r\n"
    << "Content-Type: text/plain; charset=UTF-8\r\n"
    << "Content-Transfer-Encoding: 8bit\r\n"
    << "\r\n";
  std::istringstream lines(body.str());
  std::string line;
  while (std::getline(lines, line)) m << (line == "." ? ".." : line) << "\r\n";
  return m.str();
}

std::string Notifier::curlConfig(const std::string& bodyPath) const {
  std::ostringstream c;
  c << "url = " << curlQuote(smtp_.url) << "\n"
    << "user = " << curlQuote(smtp_.user + ":" + smtp_.password) << "\n"
    << "mail-from = " << curlQuote(bareAddress(smtp_.from)) << "\n";
  for (const auto& a : smtp_.to) c << "mail-rcpt = " << curlQuote(bareAddress(a)) << "\n";
  c << "upload-file = " << curlQuote(bodyPath) << "\n"
    << "max-time = " << cfg_.maxTimeS << "\n"
    << "ssl-reqd\n"
    << "silent\n"
    << "show-error\n";
  return c.str();
}

bool Notifier::submit(NotifyEvent e) {
  auto drop = [&](const std::string& why) {
    ++suppressed_;
    logDebug("notificador: evento ", e.kind, " descartado (", why, ")");
    if (e.reply) e.reply->set_value({{"ok", false}, {"error", why}});
    publish();
    return false;
  };
  if (!cfg_.enabled && !e.test) return drop("notificador deshabilitado");
  if (!e.test) {
    if (e.kind == "trip" && !cfg_.notifyEachTrip) return drop("notify_each_trip = false");
    if (e.kind == "recovery" && !cfg_.notifyRecovery) return drop("notify_recovery = false");
    if (e.kind == "writer_error" && !cfg_.notifyWriterErrors) return drop("notify_writer_errors = false");
    if (e.kind == "connection_lost" && !cfg_.notifyConnectionLost) return drop("notify_connection_lost = false");
    if (cfg_.maxMessages > 0 && accepted_ >= static_cast<std::uint64_t>(cfg_.maxMessages))
      return drop("tope de max_messages = " + std::to_string(cfg_.maxMessages) + " alcanzado");
    if (e.kind == "trip" || e.kind == "recovery" || e.kind == "gave_up") {
      const std::string key = e.kind + ":" + std::to_string(e.ch) + ":" + std::to_string(e.episode);
      if (sentKeys_.count(key)) return drop("ya notificado en este episodio (" + key + ")");
      sentKeys_[key] = true;
    } else {
      const SteadyTime now = clock_.now_steady();
      auto it = lastProblem_.find(e.kind);
      const auto interval = std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::duration<double>(cfg_.problemIntervalMin * 60.0));
      if (it != lastProblem_.end() && now - it->second < interval) return drop("límite de 1 cada " +
                                                                               std::to_string(cfg_.problemIntervalMin) + " min");
      lastProblem_[e.kind] = now;
    }
  }
  if (!e.test) ++accepted_;
  ++seq_;
  Pending p;
  p.ev = std::move(e);
  p.queuedUtcNs = p.dueUtcNs = clock_.now_utc_ns();
  if (!p.ev.test && !cfg_.outboxDir.empty()) {
    PosixFileSystem().makeDirectories(cfg_.outboxDir, 0700);
    p.path = cfg_.outboxDir + "/" + std::to_string(p.queuedUtcNs) + "_" + std::to_string(::getpid()) + "_" +
             std::to_string(seq_) + ".json";
    persist(p);
  }
  pending_.push_back(std::move(p));
  publish();
  return true;
}

bool Notifier::sendOnce(Pending& p, std::string& err, int& exitCode) {
  // pid en el nombre: dos procesos (p. ej. tests en paralelo con el mismo reloj simulado) no chocan.
  const std::string path = tmpDir_ + "/mail_" + std::to_string(::getpid()) + "_" + std::to_string(clock_.now_utc_ns()) +
                           "_" + std::to_string(seq_) + ".eml";
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
  if (fd < 0) {
    err = "no se pudo crear " + path;
    return false;
  }
  const std::string message = composeMessage(p.ev);  // al enviar: incluye el retraso de entrega
  const char* d = message.data();
  std::size_t left = message.size();
  while (left > 0) {
    const ssize_t n = ::write(fd, d, left);
    if (n <= 0) break;
    d += n;
    left -= static_cast<std::size_t>(n);
  }
  ::close(fd);
  const std::vector<std::string> argv{cfg_.curlPath, "--config", "-"};
  auto r = spawner_.run(argv, curlConfig(path), std::chrono::seconds(cfg_.maxTimeS + 15));
  ::unlink(path.c_str());  // archivo temporal propio (sin secretos)
  if (!r) {
    err = redact(r.error().message);
    exitCode = -1;
    return false;
  }
  exitCode = r.value().exitCode;
  if (r.value().exitCode == 0 && !r.value().timedOut) return true;
  err = "curl salió con " + std::to_string(r.value().exitCode) + (r.value().timedOut ? " (timeout)" : "") + ": " +
        redact(r.value().stderrText);
  return false;
}

void Notifier::processDue() {
  const UtcNs now = clock_.now_utc_ns();
  for (auto it = pending_.begin(); it != pending_.end();) {
    if (it->dueUtcNs > now) {
      ++it;
      continue;
    }
    std::string err;
    int code = -1;
    ++it->attempts;
    const bool ok = sendOnce(*it, err, code);
    lastResult_ = {{"kind", it->ev.kind}, {"ch", it->ev.ch}, {"test", it->ev.test}, {"attempt", it->attempts},
                   {"ok", ok}, {"curl_exit", code}, {"error", err}, {"t_utc", formatUtcIso(clock_.now_utc_ns())},
                   {"delay_s", static_cast<double>(clock_.now_utc_ns() - it->ev.tUtcNs) / 1e9}};
    if (ok) {
      ++sent_;
      logInfo("correo enviado: ", it->ev.kind, " CH", it->ev.ch, it->ev.test ? " (PRUEBA)" : "", ", curl salió con 0",
              ", intento ", it->attempts, ", retraso ", static_cast<double>(clock_.now_utc_ns() - it->ev.tUtcNs) / 1e9, " s");
      if (it->ev.reply) it->ev.reply->set_value({{"ok", true}, {"curl_exit", code}, {"attempts", it->attempts}});
      if (!it->path.empty()) ::unlink(it->path.c_str());  // entregado: sale del outbox
      it = pending_.erase(it);
      continue;
    }
    const double ageH = static_cast<double>(now - it->ev.tUtcNs) / 3.6e12;
    if (!it->ev.test && ageH < cfg_.maxAgeH) {
      const std::size_t i = static_cast<std::size_t>(it->attempts - 1);
      const double delay = i < cfg_.retryDelaysMin.size() ? cfg_.retryDelaysMin[i] : cfg_.retryHourlyMin;
      logWarn("correo ", it->ev.kind, " CH", it->ev.ch, " falló (intento ", it->attempts, "): ", err,
              "; reintento en ", delay, " min");
      it->dueUtcNs = now + static_cast<UtcNs>(delay * 60e9);
      persist(*it);
      ++it;
      continue;
    }
    ++failed_;
    if (!it->ev.test) ++abandoned_;
    logError("correo ", it->ev.kind, " CH", it->ev.ch, " no se pudo enviar tras ", it->attempts, " intentos",
             it->ev.test ? "" : " (más de max_age_h desde el evento)", ": ", err);
    if (it->ev.reply) it->ev.reply->set_value({{"ok", false}, {"curl_exit", code}, {"error", err}, {"attempts", it->attempts}});
    if (!it->path.empty()) {  // abandonado: se aparta, no se borra
      const std::string failedDir = cfg_.outboxDir + "/failed";
      PosixFileSystem().makeDirectories(failedDir, 0700);
      ::rename(it->path.c_str(), (failedDir + it->path.substr(it->path.find_last_of('/'))).c_str());
    }
    it = pending_.erase(it);
  }
  publish();
}

void Notifier::checkProblems() {
  const auto w = hub_.writer.read().first;
  if (w.errors > lastWriterErrors_) {
    lastWriterErrors_ = w.errors;
    NotifyEvent e;
    e.kind = "writer_error";
    e.reason = "errores de escritura: " + std::to_string(w.errors) + "; último: " + w.lastError;
    e.tUtcNs = w.lastErrorUtcNs ? w.lastErrorUtcNs : clock_.now_utc_ns();
    submit(std::move(e));
  }
  const auto p = hub_.poller.read().first;
  const UtcNs now = clock_.now_utc_ns();
  const UtcNs ref = p.lastGoodReadUtcNs ? p.lastGoodReadUtcNs : hub_.startedUtcNs;
  const bool lost = !p.connected && ref > 0 &&
                    static_cast<double>(now - ref) > cfg_.connectionLostMin * 60e9;
  if (lost && !connectionAlerted_) {
    connectionAlerted_ = true;
    NotifyEvent e;
    e.kind = "connection_lost";
    e.reason = "sin lecturas del SY5527 desde " + formatUtcIso(ref) + "; último error: " + p.lastError;
    e.tUtcNs = now;
    submit(std::move(e));
  }
  if (p.connected) connectionAlerted_ = false;
}

nlohmann::json Notifier::status() const {
  return nlohmann::json{{"enabled", cfg_.enabled}, {"sent", sent_},          {"failed", failed_},
                        {"restored_from_outbox", restored_}, {"abandoned", abandoned_}, {"outbox", cfg_.outboxDir},
                        {"suppressed", suppressed_}, {"accepted", accepted_}, {"max_messages", cfg_.maxMessages}, {"pending", pending_.size()}, {"last_result", lastResult_},
                        {"recipients", smtp_.to.size()}};
}

void Notifier::publish() { hub_.notifier.publish(status()); }

void Notifier::run(BoundedQueue<NotifyEvent>& queue, const StopToken& st) {
  Logger::setThreadName("notify");
  SteadyTime nextCheck = clock_.now_steady();
  while (!st.stop_requested()) {
    auto e = queue.pop_until(std::chrono::steady_clock::now() + std::chrono::milliseconds(200));
    if (e) submit(std::move(*e));
    if (clock_.now_steady() >= nextCheck) {
      checkProblems();
      nextCheck = clock_.now_steady() + std::chrono::seconds(10);
    }
    processDue();
    if (queue.closed() && !e && pending_.empty()) break;
  }
  for (auto& p : pending_)
    if (p.ev.reply) p.ev.reply->set_value({{"ok", false}, {"error", "hvmon se está apagando"}});
  if (!pending_.empty()) logWarn("notificador: ", pending_.size(), " correos sin enviar al apagar");
}

}  // namespace daq::hvmon
