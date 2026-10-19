// sc_reader.cpp
#include "daq/io/sc_reader.h"

#include <sys/stat.h>

#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/gzip.h"

namespace daq::io {

namespace {

bool exists(const std::string& p) {
  struct stat st {};
  return ::stat(p.c_str(), &st) == 0;
}

std::vector<std::string> splitCsv(const std::string& s) {
  std::vector<std::string> out;
  std::size_t a = 0;
  for (;;) {
    const auto b = s.find(',', a);
    out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
    if (b == std::string::npos) break;
    a = b + 1;
  }
  return out;
}

// Días YYYYMMDD de from a to.
std::vector<std::string> days(const std::string& from, const std::string& to) {
  std::vector<std::string> out;
  std::tm tm{};
  tm.tm_year = std::atoi(from.substr(0, 4).c_str()) - 1900;
  tm.tm_mon = std::atoi(from.substr(4, 2).c_str()) - 1;
  tm.tm_mday = std::atoi(from.substr(6, 2).c_str());
  std::time_t t = timegm(&tm);
  for (int i = 0; i < 3660; ++i) {
    std::tm d{};
    gmtime_r(&t, &d);
    char b[16];
    std::snprintf(b, sizeof b, "%04d%02d%02d", (d.tm_year + 1900) % 10000, (d.tm_mon + 1) % 100, d.tm_mday % 100);
    if (std::string(b) > to) break;
    out.push_back(b);
    t += 86400;
  }
  return out;
}

struct Part {
  std::string path;  // archivo de datos (plano o .gz)
  std::string schemaPath;
  int part = 1;
};

// Partes existentes de un día: <dir>/YYYY/<kind>_<day>[.N].<ext>[.gz].
std::vector<Part> partsOf(const std::string& dir, const std::string& kind, const std::string& day, const std::string& ext) {
  std::vector<Part> out;
  const std::string base = dir + "/" + day.substr(0, 4) + "/" + kind + "_" + day;
  for (int n = 1; n < 100; ++n) {
    const std::string stem = n == 1 ? base : base + "." + std::to_string(n);
    const std::string plain = stem + "." + ext;
    std::string found;
    if (exists(plain)) found = plain;
    else if (exists(plain + ".gz")) found = plain + ".gz";
    if (found.empty()) {
      if (n > 1) break;
      continue;
    }
    out.push_back(Part{found, plain + ".schema.json", n});
  }
  return out;
}

char leafType(const std::string& t) {
  if (t == "int64") return 'L';
  if (t == "uint64") return 'l';
  if (t == "int") return 'I';
  if (t == "uint32") return 'i';
  return 'D';
}

struct Column {
  char type = 'D';
  double d = 0;
  Long64_t L = 0;
  ULong64_t l = 0;
  Int_t I = 0;
  UInt_t i = 0;
  void setMissing() {
    d = std::numeric_limits<double>::quiet_NaN();
    L = 0;
    l = 0;
    I = 0;
    i = 0;
  }
  void set(const std::string& s) {
    switch (type) {
      case 'L': L = std::strtoll(s.c_str(), nullptr, 10); break;
      case 'l': l = std::strtoull(s.c_str(), nullptr, 10); break;
      case 'I': I = static_cast<Int_t>(std::strtol(s.c_str(), nullptr, 10)); break;
      case 'i': i = static_cast<UInt_t>(std::strtoul(s.c_str(), nullptr, 10)); break;
      default: d = s.empty() ? std::numeric_limits<double>::quiet_NaN() : std::strtod(s.c_str(), nullptr);
    }
  }
  void* addr() {
    switch (type) {
      case 'L': return &L;
      case 'l': return &l;
      case 'I': return &I;
      case 'i': return &i;
      default: return &d;
    }
  }
};

struct OpenPart {
  Part p;
  std::unique_ptr<LineReader> lr;
  std::vector<std::string> header;
  std::vector<int> colIndex;  // columna de la parte -> columna del árbol
  int keyCol = -1;
  std::optional<std::string> line;
  long long key = 0;
  bool advance() {
    for (;;) {
      line = lr->next();
      if (!line) return false;
      if (line->empty()) continue;
      if (keyCol < 0) return true;
      const auto f = splitCsv(*line);
      key = keyCol < static_cast<int>(f.size()) ? std::atoll(f[static_cast<std::size_t>(keyCol)].c_str()) : 0;
      return true;
    }
  }
};

TTree* readCsv(const std::string& dir, const std::string& kind, const std::vector<std::string>& ds, std::string* err) {
  // Columnas: unión en orden de aparición; tipo según el esquema de la primera parte que la tiene.
  std::vector<std::string> names;
  std::map<std::string, char> types;
  std::vector<std::vector<Part>> perDay;
  for (const auto& d : ds) {
    perDay.push_back(partsOf(dir, kind, d, "csv"));
    for (const auto& p : perDay.back()) {
      auto lr = LineReader::open(p.path);
      if (!lr) continue;
      auto h = lr.value()->next();
      if (!h) continue;
      std::map<std::string, std::string> schemaTypes;
      if (auto sl = LineReader::open(p.schemaPath)) {
        std::string text;
        while (auto l = sl.value()->next()) text += *l + "\n";
        try {
          const nlohmann::json sj = nlohmann::json::parse(text);  // vive durante el bucle (no un temporal)
          for (const auto& c : sj.at("columns")) schemaTypes[c.at("name").get<std::string>()] = c.value("type", "float64");
        } catch (const std::exception&) {
        }
      }
      for (const auto& n : splitCsv(*h)) {
        if (types.count(n)) continue;
        names.push_back(n);
        types[n] = leafType(schemaTypes.count(n) ? schemaTypes[n] : "float64");
      }
    }
  }
  if (names.empty()) {
    if (err) *err = "no hay archivos " + kind + " en " + dir + " para " + ds.front() + "-" + ds.back();
    return nullptr;
  }
  auto* tree = new TTree(kind.c_str(), ("control lento " + kind + " (daqv2)").c_str());
  tree->SetDirectory(nullptr);
  std::vector<Column> cols(names.size());
  Int_t part = 0;
  for (std::size_t c = 0; c < names.size(); ++c) {
    cols[c].type = types[names[c]];
    tree->Branch(names[c].c_str(), cols[c].addr(), (names[c] + "/" + cols[c].type).c_str());
  }
  tree->Branch("part", &part, "part/I");
  for (const auto& parts : perDay) {
    std::vector<OpenPart> open;
    for (const auto& p : parts) {
      auto lr = LineReader::open(p.path);
      if (!lr) continue;
      OpenPart op;
      op.p = p;
      op.lr = std::move(lr).value();
      auto h = op.lr->next();
      if (!h) continue;
      op.header = splitCsv(*h);
      for (std::size_t i = 0; i < op.header.size(); ++i) {
        int idx = -1;
        for (std::size_t c = 0; c < names.size(); ++c)
          if (names[c] == op.header[i]) idx = static_cast<int>(c);
        op.colIndex.push_back(idx);
        if (op.header[i] == "t_resp_utc_ns" || (op.keyCol < 0 && op.header[i] == "t_utc_ns")) op.keyCol = static_cast<int>(i);
      }
      if (op.advance()) open.push_back(std::move(op));
    }
    // Intercalado por la marca de tiempo.
    while (!open.empty()) {
      std::size_t best = 0;
      for (std::size_t k = 1; k < open.size(); ++k)
        if (open[k].key < open[best].key) best = k;
      OpenPart& op = open[best];
      for (auto& c : cols) c.setMissing();
      const auto f = splitCsv(*op.line);
      for (std::size_t i = 0; i < f.size() && i < op.colIndex.size(); ++i)
        if (op.colIndex[i] >= 0) cols[static_cast<std::size_t>(op.colIndex[i])].set(f[i]);
      part = op.p.part;
      tree->Fill();
      if (!op.advance()) open.erase(open.begin() + static_cast<std::ptrdiff_t>(best));
    }
  }
  return tree;
}

TTree* readJsonl(const std::string& dir, const std::string& kind, const std::vector<std::string>& ds, std::string* err) {
  auto* tree = new TTree(kind.c_str(), ("control lento " + kind + " (daqv2)").c_str());
  tree->SetDirectory(nullptr);
  Long64_t t = 0, tPrev = 0, id = 0, cmd = 0;
  Int_t ch = -1, part = 0;
  std::string type, param, oldV, newV, origin, user, host, text;
  tree->Branch("t_utc_ns", &t, "t_utc_ns/L");
  tree->Branch("part", &part, "part/I");
  if (kind == "notes") {
    tree->Branch("user", &user);
    tree->Branch("host", &host);
    tree->Branch("text", &text);
  } else {
    tree->Branch("t_prev_read_utc_ns", &tPrev, "t_prev_read_utc_ns/L");
    tree->Branch("id", &id, "id/L");
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("command_id", &cmd, "command_id/L");
    tree->Branch("type", &type);
    tree->Branch("param", &param);
    tree->Branch("old_value", &oldV);
    tree->Branch("new_value", &newV);
    tree->Branch("origin", &origin);
  }
  bool any = false;
  for (const auto& d : ds)
    for (const auto& p : partsOf(dir, kind, d, "jsonl")) {
      auto lr = LineReader::open(p.path);
      if (!lr) continue;
      any = true;
      while (auto l = lr.value()->next()) {
        try {
          const auto j = nlohmann::json::parse(*l);
          t = j.value("t_utc_ns", Long64_t{0});
          part = p.part;
          if (kind == "notes") {
            user = j.value("user", "");
            host = j.value("host", "");
            text = j.value("text", "");
          } else {
            tPrev = j.value("t_prev_read_utc_ns", Long64_t{0});
            id = j.contains("id") && j["id"].is_number() ? j["id"].get<Long64_t>() : -1;
            ch = j.value("ch", -1);
            cmd = j.contains("command_id") && j["command_id"].is_number() ? j["command_id"].get<Long64_t>() : -1;
            type = j.value("type", "change");
            param = j.value("param", "");
            oldV = j.contains("old") ? j["old"].dump() : "";
            newV = j.contains("new") ? j["new"].dump() : "";
            origin = j.value("origin", "");
          }
          tree->Fill();
        } catch (const std::exception&) {
        }
      }
    }
  if (!any) {
    delete tree;
    if (err) *err = "no hay archivos " + kind + " en " + dir + " para " + ds.front() + "-" + ds.back();
    return nullptr;
  }
  return tree;
}

}  // namespace

TTree* readSlowControl(const std::string& dir, const std::string& kind, const std::string& fromDate,
                       const std::string& toDate, std::string* error) {
  if (fromDate.size() != 8 || toDate.size() != 8 || toDate < fromDate) {
    if (error) *error = "fechas YYYYMMDD inválidas";
    return nullptr;
  }
  const auto ds = days(fromDate, toDate);
  if (kind == "hv" || kind == "hv_full" || kind == "temp") return readCsv(dir, kind, ds, error);
  if (kind == "hv_changes" || kind == "notes") return readJsonl(dir, kind, ds, error);
  if (error) *error = "tipo desconocido: " + kind;
  return nullptr;
}

}  // namespace daq::io
