// http_page.cpp
#include "http_page.h"

#include <cstdio>
#include <sstream>

#include <TAxis.h>
#include <TGraph.h>
#include <THttpServer.h>
#include <TNamed.h>

#include "daq/hv/reading_json.h"

namespace daq::hvmon {

void HttpPage::Ring::push(double tt, double yy) {
  t[head] = tt;
  y[head] = yy;
  head = (head + 1) % t.size();
  if (n < t.size()) ++n;
}

Result<std::unique_ptr<HttpPage>> HttpPage::create(const HvmonConfig& cfg, const std::string& bindAddress) {
  using R = Result<std::unique_ptr<HttpPage>>;
  std::unique_ptr<HttpPage> page(new HttpPage());
  page->cfg_ = &cfg;
  const std::string engine = "http:" + std::to_string(cfg.httpPort) + "?bind=" + bindAddress + "&thrds=2";
  page->server_ = std::make_unique<THttpServer>(engine.c_str());
  if (!page->server_->IsAnyEngine()) return R::failure(-1, "THttpServer no pudo iniciar " + engine);
  page->server_->SetReadOnly(kTRUE);
  page->server_->SetTimer(0, kTRUE);  // solo atiende en processRequests() (hilo principal)
  page->url_ = "http://" + bindAddress + ":" + std::to_string(cfg.httpPort) + "/";
  const auto capacity = static_cast<std::size_t>(cfg.graphHours * 3600.0 / cfg.graphPeriodS) + 1;
  for (int ch : cfg.httpChannels) {
    char name[32];
    Channel c;
    c.vRing.t.assign(capacity, 0.0);
    c.vRing.y.assign(capacity, 0.0);
    c.iRing = c.vRing;
    std::snprintf(name, sizeof name, "CH%02d", ch);
    const std::string folder = std::string("/") + name;
    auto mk = [&](const char* what, const char* unit) {
      auto g = std::make_unique<TGraph>();
      const std::string gname = std::string(what) + "_" + name;
      g->SetName(gname.c_str());
      g->SetTitle((gname + ";tiempo (UTC);" + what + " [" + unit + "]").c_str());
      g->GetXaxis()->SetTimeDisplay(1);
      g->GetXaxis()->SetTimeFormat("%H:%M%F1970-01-01 00:00:00s0 GMT");
      TGraph* raw = g.get();
      page->server_->Register(folder.c_str(), raw);
      page->graphs_.push_back(std::move(g));
      return raw;
    };
    c.vmon = mk("VMon", "V");
    c.imon = mk("IMon", "A");
    page->channels_[ch] = std::move(c);
  }
  {
    auto g = std::make_unique<TGraph>();
    g->SetName("BoardTemp_degC");
    g->SetTitle("Temperatura de la placa A7030DP (no del detector);tiempo (UTC);T [degC]");
    g->GetXaxis()->SetTimeDisplay(1);
    g->GetXaxis()->SetTimeFormat("%H:%M%F1970-01-01 00:00:00s0 GMT");
    page->boardTemp_ = g.get();
    page->boardRing_.t.assign(static_cast<std::size_t>(cfg.graphHours * 3600.0 / 60.0) + 1, 0.0);
    page->boardRing_.y = page->boardRing_.t;
    page->server_->Register("/board", g.get());
    page->graphs_.push_back(std::move(g));
  }
  page->table_ = std::make_unique<TNamed>("table", "sin datos");
  page->json_ = std::make_unique<TNamed>("state_json", "{}");
  page->server_->Register("/status", page->table_.get());
  page->server_->Register("/status", page->json_.get());
  return page;
}

HttpPage::~HttpPage() {
  if (server_) {
    for (auto& g : graphs_) server_->Unregister(g.get());
    if (table_) server_->Unregister(table_.get());
    if (json_) server_->Unregister(json_.get());
  }
  server_.reset();
}

void HttpPage::refill(TGraph* g, const Ring& r) {
  g->Set(static_cast<Int_t>(r.n));
  const std::size_t cap = r.t.size();
  const std::size_t first = (r.head + cap - r.n) % cap;
  for (std::size_t i = 0; i < r.n; ++i) {
    const std::size_t k = (first + i) % cap;
    g->SetPoint(static_cast<Int_t>(i), r.t[k], r.y[k]);
  }
}

void HttpPage::update(const PollerState& p, const nlohmann::json& state, UtcNs now) {
  const double periodNs = cfg_->graphPeriodS * 1e9;
  if (p.lastGoodReadUtcNs > 0 && static_cast<double>(now - lastGraphUtc_) >= periodNs) {
    lastGraphUtc_ = now;
    const double t = static_cast<double>(p.lastGoodReadUtcNs) / 1e9;
    for (auto& kv : channels_) {
      auto it = p.values.find(kv.first);
      if (it == p.values.end()) continue;
      auto v = it->second.find("VMon");
      auto i = it->second.find("IMon");
      if (v != it->second.end()) kv.second.vRing.push(t, v->second.si);
      if (i != it->second.end()) kv.second.iRing.push(t, i->second.si);
      refill(kv.second.vmon, kv.second.vRing);
      refill(kv.second.imon, kv.second.iRing);
    }
  }
  if (p.board && p.board->tRespUtcNs != lastBoardUtc_) {  // un punto por lectura completa (60 s)
    lastBoardUtc_ = p.board->tRespUtcNs;
    if (const auto* t = p.board->get("Temp")) {
      boardRing_.push(static_cast<double>(p.board->tRespUtcNs) / 1e9, t->si);
      refill(boardTemp_, boardRing_);
    }
  }
  std::ostringstream os;
  const auto& conn = state["connection"];
  os << "hvmon  " << state.value("t_utc", "") << "\n"
     << "conectado=" << conn.value("connected", false) << "  stale=" << conn.value("stale", true)
     << "  edad=" << conn.value("age_s", -1.0) << " s  GECO=" << state.value("geco_running", false)
     << "  allow_commands=" << state.value("allow_commands", false) << "\n"
     << "ultimo error: " << conn.value("last_error", "") << "\n";
  if (state["board_temperature_degC"].is_number())
    os << "temperatura de la placa A7030DP (board temperature, no del detector): "
       << state["board_temperature_degC"].get<double>() << " degC\n";
  os << "\n";
  os << "CH  V0Set[V]  I0Set[A]     VMon[V]  IMon[A]      Pw  Status\n";
  for (const auto& c : state["channels"]) {
    const auto& v = c["values"];
    char line[256];
    std::string flags;
    for (const auto& f : c["status_flags"]) flags += (flags.empty() ? "" : ",") + f.get<std::string>();
    std::snprintf(line, sizeof line, "%02d  %8.3f  %10.3e  %7.3f  %10.3e  %2u  %s\n", c.value("ch", -1),
                  v.value("V0Set_V", 0.0), v.value("I0Set_A", 0.0), v.value("VMon_V", 0.0), v.value("IMon_A", 0.0),
                  v.value("Pw", 0u), flags.empty() ? "-" : flags.c_str());
    os << line;
  }
  os << "\nAuto-retry\nCH  estado          k  trip (UTC)                       prox. intento  ultimo resultado\n";
  for (const auto& c : state["channels"]) {
    const auto& r = c["auto_retry"];
    char line[512];
    const std::string trip = r.contains("trip_utc") && r["trip_utc"].is_string() ? r["trip_utc"].get<std::string>() : "-";
    const std::string next = r.contains("next_retry_in_s") && r["next_retry_in_s"].is_number()
                                 ? std::to_string(static_cast<long>(r["next_retry_in_s"].get<double>())) + " s"
                                 : "-";
    const std::string block = r.value("block_reason", "").empty() ? "" : " [" + r.value("block_reason", "") + "]";
    std::snprintf(line, sizeof line, "%02d  %-14s %2d  %-32s %-13s  %s%s\n", c.value("ch", -1),
                  r.value("label", r.value("state", "?")).c_str(), r.value("k", 0), trip.c_str(), next.c_str(),
                  r.value("last_outcome", "").c_str(), block.c_str());
    os << line;
  }
  if (state.contains("notifier") && state["notifier"].is_object()) {
    const auto& n = state["notifier"];
    os << "\nnotificador: habilitado=" << n.value("enabled", false) << " enviados=" << n.value("sent", 0)
       << " fallidos=" << n.value("failed", 0) << " pendientes=" << n.value("pending", 0) << "\n";
  }
  os << "\nultimos cambios:\n";
  const auto& ch = state["changes"];
  const std::size_t n = ch.size();
  for (std::size_t k = n > 10 ? n - 10 : 0; k < n; ++k) {
    const auto& c = ch[k];
    os << "  " << c.value("t_utc", "") << " CH" << c.value("ch", -1) << " " << c.value("param", "") << " "
       << c["old"].dump() << " -> " << c["new"].dump() << " (" << c.value("origin", "") << ")\n";
  }
  table_->SetTitle(os.str().c_str());
  json_->SetTitle(state.dump().c_str());
}

void HttpPage::processRequests() {
  if (server_) server_->ProcessRequests();
}

}  // namespace daq::hvmon
