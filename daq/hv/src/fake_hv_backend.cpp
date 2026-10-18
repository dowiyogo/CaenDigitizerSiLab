// fake_hv_backend.cpp
#include "daq/hv/fake_hv_backend.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace daq::hv {

namespace {

// Bits de Status (convención del manual de la Wrapper, Rev. 13).
constexpr std::uint32_t kOn = 1u << 0;
constexpr std::uint32_t kRampUp = 1u << 1;
constexpr std::uint32_t kRampDown = 1u << 2;
constexpr std::uint32_t kOvercurrent = 1u << 3;
constexpr std::uint32_t kInternalTrip = 1u << 9;

constexpr double kMaxStep = 0.01;  // s, paso máximo de integración
constexpr double kEps = 1e-9;

ParamInfo numeric(const char* name, ParamMode mode, int unit, int exp, double mn, double mx) {
  ParamInfo p;
  p.name = name;
  p.type = ParamType::Numeric;
  p.mode = mode;
  p.unitCode = unit;
  p.exponent = exp;
  p.unitFromBoard = true;
  p.minBoard = mn;
  p.maxBoard = mx;
  return p;
}

ParamInfo other(const char* name, ParamType type, ParamMode mode, const char* on = "", const char* off = "") {
  ParamInfo p;
  p.name = name;
  p.type = type;
  p.mode = mode;
  p.onState = on;
  p.offState = off;
  return p;
}

// Valor escalar o arreglo por canal.
double pick(const nlohmann::json& v, std::size_t ch, const std::string& what) {
  if (v.is_number()) return v.get<double>();
  if (v.is_array()) {
    if (ch >= v.size()) throw std::invalid_argument("guion: " + what + " tiene menos valores que canales");
    return v.at(ch).get<double>();
  }
  throw std::invalid_argument("guion: " + what + " debe ser número o arreglo");
}

double quantize(double x, double step) { return step > 0 ? std::round(x / step) * step : x; }

}  // namespace

std::vector<ParamInfo> FakeHVBackend::a7030dpParams() {
  // Copia de lo que informa la placa real (hvctl params --slot 3, canales 0 y 8,
  // 2026-10-18T22:09Z): tipos, modos, unidad, exponente, rangos y textos ON/OFF.
  const auto RW = ParamMode::ReadWrite;
  const auto RO = ParamMode::ReadOnly;
  return {
      numeric("V0Set", RW, 2, 0, 0.0, 3000.0),
      numeric("I0Set", RW, 1, -6, 0.0, 1010.0),
      numeric("V1Set", RW, 2, 0, 0.0, 3000.0),
      numeric("I1Set", RW, 1, -6, 0.0, 1010.0),
      numeric("RUp", RW, 7, 0, 1.0, 500.0),
      numeric("RDWn", RW, 7, 0, 1.0, 500.0),
      numeric("Trip", RW, 8, 0, 0.0, 1000.0),
      numeric("SVMax", RW, 2, 0, 0.0, 3000.0),
      numeric("VMon", RO, 2, 0, 0.0, 3000.0),
      numeric("IMon", RO, 1, -6, 0.0, 1000.0),
      other("Status", ParamType::ChStatus, RO),
      other("Pw", ParamType::OnOff, RW, "On", "Off"),
      other("POn", ParamType::OnOff, RW, "En", "Dis"),
      other("PDwn", ParamType::OnOff, RW, "Ramp", "Kill"),
      other("TripInt", ParamType::Binary, RW),
      other("TripExt", ParamType::Binary, RW),
      other("ZCDetect", ParamType::OnOff, RW, "On", "Off"),
      other("ZCAdjust", ParamType::OnOff, RW, "En", "Dis"),
  };
}

nlohmann::json FakeHVBackend::defaultScenario() {
  // Valores de hv_readall en el slot 3 el 2026-10-18T21:56:18Z, HV apagada (CH00-11).
  return nlohmann::json::parse(R"json({
    "seed": 20261018,
    "crate_slots": 4,
    "slots": [{
      "slot": 3, "model": "A7030DP", "description": " 12 Ch Pos 3KV 1mA (1.5W) ", "serial": 58015,
      "fw": [1, 15], "channels": 12,
      "set": {
        "V0Set": [65.75, 65.75, 65.90, 65.85, 65.75, 65.70, 66.15, 65.90, 0, 0, 0, 0],
        "I0Set": [15, 15, 15, 15, 15, 15, 15, 15, 5, 5, 5, 5],
        "V1Set": 0,
        "I1Set": [20, 20, 20, 20, 20, 20, 20, 20, 101, 101, 101, 101],
        "RUp":   [1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5],
        "RDWn": 50,
        "Trip":  [5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10],
        "SVMax": 70,
        "Pw": 0, "POn": 0, "PDwn": 0, "TripInt": 0, "TripExt": 16, "ZCDetect": 0, "ZCAdjust": 0
      },
      "load_uA": [2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0],
      "imon_offset_uA": [-0.048, -0.036, -0.046, -0.036, -0.044, -0.046, -0.046, -0.074, -0.118, -0.094, -0.074, -0.084],
      "vmon_offset_V": [0.19, 0.92, 0.67, 0.78, 0.38, 0.0, 0.0, 0.62, 1.26, 1.18, 1.15, 1.47],
      "imon_noise_uA": 0.004, "vmon_noise_V": 0.02,
      "imon_step_uA": 0.002, "vmon_step_V": 0.01
    }],
    "events": []
  })json");
}

FakeHVBackend::FakeHVBackend(const nlohmann::json& sc, const Clock& clock)
    : clock_(clock), t0_(clock.now_steady()), params_(a7030dpParams()) {
  rng_.seed(sc.value("seed", static_cast<std::uint64_t>(1)));
  if (!sc.contains("slots") || !sc["slots"].is_array() || sc["slots"].empty())
    throw std::invalid_argument("guion: falta 'slots'");
  int maxSlot = 0;
  for (const auto& js : sc["slots"]) {
    Slot s;
    s.info.slot = js.at("slot").get<int>();
    s.info.model = js.value("model", "A7030DP");
    s.info.description = js.value("description", "");
    s.info.serial = js.value("serial", 0u);
    if (js.contains("fw")) {
      s.info.fwMajor = js["fw"].at(0).get<int>();
      s.info.fwMinor = js["fw"].at(1).get<int>();
    }
    s.info.channels = js.value("channels", 12);
    if (s.info.slot < 0 || s.info.channels <= 0) throw std::invalid_argument("guion: slot o channels inválido");
    maxSlot = std::max(maxSlot, s.info.slot);
    s.imonNoiseUA = js.value("imon_noise_uA", 0.0);
    s.vmonNoiseV = js.value("vmon_noise_V", 0.0);
    s.imonStepUA = js.value("imon_step_uA", 0.002);
    s.vmonStepV = js.value("vmon_step_V", 0.01);
    const auto& set = js.contains("set") ? js["set"] : nlohmann::json::object();
    for (int c = 0; c < s.info.channels; ++c) {
      const auto ch = static_cast<std::size_t>(c);
      Channel k;
      for (const auto& p : params_) {
        if (p.name == "VMon" || p.name == "IMon" || p.name == "Status") continue;
        const double v = set.contains(p.name) ? pick(set[p.name], ch, p.name) : 0.0;
        if (p.type == ParamType::Numeric)
          k.num[p.name] = v;
        else
          k.ints[p.name] = static_cast<std::uint32_t>(v);
      }
      if (js.contains("load_uA")) k.loadUA = pick(js["load_uA"], ch, "load_uA");
      if (js.contains("imon_offset_uA")) k.imonOffsetUA = pick(js["imon_offset_uA"], ch, "imon_offset_uA");
      if (js.contains("vmon_offset_V")) k.vmonOffsetV = pick(js["vmon_offset_V"], ch, "vmon_offset_V");
      if (js.contains("initial_v_true")) k.vTrue = pick(js["initial_v_true"], ch, "initial_v_true");
      else if (k.ints["Pw"] == 1) k.vTrue = k.num["V0Set"];  // encendido y estable al inicio
      s.channels.push_back(k);
    }
    slots_.push_back(std::move(s));
  }
  crateSlots_ = std::max(sc.value("crate_slots", maxSlot + 1), maxSlot + 1);

  for (const auto& je : sc.value("events", nlohmann::json::array())) {
    Event e;
    e.t = je.at("t_s").get<double>();
    e.type = je.at("type").get<std::string>();
    e.slot = je.value("slot", slots_.front().info.slot);
    e.channel = je.value("ch", -1);
    e.duration = je.value("duration_s", 0.0);
    if (e.type == "set") {
      e.param = je.at("param").get<std::string>();
      e.value = je.at("value").get<double>();
      events_.push_back(e);
    } else if (e.type == "current") {
      e.value = je.at("extra_uA").get<double>();
      events_.push_back(e);
      Event end = e;
      end.t = e.t + e.duration;
      end.value = -e.value;
      events_.push_back(end);
    } else if (e.type == "connection_loss") {
      losses_.emplace_back(e.t, e.t + e.duration);
    } else {
      throw std::invalid_argument("guion: tipo de evento desconocido: " + e.type);
    }
    if (e.type != "connection_loss" && (e.channel < 0 || !findSlot(e.slot) ||
                                        e.channel >= findSlot(e.slot)->info.channels))
      throw std::invalid_argument("guion: evento con slot/canal inválido");
  }
  std::stable_sort(events_.begin(), events_.end(), [](const Event& a, const Event& b) { return a.t < b.t; });
}

Result<std::unique_ptr<FakeHVBackend>> FakeHVBackend::fromFile(const std::string& path, const Clock& clock) {
  using R = Result<std::unique_ptr<FakeHVBackend>>;
  try {
    if (path == "default") return std::make_unique<FakeHVBackend>(defaultScenario(), clock);
    std::ifstream in(path);
    if (!in) return R::failure(-1, "no se pudo abrir el guion " + path);
    nlohmann::json j;
    in >> j;
    // "base": "default" toma los slots del guion por defecto; el archivo aporta eventos y semilla.
    if (j.value("base", "") == "default") {
      auto merged = defaultScenario();
      for (auto it = j.begin(); it != j.end(); ++it)
        if (it.key() != "base") merged[it.key()] = it.value();
      j = merged;
    }
    return std::make_unique<FakeHVBackend>(j, clock);
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("guion inválido ") + path + ": " + e.what());
  }
}

double FakeHVBackend::seconds(SteadyTime t) const {
  return std::chrono::duration<double>(t - t0_).count();
}

FakeHVBackend::Slot* FakeHVBackend::findSlot(int slot) {
  for (auto& s : slots_)
    if (s.info.slot == slot) return &s;
  return nullptr;
}

bool FakeHVBackend::inConnectionLoss(double t) const {
  return std::any_of(losses_.begin(), losses_.end(), [t](const auto& l) { return t >= l.first && t < l.second; });
}

void FakeHVBackend::apply(const Event& e) {
  Slot* s = findSlot(e.slot);
  Channel& c = s->channels[static_cast<std::size_t>(e.channel)];
  if (e.type == "current") {
    c.extraUA += e.value;
    if (std::fabs(c.extraUA) < kEps) c.extraUA = 0.0;
    return;
  }
  // "set": cambio externo (p. ej. hecho en GECO).
  if (c.num.count(e.param)) {
    c.num[e.param] = e.value;
  } else if (c.ints.count(e.param)) {
    const auto v = static_cast<std::uint32_t>(e.value);
    if (e.param == "Pw" && v == 1u) c.tripped = false;  // encender borra el trip
    c.ints[e.param] = v;
  }
}

void FakeHVBackend::integrate(double dt) {
  for (auto& s : slots_) {
    for (auto& c : s.channels) {
      const bool pw = c.ints["Pw"] == 1u;
      const double v0 = c.num["V0Set"];
      const double target = pw ? v0 : 0.0;
      if (c.vTrue < target) {
        c.vTrue = std::min(target, c.vTrue + c.num["RUp"] * dt);
      } else if (c.vTrue > target) {
        // PDwn = 1 (Ramp): baja a RDWn; PDwn = 0 (Kill): cae de inmediato al apagar [I].
        if (!pw && c.ints["PDwn"] == 0u)
          c.vTrue = target;
        else
          c.vTrue = std::max(target, c.vTrue - c.num["RDWn"] * dt);
      }
      const double frac = v0 > 0 ? c.vTrue / v0 : 0.0;
      const double iTrue = (c.loadUA + c.extraUA) * frac;
      c.overcurrent = pw && iTrue > c.num["I0Set"];
      if (c.overcurrent) {
        c.ovcTime += dt;
        if (c.ovcTime >= c.num["Trip"] - kEps) {
          c.tripped = true;
          c.ints["Pw"] = 0u;
          c.overcurrent = false;
          c.ovcTime = 0.0;
        }
      } else {
        c.ovcTime = 0.0;
      }
    }
  }
}

void FakeHVBackend::advanceToNow() {
  // Integra por tramos entre eventos; cada evento se aplica exactamente en su instante.
  const double tNow = seconds(clock_.now_steady());
  for (;;) {
    while (!events_.empty() && events_.front().t <= simTime_ + kEps) {
      apply(events_.front());
      events_.erase(events_.begin());
    }
    if (simTime_ >= tNow - kEps) break;
    const double tNext = events_.empty() ? tNow : std::min(tNow, events_.front().t);
    while (simTime_ < tNext - kEps) {
      const double dt = std::min(kMaxStep, tNext - simTime_);
      integrate(dt);
      simTime_ += dt;
    }
    simTime_ = std::max(simTime_, tNext);
  }
}

Status FakeHVBackend::checkLink() {
  ++calls_;
  advanceToNow();
  if (inConnectionLoss(simTime_)) {
    connected_ = false;
    return Status::failure(-2, "fake: conexión perdida con el mainframe");
  }
  if (!connected_) return Status::failure(-1, "no conectado");
  return okStatus();
}

Status FakeHVBackend::connect() {
  ++calls_;
  advanceToNow();
  if (inConnectionLoss(simTime_)) return Status::failure(-2, "fake: InitSystem falló, mainframe inalcanzable");
  connected_ = true;
  return okStatus();
}

void FakeHVBackend::disconnect() { connected_ = false; }

Result<CrateMap> FakeHVBackend::crateMap() {
  if (auto st = checkLink(); !st) return Result<CrateMap>(st.error());
  CrateMap m;
  for (int i = 0; i < crateSlots_; ++i) {
    SlotInfo si;
    si.slot = i;
    if (Slot* s = findSlot(i)) si = s->info;
    m.slots.push_back(si);
  }
  return m;
}

Result<std::vector<ParamInfo>> FakeHVBackend::channelParamNames(int slot, int channel) {
  using R = Result<std::vector<ParamInfo>>;
  if (auto st = checkLink(); !st) return R(st.error());
  Slot* s = findSlot(slot);
  if (!s || channel < 0 || channel >= s->info.channels) return R::failure(-3, "fake: slot o canal inexistente");
  return params_;
}

std::uint32_t FakeHVBackend::statusWord(const Channel& c) const {
  std::uint32_t w = 0;
  const bool pw = c.ints.at("Pw") == 1u;
  const double target = pw ? c.num.at("V0Set") : 0.0;
  if (pw) w |= kOn;
  if (pw && c.vTrue < target - kEps) w |= kRampUp;
  if (c.vTrue > target + kEps) w |= kRampDown;
  if (c.overcurrent) w |= kOvercurrent;
  if (c.tripped) w |= kInternalTrip;
  return w;
}

ParamValue FakeHVBackend::value(const Slot& s, const Channel& c, const ParamInfo& info) {
  if (info.name == "VMon") {
    std::normal_distribution<double> n(0.0, s.vmonNoiseV);
    const double v = c.vTrue + c.vmonOffsetV + (s.vmonNoiseV > 0 ? n(rng_) : 0.0);
    // La placa real no informa VMon negativo (mínimo 0.000 V con HV apagada).
    return numericValue(info, std::max(0.0, quantize(v, s.vmonStepV)) + 0.0);
  }
  if (info.name == "IMon") {
    const double v0 = c.num.at("V0Set");
    const double frac = v0 > 0 ? c.vTrue / v0 : 0.0;
    double iTrue = (c.loadUA + c.extraUA) * frac;
    if (c.overcurrent) iTrue = c.num.at("I0Set");  // la placa limita la corriente
    std::normal_distribution<double> n(0.0, s.imonNoiseUA);
    const double i = iTrue + c.imonOffsetUA + (s.imonNoiseUA > 0 ? n(rng_) : 0.0);
    return numericValue(info, quantize(i, s.imonStepUA) + 0.0);  // + 0.0 evita "-0.000"
  }
  if (info.name == "Status") return integerValue(info, statusWord(c));
  if (info.type == ParamType::Numeric) return numericValue(info, c.num.at(info.name));
  return integerValue(info, c.ints.at(info.name));
}

Result<ChannelReading> FakeHVBackend::readChannels(int slot, const std::vector<int>& channels,
                                                   const std::vector<std::string>& params) {
  using R = Result<ChannelReading>;
  if (auto st = checkLink(); !st) return R(st.error());
  Slot* s = findSlot(slot);
  if (!s) return R::failure(-3, "fake: slot inexistente");
  if (channels.empty() || params.empty()) return R::failure(-1, "lista de canales o parámetros vacía");
  ChannelReading out;
  out.slot = slot;
  out.channels = channels;
  for (const auto& name : params) {
    auto it = std::find_if(params_.begin(), params_.end(), [&](const ParamInfo& p) { return p.name == name; });
    if (it == params_.end()) return R::failure(-4, "fake: parámetro inexistente " + name);
    out.params.push_back(*it);
  }
  out.tReqUtcNs = clock_.now_utc_ns();
  const SteadyTime t0 = clock_.now_steady();
  for (const auto& info : out.params) {
    std::vector<ParamValue> row;
    for (int ch : channels) {
      if (ch < 0 || ch >= s->info.channels) return R::failure(-3, "fake: canal inexistente " + std::to_string(ch));
      row.push_back(value(*s, s->channels[static_cast<std::size_t>(ch)], info));
    }
    out.values.push_back(std::move(row));
  }
  out.durationNs = std::chrono::duration_cast<std::chrono::nanoseconds>(clock_.now_steady() - t0).count();
  out.tRespUtcNs = clock_.now_utc_ns();
  return out;
}

double FakeHVBackend::trueVoltage(int slot, int channel) {
  advanceToNow();
  Slot* s = findSlot(slot);
  if (!s || channel < 0 || channel >= s->info.channels) return NAN;
  return s->channels[static_cast<std::size_t>(channel)].vTrue;
}

}  // namespace daq::hv
