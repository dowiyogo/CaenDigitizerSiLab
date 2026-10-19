// caen_hv_backend.cpp
// Implementación sobre la Wrapper: InitSystem, DeinitSystem, GetCrateMap, GetChParamInfo,
// GetChParamProp, GetChParam, GetBd*, GetError, Free, CAENHVLibSwRel y dos escrituras:
// CAENHV_SetChParam("Pw") en setChannelPw() y CAENHV_SetBdParam("ClrAlarm") en clearBoardAlarm().
#include "daq/hv/caen_hv_backend.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// El header de la Wrapper define macros uchar/ushort/ulong (ulong = unsigned int). Se incluye
// al final y se anulan las macros para que no contaminen el resto del archivo.
#include <CAENHVWrapper.h>
#undef uchar
#undef ushort
#undef ulong

namespace daq::hv {

namespace {

struct CaenFree {
  void operator()(void* p) const noexcept {
    if (p) CAENHV_Free(p);
  }
};
template <class T>
using CaenPtr = std::unique_ptr<T, CaenFree>;

Result<CAENHV_SYSTEM_TYPE_t> systemType(const std::string& name) {
  if (name == "SY5527") return SY5527;
  if (name == "SY4527") return SY4527;
  return Result<CAENHV_SYSTEM_TYPE_t>::failure(-1, "sistema HV no soportado: " + name);
}

// Unidad de respaldo si la placa no informa Unit/Exp: (código PARAM_UN_*, exponente).
std::pair<int, int> fallbackUnit(const std::string& name) {
  if (name == "IMon" || name == "I0Set" || name == "I1Set") return {PARAM_UN_AMPERE, -6};
  if (name == "RUp" || name == "RDWn") return {PARAM_UN_VPS, 0};
  if (name == "Trip") return {PARAM_UN_SECOND, 0};
  if (!name.empty() && name[0] == 'V') return {PARAM_UN_VOLT, 0};
  if (name == "SVMax") return {PARAM_UN_VOLT, 0};
  return {PARAM_UN_NONE, 0};
}

// Buffer generoso y alineado para GetChParamProp: el tamaño del dato depende de la propiedad
// (Type/Mode: ulong = unsigned int; Unit: unsigned short; Exp: short; Minval/Maxval: float;
// Onstate/Offstate: texto).
struct PropBuffer {
  alignas(8) std::array<unsigned char, 64> bytes{};
  template <class T>
  T as() const {
    T v{};
    std::memcpy(&v, bytes.data(), sizeof(T));
    return v;
  }
  std::string text() const {
    const char* p = reinterpret_cast<const char*>(bytes.data());
    return std::string(p, strnlen(p, bytes.size()));
  }
};

// La Wrapper entrega float (~7 cifras significativas). Se pasa a double por su representación
// decimal más corta de 7 cifras, para no inventar dígitos (0.048f -> 0.048, no 0.0480000004).
double floatToDouble(float f) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.7g", static_cast<double>(f));
  return std::strtod(buf, nullptr);
}

}  // namespace

CaenHVBackend::CaenHVBackend(ConnectionConfig config, const Clock* clock)
    : config_(std::move(config)) {
  if (clock == nullptr) {
    ownClock_ = std::make_unique<SystemClock>();
    clock_ = ownClock_.get();
  } else {
    clock_ = clock;
  }
}

CaenHVBackend::~CaenHVBackend() { disconnect(); }

std::string CaenHVBackend::librarySoftwareRelease() {
  const char* r = CAENHVLibSwRel();
  return r ? std::string(r) : std::string();
}

Error CaenHVBackend::wrapperError(int code, const std::string& where) const {
  const char* msg = handle_ >= 0 ? CAENHV_GetError(handle_) : nullptr;
  char hex[16];
  std::snprintf(hex, sizeof hex, "0x%X", static_cast<unsigned>(code));
  return Error{code, where + " falló (" + hex + "): " + (msg ? msg : "sin texto de la Wrapper")};
}

Status CaenHVBackend::connect() {
  if (isConnected()) return okStatus();
  auto sys = systemType(config_.system);
  if (!sys) return Status(sys.error());
  std::vector<char> arg(config_.host.begin(), config_.host.end());
  arg.push_back('\0');
  int h = -1;
  const CAENHVRESULT r = CAENHV_InitSystem(sys.value(), LINKTYPE_TCPIP, arg.data(),
                                           config_.credentials.user.c_str(),
                                           config_.credentials.password.c_str(), &h);
  if (r != CAENHV_OK) {
    const char* msg = CAENHV_GetError(h);
    char hex[16];
    std::snprintf(hex, sizeof hex, "0x%X", static_cast<unsigned>(r));
    return Status::failure(r, std::string("InitSystem ") + config_.host + " falló (" + hex +
                                  "): " + (msg ? msg : "sin texto de la Wrapper"));
  }
  handle_ = h;
  infoCache_.clear();
  return okStatus();
}

void CaenHVBackend::disconnect() {
  if (handle_ >= 0) {
    CAENHV_DeinitSystem(handle_);
    handle_ = -1;
  }
}

std::string CaenHVBackend::describe() const {
  return "CaenHVBackend(" + config_.system + "@" + config_.host + ", Wrapper " + librarySoftwareRelease() + ")";
}

Result<CrateMap> CaenHVBackend::crateMap() {
  if (!isConnected()) return Result<CrateMap>::failure(-1, "no conectado");
  unsigned short nSlots = 0;
  unsigned short* nCh = nullptr;
  unsigned short* serial = nullptr;
  unsigned char* fwMin = nullptr;
  unsigned char* fwMax = nullptr;
  char* models = nullptr;
  char* descr = nullptr;
  const CAENHVRESULT r = CAENHV_GetCrateMap(handle_, &nSlots, &nCh, &models, &descr, &serial, &fwMin, &fwMax);
  CaenPtr<unsigned short> pNch(nCh), pSerial(serial);
  CaenPtr<unsigned char> pMin(fwMin), pMax(fwMax);
  CaenPtr<char> pModels(models), pDescr(descr);
  if (r != CAENHV_OK) return wrapperError(r, "GetCrateMap");
  CrateMap map;
  const char* m = models;
  const char* d = descr;
  for (unsigned short s = 0; s < nSlots; ++s) {
    SlotInfo si;
    si.slot = s;
    si.model = m ? m : "";
    si.description = d ? d : "";
    si.channels = nCh ? nCh[s] : 0;
    si.serial = serial ? serial[s] : 0u;
    si.fwMajor = fwMax ? fwMax[s] : 0;
    si.fwMinor = fwMin ? fwMin[s] : 0;
    map.slots.push_back(si);
    if (m) m += std::strlen(m) + 1;
    if (d) d += std::strlen(d) + 1;
  }
  return map;
}

Result<ParamInfo> CaenHVBackend::paramInfo(int slot, int channel, const std::string& name) {
  auto key = std::make_pair(slot, name);
  if (auto it = infoCache_.find(key); it != infoCache_.end()) return it->second;
  const auto s = static_cast<unsigned short>(slot);
  const auto c = static_cast<unsigned short>(channel);
  ParamInfo info;
  info.name = name;
  PropBuffer b;
  CAENHVRESULT r = CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Type", b.bytes.data());
  if (r != CAENHV_OK) return wrapperError(r, "GetChParamProp(" + name + ", Type)");
  info.type = paramTypeFromInt(static_cast<int>(b.as<unsigned int>()));
  b = PropBuffer{};
  if (CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Mode", b.bytes.data()) == CAENHV_OK)
    info.mode = paramModeFromInt(static_cast<int>(b.as<unsigned int>()));
  if (info.type == ParamType::Numeric) {
    PropBuffer bu, be;
    const bool okUnit = CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Unit", bu.bytes.data()) == CAENHV_OK;
    const bool okExp = CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Exp", be.bytes.data()) == CAENHV_OK;
    if (okUnit && okExp) {
      info.unitCode = bu.as<unsigned short>();
      info.exponent = be.as<short>();
      info.unitFromBoard = true;
    } else {
      const auto fb = fallbackUnit(name);
      info.unitCode = fb.first;
      info.exponent = fb.second;
      info.unitFromBoard = false;
    }
    PropBuffer bmin, bmax;
    if (CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Minval", bmin.bytes.data()) == CAENHV_OK)
      info.minBoard = static_cast<double>(bmin.as<float>());
    if (CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Maxval", bmax.bytes.data()) == CAENHV_OK)
      info.maxBoard = static_cast<double>(bmax.as<float>());
  } else if (info.type == ParamType::OnOff) {
    PropBuffer bon, boff;
    if (CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Onstate", bon.bytes.data()) == CAENHV_OK)
      info.onState = bon.text();
    if (CAENHV_GetChParamProp(handle_, s, c, name.c_str(), "Offstate", boff.bytes.data()) == CAENHV_OK)
      info.offState = boff.text();
  }
  infoCache_[key] = info;
  return info;
}

Result<std::vector<ParamInfo>> CaenHVBackend::channelParamNames(int slot, int channel) {
  using R = Result<std::vector<ParamInfo>>;
  if (!isConnected()) return R::failure(-1, "no conectado");
  char* list = nullptr;
  int n = 0;
  const CAENHVRESULT r = CAENHV_GetChParamInfo(handle_, static_cast<unsigned short>(slot),
                                               static_cast<unsigned short>(channel), &list, &n);
  CaenPtr<char> pList(list);
  if (r != CAENHV_OK) return wrapperError(r, "GetChParamInfo");
  std::vector<ParamInfo> out;
  for (int p = 0; p < n; ++p) {
    const char* name = list + static_cast<std::ptrdiff_t>(p) * MAX_PARAM_NAME;
    auto info = paramInfo(slot, channel, std::string(name, strnlen(name, MAX_PARAM_NAME)));
    if (!info) return R(info.error());
    out.push_back(info.value());
  }
  return out;
}

Result<ChannelReading> CaenHVBackend::readChannels(int slot, const std::vector<int>& channels,
                                                   const std::vector<std::string>& params) {
  using R = Result<ChannelReading>;
  if (!isConnected()) return R::failure(-1, "no conectado");
  if (channels.empty() || params.empty()) return R::failure(-1, "lista de canales o parámetros vacía");
  ChannelReading out;
  out.slot = slot;
  out.channels = channels;
  // Las propiedades se consultan (y quedan en caché) antes de iniciar la medida de tiempo.
  for (const auto& name : params) {
    auto info = paramInfo(slot, channels.front(), name);
    if (!info) return R(info.error());
    out.params.push_back(info.value());
  }
  std::vector<unsigned short> list(channels.size());
  for (std::size_t i = 0; i < channels.size(); ++i) list[i] = static_cast<unsigned short>(channels[i]);
  const auto nCh = static_cast<unsigned short>(channels.size());
  std::vector<float> fbuf(channels.size());
  std::vector<unsigned int> ubuf(channels.size());

  out.tReqUtcNs = clock_->now_utc_ns();
  const SteadyTime t0 = clock_->now_steady();
  for (const auto& info : out.params) {
    std::vector<ParamValue> row;
    row.reserve(channels.size());
    if (info.type == ParamType::Numeric) {
      const CAENHVRESULT r = CAENHV_GetChParam(handle_, static_cast<unsigned short>(slot), info.name.c_str(),
                                               nCh, list.data(), fbuf.data());
      if (r != CAENHV_OK) return wrapperError(r, "GetChParam(" + info.name + ")");
      for (float f : fbuf) row.push_back(numericValue(info, floatToDouble(f)));
    } else {
      const CAENHVRESULT r = CAENHV_GetChParam(handle_, static_cast<unsigned short>(slot), info.name.c_str(),
                                               nCh, list.data(), ubuf.data());
      if (r != CAENHV_OK) return wrapperError(r, "GetChParam(" + info.name + ")");
      for (unsigned int u : ubuf) row.push_back(integerValue(info, u));
    }
    out.values.push_back(std::move(row));
  }
  out.durationNs = std::chrono::duration_cast<std::chrono::nanoseconds>(clock_->now_steady() - t0).count();
  out.tRespUtcNs = clock_->now_utc_ns();
  return out;
}

Result<ParamInfo> CaenHVBackend::boardParamInfo(int slot, const std::string& name) {
  auto key = std::make_pair(slot, name);
  if (auto it = bdInfoCache_.find(key); it != bdInfoCache_.end()) return it->second;
  const auto s = static_cast<unsigned short>(slot);
  ParamInfo info;
  info.name = name;
  PropBuffer b;
  CAENHVRESULT r = CAENHV_GetBdParamProp(handle_, s, name.c_str(), "Type", b.bytes.data());
  if (r != CAENHV_OK) return wrapperError(r, "GetBdParamProp(" + name + ", Type)");
  info.type = paramTypeFromInt(static_cast<int>(b.as<unsigned int>()));
  b = PropBuffer{};
  if (CAENHV_GetBdParamProp(handle_, s, name.c_str(), "Mode", b.bytes.data()) == CAENHV_OK)
    info.mode = paramModeFromInt(static_cast<int>(b.as<unsigned int>()));
  if (info.type == ParamType::Numeric) {
    PropBuffer bu, be, bmin, bmax;
    const bool okUnit = CAENHV_GetBdParamProp(handle_, s, name.c_str(), "Unit", bu.bytes.data()) == CAENHV_OK;
    const bool okExp = CAENHV_GetBdParamProp(handle_, s, name.c_str(), "Exp", be.bytes.data()) == CAENHV_OK;
    if (okUnit && okExp) {
      info.unitCode = bu.as<unsigned short>();
      info.exponent = be.as<short>();
      info.unitFromBoard = true;
    } else {
      info.unitCode = name == "Temp" ? PARAM_UN_CELSIUS : PARAM_UN_NONE;
      info.exponent = 0;
    }
    if (CAENHV_GetBdParamProp(handle_, s, name.c_str(), "Minval", bmin.bytes.data()) == CAENHV_OK)
      info.minBoard = static_cast<double>(bmin.as<float>());
    if (CAENHV_GetBdParamProp(handle_, s, name.c_str(), "Maxval", bmax.bytes.data()) == CAENHV_OK)
      info.maxBoard = static_cast<double>(bmax.as<float>());
  }
  bdInfoCache_[key] = info;
  return info;
}

Result<std::vector<ParamInfo>> CaenHVBackend::boardParamNames(int slot) {
  using R = Result<std::vector<ParamInfo>>;
  if (!isConnected()) return R::failure(-1, "no conectado");
  char* list = nullptr;
  const CAENHVRESULT r = CAENHV_GetBdParamInfo(handle_, static_cast<unsigned short>(slot), &list);
  CaenPtr<char> pList(list);
  if (r != CAENHV_OK) return wrapperError(r, "GetBdParamInfo");
  if (!list) return R::failure(-1, "GetBdParamInfo: lista vacía");
  // La Wrapper no da la cantidad: bloques de MAX_PARAM_NAME bytes terminados por un nombre vacío.
  std::vector<ParamInfo> out;
  for (int p = 0; p < 64; ++p) {
    const char* name = list + static_cast<std::ptrdiff_t>(p) * MAX_PARAM_NAME;
    if (name[0] == '\0') break;
    auto info = boardParamInfo(slot, std::string(name, strnlen(name, MAX_PARAM_NAME)));
    if (!info) return R(info.error());
    out.push_back(info.value());
  }
  return out;
}

Result<BoardReading> CaenHVBackend::readBoard(int slot, const std::vector<std::string>& params) {
  using R = Result<BoardReading>;
  if (!isConnected()) return R::failure(-1, "no conectado");
  BoardReading out;
  out.slot = slot;
  for (const auto& name : params) {
    auto info = boardParamInfo(slot, name);
    if (!info) return R(info.error());
    out.params.push_back(info.value());
  }
  const unsigned short s = static_cast<unsigned short>(slot);
  out.tReqUtcNs = clock_->now_utc_ns();
  const SteadyTime t0 = clock_->now_steady();
  for (const auto& info : out.params) {
    if (info.type == ParamType::Numeric) {
      float f = 0;
      const CAENHVRESULT r = CAENHV_GetBdParam(handle_, 1, &s, info.name.c_str(), &f);
      if (r != CAENHV_OK) return wrapperError(r, "GetBdParam(" + info.name + ")");
      out.values.push_back(numericValue(info, floatToDouble(f)));
    } else {
      unsigned int u = 0;
      const CAENHVRESULT r = CAENHV_GetBdParam(handle_, 1, &s, info.name.c_str(), &u);
      if (r != CAENHV_OK) return wrapperError(r, "GetBdParam(" + info.name + ")");
      out.values.push_back(integerValue(info, u));
    }
  }
  out.durationNs = std::chrono::duration_cast<std::chrono::nanoseconds>(clock_->now_steady() - t0).count();
  out.tRespUtcNs = clock_->now_utc_ns();
  return out;
}

// ÚNICA llamada de escritura de todo el código base: Pw de un canal.
Status CaenHVBackend::setChannelPw(int slot, int channel, bool on) {
  if (!isConnected()) return Status::failure(-1, "no conectado");
  unsigned short ch = static_cast<unsigned short>(channel);
  unsigned int value = on ? 1u : 0u;
  const CAENHVRESULT r = CAENHV_SetChParam(handle_, static_cast<unsigned short>(slot), "Pw", 1, &ch, &value);
  if (r != CAENHV_OK) return Status(wrapperError(r, "SetChParam(Pw)"));
  return okStatus();
}

// Segunda (y última) escritura del código base: ClrAlarm de la placa.
Status CaenHVBackend::clearBoardAlarm(int slot) {
  if (!isConnected()) return Status::failure(-1, "no conectado");
  unsigned short s = static_cast<unsigned short>(slot);
  unsigned int value = 1u;
  const CAENHVRESULT r = CAENHV_SetBdParam(handle_, 1, &s, "ClrAlarm", &value);
  if (r != CAENHV_OK) return Status(wrapperError(r, "SetBdParam(ClrAlarm)"));
  return okStatus();
}

}  // namespace daq::hv
