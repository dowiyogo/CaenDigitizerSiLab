// calibration.cpp
#include "daq/daqrun/calibration.h"

#include <cmath>
#include <fstream>

namespace daq::daqrun {

Result<PeCalibration> PeCalibration::fromJson(const nlohmann::json& j) {
  using R = Result<PeCalibration>;
  PeCalibration c;
  try {
    if (j.value("schema", "") != "daqv2.pe_calibration/1") return R::failure(-1, "calibración: schema distinto de daqv2.pe_calibration/1");
    c.version = j.at("version").get<std::string>();
    if (c.version.empty()) return R::failure(-1, "calibración sin versión");
    c.provisional.fill(true);
    c.gainMv.fill(NAN);
    c.pedestalMv.fill(NAN);
    int n = 0;
    for (const auto& e : j.at("channels")) {
      const int ch = e.at("ch").get<int>();
      if (ch < 0 || ch >= 64) return R::failure(-1, "calibración: canal fuera de rango");
      const auto i = static_cast<std::size_t>(ch);
      if (e.contains("G_mV") && e["G_mV"].is_number()) c.gainMv[i] = e["G_mV"].get<double>();
      if (e.contains("pedestal_mV") && e["pedestal_mV"].is_number()) c.pedestalMv[i] = e["pedestal_mV"].get<double>();
      c.provisional[i] = e.value("provisional", true) || !std::isfinite(c.gainMv[i]);
      ++n;
    }
    if (n != 64) return R::failure(-1, "calibración: deben estar los 64 canales");
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("calibración inválida: ") + e.what());
  }
  c.raw = j;
  return c;
}

Result<PeCalibration> PeCalibration::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) return Result<PeCalibration>::failure(-1, "no se pudo abrir la calibración " + path);
  nlohmann::json j;
  try {
    j = nlohmann::json::parse(in);
  } catch (const std::exception& e) {
    return Result<PeCalibration>::failure(-1, "calibración " + path + ": " + e.what());
  }
  auto r = fromJson(j);
  if (r) {
    auto c = r.value();
    c.path = path;
    return c;
  }
  return r;
}

}  // namespace daq::daqrun
