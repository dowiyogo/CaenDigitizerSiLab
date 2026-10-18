// temp_backend.cpp
#include "daq/temp/temp_backend.h"

#include <cmath>

namespace daq::temp {

FakeTempBackend::FakeTempBackend(const Clock& clock, const nlohmann::json& config)
    : clock_(clock), t0_(clock.now_steady()) {
  fail_ = config.value("fail", false);
  for (const auto& s : config.value("sensors", nlohmann::json::array())) {
    Sensor x;
    x.id = s.value("id", "fake");
    x.base = s.value("base_degC", 22.0);
    x.amplitude = s.value("amplitude_degC", 0.0);
    x.period = s.value("period_s", 86400.0);
    sensors_.push_back(x);
  }
  if (sensors_.empty()) sensors_.push_back(Sensor{"fake-1", 22.0, 0.0, 86400.0});
}

Result<std::vector<TempReading>> FakeTempBackend::read() {
  if (fail_) return Result<std::vector<TempReading>>::failure(-1, "fake: sensor sin respuesta");
  const double t = std::chrono::duration<double>(clock_.now_steady() - t0_).count();
  const UtcNs now = clock_.now_utc_ns();
  std::vector<TempReading> out;
  constexpr double kTwoPi = 6.283185307179586;
  for (const auto& s : sensors_)
    out.push_back(TempReading{s.id, s.base + s.amplitude * std::sin(kTwoPi * t / s.period), now});
  return out;
}

Result<std::unique_ptr<TempBackend>> makeTempBackend(const std::string& kind, const Clock& clock,
                                                     const nlohmann::json& config) {
  using R = Result<std::unique_ptr<TempBackend>>;
  if (kind == "none" || kind.empty()) return std::unique_ptr<TempBackend>();
  if (kind == "fake") return std::unique_ptr<TempBackend>(std::make_unique<FakeTempBackend>(clock, config));
  return R::failure(-1, "temperatura: backend desconocido '" + kind + "' (en la fase 1 solo 'none' o 'fake')");
}

}  // namespace daq::temp
