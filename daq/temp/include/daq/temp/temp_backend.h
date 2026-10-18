// temp_backend.h
// Interfaz de lectura de temperatura (DS18B20 por adaptador 1-Wire USB o por microcontrolador
// serie; hardware aún no elegido, docs/daq_v2_design.md, sección 12). En la fase 1 solo existe
// el fake. Unidades: °C (única magnitud fuera de SI, con la unidad en el nombre) y UTC ns.
// Seguridad entre hilos: una instancia pertenece a un solo hilo.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/core/result.h"

namespace daq::temp {

struct TempReading {
  std::string sensorId;
  double tDegC = 0.0;
  UtcNs tUtcNs = 0;
};

class TempBackend {
 public:
  virtual ~TempBackend() = default;
  virtual Result<std::vector<TempReading>> read() = 0;
  virtual std::string describe() const = 0;
};

// Sensores simulados: T = base + amplitud·sin(2π·t/período), t en s desde la construcción.
// Configuración: {"sensors": [{"id": "fake-1", "base_degC": 22, "amplitude_degC": 3,
//                             "period_s": 86400}], "fail": false}
class FakeTempBackend final : public TempBackend {
 public:
  FakeTempBackend(const Clock& clock, const nlohmann::json& config);
  Result<std::vector<TempReading>> read() override;
  std::string describe() const override { return "FakeTempBackend"; }

 private:
  struct Sensor {
    std::string id;
    double base = 22.0, amplitude = 0.0, period = 86400.0;
  };
  const Clock& clock_;
  SteadyTime t0_;
  std::vector<Sensor> sensors_;
  bool fail_ = false;
};

// "none" -> nullptr; "fake" -> FakeTempBackend; otro -> error.
Result<std::unique_ptr<TempBackend>> makeTempBackend(const std::string& kind, const Clock& clock,
                                                     const nlohmann::json& config);

}  // namespace daq::temp
