// replay_digitizer.h
// FakeDigitizer en modo replay (fase 3): reproduce lectura_vx1740_ext_nim_1M*.root (v8) como si
// viniera de la placa con la ventana B.
//  - El archivo se resuelve con hodo::resolveRawFile() (analysis/src/hodo_paths.h) y se abre en
//    modo READ; se lee entrada por entrada con solo las ramas Ch0..Ch63, event y timestamp
//    habilitadas (nunca se carga el archivo entero).
//  - Las muestras en mV del DAQ antiguo vuelven a cuentas ADC con la inversa de adc2mV
//    (AdcConversion, dc_offset 0x7FFF de polarización 0) y se redondean; el residuo máximo de
//    redondeo se publica (max_rounding_residual_adc).
//  - Recorte: muestras [old_start_sample, old_start_sample + N) del registro de 624, con
//    old_start_sample = ref_old_record_start_ns / 16 = 209 para la ventana B.
//  - Ritmo: speed = 1 respeta los intervalos de `timestamp` (µs del steady_clock del DAQ antiguo,
//    uno por lectura), speed > 1 acelera, speed = 0 entrega todo lo antes posible.
//  - TTT sintetizado al ritmo de entrega: intervalos de `timestamp` divididos por speed (speed = 0:
//    hora de lectura), con ttt_tick_ns y vuelta en ttt_wrap_bits, para que t_utc sea coherente con
//    t_read; contador = `event` del original (identifica la entrada de origen).
// Seguridad entre hilos: pertenece al hilo de lectura; los objetos ROOT (TFile, TTree) son suyos.
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "daq/core/clock.h"
#include "daq/dig/digitizer.h"

class TFile;
class TTree;

namespace daq::dig {

struct ReplayParams {
  std::string file = "lectura_vx1740_ext_nim_1M_v8.root";
  std::string tree = "data_dig";
  long long firstEntry = 0;
  long long maxEntries = 0;  // 0 = hasta el final
  double speed = 1.0;
  int oldStartSample = 209;
  int dcOffset = 0x7FFF;
  int tttWrapBits = 31;
  static Result<ReplayParams> fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
};

class ReplayDigitizer final : public DigitizerBackend {
 public:
  ReplayDigitizer(ReplayParams p, Clock& clock);
  ~ReplayDigitizer() override;

  Status configure(const DigitizerConfig& cfg) override;
  Result<AcqOrigin> start() override;
  Status stop() override;
  Result<std::size_t> readBlock(EventBlock& out, std::size_t maxEvents) override;
  nlohmann::json info() const override;
  bool simulated() const override { return true; }
  bool exhausted() const override { return done_; }

  const std::string& resolvedPath() const { return path_; }
  double maxRoundingResidualAdc() const { return maxResidual_; }
  long long entriesDelivered() const { return delivered_; }

 private:
  Status open();
  ReplayParams p_;
  Clock& clock_;
  DigitizerConfig cfg_;
  std::string path_;
  std::unique_ptr<TFile> file_;
  TTree* tree_ = nullptr;  // pertenece a file_
  std::vector<std::vector<float>*> ch_;
  int event_ = 0;
  long timestamp_ = 0;  // Long_t en el árbol antiguo
  long long next_ = 0, end_ = 0;
  long long ts0_ = -1;
  UtcNs t0_ = 0;
  bool running_ = false, done_ = false, configured_ = false;
  bool havePending_ = false;
  double maxResidual_ = 0;
  long long delivered_ = 0;
};

}  // namespace daq::dig
