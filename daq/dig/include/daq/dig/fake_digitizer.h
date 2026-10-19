// fake_digitizer.h
// FakeDigitizer sintético (fase 3): genera eventos del hodoscopio con un escenario JSON y un reloj
// inyectable (SystemClock en tiempo real, ManualClock en tests acelerados).
//  - Triggers de Poisson a trigger_rate_Hz; huecos (trigger_gap) y ráfagas (trigger_burst).
//  - Por canal: línea base (baseline_adc), ruido gaussiano (noise_mV: hod2018 ~1.9, hod2019 ~0.25 mV,
//    estudio de ventana de DAQ2-1), pulsos oscuros de Poisson (dark_rate_Hz sobre 0.5 p.e. en
//    hod2018, ~2.3 MHz; hod2019 ~24 kHz) con ganancia gain_mV_per_pe y cola de crosstalk geométrica.
//  - Muón: traza recta por las dos estaciones (posición uniforme en las 16 tiras, pendiente
//    gaussiana), pulso en signal_sample en las tiras X e Y cruzadas.
//  - Guion: hv_group_off (el grupo pierde pulsos oscuros y de muón), light_leak (tasa extra en una
//    cara física), gain_drift / gain_step, trigger_gap, trigger_burst.
//  - TTT con tick ttt_tick_ns y vuelta en ttt_wrap_bits; contador de 24 bit; memoria de la placa
//    board_memory_events (si se llena se pierden eventos y se marca board_full).
// Seguridad entre hilos: pertenece al hilo de lectura (como el backend real).
#pragma once

#include <deque>
#include <random>
#include <string>
#include <vector>

#include "daq/core/clock.h"
#include "daq/dig/digitizer.h"

namespace daq::dig {

// Geometría y cableado de un canal del digitizador (del mapeo versionado).
struct ChannelInfo {
  int station = 0;  // 0 = hod2018, 1 = hod2019
  char plane = 'X';
  int strip = 0;
  int hvGroup = 0;  // 1-8; 0 = desconocido
  std::string face;
};
// Mapeo por omisión sin grupos: hod2018 X = ch 0-15, Y = 16-31; hod2019 X = 32-47, Y = 48-63.
std::vector<ChannelInfo> defaultChannelInfo();

struct SyntheticParams {
  std::uint64_t seed = 1;
  double triggerRateHz = 0.6;
  int signalSample = 51;
  double baselineAdc[2] = {2091, 2081};
  double baselineSpreadAdc = 4;  // dispersión por canal (uniforme ±)
  double noiseMv[2] = {1.9, 0.25};
  double darkRateHz[2] = {2.3e6, 2.4e4};
  double gainMinMv[2] = {4.0, 4.0};
  double gainMaxMv[2] = {5.5, 4.0};
  double crosstalkProb[2] = {0.2, 0.05};
  double decaySamples = 2.0;
  double muonPeMean[2] = {20, 12};
  double muonEfficiency = 0.95;
  double slopeSigmaStrips = 1.0;
  int tttWrapBits = 31;
  int counterBits = 24;
  std::size_t boardMemoryEvents = 1024;
  nlohmann::json events = nlohmann::json::array();
  static Result<SyntheticParams> fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
};

class FakeDigitizer final : public DigitizerBackend {
 public:
  FakeDigitizer(SyntheticParams p, Clock& clock, std::vector<ChannelInfo> channels = defaultChannelInfo());

  Status configure(const DigitizerConfig& cfg) override;
  Result<AcqOrigin> start() override;
  Status stop() override;
  Result<std::size_t> readBlock(EventBlock& out, std::size_t maxEvents) override;
  Status sendSoftwareTrigger() override;
  nlohmann::json info() const override;
  bool simulated() const override { return true; }
  int recordLength() const override { return cfg_.recordLength; }

  // Para tests: verdad del generador.
  std::uint64_t triggersGenerated() const { return generated_; }
  std::uint64_t eventsLostBoardFull() const { return lostFull_; }
  double gainMv(int ch, UtcNs t) const;
  bool channelOff(int ch, UtcNs t) const;  // HV de su grupo apagada por el guion
  const std::vector<ChannelInfo>& channels() const { return ch_; }

 private:
  struct Pending {
    UtcNs t;
    std::uint8_t kind;      // 0 externo, 1 software
    std::uint32_t counter;  // contador de la placa al aceptar el trigger
  };
  void generateUntil(UtcNs now);
  UtcNs drawNextPoisson(UtcNs from);
  Event synthesize(const Pending& p);
  double lightLeakExtraHz(int ch, UtcNs t) const;
  bool inGap(UtcNs t, UtcNs& gapEnd) const;

  SyntheticParams p_;
  Clock& clock_;
  std::vector<ChannelInfo> ch_;
  DigitizerConfig cfg_;
  bool configured_ = false, running_ = false;
  std::mt19937_64 rng_;
  UtcNs t0_ = 0;          // inicio de adquisición (origen del TTT)
  UtcNs nextPoisson_ = 0;
  UtcNs generatedUntil_ = 0;
  std::vector<UtcNs> injected_;  // ráfagas y SW triggers pendientes de entrar a la memoria
  std::deque<Pending> memory_;
  std::uint32_t counter_ = 0;
  std::uint64_t generated_ = 0, lostFull_ = 0;
  bool fullSinceLastRead_ = false;
  std::vector<double> baseline_, gainMv_;
};

}  // namespace daq::dig
