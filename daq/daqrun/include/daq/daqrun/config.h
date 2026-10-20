// config.h
// Configuración de daqrun (fase 3), JSON (diseño 5). Ejemplo: daq/config/daqrun.example.json.
// Regla de simulación: con un digitizador simulado (fake o replay) raw_root y manifests_dir no
// pueden estar dentro de /data/daq/raw ni de /data/daq/manifests; si faltan en el JSON toman
// /data/daq/sim/raw y /data/daq/sim/manifests.
// Tipo de valor; se lee una vez al arrancar.
#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/result.h"
#include "daq/dig/digitizer.h"
#include "daq/daqrun/panels.h"
#include "daq/hv/status_bits.h"

namespace daq::daqrun {

inline constexpr const char* kRealRawRoot = "/data/daq/raw";
inline constexpr const char* kRealManifests = "/data/daq/manifests";
inline constexpr const char* kSimRawRoot = "/data/daq/sim/raw";
inline constexpr const char* kSimManifests = "/data/daq/sim/manifests";

// [A-Za-z0-9][A-Za-z0-9_.-]{0,63}, sin "..": nombre de campaña (directorio y prefijo de archivo).
bool validCampaignName(const std::string& s);

struct WindowLayout {
  int refOldRecordStartNs = 3344;
  std::pair<int, int> pedestal{0, 40}, control{40, 49}, signal{49, 58}, tail{58, 102}, post{102, 144};
  nlohmann::json toJson() const;
};

// Umbrales y ventanas de los bits de datos (data_flags.h, DAQ2-9).
struct DataFlagConfig {
  bool enabled = true;
  double windowS = 600;          // bit 7 (caras)
  double blackoutWindowS = 300;  // bit 6 (grupos): μ ≥ 20 en el grupo más chico de hod2019
  double imonWindowS = 120;      // aviso de IMon
  double baselineHours = 3;
  double baselineMinS = 1800;
  double minExpected = 20;
  // Bit 7 también decide con μ < min_expected si n ≥ min_expected y la base reúne al menos
  // min_baseline_counts cuentas (μ × historia / ventana): runs 7-8, base de la tarde ~1.4 kHz por
  // canal en la cara Trasera, μ ≈ 9 en 10 min con n = 250 y razón 20 (DAQ2-9).
  double minBaselineCounts = 30;  // incertidumbre relativa de la base ~18 %, frente a la razón 2
  double significanceP = 1e-6;
  double blackoutRatio = 0.1;
  double lightLeakRatio = 2.0;
  std::array<double, 2> blackoutThresholdMv{3.0, 1.5};
  std::array<double, 2> lightLeakThresholdMv{20.0, 3.0};
  int clearEvaluations = 3;
  double imonK = 1.5;
  double imonFracI0Set = 0.6;
  std::map<std::string, double> baselines;  // "face:hod2019 Trasera" o "group:3" (Hz por canal), "imon:CH00" (A)
  static DataFlagConfig fromJson(const nlohmann::json& j);
  nlohmann::json toJson() const;
};

struct DaqrunConfig {
  // campaign
  std::string campaign;  // también --campaign
  std::string rawRoot = kRealRawRoot;
  std::string manifestsDir = kRealManifests;
  std::string operatorName;
  // digitizer
  dig::DigitizerConfig digitizer;
  int tttWrapBits = 31;
  double tttPcToleranceS = 1.0;  // bit 4 de evt_flags
  // digitizer.ttt_clock (DAQ2-9): t_utc con el tick ajustado contra el PC (TttClock); ttt_tick_ns es
  // el tick provisional hasta el primer ajuste.
  bool tttClockFit = true;
  double tttClockWindowS = 300, tttClockReanchorS = 60, tttClockMinSpanS = 120, tttClockFirstFitS = 600;
  int tttClockMinPoints = 50;
  WindowLayout window;
  // run
  std::string mode = "physics";
  std::uint64_t maxEvents = 0;
  double maxDurationS = 0;
  double segmentMaxS = 3600;
  std::uint64_t segmentMaxBytes = 2000000000ULL;
  double autosaveS = 300;
  double flagDelayS = 5.0;
  int compression = 505;  // ROOT: algoritmo 5 (ZSTD) nivel 5
  // hv
  std::string hvStateFile = "/dev/shm/daqv2/hvmon_state.json";
  double hvPollPeriodS = 1.0;
  double hvStaleAfterS = 3.0;
  double hvSnapshotPeriodS = 10.0;
  std::vector<int> expectedOn{0, 1, 2, 3, 4, 5, 6, 7};
  hv::StatusBitMap statusBits = hv::StatusBitMap::wrapperRev13();
  std::vector<std::string> tripFlags{"internal_trip", "external_trip"};
  double vmonToleranceV = 0.5;
  double ovcKOverBaseline = 2.0;
  double imonBaselineWindowS = 3600;
  double tripBackdateS = 5.0;
  // notas y mapeo
  std::string notesDir = "/data/daq/slowcontrol";
  std::string mappingFile = "daq/config/mapping/hodo_mapping_2026_10.json";
  std::string calibrationFile;  // vacío = "none"
  // búsqueda de pulsos / monitor
  int minSeparationSamples = 4;
  double provisionalThresholdMv[2] = {20.0, 3.0};
  double trimMv = 2.0;
  double ampHistMaxMv = 100.0;
  double summaryPeriodS = 60;
  double gapProbability = 1e-6;
  double burstProbability = 1e-6;
  // monitoreo y colas
  bool httpEnabled = true;
  std::string httpBind = "127.0.0.1";
  int httpPort = 8080;
  double processRequestsHz = 15;
  std::size_t monitorQueueCapacity = 64;
  std::size_t writerQueueCapacity = 1024;  // bloques; ver 15.5 para la memoria que implica
  std::size_t hvQueueCapacity = 4096;
  double writerShutdownTimeoutS = 120;
  // Aviso por correo si el estado de hvmon está viejo/no vivo más de hvmon_stale_alert_min durante
  // un run (biblioteca del notificador de hvmon; outbox propio).
  nlohmann::json notifier = nlohmann::json::object();  // claves de NotifierConfig
  double hvmonStaleAlertMin = 5.0;
  double noTriggerAlertMin = 5.0;  // runs con trigger externo: aviso "sin triggers"
  double minTriggerRateHz = 0.1;   // ... o tasa bajo este valor en la ventana
  double triggerRateWindowMin = 10.0;
  DataFlagConfig dataFlags;  // data_flags
  PanelConfig panels;        // panels (paneles en tiempo, fase 5)

  bool simulated() const { return digitizer.backend != "caen"; }
  std::string campaignDir() const { return rawRoot + "/" + campaign; }

  static Result<DaqrunConfig> fromJson(const nlohmann::json& j);
  static Result<DaqrunConfig> fromFile(const std::string& path);
  nlohmann::json toJson() const;
};

}  // namespace daq::daqrun
