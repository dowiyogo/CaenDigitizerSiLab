// run_writer.h
// Escritor ROOT de daqrun (diseño 4.1, 4.2, 3.1). Pertenece al hilo escritor: es el único dueño
// del TFile y de sus TTree.
//  - Eventos: línea de retardo flag_delay_s (5 s) antes de escribir, para fijar los bits HV con la
//    lectura posterior (FlagEngine). Flags de evento: EventFlagger. TTT desenrollado con la hora del
//    PC; t_utc_ns = origen + ttt_ns, con el tick ajustado contra el PC por run (TttClock, DAQ2-9;
//    método y anclas por segmento en run_info_json.ttt_clock).
//  - Árboles: events (adc UShort_t[64][N], ch_flags[64], evt_flags, ...), hv (cada
//    hv_snapshot_period_s y en cada cambio, más uno al abrir cada segmento), hv_changes (dedup por
//    id), notes (las del intervalo del segmento, al cerrarlo), monitor (resúmenes por período).
//  - Segmentos: <base>_segMMM.root; cierre por segment_max_s (reloj de pared) o segment_max_bytes;
//    AutoSave cada autosave_s. Al cerrar: cabecera (config_json, mapping_json, mapping_version,
//    mapping_sha256, calibration_json/version, run_info_json, hvmon_state_at_start), Write+Close,
//    fsync, chmod 0444, SHA256 en SHA256SUMS, entrada en el índice (reemplazo atómico).
//  - finish(): libera la línea de retardo, cierra el segmento, índice con complete = true, 0444,
//    su SHA256 en SHA256SUMS, y el run en el manifiesto de la campaña.
// Seguridad entre hilos: todos los métodos en el hilo escritor; WriterStats es atómico y se lee
// desde cualquier hilo.
#pragma once

#include <array>
#include <atomic>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <RtypesCore.h>
#include <nlohmann/json.hpp>

#include "daq/core/clock.h"
#include "daq/core/result.h"
#include "daq/daqrun/campaign.h"
#include "daq/daqrun/config.h"
#include "daq/daqrun/hv_client.h"
#include "daq/daqrun/monitor.h"
#include "daq/dig/digitizer.h"
#include "daq/hvmon/mapping.h"

class TFile;
class TTree;

namespace daq::io {

struct RunStatic {
  nlohmann::json configJson;
  std::string configSha256;
  std::string mappingText;  // contenido del archivo de mapeo
  nlohmann::json digitizerInfo;       // al arrancar
  nlohmann::json digitizerInfoFinal;  // al terminar la lectura (p. ej. residuo de redondeo del replay)
  dig::AcqOrigin origin;
  nlohmann::json hvAtStart;  // estado de hvmon al arrancar (o null)
  std::string gitCommit = "unknown";
  bool gitDirty = false;
  std::string host;
  std::string rootVersion;
};

struct WriterStats {
  std::atomic<std::uint64_t> eventsWritten{0}, segmentsClosed{0}, pendingDelay{0}, hvRows{0}, hvChangeRows{0};
  std::atomic<std::uint64_t> bytesClosed{0}, autosaves{0}, maxPendingDelay{0};
  std::array<std::atomic<std::uint64_t>, 16> chFlagEvents{};   // eventos con el bit en algún canal
  std::array<std::atomic<std::uint64_t>, 9> evtFlagEvents{};
  std::atomic<std::uint64_t> lastEvtSeq{0};
  std::string currentSegment() const;
  void setCurrentSegment(const std::string& s);
  nlohmann::json toJson() const;

 private:
  mutable std::mutex m_;
  std::string seg_;
};

class RunWriter {
 public:
  RunWriter(const daqrun::DaqrunConfig& cfg, const daqrun::RunReservation& run, const hvmon::Mapping& mapping,
            RunStatic st, Clock& clock, WriterStats& stats);
  ~RunWriter();
  Status open();
  Status addBlock(const dig::EventBlockPtr& b);
  void addHv(const daqrun::HvSample& s);
  void addSummary(const daqrun::MonitorSummary& m);
  // Libera eventos vencidos, AutoSave y (si allowRotate) rotación por tiempo/tamaño. Durante el
  // cierre del run se llama con allowRotate = false: no se abre un segmento vacío al final.
  Status tick(bool allowRotate = true);
  Status finish(bool cleanShutdown);  // cierre del run
  // Información del digitizador al terminar la lectura; va en run_info_json del último segmento.
  // Se llama desde el hilo escritor después de que la cola de eventos se cerró (el hilo de lectura
  // ya no usa el backend).
  void setFinalDigitizerInfo(nlohmann::json j) { st_.digitizerInfoFinal = std::move(j); }
  const nlohmann::json& index() const { return index_; }
  std::string segmentPath(int seg) const;
  std::string indexPath() const;

 private:
  struct Pending {
    dig::EventBlockPtr blk;
    std::size_t idx;
    std::uint64_t seq;
    std::int64_t tttNs;
    UtcNs tUtc;
    std::int64_t dtPrev;
    std::uint32_t evtFlags;
  };
  Status openSegment();
  Status closeSegment(bool clean);
  Status release(bool all);
  Status writeEvent(const Pending& p);
  void fillHv(const daqrun::HvSample& s);
  bool hvChanged(const daqrun::HvSample& s) const;
  nlohmann::json runInfo(bool final) const;
  Status writeIndex(bool complete);

  const daqrun::DaqrunConfig& cfg_;
  const daqrun::RunReservation& run_;
  const hvmon::Mapping& mapping_;
  RunStatic st_;
  Clock& clock_;
  WriterStats& stats_;
  daqrun::FlagEngine flags_;
  daqrun::EventFlagger evtFlagger_;
  daqrun::TttClock ttt_;
  std::size_t segAnchorFirst_ = 0;
  std::deque<Pending> pending_;
  std::uint64_t nextSeq_ = 0;
  std::int64_t lastTttNs_ = -1;
  int N_;

  // Segmento actual.
  int seg_ = 0;
  std::unique_ptr<TFile> file_;
  TTree *events_ = nullptr, *hv_ = nullptr, *hvChanges_ = nullptr, *notes_ = nullptr, *monitor_ = nullptr;
  UtcNs segOpenUtc_ = 0, lastAutosaveUtc_ = 0;
  UtcNs segFirstUtc_ = 0, segLastUtc_ = 0;
  std::uint64_t segEvents_ = 0, segFirstSeq_ = 0, segLastSeq_ = 0;
  UtcNs runStartUtc_ = 0, runFirstEventUtc_ = 0, runLastEventUtc_ = 0;
  std::uint64_t runEvents_ = 0;
  nlohmann::json index_;
  bool finished_ = false;

  // Buffers de ramas.
  struct EvBuf {
    ULong64_t evtSeq = 0, tttNs = 0, hvSeq = 0;
    UInt_t counter = 0, tttRaw = 0, evtFlags = 0;
    Long64_t tUtc = 0, tRead = 0, tMono = 0, dtPrev = 0;
    UChar_t kind = 0;
    std::vector<UShort_t> adc;
    UShort_t chFlags[64] = {0};
  } ev_;
  struct HvBuf {
    Long64_t received = 0, tReq = 0, tResp = 0;
    ULong64_t stateSeq = 0, readSeq = 0;
    Bool_t fileOk = false, alive = false, stale = true, connected = false;
    Double_t v0set[12] = {0}, i0set[12] = {0}, vmon[12] = {0}, imon[12] = {0};
    UInt_t status[12] = {0}, pw[12] = {0};
    std::vector<std::string> retry;
  } hvb_;
  daqrun::HvSample lastHvWritten_;
  bool haveHvWritten_ = false;
  daqrun::HvSample lastHv_;
  bool haveHv_ = false;
  struct ChgBuf {
    Long64_t tDetect = 0, tPrev = 0, id = 0;
    Int_t hvCh = 0;
    std::string param, oldV, newV, origin;
  } chg_;
  struct NoteBuf {
    Long64_t t = 0;
    std::string user, host, text;
  } note_;
  struct MonBuf {
    Long64_t tStart = 0, tEnd = 0;
    ULong64_t nEvents = 0, monitorDropped = 0;
    Double_t rate = 0;
    UInt_t oot[64] = {0}, inTime[64] = {0};
    Double_t ootTime[64] = {0}, ootRate[64] = {0};
    Long64_t writerQ = 0, monitorQ = 0;
    ULong64_t chFlagEvents[16] = {0}, evtFlagEvents[9] = {0};
  } mon_;
  std::array<std::uint64_t, 16> periodCh_{};
  std::array<std::uint64_t, 9> periodEvt_{};
};

}  // namespace daq::io
