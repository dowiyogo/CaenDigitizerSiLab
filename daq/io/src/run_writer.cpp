// run_writer.cpp
#include "daq/io/run_writer.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <TFile.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>

#include <cerrno>
#include <cstring>

#include "daq/core/file_system.h"
#include "daq/core/logger.h"
#include "daq/core/sha256.h"
#include "daq/daqrun/notes.h"

namespace daq::io {

namespace {

constexpr double kNs = 1e9;

Status fsyncPath(const std::string& path) {
  const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
  if (fd < 0) return Status::failure(errno, "open " + path + ": " + std::strerror(errno));
  const int r = ::fsync(fd);
  const int e = errno;
  ::close(fd);
  if (r != 0) return Status::failure(e, "fsync " + path + ": " + std::strerror(e));
  return okStatus();
}

}  // namespace

std::string WriterStats::currentSegment() const {
  std::lock_guard<std::mutex> lk(m_);
  return seg_;
}
void WriterStats::setCurrentSegment(const std::string& s) {
  std::lock_guard<std::mutex> lk(m_);
  seg_ = s;
}

nlohmann::json WriterStats::toJson() const {
  nlohmann::json ch = nlohmann::json::array(), ev = nlohmann::json::array();
  for (const auto& a : chFlagEvents) ch.push_back(a.load());
  for (const auto& a : evtFlagEvents) ev.push_back(a.load());
  return {{"events_written", eventsWritten.load()},
          {"segments_closed", segmentsClosed.load()},
          {"pending_delay", pendingDelay.load()},
          {"max_pending_delay", maxPendingDelay.load()},
          {"hv_rows", hvRows.load()},
          {"hv_change_rows", hvChangeRows.load()},
          {"bytes_closed", bytesClosed.load()},
          {"autosaves", autosaves.load()},
          {"last_evt_seq", lastEvtSeq.load()},
          {"current_segment", currentSegment()},
          {"ch_flag_events", ch},
          {"evt_flag_events", ev}};
}

RunWriter::RunWriter(const daqrun::DaqrunConfig& cfg, const daqrun::RunReservation& run, const hvmon::Mapping& mapping,
                     RunStatic st, Clock& clock, WriterStats& stats)
    : cfg_(cfg),
      run_(run),
      mapping_(mapping),
      st_(std::move(st)),
      clock_(clock),
      stats_(stats),
      flags_(cfg, mapping),
      evtFlagger_(cfg.gapProbability, cfg.burstProbability, cfg.tttPcToleranceS),
      ttt_(st_.origin.utcBeforeStart, cfg.digitizer.tttTickNs, cfg.tttWrapBits),
      N_(cfg.digitizer.recordLength) {
  ev_.adc.resize(static_cast<std::size_t>(64 * N_));
  runStartUtc_ = st_.origin.utcBeforeStart;
  index_ = {{"schema", "daqv2.run_index/1"},
            {"campaign", run_.location().name},
            {"run", run_.run()},
            {"base_name", run_.baseName()},
            {"mode", cfg_.mode},
            {"simulated", cfg_.simulated()},
            {"digitizer_backend", cfg_.digitizer.backend},
            {"start_utc_ns", runStartUtc_},
            {"start_utc", formatUtcIso(runStartUtc_)},
            {"end_utc_ns", 0},
            {"complete", false},
            {"config_sha256", st_.configSha256},
            {"mapping_version", mapping_.version},
            {"mapping_sha256", mapping_.sha256},
            {"calibration_version", cfg_.calibrationFile.empty() ? "none" : cfg_.calibrationFile},
            {"record_length", N_},
            {"segments", nlohmann::json::array()}};
}

RunWriter::~RunWriter() {
  if (file_ && !finished_) {
    logError("escritor destruido sin finish(): se cierra el segmento ", seg_, " como no limpio");
    closeSegment(false);
  }
}

std::string RunWriter::segmentPath(int seg) const { return run_.location().dir() + "/" + run_.segmentFile(seg); }
std::string RunWriter::indexPath() const { return run_.location().dir() + "/" + run_.indexFile(); }

Status RunWriter::open() { return openSegment(); }

Status RunWriter::openSegment() {
  ++seg_;
  const std::string path = segmentPath(seg_);
  if (daqrun::fileNameUsedInAnyCampaign(run_.location().rawRoot, run_.segmentFile(seg_)))
    return Status::failure(-1, "el segmento " + run_.segmentFile(seg_) + " ya existe");
  file_.reset(TFile::Open(path.c_str(), "CREATE", "", cfg_.compression));
  if (!file_ || file_->IsZombie()) return Status::failure(-1, "no se pudo crear " + path);
  file_->cd();
  events_ = new TTree("events", "eventos del digitizador (adc en cuentas ADC)");
  hv_ = new TTree("hv", "estado de hvmon (cada hv_snapshot_period_s y en cada cambio); SI: V, A");
  hvChanges_ = new TTree("hv_changes", "cambios detectados por hvmon (dedup por id)");
  notes_ = new TTree("notes", "notas del operador del intervalo del segmento");
  monitor_ = new TTree("monitor", "resúmenes del monitor por período (umbrales provisionales en mV)");
  for (TTree* t : {events_, hv_, hvChanges_, notes_, monitor_}) {
    t->SetDirectory(file_.get());
    t->SetAutoSave(0);  // AutoSave por tiempo (autosave_s), no por bytes
  }
  const std::string adcLeaf = "adc[64][" + std::to_string(N_) + "]/s";
  events_->Branch("evt_seq", &ev_.evtSeq, "evt_seq/l");
  events_->Branch("evt_counter_hw", &ev_.counter, "evt_counter_hw/i");
  events_->Branch("ttt_raw", &ev_.tttRaw, "ttt_raw/i");
  events_->Branch("ttt_ns", &ev_.tttNs, "ttt_ns/l");
  events_->Branch("t_utc_ns", &ev_.tUtc, "t_utc_ns/L");
  events_->Branch("t_read_utc_ns", &ev_.tRead, "t_read_utc_ns/L");
  events_->Branch("t_read_mono_ns", &ev_.tMono, "t_read_mono_ns/L");
  events_->Branch("dt_prev_ns", &ev_.dtPrev, "dt_prev_ns/L");
  events_->Branch("trigger_kind", &ev_.kind, "trigger_kind/b");
  events_->Branch("adc", ev_.adc.data(), adcLeaf.c_str());
  events_->Branch("ch_flags", ev_.chFlags, "ch_flags[64]/s");
  events_->Branch("evt_flags", &ev_.evtFlags, "evt_flags/i");
  events_->Branch("hv_seq", &ev_.hvSeq, "hv_seq/l");

  hv_->Branch("received_utc_ns", &hvb_.received, "received_utc_ns/L");
  hv_->Branch("state_seq", &hvb_.stateSeq, "state_seq/l");
  hv_->Branch("read_seq", &hvb_.readSeq, "read_seq/l");
  hv_->Branch("t_req_utc_ns", &hvb_.tReq, "t_req_utc_ns/L");
  hv_->Branch("t_resp_utc_ns", &hvb_.tResp, "t_resp_utc_ns/L");
  hv_->Branch("file_ok", &hvb_.fileOk, "file_ok/O");
  hv_->Branch("alive", &hvb_.alive, "alive/O");
  hv_->Branch("stale", &hvb_.stale, "stale/O");
  hv_->Branch("connected", &hvb_.connected, "connected/O");
  hv_->Branch("v0set_V", hvb_.v0set, "v0set_V[12]/D");
  hv_->Branch("i0set_A", hvb_.i0set, "i0set_A[12]/D");
  hv_->Branch("vmon_V", hvb_.vmon, "vmon_V[12]/D");
  hv_->Branch("imon_A", hvb_.imon, "imon_A[12]/D");
  hv_->Branch("status", hvb_.status, "status[12]/i");
  hv_->Branch("pw", hvb_.pw, "pw[12]/i");
  hv_->Branch("retry_state", &hvb_.retry);

  hvChanges_->Branch("t_detect_utc_ns", &chg_.tDetect, "t_detect_utc_ns/L");
  hvChanges_->Branch("t_prev_read_utc_ns", &chg_.tPrev, "t_prev_read_utc_ns/L");
  hvChanges_->Branch("change_id", &chg_.id, "change_id/L");
  hvChanges_->Branch("hv_ch", &chg_.hvCh, "hv_ch/I");
  hvChanges_->Branch("param", &chg_.param);
  hvChanges_->Branch("old_value", &chg_.oldV);
  hvChanges_->Branch("new_value", &chg_.newV);
  hvChanges_->Branch("origin", &chg_.origin);

  notes_->Branch("t_utc_ns", &note_.t, "t_utc_ns/L");
  notes_->Branch("user", &note_.user);
  notes_->Branch("host", &note_.host);
  notes_->Branch("text", &note_.text);

  monitor_->Branch("t_start_utc_ns", &mon_.tStart, "t_start_utc_ns/L");
  monitor_->Branch("t_end_utc_ns", &mon_.tEnd, "t_end_utc_ns/L");
  monitor_->Branch("n_events", &mon_.nEvents, "n_events/l");
  monitor_->Branch("trigger_rate_Hz", &mon_.rate, "trigger_rate_Hz/D");
  monitor_->Branch("oot_counts", mon_.oot, "oot_counts[64]/i");
  monitor_->Branch("oot_time_s", mon_.ootTime, "oot_time_s[64]/D");
  monitor_->Branch("oot_rate_Hz", mon_.ootRate, "oot_rate_Hz[64]/D");
  monitor_->Branch("in_time_counts", mon_.inTime, "in_time_counts[64]/i");
  monitor_->Branch("writer_queue", &mon_.writerQ, "writer_queue/L");
  monitor_->Branch("monitor_queue", &mon_.monitorQ, "monitor_queue/L");
  monitor_->Branch("monitor_dropped", &mon_.monitorDropped, "monitor_dropped/l");
  monitor_->Branch("ch_flag_events", mon_.chFlagEvents, "ch_flag_events[16]/l");
  monitor_->Branch("evt_flag_events", mon_.evtFlagEvents, "evt_flag_events[9]/l");

  segOpenUtc_ = lastAutosaveUtc_ = clock_.now_utc_ns();
  segEvents_ = 0;
  segFirstUtc_ = segLastUtc_ = 0;
  if (haveHv_) fillHv(lastHv_);  // estado HV al abrir el segmento
  stats_.setCurrentSegment(run_.segmentFile(seg_));
  logInfo("segmento abierto: ", path);
  return okStatus();
}

nlohmann::json RunWriter::runInfo(bool final) const {
  const dig::AdcConversion conv{cfg_.digitizer.nBit, cfg_.digitizer.maxVpp, cfg_.digitizer.groupDcOffset[0]};
  nlohmann::json adc = conv.toJson();
  adc["group_dc_offset_written"] = cfg_.digitizer.groupDcOffset;
  adc["group_dc_offset_read"] = cfg_.digitizer.groupDcOffset;  // fake: sin lectura de la placa
  return {{"schema", "daqv2.run_info/1"},
          {"campaign", run_.location().name},
          {"run", run_.run()},
          {"segment", seg_},
          {"mode", cfg_.mode},
          {"simulated", cfg_.simulated()},
          {"host", st_.host},
          {"git_commit", st_.gitCommit},
          {"git_dirty", st_.gitDirty},
          {"versions",
           {{"root", st_.rootVersion},
            {"caen_digitizer", st_.digitizerInfo.value("caen_digitizer_version", "unknown")},
            {"caen_comm", st_.digitizerInfo.value("caen_comm_version", "unknown")},
            {"caen_vme", st_.digitizerInfo.value("caen_vme_version", "unknown")}}},
          {"digitizer", st_.digitizerInfo},
          {"digitizer_final", st_.digitizerInfoFinal},
          {"ttt_origin", {{"utc_before_start_ns", st_.origin.utcBeforeStart}, {"utc_after_start_ns", st_.origin.utcAfterStart}}},
          {"ttt_tick_ns", cfg_.digitizer.tttTickNs},
          {"ttt_wrap_bits", cfg_.tttWrapBits},
          {"run_start_utc_ns", runStartUtc_},
          {"segment_open_utc_ns", segOpenUtc_},
          {"segment_first_event_utc_ns", segFirstUtc_},
          {"segment_last_event_utc_ns", segLastUtc_},
          {"segment_events", segEvents_},
          {"run_events_so_far", runEvents_},
          {"final_segment", final},
          {"record_length", N_},
          {"post_trigger_register_written", cfg_.digitizer.postTriggerRegister},
          {"post_trigger_register_read", cfg_.digitizer.postTriggerRegister},
          {"window_layout", cfg_.window.toJson()},
          {"ref_old_record_start_ns", cfg_.window.refOldRecordStartNs},
          {"adc_conversion", adc},
          {"calibration_state", cfg_.calibrationFile.empty() ? "none" : "loaded"},
          {"flag_delay_s", cfg_.flagDelayS},
          {"mapping_version", mapping_.version},
          {"mapping_sha256", mapping_.sha256}};
}

Status RunWriter::closeSegment(bool clean) {
  if (!file_) return okStatus();
  const UtcNs closeUtc = clock_.now_utc_ns();
  // Notas del intervalo del segmento.
  for (const auto& n : daqrun::readNotes(cfg_.notesDir, segOpenUtc_, closeUtc + 1)) {
    note_.t = n.tUtcNs;
    note_.user = n.user;
    note_.host = n.host;
    note_.text = n.text;
    notes_->Fill();
  }
  file_->cd();
  const bool final = finished_;
  TNamed("config_json", st_.configJson.dump().c_str()).Write();
  TNamed("config_sha256", st_.configSha256.c_str()).Write();
  TNamed("mapping_json", st_.mappingText.c_str()).Write();
  TNamed("mapping_version", mapping_.version.c_str()).Write();
  TNamed("mapping_sha256", mapping_.sha256.c_str()).Write();
  TNamed("calibration_json", cfg_.calibrationFile.empty() ? "none" : cfg_.calibrationFile.c_str()).Write();
  TNamed("calibration_version", cfg_.calibrationFile.empty() ? "none" : cfg_.calibrationFile.c_str()).Write();
  TNamed("run_info_json", runInfo(final).dump().c_str()).Write();
  TNamed("hvmon_state_at_start", st_.hvAtStart.dump().c_str()).Write();
  for (TTree* t : {events_, hv_, hvChanges_, notes_, monitor_}) t->Write("", TObject::kOverwrite);
  file_->Close();
  file_.reset();  // borra los árboles (pertenecen al archivo)
  events_ = hv_ = hvChanges_ = notes_ = monitor_ = nullptr;
  const std::string path = segmentPath(seg_);
  if (auto s = fsyncPath(path); !s) logError(s.error().message);
  if (::chmod(path.c_str(), 0444) != 0) logError("chmod a-w ", path, ": ", std::strerror(errno));
  auto sha = sha256File(path);
  if (!sha) return Status(sha.error());
  if (auto s = daqrun::appendSha256Sum(run_.location().dir(), run_.segmentFile(seg_), sha.value()); !s) return s;
  struct stat stt {};
  ::stat(path.c_str(), &stt);
  stats_.bytesClosed += static_cast<std::uint64_t>(stt.st_size);
  ++stats_.segmentsClosed;
  index_["segments"].push_back({{"seg", seg_},
                                {"file", run_.segmentFile(seg_)},
                                {"sha256", sha.value()},
                                {"bytes", static_cast<std::uint64_t>(stt.st_size)},
                                {"open_utc_ns", segOpenUtc_},
                                {"close_utc_ns", closeUtc},
                                {"first_utc_ns", segFirstUtc_},
                                {"last_utc_ns", segLastUtc_},
                                {"n_events", segEvents_},
                                {"first_evt_seq", segEvents_ ? segFirstSeq_ : 0},
                                {"last_evt_seq", segEvents_ ? segLastSeq_ : 0},
                                {"closed_cleanly", clean}});
  logInfo("segmento cerrado: ", path, " (", segEvents_, " eventos, ", stt.st_size, " bytes, sha256 ", sha.value(), ")");
  return writeIndex(false);
}

Status RunWriter::writeIndex(bool complete) {
  index_["complete"] = complete;
  index_["end_utc_ns"] = runLastEventUtc_;
  index_["end_utc"] = runLastEventUtc_ ? formatUtcIso(runLastEventUtc_) : "";
  index_["first_event_utc_ns"] = runFirstEventUtc_;
  index_["n_events"] = runEvents_;
  index_["updated_utc_ns"] = clock_.now_utc_ns();
  PosixFileSystem fs;
  return fs.writeAtomic(indexPath(), index_.dump(2) + "\n", 0644);
}

Status RunWriter::addBlock(const dig::EventBlockPtr& b) {
  for (std::size_t i = 0; i < b->events.size(); ++i) {
    const auto& e = b->events[i];
    Pending p;
    p.blk = b;
    p.idx = i;
    p.seq = nextSeq_++;
    p.tttNs = ttt_.tttNs(e.tttRaw, b->tReadUtcNs);
    p.tUtc = st_.origin.utcBeforeStart + p.tttNs;
    p.dtPrev = lastTttNs_ >= 0 ? p.tttNs - lastTttNs_ : 0;
    lastTttNs_ = p.tttNs;
    p.evtFlags = evtFlagger_.flags(p.tUtc, b->tReadUtcNs, e.counterHw, b->boardFull && i == 0, p.dtPrev);
    if (cfg_.calibrationFile.empty()) p.evtFlags |= daqrun::evtflag::NoCalibration;
    pending_.push_back(std::move(p));
  }
  stats_.pendingDelay = pending_.size();
  if (pending_.size() > stats_.maxPendingDelay.load()) stats_.maxPendingDelay = pending_.size();
  return okStatus();
}

bool RunWriter::hvChanged(const daqrun::HvSample& s) const {
  if (!haveHvWritten_) return true;
  const auto& a = lastHvWritten_;
  if (a.alive != s.alive || a.stale != s.stale || a.connected != s.connected || a.fileOk != s.fileOk) return true;
  for (int h = 0; h < daqrun::kHvChannels; ++h) {
    const auto& x = a.ch[static_cast<std::size_t>(h)];
    const auto& y = s.ch[static_cast<std::size_t>(h)];
    if (x.pw != y.pw || x.status != y.status || x.v0setV != y.v0setV || x.i0setA != y.i0setA || x.retryState != y.retryState)
      return true;
  }
  return false;
}

void RunWriter::fillHv(const daqrun::HvSample& s) {
  if (!hv_) return;
  hvb_.received = s.receivedUtcNs;
  hvb_.stateSeq = s.stateSeq;
  hvb_.readSeq = s.readSeq;
  hvb_.tReq = s.tReqUtcNs;
  hvb_.tResp = s.tRespUtcNs;
  hvb_.fileOk = s.fileOk;
  hvb_.alive = s.alive;
  hvb_.stale = s.stale;
  hvb_.connected = s.connected;
  hvb_.retry.assign(daqrun::kHvChannels, "");
  for (int h = 0; h < daqrun::kHvChannels; ++h) {
    const auto& c = s.ch[static_cast<std::size_t>(h)];
    hvb_.v0set[h] = c.v0setV;
    hvb_.i0set[h] = c.i0setA;
    hvb_.vmon[h] = c.vmonV;
    hvb_.imon[h] = c.imonA;
    hvb_.status[h] = c.status;
    hvb_.pw[h] = c.pw;
    hvb_.retry[static_cast<std::size_t>(h)] = c.retryState;
  }
  hv_->Fill();
  ++stats_.hvRows;
  lastHvWritten_ = s;
  haveHvWritten_ = true;
}

void RunWriter::addHv(const daqrun::HvSample& s) {
  flags_.addSample(s);
  if (!haveHvWritten_ || hvChanged(s) ||
      static_cast<double>(s.receivedUtcNs - lastHvWritten_.receivedUtcNs) / kNs >= cfg_.hvSnapshotPeriodS)
    fillHv(s);
  for (const auto& c : s.newChanges) {
    if (!hvChanges_) break;
    chg_.tDetect = c.value("t_utc_ns", Long64_t{0});
    chg_.tPrev = c.value("t_prev_read_utc_ns", Long64_t{0});
    chg_.id = c.value("id", Long64_t{-1});
    chg_.hvCh = c.value("ch", -1);
    chg_.param = c.value("param", "");
    chg_.oldV = c.contains("old") ? c["old"].dump() : "";
    chg_.newV = c.contains("new") ? c["new"].dump() : "";
    chg_.origin = c.value("origin", "");
    hvChanges_->Fill();
    ++stats_.hvChangeRows;
  }
  lastHv_ = s;
  lastHv_.newChanges = nlohmann::json::array();
  haveHv_ = true;
}

void RunWriter::addSummary(const daqrun::MonitorSummary& m) {
  if (!monitor_) return;
  mon_.tStart = m.tStartUtcNs;
  mon_.tEnd = m.tEndUtcNs;
  mon_.nEvents = m.nEvents;
  mon_.rate = m.triggerRateHz;
  for (int c = 0; c < 64; ++c) {
    mon_.oot[c] = m.ootCounts[static_cast<std::size_t>(c)];
    mon_.inTime[c] = m.inTimeCounts[static_cast<std::size_t>(c)];
    mon_.ootTime[c] = m.ootTimeS[static_cast<std::size_t>(c)];
    mon_.ootRate[c] = m.ootRateHz[static_cast<std::size_t>(c)];
  }
  mon_.writerQ = m.queues.value("writer", -1);
  mon_.monitorQ = m.queues.value("monitor", -1);
  mon_.monitorDropped = m.monitorDropped;
  for (int b = 0; b < 16; ++b) mon_.chFlagEvents[b] = periodCh_[static_cast<std::size_t>(b)];
  for (int b = 0; b < 9; ++b) mon_.evtFlagEvents[b] = periodEvt_[static_cast<std::size_t>(b)];
  periodCh_.fill(0);
  periodEvt_.fill(0);
  monitor_->Fill();
}

Status RunWriter::writeEvent(const Pending& p) {
  const dig::Event& e = p.blk->events[p.idx];
  const daqrun::ChannelFlags cf = flags_.flagsFor(p.tUtc);
  ev_.evtSeq = p.seq;
  ev_.counter = e.counterHw;
  ev_.tttRaw = e.tttRaw;
  ev_.tttNs = static_cast<ULong64_t>(p.tttNs);
  ev_.tUtc = p.tUtc;
  ev_.tRead = p.blk->tReadUtcNs;
  ev_.tMono = p.blk->tReadMonoNs;
  ev_.dtPrev = p.dtPrev;
  ev_.kind = e.triggerKind;
  if (static_cast<int>(e.adc.size()) != 64 * N_) return Status::failure(-1, "evento con tamaño de registro distinto de N");
  std::copy(e.adc.begin(), e.adc.end(), ev_.adc.begin());
  std::uint16_t anyCh = 0;
  for (int c = 0; c < 64; ++c) {
    ev_.chFlags[c] = cf.ch[static_cast<std::size_t>(c)];
    anyCh = static_cast<std::uint16_t>(anyCh | cf.ch[static_cast<std::size_t>(c)]);
  }
  ev_.evtFlags = p.evtFlags | (cf.anyHv ? daqrun::evtflag::HvAnyFlag : 0u);
  ev_.hvSeq = cf.hvSeq;
  events_->Fill();
  for (int b = 0; b < 16; ++b)
    if (anyCh & (1u << b)) {
      ++stats_.chFlagEvents[static_cast<std::size_t>(b)];
      ++periodCh_[static_cast<std::size_t>(b)];
    }
  for (int b = 0; b < 9; ++b)
    if (ev_.evtFlags & (1u << b)) {
      ++stats_.evtFlagEvents[static_cast<std::size_t>(b)];
      ++periodEvt_[static_cast<std::size_t>(b)];
    }
  if (segEvents_ == 0) {
    segFirstUtc_ = p.tUtc;
    segFirstSeq_ = p.seq;
  }
  segLastUtc_ = p.tUtc;
  segLastSeq_ = p.seq;
  ++segEvents_;
  if (runEvents_ == 0) runFirstEventUtc_ = p.tUtc;
  runLastEventUtc_ = p.tUtc;
  ++runEvents_;
  ++stats_.eventsWritten;
  stats_.lastEvtSeq = p.seq;
  return okStatus();
}

Status RunWriter::release(bool all) {
  const UtcNs now = clock_.now_utc_ns();
  const auto delay = static_cast<UtcNs>(cfg_.flagDelayS * kNs);
  while (!pending_.empty() && (all || pending_.front().tUtc + delay <= now)) {
    if (auto s = writeEvent(pending_.front()); !s) return s;
    pending_.pop_front();
  }
  stats_.pendingDelay = pending_.size();
  const UtcNs oldest = pending_.empty() ? now : pending_.front().tUtc;
  flags_.prune(oldest - static_cast<UtcNs>((cfg_.tripBackdateS + 2 * cfg_.hvStaleAfterS + 60.0) * kNs));
  return okStatus();
}

Status RunWriter::tick(bool allowRotate) {
  if (!file_) return Status::failure(-1, "escritor sin segmento abierto");
  if (auto s = release(false); !s) return s;
  const UtcNs now = clock_.now_utc_ns();
  const bool byTime = static_cast<double>(now - segOpenUtc_) / kNs >= cfg_.segmentMaxS;
  const bool byBytes = file_->GetBytesWritten() >= static_cast<Long64_t>(cfg_.segmentMaxBytes);
  if (allowRotate && (byTime || byBytes)) {
    if (auto s = closeSegment(true); !s) return s;
    return openSegment();
  }
  if (static_cast<double>(now - lastAutosaveUtc_) / kNs >= cfg_.autosaveS) {
    for (TTree* t : {events_, hv_, hvChanges_, notes_, monitor_}) t->AutoSave("SaveSelf");
    lastAutosaveUtc_ = now;
    ++stats_.autosaves;
  }
  return okStatus();
}

Status RunWriter::finish(bool cleanShutdown) {
  if (finished_) return okStatus();
  Status st = release(true);
  finished_ = true;
  if (auto s = closeSegment(cleanShutdown); !s && st) st = s;
  index_["closed_cleanly"] = cleanShutdown;
  if (auto s = writeIndex(true); !s && st) st = s;
  const std::string ip = indexPath();
  if (::chmod(ip.c_str(), 0444) != 0) logError("chmod a-w ", ip, ": ", std::strerror(errno));
  auto sha = sha256File(ip);
  if (sha) {
    if (auto s = daqrun::appendSha256Sum(run_.location().dir(), run_.indexFile(), sha.value()); !s && st) st = s;
  }
  const nlohmann::json entry{{"run", run_.run()},
                             {"index", run_.indexFile()},
                             {"index_sha256", sha ? sha.value() : ""},
                             {"start_utc_ns", runStartUtc_},
                             {"end_utc_ns", runLastEventUtc_},
                             {"n_events", runEvents_},
                             {"segments", index_["segments"].size()},
                             {"complete", true},
                             {"closed_cleanly", cleanShutdown},
                             {"mapping_version", mapping_.version},
                             {"calibration_version", cfg_.calibrationFile.empty() ? "none" : cfg_.calibrationFile}};
  if (auto s = daqrun::updateManifestRun(run_.location(), entry); !s && st) st = s;
  stats_.setCurrentSegment("");
  return st;
}

}  // namespace daq::io
