// config.cpp
#include "daq/daqrun/config.h"

#include <cctype>
#include <fstream>

#include "daq/core/file_system.h"

namespace daq::daqrun {

namespace {

std::pair<int, int> range(const nlohmann::json& j, const char* key, std::pair<int, int> def) {
  if (!j.contains(key)) return def;
  return {j[key].at(0).get<int>(), j[key].at(1).get<int>()};
}

}  // namespace

bool validCampaignName(const std::string& s) {
  if (s.empty() || s.size() > 64 || !std::isalnum(static_cast<unsigned char>(s[0]))) return false;
  for (char ch : s) {
    const auto u = static_cast<unsigned char>(ch);
    if (!std::isalnum(u) && ch != '_' && ch != '.' && ch != '-') return false;
  }
  return s.find("..") == std::string::npos;
}

nlohmann::json WindowLayout::toJson() const {
  auto r = [](std::pair<int, int> p) { return nlohmann::json::array({p.first, p.second}); };
  return {{"ref_old_record_start_ns", refOldRecordStartNs},
          {"pedestal", r(pedestal)},
          {"control", r(control)},
          {"signal", r(signal)},
          {"tail", r(tail)},
          {"post", r(post)},
          {"sample_ns", 16}};
}

Result<DaqrunConfig> DaqrunConfig::fromJson(const nlohmann::json& jIn) {
  using R = Result<DaqrunConfig>;
  DaqrunConfig c;
  nlohmann::json j = jIn.is_object() ? jIn : nlohmann::json::object();
  try {
    const nlohmann::json dj = j.value("digitizer", nlohmann::json::object());
    auto d = dig::DigitizerConfig::fromJson(dj);
    if (!d) return R(d.error());
    c.digitizer = d.value();
    c.tttWrapBits = dj.value("ttt_wrap_bits", c.tttWrapBits);
    c.tttPcToleranceS = dj.value("ttt_pc_tolerance_s", c.tttPcToleranceS);
    c.dataFlags = DataFlagConfig::fromJson(j.value("data_flags", nlohmann::json::object()));
    c.panels = PanelConfig::fromJson(j.value("panels", nlohmann::json::object()));
    if (dj.contains("ttt_clock")) {
      const auto& k = dj["ttt_clock"];
      c.tttClockFit = k.value("fit", c.tttClockFit);
      c.tttClockWindowS = k.value("window_s", c.tttClockWindowS);
      c.tttClockReanchorS = k.value("reanchor_s", c.tttClockReanchorS);
      c.tttClockMinSpanS = k.value("min_span_s", c.tttClockMinSpanS);
      c.tttClockMinPoints = k.value("min_points", c.tttClockMinPoints);
      c.tttClockFirstFitS = k.value("first_fit_s", c.tttClockFirstFitS);
    }
    if (dj.contains("window_layout")) {
      const auto& w = dj["window_layout"];
      c.window.refOldRecordStartNs = w.value("ref_old_record_start_ns", c.window.refOldRecordStartNs);
      c.window.pedestal = range(w, "pedestal", c.window.pedestal);
      c.window.control = range(w, "control", c.window.control);
      c.window.signal = range(w, "signal", c.window.signal);
      c.window.tail = range(w, "tail", c.window.tail);
      c.window.post = range(w, "post", c.window.post);
    }
    // Campaña: con digitizador simulado, los valores por omisión son los de simulación.
    nlohmann::json cj = j.value("campaign", nlohmann::json::object());
    const bool sim = c.digitizer.backend != "caen";
    if (sim && !cj.contains("raw_root")) cj["raw_root"] = kSimRawRoot;
    if (sim && !cj.contains("manifests_dir")) cj["manifests_dir"] = kSimManifests;
    c.campaign = cj.value("name", c.campaign);
    c.rawRoot = cj.value("raw_root", c.rawRoot);
    c.manifestsDir = cj.value("manifests_dir", c.manifestsDir);
    c.operatorName = cj.value("operator", c.operatorName);
    const nlohmann::json rj = j.value("run", nlohmann::json::object());
    c.mode = rj.value("mode", c.mode);
    c.maxEvents = rj.value("max_events", c.maxEvents);
    c.maxDurationS = rj.value("max_duration_s", c.maxDurationS);
    c.segmentMaxS = rj.value("segment_max_s", c.segmentMaxS);
    c.segmentMaxBytes = rj.value("segment_max_bytes", c.segmentMaxBytes);
    c.autosaveS = rj.value("autosave_s", c.autosaveS);
    c.flagDelayS = rj.value("flag_delay_s", c.flagDelayS);
    c.compression = rj.value("compression", c.compression);
    const nlohmann::json hj = j.value("hv", nlohmann::json::object());
    c.hvStateFile = hj.value("state_file", c.hvStateFile);
    c.hvPollPeriodS = hj.value("poll_period_s", c.hvPollPeriodS);
    c.hvStaleAfterS = hj.value("stale_after_s", c.hvStaleAfterS);
    c.hvSnapshotPeriodS = hj.value("snapshot_period_s", c.hvSnapshotPeriodS);
    c.expectedOn = hj.value("expected_on", c.expectedOn);
    if (hj.contains("status_bits")) {
      auto sb = hv::StatusBitMap::fromJson(hj["status_bits"]);
      if (!sb) return R(sb.error());
      c.statusBits = sb.value();
    }
    c.tripFlags = hj.value("trip_flags", c.tripFlags);
    c.vmonToleranceV = hj.value("vmon_tolerance_V", c.vmonToleranceV);
    c.ovcKOverBaseline = hj.value("ovc_k_over_baseline", c.ovcKOverBaseline);
    c.imonBaselineWindowS = hj.value("imon_baseline_window_s", c.imonBaselineWindowS);
    c.tripBackdateS = hj.value("trip_backdate_s", c.tripBackdateS);
    c.notesDir = j.value("notes", nlohmann::json::object()).value("slowcontrol_dir", c.notesDir);
    c.mappingFile = j.value("mapping", nlohmann::json::object()).value("file", c.mappingFile);
    if (j.contains("calibration") && j["calibration"].contains("file") && j["calibration"]["file"].is_string())
      c.calibrationFile = j["calibration"]["file"].get<std::string>();
    const nlohmann::json pj = j.value("pulse_finder", nlohmann::json::object());
    c.minSeparationSamples = pj.value("min_separation_samples", c.minSeparationSamples);
    if (pj.contains("provisional_thresholds_mV")) {
      c.provisionalThresholdMv[0] = pj["provisional_thresholds_mV"].value("hod2018", c.provisionalThresholdMv[0]);
      c.provisionalThresholdMv[1] = pj["provisional_thresholds_mV"].value("hod2019", c.provisionalThresholdMv[1]);
    }
    c.trimMv = pj.value("trim_mV", c.trimMv);
    c.ampHistMaxMv = pj.value("amp_hist_max_mV", c.ampHistMaxMv);
    const nlohmann::json mj = j.value("monitoring", nlohmann::json::object());
    c.summaryPeriodS = mj.value("summary_period_s", c.summaryPeriodS);
    c.gapProbability = mj.value("gap_probability", c.gapProbability);
    c.burstProbability = mj.value("burst_probability", c.burstProbability);
    if (mj.contains("http")) {
      c.httpEnabled = mj["http"].value("enabled", c.httpEnabled);
      c.httpBind = mj["http"].value("bind", c.httpBind);
      c.httpPort = mj["http"].value("port", c.httpPort);
    }
    c.processRequestsHz = mj.value("process_requests_hz", c.processRequestsHz);
    c.monitorQueueCapacity = mj.value("monitor_queue_capacity", c.monitorQueueCapacity);
    const nlohmann::json qj = j.value("queues", nlohmann::json::object());
    c.writerQueueCapacity = qj.value("writer_capacity", c.writerQueueCapacity);
    c.hvQueueCapacity = qj.value("hv_capacity", c.hvQueueCapacity);
    c.writerShutdownTimeoutS = qj.value("writer_shutdown_timeout_s", c.writerShutdownTimeoutS);
    c.notifier = j.value("notifier", nlohmann::json::object());
    c.hvmonStaleAlertMin = c.notifier.value("hvmon_stale_alert_min", c.hvmonStaleAlertMin);
    c.noTriggerAlertMin = c.notifier.value("no_trigger_alert_min", c.noTriggerAlertMin);
    c.minTriggerRateHz = c.notifier.value("min_trigger_rate_Hz", c.minTriggerRateHz);
    c.triggerRateWindowMin = c.notifier.value("trigger_rate_window_min", c.triggerRateWindowMin);
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("configuración de daqrun inválida: ") + e.what());
  }
  // Validación.
  const int n = c.digitizer.recordLength;
  const std::pair<int, int>* regions[] = {&c.window.pedestal, &c.window.control, &c.window.signal, &c.window.tail,
                                          &c.window.post};
  int prevEnd = 0;
  for (const auto* r : regions) {
    if (r->first < prevEnd || r->second <= r->first || r->second > n)
      return R::failure(-1, "window_layout: regiones fuera de [0, N) o solapadas (N = " + std::to_string(n) + ")");
    prevEnd = r->second;
  }
  if (c.mode != "physics" && c.mode != "calibration") return R::failure(-1, "run.mode debe ser physics o calibration");
  if (c.writerQueueCapacity == 0 || c.monitorQueueCapacity == 0 || c.hvQueueCapacity == 0)
    return R::failure(-1, "queues: las capacidades deben ser > 0");
  if (c.hvmonStaleAlertMin <= 0 || c.noTriggerAlertMin <= 0 || c.minTriggerRateHz < 0 || c.triggerRateWindowMin <= 0)
    return R::failure(-1, "notifier.hvmon_stale_alert_min y no_trigger_alert_min deben ser > 0");
  if (c.segmentMaxS <= 0 || c.autosaveS <= 0 || c.flagDelayS < 0 || c.hvPollPeriodS <= 0 || c.hvStaleAfterS <= 0 ||
      c.hvSnapshotPeriodS <= 0 || c.summaryPeriodS <= 0)
    return R::failure(-1, "daqrun: períodos inválidos");
  if (!c.campaign.empty() && !validCampaignName(c.campaign))
    return R::failure(-1, "campaign.name inválido: " + c.campaign);
  for (const auto& f : c.tripFlags)
    if (c.statusBits.bitOf(f) < 0) return R::failure(-1, "hv.trip_flags: bit desconocido " + f);
  if (c.httpBind != "127.0.0.1" && c.httpBind != "localhost" && c.httpBind != "tailscale")
    return R::failure(-1, "monitoring.http.bind debe ser 127.0.0.1, localhost o tailscale");
  if (c.simulated()) {
    if (pathIsWithin(c.rawRoot, kRealRawRoot))
      return R::failure(-1, "campaign.raw_root = " + c.rawRoot + ": un run simulado (digitizador " + c.digitizer.backend +
                                ") no puede escribir en " + kRealRawRoot + " (usar " + kSimRawRoot + ")");
    if (pathIsWithin(c.manifestsDir, kRealManifests))
      return R::failure(-1, "campaign.manifests_dir = " + c.manifestsDir + ": un run simulado no puede escribir en " +
                                kRealManifests + " (usar " + kSimManifests + ")");
  }
  return c;
}

Result<DaqrunConfig> DaqrunConfig::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) return Result<DaqrunConfig>::failure(-1, "no se pudo abrir la configuración " + path);
  nlohmann::json j;
  try {
    in >> j;
  } catch (const std::exception& e) {
    return Result<DaqrunConfig>::failure(-1, "JSON inválido en " + path + ": " + e.what());
  }
  return fromJson(j);
}

nlohmann::json DaqrunConfig::toJson() const {
  nlohmann::json bits = nlohmann::json::object();
  for (const auto& b : statusBits.bits()) bits[b.first] = b.second;
  nlohmann::json dj = digitizer.toJson();
  dj["ttt_wrap_bits"] = tttWrapBits;
  dj["ttt_pc_tolerance_s"] = tttPcToleranceS;
  dj["ttt_clock"] = {{"fit", tttClockFit}, {"window_s", tttClockWindowS}, {"reanchor_s", tttClockReanchorS},
                     {"min_span_s", tttClockMinSpanS}, {"min_points", tttClockMinPoints}, {"first_fit_s", tttClockFirstFitS}};
  dj["window_layout"] = window.toJson();
  return {{"campaign", {{"name", campaign}, {"raw_root", rawRoot}, {"manifests_dir", manifestsDir}, {"operator", operatorName}}},
          {"digitizer", dj},
          {"run",
           {{"mode", mode},
            {"max_events", maxEvents},
            {"max_duration_s", maxDurationS},
            {"segment_max_s", segmentMaxS},
            {"segment_max_bytes", segmentMaxBytes},
            {"autosave_s", autosaveS},
            {"flag_delay_s", flagDelayS},
            {"compression", compression}}},
          {"hv",
           {{"state_file", hvStateFile},
            {"poll_period_s", hvPollPeriodS},
            {"stale_after_s", hvStaleAfterS},
            {"snapshot_period_s", hvSnapshotPeriodS},
            {"expected_on", expectedOn},
            {"status_bits", bits},
            {"trip_flags", tripFlags},
            {"vmon_tolerance_V", vmonToleranceV},
            {"ovc_k_over_baseline", ovcKOverBaseline},
            {"imon_baseline_window_s", imonBaselineWindowS},
            {"trip_backdate_s", tripBackdateS}}},
          {"notes", {{"slowcontrol_dir", notesDir}}},
          {"mapping", {{"file", mappingFile}}},
          {"calibration", {{"file", calibrationFile.empty() ? nlohmann::json() : nlohmann::json(calibrationFile)}}},
          {"pulse_finder",
           {{"min_separation_samples", minSeparationSamples},
            {"provisional_thresholds_mV", {{"hod2018", provisionalThresholdMv[0]}, {"hod2019", provisionalThresholdMv[1]}}},
            {"trim_mV", trimMv},
            {"amp_hist_max_mV", ampHistMaxMv}}},
          {"monitoring",
           {{"summary_period_s", summaryPeriodS},
            {"gap_probability", gapProbability},
            {"burst_probability", burstProbability},
            {"http", {{"enabled", httpEnabled}, {"bind", httpBind}, {"port", httpPort}}},
            {"process_requests_hz", processRequestsHz},
            {"monitor_queue_capacity", monitorQueueCapacity}}},
          {"queues",
           {{"writer_capacity", writerQueueCapacity},
            {"hv_capacity", hvQueueCapacity},
            {"writer_shutdown_timeout_s", writerShutdownTimeoutS}}},
          {"notifier", [&] {
             nlohmann::json n = notifier;
             n["hvmon_stale_alert_min"] = hvmonStaleAlertMin;
             n["no_trigger_alert_min"] = noTriggerAlertMin;
             n["min_trigger_rate_Hz"] = minTriggerRateHz;
             n["trigger_rate_window_min"] = triggerRateWindowMin;
             return n;
           }()},
          {"data_flags", dataFlags.toJson()},
          {"panels", panels.toJson()}};
}

}  // namespace daq::daqrun
