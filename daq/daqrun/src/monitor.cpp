// monitor.cpp
#include "daq/daqrun/monitor.h"

#include <cmath>

#include "daq/pulse_finder.h"

namespace daq::daqrun {

namespace {
constexpr double kSampleS = 16e-9;
pulse::Region region(std::pair<int, int> p) { return pulse::Region{p.first, p.second}; }
}  // namespace

Monitor::Monitor(const DaqrunConfig& cfg, UtcNs originUtc, std::function<nlohmann::json()> queueLevels)
    : cfg_(cfg),
      ttt_(originUtc, cfg.digitizer.tttTickNs, cfg.tttWrapBits),
      queues_(std::move(queueLevels)),
      lsb_(dig::AdcConversion{cfg.digitizer.nBit, cfg.digitizer.maxVpp, 0x7FFF}.lsbMv()),
      origin_(originUtc) {
  v_.ampBinMv = lsb_;
  v_.ampBins = std::max(1, static_cast<int>(cfg.ampHistMaxMv / lsb_));
  v_.amp.assign(static_cast<std::size_t>(64 * v_.ampBins), 0);
  cur_.tStartUtcNs = originUtc;
}

void Monitor::process(const dig::EventBlock& blk) {
  const int n = blk.recordLength;
  const auto& w = cfg_.window;
  const std::vector<pulse::Region> oot{region(w.pedestal), region(w.control), region(w.post)};
  const std::vector<pulse::Region> inTime{region(w.signal)};
  const std::vector<pulse::Region> tail{region(w.tail)};
  const double ootBaseS = (oot[0].size() + oot[1].size() + oot[2].size()) * kSampleS;
  const double tailS = tail[0].size() * kSampleS;
  const double trimAdc = cfg_.trimMv / lsb_;
  const double lowest = 2.0 * lsb_;  // espectro: todo pulso sobre 2 LSB
  std::vector<pulse::Pulse> pulses;
  for (const auto& ev : blk.events) {
    const std::int64_t t = ttt_.tttNs(ev.tttRaw, blk.tReadUtcNs);
    const UtcNs tUtc = origin_ + t;
    if (lastTttNs_ >= 0) {
      const double dt = static_cast<double>(t - lastTttNs_) * 1e-9;
      if (dt > 0) {
        const double x = std::log10(dt);
        int b = static_cast<int>((x - MonitorView::kDtLogMin) / (MonitorView::kDtLogMax - MonitorView::kDtLogMin) * MonitorView::kDtBins);
        b = std::min(MonitorView::kDtBins - 1, std::max(0, b));
        ++v_.dtHist[static_cast<std::size_t>(b)];
      }
    }
    lastTttNs_ = t;
    if (v_.events == 0) v_.firstEventUtcNs = tUtc;
    v_.lastEventUtcNs = tUtc;
    ++v_.events;
    ++cur_.nEvents;
    for (int c = 0; c < dig::kChannels; ++c) {
      const std::uint16_t* x = ev.adc.data() + static_cast<std::size_t>(c * n);
      const int st = c < 32 ? 0 : 1;
      const double thr = cfg_.provisionalThresholdMv[st];
      const double bl = pulse::trimmedMean(x, region(w.pedestal), trimAdc);
      const int nIn = pulse::countPulses(x, n, bl, lsb_, thr, cfg_.minSeparationSamples, inTime);
      pulses.clear();
      pulse::findPulses(x, n, bl, lsb_, lowest, cfg_.minSeparationSamples, oot, pulses);
      double exposed = ootBaseS;
      if (nIn == 0) {
        pulse::findPulses(x, n, bl, lsb_, lowest, cfg_.minSeparationSamples, tail, pulses);
        exposed += tailS;
      }
      std::uint32_t above = 0;
      for (const auto& p : pulses) {
        if (p.amplitudeMv > thr) ++above;
        const int b = static_cast<int>(p.amplitudeMv / lsb_);
        if (b >= 0 && b < v_.ampBins) ++v_.amp[static_cast<std::size_t>(c * v_.ampBins + b)];
      }
      const auto ci = static_cast<std::size_t>(c);
      cur_.ootCounts[ci] += above;
      cur_.ootTimeS[ci] += exposed;
      cur_.lowCounts[ci] += static_cast<std::uint32_t>(pulses.size());
      cur_.inTimeCounts[ci] += static_cast<std::uint32_t>(nIn);
      v_.ootCounts[ci] += above;
      v_.ootTimeS[ci] += exposed;
      v_.inTimeCounts[ci] += static_cast<std::uint64_t>(nIn);
    }
  }
}

bool Monitor::takeSummary(UtcNs now, MonitorSummary& out, std::uint64_t monitorDropped) {
  if (static_cast<double>(now - cur_.tStartUtcNs) / 1e9 < cfg_.summaryPeriodS) return false;
  cur_.tEndUtcNs = now;
  const double span = static_cast<double>(now - cur_.tStartUtcNs) / 1e9;
  cur_.triggerRateHz = span > 0 ? static_cast<double>(cur_.nEvents) / span : 0;
  for (std::size_t c = 0; c < 64; ++c) cur_.ootRateHz[c] = cur_.ootTimeS[c] > 0 ? cur_.ootCounts[c] / cur_.ootTimeS[c] : 0;
  if (queues_) cur_.queues = queues_();
  cur_.monitorDropped = monitorDropped;
  out = cur_;
  v_.last = cur_;
  v_.rateHistory.emplace_back(now, cur_.triggerRateHz);
  if (v_.rateHistory.size() > 1440) v_.rateHistory.erase(v_.rateHistory.begin());
  v_.recent.push_back(cur_);
  if (v_.recent.size() > MonitorView::kRecent) v_.recent.erase(v_.recent.begin());
  cur_ = MonitorSummary{};
  cur_.tStartUtcNs = now;
  return true;
}

void Monitor::publishView() { view_.publish(v_); }

}  // namespace daq::daqrun
