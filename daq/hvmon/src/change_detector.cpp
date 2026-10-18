// change_detector.cpp
#include "daq/hvmon/change_detector.h"

#include <algorithm>

namespace daq::hvmon {

ChangeDetector::ChangeDetector(hv::StatusBitMap bits, std::vector<std::string> tripFlags)
    : bits_(std::move(bits)), tripFlags_(std::move(tripFlags)) {}

bool ChangeDetector::hasTrip(std::uint32_t status) const {
  return std::any_of(tripFlags_.begin(), tripFlags_.end(), [&](const std::string& f) { return bits_.has(status, f); });
}

void ChangeDetector::addPending(const PendingCommand& p) { pending_.push_back(p); }

void ChangeDetector::setPendingEnd(std::uint64_t id, SteadyTime tEnd) {
  for (auto& p : pending_)
    if (p.id == id) p.tEnd = tEnd;
}

void ChangeDetector::expire(SteadyTime now) {
  pending_.erase(std::remove_if(pending_.begin(), pending_.end(), [&](const PendingCommand& p) { return p.tEnd < now; }),
                 pending_.end());
}

namespace {
nlohmann::json valueJson(const hv::ParamValue& v) {
  return v.isNumeric() ? nlohmann::json(v.si) : nlohmann::json(v.raw);
}
bool sameValue(const hv::ParamValue& a, const hv::ParamValue& b) {
  return a.isNumeric() ? a.boardValue == b.boardValue : a.raw == b.raw;
}
}  // namespace

std::vector<ChangeRecord> ChangeDetector::compare(const hv::ChannelReading& prev, const hv::ChannelReading& cur,
                                                  const std::vector<std::string>& tracked, const std::string& readKind,
                                                  bool duringDisconnect, SteadyTime now,
                                                  std::uint64_t& nextId) const {
  std::vector<ChangeRecord> out;
  if (prev.params.empty() || prev.channels.empty()) return out;
  for (int ch : cur.channels) {
    const hv::ParamValue* statusNew = cur.get("Status", ch);
    const bool trip = statusNew && hasTrip(statusNew->raw);
    for (const auto& name : tracked) {
      const hv::ParamValue* a = prev.get(name, ch);
      const hv::ParamValue* b = cur.get(name, ch);
      if (!a || !b || sameValue(*a, *b)) continue;
      ChangeRecord r;
      r.id = nextId++;
      r.tUtcNs = cur.tRespUtcNs;
      r.tPrevReadUtcNs = prev.tReqUtcNs;
      r.slot = cur.slot;
      r.ch = ch;
      r.param = name;
      for (const auto& p : cur.params)
        if (p.name == name) r.unit = p.siUnit();
      r.oldValue = valueJson(*a);
      r.newValue = valueJson(*b);
      r.duringDisconnect = duringDisconnect;
      r.readKind = readKind;
      const bool pwOrStatus = name == "Pw" || name == "Status";
      r.origin = "external";
      if (pwOrStatus && trip) {
        r.origin = "trip";
      } else if (pwOrStatus) {
        for (const auto& p : pending_) {
          if (p.ch == ch && now >= p.tStart && now <= p.tEnd) {
            r.origin = "hvmon_command";
            r.commandId = p.id;
          }
        }
      }
      out.push_back(r);
    }
  }
  return out;
}

}  // namespace daq::hvmon
