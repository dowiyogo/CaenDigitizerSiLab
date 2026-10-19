// replay_digitizer.cpp
#include "daq/dig/replay_digitizer.h"

#include <TFile.h>
#include <TTree.h>

#include <cmath>

#include "hodo_paths.h"

namespace daq::dig {

Result<ReplayParams> ReplayParams::fromJson(const nlohmann::json& j) {
  using R = Result<ReplayParams>;
  ReplayParams p;
  try {
    p.file = j.value("file", p.file);
    p.tree = j.value("tree", p.tree);
    p.firstEntry = j.value("first_entry", p.firstEntry);
    p.maxEntries = j.value("max_entries", p.maxEntries);
    p.speed = j.value("speed", p.speed);
    p.oldStartSample = j.value("old_start_sample", p.oldStartSample);
    p.dcOffset = j.value("dc_offset", p.dcOffset);
    p.tttWrapBits = j.value("ttt_wrap_bits", p.tttWrapBits);
  } catch (const std::exception& e) {
    return R::failure(-1, std::string("replay inválido: ") + e.what());
  }
  if (p.firstEntry < 0 || p.maxEntries < 0 || p.speed < 0 || p.oldStartSample < 0)
    return R::failure(-1, "replay: parámetros fuera de rango");
  return p;
}

nlohmann::json ReplayParams::toJson() const {
  return {{"file", file},           {"tree", tree},     {"first_entry", firstEntry},
          {"max_entries", maxEntries}, {"speed", speed}, {"old_start_sample", oldStartSample},
          {"dc_offset", dcOffset},  {"ttt_wrap_bits", tttWrapBits}};
}

ReplayDigitizer::ReplayDigitizer(ReplayParams p, Clock& clock) : p_(std::move(p)), clock_(clock) {}

ReplayDigitizer::~ReplayDigitizer() {
  if (file_) file_->Close();
  file_.reset();  // borra el árbol, que apunta a ch_
  // SetBranchAddress con puntero nulo: ROOT crea los std::vector<float> y el llamador los libera.
  for (auto* v : ch_) delete v;
}

Status ReplayDigitizer::configure(const DigitizerConfig& cfg) {
  cfg_ = cfg;
  configured_ = true;
  return open();
}

Status ReplayDigitizer::open() {
  path_ = hodo::resolveRawFile(p_.file);
  if (path_.empty()) return Status::failure(-1, "replay: nombre ambiguo " + p_.file);
  file_.reset(TFile::Open(path_.c_str(), "READ"));
  if (!file_ || file_->IsZombie()) return Status::failure(-1, "replay: no se pudo abrir " + path_ + " en modo READ");
  tree_ = dynamic_cast<TTree*>(file_->Get(p_.tree.c_str()));
  if (!tree_) return Status::failure(-1, "replay: no hay árbol " + p_.tree + " en " + path_);
  tree_->SetBranchStatus("*", 0);
  ch_.assign(kChannels, nullptr);
  for (int c = 0; c < kChannels; ++c) {
    const std::string b = "Ch" + std::to_string(c);
    tree_->SetBranchStatus(b.c_str(), 1);
    tree_->SetBranchAddress(b.c_str(), &ch_[static_cast<std::size_t>(c)]);
  }
  tree_->SetBranchStatus("event", 1);
  tree_->SetBranchAddress("event", &event_);
  tree_->SetBranchStatus("timestamp", 1);
  tree_->SetBranchAddress("timestamp", &timestamp_);
  tree_->SetCacheSize(64LL << 20);
  const long long n = tree_->GetEntries();
  next_ = std::min(p_.firstEntry, n);
  end_ = p_.maxEntries > 0 ? std::min(n, next_ + p_.maxEntries) : n;
  return okStatus();
}

Result<AcqOrigin> ReplayDigitizer::start() {
  if (!configured_ || !tree_) return Result<AcqOrigin>::failure(-1, "replay: start() sin configure()");
  if (p_.oldStartSample + cfg_.recordLength > 624)
    return Result<AcqOrigin>::failure(-1, "replay: la ventana sale del registro antiguo de 624 muestras");
  AcqOrigin o;
  o.utcBeforeStart = clock_.now_utc_ns();
  t0_ = o.utcBeforeStart;
  o.utcAfterStart = clock_.now_utc_ns();
  running_ = true;
  done_ = next_ >= end_;
  return o;
}

Status ReplayDigitizer::stop() {
  running_ = false;
  return okStatus();
}

Result<std::size_t> ReplayDigitizer::readBlock(EventBlock& out, std::size_t maxEvents) {
  out.tReadUtcNs = clock_.now_utc_ns();
  out.tReadMonoNs = monotonicNowNs();
  out.recordLength = cfg_.recordLength;
  if (!running_ || done_) return std::size_t{0};
  const AdcConversion conv{cfg_.nBit, cfg_.maxVpp, p_.dcOffset};
  const double maxAdc = static_cast<double>((1 << cfg_.nBit) - 1);
  const int n = cfg_.recordLength;
  std::size_t got = 0;
  while (got < maxEvents && next_ < end_) {
    // El evento ya cargado (havePending_) espera su momento si speed > 0.
    if (!havePending_) {
      if (tree_->GetEntry(next_) <= 0) return Result<std::size_t>::failure(-1, "replay: GetEntry falló en " + std::to_string(next_));
      havePending_ = true;
      if (ts0_ < 0) ts0_ = timestamp_;
    }
    const double elapsedUs = static_cast<double>(timestamp_ - ts0_);
    if (p_.speed > 0) {
      const UtcNs due = t0_ + static_cast<UtcNs>(elapsedUs * 1000.0 / p_.speed);
      if (due > out.tReadUtcNs) break;
    }
    Event ev;
    ev.counterHw = static_cast<std::uint32_t>(event_);
    // TTT al ritmo de entrega (intervalos originales / speed; con speed = 0, la hora de lectura):
    // así t_utc = origen + ttt queda coherente con t_read y el desenrollado con la hora del PC.
    const double deliveredNs =
        p_.speed > 0 ? elapsedUs * 1000.0 / p_.speed : static_cast<double>(out.tReadUtcNs - t0_);
    const std::uint64_t ticks = static_cast<std::uint64_t>(std::max(0.0, deliveredNs) / cfg_.tttTickNs);
    const std::uint64_t wrap = p_.tttWrapBits >= 32 ? (1ull << 32) : (1ull << p_.tttWrapBits);
    ev.tttRaw = static_cast<std::uint32_t>(ticks % wrap);
    ev.triggerKind = 0;
    ev.adc.resize(static_cast<std::size_t>(kChannels * n));
    for (int c = 0; c < kChannels; ++c) {
      const std::vector<float>& src = *ch_[static_cast<std::size_t>(c)];
      if (static_cast<int>(src.size()) < p_.oldStartSample + n)
        return Result<std::size_t>::failure(-1, "replay: entrada " + std::to_string(next_) + " con registro corto");
      for (int k = 0; k < n; ++k) {
        const double a = conv.toAdc(static_cast<double>(src[static_cast<std::size_t>(p_.oldStartSample + k)]));
        const double r = std::round(a);
        maxResidual_ = std::max(maxResidual_, std::fabs(a - r));
        ev.adc[static_cast<std::size_t>(c * n + k)] = static_cast<std::uint16_t>(std::min(maxAdc, std::max(0.0, r)));
      }
    }
    out.events.push_back(std::move(ev));
    havePending_ = false;
    ++next_;
    ++delivered_;
    ++got;
  }
  if (next_ >= end_) done_ = true;
  return got;
}

nlohmann::json ReplayDigitizer::info() const {
  return {{"model", "FakeDigitizer (replay)"},
          {"serial", "none"},
          {"roc_firmware", "unknown"},
          {"amc_firmware", "unknown"},
          {"caen_digitizer_version", "unknown"},
          {"caen_comm_version", "unknown"},
          {"caen_vme_version", "unknown"},
          {"replay", p_.toJson()},
          {"source_path", path_},
          {"max_rounding_residual_adc", maxResidual_},
          {"entries_delivered", delivered_}};
}

}  // namespace daq::dig
