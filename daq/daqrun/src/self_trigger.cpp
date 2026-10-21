// self_trigger.cpp
#include "daq/daqrun/self_trigger.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace daq::daqrun {

namespace {
double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v.empty() ? NAN : v[v.size() / 2];
}
}  // namespace

Status applySelfTriggerThresholds(dig::DigitizerConfig& cfg, const PeCalibration& cal) {
  auto& st = cfg.self;
  if (!st.enabled) return okStatus();
  if (static_cast<int>(st.thresholdPe.size()) != dig::kGroups) return Status::failure(-1, "auto-trigger: threshold_pe debe tener 8 valores");
  st.thresholdAdc.assign(dig::kGroups, 0);
  nlohmann::json conv = nlohmann::json::object();
  conv["method"] = "umbral ADC = round(mediana de las líneas base − T·mediana de G / LSB) por grupo; línea base = pedestal de la calibración con el DC offset del grupo";
  conv["calibration_version"] = cal.version;
  for (int g = 0; g < dig::kGroups; ++g) {
    if (!((st.groupMask >> g) & 1u)) continue;
    dig::AdcConversion ac{cfg.nBit, cfg.maxVpp, cfg.groupDcOffset[static_cast<std::size_t>(g)]};
    const std::uint32_t mask = st.channelMask(g);
    std::vector<double> bl, gain, stationGain;
    for (int c = 8 * g; c < 8 * g + 8; ++c) {
      const double gc = cal.gainMv[static_cast<std::size_t>(c)];
      if (!((mask >> (c % 8)) & 1u)) continue;
      bl.push_back(ac.toAdc(cal.pedestalMv[static_cast<std::size_t>(c)]));
      if (gc > 0) gain.push_back(gc);
    }
    for (int c = (g < 4 ? 0 : 32); c < (g < 4 ? 32 : 64); ++c)
      if (cal.gainMv[static_cast<std::size_t>(c)] > 0) stationGain.push_back(cal.gainMv[static_cast<std::size_t>(c)]);
    if (bl.empty()) return Status::failure(-1, "auto-trigger: el grupo " + std::to_string(g) + " no tiene canales habilitados");
    const double gMed = !gain.empty() ? median(gain) : median(stationGain);
    if (!(gMed > 0)) return Status::failure(-1, "auto-trigger: sin G para el grupo " + std::to_string(g) + " en la calibración");
    const double blMed = median(bl);
    const double T = st.thresholdPe[static_cast<std::size_t>(g)];
    const int thr = static_cast<int>(std::lround(blMed - T * gMed / ac.lsbMv()));
    if (thr < 0 || thr > 4095) return Status::failure(-1, "auto-trigger: umbral fuera de 0-4095 en el grupo " + std::to_string(g));
    st.thresholdAdc[static_cast<std::size_t>(g)] = thr;
    nlohmann::json eff = nlohmann::json::object();
    for (int c = 8 * g; c < 8 * g + 8; ++c) {
      if (!((mask >> (c % 8)) & 1u)) continue;
      const double gc = cal.gainMv[static_cast<std::size_t>(c)] > 0 ? cal.gainMv[static_cast<std::size_t>(c)] : gMed;
      eff[std::to_string(c)] = (ac.toAdc(cal.pedestalMv[static_cast<std::size_t>(c)]) - thr) * ac.lsbMv() / gc;
    }
    conv[std::to_string(g)] = {{"threshold_pe", T}, {"threshold_adc", thr}, {"baseline_median_adc", blMed}, {"G_median_mV", gMed},
                               {"channel_mask", mask}, {"effective_threshold_pe", eff}};
  }
  st.conversion = conv;
  return okStatus();
}

}  // namespace daq::daqrun
