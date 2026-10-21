// panels.cpp
#include "daq/daqrun/panels.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "daq/daqrun/config.h"
#include "daq/daqrun/monitor.h"
#include "daq/gain_fit.h"

namespace daq::daqrun {

namespace {
constexpr double kSampleNs = 16.0;
const char* kPlaneName[4] = {"hod2018_X", "hod2018_Y", "hod2019_X", "hod2019_Y"};
}  // namespace

void RollingHist::init(std::string name, std::string title, int nx, double x0, double x1, int ny, double y0, double y1,
                       std::size_t slices) {
  h_ = PanelHist{std::move(name), std::move(title), nx, ny, x0, x1, y0, y1, {}, {}};
  const auto n = static_cast<std::size_t>(nx * ny);
  h_.cum.assign(n, 0);
  h_.hour.assign(n, 0);
  cur_.assign(n, 0);
  slices_ = slices;
  old_.clear();
}

void RollingHist::fillBin(int ix, int iy) {
  if (ix < 0 || ix >= h_.nx || iy < 0 || iy >= h_.ny) return;
  const auto i = static_cast<std::size_t>(iy * h_.nx + ix);
  ++h_.cum[i];
  ++h_.hour[i];
  ++cur_[i];
}

void RollingHist::fill(double x, double y) {
  const int ix = static_cast<int>(std::floor((x - h_.x0) / (h_.x1 - h_.x0) * h_.nx));
  const int iy = h_.ny == 1 ? 0 : static_cast<int>(std::floor((y - h_.y0) / (h_.y1 - h_.y0) * h_.ny));
  fillBin(std::min(ix, h_.nx - 1), std::min(iy, h_.ny - 1));  // el último bin acumula el desborde
}

void RollingHist::closeSlice() {
  old_.push_back(cur_);
  std::fill(cur_.begin(), cur_.end(), 0u);
  if (old_.size() >= slices_) {
    const auto& o = old_.front();
    for (std::size_t i = 0; i < o.size(); ++i) h_.hour[i] -= o[i];
    old_.pop_front();
  }
}

nlohmann::json PanelView::gainSummary() const {
  // Última hora con resultados.
  UtcNs last = 0;
  for (const auto& g : gains) last = std::max(last, g.hourStartUtcNs);
  nlohmann::json per = nlohmann::json::array();
  std::map<std::string, int> count;
  for (const auto& g : gains)
    if (g.hourStartUtcNs == last && g.status != "dark_rate") {
      per.push_back({{"ch", g.ch}, {"G_mV", g.value}, {"status", g.status}, {"peaks", g.peaks}});
      ++count[g.status];
    }
  return {{"hour_start_utc_ns", last}, {"hod2018", per}, {"status_count", count}, {"provisional", true}};
}

PanelConfig PanelConfig::fromJson(const nlohmann::json& j) {
  PanelConfig c;
  if (!j.is_object()) return c;
  c.inTimeAmpMaxMv = j.value("intime_amp_max_mV", c.inTimeAmpMaxMv);
  c.inTimeAmpBinMv = j.value("intime_amp_bin_mV", c.inTimeAmpBinMv);
  c.integralMaxMvNs = j.value("integral_max_mVns", c.integralMaxMvNs);
  c.integralBinMvNs = j.value("integral_bin_mVns", c.integralBinMvNs);
  if (j.contains("hit_threshold_mV")) c.hitThresholdMv = {j["hit_threshold_mV"].at(0).get<double>(), j["hit_threshold_mV"].at(1).get<double>()};
  if (j.contains("mult_thresholds_mV"))
    for (std::size_t s = 0; s < 2; ++s)
      for (std::size_t k = 0; k < 2; ++k) c.multThresholdMv[s][k] = j["mult_thresholds_mV"].at(s).at(k).get<double>();
  c.persistenceEvents = j.value("persistence_events", c.persistenceEvents);
  c.gainPeriodS = j.value("gain_period_s", c.gainPeriodS);
  c.gainSeedMv = j.value("gain_seed_mV", c.gainSeedMv);
  return c;
}

nlohmann::json PanelConfig::toJson() const {
  return {{"intime_amp_max_mV", inTimeAmpMaxMv},
          {"intime_amp_bin_mV", inTimeAmpBinMv},
          {"integral_max_mVns", integralMaxMvNs},
          {"integral_bin_mVns", integralBinMvNs},
          {"hit_threshold_mV", hitThresholdMv},
          {"mult_thresholds_mV", multThresholdMv},
          {"persistence_events", persistenceEvents},
          {"gain_period_s", gainPeriodS},
          {"gain_seed_mV", gainSeedMv}};
}

Panels::Panels(const DaqrunConfig& cfg, const hvmon::Mapping* mapping, int recordLength, double lsbMv)
    : cfg_(cfg), N_(recordLength), lsb_(lsbMv) {
  station_.fill(-1);
  plane_.fill(-1);
  strip_.fill(-1);
  stripSheet_.fill(-1);
  if (mapping)
    for (const auto& e : mapping->channels) {
      if (e.digCh < 0 || e.digCh >= 64) continue;
      const int st = e.station == "hod2019" ? 1 : 0;
      station_[static_cast<std::size_t>(e.digCh)] = st;
      plane_[static_cast<std::size_t>(e.digCh)] = 2 * st + (e.plane == "Y" ? 1 : 0);
      strip_[static_cast<std::size_t>(e.digCh)] = e.strip;
      stripSheet_[static_cast<std::size_t>(e.digCh)] = e.stripSheet;
    }
  const auto& p = cfg.panels;
  const int na = std::max(1, static_cast<int>(p.inTimeAmpMaxMv / p.inTimeAmpBinMv));
  amp_.init("intime_amp", "Amplitud en tiempo (línea base − mínimo en signal);canal;amplitud [mV]", 64, -0.5, 63.5, na, 0,
            na * p.inTimeAmpBinMv);
  const int ni = std::max(1, static_cast<int>(p.integralMaxMvNs / p.integralBinMvNs));
  integ_.init("intime_integral", "Integral en tiempo (QDC);canal;integral [mV#upoint ns]", 64, -0.5, 63.5, ni, 0,
              ni * p.integralBinMvNs);
  minS_.init("min_sample", "Muestra del mínimo (amplitud > umbral de impacto);canal;muestra", 64, -0.5, 63.5, N_, -0.5, N_ - 0.5);
  // Títulos con una segunda línea: qué muestra y cómo se ve normal (DAQ2-11).
  const char* orderTag[2] = {"", "_hoja"};
  const char* orderName[2] = {"orden MEDIDO (mapeo r3)", "orden de la HOJA de octubre"};
  for (int o = 0; o < 2; ++o)
    for (int s = 0; s < 2; ++s)
      hit_[static_cast<std::size_t>(o)][static_cast<std::size_t>(s)].init(
          std::string("hitmap_") + (s ? "hod2019" : "hod2018") + orderTag[o],
          std::string("#splitline{Impactos X #times Y ") + (s ? "hod2019" : "hod2018") + ", " + orderName[o] + " (en tiempo > " +
              std::to_string(static_cast<int>(p.hitThresholdMv[static_cast<std::size_t>(s)])) +
              " mV)}{#scale[0.7]{Pares X-Y del mismo evento; normal: una banda continua cerca de la diagonal o una mancha compacta, no puntos dispersos}};tira X;tira Y",
          16, -0.5, 15.5, 16, -0.5, 15.5);
  for (int q = 0; q < 4; ++q) {
    for (int o = 0; o < 2; ++o)
      coinc_[static_cast<std::size_t>(o)][static_cast<std::size_t>(q)].init(
          std::string("coinc_") + kPlaneName[q] + orderTag[o],
          std::string("#splitline{Coincidencias entre tiras ") + kPlaneName[q] + ", " + orderName[o] +
              "}{#scale[0.7]{Tiras encendidas juntas; normal (orden correcto): cuentas pegadas a la diagonal, en las vecinas}};tira;tira",
          16, -0.5, 15.5, 16, -0.5, 15.5);
    for (int o = 0; o < 2; ++o)
      coincMix_[static_cast<std::size_t>(o)][static_cast<std::size_t>(q)].init(
          std::string("coinc_mix_") + kPlaneName[q] + orderTag[o],
          std::string("#splitline{Coincidencias al azar entre tiras ") + kPlaneName[q] + ", " + orderName[o] +
              "}{#scale[0.7]{Tiras del evento actual con las del anterior (eventos mezclados): lo que da el azar; sin estructura en la diagonal}};tira;tira",
          16, -0.5, 15.5, 16, -0.5, 15.5);
    for (int k = 0; k < 2; ++k) {
      const double thr = p.multThresholdMv[static_cast<std::size_t>(q / 2)][static_cast<std::size_t>(k)];
      mult_[static_cast<std::size_t>(q)][static_cast<std::size_t>(k)].init(
          std::string("mult_") + kPlaneName[q] + "_" + std::to_string(static_cast<int>(thr)) + "mV",
          std::string("#splitline{Multiplicidad ") + kPlaneName[q] + " (> " + std::to_string(static_cast<int>(thr)) +
              " mV)}{#scale[0.7]{Tiras encendidas por evento; normal: casi todo en 0, 1 y 2}};tiras;eventos", 17,
          -0.5, 16.5);
    }
  }
  for (std::size_t g = 0; g < 8; ++g) {
    persPrev_[g].assign(static_cast<std::size_t>(N_ * v_.persistAdcBins), 0);
    persCur_[g] = persPrev_[g];
    v_.persistence[g] = persPrev_[g];
  }
  v_.recordLength = N_;
  v_.lastEvent.assign(static_cast<std::size_t>(64 * N_), 0);
  seed_.fill(p.gainSeedMv);
  publishHists();
}

void Panels::event(const dig::Event& ev, const std::array<double, 64>& baselineAdc, UtcNs tUtc) {
  const auto& w = cfg_.window;
  const int s0 = std::max(0, w.signal.first), s1 = std::min(N_, w.signal.second);
  std::array<std::array<std::uint32_t, 4>, 2> hitMask{};  // [orden][plano]
  std::array<std::array<int, 2>, 4> mult{};
  for (int c = 0; c < 64; ++c) {
    const std::uint16_t* x = ev.adc.data() + static_cast<std::size_t>(c * N_);
    const double bl = baselineAdc[static_cast<std::size_t>(c)];
    std::uint16_t mn = x[s0];
    double integ = 0;
    for (int k = s0; k < s1; ++k) {
      mn = std::min(mn, x[k]);
      integ += (bl - x[k]) * lsb_ * kSampleNs;
    }
    const double amp = (bl - mn) * lsb_;
    amp_.fill(c, amp);
    integ_.fill(c, integ);
    const int st = station_[static_cast<std::size_t>(c)] < 0 ? (c < 32 ? 0 : 1) : station_[static_cast<std::size_t>(c)];
    const double hitThr = cfg_.panels.hitThresholdMv[static_cast<std::size_t>(st)];
    int kmin = 0;
    for (int k = 1; k < N_; ++k)
      if (x[k] < x[kmin]) kmin = k;
    if ((bl - x[kmin]) * lsb_ > hitThr) minS_.fillBin(c, kmin);
    const int q = plane_[static_cast<std::size_t>(c)], strip = strip_[static_cast<std::size_t>(c)];
    const int stripH = stripSheet_[static_cast<std::size_t>(c)];
    if (q >= 0 && strip >= 0 && strip < 16) {
      if (amp > hitThr) {
        hitMask[0][static_cast<std::size_t>(q)] |= 1u << strip;
        if (stripH >= 0 && stripH < 16) hitMask[1][static_cast<std::size_t>(q)] |= 1u << stripH;
      }
      for (std::size_t k = 0; k < 2; ++k)
        if (amp > cfg_.panels.multThresholdMv[static_cast<std::size_t>(q / 2)][k]) ++mult[static_cast<std::size_t>(q)][k];
    }
    const auto g = static_cast<std::size_t>(c / 8);
    for (int k = 0; k < N_; ++k) ++persCur_[g][static_cast<std::size_t>(k * v_.persistAdcBins + std::min(v_.persistAdcBins - 1, x[k] >> 4))];
  }
  for (std::size_t o = 0; o < 2; ++o) {
    for (int s = 0; s < 2; ++s) {
      const std::uint32_t mx = hitMask[o][static_cast<std::size_t>(2 * s)], my = hitMask[o][static_cast<std::size_t>(2 * s + 1)];
      for (int i = 0; i < 16; ++i)
        if (mx & (1u << i))
          for (int j = 0; j < 16; ++j)
            if (my & (1u << j)) hit_[o][static_cast<std::size_t>(s)].fillBin(i, j);
    }
    for (std::size_t q = 0; q < 4; ++q) {
      const std::uint32_t m = hitMask[o][q];
      for (int i = 0; i < 16; ++i)
        if (m & (1u << i))
          for (int j = 0; j < 16; ++j)
            if (m & (1u << j)) coinc_[o][q].fillBin(i, j);
      // Eventos mezclados: (i en este evento, j en el anterior), i ≠ j; en esperanza simétrica como coinc_.
      const std::uint32_t pm = prevMask_[o][q];
      if (v_.events > 0)
        for (int i = 0; i < 16; ++i)
          if (m & (1u << i))
            for (int j = 0; j < 16; ++j)
              if (j != i && (pm & (1u << j))) coincMix_[o][q].fillBin(i, j);
      prevMask_[o][q] = m;
    }
  }
  for (std::size_t q = 0; q < 4; ++q)
    for (std::size_t k = 0; k < 2; ++k) mult_[q][k].fillBin(mult[q][k]);
  if (++persCount_ >= std::max(1, cfg_.panels.persistenceEvents / 2)) {  // bloque anterior + en curso ≈ últimos ~1000
    persPrev_ = persCur_;
    for (auto& v : persCur_) std::fill(v.begin(), v.end(), 0u);
    persCount_ = 0;
  }
  std::copy(ev.adc.begin(), ev.adc.end(), v_.lastEvent.begin());
  v_.lastEventUtcNs = tUtc;
  ++v_.events;
}

void Panels::publishHists() {
  v_.hists.clear();
  for (const RollingHist* h : {&amp_, &integ_, &minS_}) v_.hists.push_back(h->view());
  for (const auto& o : hit_)
    for (const auto& h : o) v_.hists.push_back(h.view());
  for (std::size_t o = 0; o < 2; ++o)
    for (std::size_t q = 0; q < 4; ++q) {
      const auto& h = coinc_[o][q].view();
      v_.hists.push_back(h);
      double all = 0, nb = 0;
      for (int i = 0; i < 16; ++i)
        for (int j = 0; j < 16; ++j)
          if (i != j) {
            const double x = static_cast<double>(h.cum[static_cast<std::size_t>(j * 16 + i)]);
            all += x;
            if (std::abs(i - j) == 1) nb += x;
          }
      v_.neighbourFrac[o][q] = all > 0 ? nb / all : -1;
      v_.neighbourPairs[o][q] = all / 2;
    }
  for (const auto& o : coincMix_)
    for (const auto& h : o) v_.hists.push_back(h.view());
  for (const auto& q : mult_)
    for (const auto& h : q) v_.hists.push_back(h.view());
  for (std::size_t g = 0; g < 8; ++g) {
    auto& out = v_.persistence[g];
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = persPrev_[g][i] + persCur_[g][i];
  }
}

void Panels::period(const MonitorSummary& s, const std::vector<std::uint32_t>& ampCum, int ampBins, double ampBinMv) {
  for (RollingHist* h : {&amp_, &integ_, &minS_}) h->closeSlice();
  for (auto& o : hit_)
    for (auto& h : o) h.closeSlice();
  for (auto& o : coinc_)
    for (auto& h : o) h.closeSlice();
  for (auto& o : coincMix_)
    for (auto& h : o) h.closeSlice();
  for (auto& q : mult_)
    for (auto& h : q) h.closeSlice();
  // Mapa canal × tiempo y vivos por plano.
  std::array<float, 64> r{};
  std::array<int, 4> live{};
  for (std::size_t c = 0; c < 64; ++c) {
    r[c] = static_cast<float>(s.ootTimeS[c] > 0 ? s.ootCounts[c] / s.ootTimeS[c] : 0.0);
    const int q = plane_[c] >= 0 ? plane_[c] : static_cast<int>(c / 16);
    if (s.lowCounts[c] > 0) ++live[static_cast<std::size_t>(q)];
  }
  v_.rateMapT.push_back(s.tEndUtcNs);
  v_.rateMap.push_back(r);
  v_.livePerPlane.push_back(live);
  if (v_.rateMapT.size() > 1440) {
    v_.rateMapT.erase(v_.rateMapT.begin());
    v_.rateMap.erase(v_.rateMap.begin());
    v_.livePerPlane.erase(v_.livePerPlane.begin());
  }
  // hod2019: tasa > umbral provisional (3 mV) acumulada en la hora.
  for (std::size_t c = 32; c < 64; ++c) {
    darkCounts_[c] += s.ootCounts[c];
    darkTime_[c] += s.ootTimeS[c];
  }
  // Ganancia por hora (hod2018) desde el espectro fuera de tiempo acumulado del monitor.
  if (hourStart_ == 0) {
    hourStart_ = s.tStartUtcNs;
    ampAtHour_ = ampCum;
  } else if (static_cast<double>(s.tEndUtcNs - hourStart_) / 1e9 >= cfg_.panels.gainPeriodS && !ampCum.empty()) {
    for (int c = 0; c < 32; ++c) {
      std::vector<std::uint32_t> spec(static_cast<std::size_t>(ampBins));
      for (int b = 0; b < ampBins; ++b) {
        const auto i = static_cast<std::size_t>(c * ampBins + b);
        spec[static_cast<std::size_t>(b)] = ampCum[i] - (i < ampAtHour_.size() ? ampAtHour_[i] : 0u);
      }
      auto& seed = seed_[static_cast<std::size_t>(c)];
      const auto res = gain::fitGain(spec, ampBinMv, seed);
      v_.gains.push_back(GainPoint{hourStart_, c, res.gainMv, res.status, res.peaks});
      if (res.accepted()) seed = res.gainMv;  // semilla de la hora siguiente
    }
    for (std::size_t c = 32; c < 64; ++c) {
      v_.gains.push_back(GainPoint{hourStart_, static_cast<int>(c), darkTime_[c] > 0 ? darkCounts_[c] / darkTime_[c] : 0.0, "dark_rate", 0});
      darkCounts_[c] = darkTime_[c] = 0;
    }
    hourStart_ = s.tEndUtcNs;
    ampAtHour_ = ampCum;
  }
  v_.tUtcNs = s.tEndUtcNs;
  publishHists();
}

}  // namespace daq::daqrun
