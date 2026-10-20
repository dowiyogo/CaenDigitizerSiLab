// Tests de los paneles en tiempo (DAQ2-9): ajuste de ganancia por separación de picos (G, G/2 y 2G
// ambiguos, espaciado irregular, insuficiente), última hora de RollingHist y llenado de mapas de
// impactos, coincidencias y muestra del mínimo con un muón sintético.
#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "daq/daqrun/config.h"
#include "daq/daqrun/monitor.h"
#include "daq/daqrun/panels.h"
#include "daq/gain_fit.h"

using namespace daq;
using namespace daq::daqrun;

namespace {
constexpr double kLsb = 2000.0 / 4095.0;

// Espectro de picos gaussianos en p0 + n·G (n = 1..6), anchos crecientes, más un fondo bajo.
std::vector<std::uint32_t> spectrum(double p0, double g, double sigma0 = 0.6, unsigned seed = 3) {
  std::vector<double> mu(205, 0.0);
  for (int n = 1; n <= 6; ++n) {
    const double pos = p0 + n * g, s = std::sqrt(sigma0 * sigma0 + 0.04 * n), a = 4000.0 * std::exp(-0.45 * n);
    for (std::size_t b = 0; b < mu.size(); ++b) {
      const double x = (static_cast<double>(b) + 0.5) * kLsb;
      mu[b] += a * kLsb / (s * std::sqrt(2 * M_PI)) * std::exp(-0.5 * (x - pos) * (x - pos) / (s * s));
    }
  }
  std::mt19937 rng(seed);
  std::vector<std::uint32_t> out(mu.size());
  for (std::size_t b = 0; b < mu.size(); ++b) out[b] = static_cast<std::uint32_t>(std::poisson_distribution<int>(mu[b] + 2.0)(rng));
  return out;
}
}  // namespace

TEST(GainFit, RecoversGainAndFlagsHalfAndDoubleSeedsAsAmbiguous) {
  const auto sp = spectrum(-0.5, 4.0);
  const auto ok = gain::fitGain(sp, kLsb, 4.0);
  EXPECT_EQ(ok.status, "ok");
  EXPECT_NEAR(ok.gainMv, 4.0, 0.2);
  EXPECT_GE(ok.peaks, 3);
  EXPECT_TRUE(ok.accepted());
  EXPECT_EQ(gain::fitGain(sp, kLsb, 2.0).status, "ambiguous");  // semilla G/2: la medida es 2× la semilla
  EXPECT_EQ(gain::fitGain(sp, kLsb, 8.0).status, "ambiguous");  // semilla 2G
  const auto ch = gain::fitGain(sp, kLsb, 5.2);
  EXPECT_EQ(ch.status, "changed");  // cambio > 20 %: se acepta y se marca
  EXPECT_TRUE(ch.accepted());
}

TEST(GainFit, IrregularSpacingIsRejectedAndFewPeaksInsufficient) {
  // Picos en 6, 10, 17, 21 mV: diferencias 4, 7, 4 (no compatibles con una sola G).
  std::vector<std::uint32_t> sp(205, 2);
  for (double pos : {6.0, 10.0, 17.0, 21.0})
    for (std::size_t b = 0; b < sp.size(); ++b) {
      const double x = (static_cast<double>(b) + 0.5) * kLsb;
      sp[b] += static_cast<std::uint32_t>(600 * std::exp(-0.5 * (x - pos) * (x - pos) / 0.36));
    }
  EXPECT_EQ(gain::fitGain(sp, kLsb, 4.0).status, "rejected");
  std::vector<std::uint32_t> flat(205, 3);
  EXPECT_EQ(gain::fitGain(flat, kLsb, 4.0).status, "insufficient");
}

TEST(Panels, RollingHistLastHourDropsOldSlices) {
  RollingHist h;
  h.init("t", "t", 4, 0, 4, 1, 0, 1, 3);  // 3 rebanadas
  h.fill(0.5);
  h.closeSlice();
  h.fill(1.5);
  h.closeSlice();
  h.fill(1.5);
  EXPECT_EQ(h.view().hour[0], 1u);
  h.closeSlice();  // la primera rebanada sale de la "hora"
  EXPECT_EQ(h.view().hour[0], 0u);
  EXPECT_EQ(h.view().hour[1], 2u);
  EXPECT_EQ(h.view().cum[0], 1u);  // el acumulado no pierde nada
}

TEST(Panels, MuonFillsHitMapCoincidenceAndMinimumSample) {
  auto cfg = DaqrunConfig::fromJson({{"campaign", {{"name", "t"}}}}).value();
  auto m = hvmon::Mapping::fromFile(DAQ_SOURCE_DIR "/config/mapping/hodo_mapping_2026_10.json").value();
  const int N = 144;
  Panels p(cfg, &m, N, kLsb);
  dig::Event ev;
  ev.adc.assign(64 * N, 2100);
  std::array<double, 64> bl{};
  bl.fill(2100);
  // hod2018: X tiras 3 y 4 (ch 3, 4), Y tira 10 (ch 26); mínimo en la muestra 53, 120 mV.
  for (int c : {3, 4, 26}) {
    ev.adc[static_cast<std::size_t>(c * N + 53)] = static_cast<std::uint16_t>(2100 - 120 / kLsb);
    ev.adc[static_cast<std::size_t>(c * N + 54)] = static_cast<std::uint16_t>(2100 - 80 / kLsb);
  }
  p.event(ev, bl, 1);
  p.refresh();
  const PanelView& v = p.view();
  auto find = [&](const std::string& n) -> const PanelHist& {
    for (const auto& h : v.hists)
      if (h.name == n) return h;
    ADD_FAILURE() << n;
    return v.hists.front();
  };
  const auto& hm = find("hitmap_hod2018");
  EXPECT_EQ(hm.cum[static_cast<std::size_t>(10 * 16 + 3)], 1u);
  EXPECT_EQ(hm.cum[static_cast<std::size_t>(10 * 16 + 4)], 1u);
  std::uint64_t total = 0;
  for (auto x : hm.cum) total += x;
  EXPECT_EQ(total, 2u);
  const auto& cx = find("coinc_hod2018_X");
  EXPECT_EQ(cx.cum[static_cast<std::size_t>(3 * 16 + 4)], 1u);  // vecinas 3-4
  EXPECT_EQ(cx.cum[static_cast<std::size_t>(4 * 16 + 4)], 1u);  // diagonal: tira sola
  const auto& ms = find("min_sample");
  EXPECT_EQ(ms.cum[static_cast<std::size_t>(53 * 64 + 3)], 1u);
  const auto& mult = find("mult_hod2018_X_20mV");
  EXPECT_EQ(mult.cum[2], 1u);
  EXPECT_EQ(v.lastEvent[static_cast<std::size_t>(3 * N + 53)], ev.adc[static_cast<std::size_t>(3 * N + 53)]);
}
