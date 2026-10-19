// Tests de las piezas compartidas con el análisis: SHA-256, búsqueda de pulsos y estadística de Δt.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <vector>

#include "daq/core/sha256.h"
#include "daq/pulse_finder.h"
#include "daq/trigger_interval_stats.h"

using namespace daq;

TEST(Sha256, KnownVectorsAndStreaming) {
  EXPECT_EQ(sha256Hex(""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  EXPECT_EQ(sha256Hex("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  EXPECT_EQ(sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
  const std::string million(1000000, 'a');
  EXPECT_EQ(sha256Hex(million), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
  // Por partes de tamaños irregulares = de una vez; y el archivo leído por bloques.
  Sha256 h;
  for (std::size_t i = 0; i < million.size();) {
    const std::size_t n = std::min<std::size_t>(1 + (i % 97), million.size() - i);
    h.update(million.data() + i, n);
    i += n;
  }
  EXPECT_EQ(h.hexDigest(), sha256Hex(million));
  const std::string f = ::testing::TempDir() + "daq_sha_test.bin";
  { std::ofstream(f, std::ios::binary) << million; }
  auto r = sha256File(f);
  ASSERT_TRUE(r.ok());
  EXPECT_EQ(r.value(), sha256Hex(million));
  std::remove(f.c_str());
  EXPECT_FALSE(sha256File("/no/existe").ok());
}

namespace {

// Forma de onda sintética en cuentas ADC: línea base + ruido gaussiano + pulsos negativos
// (decaimiento exponencial de 2 muestras) en posiciones conocidas.
std::vector<std::uint16_t> waveform(std::mt19937_64& rng, int n, double baseline, double noiseMv,
                                    const std::vector<std::pair<int, double>>& pulses) {
  const double lsb = pulse::lsbMv();
  std::vector<double> w(static_cast<std::size_t>(n), baseline);
  for (const auto& p : pulses)
    for (int k = p.first; k < n && k < p.first + 16; ++k) w[static_cast<std::size_t>(k)] -= p.second / lsb * std::exp(-(k - p.first) / 2.0);
  std::normal_distribution<double> g(0.0, noiseMv / lsb);
  std::vector<std::uint16_t> out(static_cast<std::size_t>(n));
  for (int k = 0; k < n; ++k) out[static_cast<std::size_t>(k)] = static_cast<std::uint16_t>(std::lround(w[static_cast<std::size_t>(k)] + g(rng)));
  return out;
}

struct Score {
  double efficiency, falsePerSample;
};

// Pulsos en [0,140) con separación > 8 muestras; cuenta los encontrados a ±1 muestra.
Score score(double noiseMv, double ampMv, double thrMv, int trials) {
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> pos(2, 130);
  int injected = 0, found = 0, fakes = 0;
  const double lsb = pulse::lsbMv();
  const std::vector<pulse::Region> all{{0, 144}};
  for (int t = 0; t < trials; ++t) {
    std::vector<std::pair<int, double>> ps;
    if (t % 2 == 0) {
      const int a = pos(rng);
      ps.push_back({a, ampMv});
      ++injected;
    }
    auto w = waveform(rng, 144, 2090, noiseMv, ps);
    const double bl = pulse::trimmedMean(w.data(), pulse::Region{0, 144}, 2.0 / lsb);
    std::vector<pulse::Pulse> out;
    pulse::findPulses(w.data(), 144, bl, lsb, thrMv, 4, all, out);
    for (const auto& p : out) {
      bool match = false;
      for (const auto& q : ps) match = match || std::abs(p.sample - q.first) <= 1;
      if (match) ++found;
      else ++fakes;
    }
  }
  return {static_cast<double>(found) / injected, static_cast<double>(fakes) / (trials * 144.0)};
}

}  // namespace

TEST(PulseFinder, EfficiencyAndFalseRateOnSyntheticWaveforms) {
  // hod2018: ruido 1.9 mV, pulso de ~5 p.e. (24 mV) con el umbral provisional de 20 mV; y pulso de
  // ~3 p.e. (14 mV) con umbral de 10 mV (5.3 σ del ruido): eficiencia alta, falsos raros.
  const Score a = score(1.9, 24.0, 20.0, 4000);
  EXPECT_GT(a.efficiency, 0.97);
  EXPECT_LT(a.falsePerSample, 1e-5);
  const Score b = score(1.9, 14.0, 10.0, 4000);
  EXPECT_GT(b.efficiency, 0.97);
  EXPECT_LT(b.falsePerSample, 1e-4);
  // hod2019: ruido 0.25 mV, pulso de 1 p.e. (~4 mV), umbral 3 mV (12 σ): sin falsos.
  const Score c = score(0.25, 4.0, 3.0, 4000);
  EXPECT_GT(c.efficiency, 0.97);
  EXPECT_EQ(c.falsePerSample, 0.0);
  std::printf("eficiencia/falsos por muestra: hod2018 24mV@20mV %.4f/%.2e, 14mV@10mV %.4f/%.2e, hod2019 4mV@3mV %.4f/%.2e\n",
              a.efficiency, a.falsePerSample, b.efficiency, b.falsePerSample, c.efficiency, c.falsePerSample);
}

TEST(PulseFinder, MinimumSeparationKeepsLargestAndRegionsAndUnits) {
  std::mt19937_64 rng(1);
  // Dos pulsos a 2 muestras: queda el mayor. Otro a 30 muestras: aparte.
  auto w = waveform(rng, 144, 2090, 0.0, {{50, 10.0}, {52, 30.0}, {82, 15.0}});
  std::vector<pulse::Pulse> out;
  pulse::findPulses(w.data(), 144, 2090.0, pulse::lsbMv(), 5.0, 4, {{0, 144}}, out);
  ASSERT_EQ(out.size(), 2u);
  EXPECT_EQ(out[0].sample, 52);
  EXPECT_NEAR(out[0].amplitudeMv, 30.0 + 10.0 * std::exp(-1.0), 0.5);
  EXPECT_EQ(out[1].sample, 82);
  // Solo la región [60, 144).
  out.clear();
  pulse::findPulses(w.data(), 144, 2090.0, pulse::lsbMv(), 5.0, 4, {{60, 144}}, out);
  ASSERT_EQ(out.size(), 1u);
  EXPECT_EQ(out[0].sample, 82);
  // Muestras en mV (archivos antiguos): escala 1.
  std::vector<float> mv(144, 20.0f);
  mv[70] = 5.0f;
  EXPECT_EQ(pulse::countPulses(mv.data(), 144, pulse::trimmedMean(mv.data(), {0, 40}, 2.0), 1.0, 10.0, 4, {{0, 144}}), 1);
  EXPECT_DOUBLE_EQ(pulse::thresholdFromPe(3.0, 0.5, 4.5), 14.0);
  EXPECT_NEAR(pulse::lsbMv(), 0.48840048840, 1e-10);
}

TEST(TriggerIntervalStats, ExponentialFitGapsAndBursts) {
  std::mt19937_64 rng(7);
  std::exponential_distribution<double> ex(0.6);
  std::vector<std::int64_t> t;
  double now = 0;
  for (int i = 0; i < 20000; ++i) {
    now += ex(rng);
    t.push_back(static_cast<std::int64_t>(now * 1e9));
  }
  auto s = tis::summarize(t);
  EXPECT_NEAR(s.rateHz, 0.6, 3 * s.rateErrHz);
  EXPECT_GT(s.ksP, 1e-3);
  EXPECT_GT(s.chi2P, 1e-3);
  EXPECT_EQ(s.nGaps, 0u);
  EXPECT_EQ(s.nBursts, 0u);
  // Hueco de 120 s (sin triggers) a mitad del run y ráfaga de 30 triggers en 1 s.
  std::vector<std::int64_t> g;
  const std::int64_t mid = t[10000];
  for (auto x : t) {
    if (x > mid && x < mid + 120000000000LL) continue;
    g.push_back(x);
  }
  const std::int64_t b0 = t[15000] + 1000000;
  for (int i = 0; i < 30; ++i) g.push_back(b0 + i * 30000000LL);
  std::sort(g.begin(), g.end());
  auto s2 = tis::summarize(g);
  EXPECT_EQ(s2.nGaps, 1u);
  ASSERT_EQ(s2.gapIndices.size(), 1u);
  EXPECT_GT(static_cast<double>(g[s2.gapIndices[0] + 1] - g[s2.gapIndices[0]]) * 1e-9, 100.0);
  EXPECT_NEAR(s2.gapThresholdS, -std::log(1e-6) / s2.rateHz, 1e-9);
  EXPECT_GE(s2.nBursts, 1u);
  EXPECT_GE(s2.burstEvents, 30u);
  std::printf("tasa %.4f ± %.4f Hz, KS p %.3f, chi2 p %.3f; con hueco y ráfaga: huecos %zu (T = %.1f s), ráfagas %zu\n",
              s.rateHz, s.rateErrHz, s.ksP, s.chi2P, s2.nGaps, s2.gapThresholdS, s2.nBursts);
  // Resolución declarada: eventos leídos juntos comparten marca.
  std::vector<std::int64_t> q;
  for (auto x : t) q.push_back(x / 100000000 * 100000000);  // 100 ms
  tis::Options o;
  o.resolutionNs = 1e8;
  auto s3 = tis::summarize(q, o);
  EXPECT_GT(s3.nUnresolved, 0u);
  EXPECT_NEAR(s3.rateHz, 0.6, 0.03);
  // Colas de Poisson y p-valor de χ².
  EXPECT_NEAR(tis::poissonUpperTail(1, 2.0), 1 - std::exp(-2.0), 1e-12);
  EXPECT_NEAR(tis::poissonUpperTail(3, 1.0), 1 - std::exp(-1.0) * (1 + 1 + 0.5), 1e-12);
  EXPECT_NEAR(tis::chi2Pvalue(2.0, 2), std::exp(-1.0), 1e-12);
  EXPECT_LT(tis::poissonUpperTail(30, 5.0), 1e-9);
}
