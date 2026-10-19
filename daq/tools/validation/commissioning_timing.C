// commissioning_timing.C — DAQ2-6 (solo lectura). TTT, contador y Δt de un run de daqrun.
//  - TTT crudo: máximo observado y caídas (vueltas); se desenrolla con `bits` y la hora del PC.
//  - Paso del TTT: ajuste lineal de TTT desenrollado (ticks) contra t_read_utc_ns (ns del PC):
//    pendiente = ticks/ns -> tick = 1/pendiente; residuos RMS (incluyen la latencia de lectura).
//  - Período de vuelta = 2^bits × tick.
//  - Contador de la placa: máximo, saltos (≠ +1 módulo 2^24).
//  - Δt (con el TTT y el tick medido): trigger_interval_stats.h (tasa, KS, χ², huecos, ráfagas).
// Uso: root -l -b -q 'daq/tools/validation/commissioning_timing.C("<dir>/<base>_seg*.root", 31)'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "daq/trigger_interval_stats.h"

void commissioning_timing(const char* pattern, int bits = 31) {
  TChain t("events");
  t.Add(pattern);
  UInt_t raw = 0, cnt = 0;
  Long64_t tread = 0;
  t.SetBranchStatus("*", 0);
  for (const char* b : {"ttt_raw", "evt_counter_hw", "t_read_utc_ns"}) t.SetBranchStatus(b, 1);
  t.SetBranchAddress("ttt_raw", &raw);
  t.SetBranchAddress("evt_counter_hw", &cnt);
  t.SetBranchAddress("t_read_utc_ns", &tread);
  const Long64_t n = t.GetEntries();
  std::vector<double> x, y;  // x: t_read (ns desde el primero), y: ticks desenrollados
  UInt_t maxRaw = 0, maxCnt = 0, prevRaw = 0, prevCnt = 0;
  long drops = 0, jumps = 0, above31 = 0;
  long double wraps = 0;
  const long double period = std::ldexp(1.0L, bits);
  Long64_t t0 = 0;
  for (Long64_t i = 0; i < n; ++i) {
    t.GetEntry(i);
    if (i == 0) t0 = tread;
    maxRaw = std::max(maxRaw, raw);
    maxCnt = std::max(maxCnt, cnt);
    if (raw >> 31) ++above31;
    const UInt_t r = bits >= 32 ? raw : (raw & ((1u << bits) - 1));
    if (i > 0 && r < prevRaw) {
      ++drops;
      wraps += 1;
    }
    if (i > 0 && ((cnt - prevCnt) & 0xFFFFFF) != 1) ++jumps;
    prevRaw = r;
    prevCnt = cnt;
    x.push_back(static_cast<double>(tread - t0));
    y.push_back(static_cast<double>(wraps * period + r));
  }
  // Regresión lineal y = a + b x.
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    sx += x[i];
    sy += y[i];
    sxx += x[i] * x[i];
    sxy += x[i] * y[i];
  }
  const double m = static_cast<double>(x.size());
  const double b = (m * sxy - sx * sy) / (m * sxx - sx * sx), a = (sy - b * sx) / m;
  double res2 = 0, resMax = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    const double r = (y[i] - a - b * x[i]) / b;  // ns
    res2 += r * r;
    resMax = std::max(resMax, std::fabs(r));
  }
  const double tick = 1.0 / b;
  printf("eventos %lld; TTT crudo máximo %u (0x%08X), eventos con bit 31: %ld; caídas con %d bits: %ld\n", n, maxRaw, maxRaw, above31,
         bits, drops);
  printf("ajuste TTT vs PC: tick = %.6f ns, residuo RMS %.1f µs, máx %.1f µs; período de vuelta 2^%d × tick = %.3f s\n", tick,
         std::sqrt(res2 / m) / 1e3, resMax / 1e3, bits, static_cast<double>(period) * tick / 1e9);
  printf("contador: máximo %u, saltos ≠ +1 (mód 2^24): %ld\n", maxCnt, jumps);
  std::vector<std::int64_t> tt;
  for (double v : y) tt.push_back(static_cast<std::int64_t>(std::llround(v * tick)));
  const auto s = daq::tis::summarize(tt);
  printf("Δt: tasa %.4f ± %.4f Hz, Δt medio %.4f s, máx %.3f s, KS D %.4f p %.3g, χ² %.1f/%d p %.3g, huecos %zu (T = %.1f s), "
         "ráfagas %zu (%zu eventos)\n",
         s.rateHz, s.rateErrHz, s.meanDtS, s.maxDtS, s.ksD, s.ksP, s.chi2, s.chi2Ndf, s.chi2P, s.nGaps, s.gapThresholdS, s.nBursts,
         s.burstEvents);
}
