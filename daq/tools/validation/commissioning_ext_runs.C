// commissioning_ext_runs.C — DAQ2-6, paso 4.1 (solo lectura). Varios runs de daqrun con trigger
// externo, ventana B y HV encendida, analizados juntos (runs 7 y 8 de commissioning_2026_10).
//  1. Δt: t_utc_ns de todos los runs; el intervalo entre el último evento de un run y el primero del
//     siguiente (tiempo muerto entre runs) se quita desplazando los runs siguientes, así que hay
//     N − 1 − (runs − 1) intervalos reales. trigger_interval_stats.h (tasa, KS, χ², huecos, ráfagas)
//     y ajuste exponencial (log-verosimilitud, TH1 de Δt en [0, 15) s); tasa por bin de bin_s.
//     Contador de la placa: saltos ≠ +1 dentro de cada run; flags de evento Gap/Burst/BoardFull.
//  2. Muón: por estación, histograma de la muestra del mínimo en eventos con amplitud > umbral
//     (hod2018 amp18, hod2019 amp19); se imprimen las cuentas por muestra en [40, 70), el fondo
//     medio por muestra en [10,40)+[102,140) y el exceso en [sig_lo, sig_hi).
//  3. Fuera de tiempo: zonas [0,49)+[102,144) (91 muestras), separación 4; hod2018 > 3·G (G de v8,
//     mayo: solo referencia, no es la calibración de octubre), hod2019 > 3 mV y > 5 mV. Por canal en
//     todo el período y por bin de bin_s (CSV <out>.oot_bins.csv); medianas por estación por bin.
// Salidas: <out>.root (histogramas), <out>.oot_channels.csv, <out>.oot_bins.csv, <out>.rate_bins.csv.
// Uso: root -l -b -q 'daq/tools/validation/commissioning_ext_runs.C("<dir>/<base>T17*_seg*.root", "out")'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"
#include "daq/trigger_interval_stats.h"

namespace {
double median(std::vector<double> v) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
}  // namespace

void commissioning_ext_runs(const char* pattern, const char* out, double binS = 600.0, int sigLo = 49, int sigHi = 58,
                            double amp18 = 50.0, double amp19 = 5.0, const char* gains = "outputs/pe_calibration_v8.txt") {
  std::vector<double> G(64, NAN);
  {
    std::ifstream in(gains);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int ch;
      std::string pl;
      double p0, g;
      if (ss >> ch >> pl >> p0 >> g && ch >= 0 && ch < 64) G[ch] = g;
    }
  }
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  if (n != 144) {
    printf("registro de %d muestras: esta macro es para la ventana B (144)\n", n);
    return;
  }
  std::vector<UShort_t> adc(64 * n);
  ULong64_t seq = 0;
  UInt_t cnt = 0, flags = 0;
  Long64_t tutc = 0;
  t.SetBranchStatus("*", 0);
  for (const char* b : {"adc", "evt_seq", "evt_counter_hw", "evt_flags", "t_utc_ns"}) t.SetBranchStatus(b, 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("evt_seq", &seq);
  t.SetBranchAddress("evt_counter_hw", &cnt);
  t.SetBranchAddress("evt_flags", &flags);
  t.SetBranchAddress("t_utc_ns", &tutc);
  const Long64_t ne = t.GetEntries();
  const double lsb = daq::pulse::lsbMv();
  const std::vector<daq::pulse::Region> zones{{0, 49}, {102, 144}};
  const double exposed = 91 * 16e-9;  // s por evento y canal
  t.GetEntry(0);
  const Long64_t tFirst = tutc;
  t.GetEntry(ne - 1);
  const int nbins = static_cast<int>((tutc - tFirst) / 1e9 / binS) + 1;

  TFile fo((std::string(out) + ".root").c_str(), "RECREATE");
  TH1D hPos18("min_sample_hod2018", Form("hod2018: muestra del mínimo (> %.0f mV);muestra;canales×eventos", amp18), n, -0.5, n - 0.5);
  TH1D hPos19("min_sample_hod2019", Form("hod2019: muestra del mínimo (> %.0f mV);muestra;canales×eventos", amp19), n, -0.5, n - 0.5);
  TH1D hDt("dt", "#Deltat entre triggers;#Deltat [s];eventos", 150, 0, 15);

  std::vector<std::int64_t> tt;  // t_utc_ns con el tiempo muerto entre runs quitado
  std::vector<std::pair<Long64_t, Long64_t>> runs;  // (primer, último) t_utc_ns por run
  std::vector<long> evBin(nbins, 0);
  std::vector<long> ootAll18(64, 0), ootAll19a(64, 0), ootAll19b(64, 0);
  // conteos por bin y canal: hod2018 > 3G (o hod2019 > 3 mV), hod2019 > 5 mV
  std::vector<std::vector<long>> binA(nbins, std::vector<long>(64, 0)), binB(nbins, std::vector<long>(64, 0));
  long jumps = 0, fGap = 0, fBurst = 0, fFull = 0;
  ULong64_t prevSeq = 0;
  UInt_t prevCnt = 0;
  Long64_t shift = 0, prevT = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const bool newRun = e == 0 || seq < prevSeq;
    if (newRun) {
      if (e > 0) shift += tutc - prevT;  // se quita el intervalo entre runs
      runs.push_back({tutc, tutc});
    } else if (((cnt - prevCnt) & 0xFFFFFF) != 1) {
      ++jumps;
    }
    runs.back().second = tutc;
    if (!newRun) hDt.Fill((tutc - prevT) / 1e9);
    tt.push_back(tutc - shift);
    prevSeq = seq;
    prevCnt = cnt;
    prevT = tutc;
    if (flags & 1u) ++fGap;
    if (flags & 2u) ++fBurst;
    if (flags & 8u) ++fFull;
    const int b = std::min(nbins - 1, static_cast<int>((tutc - tFirst) / 1e9 / binS));
    ++evBin[b];
    for (int c = 0; c < 64; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      int kmin = 0;
      for (int k = 1; k < n; ++k)
        if (x[k] < x[kmin]) kmin = k;
      const double amp = (bl - x[kmin]) * lsb;
      if (amp > (c < 32 ? amp18 : amp19)) (c < 32 ? hPos18 : hPos19).Fill(kmin);
      std::vector<daq::pulse::Pulse> p;
      daq::pulse::findPulses(x, n, bl, lsb, 2 * lsb, 4, zones, p);
      for (const auto& q : p) {
        if (c < 32) {
          if (std::isfinite(G[c]) && q.amplitudeMv > 3.0 * G[c]) {
            ++ootAll18[c];
            ++binA[b][c];
          }
        } else {
          if (q.amplitudeMv > 3.0) {
            ++ootAll19a[c];
            ++binA[b][c];
          }
          if (q.amplitudeMv > 5.0) {
            ++ootAll19b[c];
            ++binB[b][c];
          }
        }
      }
    }
  }
  hPos18.Write();
  hPos19.Write();
  hDt.Write();

  // 1. Δt
  printf("eventos %lld en %zu runs; registro %d muestras\n", ne, runs.size(), n);
  double live = 0;
  for (size_t r = 0; r < runs.size(); ++r) {
    printf("  run %zu: primer evento %lld, último %lld, %.1f s\n", r + 1, runs[r].first, runs[r].second,
           (runs[r].second - runs[r].first) / 1e9);
    live += (runs[r].second - runs[r].first) / 1e9;
    if (r > 0) printf("  intervalo entre runs quitado: %.3f s\n", (runs[r].first - runs[r - 1].second) / 1e9);
  }
  const auto s = daq::tis::summarize(tt);
  printf("Δt (sin los intervalos entre runs): tasa %.4f ± %.4f Hz, Δt medio %.4f s, máx %.3f s, KS D %.4f p %.3g, χ² %.1f/%d p %.3g, "
         "huecos %zu (T = %.1f s), ráfagas %zu (%zu eventos)\n",
         s.rateHz, s.rateErrHz, s.meanDtS, s.maxDtS, s.ksD, s.ksP, s.chi2, s.chi2Ndf, s.chi2P, s.nGaps, s.gapThresholdS, s.nBursts,
         s.burstEvents);
  TF1 fe("fe", "[0]*exp(-[1]*x)", 0, 15);
  fe.SetParameters(hDt.GetEntries() * 0.1 * 0.65, 0.65);
  hDt.Fit(&fe, "QLN0");
  printf("ajuste exponencial de Δt (L, [0,15) s, %.0f entradas): pendiente %.4f ± %.4f Hz, χ²/ndf %.1f/%d\n", hDt.GetEntries(),
         fe.GetParameter(1), fe.GetParError(1), fe.GetChisquare(), fe.GetNDF());
  printf("contador de la placa: saltos ≠ +1 dentro de los runs %ld; flags de evento: Gap %ld, Burst %ld, BoardFull %ld\n", jumps, fGap,
         fBurst, fFull);
  FILE* rb = fopen((std::string(out) + ".rate_bins.csv").c_str(), "w");
  fprintf(rb, "bin,t_start_utc_ns,events,live_s,rate_Hz\n");
  printf("tasa por bin de %.0f s (tiempo vivo = parte del bin dentro de algún run):\n", binS);
  for (int b = 0; b < nbins; ++b) {
    const Long64_t b0 = tFirst + static_cast<Long64_t>(b * binS * 1e9), b1 = tFirst + static_cast<Long64_t>((b + 1) * binS * 1e9);
    double lv = 0;
    for (const auto& r : runs) lv += std::max<Long64_t>(0, std::min(b1, r.second) - std::max(b0, r.first)) / 1e9;
    fprintf(rb, "%d,%lld,%ld,%.1f,%.4f\n", b, b0, evBin[b], lv, lv > 0 ? evBin[b] / lv : NAN);
    printf("  bin %2d  %5ld eventos  %6.1f s  %.4f ± %.4f Hz\n", b, evBin[b], lv, lv > 0 ? evBin[b] / lv : NAN,
           lv > 0 ? std::sqrt(static_cast<double>(evBin[b])) / lv : NAN);
  }
  fclose(rb);

  // 2. Muón
  auto muon = [&](TH1D& h, const char* name) {
    double bg = 0;
    int nbg = 0;
    for (int k = 10; k < 40; ++k, ++nbg) bg += h.GetBinContent(k + 1);
    for (int k = 102; k < 140; ++k, ++nbg) bg += h.GetBinContent(k + 1);
    bg /= nbg;
    const double inWin = h.Integral(sigLo + 1, sigHi);
    printf("%s: %.0f canales×eventos; fondo %.2f por muestra; en [%d,%d): %.0f (exceso %.1f); moda %d\n", name, h.Integral(), bg, sigLo,
           sigHi, inWin, inWin - bg * (sigHi - sigLo), h.GetMaximumBin() - 1);
    printf("  cuentas por muestra 40-69:");
    for (int k = 40; k < 70; ++k) printf("%s%d:%.0f", k % 10 == 0 ? "\n   " : " ", k, h.GetBinContent(k + 1));
    printf("\n");
  };
  muon(hPos18, Form("hod2018 (> %.0f mV)", amp18));
  muon(hPos19, Form("hod2019 (> %.0f mV)", amp19));

  // 3. Fuera de tiempo
  const double T = ne * exposed;
  FILE* cc = fopen((std::string(out) + ".oot_channels.csv").c_str(), "w");
  fprintf(cc, "ch,station,thr_A_mV,count_A,rate_A_kHz,thr_B_mV,count_B,rate_B_kHz\n");
  std::vector<double> r18, r19a, r19b;
  for (int c = 0; c < 64; ++c) {
    if (c < 32) {
      if (!std::isfinite(G[c])) continue;
      const double r = ootAll18[c] / T / 1e3;
      fprintf(cc, "%d,hod2018,%.3f,%ld,%.3f,,,\n", c, 3 * G[c], ootAll18[c], r);
      if (c != 3) r18.push_back(r);
    } else {
      const double ra = ootAll19a[c] / T / 1e3, rb2 = ootAll19b[c] / T / 1e3;
      fprintf(cc, "%d,hod2019,3.000,%ld,%.3f,5.000,%ld,%.3f\n", c, ootAll19a[c], ra, ootAll19b[c], rb2);
      r19a.push_back(ra);
      r19b.push_back(rb2);
    }
  }
  fclose(cc);
  auto mm = [](std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return std::vector<double>{v[v.size() / 2], v.front(), v.back()};
  };
  const auto a = mm(r18), b19a = mm(r19a), b19b = mm(r19b);
  printf("fuera de tiempo (91 muestras, %.4f s por canal): hod2018 > 3 p.e. (G de v8) mediana %.2f kHz (mín %.2f, máx %.2f, %zu canales, "
         "sin ch 3); hod2019 > 3 mV mediana %.3f kHz (mín %.3f, máx %.3f); hod2019 > 5 mV mediana %.3f kHz (mín %.3f, máx %.3f)\n",
         T, a[0], a[1], a[2], r18.size(), b19a[0], b19a[1], b19a[2], b19b[0], b19b[1], b19b[2]);
  long tot19a = 0, tot19b = 0;
  for (int c = 32; c < 64; ++c) {
    tot19a += ootAll19a[c];
    tot19b += ootAll19b[c];
  }
  printf("hod2019 total (32 canales): > 3 mV %ld pulsos = %.3f kHz por canal; > 5 mV %ld pulsos = %.3f kHz por canal\n", tot19a,
         tot19a / (32 * T) / 1e3, tot19b, tot19b / (32 * T) / 1e3);
  FILE* bc = fopen((std::string(out) + ".oot_bins.csv").c_str(), "w");
  fprintf(bc, "bin,t_start_utc_ns,events,ch,count_A,rate_A_kHz,count_B,rate_B_kHz\n");
  printf("por bin: mediana hod2018 > 3 p.e. [kHz] | hod2019 media por canal > 3 mV, > 5 mV [kHz]\n");
  for (int b = 0; b < nbins; ++b) {
    if (evBin[b] == 0) continue;
    const double Tb = evBin[b] * exposed;
    std::vector<double> m18;
    long s19a = 0, s19b = 0;
    for (int c = 0; c < 64; ++c) {
      if (c < 32 && !std::isfinite(G[c])) continue;
      fprintf(bc, "%d,%lld,%ld,%d,%ld,%.3f,%ld,%.3f\n", b, tFirst + static_cast<Long64_t>(b * binS * 1e9), evBin[b], c, binA[b][c],
              binA[b][c] / Tb / 1e3, binB[b][c], binB[b][c] / Tb / 1e3);
      if (c < 32 && c != 3) m18.push_back(binA[b][c] / Tb / 1e3);
      if (c >= 32) {
        s19a += binA[b][c];
        s19b += binB[b][c];
      }
    }
    printf("  bin %2d  %5ld ev  %8.2f | %7.3f %7.3f\n", b, evBin[b], median(m18), s19a / (32 * Tb) / 1e3, s19b / (32 * Tb) / 1e3);
  }
  fclose(bc);
  fo.Close();
}
