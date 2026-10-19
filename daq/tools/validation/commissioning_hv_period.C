// commissioning_hv_period.C — DAQ2-6, paso 4 (solo lectura). IMon y VMon por canal HV (CH00-07) en
// bins de bin_s (600 s) sobre un período [from, to] (UTC, "YYYY-MM-DDTHH:MM:SSZ") desde los
// archivos de control lento de hvmon (hv y hv_full, lector libdaq_screader.so), con la temperatura
// de la placa A7030DP (bd_Temp_degC de hv_full; sigue la del rack, no la del detector).
// Si se da oot_bins (CSV <out>.oot_bins.csv de commissioning_ext_runs.C con los mismos bin_s y
// origen = primer evento), calcula por grupo HV la tasa fuera de tiempo por canal en cada bin
// (hod2018 > 3·G de v8 y hod2019 > 3 mV: columna A; hod2019 > 5 mV: columna B) y el coeficiente de
// Pearson con el IMon del canal HV en los mismos bins. Canales de grupo desconocido no entran.
// Uso: root -l -b -q 'daq/tools/validation/commissioning_hv_period.C("20261019", "2026-10-19T16:30:00Z",
//        "2026-10-19T21:00:00Z", "out", 600, "outputs/daq2-6/s4a_runs78.oot_bins.csv")'
R__ADD_INCLUDE_PATH(daq/io/include)
R__LOAD_LIBRARY(daq/build/io/libdaq_screader.so)
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/io/sc_reader.h"

namespace {
Long64_t parseUtc(const char* s) {
  std::tm tm{};
  if (!strptime(s, "%Y-%m-%dT%H:%M:%SZ", &tm)) return 0;
  return static_cast<Long64_t>(timegm(&tm)) * 1000000000LL;
}
std::string fmtUtc(Long64_t ns) {
  const time_t s = ns / 1000000000LL;
  char b[32];
  strftime(b, sizeof b, "%H:%MZ", gmtime(&s));
  return b;
}
double pearson(const std::vector<double>& a, const std::vector<double>& b) {
  const size_t n = a.size();
  if (n < 3) return NAN;
  double ma = 0, mb = 0;
  for (size_t i = 0; i < n; ++i) {
    ma += a[i];
    mb += b[i];
  }
  ma /= n;
  mb /= n;
  double sab = 0, saa = 0, sbb = 0;
  for (size_t i = 0; i < n; ++i) {
    sab += (a[i] - ma) * (b[i] - mb);
    saa += (a[i] - ma) * (a[i] - ma);
    sbb += (b[i] - mb) * (b[i] - mb);
  }
  return sab / std::sqrt(saa * sbb);
}
}  // namespace

void commissioning_hv_period(const char* day, const char* from, const char* to, const char* out, double binS = 600.0,
                             const char* ootBins = "", const char* mappingFile = "daq/config/mapping/hodo_mapping_2026_10.json",
                             const char* scDir = "/data/daq/slowcontrol") {
  const Long64_t t0 = parseUtc(from), t1 = parseUtc(to);
  const int nbins = static_cast<int>((t1 - t0) / 1e9 / binS) + 1;
  std::vector<std::vector<double>> iS(8, std::vector<double>(nbins, 0)), vS = iS, nS = iS;
  std::vector<double> tS(nbins, 0), tN(nbins, 0), v0set(8, NAN), i0set(8, NAN);
  std::vector<double> iMin(8, 1e9), iMax(8, -1e9);
  std::string err;
  TTree* sc = daq::io::readSlowControl(scDir, "hv", day, day, &err);
  if (!sc) {
    printf("hv: %s\n", err.c_str());
    return;
  }
  Long64_t tr = 0;
  Int_t ch = 0;
  Double_t imon = 0, vmon = 0, v0 = 0, i0 = 0;
  UInt_t pw = 0;
  sc->SetBranchAddress("t_resp_utc_ns", &tr);
  sc->SetBranchAddress("ch", &ch);
  sc->SetBranchAddress("IMon_A", &imon);
  sc->SetBranchAddress("VMon_V", &vmon);
  sc->SetBranchAddress("V0Set_V", &v0);
  sc->SetBranchAddress("I0Set_A", &i0);
  sc->SetBranchAddress("Pw", &pw);
  long offRows = 0;
  for (Long64_t i = 0; i < sc->GetEntries(); ++i) {
    sc->GetEntry(i);
    if (ch < 0 || ch >= 8 || tr < t0 || tr >= t1) continue;
    if (pw != 1) ++offRows;
    const int b = static_cast<int>((tr - t0) / 1e9 / binS);
    iS[ch][b] += imon;
    vS[ch][b] += vmon;
    nS[ch][b] += 1;
    v0set[ch] = v0;
    i0set[ch] = i0;
    iMin[ch] = std::min(iMin[ch], imon);
    iMax[ch] = std::max(iMax[ch], imon);
  }
  delete sc;
  TTree* full = daq::io::readSlowControl(scDir, "hv_full", day, day, &err);
  double tempMin = 1e9, tempMax = -1e9;
  if (full) {
    Double_t temp = NAN;
    full->SetBranchAddress("t_resp_utc_ns", &tr);
    full->SetBranchAddress("ch", &ch);
    full->SetBranchAddress("bd_Temp_degC", &temp);
    for (Long64_t i = 0; i < full->GetEntries(); ++i) {
      full->GetEntry(i);
      if (ch != 0 || tr < t0 || tr >= t1 || !std::isfinite(temp)) continue;
      const int b = static_cast<int>((tr - t0) / 1e9 / binS);
      tS[b] += temp;
      tN[b] += 1;
      tempMin = std::min(tempMin, temp);
      tempMax = std::max(tempMax, temp);
    }
    delete full;
  } else {
    printf("hv_full: %s\n", err.c_str());
  }
  printf("período %s -> %s, %d bins de %.0f s; filas con Pw != 1 en CH00-07: %ld\n", from, to, nbins, binS, offRows);
  printf("temperatura de la placa A7030DP: %.1f-%.1f degC\n", tempMin, tempMax);
  printf("CH  V0Set_V  I0Set_uA  IMon_uA mín-máx (lecturas)   margen I0Set/IMon_máx   VMon-V0Set_V mín-máx (bins)\n");
  for (int h = 0; h < 8; ++h) {
    double dmin = 1e9, dmax = -1e9;
    for (int b = 0; b < nbins; ++b)
      if (nS[h][b] > 0) {
        const double d = vS[h][b] / nS[h][b] - v0set[h];
        dmin = std::min(dmin, d);
        dmax = std::max(dmax, d);
      }
    printf("%02d  %7.3f  %8.3f  %6.3f-%6.3f            %5.2f                    %6.3f-%6.3f\n", h, v0set[h], i0set[h] * 1e6, iMin[h] * 1e6,
           iMax[h] * 1e6, i0set[h] / iMax[h], dmin, dmax);
  }
  FILE* csv = fopen((std::string(out) + ".csv").c_str(), "w");
  fprintf(csv, "bin,t_start_utc_ns,board_temp_degC");
  for (int h = 0; h < 8; ++h) fprintf(csv, ",imon_ch%02d_uA,vmon_ch%02d_V", h, h);
  fprintf(csv, "\n");
  printf("bin  inicio  T_placa  IMon_uA CH00..CH07\n");
  for (int b = 0; b < nbins; ++b) {
    fprintf(csv, "%d,%lld,%.2f", b, t0 + static_cast<Long64_t>(b * binS * 1e9), tN[b] > 0 ? tS[b] / tN[b] : NAN);
    printf("%3d  %s  %5.2f  ", b, fmtUtc(t0 + static_cast<Long64_t>(b * binS * 1e9)).c_str(), tN[b] > 0 ? tS[b] / tN[b] : NAN);
    for (int h = 0; h < 8; ++h) {
      const double im = nS[h][b] > 0 ? iS[h][b] / nS[h][b] * 1e6 : NAN, vm = nS[h][b] > 0 ? vS[h][b] / nS[h][b] : NAN;
      fprintf(csv, ",%.4f,%.3f", im, vm);
      printf(" %6.3f", im);
    }
    printf("\n");
    fprintf(csv, "\n");
  }
  fclose(csv);
  if (!ootBins || !*ootBins) return;
  // Correlación con las tasas fuera de tiempo por grupo.
  std::ifstream mi(mappingFile);
  const nlohmann::json mj = nlohmann::json::parse(mi);
  std::vector<int> hvOf(64, -1);
  for (const auto& c : mj["channels"])
    if (!c["hv_group"].is_null()) hvOf[c["dig_ch"].get<int>()] = c["hv_group"].get<int>() - 1;
  std::ifstream ob(ootBins);
  std::string line;
  std::getline(ob, line);
  // (bin de ext_runs) -> t_start, eventos; conteos por grupo y estación
  std::map<int, Long64_t> binStart;
  std::map<int, long> binEv;
  std::map<int, std::vector<double>> a18, a19, b19;  // conteos por grupo
  std::vector<int> n18(8, 0), n19(8, 0);
  std::vector<bool> seen(64, false);
  while (std::getline(ob, line)) {
    std::istringstream ss(line);
    std::string f[8];
    for (auto& x : f) std::getline(ss, x, ',');
    const int b = std::stoi(f[0]), c = std::stoi(f[3]);
    binStart[b] = std::stoll(f[1]);
    binEv[b] = std::stol(f[2]);
    const int h = hvOf[c];
    if (h < 0 || h >= 8) continue;
    if (!seen[c]) {
      seen[c] = true;
      (c < 32 ? n18 : n19)[h]++;
    }
    for (auto* m : {&a18, &a19, &b19})
      if (!m->count(b)) (*m)[b] = std::vector<double>(8, 0);
    if (c < 32)
      a18[b][h] += std::stod(f[4]);
    else {
      a19[b][h] += std::stod(f[4]);
      b19[b][h] += std::stod(f[6]);
    }
  }
  printf("\ncorrelación por grupo HV (bins de commissioning_ext_runs; IMon del mismo intervalo):\n");
  printf("grupo  canales18 canales19   tasa18_kHz ini/máx   r18    tasa19>3mV_kHz ini/máx  r19_3   tasa19>5mV máx  r19_5\n");
  for (int h = 0; h < 8; ++h) {
    std::vector<double> I, R18, R19a, R19b;
    for (const auto& [b, ts] : binStart) {
      if (binEv[b] < 100) continue;  // bin parcial al final
      double is = 0, in = 0;
      // IMon medio en [ts, ts + bin_s) con las filas ya agregadas por bin del período: se recalcula
      // con el bin del período que contiene el centro del bin de ext_runs.
      const Long64_t mid = ts + static_cast<Long64_t>(binS * 0.5e9);
      const int pb = static_cast<int>((mid - t0) / 1e9 / binS);
      if (pb < 0 || pb >= nbins || nS[h][pb] == 0) continue;
      is = iS[h][pb];
      in = nS[h][pb];
      const double T = binEv[b] * 91 * 16e-9;
      I.push_back(is / in * 1e6);
      R18.push_back(n18[h] ? a18[b][h] / (n18[h] * T) / 1e3 : NAN);
      R19a.push_back(n19[h] ? a19[b][h] / (n19[h] * T) / 1e3 : NAN);
      R19b.push_back(n19[h] ? b19[b][h] / (n19[h] * T) / 1e3 : NAN);
    }
    auto mx = [](const std::vector<double>& v) { return v.empty() ? NAN : *std::max_element(v.begin(), v.end()); };
    printf("%d (CH%02d)  %3d       %3d        %8.2f / %8.2f   %6.3f   %7.3f / %7.3f        %6.3f   %7.3f        %6.3f\n", h + 1, h, n18[h],
           n19[h], R18.empty() ? NAN : R18[0], mx(R18), n18[h] ? pearson(I, R18) : NAN, R19a.empty() ? NAN : R19a[0], mx(R19a),
           n19[h] ? pearson(I, R19a) : NAN, mx(R19b), n19[h] ? pearson(I, R19b) : NAN);
  }
}
