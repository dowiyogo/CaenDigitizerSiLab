// commissioning_hv_vs_rate.C — DAQ2-6, paso 4 (solo lectura). IMon y VMon por canal HV contra el
// tiempo (archivos de control lento de hvmon, lector libdaq_screader.so) durante un run, y su
// correlación con la tasa de pulsos fuera de tiempo por grupo HV (canales del digitizador del grupo
// según el mapeo JSON; los de grupo desconocido no entran).
//  - Bins de bin_s (60 s) entre el primer y el último evento del run.
//  - Tasa fuera de tiempo por canal: zonas de la ventana B [0,49)+[102,144) (91 muestras) en
//    registros de 144; en registros completos [20,250)+[280,N−2). Umbral hod2018 3·G (G de v8, solo
//    referencia), hod2019 3 mV. Por grupo: suma de pulsos / (eventos × tiempo expuesto × canales).
//  - Por canal HV: IMon y VMon medios por bin, V0Set; resumen (primero/último bin, mín, máx,
//    VMon − V0Set) y coeficiente de correlación de Pearson entre IMon y la tasa del grupo en los bins.
// Salidas: CSV por bin (out.csv) y un ROOT con TGraph (out.root).
// Uso: root -l -b -q 'daq/tools/validation/commissioning_hv_vs_rate.C("<run>_seg*.root", "20261019", "out")'
R__ADD_INCLUDE_PATH(daq/core/include)
R__ADD_INCLUDE_PATH(daq/io/include)
R__LOAD_LIBRARY(daq/build/io/libdaq_screader.so)
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/io/sc_reader.h"
#include "daq/pulse_finder.h"

namespace {
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

void commissioning_hv_vs_rate(const char* pattern, const char* day, const char* out, double binS = 60.0,
                              const char* mappingFile = "daq/config/mapping/hodo_mapping_2026_10.json",
                              const char* gains = "outputs/pe_calibration_v8.txt", const char* scDir = "/data/daq/slowcontrol") {
  // Mapeo: canal del digitizador -> canal HV.
  std::ifstream mi(mappingFile);
  const nlohmann::json mj = nlohmann::json::parse(mi);
  std::vector<int> hvOf(64, -1);
  for (const auto& c : mj["channels"])
    if (!c["hv_group"].is_null()) hvOf[c["dig_ch"].get<int>()] = c["hv_group"].get<int>() - 1;
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
  std::vector<UShort_t> adc(64 * n);
  Long64_t tutc = 0;
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchStatus("t_utc_ns", 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("t_utc_ns", &tutc);
  const Long64_t ne = t.GetEntries();
  t.GetEntry(0);
  const Long64_t t0 = tutc;
  t.GetEntry(ne - 1);
  const Long64_t t1 = tutc;
  const int nbins = std::max(1, static_cast<int>((t1 - t0) / 1e9 / binS) + 1);
  const bool windowB = n == 144;
  const std::vector<daq::pulse::Region> zones = windowB ? std::vector<daq::pulse::Region>{{0, 49}, {102, 144}}
                                                        : std::vector<daq::pulse::Region>{{20, 250}, {280, n - 2}};
  double exposedPerEvent = 0;
  for (const auto& z : zones) exposedPerEvent += (z.end - z.begin) * 16e-9;
  const daq::pulse::Region ped = windowB ? daq::pulse::Region{0, 40} : daq::pulse::Region{20, 200};
  const double lsb = daq::pulse::lsbMv();
  std::vector<std::vector<double>> cnt(8, std::vector<double>(nbins, 0));
  std::vector<int> nch(8, 0);
  for (int c = 0; c < 64; ++c)
    if (hvOf[c] >= 0 && hvOf[c] < 8 && (c >= 32 || std::isfinite(G[c]))) ++nch[hvOf[c]];
  std::vector<double> evBin(nbins, 0);
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    const int b = std::min(nbins - 1, static_cast<int>((tutc - t0) / 1e9 / binS));
    evBin[b] += 1;
    for (int c = 0; c < 64; ++c) {
      const int h = hvOf[c];
      if (h < 0 || h >= 8) continue;
      const double thr = c < 32 ? 3.0 * G[c] : 3.0;
      if (!std::isfinite(thr)) continue;
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, ped, 2.0 / lsb);
      cnt[h][b] += daq::pulse::countPulses(x, n, bl, lsb, thr, 4, zones);
    }
  }
  // Control lento: hv_<día>.
  std::string err;
  TTree* sc = daq::io::readSlowControl(scDir, "hv", day, day, &err);
  if (!sc) {
    printf("control lento: %s\n", err.c_str());
    return;
  }
  Long64_t tr = 0;
  Int_t ch = 0;
  Double_t imon = 0, vmon = 0, v0 = 0;
  sc->SetBranchAddress("t_resp_utc_ns", &tr);
  sc->SetBranchAddress("ch", &ch);
  sc->SetBranchAddress("IMon_A", &imon);
  sc->SetBranchAddress("VMon_V", &vmon);
  sc->SetBranchAddress("V0Set_V", &v0);
  std::vector<std::vector<double>> iS(8, std::vector<double>(nbins, 0)), vS(8, std::vector<double>(nbins, 0)), nS(8, std::vector<double>(nbins, 0));
  std::vector<double> v0set(8, 0);
  TFile fo((std::string(out) + ".root").c_str(), "RECREATE");
  std::vector<TGraph> gi(8), gv(8);
  for (Long64_t i = 0; i < sc->GetEntries(); ++i) {
    sc->GetEntry(i);
    if (ch < 0 || ch >= 8 || tr < t0 || tr > t1 + static_cast<Long64_t>(binS * 1e9)) continue;
    const int b = std::min(nbins - 1, static_cast<int>((tr - t0) / 1e9 / binS));
    iS[ch][b] += imon;
    vS[ch][b] += vmon;
    nS[ch][b] += 1;
    v0set[ch] = v0;
    gi[ch].AddPoint((tr - t0) / 60e9, imon * 1e6);
    gv[ch].AddPoint((tr - t0) / 60e9, vmon - v0);
  }
  for (int h = 0; h < 8; ++h) {
    gi[h].SetName(Form("imon_ch%02d", h));
    gi[h].SetTitle(Form("CH%02d IMon;min desde el primer evento;IMon [#muA]", h));
    gi[h].Write();
    gv[h].SetName(Form("vmon_minus_v0set_ch%02d", h));
    gv[h].SetTitle(Form("CH%02d VMon - V0Set;min desde el primer evento;V", h));
    gv[h].Write();
  }
  FILE* csv = fopen((std::string(out) + ".csv").c_str(), "w");
  fprintf(csv, "bin,t_start_utc_ns,events,hv_ch,imon_uA,vmon_V,v0set_V,oot_rate_kHz_per_channel\n");
  printf("run: %lld eventos, %d bins de %.0f s, registro de %d muestras (%s), %.0f ns expuestos por evento\n", ne, nbins, binS, n,
         windowB ? "ventana B" : "completo", exposedPerEvent * 1e9);
  printf("CH  canales  IMon_uA ini->fin (mín-máx)   VMon-V0Set_V ini->fin   tasa_kHz ini->fin (mediana)   r(IMon,tasa)\n");
  for (int h = 0; h < 8; ++h) {
    std::vector<double> I, R;
    double imin = 1e9, imax = -1e9;
    std::vector<double> rates;
    double iFirst = NAN, iLast = NAN, dFirst = NAN, dLast = NAN, rFirst = NAN, rLast = NAN;
    for (int b = 0; b < nbins; ++b) {
      if (nS[h][b] == 0 || evBin[b] == 0 || nch[h] == 0) continue;
      const double im = iS[h][b] / nS[h][b] * 1e6, vm = vS[h][b] / nS[h][b];
      const double rate = cnt[h][b] / (evBin[b] * exposedPerEvent * nch[h]) / 1e3;
      fprintf(csv, "%d,%lld,%.0f,%d,%.4f,%.3f,%.3f,%.3f\n", b, t0 + static_cast<Long64_t>(b * binS * 1e9), evBin[b], h, im, vm, v0set[h], rate);
      I.push_back(im);
      R.push_back(rate);
      rates.push_back(rate);
      imin = std::min(imin, im);
      imax = std::max(imax, im);
      if (std::isnan(iFirst)) {
        iFirst = im;
        dFirst = vm - v0set[h];
        rFirst = rate;
      }
      iLast = im;
      dLast = vm - v0set[h];
      rLast = rate;
    }
    std::sort(rates.begin(), rates.end());
    printf("%02d  %5d   %6.3f -> %6.3f (%6.3f-%6.3f)   %6.3f -> %6.3f           %7.2f -> %7.2f (%7.2f)       %6.3f\n", h, nch[h], iFirst,
           iLast, imin, imax, dFirst, dLast, rFirst, rLast, rates.empty() ? NAN : rates[rates.size() / 2], pearson(I, R));
  }
  fclose(csv);
  fo.Close();
  delete sc;
}
