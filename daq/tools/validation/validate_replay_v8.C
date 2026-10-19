// validate_replay_v8.C — validación 3 de DAQ2-5 (solo lectura).
// Lee un segmento de daqrun hecho con el replay de v8 y el archivo de origen (modo READ) y comprueba:
//  1. para 100 eventos al azar, que las 64×144 muestras recortadas coinciden con el origen dentro del
//     redondeo del ADC de la conversión (|mV_nuevo − mV_origen| ≤ LSB/2) y que el contador = `event`;
//  2. las tasas de pulsos fuera de tiempo por canal con el buscador nuevo (daq/core/include/daq/
//     pulse_finder.h: línea base = media recortada del pedestal, zonas pedestal+control+post = 91
//     muestras = 1456 ns, separación 4), hod2018 > 3 p.e. con las ganancias de v8 (solo referencia)
//     y hod2019 > 3 mV; y, para los mismos eventos, con el método del estudio de DAQ2-1 sobre el
//     registro antiguo (referencia = mediana de [20,200), zonas [209,258)+[311,353)).
// Uso (desde la raíz del repo):
//   root -l -b -q 'daq/tools/validation/validate_replay_v8.C("<segmento>.root", 30000, "outputs/x.csv")'
R__ADD_INCLUDE_PATH(daq/core/include)
R__ADD_INCLUDE_PATH(analysis/src)
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#include "daq/pulse_finder.h"
#include "hodo_paths.h"

namespace {
double adcToMv(double adc) { return (adc / 4095.0 - 1.0 + 32767.0 / 65535.0) * 2000.0; }
double medianOf(std::vector<double> v) {
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}
}  // namespace

void validate_replay_v8(const char* segment, long long firstEntry = 30000, const char* outCsv = "",
                        const char* source = "lectura_vx1740_ext_nim_1M_v8.root",
                        const char* gains = "outputs/pe_calibration_v8.txt") {
  const int N = 144, kOld0 = 209, nCh = 64;
  const double lsb = daq::pulse::lsbMv();
  // Ganancias de v8 (solo como referencia de umbral para hod2018).
  std::vector<double> G(nCh, NAN);
  {
    std::ifstream in(gains);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int ch;
      std::string pl;
      double p0, g;
      if (ss >> ch >> pl >> p0 >> g && ch >= 0 && ch < nCh) G[ch] = g;
    }
  }
  TFile* fn = TFile::Open(segment, "READ");
  TTree* tn = (TTree*)fn->Get("events");
  ULong64_t seq = 0;
  UInt_t counter = 0;
  static UShort_t adc[64][144];
  tn->SetBranchAddress("evt_seq", &seq);
  tn->SetBranchAddress("evt_counter_hw", &counter);
  tn->SetBranchAddress("adc", adc);
  const std::string srcPath = hodo::resolveRawFile(source);
  TFile* fs = TFile::Open(srcPath.c_str(), "READ");
  TTree* ts = (TTree*)fs->Get("data_dig");
  ts->SetBranchStatus("*", 0);
  std::vector<std::vector<float>*> ch(nCh, nullptr);
  int event = 0;
  for (int c = 0; c < nCh; ++c) {
    ts->SetBranchStatus(Form("Ch%d", c), 1);
    ts->SetBranchAddress(Form("Ch%d", c), &ch[c]);
  }
  ts->SetBranchStatus("event", 1);
  ts->SetBranchAddress("event", &event);
  const Long64_t n = tn->GetEntries();
  printf("segmento %s: %lld eventos; origen %s desde la entrada %lld\n", segment, n, srcPath.c_str(), firstEntry);

  // 1. Muestras de 100 eventos al azar.
  TRandom3 rnd(12345);
  long bad = 0, compared = 0, badCounter = 0;
  double maxDiff = 0;
  for (int i = 0; i < 100; ++i) {
    const Long64_t e = static_cast<Long64_t>(rnd.Uniform(0, n));
    tn->GetEntry(e);
    ts->GetEntry(firstEntry + static_cast<Long64_t>(seq));
    if (counter != static_cast<UInt_t>(event)) ++badCounter;
    for (int c = 0; c < nCh; ++c)
      for (int k = 0; k < N; ++k) {
        const double d = std::fabs(adcToMv(adc[c][k]) - (*ch[c])[kOld0 + k]);
        maxDiff = std::max(maxDiff, d);
        ++compared;
        if (d > lsb / 2 + 1e-4) ++bad;
      }
  }
  printf("1) 100 eventos al azar: %ld muestras comparadas, %ld fuera de LSB/2 (%.4f mV), diferencia máxima %.6f mV; "
         "contador distinto de event: %ld\n",
         compared, bad, lsb / 2, maxDiff, badCounter);

  // 2. Tasas fuera de tiempo.
  const std::vector<daq::pulse::Region> zonesNew{{0, 49}, {102, 144}};
  std::vector<long> cNew(nCh, 0), cOld(nCh, 0);
  for (Long64_t e = 0; e < n; ++e) {
    tn->GetEntry(e);
    ts->GetEntry(firstEntry + static_cast<Long64_t>(seq));
    for (int c = 0; c < nCh; ++c) {
      const bool h18 = c < 32;
      if (h18 && !std::isfinite(G[c])) continue;
      const double thr = h18 ? 3.0 * G[c] : 3.0;
      const double bl = daq::pulse::trimmedMean(adc[c], daq::pulse::Region{0, 40}, 2.0 / lsb);
      cNew[c] += daq::pulse::countPulses(adc[c], N, bl, lsb, thr, 4, zonesNew);
      const std::vector<float>& w = *ch[c];
      std::vector<double> ref(w.begin() + 20, w.begin() + 200);
      const double r = medianOf(ref);
      cOld[c] += daq::pulse::countPulses(w.data(), static_cast<int>(w.size()), r, 1.0, thr, 4, {{209, 258}, {311, 353}});
    }
  }
  const double T = n * 91 * 16e-9;  // s expuestos por canal
  std::vector<double> rN18, rO18, rN19, rO19;
  FILE* out = outCsv && *outCsv ? fopen(outCsv, "w") : nullptr;
  if (out) fprintf(out, "ch,station,threshold_mV,count_new,count_daq2_1_method,rate_new_kHz,rate_daq2_1_method_kHz\n");
  for (int c = 0; c < nCh; ++c) {
    const bool h18 = c < 32;
    if (h18 && (!std::isfinite(G[c]) || c == 3)) continue;  // ch3 muerto (como en DAQ2-1)
    const double a = cNew[c] / T / 1e3, b = cOld[c] / T / 1e3;
    (h18 ? rN18 : rN19).push_back(a);
    (h18 ? rO18 : rO19).push_back(b);
    if (out) fprintf(out, "%d,%s,%.4f,%ld,%ld,%.3f,%.3f\n", c, h18 ? "hod2018" : "hod2019", h18 ? 3 * G[c] : 3.0, cNew[c], cOld[c], a, b);
  }
  if (out) fclose(out);
  auto summary = [](const char* name, std::vector<double> v) {
    std::sort(v.begin(), v.end());
    printf("   %s: mediana %.2f kHz (mín %.2f, máx %.2f), %zu canales\n", name, v[v.size() / 2], v.front(), v.back(), v.size());
  };
  printf("2) tasas fuera de tiempo, %lld eventos, %.4f s expuestos por canal (91 muestras/evento):\n", n, T);
  summary("hod2018 > 3 p.e., buscador nuevo (base = media recortada del pedestal)", rN18);
  summary("hod2018 > 3 p.e., método DAQ2-1 (base = mediana [20,200) del registro antiguo)", rO18);
  summary("hod2019 > 3 mV, buscador nuevo", rN19);
  summary("hod2019 > 3 mV, método DAQ2-1", rO19);
  for (auto* v : ch) delete v;
}
