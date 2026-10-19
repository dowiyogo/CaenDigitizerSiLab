// commissioning_pedestal.C — DAQ2-6 (solo lectura). Pedestal por canal de un run de daqrun.
// Para cada canal, sobre las muestras [lo, hi) de cada evento:
//   mean_mV      media de todas las muestras (conversión de la cabecera: 0x7FFF, 12 bit, 2 Vpp)
//   rms_mV       RMS por muestra alrededor de la media del propio evento (ruido)
//   sigma_tm_mV  desviación entre eventos de la media recortada (±2 mV de la mediana) de [lo, hi):
//                la precisión del estimador de pedestal (con HV apagada no hay pulsos)
// Escribe un CSV y un resumen: medianas por estación y canales atípicos (rms fuera de
// [0.5, 2] × la mediana de su estación).
// Uso: root -l -b -q 'daq/tools/validation/commissioning_pedestal.C("<dir>/<base>_seg*.root", 0, 40, "out.csv")'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
double adcToMv(double a) { return (a / 4095.0 - 1.0 + 32767.0 / 65535.0) * 2000.0; }
double med(std::vector<double> v) {
  if (v.empty()) return NAN;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}
}  // namespace

void commissioning_pedestal(const char* pattern, int lo, int hi, const char* csv = "", Long64_t maxEvents = 0) {
  TChain t("events");
  t.Add(pattern);
  // Longitud del registro desde la hoja adc[64][N].
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchAddress("adc", adc.data());
  const Long64_t ne = maxEvents > 0 ? std::min(maxEvents, t.GetEntries()) : t.GetEntries();
  std::vector<double> sum(64, 0), sumVar(64, 0), tmSum(64, 0), tmSum2(64, 0);
  const double lsb = daq::pulse::lsbMv();
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    for (int c = 0; c < 64; ++c) {
      const UShort_t* x = adc.data() + c * n;
      double s = 0, s2 = 0;
      for (int k = lo; k < hi; ++k) {
        s += x[k];
        s2 += double(x[k]) * x[k];
      }
      const double m = s / (hi - lo);
      sum[c] += m;
      sumVar[c] += s2 / (hi - lo) - m * m;
      const double tm = daq::pulse::trimmedMean(x, daq::pulse::Region{lo, hi}, 2.0 / lsb);
      tmSum[c] += tm;
      tmSum2[c] += tm * tm;
    }
  }
  FILE* out = csv && *csv ? fopen(csv, "w") : nullptr;
  if (out) fprintf(out, "ch,station,mean_mV,rms_mV,sigma_tm_mV\n");
  std::vector<double> rms(64), mean(64), stm(64);
  for (int c = 0; c < 64; ++c) {
    mean[c] = adcToMv(sum[c] / ne);
    rms[c] = std::sqrt(sumVar[c] / ne) * lsb;
    const double m = tmSum[c] / ne;
    stm[c] = std::sqrt(std::max(0.0, tmSum2[c] / ne - m * m)) * lsb;
    if (out) fprintf(out, "%d,%s,%.4f,%.4f,%.4f\n", c, c < 32 ? "hod2018" : "hod2019", mean[c], rms[c], stm[c]);
  }
  if (out) fclose(out);
  printf("eventos %lld, N = %d muestras, región [%d, %d)\n", ne, n, lo, hi);
  for (int st = 0; st < 2; ++st) {
    std::vector<double> r(rms.begin() + 32 * st, rms.begin() + 32 * st + 32), mn(mean.begin() + 32 * st, mean.begin() + 32 * st + 32),
        s(stm.begin() + 32 * st, stm.begin() + 32 * st + 32);
    const double mr = med(r);
    printf("%s: media mediana %.3f mV (mín %.3f, máx %.3f); RMS por muestra mediana %.4f mV (mín %.4f, máx %.4f); "
           "sigma media recortada mediana %.4f mV (mín %.4f, máx %.4f)\n",
           st ? "hod2019" : "hod2018", med(mn), *std::min_element(mn.begin(), mn.end()), *std::max_element(mn.begin(), mn.end()), mr,
           *std::min_element(r.begin(), r.end()), *std::max_element(r.begin(), r.end()), med(s), *std::min_element(s.begin(), s.end()),
           *std::max_element(s.begin(), s.end()));
    printf("   atípicos (RMS fuera de [0.5, 2] × mediana):");
    for (int c = 0; c < 32; ++c)
      if (r[c] < 0.5 * mr || r[c] > 2 * mr) printf(" ch%d=%.4f", 32 * st + c, r[c]);
    printf("\n");
  }
}
