// hod19_ped.C — seguimiento 2026-10-31 (5), punto 1 (solo lectura): pedestal y ruido por canal (64 canales), para ver si la lectura de hod2019
// funcionaba en octubre.
//  (a) Sobre la ntupla de ntuple_extract.C (todos los eventos, o los de tiempo bueno con un archivo de intervalos UTC excluidos): por canal, la
//      media y los cuantiles 16/50/84 de ped (media recortada ±2 mV de las muestras [0,40) del tramo, pre-trigger) y de rms (RMS recortado ±2 mV).
//  (b) Sobre el crudo, en nev eventos a partir de la entrada start que estén en tiempo bueno: espectro de la desviación muestra a muestra respecto
//      del pedestal del evento (x[s] − media recortada de [0,40) del tramo, s en [0,40) del tramo), bins de 0.25 mV en [−20, +20) mV, más el
//      RMS completo (sin recortar) de esas muestras; y la diferencia entre muestras consecutivas x[s+1] − x[s] (mismo rango).
// Formatos: "daqv2" (events, adc[64][n], tramo = registro desde la muestra 0) y "may" (data_dig, Ch0..Ch63 en mV, tramo desde la muestra 209).
// Salida: <prefijo>_ped.csv (ch, n, ped_mean, ped_q16, ped_q50, ped_q84, rms_q16, rms_q50, rms_q84, rms_full_mean, dev_frac_lt_m3mV, dev_frac_gt_p3mV)
//         y <prefijo>_dev.csv (ch, bin_lo_mV, n_dev, n_diff).
// Uso: root -l -b -q 'hod19_ped.C+("<ntupla.root>", "<crudo(s)>", "daqv2"|"may", "<excluidos.txt o vacío>", start, nev, "<prefijo>")'
#include <TChain.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <utility>
#include <vector>

namespace {
std::vector<std::pair<double, double>> readExcl(const char* fn) {
  std::vector<std::pair<double, double>> v;
  if (!fn || !*fn) return v;
  std::ifstream f(fn);
  double a, b;
  while (f >> a >> b) v.push_back({a, b});
  return v;
}
bool good(double t, const std::vector<std::pair<double, double>>& ex) {
  for (auto& [a, b] : ex)
    if (t >= a && t < b) return false;
  return true;
}
double q(std::vector<float>& v, double p) {
  if (v.empty()) return NAN;
  const size_t k = std::min(v.size() - 1, static_cast<size_t>(p * v.size()));
  std::nth_element(v.begin(), v.begin() + k, v.end());
  return v[k];
}
double trimmedMean(const double* x, int n, double trim) {
  std::vector<double> s(x, x + n);
  std::nth_element(s.begin(), s.begin() + n / 2, s.end());
  const double m = s[n / 2];
  double sum = 0;
  int k = 0;
  for (int i = 0; i < n; ++i)
    if (std::fabs(x[i] - m) <= trim) { sum += x[i]; ++k; }
  return k ? sum / k : m;
}
}  // namespace

void hod19_ped(const char* ntuple, const char* raw, const char* format, const char* excl, Long64_t start, int nev, const char* prefix) {
  const auto ex = readExcl(excl);
  const bool may = std::string(format) == "may";
  const int off = may ? 209 : 0;
  // (a) ntupla
  std::vector<std::vector<float>> P(64), R(64);
  std::set<Long64_t> goodEntries;
  {
    TFile f(ntuple, "READ");
    TTree* t = (TTree*)f.Get("ev");
    double tu;
    Long64_t evn;
    float ped[64], rms[64];
    t->SetBranchStatus("*", 0);
    for (const char* b : {"t_utc_s", "ev", "ped", "rms"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("t_utc_s", &tu);
    t->SetBranchAddress("ev", &evn);
    t->SetBranchAddress("ped", ped);
    t->SetBranchAddress("rms", rms);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!good(tu, ex)) continue;
      goodEntries.insert(evn);
      for (int c = 0; c < 64; ++c) {
        P[c].push_back(ped[c]);
        R[c].push_back(rms[c]);
      }
    }
  }
  // (b) crudo
  const int NB = 160;  // 0.25 mV en [−20, 20)
  std::vector<std::vector<long>> hdev(64, std::vector<long>(NB, 0)), hdif(64, std::vector<long>(NB, 0));
  std::vector<double> s2(64, 0), lt(64, 0), gt(64, 0);
  std::vector<long> ns(64, 0);
  TChain ch(may ? "data_dig" : "events");
  ch.Add(raw);
  std::vector<std::vector<float>*> cw(64, nullptr);
  std::vector<UShort_t> adc;
  int len = 0;
  if (may) {
    ch.SetBranchStatus("*", 0);
    for (int c = 0; c < 64; ++c) {
      ch.SetBranchStatus(Form("Ch%d", c), 1);
      ch.SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
  } else {
    len = ch.GetLeaf("adc")->GetLenStatic() / 64;
    adc.resize(64 * len);
    ch.SetBranchStatus("*", 0);
    ch.SetBranchStatus("adc", 1);
    ch.SetBranchAddress("adc", adc.data());
  }
  const double lsb = 2000.0 / 4095.0;
  int used = 0;
  double x[40];
  for (Long64_t e = start; e < ch.GetEntries() && used < nev; ++e) {
    if (!goodEntries.count(e)) continue;
    ch.GetEntry(e);
    ++used;
    for (int c = 0; c < 64; ++c) {
      for (int s = 0; s < 40; ++s) x[s] = may ? (*cw[c])[off + s] : adc[c * len + s] * lsb;
      const double bl = trimmedMean(x, 40, 2.0);
      for (int s = 0; s < 40; ++s) {
        const double d = x[s] - bl;
        s2[c] += d * d;
        ++ns[c];
        if (d < -3) ++lt[c];
        if (d > 3) ++gt[c];
        const int b = static_cast<int>(std::floor((d + 20) / 0.25));
        if (b >= 0 && b < NB) ++hdev[c][b];
        if (s < 39) {
          const double dd = x[s + 1] - x[s];
          const int b2 = static_cast<int>(std::floor((dd + 20) / 0.25));
          if (b2 >= 0 && b2 < NB) ++hdif[c][b2];
        }
      }
    }
  }
  FILE* fo = fopen((std::string(prefix) + "_ped.csv").c_str(), "w");
  fprintf(fo, "ch,n_events_ntuple,ped_mean,ped_q16,ped_q50,ped_q84,rms_q16,rms_q50,rms_q84,n_events_raw,rms_full,dev_frac_lt_m3mV,dev_frac_gt_p3mV\n");
  for (int c = 0; c < 64; ++c) {
    double m = 0;
    for (float v : P[c]) m += v;
    m /= std::max<size_t>(1, P[c].size());
    const size_t n = P[c].size();
    fprintf(fo, "%d,%zu,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.6f,%.6f\n", c, n, m, q(P[c], 0.16), q(P[c], 0.5), q(P[c], 0.84), q(R[c], 0.16), q(R[c], 0.5),
            q(R[c], 0.84), used, ns[c] ? std::sqrt(s2[c] / ns[c]) : NAN, ns[c] ? lt[c] / ns[c] : NAN, ns[c] ? gt[c] / ns[c] : NAN);
  }
  fclose(fo);
  FILE* fd = fopen((std::string(prefix) + "_dev.csv").c_str(), "w");
  fprintf(fd, "ch,bin_lo_mV,n_dev,n_diff\n");
  for (int c = 0; c < 64; ++c)
    for (int b = 0; b < NB; ++b) fprintf(fd, "%d,%.2f,%ld,%ld\n", c, -20 + 0.25 * b, hdev[c][b], hdif[c][b]);
  fclose(fd);
  printf("ntupla: %zu eventos de tiempo bueno; crudo: %d eventos desde la entrada %lld\n", P[0].size(), used, start);
}
