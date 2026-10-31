// vbd_events.C — seguimiento 2026-10-31 (9), punto 1 (solo lectura de un CSV): cuántos eventos SW hacen falta para que el ajuste de G de c1
// (daq/core/include/daq/gain_fit.h, opciones por omisión: picos entre 5 y 26 mV, ≥ 3 picos de ≥ 20 cuentas en el espectro suavizado) converja
// en un escalón de tensión baja. Toma un espectro medido de pulsos oscuros (bins de 1 LSB, salida de gain_spectra.C/vbd_gain.C: ch,bin_lo_mV,n)
// con N0 eventos, lo escala a f·N0 eventos (cada bin ~ Poisson(f·n), réplicas con semillas fijas) y cuenta los canales con el ajuste aceptado
// (ok/changed) y la dispersión de G entre réplicas.
// Semilla del ajuste: G de c1 × seed_scale (la G esperada en ese escalón).
// Salida: CSV a stdout (factor, events, replica, n_channels, n_accepted, G_median_mV, n_insufficient, n_rejected, n_ambiguous) y, por canal,
// <prefijo>_perch.csv.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'vbd_events.C+("<spec.csv>", N0, "<c1.json>", seed_scale, "<factores coma>", replicas, "<prefijo>"[, "lo,hi,tol,min_peaks"])'
// opts vacío = opciones de c1; si no, la variante para G baja (p. ej. "2.0,14,0.3,3").
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/gain_fit.h"

void vbd_events(const char* specCsv, long n0, const char* c1Json, double seedScale, const char* factors, int replicas, const char* prefix,
                const char* opts = "") {
  daq::gain::Options o;  // opciones de c1 salvo que opts = "lo,hi,tol,min_peaks" (variante para G baja)
  if (opts && *opts) {
    double lo, hi, tol;
    int mp;
    if (sscanf(opts, "%lf,%lf,%lf,%d", &lo, &hi, &tol, &mp) == 4) {
      o.loMv = lo;
      o.hiMv = hi;
      o.spacingTolerance = tol;
      o.minPeaks = mp;
    }
  }
  std::map<int, std::vector<double>> sp;
  double bin = 0;
  {
    std::ifstream f(specCsv);
    std::string line;
    std::getline(f, line);
    double prev = -1;
    while (std::getline(f, line)) {
      int c;
      double lo, n;
      if (sscanf(line.c_str(), "%d,%lf,%lf", &c, &lo, &n) != 3) continue;
      if (c == 0 && prev >= 0 && bin == 0) bin = lo - prev;
      if (c == 0) prev = lo;
      sp[c].push_back(n);
    }
  }
  std::vector<double> G(32, 0);
  {
    std::ifstream in(c1Json);
    auto j = nlohmann::json::parse(in);
    std::vector<double> ok;
    for (const auto& e : j["channels"])
      if (e["ch"].get<int>() < 32 && e.contains("G_mV") && e["G_mV"].is_number()) {
        G[e["ch"].get<int>()] = e["G_mV"].get<double>();
        ok.push_back(e["G_mV"].get<double>());
      }
    std::sort(ok.begin(), ok.end());
    for (auto& g : G)
      if (g <= 0) g = ok[ok.size() / 2];
  }
  std::vector<double> fs;
  {
    std::stringstream ss(factors);
    std::string t;
    while (std::getline(ss, t, ',')) fs.push_back(std::stod(t));
  }
  TRandom3 rng(20261031);
  FILE* fc = fopen((std::string(prefix) + "_perch.csv").c_str(), "w");
  fprintf(fc, "factor,events,ch,accepted_fraction,G_median_mV,G_sd_mV\n");
  printf("factor,events,replica,n_channels,n_accepted,G_median_mV,n_insufficient,n_rejected,n_ambiguous\n");
  for (double f : fs) {
    std::map<int, std::vector<double>> gacc;
    std::map<int, int> nacc;
    for (int r = 0; r < replicas; ++r) {
      int na = 0, nc = 0, ni = 0, nr = 0, nm = 0;
      std::vector<double> gs;
      for (auto& [c, h] : sp) {
        if (c == 3 || c >= 32) continue;  // ch3 muerto
        ++nc;
        std::vector<double> s(h.size());
        for (size_t b = 0; b + 1 < h.size(); ++b) s[b] = (f == 1.0 && r == 0) ? h[b] : rng.Poisson(f * h[b]);
        s.back() = 0;  // desborde fuera
        const auto res = daq::gain::fitGain(s, bin, G[c] * seedScale, o);
        if (res.status == "insufficient") ++ni;
        if (res.status == "rejected") ++nr;
        if (res.status == "ambiguous") ++nm;
        if (res.accepted()) {
          ++na;
          gs.push_back(res.gainMv);
          gacc[c].push_back(res.gainMv);
          ++nacc[c];
        }
      }
      std::sort(gs.begin(), gs.end());
      printf("%.2f,%ld,%d,%d,%d,%.4f,%d,%d,%d\n", f, static_cast<long>(f * n0), r, nc, na, gs.empty() ? 0.0 : gs[gs.size() / 2], ni, nr, nm);
    }
    for (auto& [c, h] : sp) {
      if (c == 3 || c >= 32) continue;
      auto v = gacc[c];
      double m = 0, sd = 0;
      if (!v.empty()) {
        std::sort(v.begin(), v.end());
        m = v[v.size() / 2];
        double mu = 0;
        for (double x : v) mu += x;
        mu /= v.size();
        for (double x : v) sd += (x - mu) * (x - mu);
        sd = v.size() > 1 ? std::sqrt(sd / (v.size() - 1)) : 0;
      }
      fprintf(fc, "%.2f,%ld,%d,%.3f,%.4f,%.4f\n", f, static_cast<long>(f * n0), c, static_cast<double>(nacc[c]) / replicas, m, sd);
    }
  }
  fclose(fc);
}
