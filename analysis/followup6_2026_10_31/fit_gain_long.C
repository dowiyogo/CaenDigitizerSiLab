// fit_gain_long.C — seguimiento 2026-10-31 (6), punto 3 (solo lectura): G por canal desde el espectro oscuro de dark_spectrum_long.C con el
// método de c1/PROC (daq/core/include/daq/gain_fit.h: separación de picos sobre el residuo del espectro suavizado), más un ajuste de ROOT
// de los picos para dar χ²/ndf.
// Configuraciones de gain_fit.h:
//   "proc":  opciones por defecto (picos entre 5 y 26 mV, mínimo 3 picos, tolerancia de espaciado 0.2), bins de 0.25 mV (como hitdef_extract.C);
//   "low":   variante para una G pequeña (picos entre 1.5 y 12 mV), bins reagrupados a 0.5 mV para no ver la cuantización del ADC (LSB 0.488 mV,
//            que en bins de 0.25 mV da una oscilación par/impar);
// cada una con una rejilla de semillas (1.0 … 6.0 mV cada 0.5) y, si se da, la semilla propia del canal (hod2018: reajuste de PROC).
// Ajuste de ROOT (solo si gain_fit.h acepta con la semilla propia o, sin ella, con la mediana de las G aceptadas de la rejilla): suma de
// gaussianas en μ_n = p0 + n·G (n = los índices de p.e. de los picos que caen en la ventana), σ_n² = s0² + n·s1², amplitudes libres, fondo
// exponencial, en bins de 0.5 mV (en dos pasos: amplitudes y fondo con G, s0, s1 fijos, luego todo libre); χ²/ndf con errores de
// Poisson en la ventana de la configuración.
// Salida: <prefijo>_fits.csv (ch, config, seed_mV, status, peaks, G_mV, positions) y <prefijo>_rootfit.csv (ch, config, G_mV, G_err, p0_mV,
//         s0, s1, chi2, ndf, status).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'fit_gain_long.C+("<spec.csv>", ch_lo, ch_hi, "<semillas propias: 64 valores con coma o vacío>", "<prefijo>")'
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "daq/gain_fit.h"

namespace {
struct Cfg {
  const char* name;
  int rebin;
  double lo, hi;
};

void rootFit(const std::vector<double>& raw, const Cfg& cfg, double G, int ch, FILE* fo) {
  // siempre en bins de 0.5 mV (sin la oscilación par/impar de la cuantización); picos con μ = n·G dentro de la ventana
  std::vector<double> spec;
  for (size_t b = 0; b + 2 <= raw.size() - 1; b += 2) spec.push_back(raw[b] + raw[b + 1]);
  const double bin = 0.5;
  const int n0 = std::max(1, static_cast<int>(std::ceil(cfg.lo / G)));
  const int n1 = static_cast<int>(std::floor(cfg.hi / G));
  const int np = n1 - n0 + 1;
  if (np < 2) return;
  TH1D h(Form("h%d_%s", ch, cfg.name), "", static_cast<int>(spec.size()), 0, spec.size() * bin);
  for (size_t b = 0; b < spec.size(); ++b) {
    h.SetBinContent(b + 1, spec[b]);
    h.SetBinError(b + 1, std::sqrt(std::max(1.0, spec[b])));
  }
  // parámetros: 0 p0, 1 G, 2 s0, 3 s1, 4 fondo A, 5 fondo pendiente, 6.. amplitudes
  auto fn = [n0, np](double* x, double* p) {
    double v = std::exp(p[4] + p[5] * x[0]);
    for (int k = 0; k < np; ++k) {
      const int n = n0 + k;
      const double mu = p[0] + n * p[1];
      const double s2 = p[2] * p[2] + n * p[3] * p[3];
      v += p[6 + k] * std::exp(-0.5 * (x[0] - mu) * (x[0] - mu) / s2);
    }
    return v;
  };
  TF1 f(Form("f%d_%s", ch, cfg.name), fn, cfg.lo, cfg.hi, 6 + np);
  f.SetParameter(0, 0.0);
  f.SetParameter(1, G);
  f.SetParLimits(1, 0.5 * G, 1.5 * G);
  f.SetParameter(2, 0.3 * G);
  f.SetParLimits(2, 0.02, 2 * G);
  f.SetParameter(3, 0.1 * G);
  f.SetParLimits(3, 0.0, G);
  const int blo = h.FindBin(cfg.lo), bhi = h.FindBin(cfg.hi);
  const double ylo = std::max(1.0, h.GetBinContent(blo)), yhi = std::max(1.0, h.GetBinContent(bhi));
  const double slope = std::log(yhi / ylo) / (cfg.hi - cfg.lo);
  f.SetParameter(4, std::log(ylo) - slope * cfg.lo);
  f.SetParameter(5, slope);
  for (int k = 0; k < np; ++k) {
    const double mu = (n0 + k) * G;
    const double y = std::max(1.0, h.GetBinContent(h.FindBin(mu)));
    const double bk = std::exp(f.GetParameter(4) + slope * mu);
    f.SetParameter(6 + k, std::max(0.05 * y, y - bk));
    f.SetParLimits(6 + k, 0.0, 10 * std::max(10.0, y));
  }
  // dos pasos: primero amplitudes y fondo con G, s0, s1 fijos; luego todo libre
  for (int k = 0; k < 4; ++k) f.FixParameter(k, f.GetParameter(k));
  h.Fit(&f, "QRN0");
  f.ReleaseParameter(0);
  f.ReleaseParameter(1);
  f.SetParLimits(1, 0.5 * G, 1.5 * G);
  f.ReleaseParameter(2);
  f.SetParLimits(2, 0.02, 2 * G);
  f.ReleaseParameter(3);
  f.SetParLimits(3, 0.0, G);
  auto r = h.Fit(&f, "QSRN0");
  fprintf(fo, "%d,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%d,%d,%d,%d\n", ch, cfg.name, f.GetParameter(1), f.GetParError(1), f.GetParameter(0), f.GetParameter(2),
          f.GetParameter(3), r->Chi2(), r->Ndf(), static_cast<int>(r), n0, n1);
}
}  // namespace

void fit_gain_long(const char* specCsv, int chLo, int chHi, const char* ownSeeds, const char* prefix) {
  std::map<int, std::vector<double>> sp;
  {
    std::ifstream f(specCsv);
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
      int c;
      double lo, n;
      if (sscanf(line.c_str(), "%d,%lf,%lf", &c, &lo, &n) == 3) sp[c].push_back(n);
    }
  }
  std::vector<double> own;
  {
    std::stringstream ss(ownSeeds ? ownSeeds : "");
    std::string tok;
    while (std::getline(ss, tok, ',')) own.push_back(std::stod(tok));
  }
  const Cfg cfgs[] = {{"proc", 1, 5.0, 26.0}, {"low", 2, 1.5, 12.0}};
  FILE* fo = fopen((std::string(prefix) + "_fits.csv").c_str(), "w");
  FILE* fr = fopen((std::string(prefix) + "_rootfit.csv").c_str(), "w");
  fprintf(fo, "ch,config,seed_kind,seed_mV,status,peaks,G_mV,positions_mV\n");
  fprintf(fr, "ch,config,G_mV,G_err,p0_mV,s0_mV,s1_mV,chi2,ndf,fit_status,n_first,n_last\n");
  for (int c = chLo; c <= chHi; ++c) {
    const auto& raw = sp[c];
    for (const auto& cfg : cfgs) {
      std::vector<double> s;
      for (size_t b = 0; b + cfg.rebin <= raw.size() - 1; b += cfg.rebin) {  // sin el bin de desborde
        double v = 0;
        for (int k = 0; k < cfg.rebin; ++k) v += raw[b + k];
        s.push_back(v);
      }
      const double bin = 0.25 * cfg.rebin;
      daq::gain::Options o;
      o.loMv = cfg.lo;
      o.hiMv = cfg.hi;
      std::vector<std::pair<std::string, double>> seeds;
      if (c < static_cast<int>(own.size()) && own[c] > 0) seeds.push_back({"own", own[c]});
      for (double sd = 1.0; sd <= 6.0001; sd += 0.5) seeds.push_back({"grid", sd});
      std::vector<double> acc;
      daq::gain::Result best;
      bool haveOwn = false;
      for (auto& [kind, sd] : seeds) {
        const auto r = daq::gain::fitGain(s, bin, sd, o);
        std::string pos;
        for (double v : r.positionsMv) pos += (pos.empty() ? "" : " ") + std::to_string(v);
        fprintf(fo, "%d,%s,%s,%.2f,%s,%d,%.5f,%s\n", c, cfg.name, kind.c_str(), sd, r.status.c_str(), r.peaks, r.gainMv, pos.c_str());
        if (kind == "own" && r.accepted()) {
          best = r;
          haveOwn = true;
        }
        if (kind == "grid" && r.accepted()) acc.push_back(r.gainMv);
      }
      if (!haveOwn && !acc.empty()) {
        std::sort(acc.begin(), acc.end());
        const double med = acc[acc.size() / 2];
        best = daq::gain::fitGain(s, bin, med, o);
      }
      if (best.accepted()) rootFit(raw, cfg, best.gainMv, c, fr);
    }
  }
  fclose(fo);
  fclose(fr);
}
