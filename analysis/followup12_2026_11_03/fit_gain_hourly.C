// fit_gain_hourly.C — seguimiento 2026-11-03 (12), punto B: G por canal de hod2018 y hora UTC desde el espectro oscuro.
// Fuentes del espectro (bins de 0.25 mV de 0 a 100 mV, bin 400 = desborde, no se usa):
//   "spec": árbol "spec" de la ntupla XPLUS-4b (método del reajuste de G de hitdef_extract.C / PROC: tramo [0,49) ∪ [102,144), umbral 0.9 mV,
//           separación 4; 91 muestras por evento) — fuente central, la misma definición de G que PROC;
//   "dark": <prefijo>_dark.bin de ev_dump.C (pulsos de la región pre-trigger, buscador del DAQ con separación 3 y piso 1 mV) — variante.
// Ajuste: el de fit_gain_long.C (seguimiento 6): bins de 0.5 mV, suma de gaussianas en μ_n = p0 + n·G con σ_n² = s0² + n·s1², amplitudes libres
// y fondo exponencial, en dos pasos (amplitudes y fondo con G, p0, s0, s1 fijos; luego todo libre), en la ventana [1.5, 5.5]·G_semilla (picos de 2 a 5
// p.e.; el bulto de 1-4 mV tapa el de 1 p.e.). Primero el espectro del run entero (horas elegidas) con la semilla de PROC (G_fit del _peaks.csv);
// luego cada hora con los parámetros del run como punto de partida; central: p0 fijo al valor del run del canal (con solo los picos de 2-5 p.e.,
// (ajuste del run aceptado con estado de MIGRAD 0 o 1 y error > 0; si no, p0 libre también en las horas)
// p0 y G están correlacionados y p0 libre añade ruido a G); variante fixP0 = false: p0 libre en ±0.5·G_semilla.
// Salida CSV: hour_utc_s,ch,events,G_mV,G_err_mV,p0_mV,s0_mV,s1_mV,chi2,ndf,fit_status,counts_window   (hour_utc_s = −1: run entero)
// Uso: root -l -b -q 'analysis/followup12_2026_11_03/fit_gain_hourly.C+("spec"|"dark", "<ntupla.root o _dark.bin>", "<semillas 32 valores>", "<salida.csv>"[, fixP0 = true])'
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TTree.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Par {
  double G, p0, s0, s1;
};
struct Out {
  double G = NAN, Ge = NAN, p0 = NAN, s0 = NAN, s1 = NAN, chi2 = NAN;
  int ndf = 0, st = -1;
  double cnt = 0;
};
Out fitOne(const std::vector<double>& raw, const Par& start, double Gseed, int tag, bool fixP0 = false) {
  std::vector<double> spec;
  for (size_t b = 0; b + 1 < 400; b += 2) spec.push_back(raw[b] + raw[b + 1]);
  const double bin = 0.5, lo = 1.5 * Gseed, hi = 5.5 * Gseed;
  const int n0 = 2, n1 = 5, np = n1 - n0 + 1;
  TH1D h(Form("h_%d", tag), "", static_cast<int>(spec.size()), 0, spec.size() * bin);
  for (size_t b = 0; b < spec.size(); ++b) {
    h.SetBinContent(b + 1, spec[b]);
    h.SetBinError(b + 1, std::sqrt(std::max(1.0, spec[b])));
  }
  Out o;
  o.cnt = h.Integral(h.FindBin(lo), h.FindBin(hi));
  if (o.cnt < 50) return o;
  auto fn = [np](double* x, double* p) {
    double v = std::exp(p[4] + p[5] * x[0]);
    for (int k = 0; k < np; ++k) {
      const int n = n0 + k;
      const double mu = p[0] + n * p[1];
      const double s2 = p[2] * p[2] + n * p[3] * p[3];
      v += p[6 + k] * std::exp(-0.5 * (x[0] - mu) * (x[0] - mu) / s2);
    }
    return v;
  };
  TF1 f(Form("f_%d", tag), fn, lo, hi, 6 + np);
  const double G = start.G;
  f.SetParameter(0, start.p0);
  f.SetParameter(1, G);
  f.SetParameter(2, start.s0);
  f.SetParameter(3, start.s1);
  const int blo = h.FindBin(lo), bhi = h.FindBin(hi);
  const double ylo = std::max(1.0, h.GetBinContent(blo)), yhi = std::max(1.0, h.GetBinContent(bhi));
  const double slope = std::log(yhi / ylo) / (hi - lo);
  f.SetParameter(4, std::log(ylo) - slope * lo);
  f.SetParameter(5, slope);
  for (int k = 0; k < np; ++k) {
    const double mu = start.p0 + (n0 + k) * G;
    const double y = std::max(1.0, h.GetBinContent(h.FindBin(mu)));
    const double bk = std::exp(f.GetParameter(4) + slope * mu);
    f.SetParameter(6 + k, std::max(0.05 * y, y - bk));
    f.SetParLimits(6 + k, 0.0, 10 * std::max(10.0, y));
  }
  for (int k = 0; k < 4; ++k) f.FixParameter(k, f.GetParameter(k));
  h.Fit(&f, "QRN0");
  if (!fixP0) {
    f.ReleaseParameter(0);
    f.SetParLimits(0, -0.5 * Gseed, 0.5 * Gseed);
  }
  f.ReleaseParameter(1);
  f.SetParLimits(1, 0.6 * Gseed, 1.4 * Gseed);
  f.ReleaseParameter(2);
  f.SetParLimits(2, 0.02, 2 * Gseed);
  f.ReleaseParameter(3);
  f.SetParLimits(3, 0.0, Gseed);
  auto r = h.Fit(&f, "QSRN0");
  o.G = f.GetParameter(1);
  o.Ge = f.GetParError(1);
  o.p0 = f.GetParameter(0);
  o.s0 = f.GetParameter(2);
  o.s1 = f.GetParameter(3);
  o.chi2 = r->Chi2();
  o.ndf = r->Ndf();
  o.st = static_cast<int>(r);
  return o;
}
}  // namespace

void fit_gain_hourly(const char* source, const char* input, const char* seeds, const char* outCsv, bool fixP0 = true) {
  std::vector<double> seed;
  {
    std::stringstream ss(seeds);
    std::string tok;
    while (std::getline(ss, tok, ',')) seed.push_back(std::stod(tok));
  }
  // espectros por (hora, canal)
  std::map<std::pair<int64_t, int>, std::pair<long, std::vector<double>>> sp;
  if (std::string(source) == "spec") {
    TFile f(input, "READ");
    auto* t = static_cast<TTree*>(f.Get("spec"));
    Long64_t hr = 0, nev = 0;
    int ch = 0;
    unsigned counts[401];
    t->SetBranchAddress("hour_utc_s", &hr);
    t->SetBranchAddress("ch", &ch);
    t->SetBranchAddress("events", &nev);
    t->SetBranchAddress("counts", counts);
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      if (ch >= 32) continue;
      auto& a = sp[{hr, ch}];
      a.first = nev;
      a.second.assign(counts, counts + 401);
    }
  } else {
    FILE* fd = fopen(input, "rb");
    int64_t hr, ev;
    int32_t ch, pad;
    double expo;
    uint32_t c[401];
    while (fread(&hr, 8, 1, fd) == 1) {
      if (fread(&ch, 4, 1, fd) != 1 || fread(&pad, 4, 1, fd) != 1 || fread(&ev, 8, 1, fd) != 1 || fread(&expo, 8, 1, fd) != 1 ||
          fread(c, 4, 401, fd) != 401)
        break;
      if (ch >= 32) continue;
      auto& a = sp[{hr, ch}];
      a.first = ev;
      a.second.assign(c, c + 401);
    }
    fclose(fd);
  }
  FILE* fo = fopen(outCsv, "w");
  fprintf(fo, "# fit_gain_hourly.C fuente %s: %s; p0 por hora %s\n", source, input, fixP0 ? "fijo al del run" : "libre");
  fprintf(fo, "hour_utc_s,ch,events,G_mV,G_err_mV,p0_mV,s0_mV,s1_mV,chi2,ndf,fit_status,counts_window\n");
  int tag = 0;
  for (int c = 0; c < 32; ++c) {
    if (c == 3) continue;  // sin señal
    const double gs = c < static_cast<int>(seed.size()) ? seed[c] : 4.0;
    std::vector<double> tot(401, 0.0);
    long evt = 0;
    for (auto& [k, a] : sp)
      if (k.second == c) {
        evt += a.first;
        for (int b = 0; b < 401; ++b) tot[b] += a.second[b];
      }
    const Out R = fitOne(tot, {gs, 0.0, 0.3 * gs, 0.1 * gs}, gs, tag++);
    fprintf(fo, "-1,%d,%ld,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%d,%d,%.0f\n", c, evt, R.G, R.Ge, R.p0, R.s0, R.s1, R.chi2, R.ndf, R.st, R.cnt);
    const bool runOk = (R.st == 0 || R.st == 1) && R.Ge > 0;
    const Par st = runOk ? Par{R.G, R.p0, R.s0, R.s1} : Par{gs, 0.0, 0.3 * gs, 0.1 * gs};
    for (auto& [k, a] : sp) {
      if (k.second != c) continue;
      const Out o = fitOne(a.second, st, gs, tag++, fixP0 && runOk);
      fprintf(fo, "%lld,%d,%ld,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%d,%d,%.0f\n", static_cast<long long>(k.first), c, a.first, o.G, o.Ge, o.p0, o.s0, o.s1,
              o.chi2, o.ndf, o.st, o.cnt);
    }
    printf("ch%d: G_run %.4f ± %.4f (semilla %.4f) χ² %.1f/%d\n", c, R.G, R.Ge, gs, R.chi2, R.ndf);
  }
  fclose(fo);
}
