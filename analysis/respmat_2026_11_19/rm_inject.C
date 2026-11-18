// rm_inject.C — seguimiento 2026-11-19 (matriz de respuesta), segunda pasada: inyección de n × plantilla de 1 p.e. en segmentos reales de ruido de hod2018
// (crudos solo en lectura, READ) y estimador exacto de los espectros (hitdef_extract.C: pedestal = media recortada ±2 mV de la mediana de [0,40); pulse_finder.h,
// umbral 0.5 mV, separación 4, región signal [49,58); aS = mayor amplitud, 0 si no hay pulso) con la calibración de los espectros (A = aS / G).
// Para cada evento y canal elegibles (máscara de rm_prep.py, mismo orden de eventos que rm_extract.C) y cada n = 0 … 60 (una inyección por n, con sorteos
// independientes): muestra entera j de la distribución neta de floor(t_CFD) de los muones etiquetados de la capa, fase u uniforme en [0,1), t = j + u;
// x_i −= round(n · T_c(i − t) / lsb) · lsb (la señal inyectada se cuantiza en cuentas del ADC como los datos; T_c = plantilla de rm_prep.py en mV, τ = −1.5 … 40,
// interpolación lineal en la rejilla de 1/8 de muestra, 0 fuera). Sin diafonía ni dispersión temporal del centelleo [supuesto].
// El pedestal se recalcula solo si la inyección toca [0,40) (con τ ≥ −1.5 y t ≥ 40 no la toca: el valor es idéntico).
// Salida <prefijo>_resp.bin (float64, en este orden):
//   H[32][61][241]   cuentas de A en bins de 0.25 p.e. medidos de 0 a 60 + desborde (A ≥ 60);
//   F[32][4][600]    n = 0 … 3, bins de 0.01 p.e. de 0 a 6;
//   S[32][61][3]     número de inyecciones, Σ A, Σ A²;
//   Hc[32][241], Fc[32][600]  n = 0 con la ventana de control [99,108) (aC / G, mismo pedestal).
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/respmat_2026_11_19/rm_inject.C+("<archivo(s)>", "may"|"daqv2", "<G 32 valores>", "<work>/<camp>", <semilla>)'
//   lee <work>/<camp>_tmpl.bin, _tdist.bin, _mask.bin; argumentos opcionales: plantilla alternativa y archivo de salida
//   (variante "gnorm" de rm_gnorm.py: rm_inject.C+(..., <semilla>, "<work>/<camp>_gnorm_tmpl.bin", "<work>/<camp>_gnorm_resp.bin")).
//   Noveno argumento opcional: distribución de tiempos alternativa (variante "tcut" de rm_tcut.py, muestras enteras 48-55).
#include <TChain.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
constexpr int kSpan = 144, kN = 61, kNb = 241, kNf = 600, kNT = 333;
constexpr double kTau0 = -1.5, kDt = 0.125;
template <class T>
std::vector<T> readAll(const std::string& path) {
  std::vector<T> v;
  FILE* f = fopen(path.c_str(), "rb");
  if (!f) return v;
  fseek(f, 0, SEEK_END);
  const long sz = ftell(f);
  fseek(f, 0, SEEK_SET);
  v.resize(static_cast<std::size_t>(sz) / sizeof(T));
  if (fread(v.data(), sizeof(T), v.size(), f) != v.size()) v.clear();
  fclose(f);
  return v;
}
}  // namespace

void rm_inject(const char* pattern, const char* format, const char* gains, const char* pref, long seed, const char* tmplPath = "", const char* outPath = "",
               const char* tdistPath = "") {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  }
  const auto tmpl = readAll<double>(std::string(tmplPath).empty() ? std::string(pref) + "_tmpl.bin" : std::string(tmplPath));
  const auto tdist = readAll<double>(std::string(tdistPath).empty() ? std::string(pref) + "_tdist.bin" : std::string(tdistPath));
  const auto mask = readAll<uint8_t>(std::string(pref) + "_mask.bin");
  if (tmpl.size() != 32u * kNT || tdist.size() != 32u * kSpan || mask.empty()) {
    printf("entradas inválidas: tmpl %zu, tdist %zu, mask %zu\n", tmpl.size(), tdist.size(), mask.size());
    return;
  }
  std::vector<std::discrete_distribution<int>> dj;
  for (int c = 0; c < 32; ++c) {
    std::vector<double> w(tdist.begin() + c * kSpan, tdist.begin() + (c + 1) * kSpan);
    double s = 0;
    for (double v : w) s += v;
    if (s <= 0) w.assign(kSpan, 0.0), w[50] = 1.0;
    dj.emplace_back(w.begin(), w.end());
  }
  std::mt19937_64 rng(static_cast<uint64_t>(seed));
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  const bool may = std::string(format) == "may";
  const int off = may ? 258 - 49 : 0;
  const double lsb = daq::pulse::lsbMv();
  TChain t(may ? "data_dig" : "events");
  t.Add(pattern);
  t.SetBranchStatus("*", 0);
  std::vector<float>* cw[32] = {};
  std::vector<UShort_t> adc;
  int n = 0;
  if (may) {
    for (int c = 0; c < 32; ++c) {
      t.SetBranchStatus(Form("Ch%d", c), 1);
      t.SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
  } else {
    t.GetEntry(0);
    n = t.GetLeaf("adc")->GetLenStatic() / 64;
    adc.resize(static_cast<std::size_t>(64 * n));
    t.SetBranchStatus("adc", 1);
    t.SetBranchAddress("adc", adc.data());
  }
  t.SetCacheSize(256LL << 20);
  std::vector<double> H(32u * kN * kNb, 0.0), F(32u * 4 * kNf, 0.0), S(32u * kN * 3, 0.0), Hc(32u * kNb, 0.0), Fc(32u * kNf, 0.0);
  const std::vector<daq::pulse::Region> sig{{49, 58}}, ctl{{99, 108}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x(kSpan), xi(kSpan);
  const std::size_t nev = mask.size() / 32;
  const Long64_t ne = t.GetEntries();
  std::size_t iw = 0;
  long ninj = 0, nblre = 0;
  auto amax = [&](const std::vector<double>& v, double bl, const std::vector<daq::pulse::Region>& r) {
    p.clear();
    daq::pulse::findPulses(v.data(), kSpan, bl, 1.0, 0.5, 4, r, p);
    double a = 0;
    for (const auto& u : p) a = std::max(a, u.amplitudeMv);
    return a;
  };
  auto fill = [&](double A, double* h, double* fh) {
    h[A >= 60.0 ? kNb - 1 : static_cast<int>(A / 0.25)] += 1;
    if (fh && A < 6.0) fh[static_cast<int>(A / 0.01)] += 1;
  };
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    bool okEv = true;
    if (may) {
      for (int c = 0; c < 32; ++c)
        if (!cw[c] || static_cast<int>(cw[c]->size()) < off + kSpan) okEv = false;
    } else if (n < kSpan) {
      okEv = false;
    }
    if (!okEv) continue;  // como rm_extract.C: el evento no se escribió
    if (iw >= nev) {
      printf("más eventos que la máscara\n");
      break;
    }
    for (int c = 0; c < 32; ++c) {
      if (!mask[iw * 32 + static_cast<std::size_t>(c)]) continue;
      for (int i = 0; i < kSpan; ++i)
        x[static_cast<std::size_t>(i)] = may ? (*cw[c])[static_cast<std::size_t>(off + i)] : adc[static_cast<std::size_t>(c * n + i)] * lsb;
      const double bl0 = daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0);
      const double g = G[static_cast<std::size_t>(c)];
      const double* T = tmpl.data() + static_cast<std::size_t>(c) * kNT;
      {  // ventana de control, n = 0
        const double A = amax(x, bl0, ctl) / g;
        fill(A, Hc.data() + static_cast<std::size_t>(c) * kNb, Fc.data() + static_cast<std::size_t>(c) * kNf);
      }
      for (int k = 0; k < kN; ++k) {
        xi = x;
        double bl = bl0;
        if (k > 0) {
          const double tt = dj[static_cast<std::size_t>(c)](rng) + uni(rng);
          const int i0 = std::max(0, static_cast<int>(std::ceil(tt + kTau0))), i1 = std::min(kSpan - 1, static_cast<int>(std::floor(tt + kTau0 + kDt * (kNT - 1))));
          for (int i = i0; i <= i1; ++i) {
            const double q = (i - tt - kTau0) / kDt;
            const int m = std::min(kNT - 2, static_cast<int>(q));
            const double fr = q - m;
            const double v = k * (T[m] * (1 - fr) + T[m + 1] * fr);
            xi[static_cast<std::size_t>(i)] -= std::round(v / lsb) * lsb;
          }
          if (i0 < 40) bl = daq::pulse::trimmedMean(xi.data(), daq::pulse::Region{0, 40}, 2.0), ++nblre;
          ++ninj;
        }
        const double A = amax(xi, bl, sig) / g;
        const std::size_t cn = static_cast<std::size_t>(c) * kN + static_cast<std::size_t>(k);
        fill(A, H.data() + cn * kNb, k < 4 ? F.data() + (static_cast<std::size_t>(c) * 4 + static_cast<std::size_t>(k)) * kNf : nullptr);
        S[cn * 3] += 1, S[cn * 3 + 1] += A, S[cn * 3 + 2] += A * A;
      }
    }
    ++iw;
    if (e % 20000 == 0) printf("evento %lld de %lld\n", e, ne), fflush(stdout);
  }
  FILE* fo = fopen((std::string(outPath).empty() ? std::string(pref) + "_resp.bin" : std::string(outPath)).c_str(), "wb");
  for (const auto* v : {&H, &F, &S, &Hc, &Fc}) fwrite(v->data(), sizeof(double), v->size(), fo);
  fclose(fo);
  printf("eventos %lld, usados %zu (máscara %zu), inyecciones con n > 0: %ld, pedestal recalculado: %ld\n", ne, iw, nev, ninj, nblre);
}
