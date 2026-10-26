// dark_stage2.C — XPLUS-4b etapa 2 (sobre las ntuplas de ntuple_extract.C): pulsos oscuros por canal.
//  - Escalera: pulsos sobre T·G (hod2018, T en p.e.; G por canal como argumento) o sobre T mV (hod2019) en la región elegida.
//    Exposición por canal = Σ_eventos (longitud de la región − 3 muestras × pulsos de ese canal en la región) × 16 ns: el tiempo
//    muerto de 3 muestras de cada pulso hallado (piso 1 mV) se descuenta, así la tasa es directa.
//  - Intervalos: Δ (muestras) entre pulsos consecutivos > 0.5 p.e. del mismo canal dentro de la misma región.
//  - Pulsos > 6 p.e.: canales con pulso en la misma muestra ±1 (otros canales de la misma coordenada de hod2018 > 1 p.e., de la
//    coordenada ortogonal > 1 p.e., de hod2019 > 3 mV) y Δ al pulso anterior del mismo canal (> 0.5 p.e.) en la región.
// Regiones: "pre" = [20, s0) (física: antes de la ventana de señal), "post" = [s0+9, n), "all" = las dos (SW trigger: todo [20, n)).
// X de hod2018 = ch 0-15, Y = ch 16-31 (stage2.py), hod2019 = ch 32-63. ch3 se excluye de hod2018 (sin señal).
// Salida: <prefijo>_stair.csv, <prefijo>_dt.csv, <prefijo>_big.csv, <prefijo>_expo.csv.
// Eventos: todos, o sin los intervalos UTC de un archivo de excluidos (tiempo bueno de mayo, de export_inputs.py).
// Uso: root -l -b -q 'analysis/night_2026-10-25/dark_stage2.C+("<ntupla.root>", "<G de ch0-31, coma>", s0, "pre|post|all", "<prefijo>"[, "<excluidos.txt>"])'
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <utility>
#include <cstdio>
#include <map>
#include <sstream>
#include <string>
#include <vector>

void dark_stage2(const char* ntuple, const char* gains, int s0, const char* region, const char* prefix, const char* excluded = "") {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  }
  const std::string reg(region);
  std::vector<std::pair<double, double>> exI;  // intervalos UTC excluidos (tiempo bueno), "t0 t1" por línea
  if (excluded && *excluded) {
    std::ifstream in(excluded);
    double a, b;
    while (in >> a >> b) exI.emplace_back(a, b);
  }
  const double kTs = 16e-9;
  const std::vector<double> Tpe{0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
  const std::vector<double> Tmv{1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0, 12.0, 20.0};
  TFile f(ntuple);
  auto* t = static_cast<TTree*>(f.Get("ev"));
  std::vector<UChar_t>* pch = nullptr;
  std::vector<UShort_t>* psmp = nullptr;
  std::vector<float>* pamp = nullptr;
  int reglen = 0;
  double tU = 0;
  t->SetBranchAddress("t_utc_s", &tU);
  t->SetBranchAddress("pch", &pch);
  t->SetBranchAddress("psmp", &psmp);
  t->SetBranchAddress("pamp", &pamp);
  t->SetBranchAddress("reglen", &reglen);
  std::vector<std::vector<long>> cnt(64, std::vector<long>(9, 0));
  std::vector<double> expo(64, 0.0);
  std::vector<std::map<int, long>> dt(64);
  FILE* fbig = fopen((std::string(prefix) + "_big.csv").c_str(), "w");
  fprintf(fbig, "event,ch,sample,amp_mV,amp_pe,n_same_coord,n_orth_coord,n_h19,dt_prev_samples\n");
  const Long64_t ne = t->GetEntries();
  long nEv = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t->GetEntry(e);
    bool skip = false;
    for (auto& [a, b] : exI)
      if (tU >= a && tU < b) skip = true;
    if (skip) continue;
    ++nEv;
    auto inReg = [&](int k) {
      if (reg == "all") return k >= 20;
      if (reg == "pre") return k >= 20 && k < s0;
      return k >= s0 + 9;
    };
    const int lenPre = s0 - 20, lenPost = reglen - (s0 + 9);
    const int len = reg == "all" ? (s0 > 0 ? lenPre + std::max(0, lenPost) : reglen - 20) : (reg == "pre" ? lenPre : lenPost);
    std::vector<int> nInReg(64, 0);
    std::vector<int> lastK(64, -1000000);
    std::vector<int> lastRegPre(64, -1);
    // pulsos por canal en orden de muestra (la ntupla ya está ordenada por canal y muestra dentro de cada región)
    const std::size_t np = pch->size();
    for (std::size_t i = 0; i < np; ++i) {
      const int c = (*pch)[i], k = (*psmp)[i];
      if (!inReg(k)) continue;
      ++nInReg[static_cast<std::size_t>(c)];
      const double a = (*pamp)[i];
      const bool h18 = c < 32;
      const double g = h18 && c < static_cast<int>(G.size()) ? G[static_cast<std::size_t>(c)] : 0;
      for (int j = 0; j < 9; ++j) {
        const double thr = h18 ? Tpe[static_cast<std::size_t>(j)] * g : Tmv[static_cast<std::size_t>(j)];
        if ((h18 ? g > 0 : true) && a > thr) ++cnt[static_cast<std::size_t>(c)][static_cast<std::size_t>(j)];
      }
      const bool above05 = h18 ? (g > 0 && a > 0.5 * g) : a > 1.0;
      const int regId = k < s0 ? 0 : 1;
      if (above05) {
        if (lastRegPre[static_cast<std::size_t>(c)] == regId && lastK[static_cast<std::size_t>(c)] > -1000000)
          ++dt[static_cast<std::size_t>(c)][k - lastK[static_cast<std::size_t>(c)]];
        if (h18 && g > 0 && a > 6.0 * g && c != 3) {
          int nSame = 0, nOrth = 0, n19 = 0;
          for (std::size_t q = 0; q < np; ++q) {
            const int c2 = (*pch)[q], k2 = (*psmp)[q];
            if (c2 == c || std::abs(k2 - k) > 1) continue;
            const double a2 = (*pamp)[q];
            if (c2 < 32) {
              if (c2 == 3) continue;
              const double g2 = c2 < static_cast<int>(G.size()) ? G[static_cast<std::size_t>(c2)] : 0;
              if (g2 <= 0 || a2 <= g2) continue;
              ((c2 < 16) == (c < 16) ? nSame : nOrth)++;
            } else if (a2 > 3.0) {
              ++n19;
            }
          }
          const int dprev = lastRegPre[static_cast<std::size_t>(c)] == regId ? k - lastK[static_cast<std::size_t>(c)] : -1;
          fprintf(fbig, "%lld,%d,%d,%.3f,%.3f,%d,%d,%d,%d\n", e, c, k, a, a / g, nSame, nOrth, n19, dprev);
        }
        lastK[static_cast<std::size_t>(c)] = k;
        lastRegPre[static_cast<std::size_t>(c)] = regId;
      }
    }
    for (int c = 0; c < 64; ++c) expo[static_cast<std::size_t>(c)] += (len - 3.0 * nInReg[static_cast<std::size_t>(c)]) * kTs;
  }
  fclose(fbig);
  FILE* fs = fopen((std::string(prefix) + "_stair.csv").c_str(), "w");
  fprintf(fs, "ch,threshold,unit,counts,exposure_s,rate_cps,rate_err_cps\n");
  for (int c = 0; c < 64; ++c)
    for (int j = 0; j < 9; ++j) {
      const bool h18 = c < 32;
      const long k = cnt[static_cast<std::size_t>(c)][static_cast<std::size_t>(j)];
      const double ex = expo[static_cast<std::size_t>(c)];
      fprintf(fs, "%d,%.2f,%s,%ld,%.9f,%.6g,%.6g\n", c, h18 ? Tpe[static_cast<std::size_t>(j)] : Tmv[static_cast<std::size_t>(j)], h18 ? "pe" : "mV", k, ex,
              ex > 0 ? k / ex : 0.0, ex > 0 ? std::sqrt(static_cast<double>(k)) / ex : 0.0);
    }
  fclose(fs);
  FILE* fd = fopen((std::string(prefix) + "_dt.csv").c_str(), "w");
  fprintf(fd, "ch,dt_samples,count\n");
  for (int c = 0; c < 64; ++c)
    for (auto& [d, n] : dt[static_cast<std::size_t>(c)]) fprintf(fd, "%d,%d,%ld\n", c, d, n);
  fclose(fd);
  FILE* fe = fopen((std::string(prefix) + "_expo.csv").c_str(), "w");
  fprintf(fe, "events,region,s0\n%ld,%s,%d\n", nEv, region, s0);
  fclose(fe);
  printf("eventos %ld región %s\n", nEv, region);
}
