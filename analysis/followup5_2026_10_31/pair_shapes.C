// pair_shapes.C — seguimiento 2026-10-31 (5), punto 3 (solo lectura): pares de pulsos del mismo canal separados por pocas muestras en la
// región pre-trigger de mayo v8 (tiempo bueno nuevo), mirados en la forma de onda cruda.
// Pulsos: listas de la ntupla de ntuple_extract.C (buscador del DAQ, separación 3 muestras, piso 1 mV, pedestal local de la región
// [20, s0)), las mismas que usa dark_stage2.C. Par = dos pulsos consecutivos del canal en la región, ambos > 0.5 p.e., Δ = k2 − k1 en
// muestras (3 … 12; el punto pide 3-6). Sobre el crudo (Ch<c>, mV), con P = pedestal local = x[k1] + A1:
//   A1, A2 = amplitudes del buscador (mV sobre P); valle V = P − max(x[k1 … k2]) (la amplitud en el punto más alto entre los dos mínimos);
//   r2 = A2 − V = profundidad del segundo mínimo bajo el valle; T(d) = forma media normalizada de pulsos aislados del canal (0.5-1.5 p.e.,
//   sin otro pulso en ±12 muestras), E = A2 − A1·T(Δ) = amplitud del segundo pulso sobre la cola esperada del primero.
//   Clasificación por valle: "independiente" si r2 ≥ 0.5 G (mínimo local propio, separado del primero por una subida de al menos medio
//   p.e.), "cola" si r2 < 0.5 G. Este criterio solo no basta para A2 < 1 p.e. (el ruido llena el valle: a Δ = 12, donde la cola ya no
//   existe, la mitad de los pares cae en "cola"), así que con un archivo de forma (segunda pasada, la de la primera pasada) se usa también
//   E: "espurio" (cola partida por el buscador) si E < 0.5 G, es decir si el segundo pulso no tiene medio p.e. sobre la cola esperada.
// Escalera del canal en la misma región: N(> 0.5 p.e.), N(> 1.5 p.e.), y las mismas cuentas solo de los segundos pulsos "cola" con Δ
// 3-6 y 3-12, para corregir p_ct = N(>1.5)/N(>0.5).
// Salida: <prefijo>_summary.csv (ch, dmin, dmax, clase, n, n_A2_gt1p5), <prefijo>_hist.csv (ch, clase, Δ, bin A2/G, bin r2/G, n),
//         <prefijo>_template.csv (ch, d, T), <prefijo>_pairs.csv (los primeros npairs pares por canal).
//         con forma: <prefijo>_spurious.csv (ch, Δ, n_pares, n_espurios con A2 > 0.5 y > 1.5 p.e.) y <prefijo>_E_hist.csv (E/G por Δ).
// Formatos: "may" (data_dig, Ch<c> en mV) y "daqv2" (events, adc[64][n] en cuentas; el crudo puede ser un patrón de segmentos, la entrada
// de la ntupla es la del TChain). Región [lo, hi) en muestras del registro: mayo pre = [20, 258); octubre pre = [20, 49), post = [58, 144)
// (cada región de la ntupla tiene su pedestal local; no se mezclan pulsos de dos regiones).
// Uso: root -l -b -q 'pair_shapes.C+("<ntupla>", "<crudo>", "may"|"daqv2", "<excluidos o vacío>", lo, hi, "c1,c2,…", "G1,G2,…", npairs,
//                                       "<prefijo>"[, "<forma.csv>"])'
#include <TChain.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
std::vector<double> parse(const char* s) {
  std::vector<double> v;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ',')) v.push_back(std::stod(t));
  return v;
}
}  // namespace

void pair_shapes(const char* ntuple, const char* raw, const char* format, const char* excl, int lo, int s0, const char* chans, const char* gains, int npairs, const char* prefix, const char* templ = "") {
  const auto chv = parse(chans), gv = parse(gains);
  const int NC = static_cast<int>(chv.size());
  std::vector<std::pair<double, double>> ex;
  if (excl && *excl) {
    std::ifstream f(excl);
    double a, b;
    while (f >> a >> b) ex.push_back({a, b});
  }
  // forma de la primera pasada (ch, n_isolated, d, T), si se da
  std::map<int, std::vector<double>> T0;
  if (templ && *templ) {
    std::ifstream f(templ);
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
      int c, d;
      long n;
      double v;
      if (sscanf(line.c_str(), "%d,%ld,%d,%lf", &c, &n, &d, &v) == 4) {
        auto& tv = T0[c];
        if (static_cast<int>(tv.size()) <= d) tv.resize(d + 1, NAN);
        tv[d] = v;
      }
    }
  }
  const bool haveT = !T0.empty();
  TFile fn(ntuple, "READ");
  TTree* t = (TTree*)fn.Get("ev");
  double tu;
  Long64_t evn;
  std::vector<UChar_t>* pch = nullptr;
  std::vector<UShort_t>* psmp = nullptr;
  std::vector<float>* pamp = nullptr;
  t->SetBranchStatus("*", 0);
  for (const char* b : {"t_utc_s", "ev", "pch", "psmp", "pamp"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("t_utc_s", &tu);
  t->SetBranchAddress("ev", &evn);
  t->SetBranchAddress("pch", &pch);
  t->SetBranchAddress("psmp", &psmp);
  t->SetBranchAddress("pamp", &pamp);
  const bool may = std::string(format) == "may";
  TChain chn(may ? "data_dig" : "events");
  chn.Add(raw);
  TChain* r = &chn;
  r->SetBranchStatus("*", 0);
  std::vector<std::vector<float>*> w(NC, nullptr);
  std::vector<UShort_t> adc;
  int len = 0;
  if (may) {
    for (int i = 0; i < NC; ++i) {
      const char* bn = Form("Ch%d", static_cast<int>(chv[i]));
      r->SetBranchStatus(bn, 1);
      r->SetBranchAddress(bn, &w[i]);
    }
  } else {
    len = r->GetLeaf("adc")->GetLenStatic() / 64;
    adc.resize(64 * len);
    r->SetBranchStatus("adc", 1);
    r->SetBranchAddress("adc", adc.data());
    for (int i = 0; i < NC; ++i) w[i] = new std::vector<float>(len);
  }
  const double lsb = 2000.0 / 4095.0;
  const int DMAX = 12, TD = 13, NBA = 120, NBR = 120;  // A2/G y r2/G en bins de 0.05 p.e. hasta 6 p.e. (r2 desde −1)
  std::vector<std::vector<double>> tsum(NC, std::vector<double>(TD, 0));
  std::vector<long> tn(NC, 0), n05(NC, 0), n15(NC, 0);
  // [ch][clase 0 = independiente, 1 = cola][Δ][binA][binR]
  std::vector<std::vector<std::vector<std::vector<long>>>> hA(NC, std::vector<std::vector<std::vector<long>>>(2, std::vector<std::vector<long>>(DMAX + 1, std::vector<long>(NBA, 0))));
  std::vector<std::vector<std::vector<long>>> hR(NC, std::vector<std::vector<long>>(DMAX + 1, std::vector<long>(NBR, 0)));
  struct Pair { int d; double a1, a2, v, r2, e; Long64_t ev; int k1; };
  std::vector<std::vector<Pair>> pend(NC);   // pares a los que les falta E (necesita T); se guardan todos los campos crudos
  std::vector<std::vector<std::vector<long>>> cnt(NC, std::vector<std::vector<long>>(DMAX + 1, std::vector<long>(4, 0)));  // n_ind, n_cola, n_ind>1.5, n_cola>1.5
  std::vector<std::vector<std::vector<long>>> cola05(NC, std::vector<std::vector<long>>(DMAX + 1, std::vector<long>(2, 0)));
  const int NBE = 120;  // E/G en bins de 0.05 p.e. desde −1
  std::vector<std::vector<std::vector<long>>> hE(NC, std::vector<std::vector<long>>(DMAX + 1, std::vector<long>(NBE, 0)));
  std::vector<std::vector<std::vector<long>>> spur(NC, std::vector<std::vector<long>>(DMAX + 1, std::vector<long>(3, 0)));  // pares, esp. >0.5, esp. >1.5
  long used = 0;
  const Long64_t ne = t->GetEntries();
  for (Long64_t i = 0; i < ne; ++i) {
    t->GetEntry(i);
    bool good = true;
    for (auto& [a, b] : ex)
      if (tu >= a && tu < b) { good = false; break; }
    if (!good) continue;
    bool loaded = false;
    ++used;
    for (int ic = 0; ic < NC; ++ic) {
      const int c = static_cast<int>(chv[ic]);
      const double G = gv[ic];
      std::vector<std::pair<int, double>> pl;  // pulsos del canal en la región pre, en orden
      for (size_t q = 0; q < pch->size(); ++q)
        if ((*pch)[q] == c && (*psmp)[q] >= lo && (*psmp)[q] < s0) pl.push_back({(*psmp)[q], (*pamp)[q]});
      for (auto& p : pl) {
        if (p.second > 0.5 * G) ++n05[ic];
        if (p.second > 1.5 * G) ++n15[ic];
      }
      // pulsos > 0.5 p.e. y pares consecutivos
      std::vector<std::pair<int, double>> pg;
      for (auto& p : pl)
        if (p.second > 0.5 * G) pg.push_back(p);
      bool needRaw = false;
      for (size_t q = 0; q < pg.size(); ++q) {
        const bool iso = (q == 0 || pg[q].first - pg[q - 1].first > 12) && (q + 1 == pg.size() || pg[q + 1].first - pg[q].first > 12) &&
                         pg[q].second < 1.5 * G && pg[q].first + TD < s0;
        if (iso && tn[ic] < 200000) needRaw = true;
        if (q > 0 && pg[q].first - pg[q - 1].first <= DMAX) needRaw = true;
      }
      if (!needRaw) continue;
      if (!loaded) {
        r->GetEntry(evn);
        if (!may)
          for (int jc = 0; jc < NC; ++jc)
            for (int k = 0; k < len; ++k) (*w[jc])[k] = adc[static_cast<int>(chv[jc]) * len + k] * lsb;
        loaded = true;
      }
      const std::vector<float>& x = *w[ic];
      // también los pulsos de < 0.5 p.e. cuentan para el aislamiento
      for (size_t q = 0; q < pg.size(); ++q) {
        const int k = pg[q].first;
        bool iso = pg[q].second < 1.5 * G && k + TD < s0 && tn[ic] < 200000;
        for (auto& p : pl)
          if (p.first != k && std::abs(p.first - k) <= 12) iso = false;
        if (iso) {
          const double P = x[k] + pg[q].second;
          for (int d = 0; d < TD; ++d) tsum[ic][d] += (P - x[k + d]) / pg[q].second;
          ++tn[ic];
        }
        if (q == 0) continue;
        const int k1 = pg[q - 1].first, k2 = k, d = k2 - k1;
        if (d > DMAX) continue;
        const double a1 = pg[q - 1].second, a2 = pg[q].second, P = x[k1] + a1;
        double mx = x[k1];
        for (int s = k1; s <= k2; ++s) mx = std::max<double>(mx, x[s]);
        const double v = P - mx, r2 = a2 - v;
        const int cl = r2 >= 0.5 * G ? 0 : 1;
        const int bA = std::min(NBA - 1, std::max(0, static_cast<int>(a2 / G / 0.05)));
        const int bR = std::min(NBR - 1, std::max(0, static_cast<int>((r2 / G + 1.0) / 0.05)));
        ++hA[ic][cl][d][bA];
        ++hR[ic][d][bR];
        ++cnt[ic][d][cl];
        if (a2 > 1.5 * G) ++cnt[ic][d][2 + cl];
        if (cl == 1) ++cola05[ic][d][a2 > 1.5 * G ? 1 : 0];
        if (haveT) {
          const double e = a2 - a1 * T0[c][d];
          ++hE[ic][d][std::min(NBE - 1, std::max(0, static_cast<int>((e / G + 1.0) / 0.05)))];
          ++spur[ic][d][0];
          if (e < 0.5 * G) {
            ++spur[ic][d][1];
            if (a2 > 1.5 * G) ++spur[ic][d][2];
          }
        }
        if (static_cast<int>(pend[ic].size()) < npairs && d <= 6) pend[ic].push_back({d, a1, a2, v, r2, 0, evn, k1});
      }
    }
  }
  FILE* ft = fopen((std::string(prefix) + "_template.csv").c_str(), "w");
  fprintf(ft, "ch,n_isolated,d,T\n");
  std::vector<std::vector<double>> T(NC, std::vector<double>(TD, 0));
  for (int ic = 0; ic < NC; ++ic)
    for (int d = 0; d < TD; ++d) {
      T[ic][d] = tn[ic] ? tsum[ic][d] / tn[ic] : NAN;
      fprintf(ft, "%d,%ld,%d,%.5f\n", static_cast<int>(chv[ic]), tn[ic], d, T[ic][d]);
    }
  fclose(ft);
  FILE* fp = fopen((std::string(prefix) + "_pairs.csv").c_str(), "w");
  fprintf(fp, "ch,G_mV,event,k1,delta,A1_mV,A2_mV,valley_mV,r2_mV,E_mV,A2_pe,r2_pe,E_pe,class\n");
  for (int ic = 0; ic < NC; ++ic)
    for (auto& p : pend[ic]) {
      const double G = gv[ic], e = p.a2 - p.a1 * T[ic][p.d];
      fprintf(fp, "%d,%.4f,%lld,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.4f,%s\n", static_cast<int>(chv[ic]), G, p.ev, p.k1, p.d, p.a1, p.a2, p.v, p.r2, e,
              p.a2 / G, p.r2 / G, e / G, p.r2 >= 0.5 * G ? "independiente" : "cola");
    }
  fclose(fp);
  FILE* fs = fopen((std::string(prefix) + "_summary.csv").c_str(), "w");
  fprintf(fs, "ch,G_mV,events,N_gt0p5,N_gt1p5,delta,n_independent,n_tail,n_independent_A2gt1p5,n_tail_A2gt1p5\n");
  for (int ic = 0; ic < NC; ++ic)
    for (int d = 3; d <= DMAX; ++d)
      fprintf(fs, "%d,%.4f,%ld,%ld,%ld,%d,%ld,%ld,%ld,%ld\n", static_cast<int>(chv[ic]), gv[ic], used, n05[ic], n15[ic], d, cnt[ic][d][0], cnt[ic][d][1], cnt[ic][d][2],
              cnt[ic][d][3]);
  fclose(fs);
  FILE* fh = fopen((std::string(prefix) + "_hist.csv").c_str(), "w");
  fprintf(fh, "ch,kind,delta,bin_lo_pe,n\n");
  for (int ic = 0; ic < NC; ++ic)
    for (int d = 3; d <= DMAX; ++d) {
      for (int cl = 0; cl < 2; ++cl)
        for (int b = 0; b < NBA; ++b)
          if (hA[ic][cl][d][b]) fprintf(fh, "%d,%s,%d,%.2f,%ld\n", static_cast<int>(chv[ic]), cl ? "A2_tail" : "A2_independent", d, 0.05 * b, hA[ic][cl][d][b]);
      for (int b = 0; b < NBR; ++b)
        if (hR[ic][d][b]) fprintf(fh, "%d,r2,%d,%.2f,%ld\n", static_cast<int>(chv[ic]), d, -1.0 + 0.05 * b, hR[ic][d][b]);
    }
  fclose(fh);
  if (haveT) {
    FILE* fq = fopen((std::string(prefix) + "_spurious.csv").c_str(), "w");
    fprintf(fq, "ch,G_mV,events,N_gt0p5,N_gt1p5,delta,n_pairs,n_spurious_A2gt0p5,n_spurious_A2gt1p5\n");
    FILE* fe = fopen((std::string(prefix) + "_E_hist.csv").c_str(), "w");
    fprintf(fe, "ch,delta,bin_lo_pe,n\n");
    for (int ic = 0; ic < NC; ++ic)
      for (int d = 3; d <= DMAX; ++d) {
        fprintf(fq, "%d,%.4f,%ld,%ld,%ld,%d,%ld,%ld,%ld\n", static_cast<int>(chv[ic]), gv[ic], used, n05[ic], n15[ic], d, spur[ic][d][0], spur[ic][d][1], spur[ic][d][2]);
        for (int b = 0; b < NBE; ++b)
          if (hE[ic][d][b]) fprintf(fe, "%d,%d,%.2f,%ld\n", static_cast<int>(chv[ic]), d, -1.0 + 0.05 * b, hE[ic][d][b]);
      }
    fclose(fq);
    fclose(fe);
    for (int ic = 0; ic < NC; ++ic) {
      long s5 = 0, s15 = 0, sa5 = 0, sa15 = 0;
      for (int d = 3; d <= 6; ++d) { s5 += spur[ic][d][1]; s15 += spur[ic][d][2]; }
      for (int d = 3; d <= DMAX; ++d) { sa5 += spur[ic][d][1]; sa15 += spur[ic][d][2]; }
      printf("ch%d espurios (E<0.5G) Δ3-6: %ld (A2>1.5: %ld); Δ3-12: %ld (A2>1.5: %ld)\n", static_cast<int>(chv[ic]), s5, s15, sa5, sa15);
    }
  }
  printf("eventos de tiempo bueno: %ld\n", used);
  for (int ic = 0; ic < NC; ++ic) {
    long ind = 0, tail = 0, ind15 = 0, tail15 = 0;
    for (int d = 3; d <= 6; ++d) { ind += cnt[ic][d][0]; tail += cnt[ic][d][1]; ind15 += cnt[ic][d][2]; tail15 += cnt[ic][d][3]; }
    printf("ch%d G=%.3f N05=%ld N15=%ld | Δ3-6: independientes %ld (A2>1.5: %ld), cola %ld (A2>1.5: %ld) | T(3..6)=%.3f %.3f %.3f %.3f (n_aisl=%ld)\n",
           static_cast<int>(chv[ic]), gv[ic], n05[ic], n15[ic], ind, ind15, tail, tail15, T[ic][3], T[ic][4], T[ic][5], T[ic][6], tn[ic]);
  }
}
