// orientation_corner_test.C — seguimiento 2026-10-31, punto 5 (preparada, NO ejecutada): prueba de orientación de hod2018 con una fuente
// o un LED en la esquina (0,0) (esquina inferior izquierda, sistema de René del 2026-10-30: barras X ∥ eje Y, barras Y ∥ eje X).
// Para usar cuando haya alguien en el laboratorio y el DAQ funcione (una tarea que autorice el digitalizador):
//   1. run de referencia sin fuente (SW trigger o trigger externo, unos minutos);
//   2. run con la fuente (β, p. ej. Sr-90 con colimador) o el LED en la esquina (0,0), mismas condiciones.
// La macro lee los dos runs (solo lectura), calcula por canal de hod2018 la fracción de eventos con un pulso > thr_mV en todo el registro
// (amplitud = pedestal − mínimo; pedestal = mediana de las primeras 40 muestras), el exceso fuente − referencia y su significancia, y
// compara las dos orientaciones posibles de cada coordenada:
//   orden medido (r4, sentido absoluto libre): X índice i = 0…15 → ch 15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14;
//                                              Y índice j = 0…15 → ch 16, 17, …, 31.
//   X «+» (i crece con x): la esquina x = 0 está en i = 0, 1, 2 → ch 15, 0, 13.   X «−»: i = 15, 14, 13 → ch 14, 1, 12.
//   Y «+» (j crece con y): la esquina y = 0 está en j = 0, 1, 2 → ch 16, 17, 18.  Y «−»: j = 15, 14, 13 → ch 31, 30, 29.
// (ch3, i = 12, está muerto; no entra en ninguna de las listas de esquina.) Decisión por coordenada: la orientación cuyo grupo de esquina
// tiene el mayor exceso, si la diferencia entre los dos grupos supera 5σ; si no, «no decide».
// Además, con una fuente β desde arriba la capa superior recibe más que la inferior: se comparan las barras de la esquina de cada capa
// (X: capa A = índices pares = canales impares, capa B = canales pares; Y: A = ch pares 16…30, B = impares 17…31) para identificar la capa
// superior. Con PROC-2 (along_bar.json): si Y es «+», la capa X superior debería ser la A (canales impares); queda como comprobación cruzada.
// Uso: root -l -b -q 'orientation_corner_test.C+("<run con fuente: archivos separados por comas>", "<run de referencia: idem>", 10.0)'
#include <TFile.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
const std::array<int, 16> kX = {15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14};
const double kLsb = 2000.0 / 4095.0;  // mV por cuenta (12 bit, 2 Vpp)

struct Counts {
  long n = 0;
  std::array<long, 32> k{};
};

Counts count(const char* files, double thr) {
  Counts c;
  TObjArray* a = TString(files).Tokenize(",");
  for (int i = 0; i < a->GetEntries(); ++i) {
    TFile f(((TObjString*)a->At(i))->GetString(), "READ");
    TTree* t = (TTree*)f.Get("events");
    if (!t) continue;
    const int len = t->GetLeaf("adc")->GetLenStatic() / 64;
    std::vector<UShort_t> adc(64 * len);
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("adc", 1);
    t->SetBranchAddress("adc", adc.data());
    std::vector<double> x(len);
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      ++c.n;
      for (int ch = 0; ch < 32; ++ch) {
        for (int s = 0; s < len; ++s) x[s] = adc[ch * len + s] * kLsb;
        std::vector<double> p(x.begin(), x.begin() + std::min(40, len));
        std::nth_element(p.begin(), p.begin() + p.size() / 2, p.end());
        const double bl = p[p.size() / 2];
        if (bl - *std::min_element(x.begin(), x.end()) > thr) ++c.k[ch];
      }
    }
  }
  return c;
}

double excess(const Counts& s, const Counts& r, const std::vector<int>& chs, double& err) {
  double ex = 0, v = 0;
  for (int ch : chs) {
    const double fs = double(s.k[ch]) / s.n, fr = double(r.k[ch]) / r.n;
    ex += fs - fr;
    v += fs * (1 - fs) / s.n + fr * (1 - fr) / r.n;
  }
  err = std::sqrt(v);
  return ex;
}

void decide(const char* name, const Counts& s, const Counts& r, const std::vector<int>& plus, const std::vector<int>& minus) {
  double ep, em;
  const double xp = excess(s, r, plus, ep), xm = excess(s, r, minus, em);
  const double z = (xp - xm) / std::hypot(ep, em);
  printf("%s: exceso esquina «+» (ch", name);
  for (int ch : plus) printf(" %d", ch);
  printf(") = %.4f ± %.4f; esquina «−» (ch", xp, ep);
  for (int ch : minus) printf(" %d", ch);
  printf(") = %.4f ± %.4f; diferencia %.1fσ -> %s\n", xm, em, z, z > 5 ? "orientación «+»" : (z < -5 ? "orientación «−»" : "no decide"));
}
}  // namespace

void orientation_corner_test(const char* sourceFiles, const char* refFiles, double thrMv = 10.0) {
  const Counts s = count(sourceFiles, thrMv), r = count(refFiles, thrMv);
  printf("eventos: fuente %ld, referencia %ld; umbral %.1f mV en todo el registro\n", s.n, r.n, thrMv);
  printf("canal  índice  f_fuente  f_ref  exceso/σ\n");
  for (int ch = 0; ch < 32; ++ch) {
    int idx = -1;
    for (int i = 0; i < 16; ++i)
      if ((ch < 16 && kX[i] == ch) || (ch >= 16 && ch - 16 == i)) idx = i;
    const double fs = double(s.k[ch]) / s.n, fr = double(r.k[ch]) / r.n;
    const double e = std::sqrt(fs * (1 - fs) / s.n + fr * (1 - fr) / r.n);
    printf("ch%02d   %s%2d     %.5f   %.5f  %+.1f\n", ch, ch < 16 ? "X" : "Y", idx, fs, fr, e > 0 ? (fs - fr) / e : 0.0);
  }
  decide("X", s, r, {15, 0, 13}, {14, 1, 12});
  decide("Y", s, r, {16, 17, 18}, {31, 30, 29});
  // capa superior: barra de esquina de cada capa en el extremo con más exceso (índices 0/1 o 15/14)
  double e0, e1, e15, e14;
  const double a0 = excess(s, r, {15}, e0), b1 = excess(s, r, {0}, e1), a15 = excess(s, r, {14}, e15), b14 = excess(s, r, {1}, e14);
  printf("X, capas en la esquina: índice 0 (capa A, ch15) %.4f ± %.4f, índice 1 (capa B, ch0) %.4f ± %.4f; índice 15 (capa B, ch14) %.4f ± %.4f, "
         "índice 14 (capa A, ch1) %.4f ± %.4f\n", a0, e0, b1, e1, a15, e15, b14, e14);
  double f0, f1, f15, f14;
  const double y0 = excess(s, r, {16}, f0), y1 = excess(s, r, {17}, f1), y15 = excess(s, r, {31}, f15), y14 = excess(s, r, {30}, f14);
  printf("Y, capas en la esquina: j=0 (capa A, ch16) %.4f ± %.4f, j=1 (capa B, ch17) %.4f ± %.4f; j=15 (capa B, ch31) %.4f ± %.4f, "
         "j=14 (capa A, ch30) %.4f ± %.4f\n", y0, f0, y1, f1, y15, f15, y14, f14);
  printf("Con fuente β desde arriba, la capa con más exceso en la esquina es la superior (comprobación cruzada con along_bar).\n");
}
