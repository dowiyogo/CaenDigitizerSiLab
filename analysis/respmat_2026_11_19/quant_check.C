// quant_check.C — seguimiento 2026-11-19: ¿las muestras de mayo (mV en float) están cuantizadas en el paso del ADC (2000/4095 mV)?
// Lee los primeros nev eventos de ch0-31 (READ) y escribe el histograma de la parte fraccionaria de (x/lsb − x0/lsb) por canal.
// Uso: root -l -b -q 'analysis/respmat_2026_11_19/quant_check.C("<archivo mayo>", 200)'
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <vector>

void quant_check(const char* file, int nev) {
  const double lsb = 2000.0 / 4095.0;
  TFile f(file, "READ");
  auto* t = static_cast<TTree*>(f.Get("data_dig"));
  std::vector<float>* w[32] = {};
  t->SetBranchStatus("*", 0);
  for (int c = 0; c < 32; ++c) t->SetBranchStatus(Form("Ch%d", c), 1), t->SetBranchAddress(Form("Ch%d", c), &w[c]);
  long nall = 0, nint = 0;
  double maxdev = 0;
  for (int e = 0; e < nev; ++e) {
    t->GetEntry(e);
    for (int c = 0; c < 32; ++c) {
      const double r0 = (*w[c])[0] / lsb;
      for (float v : *w[c]) {
        const double d = v / lsb - r0;
        const double fr = std::fabs(d - std::round(d));
        ++nall;
        if (fr < 1e-3) ++nint;
        maxdev = std::max(maxdev, fr);
      }
    }
  }
  printf("muestras %ld, con (x - x0)/lsb entero a 1e-3: %ld (fracción %.6f), máx desvío %.6f\n", nall, nint, double(nint) / nall, maxdev);
  t->GetEntry(0);
  printf("ch0, primeras 8 muestras (mV) y x/lsb:");
  for (int i = 0; i < 8; ++i) printf(" %.6f/%.4f", (*w[0])[i], (*w[0])[i] / lsb);
  printf("\n");
}
