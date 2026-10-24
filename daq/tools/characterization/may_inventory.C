// may_inventory.C — PROC-1 (solo lectura). Inventario de un archivo de mayo (árbol data_dig): entradas, primer y último
// timestamp (µs del contador de la tarjeta), tiempo UTC con la época de referencia de hod2018_good_time_step2.C,
// muestras por registro y paso (vector time del primer evento), Δt mínimo entre triggers consecutivos y huecos > 20 s.
// Salida: una línea JSON por archivo en <salida> (se agrega).
// Uso: root -l -b -q 'daq/tools/characterization/may_inventory.C+("<archivo>", "<etiqueta>", "<salida.jsonl>")'
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

void may_inventory(const char* path, const char* label, const char* out) {
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  auto* t = f->Get<TTree>("data_dig");
  std::vector<float>* tw = nullptr;
  Long_t ts = 0;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("timestamp", 1);
  t->SetBranchAddress("timestamp", &ts);
  const Long64_t n = t->GetEntries();
  // primer evento: también el vector de tiempos
  t->SetBranchStatus("time", 1);
  t->SetBranchAddress("time", &tw);
  t->GetEntry(0);
  const int ns = tw ? static_cast<int>(tw->size()) : 0;
  const double dt = (tw && tw->size() > 1) ? (*tw)[1] - (*tw)[0] : NAN;
  const double tlast = (tw && !tw->empty()) ? tw->back() : NAN;
  t->SetBranchStatus("time", 0);
  const double bootRefEpoch = 1778886208.8;  // UTC de 2026-05-15T19:03:28.8-04:00 (hod2018_good_time_step2.C)
  Long_t first = 0, prev = 0, last = 0;
  double minDt = 1e30, gapSum = 0;
  long nGaps = 0, nBack = 0;
  for (Long64_t e = 0; e < n; ++e) {
    t->GetEntry(e);
    if (e == 0) first = ts;
    if (e > 0) {
      const double d = (ts - prev) * 1e-6;
      if (d < 0) ++nBack;
      if (d >= 0) minDt = std::min(minDt, d);
      if (d > 20) {
        ++nGaps;
        gapSum += d;
      }
    }
    prev = ts;
    last = ts;
  }
  FILE* o = fopen(out, "a");
  fprintf(o,
          "{\"label\": \"%s\", \"file\": \"%s\", \"entries\": %lld, \"samples\": %d, \"sample_ns\": %.3f, \"record_last_ns\": %.1f, "
          "\"ts_first_us\": %ld, \"ts_last_us\": %ld, \"utc_first_s\": %.3f, \"utc_last_s\": %.3f, \"duration_s\": %.3f, "
          "\"min_dt_s\": %.6f, \"gaps_gt20s\": %ld, \"gap_sum_s\": %.1f, \"backward_steps\": %ld}\n",
          label, path, n, ns, dt, tlast, first, last, bootRefEpoch + first * 1e-6, bootRefEpoch + last * 1e-6, (last - first) * 1e-6, minDt, nGaps,
          gapSum, nBack);
  fclose(o);
  printf("%s: %lld entradas\n", label, n);
}
