// manuscript_xy.C — PROC-1 etapa 2.5 / 2.8 (solo lectura): reproduce la definición de impacto del manuscrito (análisis de
// mayo, analysis/src/reconstruct_corrected.cxx: extract_hit_window, sin modificarlo) para un archivo de mayo:
//  - pedestal = media de las muestras con t < 2500 ns; altura = pedestal − mínimo en la ventana;
//  - ventana del muón [4000, 4400] ns; ventana de control [6000, 6400] ns (+2000 ns, kCtlShift);
//  - impacto de un plano: algún canal del plano con altura > 50 mV (hod2018) o > 35 mV (hod2019);
//  - planos (kCfg_2026_05, analysis/src/hodo_mapping.h): hod2018 X = ch0-15, Y = ch16-31; hod2019 X = ch48-63,
//    Y = ch32-47.
// Salida binaria: por evento {double t_us (timestamp del contador); uint8 bits}, con bits 0-7 = X18 en tiempo, Y18 en tiempo,
// X18 control, Y18 control, X19 en tiempo, Y19 en tiempo, X19 control, Y19 control (9 bytes por evento).
// Uso: root -l -b -q 'daq/tools/proceeding/manuscript_xy.C+("<archivo de mayo>", "<salida.bin>")'
#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <cstdio>
#include <memory>
#include <vector>

void manuscript_xy(const char* path, const char* out, double thr18 = 50.0, double thr19 = 35.0) {
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  auto* t = f->Get<TTree>("data_dig");
  std::vector<float>* tw = nullptr;
  std::vector<float>* cw[64] = {};
  Long_t ts = 0;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("time", 1);
  t->SetBranchAddress("time", &tw);
  t->SetBranchStatus("timestamp", 1);
  t->SetBranchAddress("timestamp", &ts);
  for (int c = 0; c < 64; ++c) {
    t->SetBranchStatus(Form("Ch%d", c), 1);
    t->SetBranchAddress(Form("Ch%d", c), &cw[c]);
  }
  t->SetCacheSize(256LL << 20);
  FILE* fo = fopen(out, "wb");
  const Long64_t ne = t->GetEntries();
  for (Long64_t e = 0; e < ne; ++e) {
    t->GetEntry(e);
    std::uint8_t bits = 0;
    const auto& tv = *tw;
    const int n = static_cast<int>(tv.size());
    for (int c = 0; c < 64; ++c) {
      if (!cw[c] || static_cast<int>(cw[c]->size()) != n) continue;
      const auto& wf = *cw[c];
      double s = 0;
      int nb = 0;
      for (int i = 0; i < n; ++i)
        if (tv[i] < 2500.0) {
          s += wf[i];
          ++nb;
        }
      if (nb < 2) continue;
      const double ped = s / nb;
      double vs = 1e30, vc = 1e30;
      for (int i = 0; i < n; ++i) {
        if (tv[i] >= 4000.0 && tv[i] <= 4400.0 && wf[i] < vs) vs = wf[i];
        if (tv[i] >= 6000.0 && tv[i] <= 6400.0 && wf[i] < vc) vc = wf[i];
      }
      const bool st18 = c < 32, x = c < 16 || c >= 48;
      const double thr = st18 ? thr18 : thr19;
      const int base = st18 ? 0 : 4;
      if (ped - vs > thr) bits |= static_cast<std::uint8_t>(1u << (base + (x ? 0 : 1)));
      if (ped - vc > thr) bits |= static_cast<std::uint8_t>(1u << (base + 2 + (x ? 0 : 1)));
    }
    const double tu = static_cast<double>(ts);
    fwrite(&tu, sizeof(double), 1, fo);
    fwrite(&bits, 1, 1, fo);
  }
  fclose(fo);
  printf("eventos %lld\n", ne);
}
