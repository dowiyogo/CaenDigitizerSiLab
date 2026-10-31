// h19_timeline.C — seguimiento 2026-10-31 (5), punto 1 (solo lectura de una ntupla): pulsos de hod2019 (ch 32-63) y de hod2018 (ch 0-31, sin ch3)
// fuera de la ventana de señal (listas pch/psmp/pamp de ntuple_extract.C), por minuto UTC: eventos, pulsos > 1 mV y > 3 mV por evento de cada
// estación, y el número de canales de hod2019 con algún pulso > 3 mV en ese minuto. Salida CSV en stdout.
// Uso: root -l -b -q 'h19_timeline.C+("<ntupla.root>", t0_utc_s, t1_utc_s)'
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

void h19_timeline(const char* nt, double t0, double t1) {
  TFile f(nt, "READ");
  TTree* t = (TTree*)f.Get("ev");
  std::vector<unsigned char>* pch = nullptr;
  std::vector<float>* pamp = nullptr;
  double tu;
  t->SetBranchStatus("*", 0);
  for (const char* b : {"pch", "pamp", "t_utc_s"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pch", &pch);
  t->SetBranchAddress("pamp", &pamp);
  t->SetBranchAddress("t_utc_s", &tu);
  struct M { long ev = 0, h19_1 = 0, h19_3 = 0, h18_1 = 0; std::set<int> ch3; };
  std::map<long, M> m;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (tu < t0 || tu >= t1) continue;
    auto& x = m[static_cast<long>(std::floor(tu / 60)) * 60];
    ++x.ev;
    for (size_t k = 0; k < pch->size(); ++k) {
      const int c = (*pch)[k];
      const float a = (*pamp)[k];
      if (c >= 32) {
        if (a >= 1) ++x.h19_1;
        if (a >= 3) { ++x.h19_3; x.ch3.insert(c); }
      } else if (c != 3 && a >= 1) ++x.h18_1;
    }
  }
  printf("minute_utc_s,events,h19_gt1mV_per_ev,h19_gt3mV_per_ev,h19_channels_gt3mV,h18_gt1mV_per_ev\n");
  for (auto& [k, x] : m)
    printf("%ld,%ld,%.4f,%.4f,%zu,%.2f\n", k, x.ev, double(x.h19_1) / x.ev, double(x.h19_3) / x.ev, x.ch3.size(), double(x.h18_1) / x.ev);
}
