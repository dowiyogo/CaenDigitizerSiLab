// profile_4b.C — XPLUS-4b etapa 4b: perfil temporal de impactos a lo largo del registro, por coordenada de hod2018.
// Para cada ventana de 9 muestras [w, w+9) que recorre el registro, la fracción de eventos con al menos un canal de la coordenada
// (X = ch 0-15 sin ch3, Y = ch 16-31) con un pulso del vector de la ntupla (pedestal local, separación 3, piso 1 mV) > T·G.
// Las ventanas que se solapan con la de señal [s0, s0+9) se descartan (la ntupla no guarda pulsos ahí); la ventana de señal se da
// aparte con aS > T·G (definición del estimador), y la fuera de tiempo del estimador con aC > T·G (mismo umbral).
// Eventos: todos, o los de tiempo bueno si se pasa un archivo de intervalos UTC excluidos ("t0 t1" por línea, de may_common).
// s0 < 0: SW trigger, sin ventana de señal.
// Salida: <prefijo>_profile.csv (coord, w, n_events, n_hit, frac) y <prefijo>_windows.csv (coord, signal_aS, ctl_aC con n y fracción).
// Uso: root -l -b -q 'analysis/night_2026-10-25/profile_4b.C+("<ntupla.root>", "<G ch0-31, coma>", T, s0, "<excluidos.txt o ''>", "<prefijo>")'
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

void profile_4b(const char* ntuple, const char* gains, double T, int s0, const char* excluded, const char* prefix) {
  std::vector<double> G;
  {
    std::stringstream ss(gains);
    std::string tok;
    while (std::getline(ss, tok, ',')) G.push_back(std::stod(tok));
  }
  std::vector<std::pair<double, double>> ex;
  if (excluded && *excluded) {
    std::ifstream in(excluded);
    double a, b;
    while (in >> a >> b) ex.emplace_back(a, b);
  }
  TFile f(ntuple);
  auto* t = static_cast<TTree*>(f.Get("ev"));
  std::vector<UChar_t>* pch = nullptr;
  std::vector<UShort_t>* psmp = nullptr;
  std::vector<float>* pamp = nullptr;
  int reglen = 0;
  double tU = 0;
  float aS[64], aC[64];
  t->SetBranchAddress("pch", &pch);
  t->SetBranchAddress("psmp", &psmp);
  t->SetBranchAddress("pamp", &pamp);
  t->SetBranchAddress("reglen", &reglen);
  t->SetBranchAddress("t_utc_s", &tU);
  t->SetBranchAddress("aS", aS);
  t->SetBranchAddress("aC", aC);
  const int nMax = 1024;
  std::vector<std::vector<long>> hit(2, std::vector<long>(nMax, 0));
  std::vector<long> nEvW(nMax, 0);
  long nEv = 0, sigHit[2] = {0, 0}, ctlHit[2] = {0, 0};
  std::vector<char> mark(nMax);
  for (Long64_t e = 0; e < t->GetEntries(); ++e) {
    t->GetEntry(e);
    bool bad = false;
    for (auto& [a, b] : ex)
      if (tU >= a && tU < b) bad = true;
    if (bad) continue;
    ++nEv;
    for (int w = 20; w + 9 <= reglen && w < nMax; ++w) ++nEvW[static_cast<std::size_t>(w)];
    for (int co = 0; co < 2; ++co) {
      std::fill(mark.begin(), mark.end(), 0);
      bool s = false, cc = false;
      for (int c = co * 16; c < co * 16 + 16; ++c) {
        if (c == 3) continue;
        const double thr = T * G[static_cast<std::size_t>(c)];
        if (aS[c] > thr) s = true;
        if (aC[c] > thr) cc = true;
      }
      if (s) ++sigHit[co];
      if (cc) ++ctlHit[co];
      for (std::size_t i = 0; i < pch->size(); ++i) {
        const int c = (*pch)[i];
        if (c == 3 || c / 16 != co) continue;
        if ((*pamp)[i] <= T * G[static_cast<std::size_t>(c)]) continue;
        const int k = (*psmp)[i];
        for (int w = std::max(20, k - 8); w <= k && w < nMax; ++w) mark[static_cast<std::size_t>(w)] = 1;
      }
      for (int w = 20; w + 9 <= reglen && w < nMax; ++w)
        if (mark[static_cast<std::size_t>(w)]) ++hit[static_cast<std::size_t>(co)][static_cast<std::size_t>(w)];
    }
  }
  FILE* fp = fopen((std::string(prefix) + "_profile.csv").c_str(), "w");
  fprintf(fp, "coord,w_sample,w_ns,n_events,n_hit,frac,overlaps_signal\n");
  for (int co = 0; co < 2; ++co)
    for (int w = 20; w < nMax; ++w) {
      if (nEvW[static_cast<std::size_t>(w)] == 0) continue;
      const bool ov = s0 >= 0 && w + 9 > s0 && w < s0 + 9;
      fprintf(fp, "%s,%d,%d,%ld,%ld,%.8f,%d\n", co ? "Y" : "X", w, w * 16, nEvW[static_cast<std::size_t>(w)], hit[static_cast<std::size_t>(co)][static_cast<std::size_t>(w)],
              static_cast<double>(hit[static_cast<std::size_t>(co)][static_cast<std::size_t>(w)]) / nEvW[static_cast<std::size_t>(w)], ov ? 1 : 0);
    }
  fclose(fp);
  FILE* fw = fopen((std::string(prefix) + "_windows.csv").c_str(), "w");
  fprintf(fw, "coord,n_events,signal_aS_hit,signal_frac,ctl_aC_hit,ctl_frac\n");
  for (int co = 0; co < 2; ++co)
    fprintf(fw, "%s,%ld,%ld,%.8f,%ld,%.8f\n", co ? "Y" : "X", nEv, sigHit[co], nEv ? static_cast<double>(sigHit[co]) / nEv : 0, ctlHit[co], nEv ? static_cast<double>(ctlHit[co]) / nEv : 0);
  fclose(fw);
  printf("eventos usados %ld\n", nEv);
}
