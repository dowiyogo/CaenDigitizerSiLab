// validate_hv_snapshots.C — validación 1 de DAQ2-5 (solo lectura).
// Compara cada fila del árbol hv de un run de daqrun con la lectura de hvmon del mismo read_seq en
// los archivos de control lento (hv_ y hv_full_ del día, leídos con libdaq_screader.so): Pw, Status,
// V0Set, VMon, IMon de los 12 canales.
// Uso: root -l -b -q 'daq/tools/validation/validate_hv_snapshots.C("<dir>/<base>_seg*.root", "20261019")'
R__ADD_INCLUDE_PATH(daq/io/include)
R__LOAD_LIBRARY(daq/build/io/libdaq_screader.so)
#include <TChain.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <tuple>

#include "daq/io/sc_reader.h"

void validate_hv_snapshots(const char* pattern, const char* day, const char* scDir = "/data/daq/slowcontrol") {
  struct Row {
    double v0, vmon, imon;
    UInt_t status, pw;
  };
  std::map<std::pair<ULong64_t, int>, Row> sc;  // (read_seq, ch) -> valores
  // read_seq vuelve a 1 en cada arranque de hvmon: solo filas desde el inicio del run (− 10 s).
  Long64_t tMin = 0;
  {
    TChain h0("hv");
    h0.Add(pattern);
    tMin = static_cast<Long64_t>(h0.GetMinimum("t_resp_utc_ns")) - 10000000000LL;
  }
  for (const char* kind : {"hv", "hv_full"}) {
    std::string err;
    TTree* t = daq::io::readSlowControl(scDir, kind, day, day, &err);
    if (!t) {
      printf("%s: %s\n", kind, err.c_str());
      continue;
    }
    ULong64_t seq = 0;
    Int_t ch = 0;
    Double_t v0 = 0, vm = 0, im = 0;
    UInt_t st = 0, pw = 0;
    Long64_t tr = 0;
    t->SetBranchAddress("t_resp_utc_ns", &tr);
    t->SetBranchAddress("read_seq", &seq);
    t->SetBranchAddress("ch", &ch);
    t->SetBranchAddress("V0Set_V", &v0);
    t->SetBranchAddress("VMon_V", &vm);
    t->SetBranchAddress("IMon_A", &im);
    t->SetBranchAddress("Status", &st);
    t->SetBranchAddress("Pw", &pw);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (tr >= tMin) sc[{seq, ch}] = Row{v0, vm, im, st, pw};
    }
    printf("%s_%s: %lld filas\n", kind, day, t->GetEntries());
    delete t;
  }
  TChain hv("hv");
  hv.Add(pattern);
  ULong64_t readSeq = 0;
  Bool_t alive = false;
  Double_t v0[12], vm[12], im[12];
  UInt_t st[12], pw[12];
  hv.SetBranchAddress("read_seq", &readSeq);
  hv.SetBranchAddress("alive", &alive);
  hv.SetBranchAddress("v0set_V", v0);
  hv.SetBranchAddress("vmon_V", vm);
  hv.SetBranchAddress("imon_A", im);
  hv.SetBranchAddress("status", st);
  hv.SetBranchAddress("pw", pw);
  long rows = 0, aliveRows = 0, matched = 0, missing = 0, mismatch = 0;
  double maxRelV = 0, maxAbsI = 0;
  for (Long64_t i = 0; i < hv.GetEntries(); ++i) {
    hv.GetEntry(i);
    ++rows;
    if (!alive) continue;
    ++aliveRows;
    bool rowOk = true, found = true;
    for (int c = 0; c < 12; ++c) {
      auto it = sc.find({readSeq, c});
      if (it == sc.end()) {
        found = false;
        break;
      }
      const Row& r = it->second;
      const double dv = std::max(std::fabs(r.vmon - vm[c]), std::fabs(r.v0 - v0[c]));
      const double di = std::fabs(r.imon - im[c]);
      maxRelV = std::max(maxRelV, dv);
      maxAbsI = std::max(maxAbsI, di);
      if (r.pw != pw[c] || r.status != st[c] || dv > 1e-6 * std::max(1.0, std::fabs(r.vmon)) || di > 1e-12) rowOk = false;
    }
    if (!found) {
      ++missing;
      continue;
    }
    ++matched;
    if (!rowOk) ++mismatch;
  }
  printf("árbol hv: %ld filas, %ld con hvmon vivo; con read_seq en el control lento: %ld, sin él: %ld; distintas: %ld\n", rows,
         aliveRows, matched, missing, mismatch);
  printf("diferencia máxima: V %.3g V, IMon %.3g A\n", maxRelV, maxAbsI);
}
