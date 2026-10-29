// timing_interp.C — DAQ2-22 etapa 3 (solo lectura de los crudos): tiempo interpolado de los pulsos grandes de hod2018.
// Por evento y canal de hod2018 (sin ch3): forma de onda en mV, pedestal = media recortada ±2 mV de [0,40) (como ntuple_extract.C),
// mínimo en la ventana signal [49,58); si la amplitud pasa thr_pe·G del canal, se guardan:
//   t_par = k + 0.5 (y[k−1] − y[k+1]) / (y[k−1] − 2 y[k] + y[k+1])  (parábola por el mínimo y sus vecinos, en muestras),
//   t_cfd = cruce del 50 % de la amplitud en el flanco de bajada, interpolado linealmente (en muestras),
// con ttt_raw del evento (fase del trigger respecto del reloj de muestreo).
// Salida CSV: seg, evt, ttt_raw, ch, amp_mV, k, t_par, t_cfd.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'timing_interp.C+("<archivos separados por comas>", "<G por canal 0-31, comas>", thr_pe, "<salida.csv>")'
#include <TFile.h>
#include <TObjString.h>
#include <TString.h>
#include <TTree.h>

#include <cstdio>
#include <vector>

#include "daq/pulse_finder.h"

void timing_interp(const char* files, const char* gains, double thrPe, const char* out) {
  std::vector<double> G;
  {
    TObjArray* g = TString(gains).Tokenize(",");
    for (int i = 0; i < g->GetEntries(); ++i) G.push_back(((TObjString*)g->At(i))->GetString().Atof());
  }
  const double lsb = daq::pulse::lsbMv();
  FILE* fo = fopen(out, "w");
  fprintf(fo, "seg,evt,ttt_raw,ch,amp_mV,k,t_par,t_cfd\n");
  TObjArray* a = TString(files).Tokenize(",");
  static UShort_t adc[64][144];
  double x[144];
  for (int f = 0; f < a->GetEntries(); ++f) {
    TFile fi(((TObjString*)a->At(f))->GetString(), "READ");
    TTree* t = (TTree*)fi.Get("events");
    UInt_t raw;
    ULong64_t seq;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"adc", "ttt_raw", "evt_seq"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("adc", adc);
    t->SetBranchAddress("ttt_raw", &raw);
    t->SetBranchAddress("evt_seq", &seq);
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      for (int c = 0; c < 32; ++c) {
        if (c == 3) continue;
        for (int s = 0; s < 144; ++s) x[s] = adc[c][s] * lsb;
        const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0);
        int k = 49;
        for (int s = 50; s < 58; ++s)
          if (x[s] < x[k]) k = s;
        const double amp = bl - x[k];
        if (amp < thrPe * G[c]) continue;
        const double den = x[k - 1] - 2 * x[k] + x[k + 1];
        const double tp = den != 0 ? k + 0.5 * (x[k - 1] - x[k + 1]) / den : k;
        const double lev = bl - 0.5 * amp;
        double tc = -1;
        for (int j = k - 1; j >= 40; --j)
          if (x[j] > lev) {
            tc = j + (x[j] - lev) / (x[j] - x[j + 1]);
            break;
          }
        fprintf(fo, "%d,%llu,%u,%d,%.3f,%d,%.4f,%.4f\n", f + 1, seq, raw, c, amp, k, tp, tc);
      }
    }
  }
  fclose(fo);
}
