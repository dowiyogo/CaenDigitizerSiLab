// ntuple_extract.C — XPLUS-4b etapa 1 (solo lectura de los crudos): ntupla reducida, una pasada por archivo crudo.
// Formatos: "may" (árbol data_dig, Ch0..Ch63 vector<float> en mV, timestamp en µs), "daqv2" (árbol events, adc UShort_t[64·n],
// t_utc_ns) y "sw" (daqv2 de SW trigger: sin ventana de señal, todos los pulsos desde la muestra 20).
// Definiciones de impacto idénticas a daq/tools/characterization/hitdef_extract.C (para que el fixture se reproduzca):
//   tramo de 144 muestras desde s0 − 49 (mayo s0 = 258 = 4128 ns; daqv2 s0 = 49); pedestal = media recortada ±2 mV de [0,40)
//   del tramo; signal [49,58), fuera de tiempo [99,108) del tramo; findPulses umbral 0.5 mV, separación 4; aS, aC = mayor
//   amplitud; aG/kG = mínimo del tramo; aM/aMc = definición del análisis de mayo.
// Además, por canal y evento: RMS recortado del pedestal, muestra del pulso de aS (kS, en coordenadas del registro), y todos los
// pulsos del registro fuera de la ventana de señal y de las 20 primeras muestras (sw: desde la 20, sin exclusión):
//   regiones [20, s0) y [s0 + 9, n), cada una con su pedestal local = media recortada ±2 mV de la propia región; buscador del DAQ
//   (pulse_finder.h) con separación 3 muestras (tiempo muerto de 48 ns) y piso 1.0 mV; se guardan canal, muestra y amplitud.
// Espectro oscuro por hora UTC y canal con el método del reajuste de G de hitdef_extract.C (pulsos en [0,49) ∪ [102,144) del
// tramo, pedestal del tramo, umbral 0.9 mV, separación 4; bins de 0.25 mV hasta 100 mV), para G por bloque con gain_fit.h.
// Salida: <prefijo>.root (árbol "ev": t_utc_s, ev, ped[64], rms[64], aS[64], kS[64], aC[64], npul, pch[npul], psmp[npul], pamp[npul],
//         reglen, pre_len, post_len; árbol "spec": hour_utc_s, ch, events, counts[401]); <prefijo>_hits.bin con el registro exacto de
//         hitdef_extract.C (nch 64) para mayo y daqv2.
// Uso: root -l -b -q -e 'gSystem->AddIncludePath("-I<repo>/daq/core/include")' \
//        'analysis/night_2026-10-25/ntuple_extract.C+("<archivo(s)>", "may"|"daqv2"|"sw", "<prefijo>"[, XPLUS4])'
#include <TChain.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "daq/pulse_finder.h"

namespace {
constexpr int kSpan = 144;
constexpr int kNch = 64;
constexpr int kNb = 400;  // 0.25 mV hasta 100 mV, más desborde
double trimmedRms(const double* x, daq::pulse::Region r, double mean, double trim) {
  const double m = daq::pulse::median(x, r);
  double s = 0;
  int n = 0;
  for (int i = r.begin; i < r.end; ++i)
    if (std::fabs(x[i] - m) <= trim) {
      s += (x[i] - mean) * (x[i] - mean);
      ++n;
    }
  return n > 1 ? std::sqrt(s / n) : 0.0;
}
}  // namespace

void ntuple_extract(const char* pattern, const char* format, const char* prefix, int /*marca XPLUS4*/ = 0) {
  const std::string fmt(format);
  const bool may = fmt == "may", sw = fmt == "sw";
  const int s0 = may ? 258 : 49;
  const int off = s0 - 49;
  const double lsb = daq::pulse::lsbMv();
  TChain t(may ? "data_dig" : "events");
  t.Add(pattern);
  t.SetBranchStatus("*", 0);
  std::vector<float>* cw[kNch] = {};
  std::vector<UShort_t> adc;
  Long_t ts = 0;
  Long64_t tutc = 0;
  int n = 0;
  if (may) {
    for (int c = 0; c < kNch; ++c) {
      t.SetBranchStatus(Form("Ch%d", c), 1);
      t.SetBranchAddress(Form("Ch%d", c), &cw[c]);
    }
    t.SetBranchStatus("timestamp", 1);
    t.SetBranchAddress("timestamp", &ts);
    t.SetCacheSize(256LL << 20);
  } else {
    n = t.GetLeaf("adc")->GetLenStatic() / kNch;
    adc.resize(static_cast<std::size_t>(kNch * n));
    t.SetBranchStatus("adc", 1);
    t.SetBranchStatus("t_utc_ns", 1);
    t.SetBranchAddress("adc", adc.data());
    t.SetBranchAddress("t_utc_ns", &tutc);
  }
  const double bootRefEpoch = 1778886208.8;  // como hitdef_extract.C
  TFile fo((std::string(prefix) + ".root").c_str(), "RECREATE");
  fo.SetCompressionSettings(505);  // ZSTD nivel 5
  TTree ev("ev", "XPLUS-4b ntupla reducida (amplitudes en mV, muestras de 16 ns)");
  double tU = 0;
  Long64_t evn = 0;
  float ped[kNch], rms[kNch], aS[kNch], aC[kNch];
  UShort_t kS[kNch];
  int npul = 0, reglen = 0, preLen = 0, postLen = 0;
  std::vector<UChar_t> pch;
  std::vector<UShort_t> psmp;
  std::vector<float> pamp;
  ev.Branch("t_utc_s", &tU, "t_utc_s/D");
  ev.Branch("ev", &evn, "ev/L");
  ev.Branch("ped", ped, "ped[64]/F");
  ev.Branch("rms", rms, "rms[64]/F");
  ev.Branch("aS", aS, "aS[64]/F");
  ev.Branch("kS", kS, "kS[64]/s");
  ev.Branch("aC", aC, "aC[64]/F");
  ev.Branch("reglen", &reglen, "reglen/I");
  ev.Branch("pre_len", &preLen, "pre_len/I");
  ev.Branch("post_len", &postLen, "post_len/I");
  ev.Branch("npul", &npul, "npul/I");
  ev.Branch("pch", &pch);
  ev.Branch("psmp", &psmp);
  ev.Branch("pamp", &pamp);
  // espectro oscuro por hora y canal (método de hitdef_extract.C)
  std::map<std::pair<long, int>, std::pair<long, std::vector<unsigned>>> spec;
  FILE* fb = sw ? nullptr : fopen((std::string(prefix) + "_hits.bin").c_str(), "wb");
  const std::vector<daq::pulse::Region> sig{{49, 58}}, ctl{{99, 108}}, oot{{0, 49}, {102, kSpan}};
  std::vector<daq::pulse::Pulse> p;
  std::vector<double> x, full;
  const Long64_t ne = t.GetEntries();
  Long64_t nw = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    tU = may ? bootRefEpoch + static_cast<double>(ts) * 1e-6 : static_cast<double>(tutc) * 1e-9;
    evn = e;
    bool okEv = true;
    float aG[kNch] = {}, aM[kNch] = {}, aMc[kNch] = {};
    uint8_t kG[kNch] = {};
    pch.clear();
    psmp.clear();
    pamp.clear();
    const long hour = static_cast<long>(std::floor(tU / 3600.0)) * 3600;
    for (int c = 0; c < kNch; ++c) {
      int len = may ? (cw[c] ? static_cast<int>(cw[c]->size()) : 0) : n;
      if (may && len < off + kSpan) {
        okEv = false;
        break;
      }
      full.resize(static_cast<std::size_t>(len));
      if (may)
        for (int i = 0; i < len; ++i) full[static_cast<std::size_t>(i)] = (*cw[c])[static_cast<std::size_t>(i)];
      else
        for (int i = 0; i < len; ++i) full[static_cast<std::size_t>(i)] = adc[static_cast<std::size_t>(c * n + i)] * lsb;
      reglen = len;
      aS[c] = aC[c] = 0;
      kS[c] = 0;
      ped[c] = rms[c] = 0;
      if (!sw) {
        x.assign(full.begin() + off, full.begin() + off + kSpan);
        const double bl = daq::pulse::trimmedMean(x.data(), daq::pulse::Region{0, 40}, 2.0);
        ped[c] = static_cast<float>(bl);
        rms[c] = static_cast<float>(trimmedRms(x.data(), daq::pulse::Region{0, 40}, bl, 2.0));
        p.clear();
        daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, sig, p);
        for (const auto& u : p)
          if (static_cast<float>(u.amplitudeMv) > aS[c]) {
            aS[c] = static_cast<float>(u.amplitudeMv);
            kS[c] = static_cast<UShort_t>(u.sample + off);
          }
        p.clear();
        daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.5, 4, ctl, p);
        for (const auto& u : p) aC[c] = std::max(aC[c], static_cast<float>(u.amplitudeMv));
        const int k = static_cast<int>(std::min_element(x.begin(), x.end()) - x.begin());
        aG[c] = static_cast<float>(bl - x[static_cast<std::size_t>(k)]);
        kG[c] = static_cast<uint8_t>(k);
        for (int w = 0; w < 2; ++w) {
          const int a = w == 0 ? 49 : 99;
          int km = a;
          for (int i = a; i < a + 10; ++i)
            if (x[static_cast<std::size_t>(i)] < x[static_cast<std::size_t>(km)]) km = i;
          double pd = 0;
          for (int j = -12; j <= -4; ++j) pd += x[static_cast<std::size_t>(km + j)];
          (w == 0 ? aM : aMc)[c] = static_cast<float>(pd / 9.0 - x[static_cast<std::size_t>(km)]);
        }
        p.clear();
        daq::pulse::findPulses(x.data(), kSpan, bl, 1.0, 0.9, 4, oot, p);
        auto& sp = spec[{hour, c}];
        if (sp.second.empty()) sp.second.assign(kNb + 1, 0);
        ++sp.first;
        for (const auto& u : p) ++sp.second[static_cast<std::size_t>(std::min(kNb, static_cast<int>(u.amplitudeMv / 0.25)))];
      }
      // pulsos fuera de la ventana de señal, cada región con su pedestal local
      std::vector<daq::pulse::Region> regs;
      if (sw)
        regs.push_back({20, len});
      else {
        regs.push_back({20, s0});
        regs.push_back({s0 + 9, len});
      }
      preLen = sw ? len - 20 : s0 - 20;
      postLen = sw ? 0 : len - (s0 + 9);
      for (const auto& rg : regs) {
        if (rg.size() <= 0) continue;
        const double lb = daq::pulse::trimmedMean(full.data(), rg, 2.0);
        p.clear();
        daq::pulse::findPulses(full.data(), len, lb, 1.0, 1.0, 3, std::vector<daq::pulse::Region>{rg}, p);
        for (const auto& u : p) {
          pch.push_back(static_cast<UChar_t>(c));
          psmp.push_back(static_cast<UShort_t>(u.sample));
          pamp.push_back(static_cast<float>(u.amplitudeMv));
        }
      }
    }
    if (!okEv) continue;
    npul = static_cast<int>(pch.size());
    ev.Fill();
    if (fb) {
      fwrite(&tU, sizeof(double), 1, fb);
      fwrite(aS, sizeof(float), kNch, fb);
      fwrite(aC, sizeof(float), kNch, fb);
      fwrite(aG, sizeof(float), kNch, fb);
      fwrite(kG, sizeof(uint8_t), kNch, fb);
      fwrite(aM, sizeof(float), kNch, fb);
      fwrite(aMc, sizeof(float), kNch, fb);
    }
    ++nw;
    if (nw % 20000 == 0) {
      printf("%lld/%lld\n", e, ne);
      fflush(stdout);
    }
  }
  if (fb) fclose(fb);
  TTree st("spec", "espectro oscuro por hora y canal (método de hitdef_extract.C, 0.25 mV por bin)");
  Long64_t hr = 0;
  int ch = 0;
  Long64_t nev = 0;
  unsigned counts[kNb + 1];
  st.Branch("hour_utc_s", &hr, "hour_utc_s/L");
  st.Branch("ch", &ch, "ch/I");
  st.Branch("events", &nev, "events/L");
  st.Branch("counts", counts, Form("counts[%d]/i", kNb + 1));
  for (auto& [k, v] : spec) {
    hr = k.first;
    ch = k.second;
    nev = v.first;
    std::copy(v.second.begin(), v.second.end(), counts);
    st.Fill();
  }
  fo.Write();
  fo.Close();
  printf("eventos %lld, escritos %lld, formato %s\n", ne, nw, format);
}
