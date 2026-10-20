// hist_util.h
// Utilidades para histogramas que se rellenan con SetBinContent desde contadores propios (páginas
// de daqrun). SetBinContent suma 1 a las entradas del histograma en cada llamada, así que tras N
// refrescos "Entries" vale N × bins (DAQ2-8: 515973120 = 13056 bins × 39520 refrescos en
// pulses/amp). refreshEntries() deja Entries = suma de los contenidos (sin desbordes).
#pragma once

#include <TH1.h>

namespace daq::io {

inline void refreshEntries(TH1& h) {
  double sum = 0;
  const int nx = h.GetNbinsX(), ny = h.GetDimension() > 1 ? h.GetNbinsY() : 1, nz = h.GetDimension() > 2 ? h.GetNbinsZ() : 1;
  for (int z = 1; z <= nz; ++z)
    for (int y = 1; y <= ny; ++y)
      for (int x = 1; x <= nx; ++x) sum += h.GetBinContent(h.GetBin(x, h.GetDimension() > 1 ? y : 0, h.GetDimension() > 2 ? z : 0));
  h.SetEntries(sum);
}

}  // namespace daq::io
