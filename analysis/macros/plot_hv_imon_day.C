// IMon de CH00-07 de un día de control lento de hvmon (lector daqv2: plano o .gz, partes .N).
// Desde la raíz del repo: root -l -b -q 'analysis/macros/plot_hv_imon_day.C("20261019")'
R__ADD_INCLUDE_PATH(daq/io/include)
R__LOAD_LIBRARY(daq/build/io/libdaq_screader.so)
#include "daq/io/sc_reader.h"
void plot_hv_imon_day(const char* day = "20261019", const char* dir = "/data/daq/slowcontrol",
                      const char* out = "outputs/hv_imon_day.png") {
  std::string err;
  TTree* t = daq::io::readSlowControl(dir, "hv", day, day, &err);
  if (!t) { printf("%s\n", err.c_str()); return; }
  TCanvas c("c", Form("IMon CH00-07 %s", day), 1200, 600);
  t->Draw("IMon_A*1e6:(t_resp_utc_ns%86400000000000)/3.6e12:ch", "ch<8", "colz");  // µA vs horas UTC
  c.SaveAs(out);
}
