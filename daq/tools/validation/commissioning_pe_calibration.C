// commissioning_pe_calibration.C — DAQ2-10, parte D (solo lectura de los datos). Calibración de p.e.
// de octubre desde un run de SW trigger con registro completo y todos los grupos encendidos.
//  - Por canal: espectro de amplitud de los pulsos oscuros > 2 LSB en [0,49) + [102,622) (buscador
//    del DAQ, separación 4), paso 1 LSB hasta 100 mV; pedestal = media de la línea base recortada de
//    [0,40) (mV).
//  - hod2018: G con daq/gain_fit.h (reglas de la fase 5: picos 5-26 mV, separación, G/2 y 2G
//    ambiguos contra la semilla, cambio > 20 % "changed"); semilla = G mediana del canal en la noche
//    del 2026-10-19 (gains.json de daqoffline) o 4.0 mV. Provisional si el estado no es "ok".
//  - hod2019: tasa oscura > 3 mV (indicador de ganancia); sin G, provisional.
// Salida: JSON daqv2.pe_calibration/1 con versión, run de origen, estado de HV y temperatura de la
// placa al inicio del run (hvmon_state_at_start del segmento).
// Uso: root -l -b -q 'daq/tools/validation/commissioning_pe_calibration.C+("<run>_seg*.root", "<version>", "<salida.json>")'
R__ADD_INCLUDE_PATH(daq/core/include)
#include <TChain.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TNamed.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/gain_fit.h"
#include "daq/pulse_finder.h"

void commissioning_pe_calibration(const char* pattern, const char* version, const char* outJson,
                                  const char* seeds = "/data/daq/sim/offline/night_default/gains.json") {
  std::map<int, std::vector<double>> seedVals;
  {
    std::ifstream in(seeds);
    if (in) {
      auto g = nlohmann::json::parse(in);
      for (const auto& x : g)
        if (x["status"] == "ok" || x["status"] == "changed") seedVals[x["ch"].get<int>()].push_back(x["value"].get<double>());
    }
  }
  TChain t("events");
  t.Add(pattern);
  const int n = t.GetLeaf("adc")->GetLenStatic() / 64;
  std::vector<UShort_t> adc(64 * n);
  Long64_t tread = 0;
  t.SetBranchStatus("*", 0);
  t.SetBranchStatus("adc", 1);
  t.SetBranchStatus("t_read_utc_ns", 1);
  t.SetBranchAddress("adc", adc.data());
  t.SetBranchAddress("t_read_utc_ns", &tread);
  const double lsb = daq::pulse::lsbMv();
  const int nb = static_cast<int>(100.0 / lsb);
  std::vector<std::vector<std::uint32_t>> spec(64, std::vector<std::uint32_t>(static_cast<std::size_t>(nb), 0));
  std::vector<double> pedSum(64, 0);
  std::vector<long> dark3(64, 0);
  const std::vector<daq::pulse::Region> zones{{0, 49}, {102, std::min(622, n)}};
  const double exposure = (49 + std::min(622, n) - 102) * 16e-9;
  const Long64_t ne = t.GetEntries();
  Long64_t t0 = 0, t1 = 0;
  for (Long64_t e = 0; e < ne; ++e) {
    t.GetEntry(e);
    if (e == 0) t0 = tread;
    t1 = tread;
    for (int c = 0; c < 64; ++c) {
      const UShort_t* x = adc.data() + c * n;
      const double bl = daq::pulse::trimmedMean(x, daq::pulse::Region{0, 40}, 2.0 / lsb);
      pedSum[static_cast<std::size_t>(c)] += bl;
      std::vector<daq::pulse::Pulse> p;
      daq::pulse::findPulses(x, n, bl, lsb, 2 * lsb, 4, zones, p);
      for (const auto& q : p) {
        const int b = static_cast<int>(q.amplitudeMv / lsb);
        if (b >= 0 && b < nb) ++spec[static_cast<std::size_t>(c)][static_cast<std::size_t>(b)];
        if (q.amplitudeMv > 3.0) ++dark3[static_cast<std::size_t>(c)];
      }
    }
  }
  // Cabecera del primer segmento: run, HV y temperatura de la placa al inicio.
  nlohmann::json runInfo, hvState;
  {
    std::unique_ptr<TFile> f(TFile::Open(t.GetListOfFiles()->At(0)->GetTitle(), "READ"));
    if (auto* r = f->Get<TNamed>("run_info_json")) runInfo = nlohmann::json::parse(r->GetTitle());
    if (auto* h = f->Get<TNamed>("hvmon_state_at_start")) hvState = nlohmann::json::parse(h->GetTitle());
  }
  nlohmann::json out{{"schema", "daqv2.pe_calibration/1"},
                     {"version", version},
                     {"source_run", runInfo.value("campaign", "") + " run " + std::to_string(runInfo.value("run", 0))},
                     {"source_files", pattern},
                     {"events", ne},
                     {"first_read_utc_ns", t0},
                     {"last_read_utc_ns", t1},
                     {"method", "separación de picos de pulsos oscuros (daq/gain_fit.h), reglas de la fase 5; hod2019 sin picos: tasa > 3 mV"},
                     {"board_temperature_degC", hvState.contains("board_temperature_degC") ? hvState["board_temperature_degC"] : nlohmann::json()}};
  nlohmann::json hv = nlohmann::json::array();
  for (const auto& c : hvState.value("channels", nlohmann::json::array()))
    if (c.value("ch", -1) < 8) {
      const auto& v = c["values"];
      hv.push_back({{"ch", c["ch"]}, {"Pw", v.value("Pw", 0)}, {"V0Set_V", v.value("V0Set_V", 0.0)}, {"VMon_V", v.value("VMon_V", 0.0)},
                    {"IMon_A", v.value("IMon_A", 0.0)}, {"Status", v.value("Status", 0)}});
    }
  out["hv_at_start"] = hv;
  nlohmann::json chs = nlohmann::json::array();
  std::map<std::string, int> count;
  std::vector<double> accepted;
  printf("eventos %lld, exposición %.4f s por canal\n", ne, ne * exposure);
  for (int c = 0; c < 64; ++c) {
    const auto i = static_cast<std::size_t>(c);
    const double ped = (pedSum[i] / std::max<Long64_t>(1, ne)) * lsb;  // mV sobre 0 ADC
    nlohmann::json e{{"ch", c}, {"station", c < 32 ? "hod2018" : "hod2019"}, {"pedestal_mV", ped}};
    const double dr = dark3[i] / (ne * exposure);
    e["dark_rate_above_3mV_Hz"] = dr;
    if (c < 32) {
      double seed = 4.0;
      if (seedVals.count(c)) {
        auto v = seedVals[c];
        std::sort(v.begin(), v.end());
        seed = v[v.size() / 2];
      }
      const auto r = daq::gain::fitGain(spec[i], lsb, seed);
      e["seed_mV"] = seed;
      e["status"] = r.status;
      e["peaks"] = r.peaks;
      e["peak_positions_mV"] = r.positionsMv;
      if (r.gainMv > 0) e["G_mV"] = r.gainMv;
      e["provisional"] = r.status != "ok";
      ++count[r.status];
      if (r.status == "ok") accepted.push_back(r.gainMv);
      printf("ch%02d semilla %.3f -> %-12s G %.3f mV (%d picos) ped %.2f mV\n", c, seed, r.status.c_str(), r.gainMv, r.peaks, ped);
    } else {
      e["status"] = "dark_rate_only";
      e["provisional"] = true;
      printf("ch%02d hod2019 tasa > 3 mV %.1f Hz, ped %.2f mV\n", c, dr, ped);
    }
    chs.push_back(e);
  }
  out["channels"] = chs;
  out["status_count"] = count;
  std::sort(accepted.begin(), accepted.end());
  if (!accepted.empty())
    printf("hod2018: %zu canales ok, G mediana %.3f mV (%.3f-%.3f)\n", accepted.size(), accepted[accepted.size() / 2], accepted.front(),
           accepted.back());
  for (const auto& [k, v] : count) printf("  %s: %d\n", k.c_str(), v);
  std::ofstream(outJson) << out.dump(1) << "\n";
}
