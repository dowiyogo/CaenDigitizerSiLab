#ifndef HODO_CONSTS_H
#define HODO_CONSTS_H

// Sampling period: VX1740 @ 62.5 MS/s = 16.0 ns step
constexpr double kDt_ns = 16.0;

// Physical cosmic muon time window (both stations)
#ifndef TAGPROBE_MUON_OVERRIDE
constexpr double kMuonTlo = 4000.0;
constexpr double kMuonThi = 4400.0;
#endif

// Baseline selection for Q integration (Decision Rule D1 from Step C):
// hod2018: RMS_loc/RMS_glob = 100.68/180.03 = 0.5592 < 0.7 (Rule D1 fired -> use ped_loc)
constexpr bool kQLocalPed18 = true;
// hod2019: RMS_loc/RMS_glob = 155.76/129.87 = 1.1993 >= 0.7 (Rule D1 not fired -> keep ped_glob)
constexpr bool kQLocalPed19 = false;

// Integration window parameters (Decision Rule D2 from Step B containment n95):
// kNpost18: n95 = 12 from pulse_shape.C for ped_loc (total [-3, 19] = 13761.924 mV*ns, n0=19)
constexpr int kNpre = 3;
constexpr int kNpost18 = 12;
// kNpost19: n95 = 22 from pulse_shape.C for ped_glob (total [-3, 40] = 10951.524 mV*ns, n0=40)
constexpr int kNpost19 = 22;

// Diagnostic pulse height thresholds (mV from pedestal)
constexpr double kPh18_mV = 50.0;
constexpr double kPh19_mV = 35.0;

// Q cut for hod2019 (Rule 1d: hod2019X v8 valley invalid, 5 * sigma_Q_noise = 5 * 137.4 mV·ns)
constexpr double kQcut19 = 687.0;

#endif // HODO_CONSTS_H
