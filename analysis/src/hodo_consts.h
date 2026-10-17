#ifndef HODO_CONSTS_H
#define HODO_CONSTS_H

constexpr double kDt_ns = 16.0;          // VX1740, 62.5 MS/s; equals data_dig.time step
#ifndef TAGPROBE_MUON_OVERRIDE
constexpr double kMuonTlo = 4000.0, kMuonThi = 4400.0;  // ns, both stations, all runs v1-v10 (muon_timing_scan.C)
#endif
constexpr int kNpre = 3;
constexpr int kNpost18 = 8;   // 95% of mean pulse integral, pulse_shape.C on v8
constexpr int kNpost19 = 22;  // 95% of mean pulse integral, pulse_shape.C on v8
constexpr double kPh18_mV = 50.0, kPh19_mV = 35.0;  // pulse height from pedestal

#endif // HODO_CONSTS_H

