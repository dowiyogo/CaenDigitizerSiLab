#!/usr/bin/python3 -I
# h19_fold_summary.py — seguimiento 2026-11-15 (15.1): resumen de h19_fold.json. Para la malla de Geant4 (diagonal T_X = T_Y = T) da el T en que las
# X·Y predichas valen las 3 observadas y en que dos estaciones vale 1 (interpolación lineal de ln N en T), para comparar con el T de r_X, r_Y.
# Uso: h19_fold_summary.py <h19_fold.json>
import json, math, sys

d = json.load(open(sys.argv[1]))
N = d["run_B_triggers"]


def t_for(TT, vals, target):
    for k in range(len(TT) - 1):
        a, b = vals[k], vals[k + 1]
        if a > 0 and b > 0 and (a - target) * (b - target) <= 0:
            return TT[k] + math.log(a / target) / math.log(a / b) * (TT[k + 1] - TT[k])
    return None


for c, o in d["per_conv"].items():
    for v in ("c", "lo", "hi"):
        x = o[v]
        TT = x["thresholds_pe"]
        nxy = [N * p for p in x["PXY_diag"]]
        n2 = [N * p for p in x["PXY_h18_diag"]]
        print(f"{c} [{v}]: T con N_XY = 3: {t_for(TT, nxy, 3):.2f} p.e.; T con N_dos_estaciones = 1: {t_for(TT, n2, 1):.2f} p.e.; "
              f"T(r_X) {x['T_X_pe']:.2f}, T(r_Y) {x['T_Y_pe']:.2f} (extrapolados); P_X(10) {x['PX'][-1]:.5f} = {x['PX'][-1] / 0.0144:.2f} r_X, "
              f"P_Y(10) {x['PY'][-1]:.5f} = {x['PY'][-1] / 0.0093:.2f} r_Y")
for c, o in d["per_conv"].items():
    for v in ("pc", "plo", "phi"):
        x = o[v]
        print(f"{c} [{v}] Poisson: T_X {x['T_X_pe']:.2f} p.e. G_X {x['gain_X_mV_per_pe']:.3f} mV/p.e.; T_Y {x['T_Y_pe']:.2f} G_Y {x['gain_Y_mV_per_pe']:.3f}; "
              f"N_XY {x['N_XY_pred']:.3g}; dos estaciones {x['N_two_station_pred_geom']:.3g}")
g = [o[v][k] for o in d["per_conv"].values() for v in ("c", "lo", "hi", "pc", "plo", "phi") for k in ("gain_X_mV_per_pe", "gain_Y_mV_per_pe")]
print(f"ganancia implícita, todas las variantes y convenciones: {min(g):.3f} - {max(g):.3f} mV/p.e. (cola oscura: 0.5 - 2.5)")
for v in ("c", "pc"):
    x = d["per_conv"]["A+/A+"][v]
    print(f"validación Poisson contra Geant4 (A+/A+): " if v == "pc" else "Geant4 (A+/A+): ",
          " ".join(f"T{T}: PX {x['PX'][k]:.5f}" for k, T in enumerate(x["thresholds_pe"]) if T <= 10))
