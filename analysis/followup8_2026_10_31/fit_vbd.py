#!/usr/bin/python3 -I
# fit_vbd.py — seguimiento 2026-10-31 (8), punto 1: G de hod2018 por toma y voltaje de ruptura.
# G por canal y toma con el método de c1/PROC (gain_fit.h, reglas de la fase 5: espectro fuera de tiempo en bins de 1 LSB, picos 5-26 mV, ≥ 3
# picos), solo ajustes aceptados (ok/changed):
#   nominal: runs SW scan22_dV0, scan23_dV0 y c1_source (gain_spectra.C, opciones c1, semilla c1) y física runs 4, 7, 8, 9, 10 (charact_gain.C);
#   −0.5 V: scan23_dVm0p5 con semilla 0.55 × c1 (con semilla c1 el ajuste no encuentra picos: la ventana de la tendencia y la separación
#           mínima escalan con la semilla); −1.0 y −1.5 V: ningún canal (semillas c1 y 0.15 × c1).
# Voltaje de cada canal = VMon medio de su canal SY5527 (grupo HV de r4) en la toma (runs_hv.json).
# Si un canal tiene G en ≥ 3 voltajes distintos: ajuste lineal G = k (V − V_bd). Si tiene G solo en 2 (nominal y −0.5 V): estimación de dos
# puntos [inferido, sin grados de libertad]: k = (G_nom − G_−0.5)/(V_nom − V_−0.5), V_bd = V_nom − G_nom/k, OV = G_nom/k; error de OV solo por
# la dispersión de G_nom entre las tomas SW nominales. También el número de pulsos fuera de tiempo por evento en cada toma (sustento).
# Salida: <out>/vbd.json, vbd_channels.csv. Uso: fit_vbd.py <out dir>
import csv, json, math, os, statistics as st, sys

import numpy as np

OUT = sys.argv[1]
G_ = OUT + "/gain"
R = {e["label"]: e for e in json.load(open(OUT + "/runs_hv.json"))}
MAP = {c["dig_ch"]: c for c in json.load(open("daq/config/mapping/hodo_mapping_2026_10_r4.json"))["channels"]}
CH = [c for c in range(32) if c != 3]


def spectra_fits(name, take):
    out = {}
    for r in csv.DictReader(open(f"{G_}/S_{name}_fits.csv")):
        if r["options"] == "c1" and r["status"] in ("ok", "changed"):
            out[int(r["ch"])] = float(r["G_mV"])
    return out


def charact(name):
    return {int(r["ch"]): float(r["G_mV"]) for r in csv.DictReader(open(f"{G_}/G_{name}.csv")) if r["status"] in ("ok", "changed")}


TAKES = [("scan22_dV0", spectra_fits("scan22_dV0", "scan22_dV0"), "SW"), ("scan23_dV0", spectra_fits("scan23_dV0", "scan23_dV0"), "SW"),
         ("c1_source", spectra_fits("c1_source", "c1_source"), "SW"), ("scan23_dVm0p5", spectra_fits("scan23_dVm0p5_seed0p55", "scan23_dVm0p5"), "SW"),
         ("scan23_dVm1p0", spectra_fits("scan23_dVm1p0_seed0p15", "scan23_dVm1p0"), "SW"), ("scan23_dVm1p5", spectra_fits("scan23_dVm1p5", "scan23_dVm1p5"), "SW")] + \
        [(f"phys_run{r}", charact(f"phys_run{r}"), "ext") for r in (4, 7, 8, 9, 10)]
NOMSW = ("scan22_dV0", "scan23_dV0", "c1_source")
# pulsos fuera de tiempo por evento (> 1.5 mV y > 3 mV) en los runs SW, del espectro guardado
CNT = {}
for name, f in (("scan22_dV0", "scan22_dV0"), ("scan23_dV0", "scan23_dV0"), ("c1_source", "c1_source"), ("scan23_dVm0p5", "scan23_dVm0p5"),
                ("scan23_dVm1p0", "scan23_dVm1p0"), ("scan23_dVm1p5", "scan23_dVm1p5")):
    d = np.loadtxt(f"{G_}/S_{f}_spec.csv", delimiter=",", skiprows=1)
    n = R[name]["n_events"]
    CNT[name] = {c: {"gt1p5mV": float(d[(d[:, 0] == c) & (d[:, 1] >= 1.5)][:, 2].sum() / n), "gt3mV": float(d[(d[:, 0] == c) & (d[:, 1] >= 3.0)][:, 2].sum() / n)}
                 for c in CH}
rows = []
for c in CH:
    hv = MAP[c]["hv_group"] - 1
    r = {"ch": c, "hv_group": hv + 1, "plane": MAP[c]["plane"]}
    pts = []
    for name, g, trig in TAKES:
        v = R[name]["channels"][str(hv)]["VMon_mean_V" if False else "VMon_V"]["mean"]
        r[f"G_{name}"] = g.get(c)
        r[f"V_{name}"] = v
        if c in g:
            pts.append((round(R[name]["channels"][str(hv)]["V0Set_V"]["mean"], 2), v, g[c], name))
    r["distinct_V0Set_with_G"] = sorted({p[0] for p in pts})
    r["n_distinct_V"] = len(r["distinct_V0Set_with_G"])
    nom = [r[f"G_{n}"] for n in NOMSW if r[f"G_{n}"] is not None]
    r["G_nom_mean"] = st.mean(nom) if nom else None
    r["G_nom_sd"] = st.stdev(nom) if len(nom) > 1 else None
    if r["n_distinct_V"] >= 3:
        x = np.array([p[1] for p in pts]); y = np.array([p[2] for p in pts])
        A = np.vstack([x, np.ones_like(x)]).T
        (k, b), res_, *_ = np.linalg.lstsq(A, y, rcond=None)
        r["fit"] = {"k_mV_per_V": float(k), "V_bd_V": float(-b / k)}
    elif r["G_scan23_dVm0p5"] is not None and nom:
        vn = st.mean(r[f"V_{n}"] for n in NOMSW if r[f"G_{n}"] is not None)
        v1, g1, g0 = r["V_scan23_dVm0p5"], r["G_scan23_dVm0p5"], r["G_nom_mean"]
        k = (g0 - g1) / (vn - v1)
        ov = g0 / k
        eov = (r["G_nom_sd"] or 0) * abs((vn - v1) * (-g1) / (g0 - g1) ** 2)  # d(OV)/d(G0) · sd(G0)
        r["two_point"] = {"V_nom_VMon_V": vn, "V_m0p5_VMon_V": v1, "dV_V": vn - v1, "G_nom_mV": g0, "G_m0p5_mV": g1, "k_mV_per_V": k, "V_bd_V": vn - ov,
                          "OV_V": ov, "OV_err_from_Gnom_sd_V": eov, "OV_over_V0Set": ov + (R["scan23_dV0"]["channels"][str(hv)]["V0Set_V"]["mean"] - vn)}
    r["pulses_per_event"] = {n: CNT[n][c] for n in CNT}
    rows.append(r)
summary = {"n_channels": len(rows), "n_with_3_voltages": sum(r["n_distinct_V"] >= 3 for r in rows),
           "n_by_distinct_V": {str(k): sum(r["n_distinct_V"] == k for r in rows) for k in range(0, 5)},
           "two_point": [{"ch": r["ch"], "hv_group": r["hv_group"], **r["two_point"]} for r in rows if "two_point" in r]}
tp = summary["two_point"]
if tp:
    summary["two_point_median"] = {k: st.median(t[k] for t in tp) for k in ("k_mV_per_V", "OV_V", "dV_V")}
json.dump({"summary": summary, "channels": rows, "pulses_note": "pulsos fuera de tiempo por evento en [0,49) ∪ [102,624) de los runs SW (571 muestras de 16 ns)"},
          open(OUT + "/vbd.json", "w"), indent=1, default=float)
with open(OUT + "/vbd_channels.csv", "w") as f:
    names = [t[0] for t in TAKES]
    f.write("ch,plane,hv_group," + ",".join(f"G_{n}_mV" for n in names) + ",n_distinct_V_with_G,G_nom_mean_mV,G_nom_sd_mV,two_point_k_mV_per_V,two_point_OV_V," +
            ",".join(f"pulses_gt1p5mV_per_ev_{n}" for n in CNT) + "\n")
    for r in rows:
        fm = lambda x: "" if x is None else f"{x:.4f}"
        f.write(f"{r['ch']},{r['plane']},{r['hv_group']}," + ",".join(fm(r[f'G_{n}']) for n in names) + f",{r['n_distinct_V']},{fm(r['G_nom_mean'])},{fm(r['G_nom_sd'])},"
                f"{fm(r.get('two_point', {}).get('k_mV_per_V'))},{fm(r.get('two_point', {}).get('OV_V'))}," + ",".join(f"{r['pulses_per_event'][n]['gt1p5mV']:.4f}" for n in CNT) + "\n")
print("voltajes distintos con G medida, por canal:", summary["n_by_distinct_V"], "| con ≥ 3:", summary["n_with_3_voltages"])
for t in tp:
    print(f"ch{t['ch']:02d} grupo {t['hv_group']}: G_nom {t['G_nom_mV']:.3f} mV, G(−0.5) {t['G_m0p5_mV']:.3f} mV, ΔVMon {t['dV_V']:.3f} V -> k {t['k_mV_per_V']:.2f} mV/V, "
          f"V_bd {t['V_bd_V']:.2f} V (VMon), OV {t['OV_V']:.2f} ± {t['OV_err_from_Gnom_sd_V']:.2f} V")
if tp:
    print("mediana dos puntos:", {k: round(v, 3) for k, v in summary["two_point_median"].items()})
print("pulsos fuera de tiempo > 1.5 mV por evento, mediana de los 31 canales:",
      {n: round(st.median(r["pulses_per_event"][n]["gt1p5mV"] for r in rows), 3) for n in CNT})
print("pulsos > 3 mV por evento, mediana:", {n: round(st.median(r["pulses_per_event"][n]["gt3mV"] for r in rows), 4) for n in CNT})
