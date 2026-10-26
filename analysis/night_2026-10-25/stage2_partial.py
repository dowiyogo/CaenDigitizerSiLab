#!/usr/bin/python3 -I
# stage2_partial.py — XPLUS-4 etapa 2 (parcial, sin leer crudos): tasa de pulsos oscuros contra umbral desde productos ya escritos.
#  A) Barrido de día 2026-10-23 (ΔV = 0, −0.5, −1.0, −1.5 V): espectros de amplitud de pulsos oscuros de selftrig_noise.C
#     (outputs/daq2-15/biasscan_2026_10_23_day/noise_step{0..3}.csv: buscador del DAQ, separación 4 muestras, piso 1.0 mV, todo el
#     registro de 624 muestras de runs de SW trigger; exposición por canal). R(>T) = cuentas sobre T / exposición.
#     Umbrales: (i) en mV fijos = T·G(0 V) del canal para todos los escalones (sin supuestos); (ii) en p.e. del escalón: G(0 V) y
#     G(−0.5 V) de la tabla "nuevo" de report_es.md de DAQ2-15 (medidas); sin G(−0.5 V): 0.59·G(0 V) (supuesto, mediana).
#  B) Runs de física de octubre (run 8, noche 2026-10-19): cuentas por hora y canal de charact_extract.C (outputs/
#     characterization_2026_10_23/work/<run>_hourly.csv: pedestal + control + post, 1.456 µs por evento, separación 4,
#     umbrales 0.5/1.5/3/5 p.e. de c1 en hod2018, 1/3 mV en hod2019).
#  Corrección de tiempo muerto (calculada, no paralizable): R_c = R/(1 − R_floor·τ), τ = 4 × 16 ns, R_floor = tasa sobre el piso.
# Salida: <out>/stage2_partial.json, tablas .md y figuras (matplotlib). Uso: stage2_partial.py <out>
import json, math, os, re, statistics, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
D = os.path.join(REPO, "outputs/daq2-15/biasscan_2026_10_23_day")
C = os.path.join(REPO, "outputs/characterization_2026_10_23/work")
TAU = 4 * 16e-9
TS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
STEPS = {0: 0.0, 1: -0.5, 2: -1.0, 3: -1.5}
H18 = list(range(32))

# G por canal a 0 V y −0.5 V: columnas "nuevo" (DAQ2-16) de la tabla de report_es.md de DAQ2-15
G0, G05 = {}, {}
for line in open(os.path.join(D, "report_es.md")):
    m = re.match(r"\|\s*(\d+)\s*\|\s*([^|]*)\|\s*([^|]*)\|\s*([^|]*)\|\s*([^|]*)\|", line)
    if m and line.count("|") >= 11:
        ch = int(m.group(1))
        num = lambda s: float(s.strip().rstrip("*")) if re.match(r"^\s*[0-9.]+\*?\s*$", s) else None
        star = lambda s: s.strip().endswith("*")
        g0 = num(m.group(3)) if not star(m.group(3)) else None
        g05 = num(m.group(5)) if not star(m.group(5)) else None
        if ch < 32:
            G0[ch], G05[ch] = g0, g05
ratio_meas = [G05[c] / G0[c] for c in H18 if G0.get(c) and G05.get(c)]
RATIO = statistics.median(ratio_meas)


def load_spec(step):
    expo, spec = {}, {}
    for line in open(os.path.join(D, f"noise_step{step}.csv")):
        if line.startswith("#exposure"):
            _, ch, s = line.split()
            expo[int(ch)] = float(s)
        elif line[0].isdigit():
            ch, lo, n = line.strip().split(",")
            spec.setdefault(int(ch), []).append((float(lo), int(n)))
    return expo, {c: (np.array([a for a, _ in v]), np.array([b for _, b in v])) for c, v in spec.items()}


def rate_above(spec, expo, ch, thr_mV):
    lo, n = spec.get(ch, (np.array([]), np.array([])))
    k = int(n[lo >= thr_mV - 1e-9].sum())
    return k, k / expo[ch] if expo.get(ch) else float("nan")


res = {"tau_s": TAU, "thresholds_pe": TS, "G0_mV": G0, "G_m0p5_mV": G05, "ratio_G_m0p5_over_G0_median": RATIO,
       "ratio_n_channels": len(ratio_meas), "dlnG_dV_per_V": math.log(1 / RATIO) / 0.5, "scan": {}, "physics": {}}
for step, dv in STEPS.items():
    expo, spec = load_spec(step)
    rows = {}
    for ch in H18:
        if not G0.get(ch) or ch not in expo:
            continue
        kf, rfloor = rate_above(spec, expo, ch, 1.0)
        corr = 1 / (1 - rfloor * TAU) if rfloor * TAU < 1 else float("nan")
        g_step = G0[ch] if step == 0 else (G05.get(ch) if step == 1 else None)
        g_lab = "medida" if (step == 0 or (step == 1 and G05.get(ch))) else None
        if step == 1 and not G05.get(ch):
            g_step, g_lab = RATIO * G0[ch], "supuesta 0.59·G(0)"
        r = {"exposure_s": expo[ch], "R_floor_1mV_cps": rfloor, "deadtime_corr": corr, "G_step_mV": g_step, "G_step_label": g_lab,
             "fixed_mV": {}, "pe_step": {}}
        for T in TS:
            k, rr = rate_above(spec, expo, ch, T * G0[ch])
            r["fixed_mV"][str(T)] = {"thr_mV": T * G0[ch], "counts": k, "R_cps": rr, "R_err_cps": math.sqrt(k) / expo[ch], "R_corr_cps": rr * corr}
            if g_step and T * g_step >= 1.0:
                k2, r2 = rate_above(spec, expo, ch, T * g_step)
                r["pe_step"][str(T)] = {"thr_mV": T * g_step, "counts": k2, "R_cps": r2, "R_err_cps": math.sqrt(k2) / expo[ch], "R_corr_cps": r2 * corr}
        rows[ch] = r
    res["scan"][str(dv)] = rows

for run in ("run008", "night"):
    rows, hours = {}, {}
    for line in open(os.path.join(C, f"{run}_hourly.csv")):
        if line.startswith("hour"):
            continue
        f = line.strip().split(",")
        h, ch, ev, ex = int(f[0]), int(f[1]), int(f[2]), float(f[3])
        cnt = [int(x) for x in f[4:8]]
        thr = [float(x) for x in f[8:12]]
        a = rows.setdefault(ch, {"exposure_s": 0.0, "cnt": [0, 0, 0, 0], "thr_mV": thr, "events": 0})
        a["exposure_s"] += ex
        a["events"] += ev
        a["cnt"] = [x + y for x, y in zip(a["cnt"], cnt)]
        if ch < 32 and ex > 0:
            hours.setdefault(h, {})[ch] = (cnt[0] / ex, ev)
    out = {}
    for ch, a in rows.items():
        R = [k / a["exposure_s"] if a["exposure_s"] else float("nan") for k in a["cnt"]]
        corr = 1 / (1 - R[0] * TAU) if ch < 32 and R[0] * TAU < 1 else None
        out[ch] = {"events": a["events"], "exposure_s": a["exposure_s"], "thr_mV": a["thr_mV"], "counts": a["cnt"], "R_cps": R,
                   "R_err_cps": [math.sqrt(k) / a["exposure_s"] for k in a["cnt"]], "deadtime_corr_from_R0": corr,
                   "P_ct_R1p5_over_R0p5": R[1] / R[0] if ch < 32 and R[0] else None,
                   "R3_over_R1p5": R[2] / R[1] if ch < 32 and R[1] else None}
    hourly = []
    for h in sorted(hours):
        v = [x[0] for x in hours[h].values()]
        hourly.append({"hour_start_utc_s": h, "n_ch": len(v), "events_ch0": hours[h].get(0, (0, 0))[1], "DCR_median_h18_cps": statistics.median(v)})
    res["physics"][run] = {"channels": out, "hourly_DCR_h18": hourly}

# resúmenes
def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return statistics.median(xs) if xs else None


summ = {"scan": {}, "physics": {}}
for dv, rows in res["scan"].items():
    s = {"n_ch": len(rows)}
    for T in ("0.5", "1.5", "2.5", "6.0"):
        s[f"R_fixed_mV_{T}pe_median_cps"] = med([r["fixed_mV"][T]["R_cps"] for r in rows.values()])
        s[f"R_pe_step_{T}pe_median_cps"] = med([r["pe_step"].get(T, {}).get("R_cps") for r in rows.values()])
    s["P_ct_pe_step_median"] = med([r["pe_step"]["1.5"]["R_cps"] / r["pe_step"]["0.5"]["R_cps"] for r in rows.values()
                                    if "0.5" in r["pe_step"] and r["pe_step"]["0.5"]["R_cps"]])
    s["R2p5_over_R1p5_pe_step_median"] = med([r["pe_step"]["2.5"]["R_cps"] / r["pe_step"]["1.5"]["R_cps"] for r in rows.values()
                                              if "1.5" in r["pe_step"] and r["pe_step"]["1.5"]["R_cps"]])
    s["frac_above_6pe_pe_step_median"] = med([r["pe_step"]["6.0"]["R_cps"] / r["pe_step"]["0.5"]["R_cps"] for r in rows.values()
                                              if "6.0" in r["pe_step"] and "0.5" in r["pe_step"] and r["pe_step"]["0.5"]["R_cps"]])
    s["deadtime_corr_median"] = med([r["deadtime_corr"] for r in rows.values()])
    s["R_floor_1mV_median_cps"] = med([r["R_floor_1mV_cps"] for r in rows.values()])
    summ["scan"][dv] = s
for run, d in res["physics"].items():
    ch = d["channels"]
    summ["physics"][run] = {"DCR_0p5pe_median_h18_cps": med([ch[c]["R_cps"][0] for c in ch if c < 32]),
                            "DCR_0p5pe_median_h18_corr_cps": med([ch[c]["R_cps"][0] * ch[c]["deadtime_corr_from_R0"] for c in ch if c < 32 and ch[c]["deadtime_corr_from_R0"]]),
                            "P_ct_median_h18": med([ch[c]["P_ct_R1p5_over_R0p5"] for c in ch if c < 32]),
                            "R3_over_R1p5_median_h18": med([ch[c]["R3_over_R1p5"] for c in ch if c < 32]),
                            "R_1mV_median_h19_cps": med([ch[c]["R_cps"][0] for c in ch if c >= 32]),
                            "R_3mV_median_h19_cps": med([ch[c]["R_cps"][1] for c in ch if c >= 32]),
                            "R_above_1mV_hod2018_vs_hod2019_note": "hod2018 en p.e. (c1), hod2019 en mV: no comparables sin G de hod2019",
                            "hourly_DCR_min_max_cps": [min(x["DCR_median_h18_cps"] for x in d["hourly_DCR_h18"]),
                                                       max(x["DCR_median_h18_cps"] for x in d["hourly_DCR_h18"])]}
res["summary"] = summ
json.dump(res, open(os.path.join(OUT, "stage2_partial.json"), "w"), indent=1, default=float)
print(json.dumps(summ, indent=1, default=float))
print("razón G(−0.5)/G(0) mediana", RATIO, "en", len(ratio_meas), "canales; dlnG/dV", res["dlnG_dV_per_V"])

# figuras
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    for dv, rows in res["scan"].items():
        y = [med([r["fixed_mV"][str(T)]["R_cps"] for r in rows.values()]) for T in TS]
        ax.semilogy(TS, y, "o-", ms=3, label=f"ΔV = {float(dv):+.1f} V")
    ax.set_xlabel("umbral (p.e. de G a 0 V, mV fijos)")
    ax.set_ylabel("tasa mediana hod2018 (cps)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"staircase_scan_h18.{ext}"), dpi=200)
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    for c in H18:
        ys = [res["scan"][dv].get(c, {}).get("fixed_mV", {}).get("0.5", {}).get("R_cps") for dv in res["scan"]]
        if all(y is not None for y in ys):
            ax.semilogy([float(dv) for dv in res["scan"]], ys, "-", color="0.7", lw=0.6)
    ax.semilogy([float(dv) for dv in res["scan"]], [summ["scan"][dv]["R_fixed_mV_0.5pe_median_cps"] for dv in res["scan"]], "ko-", ms=3, label="mediana")
    ax.set_xlabel("ΔV respecto del nominal (V)")
    ax.set_ylabel("R(> 0.5 p.e. a 0 V) (cps)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"dcr_vs_dV_h18.{ext}"), dpi=200)
    print("figuras escritas")
except Exception as e:  # sin matplotlib: solo tablas
    print("sin figuras:", e)
