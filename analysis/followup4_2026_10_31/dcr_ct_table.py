#!/usr/bin/python3 -I
# dcr_ct_table.py — seguimiento 2026-10-31 (4), punto 1: tasa oscura y crosstalk por canal desde las salidas de dark_stage2.C (XPLUS-4b).
# dark_stage2.C cuenta TODOS los pulsos (buscador del DAQ, mínimos locales, separación 3 muestras = 48 ns de tiempo muerto, piso 1 mV) en la
# región pre-trigger de cada registro (octubre: muestras [20, 49) = 464 ns; mayo: [20, 258) = 3.81 µs), con el pedestal local de la región, y
# descuenta el tiempo muerto de cada pulso en la exposición: la tasa es directa, sin corrección de ocupación.
#   DCR = R(> 0.5 p.e.); p_ct = N(> 1.5 p.e.) / N(> 0.5 p.e.) (misma exposición), error binomial; factor 1/(1 − p_ct).
#   p_ct incluye el apilamiento accidental de dos pulsos oscuros sin resolver: p_pile ≈ 1 − exp(−DCR·τ) con τ = 32 ns (2 muestras, supuesto);
#   se da también p_ct corregido = (p_ct − p_pile)/(1 − p_pile) (inferido, depende de τ).
#   hod2019: sin G (c1 no tiene picos para hod2019): tasas sobre 1 y 3 mV; p_ct no se calcula.
# Afterpulsing, canal por canal: n(Δ) de pulsos consecutivos > 0.5 p.e. del mismo canal en la región de L muestras. Modelo B (central): pendiente
# libre de ln[n(Δ)/(L − Δ)] ajustada por canal en la ventana larga y extrapolada; modelo A (control): Poisson con la
# tasa del canal r = R(0.5)·16 ns por muestra y el tiempo muerto de 3 muestras: m(Δ) ∝ (L − Δ)·exp(−r(Δ − 3)), normalizado a lo observado en
# Δ ∈ [40, min(150, L − 20)] (la razón de las formas en esa ventana se da como control). Exceso en Δ ∈ [3, 12] muestras (48-192 ns; < 48 ns no se
# resuelve) por pulso > 0.5 p.e. del canal = P(afterpulso en 48-192 ns). El pulso anterior es cualquier pulso > 0.5 p.e. (casi siempre 1 p.e.).
# (El ajuste con todos los canales sumados de stage2_summary.py mezcla exponenciales de tasas distintas y subestima el fondo a Δ corto.)
# Uso: dcr_ct_table.py <salida dir> <tag>=<dir>/<prefijo>:<G csv o c1>:<etiqueta> ...
import csv, json, math, os, statistics as st, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

OUT = sys.argv[1]
TAU = 32e-9
TS = 16e-9
H18 = [c for c in range(32) if c != 3]


def gains(spec):
    if spec == "c1":
        g = mc.g_c1()
        return np.where(np.isfinite(g), g, np.nanmedian(g)), [c for c in range(32) if not np.isfinite(g[c])]
    g = np.array([float(x) for x in open(spec).read().strip().split(",")])
    return g, []


def afterpulse(dtfile, L, rates, counts):
    """rates: {ch: R(0.5) en cps}; counts: {ch: N(> 0.5 p.e.)}. Devuelve por canal y la mediana."""
    n = {}
    for r in csv.DictReader(open(dtfile)):
        c = int(r["ch"])
        if c in rates:
            n.setdefault(c, {})[int(r["dt_samples"])] = int(r["count"])
    lo, hi = 40, min(150, L - 20)
    if hi - lo < 10:
        return {"L_samples": L, "note": "región demasiado corta para normalizar el modelo"}
    per = {}
    for c, nc in n.items():
        # modelo A: tasa del canal R(0.5) fija; modelo B: pendiente libre ajustada en la ventana larga (ln[n/(L − Δ)] lineal, ponderado)
        xs = np.array([d for d in range(lo, hi + 1) if nc.get(d, 0) > 0])
        if len(xs) < 8 or sum(nc.get(d, 0) for d in xs) < 50:
            continue
        ys = np.array([math.log(nc[d] / (L - d)) for d in xs])
        w = np.array([nc[d] for d in xs], float)
        b, a = np.polyfit(xs, ys, 1, w=np.sqrt(w))
        mB = lambda d: math.exp(a + b * d) * (L - d)
        r = rates[c] * TS
        mA = lambda d: (L - d) * math.exp(-r * (d - 3))
        kA = sum(nc.get(d, 0) for d in range(lo, hi + 1)) / sum(mA(d) for d in range(lo, hi + 1))
        obs = sum(nc.get(d, 0) for d in range(3, 13))
        eB = sum(mB(d) for d in range(3, 13))
        eA = kA * sum(mA(d) for d in range(3, 13))
        per[c] = {"observed": obs, "expected_free_slope": eB, "expected_R05": eA, "P_ap": (obs - eB) / counts[c], "P_ap_R05_model": (obs - eA) / counts[c],
                  "P_ap_stat": math.sqrt(obs) / counts[c], "slope_rate_cps": -b / TS, "slope_over_R05": (-b / TS) / rates[c]}
    vals = [v["P_ap"] for v in per.values()]
    valsA = [v["P_ap_R05_model"] for v in per.values()]
    ratio = [v["slope_over_R05"] for v in per.values()]
    return {"L_samples": L, "norm_range_samples": [lo, hi], "window_samples": [3, 12], "window_ns": [48, 192], "per_channel": per,
            "P_ap_median": st.median(vals), "P_ap_median_err": 1.2533 * st.pstdev(vals) / math.sqrt(len(vals)), "P_ap_min_max": [min(vals), max(vals)],
            "P_ap_R05_model_median": st.median(valsA), "slope_over_R05_median": st.median(ratio), "n_channels": len(vals)}


res = {}
for spec in sys.argv[2:]:
    tag, rest = spec.split("=", 1)
    prefix, gspec, label = rest.split(":", 2)
    G, filled = gains(gspec)
    R, E, C = {}, {}, {}
    for r in csv.DictReader(open(prefix + "_stair.csv")):
        c, T = int(r["ch"]), float(r["threshold"])
        R.setdefault(c, {})[T] = (float(r["rate_cps"]), float(r["rate_err_cps"]))
        C.setdefault(c, {})[T] = int(r["counts"])
        E[c] = float(r["exposure_s"])
    ex = list(csv.DictReader(open(prefix + "_expo.csv")))[0]
    rows = []
    for c in sorted(R):
        st_ = "hod2018" if c < 32 else "hod2019"
        if c < 32:
            if c == 3 or R[c].get(0.5, (0,))[0] <= 0:
                continue
            n05, n15 = C[c][0.5], C[c][1.5]
            p = n15 / n05
            ep = math.sqrt(p * (1 - p) / n05)
            dcr = R[c][0.5][0]
            pp = 1 - math.exp(-dcr * TAU)
            pc = (p - pp) / (1 - pp)
            rows.append({"ch": c, "station": st_, "plane": "X" if c < 16 else "Y", "G_mV": float(G[c]), "G_is_median": c in filled,
                         "DCR_kHz": dcr / 1e3, "DCR_err_kHz": R[c][0.5][1] / 1e3, "N_gt0p5": n05, "N_gt1p5": n15, "p_ct": p, "p_ct_err": ep,
                         "factor_1_over_1mp": 1 / (1 - p), "factor_err": ep / (1 - p) ** 2, "p_pile_tau32ns": pp, "p_ct_pile_corrected": pc,
                         "factor_pile_corrected": 1 / (1 - pc), "R2p5_over_R1p5": R[c][2.5][0] / R[c][1.5][0] if R[c][1.5][0] else None})
        else:
            rows.append({"ch": c, "station": st_, "R_gt1mV_kHz": R[c].get(1.0, (0, 0))[0] / 1e3, "R_gt1mV_err_kHz": R[c].get(1.0, (0, 0))[1] / 1e3,
                         "R_gt3mV_kHz": R[c].get(3.0, (0, 0))[0] / 1e3, "R_gt3mV_err_kHz": R[c].get(3.0, (0, 0))[1] / 1e3})
    h18 = [r for r in rows if r["station"] == "hod2018"]
    h19 = [r for r in rows if r["station"] == "hod2019"]
    med = lambda v: (st.median(v), 1.2533 * st.pstdev(v) / math.sqrt(len(v)))
    L = 29 if ex["region"] == "pre" and int(ex["s0"]) == 49 else (238 if ex["region"] == "pre" else (144 - (int(ex["s0"]) + 9) if int(ex["s0"]) < 100 else None))
    o = {"label": label, "events": int(ex["events"]), "region": ex["region"], "s0": int(ex["s0"]), "region_len_samples": L, "rows": rows,
         "median_hod2018": {k: med([r[k] for r in h18]) for k in ("DCR_kHz", "p_ct", "factor_1_over_1mp", "p_ct_pile_corrected", "factor_pile_corrected")},
         "median_hod2019": {k: med([r[k] for r in h19]) for k in ("R_gt1mV_kHz", "R_gt3mV_kHz")} if h19 else None}
    if L and L > 60:
        o["afterpulse"] = afterpulse(prefix + "_dt.csv", L, {r["ch"]: r["DCR_kHz"] * 1e3 for r in h18}, {r["ch"]: r["N_gt0p5"] for r in h18})
    res[tag] = o
    with open(f"{OUT}/table_{tag}.csv", "w") as f:
        f.write(f"# {label}; {o['events']} triggers; región {ex['region']} ({L} muestras de 16 ns por registro)\n")
        f.write("ch,station,plane,G_mV,DCR_kHz,DCR_err_kHz,p_ct,p_ct_err,factor_1_over_1mp,factor_err,p_ct_pile_corrected,factor_pile_corrected,R_gt1mV_kHz,R_gt3mV_kHz\n")
        for r in rows:
            g = lambda k, fmt="{:.4f}": "" if r.get(k) is None else fmt.format(r[k])
            f.write(",".join([str(r["ch"]), r["station"], r.get("plane", ""), g("G_mV"), g("DCR_kHz", "{:.1f}"), g("DCR_err_kHz", "{:.1f}"), g("p_ct"), g("p_ct_err"),
                              g("factor_1_over_1mp"), g("factor_err"), g("p_ct_pile_corrected"), g("factor_pile_corrected"), g("R_gt1mV_kHz", "{:.2f}"),
                              g("R_gt3mV_kHz", "{:.3f}")]) + "\n")
    m = o["median_hod2018"]
    print(f"{tag} ({label}): ev {o['events']} región {o['region']} L {L}")
    print(f"  hod2018 mediana: DCR {m['DCR_kHz'][0]:.0f} ± {m['DCR_kHz'][1]:.0f} kHz, p_ct {m['p_ct'][0]:.4f} ± {m['p_ct'][1]:.4f}, 1/(1−p_ct) {m['factor_1_over_1mp'][0]:.4f}, "
          f"p_ct corr {m['p_ct_pile_corrected'][0]:.4f}")
    if o["median_hod2019"]:
        print(f"  hod2019 mediana: R>1 mV {o['median_hod2019']['R_gt1mV_kHz'][0]:.2f} kHz, R>3 mV {o['median_hod2019']['R_gt3mV_kHz'][0]:.3f} kHz")
    if "afterpulse" in o:
        a = o["afterpulse"]
        print("  afterpulse:", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in a.items() if k != "per_channel"})
json.dump(res, open(f"{OUT}/dcr_ct.json", "w"), indent=1)
