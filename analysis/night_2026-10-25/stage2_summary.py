#!/usr/bin/python3 -I
# stage2_summary.py — XPLUS-4b etapa 2: tablas desde las salidas de dark_stage2.C (_stair, _dt, _big.csv.gz, _expo).
#  2.1-2.3 por canal de hod2018 (sin ch3) y mediana: R(T) en cps (tiempo muerto descontado en la exposición), DCR = R(0.5),
#          P_ct = R(1.5)/R(0.5), R(2.5)/R(1.5), fracción > 6 p.e. = R(6)/R(0.5).
#  2.4 intervalos: n(Δ) de pulsos consecutivos > 0.5 p.e. del mismo canal (todos los canales de hod2018 sumados). En una región de
#      longitud L, para un proceso de Poisson n(Δ) ∝ (L − Δ)·exp(−bΔ): se ajusta ln[n(Δ)/(L − Δ)] con una recta en Δ ∈ [40, 150] muestras
#      (0.64-2.4 µs) y se extrapola a Δ ∈ [4, 31] (64-496 ns). Exceso = (observado − extrapolado)/pulsos > 0.5 p.e. Además b contra
#      la tasa R(0.5) (para Poisson, b ≈ R(0.5)·16 ns).
#  2.7 pulsos > 6 p.e.: n_same (otros canales de la misma coordenada > 1 p.e. a ±1 muestra), n_orth, n_h19 (> 3 mV), y Δ al pulso
#      previo. Esperado por azar con las tasas medidas: p_c = 1 − exp(−R_c·3·16 ns) por canal (Poisson-binomial); para Δ_prev, la
#      fracción con pulso previo en ≤ 62 muestras (1 µs) para Poisson con R(0.5).
# Uso: stage2_summary.py <dir> <tag1> [tag2 ...]   (tag = prefijo de dark_stage2.C dentro de <dir>)
import csv, gzip, json, math, os, statistics as st, sys

import numpy as np

D, TAGS = sys.argv[1], sys.argv[2:]
H18 = [c for c in range(32) if c != 3]
TS = 16e-9


def stair(tag):
    R, E = {}, {}
    for r in csv.DictReader(open(f"{D}/{tag}_stair.csv")):
        c = int(r["ch"])
        R.setdefault(c, {})[float(r["threshold"])] = (float(r["rate_cps"]), float(r["rate_err_cps"]), int(r["counts"]))
        E[c] = float(r["exposure_s"])
    return R, E


out = {}
for tag in TAGS:
    R, E = stair(tag)
    ex = list(csv.DictReader(open(f"{D}/{tag}_expo.csv")))[0]
    nev, region, s0 = int(ex["events"]), ex["region"], int(ex["s0"])
    rows = {}
    for c in H18:
        r = R[c]
        if r[0.5][0] <= 0:
            continue
        rows[c] = {"DCR_cps": r[0.5][0], "DCR_err": r[0.5][1], "P_ct": r[1.5][0] / r[0.5][0], "R2p5_R1p5": r[2.5][0] / r[1.5][0] if r[1.5][0] else None,
                   "frac_gt6": r[6.0][0] / r[0.5][0], "R": {str(k): v[0] for k, v in r.items()}}
    med = lambda k: st.median([v[k] for v in rows.values() if v[k] is not None])
    o = {"events": nev, "region": region, "s0": s0, "n_channels": len(rows), "DCR_median_cps": med("DCR_cps"), "P_ct_median": med("P_ct"),
         "R2p5_R1p5_median": med("R2p5_R1p5"), "frac_gt6_median": med("frac_gt6"),
         "staircase_median_cps": {str(T): st.median([v["R"][str(T)] for v in rows.values()]) for T in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)},
         "per_channel": rows}
    o["DCR_median_err_cps"] = 1.2533 * st.pstdev([v["DCR_cps"] for v in rows.values()]) / math.sqrt(len(rows))
    o["P_ct_median_err"] = 1.2533 * st.pstdev([v["P_ct"] for v in rows.values()]) / math.sqrt(len(rows))
    hod19 = [R[c][1.0][0] for c in range(32, 64) if c in R]
    o["h19_R_gt1mV_median_cps"] = st.median(hod19) if hod19 else None
    o["h19_R_gt3mV_median_cps"] = st.median([R[c][3.0][0] for c in range(32, 64) if c in R])
    # 2.4 intervalos
    n = {}
    for r in csv.DictReader(open(f"{D}/{tag}_dt.csv")):
        if int(r["ch"]) in H18:
            n[int(r["dt_samples"])] = n.get(int(r["dt_samples"]), 0) + int(r["count"])
    L = (s0 - 20) if region == "pre" else ((rlen := None) or None)
    if region == "all":
        L = 604
    elif region == "post":
        L = 144 - (s0 + 9) if s0 < 100 else 624 - (s0 + 9)
    fit_lo, fit_hi = (40, min(150, L - 20))
    xs = np.array([d for d in range(fit_lo, fit_hi + 1) if n.get(d, 0) > 0])
    if len(xs) > 5:
        ys = np.array([math.log(n[d] / (L - d)) for d in xs])
        b, a = np.polyfit(xs, ys, 1)
        pred = lambda d: math.exp(a + b * d) * (L - d)
        short = range(4, 32)
        obs_s = sum(n.get(d, 0) for d in short)
        pred_s = sum(pred(d) for d in short)
        tot = sum(rr[0.5][2] for c, rr in R.items() if c in H18)
        o["intervals"] = {"L_samples": L, "fit_range_samples": [fit_lo, fit_hi], "slope_per_sample": -b, "slope_rate_cps": -b / TS,
                          "R0p5_median_cps": o["DCR_median_cps"], "short_range_samples": [4, 31], "observed_short": obs_s, "expected_short": pred_s,
                          "excess_over_pulses": (obs_s - pred_s) / tot if tot else None, "excess_rel_to_expected": obs_s / pred_s - 1,
                          "n_by_dt": {str(d): n.get(d, 0) for d in range(3, min(L, 200))}, "pred_by_dt": {str(d): pred(d) for d in range(3, min(L, 200))}}
    # 2.7 pulsos grandes
    big = f"{D}/{tag}_big.csv.gz"
    if os.path.exists(big):
        R1 = {c: (R[c][1.0][0] if c in R else 0.0) for c in range(64)}
        R3mv = {c: (R[c][3.0][0] if c in R else 0.0) for c in range(32, 64)}
        pc = lambda rate: 1 - math.exp(-rate * 3 * TS)
        cnt = {"total": 0, "coherent_h19ge2": 0, "same_ge1": 0, "same_ge3": 0, "orth_ge1": 0, "none": 0, "prev_le62": 0, "none_and_prev_le62": 0}
        exp_same, exp_orth, exp_h19ge2 = [], [], []
        ns_sum = no_sum = n19_sum = 0
        with gzip.open(big, "rt") as f:
            for r in csv.DictReader(f):
                c = int(r["ch"])
                ns, no, n19, dp = int(r["n_same_coord"]), int(r["n_orth_coord"]), int(r["n_h19"]), int(r["dt_prev_samples"])
                cnt["total"] += 1
                ns_sum += ns
                no_sum += no
                n19_sum += n19
                cnt["coherent_h19ge2"] += n19 >= 2
                cnt["same_ge1"] += ns >= 1
                cnt["same_ge3"] += ns >= 3
                cnt["orth_ge1"] += no >= 1
                cnt["none"] += ns == 0 and no == 0 and n19 == 0
                cnt["prev_le62"] += 0 < dp <= 62
                cnt["none_and_prev_le62"] += ns == 0 and no == 0 and n19 == 0 and 0 < dp <= 62
        # expectativas por azar (promedio sobre los canales de hod2018, ponderado por su R(6))
        w = {c: R[c][6.0][2] for c in H18}
        W = sum(w.values())
        e_same = sum(w[c] * sum(pc(R1[c2]) for c2 in H18 if c2 != c and (c2 < 16) == (c < 16)) for c in H18) / W
        e_orth = sum(w[c] * sum(pc(R1[c2]) for c2 in H18 if (c2 < 16) != (c < 16)) for c in H18) / W
        p19 = [pc(R3mv[c2]) for c2 in range(32, 64)]
        lam19 = sum(p19)
        e_h19ge2 = 1 - math.exp(-lam19) * (1 + lam19)
        e_prev = sum(w[c] * (1 - math.exp(-R[c][0.5][0] * 62 * TS)) for c in H18) / W
        t = cnt["total"]
        o["big_pulses"] = {"counts": cnt, "mean_n_same": ns_sum / t, "expected_n_same_chance": e_same, "mean_n_orth": no_sum / t,
                           "expected_n_orth_chance": e_orth, "mean_n_h19": n19_sum / t, "expected_n_h19_chance": lam19,
                           "frac_h19ge2": cnt["coherent_h19ge2"] / t, "expected_frac_h19ge2_chance": e_h19ge2,
                           "frac_prev_le62": cnt["prev_le62"] / t, "expected_frac_prev_le62_chance": e_prev, "frac_no_coinc": cnt["none"] / t}
    out[tag] = o
    print(f"{tag}: ev {nev} DCR {o['DCR_median_cps']:.4g} ± {o['DCR_median_err_cps']:.2g}  P_ct {o['P_ct_median']:.4f} ± {o['P_ct_median_err']:.4f}  "
          f"R2.5/R1.5 {o['R2p5_R1p5_median']:.4f}  >6 {o['frac_gt6_median']:.4f}")
    if "intervals" in o:
        i = o["intervals"]
        print(f"   intervalos: pendiente {i['slope_rate_cps']:.4g} cps (R0.5 {i['R0p5_median_cps']:.4g}); corto obs {i['observed_short']} esp {i['expected_short']:.1f} "
              f"exceso/pulsos {i['excess_over_pulses']:.4f} (rel {i['excess_rel_to_expected']:+.3f})")
    if "big_pulses" in o:
        b = o["big_pulses"]
        print(f"   >6 p.e.: n {b['counts']['total']}  n_same {b['mean_n_same']:.3f} (azar {b['expected_n_same_chance']:.3f})  n_orth {b['mean_n_orth']:.3f} "
              f"(azar {b['expected_n_orth_chance']:.3f})  n_h19 {b['mean_n_h19']:.4f} (azar {b['expected_n_h19_chance']:.4f})  h19≥2 {b['frac_h19ge2']:.4f} "
              f"(azar {b['expected_frac_h19ge2_chance']:.4f})  prev≤1µs {b['frac_prev_le62']:.3f} (azar {b['expected_frac_prev_le62_chance']:.3f})")
json.dump(out, open(f"{D}/stage2_summary.json", "w"), indent=1, default=float)
