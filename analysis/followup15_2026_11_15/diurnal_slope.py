#!/usr/bin/python3 -I
# diurnal_slope.py — seguimiento 2026-11-15 (15, punto 3b): pendiente k = d ln R(>T) / d ln G de la cola del ruido de hod2018 en mayo (run B, bins
# buenos de diurnal_tail.py), y la misma pendiente para el crosstalk p_ct.
#  - R(>T): tasa agregada de la estación por hora (Σ cuentas / Σ tiempo vivo, diurnal_tail_bins.csv), T = 0.5, 3, 5, 10, 15 p.e. medidos, umbral fijo
#    (fixedG) y umbral que sigue la G (trackedG); G = G/G_ref horaria de la estación (followup 12 B, columna G_rel_station).
#  - p_ct: mediana de la estación de N(>1.5)/N(>0.5) por hora de followup 12 B (diurnal_bins.csv; umbrales con la G de la hora), y la versión corregida
#    por apilamiento (pct_corr).
# Ajuste ponderado de mínimos cuadrados ln y = a + k ln(G/G_ref), pesos 1/σ²(ln y) (Poisson para R, error de la mediana para p_ct); como la dispersión
# entre horas supera la de Poisson, el error de k se escala por √(χ²/ndf) cuando es > 1 (se da también sin escalar).
# Expectativa si la cola escala como R0·P_Borel(N ≥ T | p) con p = 0.429: d ln R(>T)/d ln G = k(0.5) + [d ln P(N≥T|p)/d ln p − d ln P(N≥1|p)/d ln p]·k_pct
# (P(N ≥ 1) = 1); la derivada de Borel se toma numérica. Para p pequeño, P(N ≥ T) ∝ p^(T−1): k crece con T.
# Salida: <out>/diurnal_slope.json, diurnal_slope.csv. Uso: diurnal_slope.py <diurnal_tail out dir> <followup12 dir>
import csv, json, math, os, sys

import numpy as np

D, F12 = sys.argv[1], sys.argv[2]
TPE = ["0.5", "3", "5", "10", "15"]
P_MAY = 0.429


def wfit(x, y, s):
    x, y, s = map(np.asarray, (x, y, s))
    w = 1 / s ** 2
    X = np.vstack([np.ones_like(x), x]).T
    C = np.linalg.inv((X * w[:, None]).T @ X)
    p = C @ (X * w[:, None]).T @ y
    chi = float((((y - X @ p) ** 2) * w).sum())
    ndf = len(y) - 2
    sc = math.sqrt(max(1.0, chi / ndf))
    return {"n": int(len(y)), "k": float(p[1]), "k_err_scaled": float(math.sqrt(C[1, 1]) * sc), "k_err_stat": float(math.sqrt(C[1, 1])),
            "chi2": chi, "ndf": ndf, "err_scale": sc}


def borel_ge(T, p):
    T = int(math.ceil(T))
    return 1.0 - sum(math.exp(-p * n + (n - 1) * math.log(p * n) - math.lgamma(n + 1)) for n in range(1, T))


def dlnP_dlnp(T, p, h=1e-4):
    return (math.log(borel_ge(T, p * (1 + h))) - math.log(borel_ge(T, p * (1 - h)))) / (math.log(1 + h) - math.log(1 - h))


rows = [r for r in csv.DictReader(l for l in open(f"{D}/diurnal_tail_bins.csv") if not l.startswith("#"))
        if r["campaign"] == "may" and r["good_bin"] == "1"]
res = {"sample": "May run B from 5.25 h, good-time 1 h bins (diurnal_tail.py)", "p_borel": P_MAY, "rate": {}, "pct": {}}
out_rows = []
for var in ("fixedG", "trackedG"):
    res["rate"][var] = {}
    for T in TPE:
        rr = [r for r in rows if r["variant"] == var and r["T_pe"] == T]
        x = [math.log(float(r["G_rel_station"])) for r in rr]
        y = [math.log(float(r["pooled_rate_cps"])) for r in rr]
        s = [float(r["pooled_rate_err_cps"]) / float(r["pooled_rate_cps"]) for r in rr]
        f = wfit(x, y, s)
        res["rate"][var][T] = f
        out_rows.append(["R(>T)", var, T, f["n"], f"{f['k']:.4f}", f"{f['k_err_scaled']:.4f}", f"{f['k_err_stat']:.4f}", f"{f['chi2']:.1f}", f["ndf"]])
pb = [r for r in csv.DictReader(l for l in open(f"{F12}/B_diurnal/diurnal_bins.csv") if not l.startswith("#")) if r["campaign"] == "may" and r["good_bin"] == "1"]
for key in ("pct", "pct_corr"):
    x = [math.log(float(r["G_rel_median"])) for r in pb]
    y = [math.log(float(r[f"{key}_median"])) for r in pb]
    s = [float(r[f"{key}_median_err"]) / float(r[f"{key}_median"]) for r in pb]
    f = wfit(x, y, s)
    res["pct"][key] = f
    out_rows.append([key, "hourly G thresholds", "1.5/0.5", f["n"], f"{f['k']:.4f}", f"{f['k_err_scaled']:.4f}", f"{f['k_err_stat']:.4f}", f"{f['chi2']:.1f}", f["ndf"]])
# expectativa de Borel a partir de k(0.5) y k_pct
kp, kpe = res["pct"]["pct"]["k"], res["pct"]["pct"]["k_err_scaled"]
res["borel_expectation"] = {}
for var in ("fixedG", "trackedG"):
    k05, k05e = res["rate"][var]["0.5"]["k"], res["rate"][var]["0.5"]["k_err_scaled"]
    e = {}
    for T in TPE[1:]:
        d = dlnP_dlnp(float(T), P_MAY)
        kexp = k05 + d * kp
        e[T] = {"dlnP_dlnp": d, "k_expected": kexp, "k_expected_err": math.hypot(k05e, d * kpe), "k_measured": res["rate"][var][T]["k"],
                "k_measured_err": res["rate"][var][T]["k_err_scaled"]}
        out_rows.append(["Borel expectation", var, T, "", f"{kexp:.4f}", f"{math.hypot(k05e, d * kpe):.4f}", "", "", f"dlnP/dlnp={d:.3f}"])
    res["borel_expectation"][var] = e
with open(f"{D}/diurnal_slope.csv", "w") as f:
    f.write("# diurnal_slope.py (seguimiento 15, 2026-11-15): k = d ln y / d ln(G/G_ref), weighted LS over May run B good-time 1 h bins; k_err_scaled = "
            "error x sqrt(chi2/ndf) when > 1; Borel expectation k(T) = k(0.5) + dlnP_Borel(N>=T|p)/dlnp * k_pct with p = 0.429\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["quantity", "variant", "T_pe", "n_bins", "k", "k_err_scaled", "k_err_stat", "chi2", "ndf"])
    w.writerows(out_rows)
json.dump(res, open(f"{D}/diurnal_slope.json", "w"), indent=1)
for r in out_rows:
    print(" ".join(str(v) for v in r))
