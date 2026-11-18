#!/usr/bin/python3 -I
# epoxy_mixture.py — seguimiento 2026-11-18, punto 2: límite superior al número de canales de hod2018 con una pérdida relativa ≥ 26 % (respuesta de
# Geant4 FU18 de una junta de epoxi despegada, separación ≥ 0.05 mm), con un modelo de mezcla bien definido. Reemplaza el límite de k de
# iso_xtalk.py (k canales elegidos para maximizar L): esa elección tiene libertad combinatoria y −2Δln L no sigue la distribución de χ² (dio k_mejor = 12).
# Datos: razones oct/mayo por canal con la selección 29-54 mm (outputs/hod2018_ratio_reconcile/2026-11-16/r1_C.csv, variante thr3_far), 31 canales.
# Modelo: cada canal, independientemente, despegado con probabilidad f: r_i ~ (1 − f)·N(μ, σ_i² + τ²) + f·N(0.74·μ, σ_i² + 0.74²·τ²).
#   M0: τ = 0; M1: τ libre. Perfil de −2 ln L en f (μ, τ libres en rejilla); f_95 = el f con −2Δln L = 2.71 por encima del mínimo (una cola);
#   k_95 = f_95 · 31 [calculado]. También el número esperado de canales despegados a f_95 y P(≥ 1 despegado) = 1 − (1 − f_95)^31.
# f ≤ f_max (argumento; 1.0 por defecto). Con f libre hasta 1 el modelo es degenerado: «todos despegados» con μ/0.74 es un factor global (M0 dio
# f_mejor = 0.92); la pregunta es por una minoría de canales con pérdida relativa a los demás, de ahí la variante f_max = 0.5 (sufijo _v2).
# Salida: <out>/epoxy_mixture<sufijo>.json, epoxy_mixture<sufijo>.csv (perfil). Uso: epoxy_mixture.py <out dir> [f_max [sufijo]]
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
FMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
SUF = sys.argv[3] if len(sys.argv) > 3 else ""
rc = [r for r in csv.DictReader(open(REPO + "/outputs/hod2018_ratio_reconcile/2026-11-16/r1_C.csv")) if r["variant"] == "thr3_far"]
R = np.array([float(r["ratio_oct_over_may"]) for r in rc])
E = np.array([float(r["err"]) for r in rc])
LOSS = 0.74


def m2ll(f, mu, tau):
    v0 = E ** 2 + tau ** 2
    v1 = E ** 2 + (LOSS * tau) ** 2
    l0 = np.exp(-0.5 * (R - mu) ** 2 / v0) / np.sqrt(2 * np.pi * v0)
    l1 = np.exp(-0.5 * (R - LOSS * mu) ** 2 / v1) / np.sqrt(2 * np.pi * v1)
    return -2 * float(np.log((1 - f) * l0 + f * l1 + 1e-300).sum())


fgrid = np.round(np.arange(0, FMAX + 1e-9, 0.005), 4)
mus = np.linspace(0.55, 0.95, 401)
taus = {"M0_tau0": np.array([0.0]), "M1_tau_free": np.linspace(0, 0.15, 151)}
res, prow = {}, []
for name, tg in taus.items():
    prof = []
    for f in fgrid:
        best = (1e300, None, None)
        for tau in tg:
            v = [m2ll(f, m, tau) for m in mus]
            i = int(np.argmin(v))
            if v[i] < best[0]:
                best = (v[i], float(mus[i]), float(tau))
        prof.append((float(f),) + best)
    v = np.array([p[1] for p in prof])
    ib = int(np.argmin(v))
    d = v - v[ib]
    f95 = float(fgrid[ib:][np.argmax(d[ib:] >= 2.71)]) if (d[ib:] >= 2.71).any() else float("nan")  # nan: sin límite en [0, f_max]
    res[name] = {"f_best": float(fgrid[ib]), "m2ll_best": float(v[ib]), "mu_best": prof[ib][2], "tau_best": prof[ib][3], "m2ll_f0_minus_best": float(d[0]),
                 "f_95CL_upper": f95, "k_95CL_upper_eq": f95 * len(R), "P_at_least_one_at_f95": 1 - (1 - f95) ** len(R), "n_channels": len(R)}
    for p, dd in zip(prof, d):
        prow.append([name, p[0], f"{p[1]:.4f}", f"{dd:.4f}", f"{p[2]:.4f}", f"{p[3]:.4f}"])
    print(name, json.dumps(res[name]), flush=True)
json.dump(res, open(OUT + f"/epoxy_mixture{SUF}.json", "w"), indent=1)
with open(OUT + f"/epoxy_mixture{SUF}.csv", "w") as fo:
    w = csv.writer(fo, lineterminator="\n")
    w.writerow(["model", "f", "m2lnL", "delta", "mu", "tau"])
    w.writerows(prow)
