#!/usr/bin/python3 -I
# g_dead_ch3.py — 2026-11-02: g de hod2018 con la tira X12 (ch3, centro 9 mm, capa A) muerta, con el MC del banco.
#  - Muestra: la de stage5.py (f_2): bench_mc.py, cos³θ, placas del CAD, semilla 42, trozos de 2e6, triggers T1 ∧ T2; los primeros
#    30e6 rayos son exactamente los de stage5 (791485 triggers). Se sigue hasta 1e8 (la muestra de bar_mc.py / bar_mc_g.json, de la que
#    sale g = 0.980) para comprobar que el valor de 16 barras se reproduce y para reducir el error estadístico.
#  - g = P(X ∧ Y | S), S = cruza el cuadrado de 33 × 33 mm en z = 0; X, Y = camino > 0 en alguna barra X (Y). g_dead: igual con la
#    columna X 12 apagada (como H["X"][:, DEAD_X] = False en stage5). Bloque X arriba y bloque Y arriba.
#  - ε_XY = q / (f_2018 · g) con q de stage2.json (mayo v8 tiempo bueno, ch3 fuera), f_2018 del CAD; incertidumbre de aceptancia
#    = f_2018 (3.50 %) ⊕ g (3.06 %), sin el término de cobertura de ch3.
# Salida: <out>/g_dead_ch3.json. Uso: g_dead_ch3.py <dir del paquete> <out dir> [N]
import json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import bench_mc as bm  # noqa: E402

PKG, OUT = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000_000
N5 = 30_000_000
DEAD_X = 12
geo = bm.GEOMS["CAD"]
rng = np.random.default_rng(42)
keys = ("trig", "S") + tuple(f"{k}_{top}" for top in "XY" for k in ("XY_S", "XY_S_dead", "X", "X_dead", "Y", "XY", "XY_dead"))
c = dict.fromkeys(keys, 0)
snap = None
for start in range(0, N, 2_000_000):
    m = min(2_000_000, N - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, geo)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    S = bm.in_square(rt, geo["H18"])
    c["trig"] += int(t.sum())
    c["S"] += int(S.sum())
    for top in "XY":
        p = bm.bar_paths(rt, geo["H18"], top)
        hx = p["X"] > 0
        X, Y = hx.any(1), (p["Y"] > 0).any(1)
        hx[:, DEAD_X] = False
        Xd = hx.any(1)
        for k, v in (("XY_S", X & Y & S), ("XY_S_dead", Xd & Y & S), ("X", X), ("X_dead", Xd), ("Y", Y), ("XY", X & Y), ("XY_dead", Xd & Y)):
            c[f"{k}_{top}"] += int(v.sum())
    if start + m == N5:
        snap = dict(c)


def frac(k, n):
    p = k / n
    return {"value": p, "err": math.sqrt(max(p * (1 - p), 0) / n)}


def summary(c):
    d = {"counts": c}
    for top in "XY":
        n16, nd = c[f"XY_S_{top}"], c[f"XY_S_dead_{top}"]
        d[f"{top}_top"] = {"g16": frac(n16, c["S"]), "g_dead": frac(nd, c["S"]),
                           # razón g_dead/g16 con muestras anidadas: error binomial de (n16 − nd) / n16
                           "ratio": frac(nd, n16), "f_2018": frac(c["S"], c["trig"])}
    d["average"] = {k: 0.5 * (d["X_top"][k]["value"] + d["Y_top"][k]["value"]) for k in ("g16", "g_dead", "ratio")}
    return d


res = {"rays": N, "seed": 42, "model": "cos3", "geometry": "CAD", "dead_column": DEAD_X, "dead_centre_mm": -15.0 + 2.0 * DEAD_X,
       "stage5_sample": summary(snap) if snap else None, "full": summary(c)}

# ε_XY con la g nueva
s2 = json.load(open(PKG + "/work/stage2.json"))
F18 = s2["f_2018"]
F18_REL = s2["acceptance_rel_syst"]
G_SYST_REL = 0.03 / 0.98
use = res["stage5_sample"] if snap else res["full"]
g_new = use["average"]["g_dead"]
res["g_used_new"] = {"value": g_new, "ratio_to_0.98": g_new / 0.98, "X_top_ratio_to_0.98": use["X_top"]["g_dead"]["value"] / 0.98,
                     "Y_top_ratio_to_0.98": use["Y_top"]["g_dead"]["value"] / 0.98,
                     "rule": "g_dead de la muestra de stage5, promedio de los dos órdenes de bloque (como 0.980 para 16 barras)"}
res["acceptance_rel_syst_no_ch3"] = math.sqrt(F18_REL ** 2 + G_SYST_REL ** 2)
res["eps_XY"] = {}
for T in ("5", "10"):
    rr = s2["datasets"]["may8_good"]["per_T"][T]
    q, qe, qc = rr["f_mix"], rr["f_mix_err"], rr["f_ctl"]
    e = q / (F18 * g_new)
    res["eps_XY"][T] = {"q": q, "value": e, "stat": qe / (F18 * g_new), "syst_estimator": abs(q - qc) / (F18 * g_new),
                        "syst_acceptance_no_ch3": e * res["acceptance_rel_syst_no_ch3"], "old_g0.98": rr["eps_XY"]["value"]}
os.makedirs(OUT, exist_ok=True)
json.dump(res, open(OUT + "/g_dead_ch3.json", "w"), indent=1)
for lab in ("stage5_sample", "full"):
    d = res[lab]
    if d is None:
        continue
    print(lab, "triggers", d["counts"]["trig"], "S", d["counts"]["S"])
    for top in "XY":
        x = d[f"{top}_top"]
        print(f"  {top} arriba: g16 {x['g16']['value']:.4f} ± {x['g16']['err']:.4f}  g_dead {x['g_dead']['value']:.4f} ± {x['g_dead']['err']:.4f}"
              f"  razón {x['ratio']['value']:.4f} ± {x['ratio']['err']:.4f}")
    print("  promedio: g16 %.4f g_dead %.4f razón %.4f" % (d["average"]["g16"], d["average"]["g_dead"], d["average"]["ratio"]))
u = res["g_used_new"]
print("g nueva = %.4f (X arriba / 0.98 = %.4f, Y arriba / 0.98 = %.4f, promedio / 0.98 = %.4f); aceptancia sin ch3 = %.4f"
      % (g_new, u["X_top_ratio_to_0.98"], u["Y_top_ratio_to_0.98"], u["ratio_to_0.98"], res["acceptance_rel_syst_no_ch3"]))
for T, x in res["eps_XY"].items():
    print(f"T {T}: q {x['q']:.5f} ε_XY {x['value']:.4f} ± {x['stat']:.4f} (stat) ± {x['syst_estimator']:.4f} (estimador) ± "
          f"{x['syst_acceptance_no_ch3']:.4f} (aceptancia); con g = 0.98: {x['old_g0.98']:.4f}")
