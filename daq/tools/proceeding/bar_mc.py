#!/usr/bin/python3 -I
# bar_mc.py — PROC-1 etapa 2.3 (y base de la etapa 5): factores geométricos de hod2018 a nivel de barra con el MC del banco
# (daq/tools/characterization/bench_mc.py; geometría de barras del CAD de René, 2026-10-24).
# Entre los triggers (T1 ∧ T2):
#   S = cruza el cuadrado activo de 33 × 33 mm en z = 0 (la definición de f_2018 de geom_acceptance.py);
#   X, Y = cruza al menos una barra X (Y) con longitud de camino > 0;
#   g     = P(X ∧ Y | S)       (A6: fracción de las trayectorias que cruzan hod2018 que dan impacto X e Y con detección perfecta);
#   g_XY  = P(X | Y), g_YX = P(Y | X)  (para ε_X = ε(X|Y)/g_XY, ε_Y = ε(Y|X)/g_YX);
#   f_XY  = P(X ∧ Y)            (fracción de triggers con impacto X e Y geométrico, sin pasar por S).
# Se dan para el bloque X arriba y el bloque Y arriba (no hay registro de cuál está arriba), con la geometría OLD de las placas
# (la de f_2018 = 0.03646) y la del CAD.
# PROC-1 etapa 3: además, cruces por tira (barra ordenada por centro, −15 … 15 mm = coordenada 0 … 15) de hod2018 (bloque X
# arriba) y hod2019 (z = 98.4 mm), por trigger: la expectativa de impactos por tira con detección perfecta.
# Uso: bar_mc.py <salida.json> [N] [semilla]
import json, math, os, sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "characterization"))
import bench_mc as bm  # noqa: E402

OUTF = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 50_000_000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 42
CHUNK = 2_000_000


def run(geom, model="cos3"):
    g = bm.GEOMS[geom]
    rng = np.random.default_rng(SEED)
    c = {k: 0 for k in ("trig", "S", "XY_S", "X", "Y", "XY", "XY_S_topY", "X_topY", "Y_topY", "XY_topY")}
    strips = {f"{st}_{co}": np.zeros(16) for st in ("h18", "h19") for co in ("X", "Y")}  # PROC-1 etapa 3: cruces por tira
    for start in range(0, N, CHUNK):
        m = min(CHUNK, N - start)
        r = bm.generate(rng, m, model)
        t = bm.trigger(r, g)
        rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
        S = bm.in_square(rt, g["H18"])
        c["trig"] += int(t.sum())
        c["S"] += int(S.sum())
        for top, suf in (("X", ""), ("Y", "_topY")):
            p = bm.bar_paths(rt, g["H18"], top)
            X, Y = (p["X"] > 0).any(1), (p["Y"] > 0).any(1)
            c["XY_S" + suf] += int((X & Y & S).sum())
            c["X" + suf] += int(X.sum())
            c["Y" + suf] += int(Y.sum())
            c["XY" + suf] += int((X & Y).sum())
            if top == "X":
                strips["h18_X"] += (p["X"] > 0).sum(0)
                strips["h18_Y"] += (p["Y"] > 0).sum(0)
        p19 = bm.bar_paths(rt, g["H19"], "X")
        strips["h19_X"] += (p19["X"] > 0).sum(0)
        strips["h19_Y"] += (p19["Y"] > 0).sum(0)
    out = {"geometry": geom, "model": model, "rays": N, "seed": SEED, "counts": c,
           "strip_crossings": {k: v.tolist() for k, v in strips.items()},
           "strip_fraction_per_trigger": {k: (v / c["trig"]).tolist() for k, v in strips.items()},
           "strip_spread": {k: {"max_over_min": float(v.max() / v.min()), "rel_sd": float(v.std() / v.mean())} for k, v in strips.items()}}
    for suf, lab in (("", "X arriba"), ("_topY", "Y arriba")):
        def frac(k, n):
            p = k / n
            return {"value": p, "err": math.sqrt(max(p * (1 - p), 0) / n)}
        out[lab] = {"g": frac(c["XY_S" + suf], c["S"]), "g_XY": frac(c["XY" + suf], c["Y" + suf]), "g_YX": frac(c["XY" + suf], c["X" + suf]),
                    "f_XY": frac(c["XY" + suf], c["trig"]), "f_X": frac(c["X" + suf], c["trig"]), "f_Y": frac(c["Y" + suf], c["trig"])}
    out["f_2018"] = {"value": c["S"] / c["trig"], "err": math.sqrt(c["S"] / c["trig"] * (1 - c["S"] / c["trig"]) / c["trig"])}
    return out


res = {"bars": {"layer_A_centres_mm": bm.LAYER_A.tolist(), "layer_B_centres_mm": bm.LAYER_B.tolist(), "width_mm": bm.BAR_W, "thickness_mm": bm.BAR_T,
                "length_mm": bm.BAR_L, "layer_dz_mm": bm.LAYER_DZ, "block_dz_mm": bm.BLOCK_DZ, "source": "CAD (René, 2026-10-24)"},
       "runs": [run("OLD"), run("CAD")]}
json.dump(res, open(OUTF, "w"), indent=1)
for r in res["runs"]:
    print(r["geometry"], "f_2018 %.5f ± %.5f" % (r["f_2018"]["value"], r["f_2018"]["err"]), "triggers", r["counts"]["trig"])
    for lab in ("X arriba", "Y arriba"):
        d = r[lab]
        print("  ", lab, " ".join(f"{k} {v['value']:.4f} ± {v['err']:.4f}" for k, v in d.items()))
