#!/usr/bin/python3 -I
# geom_sensitivity.py — DAQ2-15: sensibilidad a la geometría de la fracción de triggers externos que cruza el área
# activa de hod2018 (f_2018), con rayos rectos (sin Geant4).
# Geometría (analysis/python/geom_acceptance.py, constantes del MC de mayo; z hacia arriba, muones bajan):
#   hod2018 z = 0, área activa ±16.5 mm (33 × 33 mm); T1 z = 124.1 mm, 150 × 150 mm, centro (−42.5, −42.5);
#   T2 z = 134.1 mm, 150 × 150 mm, centro (+42.5, +42.5) -> superposición de 65 × 65 mm centrada en (0, 0).
# Distribución angular: intensidad ∝ cos²θ por ángulo sólido; con posiciones uniformes en un plano horizontal el
# flujo lleva además el factor de proyección cosθ (dN/dΩ dA ∝ cos³θ). Se dan las dos: "cos2" reproduce
# geom_acceptance.py (sin proyección, 2.80 %); "cos3" es la que corrigió el handoff (f_both 1.27 % -> 1.67 %).
# Variaciones con los mismos rayos (diferencias correlacionadas): estación o placas desplazadas 2, 5, 10 mm en x o
# y, y rotadas 2° alrededor del eje vertical (estación: alrededor de su centro; placas: alrededor de (0, 0)).
# Error: f binomial; diferencia con la nominal = sqrt(ganados + perdidos) / triggers.
# Uso: geom_sensitivity.py <salida.json> [N]
import json, math, sys

import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_mc  # noqa: E402

Z_START, Z_T2, Z_T1, Z_H18 = 160.0, 134.1, 124.1, 0.0
HALF = 16.5
T1C, T2C, TH = (-42.5, -42.5), (42.5, 42.5), 75.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20_000_000
CHUNK = 2_000_000


def rot(x, y, cx, cy, deg):
    a = math.radians(deg)
    dx, dy = x - cx, y - cy
    return cx + dx * math.cos(a) + dy * math.sin(a), cy - dx * math.sin(a) + dy * math.cos(a)


VARS = [("nominal", {})]
for d in (2, 5, 10):
    for ax in ("x", "y"):
        VARS.append((f"estación {ax} +{d} mm", {"st_" + ax: d}))
        VARS.append((f"placas {ax} +{d} mm", {"pl_" + ax: d}))
VARS.append(("estación rotada 2°", {"st_rot": 2.0}))
VARS.append(("placas rotadas 2°", {"pl_rot": 2.0}))


def run(model, seed=42):
    rng = np.random.default_rng(seed)
    cnt = {k: 0 for k, _ in VARS}
    trig_n = {k: 0 for k, _ in VARS}
    gain = {k: 0 for k, _ in VARS}
    loss = {k: 0 for k, _ in VARS}
    tot = 0
    for start in range(0, N, CHUNK):
        m = min(CHUNK, N - start)
        rays = bench_mc.generate(rng, m, model)  # generador compartido (PROC-1), mismo orden de llamadas que antes
        p2, p1, p18 = rays.prop(Z_T2), rays.prop(Z_T1), rays.prop(Z_H18)
        nom = None
        for k, v in VARS:
            sx, sy = v.get("pl_x", 0.0), v.get("pl_y", 0.0)
            a2x, a2y = p2[0] - sx, p2[1] - sy
            a1x, a1y = p1[0] - sx, p1[1] - sy
            if "pl_rot" in v:
                a2x, a2y = rot(a2x, a2y, 0, 0, v["pl_rot"])
                a1x, a1y = rot(a1x, a1y, 0, 0, v["pl_rot"])
            trig = (np.abs(a2x - T2C[0]) <= TH) & (np.abs(a2y - T2C[1]) <= TH) & (np.abs(a1x - T1C[0]) <= TH) & (np.abs(a1y - T1C[1]) <= TH)
            hx, hy = p18[0] - v.get("st_x", 0.0), p18[1] - v.get("st_y", 0.0)
            if "st_rot" in v:
                hx, hy = rot(hx, hy, 0, 0, v["st_rot"])
            ok = trig & (np.abs(hx) <= HALF) & (np.abs(hy) <= HALF)
            trig_n[k] += int(trig.sum())
            cnt[k] += int(ok.sum())
            if k == "nominal":
                nom = (trig, ok)
            else:
                # cambio de f = ok/trig: rayos que entran o salen (en trig o en ok)
                gain[k] += int((ok & ~nom[1]).sum())
                loss[k] += int((~ok & nom[1]).sum())
        tot += m
    out = {}
    fn = cnt["nominal"] / trig_n["nominal"]
    for k, _ in VARS:
        f = cnt[k] / trig_n[k]
        e = math.sqrt(f * (1 - f) / trig_n[k])
        d = {"f": f, "f_err": e, "n_trig": trig_n[k], "n_in": cnt[k]}
        if k != "nominal":
            d["rel_change"] = f / fn - 1
            d["rel_change_err"] = math.sqrt(gain[k] + loss[k]) / cnt["nominal"]
        out[k] = d
    return out, tot


res = {"geometry": {"source": "analysis/python/geom_acceptance.py (constantes); handoff 2026-10-18 línea 125 (3.60 %) y 163 (factor de proyección)",
                    "hod2018_z_mm": Z_H18, "hod2018_active_half_mm": HALF, "T1": {"z_mm": Z_T1, "center_mm": T1C, "side_mm": 2 * TH},
                    "T2": {"z_mm": Z_T2, "center_mm": T2C, "side_mm": 2 * TH}, "overlap_mm": "65 × 65 centrada en (0, 0)", "rays": N, "seed": 42}}
for model in ("cos2", "cos3"):
    res[model], _ = run(model)
json.dump(res, open(sys.argv[1], "w"), ensure_ascii=False, indent=1)
for model in ("cos2", "cos3"):
    print(model)
    for k, d in res[model].items():
        print(f"  {k:24s} f = {d['f']:.5f} ± {d['f_err']:.5f}" + (f"  cambio {100 * d['rel_change']:+.2f} ± {100 * d['rel_change_err']:.2f} %" if "rel_change" in d else ""))
