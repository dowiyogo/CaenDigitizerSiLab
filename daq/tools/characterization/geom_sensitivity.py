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
# PROC-1 etapa 6: también f_2019 (cruza el área activa de hod2019) y f_both (las dos), las variaciones D ± 1 mm (z de hod2019)
# y placas z ± 2 mm (T1 y T2 juntas), y la geometría de las placas a elegir (bench_mc.GEOMS: OLD = la de arriba, CAD = A5).
# Con la geometría OLD y las variaciones de antes, f y rel_change de hod2018 no cambian (mismos rayos).
# En las variaciones de la estación (desplazada o rotada) hod2018 y hod2019 se mueven juntas.
# Uso: geom_sensitivity.py <salida.json> [N] [OLD|CAD]
import json, math, sys

import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_mc  # noqa: E402

Z_START, Z_H18 = 160.0, 0.0
HALF = 16.5
T1C, T2C, TH = (-42.5, -42.5), (42.5, 42.5), 75.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20_000_000
GEOM = sys.argv[3] if len(sys.argv) > 3 else "OLD"
Z_T2, Z_T1, Z_H19 = bench_mc.GEOMS[GEOM]["T2"], bench_mc.GEOMS[GEOM]["T1"], bench_mc.GEOMS[GEOM]["H19"]
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
for d in (1.0, -1.0):
    VARS.append((f"D {d:+.0f} mm", {"dD": d}))
for d in (2.0, -2.0):
    VARS.append((f"placas z {d:+.0f} mm", {"pl_z": d}))


def run(model, seed=42):
    rng = np.random.default_rng(seed)
    cnt = {k: 0 for k, _ in VARS}
    trig_n = {k: 0 for k, _ in VARS}
    gain = {k: 0 for k, _ in VARS}
    loss = {k: 0 for k, _ in VARS}
    c19 = {k: 0 for k, _ in VARS}
    cboth = {k: 0 for k, _ in VARS}
    g19 = {k: [0, 0, 0, 0] for k, _ in VARS}  # ganados/perdidos de hod2019 y de las dos
    tot = 0
    for start in range(0, N, CHUNK):
        m = min(CHUNK, N - start)
        rays = bench_mc.generate(rng, m, model)  # generador compartido (PROC-1), mismo orden de llamadas que antes
        p2, p1, p18, p19 = rays.prop(Z_T2), rays.prop(Z_T1), rays.prop(Z_H18), rays.prop(Z_H19)
        nom = None
        for k, v in VARS:
            if "pl_z" in v:
                p2, p1 = rays.prop(Z_T2 + v["pl_z"]), rays.prop(Z_T1 + v["pl_z"])
            else:
                p2, p1 = rays.prop(Z_T2), rays.prop(Z_T1)
            q19 = rays.prop(Z_H19 + v.get("dD", 0.0))
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
            ok19 = trig & (np.abs(q19[0] - v.get("st_x", 0.0)) <= HALF) & (np.abs(q19[1] - v.get("st_y", 0.0)) <= HALF) if "st_rot" not in v else \
                trig & (np.abs(rot(q19[0], q19[1], 0, 0, v["st_rot"])[0]) <= HALF) & (np.abs(rot(q19[0], q19[1], 0, 0, v["st_rot"])[1]) <= HALF)
            both = ok & ok19
            trig_n[k] += int(trig.sum())
            cnt[k] += int(ok.sum())
            c19[k] += int(ok19.sum())
            cboth[k] += int(both.sum())
            if k == "nominal":
                nom = (trig, ok, ok19, both)
            else:
                g19[k][0] += int((ok19 & ~nom[2]).sum())
                g19[k][1] += int((~ok19 & nom[2]).sum())
                g19[k][2] += int((both & ~nom[3]).sum())
                g19[k][3] += int((~both & nom[3]).sum())
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
        for lab, cc in (("f_2019", c19), ("f_both", cboth)):
            v_ = cc[k] / trig_n[k]
            d[lab], d[lab + "_err"] = v_, math.sqrt(v_ * (1 - v_) / trig_n[k])
        if k != "nominal":
            d["rel_change"] = f / fn - 1
            d["rel_change_err"] = math.sqrt(gain[k] + loss[k]) / cnt["nominal"]
            d["rel_change_2019"] = d["f_2019"] / (c19["nominal"] / trig_n["nominal"]) - 1
            d["rel_change_2019_err"] = math.sqrt(g19[k][0] + g19[k][1]) / c19["nominal"]
            d["rel_change_both"] = d["f_both"] / (cboth["nominal"] / trig_n["nominal"]) - 1
            d["rel_change_both_err"] = math.sqrt(g19[k][2] + g19[k][3]) / cboth["nominal"]
        out[k] = d
    return out, tot


res = {"geometry": {"source": "analysis/python/geom_acceptance.py (constantes); handoff 2026-10-18 línea 125 (3.60 %) y 163 (factor de proyección)",
                    "plates": GEOM, "hod2019_z_mm": Z_H19, "hod2018_z_mm": Z_H18, "hod2018_active_half_mm": HALF, "T1": {"z_mm": Z_T1, "center_mm": T1C, "side_mm": 2 * TH},
                    "T2": {"z_mm": Z_T2, "center_mm": T2C, "side_mm": 2 * TH}, "overlap_mm": "65 × 65 centrada en (0, 0)", "rays": N, "seed": 42}}
for model in ("cos2", "cos3"):
    res[model], _ = run(model)
json.dump(res, open(sys.argv[1], "w"), ensure_ascii=False, indent=1)
for model in ("cos2", "cos3"):
    print(model)
    for k, d in res[model].items():
        print(f"  {k:24s} f = {d['f']:.5f} ± {d['f_err']:.5f}" + (f"  cambio {100 * d['rel_change']:+.2f} ± {100 * d['rel_change_err']:.2f} %" if "rel_change" in d else ""))
