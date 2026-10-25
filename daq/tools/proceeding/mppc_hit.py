#!/usr/bin/python3 -I
# mppc_hit.py — XPLUS-1 parte C.1 (MC, ESTIMADO): hipótesis de impacto directo en los MPPC de hod2019.
# MC del banco (bench_mc.py, z del CAD, cos³θ, 1e8 rayos, semilla 42), trayectorias del trigger T1 ∧ T2. Supuestos (sin registro):
#  - un MPPC por barra, en el extremo: X capa A en y = +37.5 mm, X capa B en y = −37.5 mm (barras X a lo largo de y); Y capa A en
#    x = +37.5 mm, Y capa B en x = −37.5 mm; el MPPC está centrado en la barra (x de la barra para las X, y para las Y) y en el plano
#    medio de su capa (z de bench_mc.layer_z, módulo de hod2019 a 98.4 mm);
#  - área activa 3 × 3 mm (ancho de la barra × 3 mm en z), con 0.3 mm de silicio de profundidad fuera del extremo de la barra
#    (supuesto); encapsulado cerámico de 4.3 × 4.3 mm con 1.4 mm de profundidad (supuesto, como pide la tarea);
#  - "cruza" = la recta atraviesa la caja (método de las placas). Los dos órdenes de bloque X/Y; se da el promedio y la mitad de la
#    diferencia.
# Comparación con el run B (h19_xy.py): rX, rY a 35 y 3 mV y las 3 coincidencias X·Y en el mismo evento (N = 143 814).
# Salida: work/mppc_hit.json. Uso: mppc_hit.py <dir del paquete> [N]
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import bench_mc as bm  # noqa: E402

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000_000
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
G = bm.GEOMS["CAD"]
HL = bm.BAR_L / 2
SHAPES = {"sensor_3x3": {"half_w": 1.5, "half_z": 1.5, "depth": 0.3}, "package_4p3x4p3": {"half_w": 2.15, "half_z": 2.15, "depth": 1.4}}


def box_hit(r, lo, hi):
    """recta x = x0 + dx·t, …: ¿atraviesa la caja [lo, hi] (por eje)? lo/hi: listas de 3 arreglos o escalares."""
    o = (r.x0, r.y0, np.full_like(r.x0, bm.Z_START))
    d = (r.dx, r.dy, r.dz)
    tmin = np.full_like(r.x0, -np.inf)
    tmax = np.full_like(r.x0, np.inf)
    for k in range(3):
        with np.errstate(divide="ignore", invalid="ignore"):
            t1 = (lo[k] - o[k]) / d[k]
            t2 = (hi[k] - o[k]) / d[k]
        par = np.abs(d[k]) < 1e-12
        inside = (o[k] >= lo[k]) & (o[k] <= hi[k])
        t1 = np.where(par, np.where(inside, -np.inf, np.inf), t1)
        t2 = np.where(par, np.where(inside, np.inf, -np.inf), t2)
        tmin = np.maximum(tmin, np.minimum(t1, t2))
        tmax = np.minimum(tmax, np.maximum(t1, t2))
    return tmax >= tmin


def mppc_any(r, top, coord, shape):
    """¿cruza algún MPPC de la coordenada `coord` de hod2019?"""
    s = SHAPES[shape]
    zl = bm.layer_z(G["H19"], top)
    hit = np.zeros(len(r.x0), bool)
    for layer, centres, sign in (("A", bm.LAYER_A, +1), ("B", bm.LAYER_B, -1)):
        z0 = zl[(coord, layer)]
        end_lo, end_hi = (HL, HL + s["depth"]) if sign > 0 else (-HL - s["depth"], -HL)
        for c in centres:
            if coord == "X":  # barra X a lo largo de y: el MPPC en y = ±37.5, centrado en x = c
                lo, hi = (c - s["half_w"], end_lo, z0 - s["half_z"]), (c + s["half_w"], end_hi, z0 + s["half_z"])
            else:             # barra Y a lo largo de x
                lo, hi = (end_lo, c - s["half_w"], z0 - s["half_z"]), (end_hi, c + s["half_w"], z0 + s["half_z"])
            hit |= box_hit(r, lo, hi)
    return hit


cnt = {(shape, top, k): 0 for shape in SHAPES for top in ("X", "Y") for k in ("X", "Y", "XY")}
trig = 0
rng = np.random.default_rng(42)
for start in range(0, N, 2_000_000):
    m = min(2_000_000, N - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, G)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    trig += int(t.sum())
    for shape in SHAPES:
        for top in ("X", "Y"):
            hx, hy = mppc_any(rt, top, "X", shape), mppc_any(rt, top, "Y", shape)
            cnt[(shape, top, "X")] += int(hx.sum())
            cnt[(shape, top, "Y")] += int(hy.sum())
            cnt[(shape, top, "XY")] += int((hx & hy).sum())
NB = 143814
MEAS = {"35mV": {"rX": 0.01441, "rY": 0.00934, "N_XY": 3}, "3mV": {"rX": 0.08768, "rY": 0.07554, "N_XY": 4087, "qN_mix": 128.7}}
S = {"stage": "XPLUS-1 C.1", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "status": "estimated (MC con supuestos)",
     "rays": N, "triggers": trig, "geometry": "CAD", "model": "cos3", "seed": 42, "shapes": SHAPES,
     "z_used_mm": {top: {f"{c}{l}": v for (c, l), v in bm.layer_z(G["H19"], top).items()} for top in ("X", "Y")},
     "assumptions": ("un MPPC por barra en su extremo: X capa A en y = +37.5 mm y capa B en y = −37.5 mm, Y capa A en x = +37.5 mm y capa B en "
                     "x = −37.5 mm (sin registro de qué extremo); centrado en la barra y en el plano medio de su capa; sensor 3 × 3 mm con 0.3 mm de "
                     "silicio (supuesto), encapsulado cerámico 4.3 × 4.3 mm con 1.4 mm (supuesto); los dos órdenes de bloque"),
     "per_shape": {}, "measured": MEAS}
for shape in SHAPES:
    d = {}
    for k in ("X", "Y", "XY"):
        v = {top: cnt[(shape, top, k)] / trig for top in ("X", "Y")}
        avg = 0.5 * (v["X"] + v["Y"])
        n = 0.5 * (cnt[(shape, "X", k)] + cnt[(shape, "Y", k)])
        d[f"P_{k}"] = {"value": avg, "stat": math.sqrt(n) / trig if n > 0 else 1.0 / trig, "syst_order": 0.5 * abs(v["X"] - v["Y"]), "counts_avg": n}
    d["expected_N_XY_run_B"] = d["P_XY"]["value"] * NB
    d["ratio_to_measured"] = {thr: {"rX_over_PX": m["rX"] / d["P_X"]["value"] if d["P_X"]["value"] > 0 else None,
                                    "rY_over_PY": m["rY"] / d["P_Y"]["value"] if d["P_Y"]["value"] > 0 else None} for thr, m in MEAS.items()}
    S["per_shape"][shape] = d
    print(shape, "P(X) %.3e ± %.1e  P(Y) %.3e ± %.1e  P(X∧Y) %.3e  N_XY esperadas %.3f" % (d["P_X"]["value"], d["P_X"]["stat"], d["P_Y"]["value"],
          d["P_Y"]["stat"], d["P_XY"]["value"], d["expected_N_XY_run_B"]), {t: {k: (round(v, 1) if v else v) for k, v in x.items()} for t, x in d["ratio_to_measured"].items()})
json.dump(S, open(W + "/mppc_hit.json", "w"), indent=1, default=float)
