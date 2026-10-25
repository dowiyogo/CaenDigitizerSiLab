#!/usr/bin/python3 -I
# block_fractions.py — PROC-3 etapa 1 (MC, estimado): fracción de los triggers (T1 ∧ T2) cuyas trayectorias cruzan un bloque de barras
# completo de cada estación, para la línea del manuscrito "about 0.4 of the triggers cross the full 33 mm × 75 mm of the X bars"
# (sección 6.1, \verif). MC del banco (bench_mc.py, z del CAD, cos³θ, 1e8 rayos, semilla 42).
#  (a) bloque X: la trayectoria está dentro de |x| ≤ 16.5 mm (33 mm de barras) y |y| ≤ 37.5 mm (75 mm de largo) en los planos medios
#      de las dos capas X;
#  (b) bloque Y: |y| ≤ 16.5 mm y |x| ≤ 37.5 mm en las dos capas Y;
#  (c) bloque X dentro de la superposición de las placas (65 × 65 mm centrada): además |x|, |y| ≤ 32.5 mm en esas dos capas.
# Para hod2019 (z = 98.4 mm) y hod2018 (z = 0). El orden de los bloques X/Y dentro de la estación no tiene registro: se dan los dos y el
# promedio (bench_mc.layer_z, capa A arriba dentro de cada bloque; el orden de las capas no cambia el criterio de las dos capas).
# Salida: work/block_fractions.json. Uso: block_fractions.py <dir del paquete> [N]
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
HB, HL, HO = 16.5, 37.5, 32.5  # semiancho del bloque, semilargo de las barras, semilado de la superposición de las placas
cnt = {}
trig = 0
rng = np.random.default_rng(42)
for start in range(0, N, 2_000_000):
    m = min(2_000_000, N - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, G)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    trig += int(t.sum())
    for st, zm in (("hod2019", G["H19"]), ("hod2018", G["H18"])):
        for top in ("X", "Y"):
            zl = bm.layer_z(zm, top)
            px = [rt.prop(zl[("X", l)]) for l in "AB"]
            py = [rt.prop(zl[("Y", l)]) for l in "AB"]
            a = np.all([(np.abs(x) <= HB) & (np.abs(y) <= HL) for x, y in px], axis=0)
            b = np.all([(np.abs(y) <= HB) & (np.abs(x) <= HL) for x, y in py], axis=0)
            c = a & np.all([(np.abs(x) <= HO) & (np.abs(y) <= HO) for x, y in px], axis=0)
            for k, v in (("a_X_block", a), ("b_Y_block", b), ("c_X_block_in_overlap", c)):
                cnt[(st, top, k)] = cnt.get((st, top, k), 0) + int(v.sum())
res = {"stage": "PROC-3 etapa 1", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "rays": N, "triggers": trig,
       "geometry": "CAD", "model": "cos3", "seed": 42, "stations": {}}
for st in ("hod2019", "hod2018"):
    d = {}
    for k in ("a_X_block", "b_Y_block", "c_X_block_in_overlap"):
        v = {top: cnt[(st, top, k)] / trig for top in ("X", "Y")}
        avg = 0.5 * (v["X"] + v["Y"])
        d[k] = {"X_block_on_top": v["X"], "Y_block_on_top": v["Y"], "value": avg, "stat": math.sqrt(avg * (1 - avg) / trig),
                "syst_order": 0.5 * abs(v["X"] - v["Y"])}
        print(st, k, "%.5f ± %.5f (orden de bloques ± %.5f)" % (avg, d[k]["stat"], d[k]["syst_order"]))
    res["stations"][st] = d
json.dump(res, open(W + "/block_fractions.json", "w"), indent=1)
