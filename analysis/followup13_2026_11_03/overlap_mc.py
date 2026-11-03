#!/usr/bin/python3 -I
# overlap_mc.py — seguimiento 2026-11-03 (13), punto B: MC del banco (bench_mc.py, muestra de stage5: cos³θ, CAD, semilla 42, trozos de 2e6, 3e7
# rayos) para pares de barras adyacentes de hod2018 cruzados por el mismo muón.
# Selección (como overlap_ratio.py en datos): trayectoria del trigger (T1 ∧ T2); etiqueta = alguna barra de la coordenada ortogonal con camino > 0;
# sonda = exactamente dos barras cruzadas (camino > 0) en la coordenada y son vecinas en el orden de posición (columnas j, j+1 de bar_paths,
# capas alternas, solapamiento de 1 mm entre c_j + 1.5 y c_{j+1} − 1.5).
# Por par: L1 = camino en la barra de menor posición, L2 = en la otra (mm); u_mid = coordenada medida en el plano medio entre las dos capas
# menos el centro del solapamiento (c_j + c_{j+1})/2 (mm); u1, u2 = la misma en el plano medio de cada capa; j; cosθ.
# Los dos órdenes de bloque (X arriba / Y arriba) por separado; dentro de un bloque la capa A está arriba (layer_z; el cociente L1/(L1+L2) no
# depende de cuál capa está arriba salvo por el signo de u, que se da por capa).
# Salida: <out>/overlap_mc.npz (arrays por top y coordenada) y overlap_mc_meta.json. Uso: overlap_mc.py <out dir> [N]
import json, os, sys, time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import bench_mc as bm  # noqa: E402

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 30_000_000
os.makedirs(OUT, exist_ok=True)
G = bm.GEOMS["CAD"]
CENTRES = sorted([(c, "A") for c in bm.LAYER_A] + [(c, "B") for c in bm.LAYER_B])
acc = {(top, co): {k: [] for k in ("L1", "L2", "u_mid", "u1", "u2", "j", "cos")} for top in "XY" for co in "XY"}
ntrig, ntag = 0, {(top, co): 0 for top in "XY" for co in "XY"}
rng = np.random.default_rng(42)
t0 = time.time()
for start in range(0, N, 2_000_000):
    m = min(2_000_000, N - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, G)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    ntrig += int(t.sum())
    for top in "XY":
        pth = bm.bar_paths(rt, G["H18"], top)
        zl = bm.layer_z(G["H18"], top)
        for co, ref in (("X", "Y"), ("Y", "X")):
            P = pth[co]
            tag = (pth[ref] > 0).any(1)
            ntag[(top, co)] += int(tag.sum())
            hit = P > 0
            two = tag & (hit.sum(1) == 2)
            j = np.argmax(hit, axis=1)
            adj = two & (j < 15) & hit[np.arange(len(j)), np.minimum(j + 1, 15)]
            idx = np.nonzero(adj)[0]
            jj = j[idx]
            L1, L2 = P[idx, jj], P[idx, jj + 1]
            cj = np.array([c for c, _ in CENTRES])
            uc = 0.5 * (cj[jj] + cj[jj + 1])
            za = np.array([zl[(co, CENTRES[k][1])] for k in jj])
            zb = np.array([zl[(co, CENTRES[k + 1][1])] for k in jj])
            sub = bm.Rays(rt.x0[idx], rt.y0[idx], rt.dx[idx], rt.dy[idx], rt.dz[idx])

            def meas(z):
                x, y = sub.prop(z)
                return x if co == "X" else y

            a = acc[(top, co)]
            a["L1"].append(L1); a["L2"].append(L2); a["j"].append(jj)
            a["u_mid"].append(meas(0.5 * (za + zb)) - uc); a["u1"].append(meas(za) - uc); a["u2"].append(meas(zb) - uc)
            a["cos"].append(-sub.dz)
out = {}
for (top, co), a in acc.items():
    for k, v in a.items():
        out[f"{top}|{co}|{k}"] = np.concatenate(v)
np.savez_compressed(f"{OUT}/overlap_mc.npz", **out)
meta = {"rays": N, "triggers": ntrig, "seed": 42, "model": "cos3", "geometry": "CAD", "z_h18_mm": G["H18"], "seconds": time.time() - t0,
        "tag_events": {f"{t}|{c}": v for (t, c), v in ntag.items()},
        "pairs": {f"{t}|{c}": int(len(out[f"{t}|{c}|L1"])) for t, c in acc}}
json.dump(meta, open(f"{OUT}/overlap_mc_meta.json", "w"), indent=1)
print(json.dumps(meta))
