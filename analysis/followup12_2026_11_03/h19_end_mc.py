#!/usr/bin/python3 -I
# h19_end_mc.py — seguimiento 2026-11-03 (12b), punto 1-2: prueba geométrica de la hipótesis «hod2019 solo detecta cerca del extremo con MPPC».
# MC del banco de stage5 (bench_mc.py: cos³θ, placas del CAD, rayos rectos, semilla 42, trozos de 2e6; los primeros 3e7 rayos son la muestra de
# stage5, se sigue hasta N = 1e8 como sensitive_length.py). hod2019 en z = 98.4 mm (módulo), misma disposición de barras que hod2018: 75 × 3 × 1 mm,
# MPPC en un extremo, ±37.5 mm; las dos capas de una coordenada se leen por extremos opuestos.
# Extremos de lectura (registro de conexiones, kCfg_2026_05 / handoff §5): X tira k par (capa A) = cara Delantera, impar (capa B) = Trasera;
# Y par (A) = Opuesta, impar (B) = Frontal. Qué cara está en +y (X) o en +x (Y) no está registrado: se dan las dos asignaciones de cada coordenada
# ("A+" = capa A con su MPPC en +37.5, B en −37.5; "A-" al revés) y los dos órdenes de bloque (X arriba, Y arriba).
# Para cada trayectoria del trigger (T1 ∧ T2) y cada barra cruzada (camino > 0, como bar_paths): w = coordenada a lo largo de la barra en el plano
# medio de su capa; d = distancia de w al extremo con MPPC. d_X = mínimo sobre las barras X cruzadas (∞ si ninguna), igual d_Y.
#   P_X(z) = P(d_X ≤ z), P_Y(z) = P(d_Y ≤ z), P_XY(z_X, z_Y) = P(d_X ≤ z_X ∧ d_Y ≤ z_Y), P_XY18 = además cruza el cuadrado activo de hod2018
#   (|x|, |y| ≤ 16.5 en z = 0, bench_mc.in_square).
# Por barra (columnas = barras ordenadas por centro, −15 … 15): cruces con d ≤ z para z en ZB (una trayectoria puede contar en dos barras).
# Salida: <out>/h19_end_mc.npz (histogramas 2D de (d_X, d_Y) en bins de 0.25 mm hasta 75 mm + bin «sin cruce», todos y con hod2018; bars|…) y
#         <out>/h19_end_mc_meta.json. Uso: h19_end_mc.py <out dir> [N]
import json, os, sys, time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import bench_mc as bm  # noqa: E402

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000_000
os.makedirs(OUT, exist_ok=True)
G = bm.GEOMS["CAD"]
HL = bm.BAR_L / 2
STEP = 0.25
EDG = np.append(np.arange(0, HL * 2 + STEP, STEP), 1e9)  # último bin: sin cruce
NB = len(EDG) - 1
CENTRES = sorted([(c, "A") for c in bm.LAYER_A] + [(c, "B") for c in bm.LAYER_B])
H = {}
for top in "XY":
    for cx in ("A+", "A-"):
        for cy in ("A+", "A-"):
            for s in ("all", "h18"):
                H[(top, cx, cy, s)] = np.zeros((NB, NB))
ZB = np.array([2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0])
PB = {(top, coord, cv): np.zeros((len(ZB), 16)) for top in "XY" for coord in "XY" for cv in ("A+", "A-")}  # cruces por barra con d ≤ z
trig = 0
trig5 = None
rng = np.random.default_rng(42)
t0 = time.time()


def dmin(rt, zl, coord, conv, path, pb=None):
    """d (mm) al extremo con MPPC, mínimo sobre las barras cruzadas de la coordenada."""
    best = np.full(len(rt.dz), np.inf)
    for j, (c, layer) in enumerate(CENTRES):
        x, y = rt.prop(zl[(coord, layer)])
        w = y if coord == "X" else x
        plus = (layer == "A") == (conv == "A+")
        d = HL - w if plus else w + HL
        hit = path[:, j] > 0
        best = np.where(hit, np.minimum(best, d), best)
        if pb is not None:
            pb[:, j] += ((d[None, :] <= ZB[:, None]) & hit[None, :]).sum(1)
    return best


for start in range(0, N, 2_000_000):
    m = min(2_000_000, N - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, G)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    trig += int(t.sum())
    s18 = bm.in_square(rt, G["H18"])
    for top in "XY":
        pth = bm.bar_paths(rt, G["H19"], top)
        zl = bm.layer_z(G["H19"], top)
        dX = {cv: dmin(rt, zl, "X", cv, pth["X"], PB[(top, "X", cv)]) for cv in ("A+", "A-")}
        dY = {cv: dmin(rt, zl, "Y", cv, pth["Y"], PB[(top, "Y", cv)]) for cv in ("A+", "A-")}
        for cx in ("A+", "A-"):
            for cy in ("A+", "A-"):
                a, b = np.minimum(dX[cx], 1e8), np.minimum(dY[cy], 1e8)
                H[(top, cx, cy, "all")] += np.histogram2d(a, b, bins=(EDG, EDG))[0]
                H[(top, cx, cy, "h18")] += np.histogram2d(a[s18], b[s18], bins=(EDG, EDG))[0]
    if start + m == 30_000_000:
        trig5 = trig
np.savez_compressed(f"{OUT}/h19_end_mc.npz", edges=EDG, zb=ZB, **{"|".join(k): v for k, v in H.items()},
                    **{"bars|" + "|".join(k): v for k, v in PB.items()})
json.dump({"rays": N, "triggers": trig, "triggers_first_3e7_stage5": trig5, "seed": 42, "model": "cos3", "geometry": "CAD",
           "z_h19_mm": G["H19"], "bin_mm": STEP, "seconds": time.time() - t0,
           "conventions": "A+: capa A (tiras pares: X Delantera, Y Opuesta) con MPPC en +37.5 (y para X, x para Y), capa B en −37.5; A-: al revés"},
          open(f"{OUT}/h19_end_mc_meta.json", "w"), indent=1)
print(f"rayos {N}, triggers {trig} (primeros 3e7: {trig5}), {time.time() - t0:.0f} s")
