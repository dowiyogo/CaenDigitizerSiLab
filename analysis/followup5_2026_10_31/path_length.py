#!/usr/bin/python3 -I
# path_length.py — seguimiento 2026-10-31 (5), punto 2: longitud de camino en la barra sonda de hod2018 para los impactos que entran en
# el espectro de coincidencias X–Y (en el MC: la trayectoria cruza la barra sonda y al menos una barra de la coordenada ortogonal), con
# el MC geométrico del banco de PROC (daq/tools/characterization/bench_mc.py), sin cambiar nada de él: geometría CAD (placas T1/T2 de
# mayo, las "viejas", 150 × 150 mm centradas en ∓42.5 mm), ley angular "cos3" (pdf de cosθ ∝ cos³θ: intensidad cos²θ por el factor de
# proyección en un plano horizontal), puntos uniformes en ±200 mm a z = 160 mm, hod2018 en z = 0; los dos órdenes de bloque (X o Y
# arriba). Un impacto = (trayectoria, barra sonda con camino > 0) con alguna barra ortogonal con camino > 0. Se dan también la variante
# con camino mínimo 0.2 mm en ambas barras (la banda intermedia de la figura de dobles de PROC) y la versión "una barra por coordenada"
# (camino de la barra con más camino en la coordenada sonda).
# Uso: path_length.py <salida.json> [N] [semilla]
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "daq", "tools", "characterization"))
import bench_mc as bm  # noqa: E402

OUTF = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 50_000_000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 42
CHUNK = 2_000_000
BW = 0.01
NB = 1000  # 0 … 10 mm
LAYER_OF_COL = ["A" if i % 2 == 0 else "B" for i in range(16)]  # centros ordenados −15 (A), −13 (B), …


def stats(h, extra=None):
    n = h.sum()
    c = np.cumsum(h)
    x = (np.arange(NB) + 0.5) * BW

    def q(p):
        return float((np.searchsorted(c, p * n) + 0.5) * BW)

    d = {"n": int(n), "mean_mm": float((h * x).sum() / n), "median_mm": q(0.5), "q05_mm": q(0.05), "q16_mm": q(0.16), "q84_mm": q(0.84),
         "q95_mm": q(0.95), "frac_lt_0p5mm": float(h[: int(0.5 / BW)].sum() / n), "frac_lt_1mm": float(h[: int(1.0 / BW)].sum() / n)}
    d["median_over_1mm"] = d["median_mm"] / 1.0
    d["mean_over_1mm"] = d["mean_mm"] / 1.0
    if extra:
        d.update(extra)
    return d


def main():
    g = bm.GEOMS["CAD"]
    rng = np.random.default_rng(SEED)
    keys = [(top, cut, coord, sel) for top in ("X", "Y") for cut in (0.0, 0.2) for coord in ("X", "Y") for sel in ("A", "B", "all", "maxbar")]
    H = {k: np.zeros(NB, np.int64) for k in keys}
    sec = {k: [0.0, 0] for k in keys}  # suma de 1/cosθ y cuenta (impactos)
    ntrig = 0
    for start in range(0, N, CHUNK):
        m = min(CHUNK, N - start)
        r = bm.generate(rng, m, "cos3")
        t = bm.trigger(r, g)
        rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
        ntrig += int(t.sum())
        inv = 1.0 / np.abs(rt.dz)
        for top in ("X", "Y"):
            p = bm.bar_paths(rt, g["H18"], top)
            for cut in (0.0, 0.2):
                for coord, oth in (("X", "Y"), ("Y", "X")):
                    po = p[oth]
                    has_o = (po > cut).any(1) if cut > 0 else (po > 0).any(1)
                    pc = p[coord]
                    hit = (pc > cut) if cut > 0 else (pc > 0)
                    hit &= has_o[:, None]
                    for sel in ("A", "B", "all"):
                        cols = [i for i in range(16) if sel == "all" or LAYER_OF_COL[i] == sel]
                        v = pc[:, cols][hit[:, cols]]
                        w = np.broadcast_to(inv[:, None], pc.shape)[:, cols][hit[:, cols]]
                        H[(top, cut, coord, sel)] += np.bincount(np.minimum((v / BW).astype(int), NB - 1), minlength=NB)
                        sec[(top, cut, coord, sel)][0] += float(w.sum())
                        sec[(top, cut, coord, sel)][1] += int(v.size)
                    ev = hit.any(1)
                    v = pc[ev].max(1)
                    H[(top, cut, coord, "maxbar")] += np.bincount(np.minimum((v / BW).astype(int), NB - 1), minlength=NB)
                    sec[(top, cut, coord, "maxbar")][0] += float(inv[ev].sum())
                    sec[(top, cut, coord, "maxbar")][1] += int(ev.sum())
    res = {"geometry": "CAD", "model": "cos3", "rays": N, "seed": SEED, "triggers": ntrig, "bin_mm": BW,
           "angle_law": "bench_mc.generate(model='cos3'): cosθ = U^(1/4), pdf de cosθ ∝ cos³θ (I ∝ cos²θ × proyección), φ uniforme, "
                        "(x0, y0) uniformes en ±200 mm a z = 160 mm; trigger = cruce de T1 y T2 (bench_mc.trigger)",
           "results": {}}
    for k in keys:
        top, cut, coord, sel = k
        res["results"][f"top{top}_cut{cut:g}_{coord}_{sel}"] = stats(H[k], {"mean_sec_theta": sec[k][0] / max(sec[k][1], 1)})
    res["histograms"] = {f"top{t}_cut{c:g}_{co}_{s}": H[(t, c, co, s)].tolist() for (t, c, co, s) in keys if s in ("all", "maxbar")}
    with open(OUTF, "w") as f:
        json.dump(res, f, indent=1)
    print(f"triggers {ntrig} de {N} rayos")
    print("clave n mediana_mm media_mm q16 q84 frac<0.5mm frac<1mm <secθ> mediana/1mm")
    for k, d in res["results"].items():
        print(k, d["n"], f"{d['median_mm']:.3f} {d['mean_mm']:.4f} {d['q16_mm']:.3f} {d['q84_mm']:.3f} {d['frac_lt_0p5mm']:.4f} "
              f"{d['frac_lt_1mm']:.4f} {d['mean_sec_theta']:.4f} {d['median_over_1mm']:.3f}")


if __name__ == "__main__":
    main()
