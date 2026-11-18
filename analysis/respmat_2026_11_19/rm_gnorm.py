#!/usr/bin/python3 -I
# rm_gnorm.py — seguimiento 2026-11-19 (matriz de respuesta), variante "gnorm" de la plantilla: la plantilla de rm_prep.py (normalizada al pico de 1 p.e. de la
# amplitud por mínimos cuadrados, "ls") reescalada para que el estimador de pico sobre línea base plana y fase uniforme dé en promedio G (la separación de picos
# de la calibración de los espectros) por p.e.: factor = G / E_φ[max_k T(k − φ)] (k entero, φ uniforme en [0,1), 400 fases, sin ruido ni cuantización).
# Salida: <out>/work/<camp>_gnorm_tmpl.bin (float64 [32][333]), <out>/gnorm.json (E/G y factor por canal).
# Uso: rm_gnorm.py <out dir>
import json, sys

import numpy as np

OUT = sys.argv[1]
grid = -1.5 + 0.125 * np.arange(333)
res = {}
for camp in ("may", "oct"):
    G = np.array([float(v) for v in open(f"{OUT}/work/G_{camp}.txt").read().split(",")])
    T = np.fromfile(f"{OUT}/work/{camp}_tmpl.bin", "<f8").reshape(32, 333)
    ph = (np.arange(400) + 0.5) / 400
    ks = np.arange(-2, 8)[None, :] + ph[:, None]
    Tg = np.zeros_like(T)
    res[camp] = {}
    for c in range(32):
        if not T[c].any():
            continue
        e = np.interp(ks, grid, T[c], left=0, right=0).max(1)
        f = G[c] / e.mean()
        Tg[c] = T[c] * f
        res[camp][str(c)] = {"E_max_sampled_over_G": float(e.mean() / G[c]), "min_over_G": float(e.min() / G[c]), "max_over_G": float(e.max() / G[c]), "scale": float(f)}
    Tg.astype("<f8").tofile(f"{OUT}/work/{camp}_gnorm_tmpl.bin")
    sc = np.array([v["scale"] for v in res[camp].values()])
    em = np.array([v["E_max_sampled_over_G"] for v in res[camp].values()])
    print(f"{camp}: E[máx muestreado]/G mediana {np.median(em):.4f} (rango {em.min():.4f}-{em.max():.4f}); factor mediana {np.median(sc):.4f}")
json.dump(res, open(f"{OUT}/gnorm.json", "w"), indent=1)
