#!/usr/bin/python3 -I
# rm_tcut.py — seguimiento 2026-11-19 (matriz de respuesta), variante "tcut" de la distribución de tiempos de inyección: la distribución neta de floor(t_CFD) de
# los muones etiquetados (work/<camp>_tdist.bin de rm_prep.py) restringida a las muestras enteras 48-55 y renormalizada. Fuera de ese intervalo el pico de la
# plantilla cae fuera de la ventana signal [49,58) (window_edge.csv de rm_core.py), mientras que los muones etiquetados tienen su pulso de aS dentro de la
# ventana: los t_CFD de 43-47 y 56 vienen de la búsqueda hacia atrás del CFD (hasta 6 muestras) sobre una línea base que sube, no del tiempo de llegada.
# Salida: work/<camp>_tcut_tdist.bin (float64 [32][144]); imprime la fracción quitada por capa.
# Uso: rm_tcut.py <out dir>
import sys

import numpy as np

OUT = sys.argv[1]
for camp in ("may", "oct"):
    t = np.fromfile(f"{OUT}/work/{camp}_tdist.bin", "<f8").reshape(32, 144)
    u = np.zeros_like(t)
    for c in range(32):
        if t[c].sum() > 0:
            u[c, 48:56] = t[c, 48:56]
            print(f"{camp} ch{c}: fracción quitada {1 - u[c].sum() / t[c].sum():.4f}")
            u[c] /= u[c].sum()
    u.astype("<f8").tofile(f"{OUT}/work/{camp}_tcut_tdist.bin")
