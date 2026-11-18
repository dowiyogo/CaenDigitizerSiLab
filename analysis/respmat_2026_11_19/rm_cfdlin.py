#!/usr/bin/python3 -I
# rm_cfdlin.py — seguimiento 2026-11-19 (matriz de respuesta), comprobación de la alineación de la plantilla: la parte fraccionaria del tiempo CFD de los pulsos
# aislados no es uniforme (el flanco de subida dura menos de una muestra de 16 ns); se compara el promedio alineado con el t_CFD medido y con t_CFD linealizado
# (fase → su distribución acumulada, uniforme por construcción). Misma selección que rm_prep.py (mayo en tiempo bueno sin filtrar aquí: todos los eventos).
# Uso: rm_cfdlin.py <out dir>   (imprime por canal el pico, su posición, el FWHM y el área del promedio; resumen de |Δpico|/pico)
import sys

import numpy as np

OUT = sys.argv[1]
ISD = np.dtype([("ev", "<u4"), ("ch", "u1"), ("isod", "u1"), ("k", "<u2"), ("a", "<f4"), ("t", "<f4"), ("y", "<f4", 18)])
for camp, pre in (("may", "may_v8"), ("oct", "oct_run8")):
    G = np.array([float(x) for x in open(f"{OUT}/work/G_{camp}.txt").read().split(",")])
    r = np.fromfile(f"{OUT}/work/{pre}_iso.bin", dtype=ISD)
    rel, frac = [], np.zeros(8)
    for c in range(32):
        if c == 3:
            continue
        s = r[r["ch"] == c]
        dk = s["k"].astype(float) - s["t"]
        q = (dk > 0) & (dk <= 2.5) & (np.abs(s["y"][:, :6]).max(1) < 0.5 * G[c])
        s, dk = s[q], dk[q]
        u = s["a"] / G[c]
        one = (u > 0.45) & (u < 1.2)
        s = s[one]
        t = s["t"].astype(float)
        fr = t % 1
        frac += np.histogram(fr, bins=8, range=(0, 1))[0]
        frl = np.searchsorted(np.sort(fr), fr) / len(fr)
        pk = {}
        for name, tt in (("raw", t), ("lin", np.floor(t) + frl)):
            tau = (np.arange(18) - 8)[None, :] + (s["k"][:, None] - tt[:, None])
            fb = np.floor((tau + 7.5) / 0.125).astype(int)
            ok = (fb >= 0) & (fb < 140)
            m = np.bincount(fb[ok], weights=s["y"][ok], minlength=140) / np.maximum(np.bincount(fb[ok], minlength=140), 1)
            fc = -7.5 + 0.125 * (np.arange(140) + 0.5)
            sh = m - m[(fc >= -3) & (fc < -1.5)].mean()
            sm = np.convolve(sh, np.ones(3) / 3, "same")
            pk[name] = (sm.max(), fc[np.argmax(sm)], 0.125 * (sm > 0.5 * sm.max()).sum(), sh[fc >= -1.5].sum() * 0.125)
        rel.append(abs(pk["lin"][0] - pk["raw"][0]) / pk["raw"][0])
        print(f"{camp} ch{c}: crudo pico {pk['raw'][0]:.3f} mV en τ {pk['raw'][1]:.4f}, FWHM {pk['raw'][2]:.3f}, área {pk['raw'][3]:.3f}; "
              f"linealizado pico {pk['lin'][0]:.3f} en τ {pk['lin'][1]:.4f}, FWHM {pk['lin'][2]:.3f}, área {pk['lin'][3]:.3f}")
    print(f"{camp}: |Δpico|/pico mediana {np.median(rel):.5f}, máx {max(rel):.5f}; parte fraccionaria de t_CFD en 8 bins (todos los canales): {frac.astype(int).tolist()}")
