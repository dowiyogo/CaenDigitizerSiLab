#!/usr/bin/python3 -I
# ray_export.py — exporta, para Geant4, los impactos seleccionados de UNA barra sonda de hod2018 del MC geométrico del banco de PROC
# (daq/tools/characterization/bench_mc.py, sin cambios): geometría CAD (T1/T2 de mayo), cos³θ (intensidad cos²θ × proyección), puntos
# uniformes en ±200 mm a z = 160 mm, trigger T1 ∧ T2, hod2018 en z = 0. Impacto = camino > 0 en la barra sonda y > 0 en alguna barra de la
# coordenada ortogonal (la selección de analysis/followup5_2026_10_31/path_length.py). Cada fila es un rayo en el sistema de la barra
# (convenciones en la cabecera del CSV). path_mm es la cuerda de bench_mc.bar_paths(); se comprueba contra la cuerda recalculada aquí.
# Uso: ray_export.py <salida.csv> [--coord Y] [--layer A] [--centre -15] [--top X] [--rows 200000] [--seed 42]
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "characterization"))
import bench_mc as bm  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--coord", default="Y", choices=("X", "Y"))
ap.add_argument("--layer", default="A", choices=("A", "B"))
ap.add_argument("--centre", type=float, default=-15.0)
ap.add_argument("--top", default="X", choices=("X", "Y"))
ap.add_argument("--rows", type=int, default=200_000)
ap.add_argument("--seed", type=int, default=42)
a = ap.parse_args()

centres = sorted([(c, "A") for c in bm.LAYER_A] + [(c, "B") for c in bm.LAYER_B])
col = centres.index((a.centre, a.layer))  # falla si la barra no existe
oth = "X" if a.coord == "Y" else "Y"
g = bm.GEOMS["CAD"]
zc = bm.layer_z(g["H18"], a.top)[(a.coord, a.layer)]
ztop, zbot = zc + bm.BAR_T / 2, zc - bm.BAR_T / 2
rng = np.random.default_rng(a.seed)
CH = 2_000_000
rows, rays, trig, nmax_check = [], 0, 0, 0.0
while sum(len(r[0]) for r in rows) < a.rows:
    r = bm.generate(rng, CH, "cos3")
    rays += CH
    t = bm.trigger(r, g)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    trig += int(t.sum())
    p = bm.bar_paths(rt, g["H18"], a.top)
    sel = (p[a.coord][:, col] > 0) & (p[oth] > 0).any(1)
    if not sel.any():
        continue
    rs = bm.Rays(rt.x0[sel], rt.y0[sel], rt.dx[sel], rt.dy[sel], rt.dz[sel])
    path = p[a.coord][sel, col]
    # sistema de la barra: u = coordenada medida − centro (ancho de 3 mm), v = coordenada a lo largo (75 mm), w = z − zc (alto de 1 mm)
    xt, yt = rs.prop(ztop)
    xb, yb = rs.prop(zbot)
    ut, ub = (xt, xb) if a.coord == "X" else (yt, yb)
    vt, vb = (yt, yb) if a.coord == "X" else (xt, xb)
    ut, ub = ut - a.centre, ub - a.centre
    du = ub - ut
    h = bm.BAR_W / 2
    # parámetro s del tramo de arriba (s = 0, w = +0.5) hacia abajo (s = 1, w = −0.5) en que |u| ≤ 1.5
    with np.errstate(divide="ignore", invalid="ignore"):
        s1, s2 = (-h - ut) / du, (h - ut) / du
    s_in = np.where(np.abs(du) > 1e-12, np.maximum(0, np.minimum(s1, s2)), 0.0)
    s_out = np.where(np.abs(du) > 1e-12, np.minimum(1, np.maximum(s1, s2)), 1.0)
    cos = np.abs(rs.dz)
    chord = (s_out - s_in) * bm.BAR_T / cos
    nmax_check = max(nmax_check, float(np.max(np.abs(chord - path))))
    top_in = s_in <= 1e-12
    bot_out = s_out >= 1 - 1e-12
    u_e = np.where(top_in, ut, np.clip(ut + s_in * du, -h, h))
    v_e = vt + s_in * (vb - vt)
    w_e = bm.BAR_T / 2 - s_in * bm.BAR_T
    entry = np.where(top_in, "top", np.where(u_e > 0, "side+u", "side-u"))
    u_x = ut + s_out * du
    exitf = np.where(bot_out, "bottom", np.where(u_x > 0, "side+u", "side-u"))
    theta = np.degrees(np.arccos(cos))
    # dirección de propagación (hacia abajo) proyectada en el plano horizontal; acimut desde +v hacia +u
    du_dir, dv_dir = (rs.dx, rs.dy) if a.coord == "X" else (rs.dy, rs.dx)
    phi = np.degrees(np.arctan2(du_dir, dv_dir)) % 360.0
    rows.append((u_e, v_e, w_e, entry, exitf, theta, phi, path))
cat = [np.concatenate([r[i] for r in rows])[: a.rows] for i in range(8)]
n = len(cat[0])
fr_side_in = float(np.mean(cat[3] != "top"))
fr_side_out = float(np.mean(cat[4] != "bottom"))
fr_side_any = float(np.mean((cat[3] != "top") | (cat[4] != "bottom")))
v_out = int(np.sum(np.abs(cat[1]) > bm.BAR_L / 2))
with open(a.out, "w") as f:
    f.write(f"# ray_export.py: impactos seleccionados de una barra sonda de hod2018, MC del banco de PROC (bench_mc.py), para Geant4.\n"
            f"# MC: geometría CAD (T1 z = {g['T1']} mm centrada en (-42.5,-42.5), T2 z = {g['T2']} mm en (+42.5,+42.5), 150 x 150 mm), "
            f"cos^3 en cos(theta) (I ~ cos^2 theta por ángulo sólido x proyección), puntos uniformes en +-200 mm a z = 160 mm, semilla {a.seed}; "
            f"{rays} rayos, {trig} triggers; se guardan los primeros {n} impactos.\n"
            f"# Selección: trigger T1^T2, camino > 0 en la barra sonda y > 0 en alguna barra {oth} (cualquier capa). Bloque {a.top} arriba.\n"
            f"# Barra sonda: coordenada {a.coord}, capa {a.layer}, centro {a.centre:+.1f} mm en el banco; z del centro = {zc:+.2f} mm "
            f"(hod2018 en z = 0); 3 mm (u) x 75 mm (v) x 1 mm (w).\n"
            f"# Sistema de la barra ({'dextrógiro' if a.coord == 'X' else 'levógiro: (u, v, w) = (y, x, z) del banco'}): u ={a.coord.lower()}_banco - ({a.centre:+.1f}) mm (a través del ancho, [-1.5, 1.5]); "
            f"v = {'y' if a.coord == 'X' else 'x'}_banco (a lo largo, eje largo; [-37.5, 37.5]); w = z_banco - zc (vertical, hacia arriba; caras en w = +-0.5).\n"
            f"# u_mm, v_mm, w_mm: punto de entrada. entry = top (cara de arriba, w = +0.5) o side+u / side-u (cara lateral de 1 mm en u = +-1.5; "
            f"w_mm da la altura). exit = bottom (w = -0.5) o side+u / side-u.\n"
            f"# theta_deg: cenit (0 = vertical). phi_deg: acimut de la dirección de propagación (hacia abajo) proyectada en el plano u-v, "
            f"medido desde +v hacia +u, en [0, 360).\n"
            f"# path_mm: cuerda de PROC (bench_mc.bar_paths), igual a la cuerda por la caja de sección 3 x 1 mm (|diferencia| máx. {nmax_check:.2e} mm). "
            f"PROC no modela las caras de los extremos (acepta si el punto medio tiene |v| <= 37.5): {v_out} filas con |v| de entrada > 37.5.\n"
            f"# v de entrada en [{cat[1].min():.3f}, {cat[1].max():.3f}] mm: la barra ortogonal exigida limita el impacto al cuadrado activo (+-16.5 mm).\n"
            f"# Fracciones: entrada lateral {fr_side_in:.4f}, salida lateral {fr_side_out:.4f}, alguna cara lateral {fr_side_any:.4f}.\n")
    f.write("u_mm,v_mm,w_mm,entry,exit,theta_deg,phi_deg,path_mm\n")
    for i in range(n):
        f.write(f"{cat[0][i]:.5f},{cat[1][i]:.5f},{cat[2][i]:.5f},{cat[3][i]},{cat[4][i]},{cat[5][i]:.5f},{cat[6][i]:.5f},{cat[7][i]:.6f}\n")
th, pa = cat[5], cat[7]
print(f"rays {rays}, triggers {trig}, rows {n}; z_centre {zc:+.2f} mm; max |chord - path| {nmax_check:.2e} mm")
print(f"theta mean {th.mean():.2f} deg, median {np.median(th):.2f} deg; path median {np.median(pa):.4f} mm, mean {pa.mean():.4f} mm")
print(f"entry side {fr_side_in:.4f}, exit side {fr_side_out:.4f}, any side {fr_side_any:.4f}; |v| > 37.5 at entry: {v_out}")
print(f"u range [{cat[0].min():.4f}, {cat[0].max():.4f}], v range [{cat[1].min():.3f}, {cat[1].max():.3f}]")
