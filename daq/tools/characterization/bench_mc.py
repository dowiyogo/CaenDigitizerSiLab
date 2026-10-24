# bench_mc.py — PROC-1: MC geométrico del banco de rayos cósmicos (rayos rectos, sin interacciones), compartido por
# geom_sensitivity.py y los scripts de daq/tools/proceeding/ (antes el generador estaba dentro de geom_sensitivity.py).
#  - generate(): puntos uniformes en un plano horizontal a Z_START (±200 mm) y dirección con pdf de cosθ ∝ c² ("cos2", el
#    taller) o c³ ("cos3", intensidad cos²θ con el factor de proyección); mismo orden de llamadas al generador que antes.
#  - Geometrías (z hacia arriba, origen en el centro de hod2018):
#      OLD: analysis/python/geom_acceptance.py: T1 z = 124.1 mm, T2 z = 134.1 mm (planos medios), hod2019 z = 98.4 mm;
#      CAD (René, 2026-10-24): módulo de 15 mm de alto, la placa inferior 14.7 mm sobre la cara superior del módulo de arriba
#           -> T1 entre 120.6 y 130.6 mm (plano medio 125.6), T2 encima (plano medio 135.6); placas de 150 × 150 × 10 mm
#           centradas en (−42.5, −42.5) (T1) y (+42.5, +42.5) (T2).
#  - Barras (CAD, René): 75 × 3 × 1 mm, paso 4.0 mm, 8 por capa, centros de capa en −15, −11, …, 13 y −13, −9, …, 15 mm;
#    3.0 mm entre centros de las dos capas de una coordenada; 6.0 mm entre centros de los bloques X e Y; bloques centrados
#    en el plano medio del módulo. Las barras X miden x (3 mm en x, 75 mm en y); las Y miden y.
import math

import numpy as np

Z_START = 160.0
HALF = 16.5
TH = 75.0
T1C, T2C = (-42.5, -42.5), (42.5, 42.5)
GEOMS = {"OLD": {"T1": 124.1, "T2": 134.1, "H19": 98.4, "H18": 0.0},
         "CAD": {"T1": 125.6, "T2": 135.6, "H19": 98.4, "H18": 0.0}}
LAYER_A = np.arange(-15.0, 14.0, 4.0)  # −15 … 13
LAYER_B = np.arange(-13.0, 16.0, 4.0)  # −13 … 15
BAR_W, BAR_T, BAR_L = 3.0, 1.0, 75.0
LAYER_DZ, BLOCK_DZ = 3.0, 6.0


class Rays:
    def __init__(self, x0, y0, dx, dy, dz):
        self.x0, self.y0, self.dx, self.dy, self.dz = x0, y0, dx, dy, dz

    def prop(self, z):
        t = (z - Z_START) / self.dz
        return self.x0 + self.dx * t, self.y0 + self.dy * t

    def cos_theta(self):
        return -self.dz


def generate(rng, m, model="cos3"):
    U = rng.uniform(0, 1, m)
    c = U ** (1 / 3) if model == "cos2" else U ** (1 / 4)  # pdf de cosθ ∝ c² o c³
    th = np.arccos(c)
    ph = rng.uniform(0, 2 * np.pi, m)
    dx, dy, dz = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), -np.cos(th)
    x0, y0 = rng.uniform(-200, 200, m), rng.uniform(-200, 200, m)
    return Rays(x0, y0, dx, dy, dz)


def trigger(r, g):
    p2, p1 = r.prop(g["T2"]), r.prop(g["T1"])
    return (np.abs(p2[0] - T2C[0]) <= TH) & (np.abs(p2[1] - T2C[1]) <= TH) & (np.abs(p1[0] - T1C[0]) <= TH) & (np.abs(p1[1] - T1C[1]) <= TH)


def in_square(r, z):
    x, y = r.prop(z)
    return (np.abs(x) <= HALF) & (np.abs(y) <= HALF)


def layer_z(z_module, top="X"):
    """Planos medios de las 4 capas de una estación: {(coordenada, capa): z}. El bloque `top` está arriba; dentro de cada
    bloque la capa A está arriba (no hay registro de cuál capa está arriba: el efecto se da aparte)."""
    zx = z_module + (BLOCK_DZ / 2 if top == "X" else -BLOCK_DZ / 2)
    zy = z_module + (-BLOCK_DZ / 2 if top == "X" else BLOCK_DZ / 2)
    return {("X", "A"): zx + LAYER_DZ / 2, ("X", "B"): zx - LAYER_DZ / 2, ("Y", "A"): zy + LAYER_DZ / 2, ("Y", "B"): zy - LAYER_DZ / 2}


def bar_paths(r, z_module, top="X"):
    """Longitud de camino (mm) en cada barra: {"X": array (n, 16), "Y": array (n, 16)}, columnas = barras ordenadas por
    centro (−15 … 15, alternando capas). Una barra X ocupa |x − xc| ≤ 1.5, |y| ≤ 37.5 y |z − zl| ≤ 0.5."""
    out = {}
    zl = layer_z(z_module, top)
    cos = np.abs(r.dz)
    for coord in ("X", "Y"):
        cols = []
        centres = sorted([(c, "A") for c in LAYER_A] + [(c, "B") for c in LAYER_B])
        for c, layer in centres:
            z0 = zl[(coord, layer)]
            za, zb = z0 - BAR_T / 2, z0 + BAR_T / 2
            ua, va = r.prop(za)
            ub, vb = r.prop(zb)
            u_a, u_b = (ua, ub) if coord == "X" else (va, vb)   # coordenada medida
            w_a, w_b = (va, vb) if coord == "X" else (ua, ub)   # coordenada a lo largo de la barra
            lo, hi = c - BAR_W / 2, c + BAR_W / 2
            # fracción del tramo [za, zb] con la coordenada medida dentro de la barra (movimiento lineal en z)
            du = u_b - u_a
            with np.errstate(divide="ignore", invalid="ignore"):
                s1 = np.where(np.abs(du) > 1e-12, (lo - u_a) / du, np.where((u_a >= lo) & (u_a <= hi), -np.inf, np.inf))
                # du = 0 (rayo sin desplazamiento en la coordenada medida): dentro -> [0, 1]; fuera -> s1 = s2 = +inf, fracción 0
                s2 = np.where(np.abs(du) > 1e-12, (hi - u_a) / du, np.inf)
            smin, smax = np.maximum(0, np.minimum(s1, s2)), np.minimum(1, np.maximum(s1, s2))
            frac = np.clip(smax - smin, 0, 1)
            along = np.abs(0.5 * (w_a + w_b)) <= BAR_L / 2
            cols.append(np.where(along, frac * BAR_T / cos, 0.0))
        out[coord] = np.stack(cols, axis=1)
    return out


def rel_err(k, n):
    p = k / n
    return p, math.sqrt(max(p * (1 - p), 0) / n)
