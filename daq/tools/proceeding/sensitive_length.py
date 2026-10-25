#!/usr/bin/python3 -I
# sensitive_length.py — PROC-3 etapa 2 (MC, ESTIMADO con un modelo): longitud sensible de las barras de hod2019.
# Modelo (supuestos, sin registro): cada barra detecta un muón solo a menos de L de su extremo con MPPC (eficiencia 1 dentro, 0 fuera);
# las dos capas de una coordenada se leen por extremos opuestos: X capa A en el extremo +y, capa B en −y; Y capa A en +x, capa B en −x.
# Capa A = coordenadas pares (centros −15, −11, …), B = impares (bench_mc). MC del banco: z del CAD, cos³θ, 1e8 rayos, semilla 42,
# los dos órdenes de bloque (promedio; la mitad de la diferencia como rango).
# Para cada trayectoria del trigger se calcula la L mínima que la detecta en X (L_x*: 37.5 − y en la capa A cruzada, o y + 37.5 en la B;
# la menor de las capas cruzadas) y en Y (L_y*). P(X)(L) = P(L_x* ≤ L), P(Y)(L) = P(L_y* ≤ L), P(X∧Y)(L_X, L_Y) = P(L_x* ≤ L_X, L_y* ≤ L_Y).
# Comparación con el run B (mayo v8, h19_xy.py): L que reproduce rX y L que reproduce rY (interpolación lineal en la malla de 0.25 mm),
# X·Y predichas = N·P(X∧Y)(L_X, L_Y), con rango (los dos órdenes de bloque y las combinaciones (L_X, L_X), (L_Y, L_Y)); compatibilidad
# de Poisson con lo observado.
# Salida: work/sensitive_length.json, work/fig_sensitive_length.json. Uso: sensitive_length.py <dir del paquete> [N]
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
STEP = 0.25
EDGES = np.arange(0, 75 + STEP, STEP)
LGRID = [1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 75]
H = {top: np.zeros((len(EDGES) - 1, len(EDGES) - 1)) for top in ("X", "Y")}  # histograma 2D de (L_x*, L_y*)
trig = 0
rng = np.random.default_rng(42)
for start in range(0, N, 2_000_000):
    m = min(2_000_000, N - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, G)
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    trig += int(t.sum())
    for top in ("X", "Y"):
        pth = bm.bar_paths(rt, G["H19"], top)
        zl = bm.layer_z(G["H19"], top)
        need = {}
        for coord in ("X", "Y"):
            lstar = np.full(len(rt.x0), np.inf)
            for l, cols, sign in (("A", slice(0, 16, 2), +1), ("B", slice(1, 16, 2), -1)):
                hit = (pth[coord][:, cols] > 0).any(1)
                x, y = rt.prop(zl[(coord, l)])
                along = y if coord == "X" else x  # las barras X van a lo largo de y, las Y a lo largo de x
                d = HL - sign * along  # distancia al extremo con MPPC (A en +, B en −)
                lstar = np.where(hit, np.minimum(lstar, d), lstar)
            need[coord] = np.clip(lstar, 0, 1e9)
        hh, _, _ = np.histogram2d(need["X"], need["Y"], bins=[EDGES, EDGES])
        H[top] += hh
# CDF: P(L_x* ≤ L_X, L_y* ≤ L_Y) en los bordes superiores de cada bin
C = {top: np.cumsum(np.cumsum(H[top], 0), 1) / trig for top in H}
PX = {top: np.cumsum(H[top].sum(1)) / trig for top in H}  # P(L_x* ≤ borde) (las trayectorias sin cruce en X tienen L* = inf: fuera)
PY = {top: np.cumsum(H[top].sum(0)) / trig for top in H}
Lup = EDGES[1:]


def at(arr, L):
    return float(np.interp(L, Lup, arr))


def pxy(top, LX, LY):
    i = min(int(round(LX / STEP)) - 1, len(Lup) - 1)
    j = min(int(round(LY / STEP)) - 1, len(Lup) - 1)
    return float(C[top][max(i, 0), max(j, 0)])


def L_for(P, val):
    return float(np.interp(val, P, Lup)) if val <= P[-1] else None


S = {"stage": "PROC-3 etapa 2", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "status": "estimated (modelo)", "rays": N,
     "triggers": trig, "geometry": "CAD", "model": "cos3", "seed": 42, "grid_step_mm": STEP,
     "assumptions": "eficiencia 1 a menos de L del extremo con MPPC y 0 fuera; X capa A leída en +y, capa B en −y; Y capa A en +x, capa B en −x (sin registro)",
     "table": []}
for L in LGRID:
    row = {"L_mm": L}
    for top in ("X", "Y"):
        row[f"PX_{top}top"], row[f"PY_{top}top"], row[f"PXY_{top}top"] = at(PX[top], L), at(PY[top], L), pxy(top, L, L)
    for k in ("PX", "PY", "PXY"):
        row[k] = 0.5 * (row[f"{k}_Xtop"] + row[f"{k}_Ytop"])
    S["table"].append(row)
    print("L %5.1f mm: P(X) %.5f  P(Y) %.5f  P(X∧Y) %.6f" % (L, row["PX"], row["PY"], row["PXY"]))
NB = 143814
MEAS = {"35mV": {"rX": 0.01441, "rY": 0.00934, "N_XY": 3, "N_XY_note": "observadas en el mismo evento (sin restar)"},
        "3mV": {"rX": 0.08768, "rY": 0.07554, "N_XY": 4087, "N_XY_note": "observadas en el mismo evento; incluye el azar (producto 4749.5, q·N restado 129 ± 100 mezclados)",
                "qN_mix": 128.7, "qN_mix_err": 100.0}}
S["comparison"] = {}
for thr, m_ in MEAS.items():
    d = {"measured": m_}
    for top in ("X", "Y"):
        Pavg_X = 0.5 * (PX["X"] + PX["Y"])
        Pavg_Y = 0.5 * (PY["X"] + PY["Y"])
    LX, LY = L_for(Pavg_X, m_["rX"]), L_for(Pavg_Y, m_["rY"])
    d["L_X_mm"], d["L_Y_mm"] = LX, LY
    preds = {}
    for top in ("X", "Y"):
        LXt, LYt = L_for(PX[top], m_["rX"]), L_for(PY[top], m_["rY"])
        preds[f"{top}top_(LX,LY)"] = NB * pxy(top, LXt, LYt)
        preds[f"{top}top_(LX,LX)"] = NB * pxy(top, LXt, LXt)
        preds[f"{top}top_(LY,LY)"] = NB * pxy(top, LYt, LYt)
        d[f"L_X_mm_{top}top"], d[f"L_Y_mm_{top}top"] = LXt, LYt
    central = 0.5 * (preds["Xtop_(LX,LY)"] + preds["Ytop_(LX,LY)"])
    d["predicted_N_XY"] = central
    d["predicted_range"] = [min(preds.values()), max(preds.values())]
    d["predictions"] = preds
    obs = m_["N_XY"]
    lo_mu = d["predicted_range"][0]
    d["poisson_P_le_obs_at_central"] = float(sum(math.exp(-central) * central ** k / math.factorial(k) for k in range(obs + 1))) if obs < 50 else None
    d["poisson_P_le_obs_at_range_min"] = float(sum(math.exp(-lo_mu) * lo_mu ** k / math.factorial(k) for k in range(obs + 1))) if obs < 50 else None
    S["comparison"][thr] = d
    print(thr, "L_X %.2f mm, L_Y %.2f mm; X·Y predichas %.1f (rango %.1f-%.1f); observadas %s" % (LX, LY, central, *d["predicted_range"], obs),
          "P(≤obs) %s / %s" % (d["poisson_P_le_obs_at_central"], d["poisson_P_le_obs_at_range_min"]))
fig = {"L": Lup[::4].tolist(), "PX": (0.5 * (PX["X"] + PX["Y"]))[::4].tolist(), "PY": (0.5 * (PY["X"] + PY["Y"]))[::4].tolist(),
       "PXY": [0.5 * (pxy("X", L, L) + pxy("Y", L, L)) for L in Lup[::4]], "measured": {k: {"rX": v["rX"], "rY": v["rY"], "pXY": v["N_XY"] / NB} for k, v in MEAS.items()}}
json.dump(fig, open(W + "/fig_sensitive_length.json", "w"), indent=1)
json.dump(S, open(W + "/sensitive_length.json", "w"), indent=1, default=float)
print("hecho")
