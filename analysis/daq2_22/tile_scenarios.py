#!/usr/bin/python3 -I
# tile_scenarios.py — DAQ2-22 etapa 1: escenarios geométricos del trigger nuevo (dos tejas START) con el MC del banco
# (daq/tools/characterization/bench_mc.py: cos³, hod2018 del CAD, mismas trazas para todos los escenarios).
# Por escenario (T1 abajo, T2 arriba; tejas como rectángulos en sus planos medios):
#   f_2018^MC = P(cruza el cuadrado activo de 33 × 33 mm en z = 0 | trigger) (la S de bar_mc.py),
#   (rX/rX₈)^MC = P(alguna barra X | trigger) / lo mismo con las placas de run 8 (CAD: 150 mm, desplazadas ±42.5 mm),
#   R_MC/R₈    = triggers del escenario / triggers con las placas de run 8 (las dos con eficiencia 100 %).
# Con lo medido en DAQ2-21 (rX/rX₈ = 0.602 ± 0.055, R₁₀/R₈ = 0.387 ± 0.004):
#   a = 1 − 0.602 / (rX/rX₈)^MC   (fracción de triggers sin muón),
#   ε_par = 0.387 (1 − a) / (R_MC/R₈)   (eficiencia del par de tejas relativa a la del par de run 8, supuesto: run 8 = 100 %).
# Uso: tile_scenarios.py <salida dir> [N trazas]
import json, math, os, sys
import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import bench_mc as bm  # noqa: E402

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 60_000_000
NT = 6_000_000  # submuestra para contar triggers (todas las trazas)
CH = 2_000_000
RX_MEAS, RX_ERR = 0.602, 0.055
R_MEAS, R_ERR = 0.387, 0.004
ZS = {"CAD": (125.6, 135.6), "CAD+20": (145.6, 155.6)}  # planos medios de T1, T2
ZALL = sorted({z for v in ZS.values() for z in v})
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(2222)
T = {z: [] for z in ZALL}
Hp = {z: [] for z in ZALL}
HX, HS = [], []
nT = 0
for s in range(0, N, CH):
    m = min(CH, N - s)
    r = bm.generate(rng, m)
    pos = {z: r.prop(z) for z in ZALL}
    if nT < NT:
        k = min(m, NT - nT)
        for z in ZALL:
            T[z].append(np.stack([pos[z][0][:k], pos[z][1][:k]]).astype(np.float32))
        nT += k
    x0, y0 = r.prop(0.0)
    sel = (np.abs(x0) <= 45) & (np.abs(y0) <= 45)
    rs = bm.Rays(r.x0[sel], r.y0[sel], r.dx[sel], r.dy[sel], r.dz[sel])
    bp = bm.bar_paths(rs, 0.0, "X")
    ax = (bp["X"] > 0).any(1)
    sq = (np.abs(x0[sel]) <= bm.HALF) & (np.abs(y0[sel]) <= bm.HALF)
    keep = ax | sq
    for z in ZALL:
        Hp[z].append(np.stack([pos[z][0][sel][keep], pos[z][1][sel][keep]]).astype(np.float32))
    HX.append(ax[keep]); HS.append(sq[keep])
T = {z: np.concatenate(v, axis=1) for z, v in T.items()}
Hp = {z: np.concatenate(v, axis=1) for z, v in Hp.items()}
HX, HS = np.concatenate(HX), np.concatenate(HS)
SC = N / NT


def inside(p, c, h):
    return (np.abs(p[0] - c[0]) <= h[0]) & (np.abs(p[1] - c[1]) <= h[1])


def counts(t1, t2, z):
    """t1, t2 = (centro, semilados); z = clave de ZS. Triggers (escalados a N), cruces S, cruces X."""
    z1, z2 = ZS[z]
    nt = (inside(T[z1], *t1) & inside(T[z2], *t2)).sum() * SC
    th = inside(Hp[z1], *t1) & inside(Hp[z2], *t2)
    return nt, int((th & HS).sum()), int((th & HX).sum())


old = counts(((-42.5, -42.5), (75, 75)), ((42.5, 42.5), (75, 75)), "CAD")
OLD = {"T": old[0], "f": old[1] / old[0], "pX": old[2] / old[0]}
scen = []


def add(name, t1, t2, z):
    nt, ns, nx = counts(t1, t2, z)
    if nt <= 0:
        return
    f, px = ns / nt, nx / nt
    rr, rate = px / OLD["pX"], nt / OLD["T"]
    a = 1 - RX_MEAS / rr
    ea = RX_ERR / rr
    eps = R_MEAS * (1 - a) / rate
    eeps = eps * math.hypot(R_ERR / R_MEAS, ea / (1 - a)) if a < 1 else float("nan")
    flag = "imposible (a < 0)" if a < 0 else ("ε_par > 1" if eps > 1 else "")
    if a < 0 and eps > 1:
        flag = "imposible (a < 0) y ε_par > 1"
    scen.append({"scenario": name, "heights": z, "f2018_MC": f, "f2018_MC_relerr": 1 / math.sqrt(max(ns, 1)), "rX_ratio_MC": rr,
                 "R_MC_over_R8": rate, "a": a, "a_err": ea, "eps_pair": eps, "eps_pair_err": eeps, "flag": flag, "n_trig_MC": nt, "n_S_MC": ns})


for z in ZS:
    for L in (100.0, 150.0):
        h = (L / 2, L / 2)
        for d in (0, 20, 40, 60):
            add(f"{int(L)} mm, alineadas, corridas {d} mm en x", ((d, 0), h), ((d, 0), h), z)
        for d in (0, 20, 40, 60):
            add(f"{int(L)} mm, T1 corrida {d} mm, T2 {d + 20} mm en x", ((d, 0), h), ((d + 20, 0), h), z)
    # tiras 5 × 200 × 10 mm: cruzadas (T1 a lo largo de x, T2 a lo largo de y) y paralelas, centradas
    add("tiras 5 × 200, cruzadas, centradas", ((0, 0), (100, 2.5)), ((0, 0), (2.5, 100)), z)
    add("tiras 5 × 200, paralelas, centradas", ((0, 0), (100, 2.5)), ((0, 0), (100, 2.5)), z)
    # fuera de la lista pedida: desplazamiento mutuo grande (como las placas de run 8, en diagonal), punto medio sobre hod2018
    for L in (100.0, 150.0):
        h = (L / 2, L / 2)
        for D in (40, 60, 80, 100, 120):
            q = D / (2 * math.sqrt(2))
            add(f"extra: {int(L)} mm, desplazadas {D} mm entre sí en diagonal", ((-q, -q), h), ((q, q), h), z)
res = {"N": N, "NT": NT, "model": "cos3, hod2018 CAD (bench_mc), tejas como rectángulos en z medio", "z_mid_mm": ZS,
       "run8_MC": {"f2018": OLD["f"], "pX": OLD["pX"], "n_trig": OLD["T"]}, "measured": {"rX_ratio": [RX_MEAS, RX_ERR], "R_ratio": [R_MEAS, R_ERR]},
       "scenarios": scen}
json.dump(res, open(f"{OUT}/tile_scenarios.json", "w"), indent=1)
print(f"run 8 (MC): f = {OLD['f']:.5f}, pX = {OLD['pX']:.5f}")
for s in scen:
    print(f"{s['heights']:6s} {s['scenario']:44s} f {s['f2018_MC']:.4f}  rX/rX8 {s['rX_ratio_MC']:.3f}  R/R8 {s['R_MC_over_R8']:.3f}  "
          f"a {s['a']:+.2f}±{s['a_err']:.2f}  ε_par {s['eps_pair']:.2f}±{s['eps_pair_err']:.2f}  {s['flag']}")
