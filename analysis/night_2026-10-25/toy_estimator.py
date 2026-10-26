#!/usr/bin/python3 -I
# toy_estimator.py — XPLUS-4 etapa 4: toy MC del estimador X·Y. No reimplementa el estimador (D2): genera los cuatro arreglos
# booleanos por evento (Xs, Ys en tiempo; Xc, Yc fuera de tiempo) y llama a xy_estimator.xy_fraction, el mismo código de PROC.
# Modelo por evento (todas las elecciones están en el informe XPLUS-4, etapa 4):
#  - cruce del muón con prob. c = f_2018·g; si cruza, impacto verdadero en X con ε_X y en Y con ε_Y (independientes).
#    ε_X = ε_Y = sqrt(q_mu/(f·g)), q_mu = q barrido (si q_mu > f·g: ε = 1 y c = q_mu). Caso asimétrico: ε = 0.95/0.75, c = f·g.
#  - actividad verdadera no muónica en tiempo por capa, B_X, B_Y (prob. b_X, b_Y), elegida para que P(T_X) = r_X medido,
#    T_X = muón_X ∨ B_X; P(B_X ∧ B_Y) = k·b_X·b_Y (k = 1 base; k = 1.5 variante "correlación mayor que la supuesta").
#  - ruido independiente en tiempo, prob. a_X, a_Y; fuera de tiempo, ruido independiente con prob. a_X^oot = (1+δ)·a_X.
#  - verdad del estimador: q_T = P(T_X ∧ T_Y) (exacta, por enumeración); verdad muónica: q_mu = c·ε_X·ε_Y.
# Salida: JSON con, por configuración y estimador (mix, ctl): media de q̂ − q_T y de q̂ − q_mu con su error, sd(q̂)/<σ̂>, cobertura
# a 1σ (|q̂ − q_T| < σ̂), diferencia media mix − ctl. Uso: toy_estimator.py <salida.json> <grupo> <réplicas> <semilla> [XPLUS4]
import json, math, os, sys, time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import xy_estimator as xy  # noqa: E402

OUT, GROUP, NREP, SEED = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
F18, G = 0.03522, 0.980  # results.json v3: f_2018 (R1) y g (R2)
# a, r medidos (fixtures XPLUS-3 = work/stage2.json, T = 5 p.e.)
DATA = {"v8": {"N": 113152, "aX": 0.26459983031674206, "aY": 0.2711927319004525, "rX": 0.06056818728068068, "rY": 0.07474595590910194},
        "run8": {"N": 44176, "aX": 0.14430912712785224, "aY": 0.12552064469395147, "rX": 0.060421681966085536, "rY": 0.06608682146462685}}


def configs():
    base = dict(asym=False, a_scale=1.0, k=1.0, delta=0.0, zero=False)
    out = []
    for ds in DATA:
        for q in (0.015, 0.0275, 0.04):
            out += [dict(base, ds=ds, q=q, name="base"),
                    dict(base, ds=ds, q=q, name="a_doble", a_scale=2.0),
                    dict(base, ds=ds, q=q, name="corr_1.5", k=1.5),
                    dict(base, ds=ds, q=q, name="oot_+20%", delta=0.2),
                    dict(base, ds=ds, q=q, name="oot_-20%", delta=-0.2),
                    dict(base, ds=ds, q=q, name="validacion_a0_r0", zero=True)]
        out.append(dict(base, ds=ds, q=None, name="asimetrico_0.95_0.75", asym=True))
    return out


def delta_scan():
    return [dict(ds=ds, q=0.0275, name=f"delta_{d:+.2f}", asym=False, a_scale=1.0, k=1.0, delta=d, zero=False)
            for ds in DATA for d in (-0.20, -0.10, -0.05, 0.0, 0.05, 0.10, 0.20, 0.30)]


def params(cf):
    d = DATA[cf["ds"]]
    c = F18 * G
    if cf["asym"]:
        eX, eY = 0.95, 0.75
    else:
        e = math.sqrt(cf["q"] / c)
        if e > 1:
            e, c = 1.0, cf["q"]
        eX = eY = e
    aX, aY = d["aX"] * cf["a_scale"], d["aY"] * cf["a_scale"]
    if cf["zero"]:
        aX = aY = 0.0
        bX = bY = 0.0
    else:
        bX = 1 - (1 - d["rX"]) / (1 - c * eX)
        bY = 1 - (1 - d["rY"]) / (1 - c * eY)
    k = cf["k"]
    p11 = k * bX * bY
    pb = {(1, 1): p11, (1, 0): bX - p11, (0, 1): bY - p11, (0, 0): 1 - bX - bY + p11}
    # q_T = P(T_X ∧ T_Y) exacta
    qT = 0.0
    for cross, pc in ((1, c), (0, 1 - c)):
        for (bx, by), pbb in pb.items():
            px = 1.0 if bx else (eX if cross else 0.0)
            py = 1.0 if by else (eY if cross else 0.0)
            qT += pc * pbb * px * py
    return dict(N=DATA[cf["ds"]]["N"], c=c, eX=eX, eY=eY, aX=aX, aY=aY, aXo=aX * (1 + cf["delta"]), aYo=aY * (1 + cf["delta"]),
                bX=bX, bY=bY, pb=pb, q_mu=c * eX * eY, q_T=qT)


def replica(rng, p):
    n = p["N"]
    cross = rng.random(n) < p["c"]
    mx = cross & (rng.random(n) < p["eX"])
    my = cross & (rng.random(n) < p["eY"])
    u = rng.random(n)
    p11, p10, p01 = p["pb"][(1, 1)], p["pb"][(1, 0)], p["pb"][(0, 1)]
    bx = u < p11 + p10
    by = (u < p11) | ((u >= p11 + p10) & (u < p11 + p10 + p01))
    Xs = mx | bx | (rng.random(n) < p["aX"])
    Ys = my | by | (rng.random(n) < p["aY"])
    Xc = rng.random(n) < p["aXo"]
    Yc = rng.random(n) < p["aYo"]
    return xy.xy_fraction(Xs, Ys, Xc, Yc)


cfs = configs() if GROUP == "main" else delta_scan()
rng = np.random.default_rng(SEED)
res = []
t0 = time.time()
for i, cf in enumerate(cfs):
    p = params(cf)
    qs = {"mix": [], "ctl": []}
    es = {"mix": [], "ctl": []}
    for _ in range(NREP):
        r = replica(rng, p)
        for t in ("mix", "ctl"):
            qs[t].append(r[f"f_{t}"])
            es[t].append(r[f"f_{t}_err"])
    row = {"config": cf, "params": {k: v for k, v in p.items() if k != "pb"}, "replicas": NREP}
    for t in ("mix", "ctl"):
        q, e = np.array(qs[t]), np.array(es[t])
        row[t] = {"mean": float(q.mean()), "bias_vs_qT": float(q.mean() - p["q_T"]), "bias_vs_qmu": float(q.mean() - p["q_mu"]),
                  "bias_err": float(q.std(ddof=1) / math.sqrt(NREP)), "sd": float(q.std(ddof=1)), "mean_sigma": float(e.mean()),
                  "sd_over_sigma": float(q.std(ddof=1) / e.mean()), "coverage_1sigma": float(np.mean(np.abs(q - p["q_T"]) < e))}
    d = np.array(qs["mix"]) - np.array(qs["ctl"])
    row["mix_minus_ctl"] = {"mean": float(d.mean()), "err": float(d.std(ddof=1) / math.sqrt(NREP)), "sd": float(d.std(ddof=1))}
    res.append(row)
    print(f"{i + 1}/{len(cfs)} {cf['ds']} {cf['name']} q={cf['q']}: q_T {p['q_T']:.5f} mix {row['mix']['bias_vs_qT']:+.2e} "
          f"ctl {row['ctl']['bias_vs_qT']:+.2e} mix-ctl {row['mix_minus_ctl']['mean']:+.2e} ({time.time() - t0:.0f} s)", flush=True)
json.dump({"group": GROUP, "replicas": NREP, "seed": SEED, "f_2018": F18, "g": G, "data": DATA, "results": res,
           "elapsed_s": time.time() - t0}, open(OUT, "w"), indent=1)
print("fin", OUT, flush=True)
