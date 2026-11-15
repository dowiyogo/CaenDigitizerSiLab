#!/usr/bin/python3 -I
# h19_fold.py — seguimiento 2026-11-15 (15, punto 1): doblado de la probabilidad de detección de hod2019 contra la distancia al extremo de lectura
# (Geant4, modelo pintado de MSI, inputs/msi/2026-11-15_optics_followup_14/detection_probability.csv: bins de 1 mm, umbrales 1, 2, 3, 5, 10 p.e.)
# con el MC del banco de stage5 (bench_mc.py: cos³θ, placas del CAD, rayos rectos, semilla 42, trozos de 2e6, N = 1e8 como h19_end_mc.py).
# Geometría y extremos de lectura como h19_end_mc.py: hod2019 en z = 98.4 mm, barras 75 × 3 × 1 mm, MPPC en un extremo (±37.5 mm), las dos capas
# de una coordenada leídas por extremos opuestos (registro de conexiones, kCfg_2026_05 / handoff §5: X tira par = Delantera, impar = Trasera;
# Y par = Opuesta, impar = Frontal). Qué cara está en + no está registrado: las cuatro asignaciones (A+/A− por coordenada) y los dos órdenes de
# bloque (X arriba, Y arriba); el resultado central es el promedio de los dos órdenes de bloque, como sensitive_length.py y h19_end_test.py.
# Por trayectoria del trigger y barra cruzada (camino > 0): d = distancia al extremo con MPPC de la coordenada a lo largo de la barra en el plano
# medio de su capa; p_T(d) = probabilidad de detección de Geant4 (interpolación lineal entre centros de bin, extremos fijos). Barras independientes:
#   P_X,T(traza) = 1 − Π_{barras X cruzadas} (1 − p_T(d)),  igual P_Y;  P_XY(traza) = P_X,T_X · P_Y,T_Y (mismo evento, barras independientes).
# Por trigger: medias sobre las trazas. Dos estaciones: P_XY con la traza dentro del cuadrado activo de hod2018 (bench_mc.in_square), eficiencia de
# hod2018 = 1 (cota geométrica, como h19_end_test.py; el manuscrito exige > 50 mV en X e Y de hod2018).
# Umbral que reproduce r_X = 0.0144 y r_Y = 0.0093 (35 mV, v8 completo, h19_end_test.py): interpolación lineal de ln P contra T entre los umbrales de
# la malla; ganancia implícita = 35 mV / T. P_XY en (T_X, T_Y): interpolación bilineal de ln P_XY en la malla 5 × 5.
# Error de Geant4: se repite el doblado con p ± error estándar en todos los bins a la vez (correlación total, cota conservadora).
# Fuera de la malla (r_X, r_Y piden más de 10 p.e.): (a) prolongación lineal de ln P con el tramo 5-10 p.e.; (b) modelo de Poisson con el número
# medio de p.e. de mean_npe_vs_distance.csv (mismo paquete), umbrales hasta 90 p.e.; los dos son modelos, no salidas de Geant4.
# Salida: <out>/h19_fold.json, h19_fold_P.csv. Uso: h19_fold.py <detection_probability.csv> <out dir> [N]
import csv, hashlib, json, math, os, sys, time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import bench_mc as bm  # noqa: E402

DP, OUT = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000_000
os.makedirs(OUT, exist_ok=True)
G = bm.GEOMS["CAD"]
HL = bm.BAR_L / 2
CENTRES = sorted([(c, "A") for c in bm.LAYER_A] + [(c, "B") for c in bm.LAYER_B])
NRUN = 143814  # run B = mayo v8 completo (h19_xy.py, h19_end_test.py)
MEAS = {"rX": 0.0144, "rY": 0.0093, "N_XY_same": 3, "N_XY_mixed": 26, "N_two_station": 1, "thr_mV": 35.0}

# ---------------------------------------------------------------- probabilidad de detección
rows = list(csv.DictReader(open(DP)))
TS = sorted({int(r["threshold_pe"]) for r in rows})
DC = sorted({float(r["distance_mm_from_readout_bin_center"]) for r in rows})
PT = {v: np.zeros(len(DC)) for v in ("c", "lo", "hi")}
PT = {T: {v: np.zeros(len(DC)) for v in ("c", "lo", "hi")} for T in TS}
for r in rows:
    T, i = int(r["threshold_pe"]), DC.index(float(r["distance_mm_from_readout_bin_center"]))
    p, e = float(r["probability"]), float(r["standard_error"])
    PT[T]["c"][i], PT[T]["lo"][i], PT[T]["hi"][i] = p, max(p - e, 0.0), min(p + e, 1.0)
DC = np.asarray(DC)
# modelo fuera de la malla (umbrales > 10 p.e.): N ~ Poisson(μ(d)), μ(d) = mean_npe_true de mean_npe_vs_distance.csv (mismo paquete de MSI,
# interpolación lineal, extremos fijos). Sin la fluctuación de Landau del depósito: se valida contra la malla de Geant4 en 1-10 p.e.
MN = os.path.join(os.path.dirname(DP), "mean_npe_vs_distance.csv")
mrows = list(csv.DictReader(open(MN)))
MD = np.array([float(r["distance_mm_from_readout_bin_center"]) for r in mrows])
MU = np.array([float(r["mean_npe_true"]) for r in mrows])
MUE = np.array([float(r["standard_error"]) for r in mrows])
TM = [1, 2, 3, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90]
LGF = np.array([math.lgamma(k + 1) for k in range(max(TM) + 1)])


def pois_ge(mu, T):
    """P(N ≥ T | Poisson(mu)), vectorizado."""
    mu = np.maximum(mu, 1e-12)
    k = np.arange(T)
    lp = k[:, None] * np.log(mu)[None, :] - mu[None, :] - LGF[:T][:, None]
    return np.clip(1.0 - np.exp(lp).sum(0), 0.0, 1.0)


# variantes: Geant4 (c, lo, hi = p ± error estándar en todos los bins) en TS; Poisson (pc, plo, phi = μ ± error estándar) en TM
VARS = ("c", "lo", "hi", "pc", "plo", "phi")


def thr_list(v):
    return TS if v in ("c", "lo", "hi") else TM


def pdet(d, k, v):
    if v in ("c", "lo", "hi"):
        return np.interp(d, DC, PT[TS[k]][v])  # fuera de [0.5, 74.5]: valor del bin extremo
    mu = np.interp(d, MD, MU + {"pc": 0.0, "plo": -1.0, "phi": 1.0}[v] * MUE)
    return pois_ge(mu, TM[k])


CONV = [(cx, cy) for cx in ("A+", "A-") for cy in ("A+", "A-")]
# acumuladores: S1[top][conv][v]["X"][iT] = Σ P_X; S2[...][iTx, iTy] = Σ P_X P_Y; S2h = idem dentro de hod2018
acc = {(top, cx, cy, v): {"X": np.zeros(len(thr_list(v))), "Y": np.zeros(len(thr_list(v))), "XY": np.zeros((len(thr_list(v)),) * 2),
                          "XYh": np.zeros((len(thr_list(v)),) * 2)} for top in "XY" for cx, cy in CONV for v in VARS}
trig = 0
rng = np.random.default_rng(42)
t0 = time.time()


def miss_prob(rt, zl, coord, conv, path, v):
    """Π (1 − p_T(d)) sobre las barras cruzadas de la coordenada, para cada T: arreglo (nT, ntraza)."""
    q = np.ones((len(thr_list(v)), len(rt.dz)))
    for j, (c, layer) in enumerate(CENTRES):
        x, y = rt.prop(zl[(coord, layer)])
        w = y if coord == "X" else x
        plus = (layer == "A") == (conv == "A+")
        d = HL - w if plus else w + HL
        hit = path[:, j] > 0
        if not hit.any():
            continue
        for k in range(len(thr_list(v))):
            q[k, hit] *= 1.0 - pdet(d[hit], k, v)
    return q


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
        for v in VARS:
            PX = {cv: 1.0 - miss_prob(rt, zl, "X", cv, pth["X"], v) for cv in ("A+", "A-")}
            PY = {cv: 1.0 - miss_prob(rt, zl, "Y", cv, pth["Y"], v) for cv in ("A+", "A-")}
            for cx, cy in CONV:
                a = acc[(top, cx, cy, v)]
                a["X"] += PX[cx].sum(1)
                a["Y"] += PY[cy].sum(1)
                a["XY"] += PX[cx] @ PY[cy].T
                a["XYh"] += PX[cx][:, s18] @ PY[cy][:, s18].T
secs = time.time() - t0


def interp_T(P, r, TT, extrap=False):
    """T (p.e.) con P(T) = r, interpolación lineal de ln P en T; con extrap, prolonga el último tramo; (T, extrapolado?) o (None, None)."""
    lp = np.log(np.maximum(P, 1e-300))
    lr = math.log(r)
    for k in range(len(TT) - 1):
        a, b = lp[k], lp[k + 1]
        if (a - lr) * (b - lr) <= 0 and a != b:
            return TT[k] + (lr - a) / (b - a) * (TT[k + 1] - TT[k]), False
    if extrap and lr < lp[-1] and lp[-1] < lp[-2]:
        a, b = lp[-2], lp[-1]
        return TT[-2] + (lr - a) / (b - a) * (TT[-1] - TT[-2]), True
    return None, None


def interp_2d(M, tx, ty, TT):
    """ln M bilineal en (T_X, T_Y) (prolonga el último tramo fuera de la malla)."""
    L = np.log(np.maximum(M, 1e-300))
    nT = len(TT)

    def idx(tv):
        k = max(0, min(nT - 2, int(np.searchsorted(TT, tv)) - 1))
        return k, (tv - TT[k]) / (TT[k + 1] - TT[k])
    i, fi = idx(tx)
    j, fj = idx(ty)
    v = (1 - fi) * (1 - fj) * L[i, j] + fi * (1 - fj) * L[i + 1, j] + (1 - fi) * fj * L[i, j + 1] + fi * fj * L[i + 1, j + 1]
    return float(math.exp(v))


res = {"input": {"file": os.path.relpath(DP, REPO), "sha256": hashlib.sha256(open(DP, "rb").read()).hexdigest(), "thresholds_pe": TS,
                 "mean_npe_file": os.path.relpath(MN, REPO), "mean_npe_sha256": hashlib.sha256(open(MN, "rb").read()).hexdigest(), "model_thresholds_pe": TM,
                 "distance_bins_mm": [float(DC[0]), float(DC[-1]), len(DC)]},
       "mc": {"rays": N, "triggers": trig, "seed": 42, "model": "cos3", "geometry": "CAD", "z_h19_mm": G["H19"], "seconds": secs},
       "run_B_triggers": NRUN, "measured_35mV": MEAS, "per_conv": {}}
prow = []
for cx, cy in CONV:
    out = {}
    for v in VARS:
        tops = {}
        for top in "XY":
            a = acc[(top, cx, cy, v)]
            tops[top] = {"PX": a["X"] / trig, "PY": a["Y"] / trig, "PXY": a["XY"] / trig, "PXYh": a["XYh"] / trig}
        avg = {k: 0.5 * (tops["X"][k] + tops["Y"][k]) for k in tops["X"]}
        o = {"PX": avg["PX"].tolist(), "PY": avg["PY"].tolist(), "PXY_diag": np.diag(avg["PXY"]).tolist(),
             "PXY_h18_diag": np.diag(avg["PXYh"]).tolist(), "PXY_grid": avg["PXY"].tolist(), "PXY_h18_grid": avg["PXYh"].tolist(),
             "Xtop": {"PX": tops["X"]["PX"].tolist(), "PY": tops["X"]["PY"].tolist()}, "Ytop": {"PX": tops["Y"]["PX"].tolist(), "PY": tops["Y"]["PY"].tolist()}}
        TT = thr_list(v)
        o["thresholds_pe"] = TT
        (TX, ex), (TY, ey) = interp_T(avg["PX"], MEAS["rX"], TT, True), interp_T(avg["PY"], MEAS["rY"], TT, True)
        o["T_X_pe"], o["T_Y_pe"], o["T_X_extrapolated"], o["T_Y_extrapolated"] = TX, TY, ex, ey
        o["gain_X_mV_per_pe"] = MEAS["thr_mV"] / TX if TX else None
        o["gain_Y_mV_per_pe"] = MEAS["thr_mV"] / TY if TY else None
        if TX and TY:
            pxy, pxyh = interp_2d(avg["PXY"], TX, TY, TT), interp_2d(avg["PXYh"], TX, TY, TT)
            o.update(PXY_at_T=pxy, N_XY_pred=NRUN * pxy, N_XY_indep_pred=NRUN * MEAS["rX"] * MEAS["rY"], PXY_h18_at_T=pxyh,
                     N_two_station_pred_geom=NRUN * pxyh, P_Y_given_X=pxy / MEAS["rX"])
        out[v] = o
        if v in ("c", "pc"):
            for k, T in enumerate(TT):
                prow.append(["geant4" if v == "c" else "poisson_mean_npe", cx, cy, T, f"{avg['PX'][k]:.6g}", f"{avg['PY'][k]:.6g}", f"{avg['PXY'][k, k]:.6g}", f"{NRUN * avg['PXY'][k, k]:.4g}",
                             f"{avg['PXYh'][k, k]:.6g}", f"{NRUN * avg['PXYh'][k, k]:.4g}",
                             f"{tops['X']['PX'][k]:.6g}", f"{tops['Y']['PX'][k]:.6g}", f"{tops['X']['PY'][k]:.6g}", f"{tops['Y']['PY'][k]:.6g}"])
    res["per_conv"][f"{cx}/{cy}"] = out
with open(f"{OUT}/h19_fold_P.csv", "w") as f:
    f.write(f"# h19_fold.py: bench MC {N} rays ({trig} triggers), hod2019 at z = {G['H19']} mm, Geant4 detection probability vs distance from the readout "
            f"end ({os.path.relpath(DP, REPO)}, sha256 {res['input']['sha256']}); average of the two block orders; conv = readout end of layer A at "
            "+37.5 (A+) or -37.5 (A-); probabilities per trigger; N_* = run B (143814 triggers) x P; two-station = hod2018 efficiency 1\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["model", "conv_X", "conv_Y", "threshold_pe", "PX", "PY", "PXY_same_event", "N_XY_runB", "PXY_h18_square", "N_two_station_runB_geom",
                "PX_Xtop", "PX_Ytop", "PY_Xtop", "PY_Ytop"])
    w.writerows(prow)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


json.dump(clean(res), open(f"{OUT}/h19_fold.json", "w"), indent=1)
print(f"entrada {res['input']['file']} sha256 {res['input']['sha256']}")
print(f"MC: rayos {N}, triggers {trig}, {secs:.0f} s; run B N = {NRUN}; independiente N·rX·rY = {NRUN * MEAS['rX'] * MEAS['rY']:.2f}")
for c, o in res["per_conv"].items():
    for v in VARS:
        x = o[v]
        print(f"{c} [{v}]: " + " ".join(f"T{T}: PX {x['PX'][k]:.5f} PY {x['PY'][k]:.5f} PXY {x['PXY_diag'][k]:.3e} NXY {NRUN * x['PXY_diag'][k]:.3f} "
                                        f"N2st {NRUN * x['PXY_h18_diag'][k]:.4f} |" for k, T in enumerate(x["thresholds_pe"])))
        if x["T_X_pe"] and x["T_Y_pe"]:
            print(f"   -> (extrapolado X {x['T_X_extrapolated']}, Y {x['T_Y_extrapolated']}) T_X {x['T_X_pe']:.3f} p.e. (G {x['gain_X_mV_per_pe']:.2f} mV/p.e.), T_Y {x['T_Y_pe']:.3f} p.e. (G {x['gain_Y_mV_per_pe']:.2f} mV/p.e.); "
                  f"P_XY {x['PXY_at_T']:.3e}, N_XY pred {x['N_XY_pred']:.3f} (obs 3, mezclados 26); dos estaciones pred (geom.) {x['N_two_station_pred_geom']:.4f} (obs 1)")
        else:
            print(f"   -> T_X {x['T_X_pe']} T_Y {x['T_Y_pe']} (fuera de la malla)")
