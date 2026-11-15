#!/usr/bin/python3 -I
# cfd_window_eff2.py — seguimiento 2026-11-15 (ventana CFD, versión 2): ε_XY de hod2018 con ventanas de coincidencia sobre tiempos CFD, para mayo v8
# (tiempo bueno, 77 558 disparos) y octubre run 8 (44 176 disparos, todos, como stage2.py), con el mismo estimador para las dos campañas.
# Entrada: <campaña>_time2.bin de time_extract2.C (pulso más alto por canal en la ventana signal [49,58) y en la fuera de tiempo [99,108), 800 ns
# después: amplitud, k, p_loc, x[k−7 … k+2]); desplazamientos cfd50 por canal (mayo: seguimiento 13 A, A_time/time_res_offsets.csv; octubre:
# time_res2.py oct, A_time/time_res_offsets.csv del paquete de octubre), columna offset_ns, t' = t_cfd − o[c] (y − 800 ns fuera de tiempo).
# G: mayo = reajuste de PROC (hd64_may_v8_peaks.csv), octubre = c1 (mediana para los canales sin G), como stage2.py; ch3 fuera.
# ε_XY = q_mix/(f_2018·g) con f_2018 de stage2.json y g = 0.9551211 (g_dead_ch3.json) en las dos campañas; aceptancia relativa 0.046481.
# Definiciones de coincidencia (una tira cuenta si amplitud > T·G y tiempo CFD finito):
#  - 144ns: ventana signal, xy_estimator.xy_fraction (= stage2.json, se comprueba); 144ns_general: el estimador generalizado con la misma ventana.
#  - medianW: |t' − t_ref| < W, t_ref = mediana por evento de los t' de todas las tiras de hod2018 sobre el umbral (W = 8, 16 ns).
#  - pairW: alguna tira X y alguna tira Y con |t'_X − t'_Y| < W, sin referencia (W = 8, 16 ns).
# Estimador generalizado (cabecera de cfd_window_eff.py): S = N_XY/N, M = acc_mix/N con mezclados (tiras X del evento i e Y del i+1, y al revés)
# con la misma definición; q_mix = [S − M + rX·rY·(A_mm − A_X − A_Y + A_0)]/D, q_ctl = [S − rX·A_X − rY·A_Y − (1 − rX − rY)·A_0]/D,
# D = A_XY − A_X − A_Y + A_0; rX, rY de la ventana de 144 ns.
# A_* por incrustación, ruido = pulsos fuera de tiempo de cada evento sobre T·G:
#  - "measured" (versión 2): con su tiempo CFD medido, t'_C − 800 ns (correlaciones entre canales incluidas); los muones verdaderos (un canal X y uno Y
#    al azar) con el tiempo medio de un muón etiquetado al azar (X, Y > 10 p.e., |t'_X − t'_Y| < 5 ns); A_mm con t'_X, t'_Y de dos muones distintos.
#    Un pulso fuera de tiempo sin cruce del 50 % no cuenta (como en la ventana signal).
#  - "uniform" (versión 1, cfd_window_eff.py): tiempos uniformes en los 144 ns, para comparar.
# Control: en muones etiquetados, densidad de tiras sobre T (sin las etiquetadas ni sus vecinas en el plano) con t' − t_ref ∈ [24, 72) ns en datos y
# en el ruido incrustado medido.
# Errores: stat = √(N_XY + (m1 + m2)/4)/(N·D)/(f_2018·g); estimador = |q_mix − q_ctl|/(f_2018·g); aceptancia = ε·0.046481.
# Salida: <out>/cfd_window_eff2_<campaña>.csv, .json. Uso: cfd_window_eff2.py may|oct <time2.bin> <offsets.csv> <out>
import csv, json, math, os, sys, warnings

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

CAMP, TB, OFFS, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
GD = json.load(open(os.path.join(REPO, "outputs/proceeding_revision/2026-11-02_revision/g_dead_ch3/g_dead_ch3.json")))
S2 = json.load(open(PK + "/stage2.json"))
F18, G_XY, ACC_REL = S2["f_2018"], GD["g_used_new"]["value"], GD["acceptance_rel_syst_no_ch3"]
S2KEY = {"may": "may8_good", "oct": "oct8"}[CAMP]
LABEL = {"may": "may_v8_goodtime", "oct": "oct_run8"}[CAMP]
WIN = [("a", "<f4", 32), ("k", "<u2", 32), ("pl", "<f4", 32), ("w", "<f4", (32, 10))]
DT2 = np.dtype([("t", "<f8"), ("bl", "<f4", 32)] + [(n + "S",) + tuple(r) for n, *r in WIN] + [(n + "C",) + tuple(r) for n, *r in WIN])
XCH = [c for c in range(16) if c != 3]
YCH = list(range(16, 32))
TS = [3, 5, 7, 10]
WS = [8.0, 16.0]
SPAN_NS, KPRE, DS, S0, DC = 144.0, 7, 16.0, 49, 800.0
TAG = 10.0
SEED = 2026111502
fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))  # noqa: E731


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return fill(g)


def cfd50(pl, w, k):
    """cfd50 de time_res.py (ns desde el inicio del tramo), (N, 32); NaN sin cruce o sin pulso."""
    a = pl.astype(float)[:, :, None] - w.astype(float)
    A = a[:, :, KPRE]
    below = a[:, :, :KPRE] < 0.5 * A[:, :, None]
    last = np.where(below.any(2), KPRE - 1 - np.argmax(below[:, :, ::-1], 2), -1)
    li = np.clip(last, 0, KPRE - 1)
    a0 = np.take_along_axis(a, li[:, :, None], 2)[:, :, 0]
    a1 = np.take_along_axis(a, (li + 1)[:, :, None], 2)[:, :, 0]
    frac = (0.5 * A - a0) / np.where(a1 != a0, a1 - a0, 1)
    t = (k.astype(float) - KPRE + li + frac) * DS
    return np.where((last >= 0) & (k != 0xFFFF), t, np.nan)


G = g_refit() if CAMP == "may" else fill(mc.g_c1())
e = np.fromfile(TB, dtype=DT2)
e = e[np.argsort(e["t"], kind="stable")]
if CAMP == "may":
    e = e[mc.good_time(e["t"], "v8")]
N = len(e)
off = np.zeros(32)
for r in csv.DictReader(open(OFFS)):
    if r["method"] == "cfd50":
        off[int(r["ch"])] = float(r["offset_ns"])
trawS = cfd50(e["plS"], e["wS"], e["kS"])
tS = trawS - off[None, :]
tC = cfd50(e["plC"], e["wC"], e["kC"]) - off[None, :] - DC
aS, aC = e["aS"].astype(float), e["aC"].astype(float)
isX = np.zeros(32, bool)
isX[XCH] = True
isY = np.zeros(32, bool)
isY[YCH] = True
LIVE = isX | isY
XI, YI = np.array(XCH), np.array(YCH)
prv, nxt = np.arange(0, N - 1), np.arange(1, N)
ar = np.arange(N)
rng = np.random.default_rng(SEED)


def dec_median(above, t, W):
    tt = np.where(above & np.isfinite(t), t, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ref = np.nanmedian(tt, 1)
    inw = np.isfinite(tt) & (np.abs(tt - ref[:, None]) < W)
    return inw[:, isX].any(1) & inw[:, isY].any(1), ref


def dec_pair(above, t, W):
    tt = np.where(above & np.isfinite(t), t, np.inf)
    out = np.zeros(len(t), bool)
    for s in range(0, len(t), 20000):  # por bloques, para no reservar (n, 15, 16) de una vez
        tx, ty = tt[s:s + 20000][:, XI], tt[s:s + 20000][:, YI]
        with np.errstate(invalid="ignore"):
            d = np.abs(tx[:, :, None] - ty[:, None, :])
        out[s:s + 20000] = (d < W).any((1, 2))
    return out, None


DEC = {"median": dec_median, "pair": dec_pair}

# muones etiquetados
pe = np.where(LIVE[None, :], aS / G[None, :], 0.0)
ix = XI[np.argmax(pe[:, XCH], 1)]
iy = YI[np.argmax(pe[:, YCH], 1)]
tXt, tYt = tS[ar, ix], tS[ar, iy]
tag = (pe[ar, ix] > TAG) & (pe[ar, iy] > TAG) & np.isfinite(tXt) & np.isfinite(tYt)
tag_c = tag & (np.abs(tXt - tYt) < 5.0)
TMU = 0.5 * (tXt + tYt)[tag_c]
TXY = np.stack([tXt[tag_c], tYt[tag_c]], 1)
DELTA = float(np.median((e["kS"].astype(float) * DS - trawS)[(pe > 5) & np.isfinite(trawS)]))
TLO = S0 * DS - DELTA - off
MUON = {"tagged": int(tag.sum()), "tagged_dt_lt_5ns": int(tag_c.sum()), "t_mu_median_ns": float(np.median(TMU)), "t_mu_std_ns": float(TMU.std()),
        "delta_kmin_minus_cfd_ns": DELTA}


def embed(offab, dec, W, cfg, noise):
    ab = offab.copy()
    if noise == "measured":
        tt = tC.copy()
        ab &= np.isfinite(tt)
    else:
        tt = TLO[None, :] + SPAN_NS * rng.random((N, 32))
    cx, cyy = rng.choice(XCH, N), rng.choice(YCH, N)
    if cfg in ("XY", "X", "Y"):
        tx = ty = TMU[rng.integers(0, len(TMU), N)]
    elif cfg == "mm":
        tx = TXY[rng.integers(0, len(TXY), N), 0]
        ty = TXY[rng.integers(0, len(TXY), N), 1]
    if cfg in ("XY", "X", "mm"):
        ab[ar, cx], tt[ar, cx] = True, tx
    if cfg in ("XY", "Y", "mm"):
        ab[ar, cyy], tt[ar, cyy] = True, ty
    if dec is None:
        return float((ab[:, isX].any(1) & ab[:, isY].any(1)).mean())
    return float(DEC[dec](ab, tt, W)[0].mean())


fg = F18 * G_XY


def finish(o):
    ev = o["q_mix"] / fg
    o["eps_XY"] = {"value": ev, "stat": o["q_mix_err"] / fg, "syst_estimator": abs(o["q_mix"] - o["q_ctl"]) / fg, "syst_acceptance": ev * ACC_REL}
    return o


res = {"campaign": LABEL, "N": int(N), "muons": MUON, "f_2018": F18, "g": G_XY, "acceptance_rel": ACC_REL, "seed": SEED, "offsets": OFFS, "per_T": {}}
rows = []
for T in TS:
    thr = T * G
    above = (aS > thr[None, :]) & LIVE[None, :]
    offab = (aC > thr[None, :]) & LIVE[None, :]
    Xs, Ys = above[:, isX].any(1), above[:, isY].any(1)
    r = xy.xy_fraction(Xs, Ys, offab[:, isX].any(1), offab[:, isY].any(1))
    ref = S2["datasets"][S2KEY]["per_T"][str(T)]
    out = {"144ns": finish({"N_XY": r["N_XY"], "acc_mix": r["acc_mix"], "D": (1 - r["aX"]) * (1 - r["aY"]), "rX": r["rX"], "rY": r["rY"],
                            "q_mix": r["f_mix"], "q_mix_err": r["f_mix_err"], "q_ctl": r["f_ctl"], "noise": "-",
                            "check_dq_mix_vs_stage2": r["f_mix"] - ref["f_mix"], "check_dq_ctl_vs_stage2": r["f_ctl"] - ref["f_ctl"]})}
    rX, rY = r["rX"], r["rY"]
    out["_cfd_nan_fraction_above_T"] = {"signal": float((above & ~np.isfinite(tS)).sum() / above.sum()),
                                        "offtime": float((offab & ~np.isfinite(tC)).sum() / max(offab.sum(), 1))}
    for dname, W in [(None, None)] + [(d, w) for d in ("median", "pair") for w in WS]:
        key = "144ns_general" if dname is None else f"{dname}{int(W)}"
        if dname is None:
            Sv, m = r["N_XY"], [float((Xs[:-1] & Ys[1:]).sum()), float((Ys[:-1] & Xs[1:]).sum())]
        else:
            pas, tref = DEC[dname](above, tS, W)
            Sv, m = float(pas.sum()), []
            for a, b in ((prv, nxt), (nxt, prv)):
                abm = np.where(isX[None, :], above[a], above[b])
                ttm = np.where(isX[None, :], tS[a], tS[b])
                m.append(float(DEC[dname](abm, ttm, W)[0].sum()))
        for noise in (("measured", "uniform") if dname is not None else ("measured",)):
            A = {c: embed(offab, dname, W, c, noise) for c in ("XY", "X", "Y", "0", "mm")}
            D = A["XY"] - A["X"] - A["Y"] + A["0"]
            acc_mix = 0.5 * (m[0] + m[1]) * N / (N - 1)
            Sf, Mf = Sv / N, acc_mix / N
            o = {"N_XY": Sv, "m1": m[0], "m2": m[1], "acc_mix": acc_mix, "A": A, "D": D, "rX": rX, "rY": rY, "noise": noise,
                 "q_mix": (Sf - Mf + rX * rY * (A["mm"] - A["X"] - A["Y"] + A["0"])) / D,
                 "q_mix_err": math.sqrt(Sv + (m[0] + m[1]) / 4) / (N * D),
                 "q_ctl": (Sf - rX * A["X"] - rY * A["Y"] - (1 - rX - rY) * A["0"]) / D}
            if dname == "median" and noise == "measured":
                sel = np.nonzero(tag)[0]
                excl = np.zeros((len(sel), 32), bool)
                for d in (-1, 0, 1):
                    for arr in (ix, iy):
                        cc = arr[sel] + d
                        okc = (cc >= 0) & (cc < 32) & ((cc < 16) == (arr[sel] < 16))
                        excl[np.nonzero(okc)[0], cc[okc]] = True
                dd = tS[sel] - tref[sel][:, None]
                dn = tC[sel] - tref[sel][:, None]
                o["check_density_tagged"] = {"n_tagged": int(len(sel)),
                                             "data_per_ns_per_event": float((above[sel] & ~excl & (dd >= 24) & (dd < 72)).sum() / len(sel) / 48.0),
                                             "offtime_measured_per_ns_per_event": float((offab[sel] & ~excl & (dn >= 24) & (dn < 72)).sum() / len(sel) / 48.0)}
            out[key if noise == "measured" else key + "_uniform"] = finish(o)
    res["per_T"][str(T)] = out
    for k, o in out.items():
        if k.startswith("_"):
            continue
        e_ = o["eps_XY"]
        rows.append([LABEL, T, k, o["noise"], N, f"{o['N_XY']:.0f}", f"{o['acc_mix']:.2f}", f"{o['D']:.5f}", f"{o['rX']:.5f}", f"{o['rY']:.5f}",
                     f"{o['q_mix']:.6f}", f"{o['q_mix_err']:.6f}", f"{o['q_ctl']:.6f}", f"{e_['value']:.4f}", f"{e_['stat']:.4f}",
                     f"{e_['syst_estimator']:.4f}", f"{e_['syst_acceptance']:.4f}"])
        extra = (" A=" + " ".join(f"{c}:{v:.4f}" for c, v in o["A"].items())) if "A" in o else \
            f" | dq vs stage2 {o['check_dq_mix_vs_stage2']:.1e} {o['check_dq_ctl_vs_stage2']:.1e}"
        if "check_density_tagged" in o:
            cd = o["check_density_tagged"]
            extra += f" | dens data {cd['data_per_ns_per_event']:.3e} offtime {cd['offtime_measured_per_ns_per_event']:.3e}"
        print(f"T={T:2d} {k:22s} N_XY={o['N_XY']:7.0f} acc_mix={o['acc_mix']:9.2f} D={o['D']:.4f} q={o['q_mix']:.5f}±{o['q_mix_err']:.5f} "
              f"q_ctl={o['q_ctl']:.5f} eps={e_['value']:.4f} ± {e_['stat']:.4f} ± {e_['syst_estimator']:.4f} ± {e_['syst_acceptance']:.4f}" + extra,
              flush=True)
    print(f"T={T:2d} fracción sin cruce CFD: {out['_cfd_nan_fraction_above_T']}", flush=True)
with open(f"{OUT}/cfd_window_eff2_{CAMP}.csv", "w") as fh:
    fh.write("# cfd_window_eff2.py (2026-11-15, version 2): hod2018 eps_XY = q_mix/(f_2018 g), g = 0.9551; definitions 144ns (xy_estimator = stage2), "
             "144ns_general, medianW (|t'-t_ref| < W, t_ref = per-event median), pairW (an X and a Y strip with |t'_X - t'_Y| < W); N_acc = mixed events "
             "with the same definition; window-passing probabilities A from embedded off-time pulses with measured CFD times (noise = measured) or uniform "
             "times (noise = uniform, version 1 method); errors stat, estimator |q_mix - q_ctl|, acceptance eps*0.046481\n")
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["campaign", "T_pe", "definition", "noise_times", "N", "N_XY", "N_acc_mix", "D", "rX", "rY", "q_mix", "q_mix_err", "q_ctl", "eps_XY",
                "err_stat", "err_estimator", "err_acceptance"])
    w.writerows(rows)
json.dump(res, open(f"{OUT}/cfd_window_eff2_{CAMP}.json", "w"), indent=1)
print(json.dumps(MUON))
