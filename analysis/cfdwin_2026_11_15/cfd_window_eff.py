#!/usr/bin/python3 -I
# cfd_window_eff.py — seguimiento 2026-11-15 (ventana CFD): ε_XY de hod2018 con una ventana de coincidencia sobre los tiempos CFD en lugar de la
# ventana signal de 144 ns (mínimo del pulso en [49,58) del tramo). Solo lectura de salidas existentes.
# Entradas (mayo, run B en tiempo bueno de may_common.good_time, 77 558 disparos):
#  - work/may_v8_time.bin del seguimiento 13 (time_extract.C): por canal el pulso más alto con el mínimo en la ventana signal (aS, kS, p_loc y
#    x[k−7 … k+2]); tiempo = cfd50 de time_res.py (fracción 0.5, pedestal local p_loc, interpolación lineal, 16 ns/muestra), corregido por el
#    desplazamiento por canal del seguimiento 13 A (A_time/time_res_offsets.csv, método cfd50, columna offset_ns): t' = t − o[c].
#  - fuera de tiempo: aC de work/hd64_may_v8_hits.bin (hitdef_extract.C, la entrada de stage2.py), unido por t evento a evento (se comprueba que aS
#    coincide en los dos archivos).
#  - G = reajuste de PROC (hd64_may_v8_peaks.csv, ok/changed; mediana para los demás), ch3 fuera; f_2018 de stage2.json; g = 0.9551211 e
#    incertidumbre de aceptancia relativa 0.046481 (f_2018 ⊕ g) de g_dead_ch3.json, como la figura g09551 (0.81 ± 0.09 ± 0.07 ± 0.04 a 5 p.e.).
# Ventana: una tira (amplitud > T·G) cuenta solo si su tiempo CFD es finito y |t' − t_ref| < W (W = 8, 16 ns); t_ref = mediana por evento de los t'
# finitos de todas las tiras de hod2018 sobre el umbral (X e Y). [Elección: el seguimiento 13 A no tiene referencia del disparo (solo Δt = t_X − t_Y),
# y una referencia fija del disparo no sirve: el tiempo medio de los muones etiquetados tiene σ = 27 ns respecto del disparo.] Con una tira X y una Y,
# t_ref es el punto medio y la condición es |t_X − t_Y| < 2W.
# Estimador. Con t_ref sacada del propio evento, la decisión X·Y ya no es el producto de decisiones por plano y la fórmula de xy_estimator
# ((N_XY − acc_mix)/(N(1−aX)(1−aY)) + rX·rY) deja de valer (un plano con un solo pulso de ruido se fija su propia referencia). Generalización,
# con las clases por evento (muón verdadero en X y en Y; solo en X; solo en Y; ninguno) y la probabilidad A de pasar X·Y con ventana en cada clase:
#   S = N_XY/N = q·A_XY + (rX − q)·A_X + (rY − q)·A_Y + (1 − rX − rY + q)·A_0
#   M = acc_mix/N (mezclados con la misma ventana: tiras X del evento i e Y del i+1, y al revés; t_ref recalculada en el seudo-evento)
#     = rX·rY·A_mm + rX(1 − rY)·A_X + (1 − rX)rY·A_Y + (1 − rX)(1 − rY)·A_0      (A_mm: muones X e Y de eventos distintos)
#   q_mix = [S − M + rX·rY·(A_mm − A_X − A_Y + A_0)] / D,   D = A_XY − A_X − A_Y + A_0
#   q_ctl = [S − rX·A_X − rY·A_Y − (1 − rX − rY)·A_0] / D      (sin mezclados: el control para el término del estimador)
# Con la ventana de 144 ns, A_XY = 1, A_X = aY, A_Y = aX, A_0 = aX·aY, A_mm = 1 y q_mix es la fórmula de xy_estimator.
#  - rX, rY: probabilidad verdadera por plano de la ventana de 144 ns (xy_fraction); la presencia de un muón en un plano no depende de la ventana.
#  - A_*: incrustación con los datos [S]: para cada evento, el patrón real de pulsos fuera de tiempo sobre T·G (aC) como ruido, cada uno con tiempo
#    uniforme en los 144 ns de la ventana (t' ∈ [784 − δ − o[c], 928 − δ − o[c]), δ = mediana de 16·kS − t_cfd de las tiras sobre el umbral); los
#    muones verdaderos en un canal X y uno Y al azar con el tiempo medio (t'_X + t'_Y)/2 de un muón etiquetado al azar (X, Y > 10 p.e. y
#    |t'_X − t'_Y| < 5 ns), los dos en el mismo tiempo para A_XY (q cuenta los X·Y verdaderos que pasan; la pérdida por resolución queda en ε);
#    para A_mm, t'_X y t'_Y de dos muones etiquetados distintos. Una pasada por los N eventos, semilla 2026111501.
#  - Control del supuesto uniforme: en los muones etiquetados, densidad de tiras sobre T (fuera de las tiras etiquetadas y sus vecinas) con t' − t_ref
#    en [24, 72) ns, contra la del ruido incrustado en el mismo intervalo.
# Errores: stat = √(N_XY + (m1 + m2)/4)/(N·D)/(f_2018·g); estimador = |q_mix − q_ctl|/(f_2018·g); aceptancia = ε·0.046481.
# Para 144 ns se da también el resultado de xy_estimator.xy_fraction (debe coincidir con stage2.json; se comprueba).
# Salida: <out>/cfd_window_eff.csv, cfd_window_eff.json. Uso: cfd_window_eff.py may <out>
import csv, json, math, os, sys, warnings

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

CAMP, OUT = sys.argv[1], sys.argv[2]
assert CAMP == "may"
os.makedirs(OUT, exist_ok=True)
F13 = os.path.join(REPO, "outputs/hod2018_timing/2026-11-03_followup13")
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
GD = json.load(open(os.path.join(REPO, "outputs/proceeding_revision/2026-11-02_revision/g_dead_ch3/g_dead_ch3.json")))
S2 = json.load(open(PK + "/stage2.json"))
F18, G_XY, ACC_REL = S2["f_2018"], GD["g_used_new"]["value"], GD["acceptance_rel_syst_no_ch3"]
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("kS", "<u2", 32), ("bl", "<f4", 32), ("pl", "<f4", 32), ("w", "<f4", (32, 10))])
XCH = [c for c in range(16) if c != 3]
YCH = list(range(16, 32))
TS = [3, 5, 7, 10]
WS = [8.0, 16.0]
SPAN_NS, KPRE, DS, S0 = 144.0, 7, 16.0, 49
TAG = 10.0
SEED = 2026111501


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


def cfd50(e):
    """t (ns desde el inicio del tramo) del cfd50 de time_res.py, (N, 32); NaN sin cruce o sin pulso."""
    w = e["w"].astype(float)
    a = e["pl"].astype(float)[:, :, None] - w
    A = a[:, :, KPRE]
    below = a[:, :, :KPRE] < 0.5 * A[:, :, None]
    last = np.where(below.any(2), KPRE - 1 - np.argmax(below[:, :, ::-1], 2), -1)
    li = np.clip(last, 0, KPRE - 1)
    a0 = np.take_along_axis(a, li[:, :, None], 2)[:, :, 0]
    a1 = np.take_along_axis(a, (li + 1)[:, :, None], 2)[:, :, 0]
    frac = (0.5 * A - a0) / np.where(a1 != a0, a1 - a0, 1)
    t = (e["kS"].astype(float) - KPRE + li + frac) * DS
    return np.where((last >= 0) & (e["kS"] != 0xFFFF), t, np.nan)


G = g_refit()
e = np.fromfile(F13 + "/work/may_v8_time.bin", dtype=DT)
e = e[np.argsort(e["t"], kind="stable")]
e = e[mc.good_time(e["t"], "v8")]
h = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
h = h[mc.good_time(h["t"], "v8")]
assert len(e) == len(h) and np.array_equal(e["t"], h["t"]), "eventos distintos en los dos archivos"
N = len(e)
JOIN = {"N": int(N), "aS_max_abs_diff_mV": float(np.abs(e["aS"] - h["aS"][:, :32]).max())}
off = np.zeros(32)
for r in csv.DictReader(open(F13 + "/A_time/time_res_offsets.csv")):
    if r["method"] == "cfd50":
        off[int(r["ch"])] = float(r["offset_ns"])
traw = cfd50(e)
tc = traw - off[None, :]
aS, aC = e["aS"].astype(float), h["aC"][:, :32].astype(float)
isX = np.zeros(32, bool)
isX[XCH] = True
isY = np.zeros(32, bool)
isY[YCH] = True
LIVE = isX | isY
prv, nxt = np.arange(0, N - 1), np.arange(1, N)
rng = np.random.default_rng(SEED)


def window_pass(above, t, W):
    """(X·Y con ventana, tiras en la ventana, t_ref) para above, t (n, 32)."""
    tt = np.where(above & np.isfinite(t), t, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ref = np.nanmedian(tt, 1)
    inw = np.isfinite(tt) & (np.abs(tt - ref[:, None]) < W)
    return inw[:, isX].any(1) & inw[:, isY].any(1), inw, ref


# muones etiquetados (para t_μ y el control): tira X y tira Y más altas, las dos > 10 p.e.
pe = np.where(LIVE[None, :], aS / G[None, :], 0.0)
ix = np.array(XCH)[np.argmax(pe[:, XCH], 1)]
iy = np.array(YCH)[np.argmax(pe[:, YCH], 1)]
ar = np.arange(N)
tXt, tYt = tc[ar, ix], tc[ar, iy]
tag = (pe[ar, ix] > TAG) & (pe[ar, iy] > TAG) & np.isfinite(tXt) & np.isfinite(tYt)
tag_c = tag & (np.abs(tXt - tYt) < 5.0)
TMU = 0.5 * (tXt + tYt)[tag_c]
TXY = np.stack([tXt[tag_c], tYt[tag_c]], 1)
# δ: 16·kS − t_cfd crudo de las tiras sobre 5 p.e.
ab5 = (pe > 5) & np.isfinite(traw)
DELTA = float(np.median((e["kS"].astype(float) * DS - traw)[ab5]))
TLO = S0 * DS - DELTA - off  # límite inferior del t' de un pulso con el mínimo en la ventana, por canal
MUON = {"tagged": int(tag.sum()), "tagged_dt_lt_5ns": int(tag_c.sum()), "t_mu_median_ns": float(np.median(TMU)), "t_mu_std_ns": float(TMU.std()),
        "t_mu_p16_p84_ns": [float(x) for x in np.percentile(TMU, [16, 84])], "delta_kmin_minus_cfd_ns": DELTA}


def embed(offab, W, cfg):
    """Probabilidad A de pasar X·Y con ventana (W en ns; W = None: ventana de 144 ns sin corte CFD) para la clase cfg ∈ XY, X, Y, 0, mm,
    con el ruido = patrón fuera de tiempo de cada evento y tiempos uniformes."""
    ab = offab.copy()
    tt = TLO[None, :] + SPAN_NS * rng.random((N, 32))
    cx = rng.choice(XCH, N)
    cyy = rng.choice(YCH, N)
    if cfg in ("XY", "X", "Y"):
        tm = TMU[rng.integers(0, len(TMU), N)]
        tx = ty = tm
    elif cfg == "mm":
        tx = TXY[rng.integers(0, len(TXY), N), 0]
        ty = TXY[rng.integers(0, len(TXY), N), 1]
    if cfg in ("XY", "X", "mm"):
        ab[ar, cx], tt[ar, cx] = True, tx
    if cfg in ("XY", "Y", "mm"):
        ab[ar, cyy], tt[ar, cyy] = True, ty
    if W is None:
        return float((ab[:, isX].any(1) & ab[:, isY].any(1)).mean())
    return float(window_pass(ab, tt, W)[0].mean())


res = {"campaign": "may_v8_goodtime", "join": JOIN, "muons": MUON, "f_2018": F18, "g": G_XY, "acceptance_rel": ACC_REL, "seed": SEED,
       "t_ref": "per-event median of t' (CFD 50 %, followup 13 A cfd50 offsets) of all hod2018 strips above T with finite time, X and Y, ch3 out",
       "per_T": {}}
rows = []
fg = F18 * G_XY


def finish(o):
    ev = o["q_mix"] / fg
    o["eps_XY"] = {"value": ev, "stat": o["q_mix_err"] / fg, "syst_estimator": abs(o["q_mix"] - o["q_ctl"]) / fg, "syst_acceptance": ev * ACC_REL}
    return o


for T in TS:
    thr = T * G
    above = (aS > thr[None, :]) & LIVE[None, :]
    offab = (aC > thr[None, :]) & LIVE[None, :]
    out = {}
    # 144 ns, xy_fraction (la definición de stage2)
    Xs, Ys = above[:, isX].any(1), above[:, isY].any(1)
    r = xy.xy_fraction(Xs, Ys, offab[:, isX].any(1), offab[:, isY].any(1))
    ref = S2["datasets"]["may8_good"]["per_T"][str(T)]
    o = {"N": N, "N_XY": r["N_XY"], "acc_mix": r["acc_mix"], "acc_ctl": r["acc_ctl"], "aX": r["aX"], "aY": r["aY"], "rX": r["rX"], "rY": r["rY"],
         "q_mix": r["f_mix"], "q_mix_err": r["f_mix_err"], "q_ctl": r["f_ctl"], "check_dq_mix_vs_stage2": r["f_mix"] - ref["f_mix"],
         "check_dq_ctl_vs_stage2": r["f_ctl"] - ref["f_ctl"]}
    out["144ns"] = finish(o)
    rX, rY = r["rX"], r["rY"]
    nan_frac = float((above & ~np.isfinite(tc)).sum() / above.sum())
    for W in [None] + WS:
        key = "144ns_general" if W is None else f"cfd{int(W)}"
        if W is None:
            Sv = r["N_XY"]
            m = [float((Xs[:-1] & Ys[1:]).sum()), float((Ys[:-1] & Xs[1:]).sum())]
        else:
            pas, _, tref = window_pass(above, tc, W)
            Sv = float(pas.sum())
            m = []
            for a, b in ((prv, nxt), (nxt, prv)):
                abm = np.where(isX[None, :], above[a], above[b])
                ttm = np.where(isX[None, :], tc[a], tc[b])
                m.append(float(window_pass(abm, ttm, W)[0].sum()))
        A = {c: embed(offab, W, c) for c in ("XY", "X", "Y", "0", "mm")}
        D = A["XY"] - A["X"] - A["Y"] + A["0"]
        acc_mix = 0.5 * (m[0] + m[1]) * N / (N - 1)
        Sf, Mf = Sv / N, acc_mix / N
        q_mix = (Sf - Mf + rX * rY * (A["mm"] - A["X"] - A["Y"] + A["0"])) / D
        q_ctl = (Sf - rX * A["X"] - rY * A["Y"] - (1 - rX - rY) * A["0"]) / D
        o = {"N": N, "N_XY": Sv, "m1": m[0], "m2": m[1], "acc_mix": acc_mix, "acc_ctl_model": N * (rX * A["X"] + rY * A["Y"] + (1 - rX - rY) * A["0"]),
             "A": A, "D": D, "rX": rX, "rY": rY, "q_mix": q_mix, "q_mix_err": math.sqrt(Sv + (m[0] + m[1]) / 4) / (N * D), "q_ctl": q_ctl,
             "cfd_nan_fraction_above_T": nan_frac}
        if W is not None:
            # control del supuesto uniforme en muones etiquetados: tiras sobre T fuera de las etiquetadas y vecinas, t' − t_ref ∈ [24, 72) ns
            sel = np.nonzero(tag)[0]
            ab = above[sel].copy()
            for d in (-1, 0, 1):
                for arr in (ix, iy):
                    cc = arr[sel] + d
                    okc = (cc >= 0) & (cc < 32) & ((cc < 16) == (arr[sel] < 16))
                    ab[np.nonzero(okc)[0], cc[okc]] = False
            dd = tc[sel] - tref[sel][:, None]
            dens_data = float((ab & (dd >= 24) & (dd < 72)).sum() / len(sel) / 48.0)
            tt = TLO[None, :] + SPAN_NS * rng.random((len(sel), 32))
            nab = offab[sel]
            dn = tt - tref[sel][:, None]
            # misma exclusión de canales para el ruido incrustado
            excl = above[sel] & ~ab
            dens_emb = float((nab & ~excl & (dn >= 24) & (dn < 72)).sum() / len(sel) / 48.0)
            o["check_uniform_tagged"] = {"n_tagged": int(len(sel)), "data_density_per_ns_per_event": dens_data, "embedded_density_per_ns_per_event": dens_emb}
        out[key] = finish(o)
    res["per_T"][str(T)] = out
    for k, o in out.items():
        e_ = o["eps_XY"]
        rows.append(["may_v8_goodtime", T, k, N, f"{o['N_XY']:.0f}", f"{o['acc_mix']:.2f}", f"{o.get('D', (1 - o.get('aX', 0)) * (1 - o.get('aY', 0))):.5f}",
                     f"{o['rX']:.5f}", f"{o['rY']:.5f}", f"{o['q_mix']:.6f}", f"{o['q_mix_err']:.6f}", f"{o['q_ctl']:.6f}", f"{e_['value']:.4f}",
                     f"{e_['stat']:.4f}", f"{e_['syst_estimator']:.4f}", f"{e_['syst_acceptance']:.4f}"])
        extra = ""
        if "A" in o:
            extra = " A=" + " ".join(f"{c}:{v:.4f}" for c, v in o["A"].items()) + f" D={o['D']:.4f}"
        if "check_uniform_tagged" in o:
            cu = o["check_uniform_tagged"]
            extra += f" | dens data {cu['data_density_per_ns_per_event']:.3e} emb {cu['embedded_density_per_ns_per_event']:.3e} nan {o['cfd_nan_fraction_above_T']:.4f}"
        if k == "144ns":
            extra = f" | dq vs stage2 {o['check_dq_mix_vs_stage2']:.1e} {o['check_dq_ctl_vs_stage2']:.1e}"
        print(f"T={T:2d} {k:13s} N_XY={o['N_XY']:7.0f} acc_mix={o['acc_mix']:9.2f} q={o['q_mix']:.5f}±{o['q_mix_err']:.5f} q_ctl={o['q_ctl']:.5f} "
              f"eps={e_['value']:.4f} ± {e_['stat']:.4f} ± {e_['syst_estimator']:.4f} ± {e_['syst_acceptance']:.4f}" + extra, flush=True)
with open(f"{OUT}/cfd_window_eff.csv", "w") as fh:
    fh.write("# cfd_window_eff.py (2026-11-15): hod2018 eps_XY = q_mix/(f_2018 g), g = 0.9551, f_2018 from stage2.json; window 144ns = signal window "
             "(xy_estimator, = stage2), 144ns_general = same with the generalized estimator, cfdW = CFD window |t'-t_ref| < W ns, t_ref = per-event median "
             "of CFD times above T; N_acc = mixed events with the same window; D = A_XY - A_X - A_Y + A_0 (144ns: (1-aX)(1-aY)); "
             "errors: stat, estimator |q_mix - q_ctl|, acceptance eps*0.046481\n")
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["campaign", "T_pe", "window", "N", "N_XY", "N_acc_mix", "D", "rX", "rY", "q_mix", "q_mix_err", "q_ctl", "eps_XY", "err_stat",
                "err_estimator", "err_acceptance"])
    w.writerows(rows)
json.dump(res, open(f"{OUT}/cfd_window_eff.json", "w"), indent=1)
print(json.dumps(JOIN), json.dumps(MUON))
