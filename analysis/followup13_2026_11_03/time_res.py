#!/usr/bin/python3 -I
# time_res.py — seguimiento 2026-11-03 (13), punto A: resolución temporal de hod2018 en run B (mayo v8, tiempo bueno de may_common.good_time).
# Entrada: <work>/may_v8_time.bin de time_extract.C (pulso más alto de la ventana signal por canal: aS, kS, bl, p_loc, x[k−7 … k+2]).
# Muón etiquetado: tira X de mayor aS/G y tira Y de mayor aS/G, las dos > 10 p.e. medidos (G = reajuste de PROC, ch3 muerto fuera).
# Mezclados (accidentales): X del evento i con Y del evento i+1, misma selección, × N/(N−1); los histogramas son mismo − mezclados.
# Tiempos del pulso (ns desde el inicio del tramo):
#   cfd50      — fracción constante 0.5 en el frente de subida: A = p_loc − x[k]; última muestra i < k con p_loc − x[i] < 0.5 A, interpolación
#                lineal entre i e i+1 (16 ns/muestra); p_loc = media de [k−12, k−4] (pedestal local de charge_extract.C).
#   cfd50_bl   — igual con el pedestal del tramo bl (media recortada de [0,40)) en lugar de p_loc.
#   kmin       — muestra del mínimo · 16 ns.
# Correcciones de Δt = t_X − t_Y: "none"; "offsets" = desplazamiento por canal o_X[c], o_Y[c] (mínimos cuadrados alternados con mediana sobre
# todos los etiquetados, Σ o_Y = 0); "offsets+prop" = además s · (d_X − d_Y), d = distancia del cruce al extremo de lectura de la tira (la
# posición a lo largo de la barra X la da la tira Y y al revés; centros d = 22.5 + 2j mm desde el extremo del índice 0, extremo de lectura por
# capa de layer_from_channels.json, como bar_targets.py) y s ajustada junto con los desplazamientos (ns/mm).
# Por bin de la menor amplitud (10-15, 15-20, 20-30, > 30 p.e. medidos, y todos): ajuste gaussiano (LM, bins de 1 ns, o 16 ns para kmin, rango
# ±2.5 σ iterado) del histograma neto, RMS del neto dentro de ±3 σ del ajuste alrededor de μ, σ por tira = σ/√2 (contribuciones X e Y iguales).
# Salida: <out>/time_res.csv, time_res_hist.csv (histogramas netos), time_res_offsets.csv, time_res.json. Uso: time_res.py <work> <out>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

WK, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("kS", "<u2", 32), ("bl", "<f4", 32), ("pl", "<f4", 32), ("w", "<f4", (32, 10))])
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ORD = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD, TAG, KPRE, DS = 3, 10.0, 7, 16.0
BINS = [("10-15", 10, 15), ("15-20", 15, 20), ("20-30", 20, 30), (">30", 30, 1e9), ("all", 10, 1e9)]


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


G = g_refit()
e = np.fromfile(f"{WK}/may_v8_time.bin", dtype=DT)
e = e[mc.good_time(e["t"], "v8")]
N = len(e)
pe = e["aS"] / G
pe[:, DEAD] = 0
xo, yo = np.asarray(ORD["h18_X"]), np.asarray(ORD["h18_Y"])
posX = {c: i for i, c in enumerate(xo)}
posY = {c: i for i, c in enumerate(yo)}
cx = xo[np.argmax(pe[:, xo], 1)]
cy = yo[np.argmax(pe[:, yo], 1)]
ax, ay = pe[np.arange(N), cx], pe[np.arange(N), cy]


def times(method, ch):
    """Tiempo (ns) del pulso del canal ch[i] en el evento i, para todos los eventos."""
    i = np.arange(N)
    k = e["kS"][i, ch].astype(float)
    if method == "kmin":
        return k * DS
    w = e["w"][i, ch]  # (N, 10), w[:, 7] = mínimo
    ped = e["pl"][i, ch] if method == "cfd50" else e["bl"][i, ch]
    a = ped[:, None] - w
    A = a[:, KPRE]
    below = a[:, :KPRE] < 0.5 * A[:, None]
    # última muestra antes del mínimo con a < A/2
    last = np.where(below.any(1), KPRE - 1 - np.argmax(below[:, ::-1], 1), -1)
    ok = last >= 0
    li = np.clip(last, 0, KPRE - 1)
    a0, a1 = a[i, li], a[i, li + 1]
    frac = np.where(ok, (0.5 * A - a0) / np.where(a1 != a0, a1 - a0, 1), np.nan)
    return (k - KPRE + li + frac) * DS


def readout_d(probe, strip_pos_other, ch_probe):
    """Distancia (mm) del cruce al extremo de lectura de la tira ch_probe (coordenada probe), con la posición a lo largo dada por la otra."""
    layer = np.where(np.array([(posX if probe == "X" else posY)[c] % 2 for c in ch_probe]) == 0, "A", "B")
    sgn = np.array([SIGN[f"may2026_v8_goodtime/{probe}{L}"]["sign"] for L in layer])
    D = 22.5 + 2.0 * strip_pos_other
    return np.where(sgn < 0, D, 75.0 - D)


jX = np.array([posX[c] for c in cx])
jY = np.array([posY[c] for c in cy])
dX = readout_d("X", jY, cx)  # tira X: posición a lo largo = índice de la tira Y
dY = readout_d("Y", jX, cy)

sel = (ax > TAG) & (ay > TAG)
nxt = np.roll(np.arange(N), -1)
selm = (ax > TAG) & (ay[nxt] > TAG)
selm[-1] = False
scale = N / (N - 1)


def lm_gauss(xc, y, err, p0, iters=60):
    p = np.array(p0, float)
    lam = 1e-3

    def f(p):
        return p[0] * np.exp(-0.5 * ((xc - p[1]) / p[2]) ** 2)

    def chi2(p):
        return float((((y - f(p)) / err) ** 2).sum())

    c0 = chi2(p)
    for _ in range(iters):
        g = np.exp(-0.5 * ((xc - p[1]) / p[2]) ** 2)
        J = np.stack([g, p[0] * g * (xc - p[1]) / p[2] ** 2, p[0] * g * (xc - p[1]) ** 2 / p[2] ** 3], 1) / err[:, None]
        r = (y - f(p)) / err
        H = J.T @ J
        step = np.linalg.solve(H + lam * np.diag(np.diag(H)), J.T @ r)
        pn = p + step
        pn[2] = abs(pn[2])
        cn = chi2(pn)
        if cn < c0:
            p, c0, lam = pn, cn, lam * 0.3
            if np.abs(step).max() < 1e-7:
                break
        else:
            lam *= 10
    g = np.exp(-0.5 * ((xc - p[1]) / p[2]) ** 2)
    J = np.stack([g, p[0] * g * (xc - p[1]) / p[2] ** 2, p[0] * g * (xc - p[1]) ** 2 / p[2] ** 3], 1) / err[:, None]
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov = np.full((3, 3), np.nan)
    return p, np.sqrt(np.abs(np.diag(cov))), c0, len(xc) - 3


def fit_hist(d_same, d_mix, bw, lim):
    edges = np.arange(-lim, lim + bw / 2, bw) - bw / 2 if bw == DS else np.arange(-lim, lim + bw / 2, bw)
    hs = np.histogram(d_same, edges)[0].astype(float)
    hm = np.histogram(d_mix, edges)[0].astype(float) * scale
    net = hs - hm
    err = np.sqrt(np.maximum(hs + scale * hm, 1.0))
    xc = 0.5 * (edges[1:] + edges[:-1])
    mu, sg = float(np.median(d_same)), float(1.4826 * np.median(np.abs(d_same - np.median(d_same))) + 0.5)
    p = [net.max(), mu, max(sg, bw / 2)]
    res = None
    for _ in range(6):
        m = np.abs(xc - p[1]) <= 2.5 * max(p[2], bw)
        if m.sum() < 4:
            break
        res = lm_gauss(xc[m], net[m], err[m], p)
        p = list(res[0])
    if res is None:
        return None, edges, hs, hm
    (A, mu, sg), (eA, emu, esg), c2, ndf = res
    m3 = np.abs(xc - mu) <= 3 * sg
    w = np.clip(net[m3], 0, None)
    rms = float(np.sqrt((w * (xc[m3] - mu) ** 2).sum() / w.sum())) if w.sum() > 0 else float("nan")
    return dict(mu=float(mu), mu_err=float(emu), sigma=float(sg), sigma_err=float(esg), chi2=float(c2), ndf=int(ndf), rms3=rms,
                n_same=int(hs.sum()), n_mix=float(hm.sum()), n_net=float(net.sum())), edges, hs, hm


def solve_offsets(dt, cxs, cys, dd, with_prop, iters=50):
    oX, oY, s = np.zeros(32), np.zeros(32), 0.0
    ok = np.isfinite(dt)
    for _ in range(iters):
        r = dt - s * dd
        for c in np.unique(cxs[ok]):
            m = ok & (cxs == c)
            oX[c] = np.median(r[m] + oY[cys[m]])
        for c in np.unique(cys[ok]):
            m = ok & (cys == c)
            oY[c] = np.median(oX[cxs[m]] - r[m])
        sh = oY[np.unique(cys[ok])].mean()
        oY[np.unique(cys[ok])] -= sh
        oX[np.unique(cxs[ok])] -= sh
        if with_prop:
            rr = dt - oX[cxs] + oY[cys]
            mad = max(1.4826 * np.median(np.abs(rr[ok] - np.median(rr[ok]))), DS / 2)  # piso: kmin es discreto (MAD = 0)
            mm = ok & (np.abs(rr - np.median(rr[ok])) < 3 * mad)
            x = dd[mm] - dd[mm].mean()
            y = rr[mm] - rr[mm].mean()
            s = float((x * y).sum() / (x * x).sum())
            s_err = float(np.sqrt(((y - s * x) ** 2).sum() / (len(x) - 2) / (x * x).sum()))
    return (oX, oY, s, s_err) if with_prop else (oX, oY, s)


rows, hrows, offs, summary = [], [], [], {"events_good_time": int(N), "tagged": int(sel.sum()), "mixed": int(selm.sum()), "scale": scale}
for method in ("cfd50", "cfd50_bl", "kmin"):
    tX, tY = times(method, cx), times(method, cy)
    tYm = times(method, cy[nxt])[: N]  # Y del evento siguiente (tiempo del evento i+1)
    dt_all = tX - tY
    # desplazamientos con los etiquetados
    S = np.nonzero(sel & np.isfinite(dt_all))[0]
    cor = {"none": (np.zeros(32), np.zeros(32), 0.0)}
    cor["offsets"] = solve_offsets(dt_all[S], cx[S], cy[S], (dX - dY)[S], False)
    oXp, oYp, sp, sp_err = solve_offsets(dt_all[S], cx[S], cy[S], (dX - dY)[S], True)
    cor["offsets+prop"] = (oXp, oYp, sp)
    summary[f"{method}_prop_slope_ns_per_mm"] = sp
    summary[f"{method}_prop_slope_err_ns_per_mm"] = sp_err
    summary[f"{method}_dXmdY_range_mm"] = [float(np.min((dX - dY)[S])), float(np.max((dX - dY)[S])), float(np.std((dX - dY)[S]))]
    summary[f"{method}_nan_fraction_tagged"] = float(np.mean(~np.isfinite(dt_all[sel])))
    for c in range(32):
        offs.append([method, c, "X" if c in posX else "Y", cor["offsets"][0][c] if c in posX else cor["offsets"][1][c],
                     cor["offsets+prop"][0][c] if c in posX else cor["offsets+prop"][1][c]])
    for cname, (oX, oY, s) in cor.items():
        d_same = dt_all - oX[cx] + oY[cy] - s * (dX - dY)
        # mezclados: X del evento i, Y del evento i+1 (canal, tiempo y posición del i+1)
        cyn = cy[nxt]
        dXm = readout_d("X", jY[nxt], cx)
        dYm = readout_d("Y", jX, cyn)
        d_mix = tX - tYm - oX[cx] + oY[cyn] - s * (dXm - dYm)
        mn_s, mn_m = np.minimum(ax, ay), np.minimum(ax, ay[nxt])
        for bname, lo, hi in BINS:
            ms = sel & (mn_s > lo) & (mn_s <= hi) & np.isfinite(d_same)
            mm = selm & (mn_m > lo) & (mn_m <= hi) & np.isfinite(d_mix)
            bw = DS if method == "kmin" else 1.0
            r, edges, hs, hm = fit_hist(d_same[ms], d_mix[mm], bw, 80.0)
            if r is None:
                continue
            rows.append([method, cname, bname, r["n_same"], round(r["n_mix"], 2), round(r["n_net"], 2), r["mu"], r["mu_err"], r["sigma"],
                         r["sigma_err"], r["rms3"], r["sigma"] / np.sqrt(2), r["sigma_err"] / np.sqrt(2), r["rms3"] / np.sqrt(2), r["chi2"], r["ndf"]])
            if method == "kmin":  # fracciones netas en Δk = 0, ±1, ±2 (bins de 16 ns centrados en múltiplos de 16 ns)
                xc = 0.5 * (edges[1:] + edges[:-1])
                net = hs - hm
                summary[f"kmin_{cname}_{bname}_frac_dk"] = {str(int(round(v / DS))): float(net[np.abs(xc - v) < 1].sum() / net.sum())
                                                            for v in (-2 * DS, -DS, 0.0, DS, 2 * DS)}
            if cname in ("none", "offsets+prop"):
                for k in range(len(hs)):
                    hrows.append([method, cname, bname, edges[k], edges[k + 1], int(hs[k]), hm[k], hs[k] - hm[k]])
            print(f"{method:9s} {cname:13s} {bname:6s} n={r['n_same']:5d} mix={r['n_mix']:7.1f} mu={r['mu']:7.3f} "
                  f"sigma={r['sigma']:6.3f}±{r['sigma_err']:.3f} rms3={r['rms3']:6.3f} strip={r['sigma'] / np.sqrt(2):6.3f} chi2/ndf={r['chi2']:.1f}/{r['ndf']}",
                  flush=True)
with open(f"{OUT}/time_res.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["# time_res.py: Δt = t_X − t_Y de muones etiquetados (X y Y > 10 p.e.), mayo v8 tiempo bueno; mismo − mezclados; ns"])
    w.writerow(["method", "correction", "min_amp_bin_pe", "n_same", "n_mix_scaled", "n_net", "mu_ns", "mu_err_ns", "sigma_ns", "sigma_err_ns",
                "rms_3sigma_ns", "sigma_strip_ns", "sigma_strip_err_ns", "rms_strip_ns", "chi2", "ndf"])
    w.writerows(rows)
with open(f"{OUT}/time_res_hist.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["method", "correction", "min_amp_bin_pe", "lo_ns", "hi_ns", "n_same", "n_mix_scaled", "n_net"])
    w.writerows(hrows)
with open(f"{OUT}/time_res_offsets.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["method", "ch", "coord", "offset_ns", "offset_with_prop_ns"])
    w.writerows(offs)
json.dump(summary, open(f"{OUT}/time_res.json", "w"), indent=1)
print(json.dumps(summary))
