#!/usr/bin/python3 -I
# overlap_ratio.py — seguimiento 2026-11-03 (13), punto B: posición dentro del solapamiento de 1 mm entre tiras vecinas de hod2018 (run B, mayo v8,
# tiempo bueno) y comparación con el MC del banco (overlap_mc.py).
# Datos (<work>/may_v8_time.bin de time_extract.C; a = aS/G en p.e. medidos, G = reajuste de PROC, ch3 muerto fuera):
#   etiqueta = alguna tira de la coordenada ortogonal > 10 p.e. (la más alta da la posición a lo largo de la barra);
#   sonda = exactamente dos tiras de la coordenada > 5 p.e., vecinas en el orden de posición (stage3.json, mayo), ninguna es ch3;
#   R = a1/(a1 + a2), 1 = la tira de menor posición. Mezclados: sonda del evento i con etiqueta del evento i+1, × N/(N−1); neto = mismo − mezclados.
#   "within_pair": R − mediana del par + 0.5 (sin restar mezclados), para separar la forma de las diferencias de respuesta entre tiras.
#   Variante "att": a·exp(d/L_eff) con d = distancia del cruce al extremo de lectura de cada tira (las dos capas se leen por extremos opuestos;
#   posición a lo largo = tira etiqueta más alta, como time_res.py) y L_eff de la mediana por capa del seguimiento 11 (may_layer_targets.csv:
#   XA 114.5, XB 79.8, YA 79.1, YB 97.5 mm).
# MC (overlap_mc.npz: rayos del trigger con etiqueta y exactamente dos barras vecinas cruzadas; los pares que contienen la posición de ch3 se
# quitan de X): R_path = L1/(L1 + L2); u_mid = posición del cruce respecto del centro del solapamiento.
#   "pois": a_i = n_i + NegBin(n_i, 1/F) (diafonía geométrica con multiplicación media F = 1.75), n_i ~ Poisson(μ L_i cosθ / 1 mm), μ tal que la
#           media medida de un cruce completo es M = 13 p.e. (fotoestadística sola, el límite pedido); umbral 5 p.e. en las dos.
#   "spec": a_i = s_i · L_i cosθ / 1 mm, s_i independientes del espectro medido de tiras únicas etiquetadas (una sola tira de la coordenada > 5 p.e.,
#           con etiqueta, sin restar mezclados) de la misma coordenada: incluye Landau, fotoestadística y la dispersión de ganancias.
#   Resolución en posición: pendiente dR/du del MC (regresión de R_path sobre u_mid) y σ_u = σ_R / |dR/du| con σ_R de "pois", frente al RMS de
#   u_mid de los pares (lo que da saber que el cruce está en el solapamiento).
# Salida: <out>/overlap_ratio.csv (histogramas de R), overlap_pairs.csv (por par de tiras), overlap_ratio.json. Uso: overlap_ratio.py <work> <mc dir> <out>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

WK, MCD, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("kS", "<u2", 32), ("bl", "<f4", 32), ("pl", "<f4", 32), ("w", "<f4", (32, 10))])
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ORD = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
LEFF = {"XA": 114.5, "XB": 79.8, "YA": 79.1, "YB": 97.5}
DEAD, TAG, THR, F, M = 3, 10.0, 5.0, 1.75, 13.0
EDG = np.linspace(0, 1, 41)


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


def stats(R, w=None):
    w = np.ones(len(R)) if w is None else w
    if len(R) == 0 or w.sum() <= 0:
        return {}
    o = np.argsort(R)
    cs = np.cumsum(w[o]) / w.sum()
    q = lambda p: float(R[o][np.searchsorted(cs, p)])
    m = float((w * R).sum() / w.sum())
    return {"n": float(w.sum()), "mean": m, "rms": float(np.sqrt((w * (R - m) ** 2).sum() / w.sum())), "median": q(0.5),
            "iqr_sigma": (q(0.75) - q(0.25)) / 1.349, "f_0.4_0.6": float(w[(R > 0.4) & (R < 0.6)].sum() / w.sum())}


def hstats(h):
    """Estadística de un histograma neto en EDG (centros), negativos a 0."""
    xc = 0.5 * (EDG[1:] + EDG[:-1])
    return stats(xc, np.clip(h, 0, None))


G = g_refit()
e = np.fromfile(f"{WK}/may_v8_time.bin", dtype=DT)
e = e[mc.good_time(e["t"], "v8")]
N = len(e)
pe = e["aS"] / G
pe[:, DEAD] = 0
scale = N / (N - 1)
nxt = np.roll(np.arange(N), -1)
res, hrows, prow = {"events_good_time": int(N), "F": F, "M_pe": M, "L_eff_mm": LEFF}, [], []
mcz = np.load(f"{MCD}/overlap_mc.npz")
rng = np.random.default_rng(7)
for co, ref in (("X", "Y"), ("Y", "X")):
    po, ro = np.asarray(ORD[f"h18_{co}"]), np.asarray(ORD[f"h18_{ref}"])
    a = pe[:, po]  # columnas = posición
    ar = pe[:, ro]
    tag = ar.max(1) > TAG
    jt = np.argmax(ar, 1)  # posición de la etiqueta a lo largo de la barra
    over = a > THR
    two = over.sum(1) == 2
    j1 = np.argmax(over, 1)
    adj = two & (j1 < 15) & over[np.arange(N), np.minimum(j1 + 1, 15)]
    bad = np.isin(po[j1], [DEAD]) | np.isin(po[np.minimum(j1 + 1, 15)], [DEAD])
    adj &= ~bad
    a1, a2 = a[np.arange(N), j1], a[np.arange(N), np.minimum(j1 + 1, 15)]
    R = a1 / (a1 + a2)
    # atenuación: distancia de cada tira a su extremo de lectura
    D = 22.5 + 2.0 * jt
    lay1 = np.where(j1 % 2 == 0, "A", "B")
    lay2 = np.where((j1 + 1) % 2 == 0, "A", "B")
    s1 = np.array([SIGN[f"may2026_v8_goodtime/{co}{L}"]["sign"] for L in lay1])
    s2 = np.array([SIGN[f"may2026_v8_goodtime/{co}{L}"]["sign"] for L in lay2])
    d1, d2 = np.where(s1 < 0, D, 75 - D), np.where(s2 < 0, D, 75 - D)
    l1 = np.array([LEFF[f"{co}{L}"] for L in lay1])
    l2 = np.array([LEFF[f"{co}{L}"] for L in lay2])
    c1, c2 = a1 * np.exp(d1 / l1), a2 * np.exp(d2 / l2)
    Ra = c1 / (c1 + c2)
    same = adj & tag
    mix = adj & tag[nxt]
    mix[-1] = False
    # mezclados para la variante att: la posición a lo largo viene de la etiqueta del evento siguiente
    Dm = 22.5 + 2.0 * jt[nxt]
    d1m, d2m = np.where(s1 < 0, Dm, 75 - Dm), np.where(s2 < 0, Dm, 75 - Dm)
    c1m, c2m = a1 * np.exp(d1m / l1), a2 * np.exp(d2m / l2)
    Ram = c1m / (c1m + c2m)
    hs, hm = np.histogram(R[same], EDG)[0], np.histogram(R[mix], EDG)[0] * scale
    hsa, hma = np.histogram(Ra[same], EDG)[0], np.histogram(Ram[mix], EDG)[0] * scale
    r = res[co] = {"pairs_same": int(same.sum()), "pairs_mixed_scaled": float(mix.sum() * scale),
                   "data": hstats(hs - hm), "data_att": hstats(hsa - hma), "data_same_unsubtracted": stats(R[same])}
    # por par de tiras (mismo evento, sin restar)
    for j in range(15):
        if DEAD in (po[j], po[j + 1]):
            continue
        m = same & (j1 == j)
        st = stats(R[m]) if m.sum() else {}
        prow.append([co, j, int(po[j]), int(po[j + 1]), int(m.sum()), st.get("median", ""), st.get("mean", ""), st.get("rms", ""),
                     float(np.median(a1[m])) if m.sum() else "", float(np.median(a2[m])) if m.sum() else ""])
    # ancho dentro de cada par: R − mediana de su par + 0.5 (mismo evento, sin restar mezclados), quita las diferencias de respuesta entre tiras
    Rw = np.full(N, np.nan)
    for j in range(15):
        m = same & (j1 == j)
        if m.sum():
            Rw[m] = R[m] - np.median(R[m]) + 0.5
    r["data_within_pair_unsubtracted"] = stats(Rw[same & np.isfinite(Rw)])
    # espectro de tiras únicas etiquetadas (para "spec")
    one = tag & (over.sum(1) == 1)
    single = a[one][over[one]]
    r["single_strip_spectrum"] = {"n": int(len(single)), "median_pe": float(np.median(single)), "mean_pe": float(np.mean(single))}
    r["pair_amplitudes_median_pe"] = [float(np.median(a1[same])), float(np.median(a2[same]))]
    for k in range(len(EDG) - 1):
        hrows.append([co, "data", "", EDG[k], EDG[k + 1], int(hs[k]), hm[k], hs[k] - hm[k]])
        hrows.append([co, "data_att", "", EDG[k], EDG[k + 1], int(hsa[k]), hma[k], hsa[k] - hma[k]])
    # MC
    for top in "XY":
        L1, L2, u, cos, jj = (mcz[f"{top}|{co}|{k}"] for k in ("L1", "L2", "u_mid", "cos", "j"))
        keep = ~np.isin(po[jj], [DEAD]) & ~np.isin(po[np.minimum(jj + 1, 15)], [DEAD])
        L1, L2, u, cos, jj = L1[keep], L2[keep], u[keep], cos[keep], jj[keep]
        Rp = L1 / (L1 + L2)
        mu = M / F
        n1, n2 = rng.poisson(mu * L1 * cos), rng.poisson(mu * L2 * cos)
        p1 = n1 + np.where(n1 > 0, rng.negative_binomial(np.maximum(n1, 1), 1 / F), 0)
        p2 = n2 + np.where(n2 > 0, rng.negative_binomial(np.maximum(n2, 1), 1 / F), 0)
        okp = (p1 > THR) & (p2 > THR)
        Rpo = p1[okp] / (p1[okp] + p2[okp])
        s1_, s2_ = rng.choice(single, len(L1)), rng.choice(single, len(L1))
        q1, q2 = s1_ * L1 * cos, s2_ * L2 * cos
        oks = (q1 > THR) & (q2 > THR)
        Rsp = q1[oks] / (q1[oks] + q2[oks])
        x = u - u.mean()
        slope = float((x * (Rp - Rp.mean())).sum() / (x * x).sum())
        # σ_R de la fotoestadística sola en R_path = 0.5 (pares verticales): dispersión de Rpo − Rp de los mismos rayos
        sigR = float(np.std(Rpo - Rp[okp]))
        rm = res[co][f"mc_top{top}"] = {
            "pairs": int(len(L1)), "R_path": stats(Rp), "R_pois": stats(Rpo), "R_spec": stats(Rsp),
            "eff_pois_both_gt5": float(okp.mean()), "eff_spec_both_gt5": float(oks.mean()),
            "u_mid": {"mean": float(u.mean()), "rms": float(u.std()), "p05": float(np.percentile(u, 5)), "p95": float(np.percentile(u, 95)),
                      "f_inside_overlap": float(np.mean(np.abs(u) <= 0.5))},
            "f_R_path_eq_0.5": float(np.mean(np.abs(Rp - 0.5) < 1e-6)), "f_L1_eq_L2_within_1pct": float(np.mean(np.abs(L1 - L2) < 0.01 * (L1 + L2))),
            "dRpath_du_per_mm": slope, "corr_Rpath_u": float(np.corrcoef(Rp, u)[0, 1]),
            "sigma_R_photostat": sigR, "sigma_u_from_ratio_mm": float(sigR / abs(slope)) if slope else float("inf")}
        for name, RR in (("mc_path", Rp), ("mc_pois", Rpo), ("mc_spec", Rsp)):
            h = np.histogram(RR, EDG)[0]
            for k in range(len(EDG) - 1):
                hrows.append([co, name, top, EDG[k], EDG[k + 1], int(h[k]), "", ""])
        hu = np.histogram(u, np.linspace(-3, 3, 61))[0]
        for k in range(60):
            hrows.append([co, "mc_u_mid", top, -3 + 0.1 * k, -3 + 0.1 * (k + 1), int(hu[k]), "", ""])
    print(co, json.dumps({k: v for k, v in r.items()}, default=float)[:2500], flush=True)
with open(f"{OUT}/overlap_ratio.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["# overlap_ratio.py: R = a1/(a1+a2) de pares de tiras vecinas > 5 p.e. con etiqueta ortogonal > 10 p.e.; MC: R_path, pois, spec; mc_u_mid en mm"])
    w.writerow(["coord", "series", "mc_top", "lo", "hi", "n_same", "n_mix_scaled", "n_net"])
    w.writerows(hrows)
with open(f"{OUT}/overlap_pairs.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["coord", "pos_j", "ch_j", "ch_j1", "n_same", "R_median", "R_mean", "R_rms", "a1_median_pe", "a2_median_pe"])
    w.writerows(prow)
json.dump(res, open(f"{OUT}/overlap_ratio.json", "w"), indent=1)
