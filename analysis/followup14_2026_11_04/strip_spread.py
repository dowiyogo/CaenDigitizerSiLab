#!/usr/bin/python3 -I
# strip_spread.py — 2026-11-04: dispersión de luz entre tiras de hod2018, por capa y campaña (mayo v8 tiempo bueno, octubre run 8), intervalo 29-54 mm.
# Por tira: mediana etiquetada como en layer_targets.py / bar_targets.py (etiqueta = alguna tira ortogonal del intervalo > 10 p.e., ch3 fuera; sonda > 5 p.e. medidos;
# mismo − mezclados, bins de 1 p.e., netos negativos a 0, mediana interpolada). Se recalcula porque bar_targets.csv solo tiene los cuatro intervalos, no 29-54.
# Bootstrap de eventos (200, los eventos con etiqueta en el mismo evento o en el siguiente): con los mismos pesos para todas las tiras de la capa; e_i = SD de la mediana de la tira i.
# S = sd(ddof = 1)/media de las medianas de las tiras vivas.
#   err_stat = SD de S sobre las réplicas (error de medida con estas tiras).
#   err_sampling = S/√(2(n − 1)) (si las n tiras se ven como muestra de una población gaussiana de tiras).
#   S_stat = √(media e_i²)/media: lo esperado solo por estadística si todas las tiras tuvieran la misma mediana. Toy: m_i ~ N(media, e_i), 20000 veces, da la media, el
#   intervalo 16-84 % de S y p = P(S_toy ≥ S). S_intr = √(S² − S_stat²) si S > S_stat. χ² = Σ (m_i − m_pond)²/e_i², ndf = n − 1.
# Variantes: without_filled_G = en octubre sin ch22 y ch31 (sin G en c1, se usa la mediana de la estación); without_Y13 = YB sin la barra Y13 (mayo ch30, octubre ch29, el mismo
# MPPC, con la mediana más alta de la capa en las dos campañas).
# Salida: <out>/strip_spread.csv, strip_medians.csv. Uso: strip_spread.py <out>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD, DEAD_X, TAG, THR, NBOOT, LO, HI = 3, 12, 10.0, 5.0, 200, 29, 54
D = 37.5 - 15.0 + 2.0 * np.arange(16)
EDG = np.arange(0, 301, 1.0)
RNG = np.random.default_rng(1104)
FILLED = {"may": {}, "oct": {22, 31}}


def med(net):
    net = np.clip(net, 0, None)
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    return float(EDG[b] + (0.5 * n - (cs[b - 1] if b else 0.0)) / net[b])


def spread(m):
    return float(np.std(m, ddof=1) / np.mean(m))


def load(camp):
    if camp == "may":
        g = np.full(32, np.nan)
        for l in open(PK + "/hd64_may_v8_peaks.csv"):
            f = l.strip().split(",")
            if f[0] != "ch" and f[3] in ("ok", "changed"):
                g[int(f[0])] = float(f[2])
        e = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
        e = e[np.argsort(e["t"], kind="stable")]
        e = e[mc.good_time(e["t"], "v8")]
    else:
        g = mc.g_c1()
        e = xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)
        e = e[np.argsort(e["t"], kind="stable")]
    return e, np.where(np.isfinite(g), g, np.nanmedian(g))


srows, mrows = [], []
for camp in ("may", "oct"):
    e, g = load(camp)
    o, skey = ST3[camp]["order"], {"may": "may2026_v8_goodtime", "oct": "oct2026_run8"}[camp]
    N = len(e)
    s = N / (N - 1)
    for probe, ref in (("X", "Y"), ("Y", "X")):
        refch = np.asarray(o[f"h18_{ref}"])
        tagS = e["aS"][:, refch] > TAG * g[refch]
        tagS[:, refch == DEAD] = False
        for L in "AB":
            sgn = SIGN[f"{skey}/{probe}{L}"]["sign"]
            dist = D if sgn < 0 else 75.0 - D
            js = [j for j in range(16) if LO <= dist[j] < HI + 1e-9 and not (probe == "Y" and j == DEAD_X)]
            T = tagS[:, js].any(1)
            Tn = np.concatenate([T[1:], [False]])
            U = np.nonzero(T | Tn)[0]
            idx = [(i, c) for i, c in enumerate(o[f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
            nom, reps, ninfo = [], np.zeros((NBOOT, len(idx))), []
            W = RNG.multinomial(len(U), np.full(len(U), 1.0 / len(U)), size=NBOOT).astype(float)
            for q, (i, c) in enumerate(idx):
                a = e["aS"][U, c] / g[c]
                m = (a > THR) & (a < EDG[-1])
                ev, bn, iS, iM = np.nonzero(m)[0], a[m].astype(int), T[U][m], Tn[U][m]
                hs, hm = np.bincount(bn[iS], minlength=300).astype(float), np.bincount(bn[iM], minlength=300).astype(float)
                nom.append(med(hs - s * hm))
                ninfo.append((hs.sum(), s * hm.sum()))
                for b in range(NBOOT):
                    w = W[b][ev]
                    reps[b, q] = med(np.bincount(bn[iS], w[iS], 300) - s * np.bincount(bn[iM], w[iM], 300))
            nom = np.array(nom)
            ei = np.nanstd(reps, axis=0, ddof=1)
            for q, (i, c) in enumerate(idx):
                mrows.append([camp, probe + L, i, int(c), int(ninfo[q][0]), round(ninfo[q][1], 3), round(ninfo[q][0] - ninfo[q][1], 3), nom[q], ei[q], int(c in FILLED[camp])])
            for var in ("all", "without_filled_G", "without_Y13"):
                keep = {"all": np.ones(len(idx), bool), "without_filled_G": np.array([c not in FILLED[camp] for _, c in idx]),
                        "without_Y13": np.array([not (probe == "Y" and i == 13) for i, _ in idx])}[var]
                if var != "all" and keep.all():
                    continue
                m, er, rp = nom[keep], ei[keep], reps[:, keep]
                n = len(m)
                Sv, mean = spread(m), float(m.mean())
                Sb = np.array([spread(r) for r in rp if np.isfinite(r).all()])
                Sst = float(np.sqrt((er ** 2).mean()) / mean)
                toy = RNG.normal(mean, er, (20000, n))
                St = toy.std(axis=1, ddof=1) / toy.mean(axis=1)
                wm = float((m / er ** 2).sum() / (1 / er ** 2).sum())
                chi = float((((m - wm) / er) ** 2).sum())
                intr = float(np.sqrt(Sv ** 2 - Sst ** 2)) if Sv > Sst else 0.0
                srows.append([camp, probe + L, var, n, mean, Sv, float(Sb.std(ddof=1)), Sv / np.sqrt(2 * (n - 1)), Sst, float(St.mean()), float(np.percentile(St, 16)), float(np.percentile(St, 84)),
                              float((St >= Sv).mean()), intr, chi, n - 1, float(m.min()), float(m.max()), int(min(x[0] - x[1] for x, k in zip(ninfo, keep) if k)), len(Sb)])
                print(f"{camp} {probe}{L} {var}: n {n} media {mean:.2f} S {Sv:.3f} ± {Sb.std(ddof=1):.3f} (muestreo {Sv / np.sqrt(2 * (n - 1)):.3f}) S_stat {Sst:.3f} toy {St.mean():.3f} "
                      f"[{np.percentile(St, 16):.3f}, {np.percentile(St, 84):.3f}] p {(St >= Sv).mean():.4f} S_intr {intr:.3f} χ² {chi:.1f}/{n - 1} rango {m.min():.1f}-{m.max():.1f}", flush=True)
fm = lambda v: f"{v:.5g}" if isinstance(v, float) else str(v)
with open(OUT + "/strip_spread.csv", "w") as f:
    f.write("# strip_spread.py: S = sd/mean of the per-strip tagged median amplitude (measured p.e.) in the 29-54 mm bin; err_stat = SD over 200 event-bootstrap resamples; err_sampling = S/sqrt(2(n-1)); "
            "S_stat = sqrt(mean e_i^2)/mean, expected from statistics alone; toy: m_i ~ N(mean, e_i); p = P(S_toy >= S); S_intrinsic = sqrt(S^2 - S_stat^2)\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "layer", "strips", "n_strips", "mean_median_pe", "S", "err_stat", "err_sampling", "S_stat", "S_toy_mean", "S_toy_p16", "S_toy_p84", "p_toy", "S_intrinsic", "chi2", "ndf", "min_median_pe", "max_median_pe", "min_n_net", "n_boot"])
    w.writerows([[fm(x) for x in r] for r in srows])
with open(OUT + "/strip_medians.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "layer", "strip_index", "channel", "n_same", "n_mix_scaled", "n_net", "median_pe", "median_pe_err_boot", "G_filled_with_station_median"])
    w.writerows([[fm(x) for x in r] for r in mrows])
