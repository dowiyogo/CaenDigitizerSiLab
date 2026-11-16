#!/usr/bin/python3 -I
# r2_profile.py — seguimiento 2026-11-16, punto 2: perfil de amplitud a lo largo de la barra de los canales Y06-Y09 en cada campaña.
# Las capas A y B de Y se leen por extremos opuestos (layer_from_channels.json); un canal en la barra equivocada tendría la pendiente del perfil con el
# signo de la otra capa. Por canal y campaña: mediana interpolada (selección de reconcile.py: etiqueta = tira X > 10 p.e., sonda > 3 p.e. medidos,
# mismo − mezclados, bins de 0.05 p.e.) en 4 grupos de 4 tiras X por índice (0-3, 4-7, 8-11, 12-15; X12 muerta fuera), error por remuestreo
# (200); ajuste ponderado de ln(mediana) contra el índice medio del grupo: pendiente por índice de tira X (2 mm) [C].
# Salida: <out>/r2_profile.csv. Uso: r2_profile.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
EDG = np.arange(0, 300.0001, 0.05)


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))  # noqa: E731
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1()[:32])}
e8 = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
DATA = {"may": e8[mc.good_time(e8["t"], "v8")], "oct": xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)}
PE = {c: DATA[c]["aS"][:, :32] / G[c][None, :32] for c in DATA}
rng = np.random.default_rng(1602)


def med(net):
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    before = cs[b - 1] if b > 0 else 0.0
    return float(EDG[b] + (0.5 * n - before) / net[b] * 0.05)


rows = []
for camp in ("may", "oct"):
    o = ORDER[camp]
    xc = np.asarray(o["h18_X"])
    t = PE[camp][:, xc] > 10.0
    t[:, xc == 3] = False
    s = len(PE[camp]) / (len(PE[camp]) - 1)
    for j in (6, 7, 8, 9):
        ch = o["h18_Y"][j]
        layer = "A" if j % 2 == 0 else "B"
        pts = []
        for k in range(4):
            js = [x for x in range(4 * k, 4 * k + 4) if x != 12]
            T = t[:, js].any(1)
            Tn = np.concatenate([T[1:], [False]])
            ps = PE[camp][:, ch]
            sel = ps > 3.0
            hs = np.histogram(ps[sel & T], bins=EDG)[0].astype(float)
            hm = np.histogram(ps[sel & Tn], bins=EDG)[0].astype(float)
            m = med(np.clip(hs - s * hm, 0, None))
            bs = [med(np.clip(rng.poisson(hs) - s * rng.poisson(hm), 0, None)) for _ in range(200)]
            pts.append((float(np.mean(js)), m, float(np.nanstd(bs))))
        x = np.array([p[0] for p in pts])
        y = np.log([p[1] for p in pts])
        w = (np.array([p[1] for p in pts]) / np.array([p[2] for p in pts])) ** 2
        X = np.vstack([np.ones_like(x), x - 7.5]).T
        C = np.linalg.inv((X * w[:, None]).T @ X)
        b = C @ (X * w[:, None]).T @ y
        rows.append([camp, f"Y{j:02d}", layer, SIGN[f"oct2026_run8/Y{layer}"]["sign"], ch] + [f"{p[1]:.3f}±{p[2]:.3f}" for p in pts] +
                    [f"{b[1]:.4f}", f"{math.sqrt(C[1, 1]):.4f}", f"{b[1] / math.sqrt(C[1, 1]):.2f}"])
with open(OUT + "/r2_profile.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "bar", "layer", "layer_readout_sign", "channel", "median_X0-3", "median_X4-7", "median_X8-11", "median_X12-15", "slope_ln_per_X_index", "err",
                "z"])
    w.writerows(rows)
for r in rows:
    print(" | ".join(map(str, r)))
