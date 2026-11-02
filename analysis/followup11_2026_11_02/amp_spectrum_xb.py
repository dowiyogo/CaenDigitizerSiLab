#!/usr/bin/python3 -I
# amp_spectrum_xb.py — 2026-11-02: espectro de amplitud en p.e. de un canal de la capa XB de hod2018 con G cerca de la mediana de la estación,
# para el panel (a) de mip_spectra_paper_2p (mip_paper_2p.C). Mismo cálculo que daq/tools/proceeding/stage4.py (figura
# amp_spectra_hod2018_pe): mayo v8 tiempo bueno, aS/G (ventana signal) y aC/G (ventana fuera de tiempo, mismo ancho), G reajustada de PROC,
# 0-30 p.e. en bins de 0.5. Canal: el de XB con |G − mediana| mínimo (mediana sobre los canales vivos con G ajustada, como stage4.py).
# Salida: <out>/fig_amp_xb.json (formato de fig_stage4.json). Uso: amp_spectrum_xb.py <out dir>
import json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
LIVE = [c for c in range(32) if c != 3]
G0 = np.full(32, np.nan)
for l in open(W + "/hd64_may_v8_peaks.csv"):
    f = l.strip().split(",")
    if f[0] != "ch" and f[3] in ("ok", "changed"):
        G0[int(f[0])] = float(f[2])
G = np.where(np.isfinite(G0), G0, np.nanmedian(G0))
fit_ok = [c for c in LIVE if np.isfinite(G0[c])]
g_med = float(np.median(G0[fit_ok]))
order = json.load(open(W + "/stage3.json"))["strip_orders"]["may"]["order"]["h18_X"]
xb = [c for i, c in enumerate(order) if i % 2 == 1 and c in fit_ok]
ch = min(xb, key=lambda c: (abs(G0[c] - g_med), c))
print("mediana de G (canales vivos con ajuste):", g_med, "| XB:", {c: round(float(G0[c]), 5) for c in xb}, "| elegido ch", ch, "G", G0[ch])
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
eg = e8[mc.good_time(e8["t"], "v8")]
edges = np.arange(0, 30.0001, 0.5)
hi, _ = np.histogram(eg["aS"][:, ch] / G[ch], edges)
ho, _ = np.histogram(eg["aC"][:, ch] / G[ch], edges)
print("eventos", len(eg), "| en tiempo en [0,30):", int(hi.sum()), "≥ 30:", int((eg["aS"][:, ch] / G[ch] >= 30).sum()), "| fuera de tiempo:", int(ho.sum()))
os.makedirs(OUT, exist_ok=True)
json.dump({"edges": edges.tolist(), "G_median_mV": g_med, "events": int(len(eg)),
           "channels": [{"rule": "XB, gain nearest the median", "ch": int(ch), "G_mV": float(G[ch]), "in_time": hi.tolist(), "out_of_time": ho.tolist()}]},
          open(OUT + "/fig_amp_xb.json", "w"), indent=1)
