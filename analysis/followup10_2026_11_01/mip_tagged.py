#!/usr/bin/python3 -I
# mip_tagged.py — 2026-11-01: espectros MIP etiquetados de hod2018 por capa (XA, XB, YA, YB), campaña (mayo v8 tiempo bueno, octubre run 8
# completo) e intervalo de distancia al extremo de lectura (cerca 21-29 mm, lejos 29-54 mm), en amplitud y en carga, en p.e. primarios.
# Entradas: <work>/{may_v8,oct_run8}_charge.bin y _dark.csv de charge_extract.C (aS = amplitud del pulso más alto de la ventana signal,
# QS = su carga en [k−3, k+12] con pedestal local; Q1 = mediana de la carga de pulsos oscuros aislados de 1 p.e. del canal).
# Selección: etiqueta = alguna tira de la coordenada ortogonal del intervalo con aS > 10 G (ch3 muerto, fuera); sonda = cada tira viva de la
# capa con aS > 3 G (amplitud) o QS > 3 Q1 (carga); mezclados: etiqueta del evento siguiente, × N/(N − 1); neto = mismo − mezclados,
# error² = mismo + escala² · mezclados. Todas las tiras vivas de la capa sumadas. p.e. primarios = (aS/G o QS/Q1) / F, F = 1.75 mayo,
# 1.46 octubre; bins de 1 p.e. primario en [0, 150).
# G: mayo = reajuste de PROC, octubre = c1 (mediana para los canales sin G). Orden de tiras y extremo de lectura: los de bar_targets.py.
# Mediana de la selección del seguimiento 6 por capa (sin etiqueta, todas las posiciones): por canal, bins de 1 p.e. medido, > 5 p.e. medidos,
# signal − fuera de tiempo, negativos a 0, centros (n + 0.5)/F con el F de campaña; mediana ponderada de la capa.
# Salida: <out>/spectra.csv, <out>/mip_tagged.json. Uso: mip_tagged.py <out dir>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

OUT = sys.argv[1]
WK = OUT + "/work"
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD = 3
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = {"near": (21, 29), "far": (29, 54)}
TAG, PROBE_T = 10.0, 3.0
F = {"may": 1.75, "oct": 1.46}
EDG = np.arange(0, 151, 1.0)
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("aC", "<f4", 32), ("QS", "<f4", 32), ("QC", "<f4", 32)])


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1())[:32]}
Q1 = {}
for camp, pre in (("may", "may_v8"), ("oct", "oct_run8")):
    q = np.full(32, np.nan)
    for r in csv.DictReader(open(f"{WK}/{pre}_dark.csv")):
        if r["Q1_median_mVns"] and int(r["ch"]) != DEAD:
            q[int(r["ch"])] = float(r["Q1_median_mVns"])
    Q1[camp] = q
DATA = {}
for camp, pre in (("may", "may_v8"), ("oct", "oct_run8")):
    e = np.fromfile(f"{WK}/{pre}_charge.bin", dtype=DT)
    e = e[np.argsort(e["t"], kind="stable")]
    if camp == "may":
        e = e[mc.good_time(e["t"], "v8")]
    DATA[camp] = e


def readout_sign(probe, layer):
    s = SIGN[f"oct2026_run8/{probe}{layer}"]["sign"]
    assert s == SIGN[f"may2026_v8_goodtime/{probe}{layer}"]["sign"]
    return s


def bin_js(probe, layer, b):
    lo, hi = BINS[b]
    d = D if readout_sign(probe, layer) < 0 else 75.0 - D
    js = [j for j in range(16) if lo <= d[j] < hi + (1e-9 if hi == 54 else 0)]
    return js


def f6_median(e, chans, g, Fc):
    vals, ws = [], []
    for c in chans:
        s, o = e["aS"][:, c] / g[c], e["aC"][:, c] / g[c]
        net = np.clip(np.histogram(s[s > 5], bins=np.arange(0, 301, 1.0))[0] - np.histogram(o[o > 5], bins=np.arange(0, 301, 1.0))[0], 0, None)
        vals.append((np.arange(300) + 0.5) / Fc); ws.append(net)
    v, w = np.concatenate(vals), np.concatenate(ws).astype(float)
    o_ = np.argsort(v); cs = np.cumsum(w[o_])
    return float(v[o_][int(np.argmax(cs >= 0.5 * cs[-1]))])


rows, res = [], {"F": F, "tag_pe": TAG, "probe_threshold_pe": PROBE_T, "events": {k: int(len(v)) for k, v in DATA.items()},
                 "Q1_median_mVns": {k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in Q1.items()}, "layers": {}}
for camp in ("may", "oct"):
    e, o, g, q1, Fc = DATA[camp], ORDER[camp], G[camp], Q1[camp], F[camp]
    N = len(e)
    s = N / (N - 1)
    for probe, ref in (("X", "Y"), ("Y", "X")):
        refch = np.asarray(o[f"h18_{ref}"])
        tagS = e["aS"][:, refch] > TAG * g[refch]
        tagS[:, refch == DEAD] = False
        for L in ("A", "B"):
            chans = [c for i, c in enumerate(o[f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
            key = f"{camp}_{probe}{L}"
            res["layers"][key] = {"channels": chans, "f6_median_pe_prim": f6_median(e, chans, g, Fc)}
            for b in BINS:
                js = bin_js(probe, L, b)
                T = tagS[:, js].any(1)
                Tn = np.concatenate([T[1:], [False]])
                res["layers"][key][f"{b}_tag_strips"] = js
                res["layers"][key][f"{b}_tag_events"] = int(T.sum())
                for var in ("amp", "charge"):
                    hs, hm = np.zeros(len(EDG) - 1), np.zeros(len(EDG) - 1)
                    for c in chans:
                        meas = e["aS"][:, c] / g[c] if var == "amp" else e["QS"][:, c] / q1[c]
                        sel = meas > PROBE_T
                        hs += np.histogram(meas[sel & T] / Fc, bins=EDG)[0]
                        hm += np.histogram(meas[sel & Tn] / Fc, bins=EDG)[0]
                    net = hs - s * hm
                    err = np.sqrt(hs + s * s * hm)
                    res["layers"][key][f"{b}_{var}_n_net"] = float(net.sum())
                    for k in range(len(net)):
                        rows.append([camp, f"{probe}{L}", b, var, EDG[k], EDG[k + 1], int(hs[k]), s * hm[k], net[k], err[k]])
            print(key, "f6 mediana", round(res["layers"][key]["f6_median_pe_prim"], 3),
                  {k: (round(v) if isinstance(v, float) else v) for k, v in res["layers"][key].items() if k.endswith(("n_net", "tag_events"))}, flush=True)
with open(OUT + "/spectra.csv", "w") as f:
    f.write("# mip_tagged.py: espectros etiquetados (mismo evento − mezclados) en p.e. primarios; amplitud aS/G/F y carga QS/Q1/F; umbral de la "
            f"sonda {PROBE_T:g} p.e. medidos en su variable; etiqueta > {TAG:g} p.e.; F mayo {F['may']}, octubre {F['oct']}\n")
    f.write("campaign,layer,pos_bin,variable,pe_lo,pe_hi,n_same,n_mix_scaled,n_net,err\n")
    for r in rows:
        f.write(",".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
json.dump(res, open(OUT + "/mip_tagged.json", "w"), indent=1)
# comprobación: aS de este archivo = aS de los archivos de impactos de PROC
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import xy_estimator as xy  # noqa: E402
h = xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)
print("comprobación aS octubre: máx |Δ| =", float(np.abs(h["aS"][:, :32] - DATA["oct"]["aS"]).max()))
