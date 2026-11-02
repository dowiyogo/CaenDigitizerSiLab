#!/usr/bin/python3 -I
# mip_paper.py — 2026-11-02: espectros MIP etiquetados de hod2018 en mayo (v8, tiempo bueno), AMPLITUD en p.e. MEDIDOS (aS/G, sin dividir
# por F), para la figura del artículo (mip_paper.C). Misma selección que analysis/followup10_2026_11_01/mip_tagged.py:
# etiqueta = alguna tira de la coordenada ortogonal del intervalo con aS > 10 G (ch3 fuera); sonda = cada tira viva de la capa con aS > 3 G;
# mezclados = etiqueta del evento siguiente, × N/(N − 1); neto = mismo − mezclados, error² = mismo + escala² · mezclados. Solo cambia el
# histograma: bins de 1 p.e. medido en [0, 300) en lugar de 1 p.e. primario.
# Entradas: <tagged>/work/may_v8_charge.bin (charge_extract.C) y las de mip_tagged.py (G reajustada de PROC, orden de tiras, extremo de lectura).
# Comprobación: los mismos valores divididos por F = 1.75 y con los bins de mip_tagged.py reproducen <tagged>/spectra.csv (mayo, amplitud).
# Salida: <out>/spectra_meas.csv. Uso: mip_paper.py <dir de 2026-11-01_tagged> <out dir>
import json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

TAGGED, OUT = sys.argv[1], sys.argv[2]
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ORDER = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD = 3
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = {"near": (21, 29), "far": (29, 54)}
TAG, PROBE_T, F = 10.0, 3.0, 1.75
EDG = np.arange(0, 301, 1.0)
EDG_PRIM = np.arange(0, 151, 1.0)
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("aC", "<f4", 32), ("QS", "<f4", 32), ("QC", "<f4", 32)])


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


g = g_refit()
g = np.where(np.isfinite(g), g, np.nanmedian(g))
e = np.fromfile(TAGGED + "/work/may_v8_charge.bin", dtype=DT)
e = e[np.argsort(e["t"], kind="stable")]
e = e[mc.good_time(e["t"], "v8")]
N = len(e)
s = N / (N - 1)
print("eventos de mayo v8 en tiempo bueno:", N)

# referencia: spectra.csv de mip_tagged.py (mayo, amplitud, p.e. primarios)
ref = {}
for l in open(TAGGED + "/spectra.csv"):
    f = l.strip().split(",")
    if f[0] == "may" and f[3] == "amp":
        ref.setdefault((f[1], f[2]), []).append((int(f[6]), float(f[7])))

rows, maxdiff = [], 0.0
for probe, refc in (("X", "Y"), ("Y", "X")):
    refch = np.asarray(ORDER[f"h18_{refc}"])
    tagS = e["aS"][:, refch] > TAG * g[refch]
    tagS[:, refch == DEAD] = False
    for L in ("A", "B"):
        chans = [c for i, c in enumerate(ORDER[f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
        sg = SIGN[f"may2026_v8_goodtime/{probe}{L}"]["sign"]
        d = D if sg < 0 else 75.0 - D
        for b, (lo, hi) in BINS.items():
            js = [j for j in range(16) if lo <= d[j] < hi + (1e-9 if hi == 54 else 0)]
            T = tagS[:, js].any(1)
            Tn = np.concatenate([T[1:], [False]])
            hs, hm = np.zeros(len(EDG) - 1), np.zeros(len(EDG) - 1)
            ps, pm = np.zeros(len(EDG_PRIM) - 1), np.zeros(len(EDG_PRIM) - 1)
            for c in chans:
                meas = e["aS"][:, c] / g[c]
                sel = meas > PROBE_T
                hs += np.histogram(meas[sel & T], bins=EDG)[0]
                hm += np.histogram(meas[sel & Tn], bins=EDG)[0]
                ps += np.histogram(meas[sel & T] / F, bins=EDG_PRIM)[0]
                pm += np.histogram(meas[sel & Tn] / F, bins=EDG_PRIM)[0]
            r = np.array([x[0] for x in ref[(probe + L, b)]]), np.array([x[1] for x in ref[(probe + L, b)]])
            maxdiff = max(maxdiff, float(np.abs(r[0] - ps).max()), float(np.abs(r[1] - s * pm).max()))
            net, err = hs - s * hm, np.sqrt(hs + s * s * hm)
            print(f"{probe}{L} {b}: n_net {net.sum():.1f} (mismo {hs.sum():.0f}, mezclados {s * hm.sum():.1f})", flush=True)
            for k in range(len(net)):
                rows.append([f"{probe}{L}", b, EDG[k], EDG[k + 1], int(hs[k]), s * hm[k], net[k], err[k]])
print("comprobación contra spectra.csv (mayo, amplitud, p.e. primarios): máx |Δ| =", maxdiff)
os.makedirs(OUT, exist_ok=True)
with open(OUT + "/spectra_meas.csv", "w") as f:
    f.write(f"# mip_paper.py: mayo v8 tiempo bueno, amplitud aS/G en p.e. medidos (sin F), etiqueta > {TAG:g} p.e., sonda > {PROBE_T:g} p.e.; "
            "neto = mismo − mezclados\n")
    f.write("layer,pos_bin,pe_lo,pe_hi,n_same,n_mix_scaled,n_net,err\n")
    for r in rows:
        f.write(",".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
