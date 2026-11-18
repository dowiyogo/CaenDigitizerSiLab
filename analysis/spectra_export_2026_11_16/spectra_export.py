#!/usr/bin/python3 -I
# spectra_export.py — seguimiento 2026-11-16 (spectra-export), punto 4: espectros de muones etiquetados de hod2018 en p.e. medidos, por capa
# (XA, XB, YA, YB) e intervalo de distancia al extremo de lectura (21-29, 29-37, 37-45, 45-54 mm), mayo y octubre, para un ajuste de forma completa.
# Selección idéntica a la de los espectros de entrada de may_targets y oct_targets (analysis/followup14_2026_11_04/layer_targets.py, modo spectra):
#  - mayo v8 en tiempo bueno, G = reajuste de PROC (hd64_may_v8_peaks.csv, ok/changed; mediana para los demás); octubre run 8, G = c1 (mediana para
#    ch22, ch31); ch3 muerto; orden de tiras de cada campaña (stage3.json), extremo de lectura de layer_from_channels.json;
#  - etiqueta = alguna tira de la coordenada ortogonal del intervalo con aS > 10 G (X12 = ch3 fuera); sonda = cada tira viva de la capa con aS/G > 3;
#  - mismo evento − mezclados (etiqueta del evento siguiente) × N/(N − 1); error² = mismo + escala² · mezclados.
# Bins de 0.25 p.e. medido de 0 a 60; las entradas sobre 60 van en una fila de desborde (pe_lo = 60, pe_hi = inf).
# Comprobación: sumando de a 4 bins se reproduce spectra_layers.csv de outputs/hod2018_{may,oct}_targets (bins de 1 p.e. bajo 60).
# Salida: <out>/spectra_{may,oct}.csv (campaign, layer, pos_bin, pe_lo, pe_hi, n_same, n_mix_scaled, n_net, err; counts = n_net),
#         <out>/calibration_{may,oct}.csv (ch, G_mV, fuente, capa y tira en el orden de la campaña), spectra_export.json.
# Uso: spectra_export.py <out dir>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
SKEY = {"may": "may2026_v8_goodtime", "oct": "oct2026_run8"}
TARGETS = {"may": REPO + "/outputs/hod2018_may_targets/2026-11-04_bootstrap/spectra_layers.csv",
           "oct": REPO + "/outputs/hod2018_oct_targets/2026-11-04/spectra_layers.csv"}
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
EDG = np.round(np.arange(0, 60.0001, 0.25), 4)
summary = {}
for camp in ("may", "oct"):
    order = ST3[camp]["order"]
    if camp == "may":
        g = np.full(64, np.nan)
        src = ["median of refit"] * 64
        for l in open(PK + "/hd64_may_v8_peaks.csv"):
            f = l.strip().split(",")
            if f[0] != "ch" and f[3] in ("ok", "changed"):
                g[int(f[0])] = float(f[2])
                src[int(f[0])] = f"PROC refit ({f[3]})"
        gf = np.where(np.isfinite(g), g, np.nanmedian(g))
        e = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
        e = e[np.argsort(e["t"], kind="stable")]
        e = e[mc.good_time(e["t"], "v8")]
    else:
        g = mc.g_c1()
        src = ["c1" if np.isfinite(g[c]) else "median of c1" for c in range(len(g))]
        gf = np.where(np.isfinite(g), g, np.nanmedian(g))
        e = xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)
        e = e[np.argsort(e["t"], kind="stable")]
    N = len(e)
    s = N / (N - 1)
    rows, chk = [], 0.0
    ref = {}
    for r in csv.DictReader(l for l in open(TARGETS[camp]) if not l.startswith("#")):
        ref[(r["layer"], r["pos_bin"], float(r["pe_lo"]))] = (float(r["n_same"]), float(r["n_mix_scaled"]), float(r["n_net"]))
    for probe, refc in (("X", "Y"), ("Y", "X")):
        refch = np.asarray(order[f"h18_{refc}"])
        tagS = e["aS"][:, refch] > 10.0 * gf[refch]
        tagS[:, refch == DEAD] = False
        for L in ("A", "B"):
            chans = [c for i, c in enumerate(order[f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
            for lo, hi in BINS:
                dist = lambda j: D[j] if SIGN[f"{SKEY[camp]}/{probe}{L}"]["sign"] < 0 else 75.0 - D[j]  # noqa: E731
                js = [j for j in range(16) if lo <= dist(j) < hi + (1e-9 if hi == 54 else 0)]
                js = [j for j in js if not (probe == "Y" and j == DEAD_X)]
                T = tagS[:, js].any(1)
                Tn = np.concatenate([T[1:], [False]])
                hs, hm = np.zeros(len(EDG)), np.zeros(len(EDG))  # último = desborde
                for c in chans:
                    meas = e["aS"][:, c] / gf[c]
                    sel = meas > 3.0
                    for msk, h in ((sel & T, hs), (sel & Tn, hm)):
                        v = meas[msk]
                        h[:-1] += np.histogram(v[v < 60], bins=EDG)[0]
                        h[-1] += (v >= 60).sum()
                net, err = hs - s * hm, np.sqrt(hs + s * s * hm)
                for k in range(len(EDG)):
                    pl, ph = (EDG[k], EDG[k + 1]) if k < len(EDG) - 1 else (60.0, float("inf"))
                    rows.append([camp, f"{probe}{L}", f"{lo}-{hi}", f"{pl:.2f}", "inf" if ph == float("inf") else f"{ph:.2f}", int(hs[k]), f"{s * hm[k]:.4f}",
                                 f"{net[k]:.4f}", f"{err[k]:.4f}"])
                # comprobación contra spectra_layers.csv (bins de 1 p.e.)
                for b in range(60):
                    a_s, a_m = hs[4 * b:4 * b + 4].sum(), s * hm[4 * b:4 * b + 4].sum()
                    r_s, r_m, _ = ref[(f"{probe}{L}", f"{lo}-{hi}", float(b))]
                    chk = max(chk, abs(a_s - r_s), abs(a_m - r_m))
                summary.setdefault(camp, {})[f"{probe}{L}/{lo}-{hi}"] = {"n_same": float(hs.sum()), "n_mix_scaled": float(s * hm.sum()), "n_net": float(net.sum()),
                                                                          "n_net_overflow_ge60": float(net[-1]), "probe_channels": [int(c) for c in chans],
                                                                          "tag_channels": [int(refch[j]) for j in js]}
    summary[camp + "_events"] = N
    summary[camp + "_max_abs_diff_vs_targets_1pe_bins"] = float(chk)
    print(f"{camp}: eventos {N}; máx |Δ| contra spectra_layers.csv (bins de 1 p.e. bajo 60): {chk:.6f}", flush=True)
    with open(f"{OUT}/spectra_{camp}.csv", "w") as f:
        f.write(f"# hod2018 {camp}: tagged-muon spectra in MEASURED p.e. (aS/G, not divided by crosstalk F), selection of layer_targets.py (tag > 10 p.e. in "
                f"the position bin, probe > 3 p.e.), net = same - mixed x N/(N-1) = {s:.8f}; counts = n_net; err = sqrt(same + scale^2 mixed); "
                f"last row of each block = overflow (>= 60 p.e.); {N} triggers\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["campaign", "layer", "pos_bin_mm", "pe_lo", "pe_hi", "n_same", "n_mix_scaled", "n_net", "err"])
        w.writerows(rows)
    lay = {}
    for probe in ("X", "Y"):
        for i, c in enumerate(order[f"h18_{probe}"]):
            lay[int(c)] = (f"{probe}{'A' if i % 2 == 0 else 'B'}", i)
    with open(f"{OUT}/calibration_{camp}.csv", "w") as f:
        f.write(f"# hod2018 {camp}: single-p.e. gain G (mV per p.e., peak spacing) used for the measured-p.e. scale aS/G; ch3 dead\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ch", "G_mV", "source", "layer", "strip_index", "note"])
        for c in range(32):
            w.writerow([c, f"{gf[c]:.5f}", src[c], lay[c][0], lay[c][1], "dead" if c == DEAD else ""])
json.dump(summary, open(OUT + "/spectra_export.json", "w"), indent=1)
