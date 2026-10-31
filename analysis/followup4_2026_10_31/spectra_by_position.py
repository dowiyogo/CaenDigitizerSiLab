#!/usr/bin/python3 -I
# spectra_by_position.py — seguimiento 2026-10-31 (4), punto 2: espectro de amplitud de coincidencias X–Y verdaderas de hod2018, octubre run 8,
# por capa y por posición a lo largo de la barra en 4 intervalos de distancia al extremo de lectura: [21, 29), [29, 37), [37, 45), [45, 54] mm.
# Misma selección y definición de amplitud que amp_vs_position.py (seguimiento 3): amplitud de pico en la ventana de 144 ns en p.e. (G de c1;
# mediana para ch22 y ch31), sonda P > 5 p.e., alguna tira ortogonal > 5 p.e., referencia = la de mayor amplitud, posición de la referencia
# d = 37.5 + (−15 + 2j) mm desde el extremo del índice 0 (orden de r4); mezclados = referencia del evento siguiente × N/(N − 1); resta.
# Extremo de lectura de cada capa [supuesto: el extremo más brillante es el del MPPC] = el del signo medido en el seguimiento 3
# (layer_from_channels.json, octubre): XA y YB bajan con d -> lectura en el extremo del índice 0, distancia = d; XB e YA suben con d ->
# lectura en el extremo opuesto, distancia = 75 − d. Excluidos: ch3 y ch29 (octubre), como sonda y como referencia.
# Bins de 1 p.e. de 0 a 300 (el recorte de cada canal está en clip.csv del seguimiento 3, 138-290 p.e.).
# Salida: <out>/spectrum_<capa>_pos<k>.csv (k = 0..3), <out>/summary.json. Uso: spectra_by_position.py <out dir>
import json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
ST3 = json.load(open(os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work/stage3.json")))["strip_orders"]["oct"]["order"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
CLIP = {}
for l in open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/oct2026_run8/clip.csv")):
    if not l.startswith("ch,"):
        f = l.strip().split(",")
        CLIP[int(f[0])] = float(f[3])
T = 5.0
EXCL = {3, 29}
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
E1 = np.arange(0, 301, 1.0)
g = mc.g_c1()
G = np.where(np.isfinite(g), g, np.nanmedian(g))
e = xy.load_hits(os.path.join(W, "hd64_oct8_hits.bin"), 64)
N = len(e)
s = N / (N - 1)
pe = e["aS"][:, :32] / G[None, :32]
hit = pe > T
for c in EXCL:
    hit[:, c] = False
order = {"X": ST3["h18_X"], "Y": ST3["h18_Y"]}
summary = {"run": "octubre run 8", "N_triggers": int(N), "T_pe": T, "excluded": sorted(EXCL), "position_bins_mm_from_readout": BINS, "layers": {}}
for coord, other in (("X", "Y"), ("Y", "X")):
    ref = order[other]
    refpe = np.where(hit[:, ref], pe[:, ref], -1.0)
    anyref = (refpe > 0).any(1)
    jbest = np.argmax(refpe, axis=1)
    anyref_n = np.concatenate([anyref[1:], [False]])
    jbest_n = np.concatenate([jbest[1:], [0]])
    for L in ("A", "B"):
        sign = SIGN[f"oct2026_run8/{coord}{L}"]["sign"]
        dist = D if sign < 0 else 75.0 - D  # distancia al extremo de lectura de cada posición j
        chans = [c for i, c in enumerate(order[coord]) if (i % 2 == (0 if L == "A" else 1)) and c not in EXCL]
        lay = {"sign_vs_d": sign, "readout_end": "índice 0" if sign < 0 else "índice 15", "channels": chans,
               "clip_pe_min": min(CLIP[c] for c in chans), "bins": []}
        for k, (lo, hi) in enumerate(BINS):
            js = [j for j in range(16) if lo <= dist[j] < hi + (1e-9 if hi == 54 else 0) and ref[j] not in EXCL]
            hs = np.zeros(300); hm = np.zeros(300)
            for c in chans:
                ms = hit[:, c] & anyref & np.isin(jbest, js)
                mm = hit[:, c] & anyref_n & np.isin(jbest_n, js)
                hs += np.histogram(pe[ms, c], bins=E1)[0]
                hm += np.histogram(pe[mm, c], bins=E1)[0]
            net = hs - s * hm
            err = np.sqrt(hs + s * s * hm)
            n = np.clip(net, 0, None)
            cs = np.cumsum(n)
            med = float(E1[np.searchsorted(cs, 0.5 * cs[-1])] + 0.5) if cs[-1] > 0 else None
            mean = float((n * (E1[:-1] + 0.5)).sum() / n.sum()) if n.sum() > 0 else None
            with open(f"{OUT}/spectrum_{coord}{L}_pos{k}.csv", "w") as f:
                f.write(f"# octubre run 8, capa {coord}{L} (canales {chans}), distancia al extremo de lectura [{lo}, {hi}) mm, referencias j = {js} "
                        f"(ch {[ref[j] for j in js]}); mismo evento − mezclados × {s:.6f}; recorte mínimo de la capa {lay['clip_pe_min']:.0f} p.e.\n")
                f.write("pe_lo,pe_hi,n_same,n_mix_scaled,n_net,n_net_err\n")
                for i in range(300):
                    f.write(f"{i},{i + 1},{int(hs[i])},{s * hm[i]:.3f},{net[i]:.3f},{err[i]:.3f}\n")
            lay["bins"].append({"range_mm": [lo, hi], "ref_j": js, "n_same": int(hs.sum()), "n_mix_scaled": float(s * hm.sum()), "n_net": float(net.sum()),
                                "n_net_err": float(math.sqrt(hs.sum() + s * s * hm.sum())), "median_net_pe": med, "mean_net_pe": mean})
        summary["layers"][f"{coord}{L}"] = lay
        print(f"{coord}{L}: lectura en el extremo del {lay['readout_end']}; " + "; ".join(
            f"[{b['range_mm'][0]},{b['range_mm'][1]}) netos {b['n_net']:.0f} ± {b['n_net_err']:.0f}, mediana {b['median_net_pe']}, media {b['mean_net_pe']:.1f}" for b in lay["bins"]))
json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1)
