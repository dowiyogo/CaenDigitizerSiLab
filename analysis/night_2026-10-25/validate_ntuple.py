#!/usr/bin/python3 -I
# validate_ntuple.py — XPLUS-4b etapa 1: validación obligatoria de la ntupla antes de producir. Sobre el _hits.bin de ntuple_extract.C
# aplica la carga de stage2.py (G de mayo reajustada de work/hd64_may_v8_peaks.csv o c1, tiempo bueno de may_common.good_time, XCH
# sin ch3, YCH 16-31, T = 5 p.e.) y xy_fraction, y exige igualdad exacta con "expected" del fixture de XPLUS-3.
# Uso: validate_ntuple.py <hits.bin> <may8_good|oct8>
import base64, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

HITS, KEY = sys.argv[1], sys.argv[2]
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
XCH = [c for c in range(16) if c != 3]
YCH = list(range(16, 32))


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
e = xy.load_hits(HITS, 64)
if KEY == "may8_good":
    e = e[mc.good_time(e["t"], "v8")]
    G = np.concatenate([fill(g_refit()), np.full(32, np.nan)])
    fx = json.load(open(os.path.join(REPO, "docs/split/fixtures/xy_estimator_may8good_T5.json")))
else:
    G = np.concatenate([fill(mc.g_c1()), np.full(32, np.nan)])
    fx = json.load(open(os.path.join(REPO, "docs/split/fixtures/xy_estimator_run8_T5.json")))
Xs, Xc = xy.plane_hits(e, G, 5, XCH)
Ys, Yc = xy.plane_hits(e, G, 5, YCH)
r = xy.xy_fraction(Xs, Ys, Xc, Yc)
bad = {k: (r[k], fx["expected"][k]) for k in fx["expected"] if r[k] != fx["expected"][k]}
print(f"{KEY}: N {r['N']} N_XY {r['N_XY']} q_mix {r['f_mix']!r} ± {r['f_mix_err']!r}; distintas del fixture: {bad or 'ninguna'}")
sys.exit(1 if bad else 0)
