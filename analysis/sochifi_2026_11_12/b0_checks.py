#!/usr/bin/python3 -I
# b0_checks.py — SOCHIFI 2026-11-12, precondiciones B0.2 y B0.3 (solo lectura):
#  - q y ε_XY = q/(f_2018·g) con el estimador compartido (xy_estimator.xy_fraction, mezclados) en el tiempo bueno de run B, en el segmento
#    estacionario (t ≥ 5.25 h) y en su intersección, a 5 y 10 p.e.; f_2018 de geom_sens_CAD_1e8.json (cos³), g = 0.95512 (g_dead_ch3.py);
#  - exceso de coincidencias entre tiras (stage1.json, strip_order_matrices: co − mezclados, z = E/√(co + mix)) en las uniones del orden de
#    mayo que las cadenas medidas no fijan, y en sus alternativas;
#  - grupos HV de mayo contra octubre en las tiras Y 12 y 15 con el segundo tramo de Y en el orden de mayo y con ese tramo invertido.
# Uso: b0_checks.py
import json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
g = np.full(32, np.nan)
for l in open(PK + "/hd64_may_v8_peaks.csv"):
    f = l.strip().split(",")
    if f[0] != "ch" and f[3] in ("ok", "changed"):
        g[int(f[0])] = float(f[2])
g = np.where(np.isfinite(g), g, np.nanmedian(g))
G = np.concatenate([g, np.full(32, np.nan)])
e = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
F18 = json.load(open(PK + "/geom_sens_CAD_1e8.json"))["cos3"]["nominal"]["f"]
GG = 0.95512
t0 = e["t"].min()
st, gt = e["t"] >= t0 + 5.25 * 3600, mc.good_time(e["t"], "v8")
XCH, YCH = [c for c in range(16) if c != 3], list(range(16, 32))
print(f"eventos v8 {len(e)}; f_2018 = {F18:.6f}; g = {GG}")
for name, m in (("tiempo bueno", gt), ("estacionario", st), ("estacionario ∩ tiempo bueno", gt & st)):
    ee = e[m]
    for T in (5, 10):
        Xs, Xc = xy.plane_hits(ee, G, T, XCH)
        Ys, Yc = xy.plane_hits(ee, G, T, YCH)
        r = xy.xy_fraction(Xs, Ys, Xc, Yc)
        q, qe, qc = r["f_mix"], r["f_mix_err"], r["f_ctl"]
        k = F18 * GG
        print(f"{name}: N {len(ee)} T {T}: N_XY {r['N_XY']} acc_mix {r['acc_mix']:.1f} acc_ctl {r['acc_ctl']} q {q:.5f} ± {qe:.5f} (ctl {qc:.5f}) "
              f"eps_XY {q / k:.4f} ± {qe / k:.4f} (stat) ± {abs(q - qc) / k:.4f} (estimador) ± {q / k * 0.0465:.4f} (aceptancia 4.65 %)")
S = json.load(open(PK + "/stage1.json"))["strip_order_matrices"]
for plane, pairs in (("hod2018 X", [(2, 11), (2, 4), (15, 11), (15, 4), (4, 9), (4, 10), (11, 9), (11, 10), (10, 12), (10, 14), (9, 12), (9, 14)]),
                     ("hod2018 Y", [(21, 29), (21, 31), (28, 29), (28, 31)])):
    d = S[plane]
    ix = {c: i for i, c in enumerate(d["channels"])}
    co, mx = np.array(d["co"], float), np.array(d["mixed"], float)
    for a, b in pairs:
        i, j = ix[a], ix[b]
        E = co[i, j] - mx[i, j]
        print(f"{plane} ({a},{b}): co {co[i, j]:.0f} mezclados {mx[i, j]:.1f} E {E:.1f} z {E / np.sqrt(max(co[i, j] + mx[i, j], 1)):.2f}")
may_y = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]["h18_Y"]
rev = may_y[:11] + may_y[11:][::-1]
grp = {}
for l in open(os.path.join(REPO, "outputs/hod2018_recabling/2026-10-31_groups/recabling_groups.log")):
    if "'bar_may'" in l:
        d = eval(l.strip())  # líneas de dict de Python escritas por recabling_groups.py
        grp[d["ch_may"]] = d["hv_group_may"]
oct_y = json.load(open(PK + "/stage3.json"))["strip_orders"]["oct"]["order"]["h18_Y"]
oct_grp = {}
import csv  # noqa: E402
for r in csv.DictReader(open(os.path.join(REPO, "outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv"))):
    if r["plane"] == "Y":
        oct_grp[int(r["ch_oct"])] = int(r["hv_group"])
for lab, order in (("orden de mayo", may_y), ("segundo tramo invertido", rev)):
    res = []
    for pos, ch in enumerate(order):
        if ch in grp:
            res.append(f"ch{ch} pos {pos}: grupo mayo {grp[ch]}, grupo octubre en pos {pos} (ch{oct_y[pos]}) {oct_grp[oct_y[pos]]}")
    print(lab + ": " + "; ".join(res))
