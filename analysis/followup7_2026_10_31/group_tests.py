#!/usr/bin/python3 -I
# group_tests.py — seguimiento 2026-10-31 (7), punto 1: pruebas entre grupos sobre la salida de recabling_groups.py.
#  - χ² entre grupos de las medias ponderadas (con el error escalado de cada grupo) contra la media común, para la razón y ΔT(6), por agrupación,
#    con su probabilidad (chi2 con k − 1 grados de libertad);
#  - diferencia A − B y su z;
#  - correlación (Pearson y Spearman) por MPPC entre la razón, ΔT(6) y la razón de ganancias G_mayo/G_oct del canal de cada campaña.
# Exclusión: la principal de recabling_groups.py (pares con ch29 o ch30). Uso: group_tests.py <dir de recabling_groups.py>
import csv, json, math, sys

import numpy as np

D = sys.argv[1]
R = json.load(open(D + "/recabling_groups.json"))
PE = {int(r["ch"]): r for r in csv.DictReader(open("outputs/hod2018_mip/2026-10-31_pe_primary/pe_primary_channels.csv"))}


def chi2_sf(x, k):
    """cola superior de chi² con k grados de libertad: 1 − P(k/2, x/2) con la serie de la gamma incompleta inferior."""
    if k <= 0:
        return float("nan")
    if x <= 0:
        return 1.0
    a, s = k / 2, x / 2
    term = tot = 1 / a
    n = 1
    while term > 1e-16 * tot:
        term *= s / (a + n)
        tot += term
        n += 1
    return 1 - math.exp(-s + a * math.log(s) - math.lgamma(a)) * tot


def between(groups):
    m = np.array([g["mean"] for g in groups.values() if g["n"] > 0])
    e = np.array([g["err_scaled"] for g in groups.values() if g["n"] > 0])
    w = 1 / e ** 2
    mu = (w * m).sum() / w.sum()
    x = float((w * (m - mu) ** 2).sum())
    k = len(m) - 1
    return x, k, chi2_sf(x, k)


for what, src in (("razón", R["ratio"]["main"]["by_mppc"]), ("ΔT(6)", R["shape"]["main"]["dT6"])):
    for key in ("plane", "layer", "hv_group", "adc_group_oct", "pos_oct", "adc_group_may", "pos_may"):
        x, k, p = between(src[key])
        print(f"{what} entre grupos de {key}: χ² = {x:.2f} / {k}, p = {p:.3g}")
    a, b = src["layer"]["A"], src["layer"]["B"]
    print(f"{what} A − B = {a['mean'] - b['mean']:+.4f} ± {math.hypot(a['err_scaled'], b['err_scaled']):.4f} (z = {(a['mean'] - b['mean']) / math.hypot(a['err_scaled'], b['err_scaled']):.2f})")
rows = [r for r in csv.DictReader(open(D + "/ratio_by_mppc.csv")) if not (int(r["ch_may"]) in (29, 30) or int(r["ch_oct"]) in (29, 30))]
sh = {r["id"]: r for r in csv.DictReader(open(D + "/shape_by_mppc.csv"))}
ratio = np.array([float(r["ratio"]) for r in rows])
dt6 = np.array([float(sh[r["id"]]["dT6"]) for r in rows])
gr = np.array([float(PE[int(r["ch_may"])]["G_may_mV"]) / float(PE[int(r["ch_oct"])]["G_oct_mV"]) for r in rows])
rank = lambda v: np.argsort(np.argsort(v)).astype(float)
for n1, a, n2, b in (("razón", ratio, "ΔT(6)", dt6), ("razón", ratio, "G_mayo/G_oct", gr), ("ΔT(6)", dt6, "G_mayo/G_oct", gr)):
    r = np.corrcoef(a, b)[0, 1]
    rs = np.corrcoef(rank(a), rank(b))[0, 1]
    t = r * math.sqrt((len(a) - 2) / (1 - r * r))
    print(f"correlación {n1} – {n2}: Pearson {r:+.3f} (t = {t:.2f}, n = {len(a)}), Spearman {rs:+.3f}")
print("ΔT(6) > 0 en", int((dt6 > 0).sum()), "de", len(dt6), "MPPC (sin exclusión: ver shape_by_mppc.csv); mín", f"{dt6.min():.4f}", "máx", f"{dt6.max():.4f}")
allsh = [float(r["dT6"]) for r in sh.values()]
print("todos los 31 MPPC: ΔT(6) > 0 en", sum(v > 0 for v in allsh), "de", len(allsh), f"(mín {min(allsh):.4f}, máx {max(allsh):.4f})")
