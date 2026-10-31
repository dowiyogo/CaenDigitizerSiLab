#!/usr/bin/python3 -I
# ov_correlation.py — seguimiento 2026-10-31 (8), punto 3: ¿la razón mayo/octubre de p.e. primarios por MPPC (ratio_by_mppc.csv del seguimiento 7)
# correlaciona con ΔOV = OV_mayo − OV_oct?
# OV = G/k [supuesto: G = k·OV, recta por el origen en V_bd], con k común = mediana de las estimaciones de dos puntos del punto 1 (vbd.json; solo 5
# canales: no hay k por canal). Con k común, las correlaciones son las mismas que con ΔG; la pendiente escala con 1/k.
#   principal (pedido): G_mayo = reajuste de PROC del canal de mayo del MPPC (pe_primary_channels.csv), G_oct = c1 del canal de octubre;
#   variante (mismo método en las dos campañas): G_oct = reajuste de PROC sobre octubre run 8 (hd64_oct8_peaks.csv, hitdef_extract.C).
# Conjuntos: los 28 MPPC de la exclusión principal del seguimiento 7 (sin pares con ch29 o ch30) y los 31.
# Pearson (con t), Spearman, n y pendiente por mínimos cuadrados ponderados (1/σ² de la razón) y sin pesos.
# Uso: ov_correlation.py <dir de recabling_groups.py> <dir vbd> <salida.json>
import csv, json, math, statistics as st, sys

import numpy as np

RG, VB, OUTF = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(RG + "/ratio_by_mppc.csv")))
PE = {int(r["ch"]): r for r in csv.DictReader(open("outputs/hod2018_mip/2026-10-31_pe_primary/pe_primary_channels.csv"))}
OCT_PROC = {int(r["ch"]): float(r["G_fit"]) for r in csv.DictReader(open("outputs/proceeding_2026_10/work/hd64_oct8_peaks.csv")) if r["status"] in ("ok", "changed")}
k = json.load(open(VB + "/vbd.json"))["summary"]["two_point_median"]["k_mV_per_V"]
rank = lambda v: np.argsort(np.argsort(v)).astype(float)
res = {"k_common_mV_per_V": k, "sets": {}}
for variant in ("c1", "oct_proc"):
    for sname, excl in (("28", True), ("31", False)):
        sel = [r for r in rows if not (excl and (int(r["ch_may"]) in (29, 30) or int(r["ch_oct"]) in (29, 30)))]
        x, y, e, ids = [], [], [], []
        for r in sel:
            cm, co = int(r["ch_may"]), int(r["ch_oct"])
            gm = float(PE[cm]["G_may_mV"])
            go = float(PE[co]["G_oct_mV"]) if variant == "c1" else OCT_PROC.get(co)
            if go is None or (variant == "c1" and co in (22, 31)):  # c1 sin G propia (mediana) -> fuera
                continue
            x.append((gm - go) / k); y.append(float(r["ratio"])); e.append(float(r["ratio_err"])); ids.append(r["id"])
        x, y, e = np.array(x), np.array(y), np.array(e)
        n = len(x)
        pr = float(np.corrcoef(x, y)[0, 1]); sr = float(np.corrcoef(rank(x), rank(y))[0, 1])
        t = pr * math.sqrt((n - 2) / (1 - pr * pr))
        A = np.vstack([x, np.ones_like(x)]).T
        w = 1 / e ** 2
        C = np.linalg.inv(A.T @ (A * w[:, None])); bw = C @ (A.T @ (w * y))
        chi = float((w * (y - A @ bw) ** 2).sum())
        bu, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ bu
        Cu = np.linalg.inv(A.T @ A) * (resid @ resid) / (n - 2)
        d = {"variant": variant, "set": sname, "n": n, "pearson": pr, "pearson_t": t, "spearman": sr,
             "slope_weighted_per_V": float(bw[0]), "slope_weighted_err": float(math.sqrt(C[0, 0])), "slope_weighted_err_scaled": float(math.sqrt(C[0, 0] * max(1, chi / (n - 2)))),
             "chi2": chi, "ndf": n - 2, "slope_unweighted_per_V": float(bu[0]), "slope_unweighted_err": float(math.sqrt(Cu[0, 0])),
             "dOV_median_V": float(np.median(x)), "dOV_range_V": [float(x.min()), float(x.max())], "mppc": ids, "dOV_V": x.tolist(), "ratio": y.tolist()}
        res["sets"][f"{variant}_{sname}"] = d
        print(f"G_oct {variant:8s} conjunto {sname}: n {n}, ΔOV mediana {d['dOV_median_V']:+.3f} V (rango {x.min():+.3f} a {x.max():+.3f}), Pearson {pr:+.3f} (t = {t:.2f}), "
              f"Spearman {sr:+.3f}, pendiente ponderada {bw[0]:+.3f} ± {d['slope_weighted_err_scaled']:.3f} /V (χ² {chi:.1f}/{n - 2}), sin pesos {bu[0]:+.3f} ± {d['slope_unweighted_err']:.3f} /V")
print(f"k común = {k:.3f} mV/V")
json.dump(res, open(OUTF, "w"), indent=1)
