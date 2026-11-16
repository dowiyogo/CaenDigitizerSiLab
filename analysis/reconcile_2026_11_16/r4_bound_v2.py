#!/usr/bin/python3 -I
# r4_bound_v2.py — seguimiento 2026-11-16, punto 4: cota superior R_prim/R_G por definición de F, versión 2.
# reconcile.py escribió r4_bound.csv leyendo xtalk_ratio.csv con csv.DictReader; ese archivo tiene dos columnas llamadas «err» (error de F_oct/F_mayo
# y error de la razón primaria octubre/mayo... y una tercera del inverso) y DictReader tomó la última (error de mayo/octubre). Los valores centrales
# de r4_bound.csv son correctos; la columna de error no. Aquí se leen las columnas por posición: 7 = prim_ratio_oct_over_may, 8 = su error.
# R_G = mediana de G_oct/G_mayo por MPPC (r4_summary.json, R_G_median). Salida: <out>/r4_bound_v2.csv. Uso: r4_bound_v2.py <out dir>
import csv, json, os, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
s = json.load(open(OUT + "/r4_summary.json"))
RG = s["R_G_median"]
r1 = json.load(open(OUT + "/r1_summary.json"))["check_0p6477"]
rows = list(csv.reader(open(REPO + "/outputs/hod2018_xtalk/2026-11-15_borel/xtalk_ratio.csv")))
hdr = rows[0]
assert hdr[7] == "prim_ratio_oct_over_may" and hdr[8] == "err", hdr
out = []
for r in rows[1:]:
    rp, erp = float(r[7]), float(r[8])
    out.append([r[0], f"{rp:.4f}", f"{erp:.4f}", f"{RG:.4f}", f"{rp / RG:.4f}", f"{erp / RG:.4f}"])
rm = r1["oct"] / r1["may"]
rme = rm * ((0.4458 / 16.1627) ** 2 + (0.2954 / 10.4680) ** 2) ** 0.5  # errores de las medianas de ampcal [M]
rp = rm * s["F_may"] / s["F_oct"]
out.append([f"muon spectrum, fit R (p {s['p_may']:.4f} / {s['p_oct']:.4f})", f"{rp:.4f}", f"{rme * s['F_may'] / s['F_oct']:.4f}", f"{RG:.4f}", f"{rp / RG:.4f}",
            f"{rme * s['F_may'] / s['F_oct'] / RG:.4f}"])
with open(OUT + "/r4_bound_v2.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["F_definition", "R_prim_oct_over_may", "err", "R_G_oct_over_may_median", "R_prim_over_R_G_upper_bound", "err"])
    w.writerows(out)
for r in out:
    print(" | ".join(r))
