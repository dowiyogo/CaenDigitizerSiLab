#!/usr/bin/python3 -I
# compare_offsets.py — seguimiento 2026-11-15: desplazamientos cfd50 por canal de octubre run 8 (time_res2.py oct) contra los de mayo (seguimiento
# 13 A), emparejados por barra (posición de la tira en su plano, stage3.json strip_orders; el Y de hod2018 se recableó, así que el mismo canal no es
# la misma MPPC), y la resolución σ(t_X − t_Y) por bin de amplitud (cfd50, correcciones none / offsets / offsets+prop) de las dos campañas.
# Salida: <out>/offsets_may_vs_oct.csv, time_res_may_vs_oct.csv; imprime el resumen. Uso: compare_offsets.py <may A_time> <oct A_time> <out>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MA, OA, OUT = sys.argv[1:4]
ORD = json.load(open(os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work/stage3.json")))["strip_orders"]
pos = {k: {"X": ORD[k]["order"]["h18_X"], "Y": ORD[k]["order"]["h18_Y"]} for k in ("may", "oct")}


def offs(d):
    return {int(r["ch"]): float(r["offset_ns"]) for r in csv.DictReader(open(d + "/time_res_offsets.csv")) if r["method"] == "cfd50"}


om, oo = offs(MA), offs(OA)
rows, dm, do = [], [], []
for coord in ("X", "Y"):
    for j in range(16):
        cm, co = pos["may"][coord][j], pos["oct"][coord][j]
        if 3 in (cm, co):
            continue
        rows.append([coord, j, cm, co, f"{om[cm]:.4f}", f"{oo[co]:.4f}", f"{oo[co] - om[cm]:.4f}"])
        dm.append(om[cm])
        do.append(oo[co])
dm, do = np.array(dm), np.array(do)
with open(f"{OUT}/offsets_may_vs_oct.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["coord", "strip_position", "ch_may", "ch_oct", "offset_may_ns", "offset_oct_ns", "oct_minus_may_ns"])
    w.writerows(rows)
print(f"barras {len(dm)}: rms desplazamiento mayo {dm.std():.3f} ns, octubre {do.std():.3f} ns, rms(oct − may) {(do - dm).std():.3f} ns, "
      f"correlación {np.corrcoef(dm, do)[0, 1]:.3f}")
tr = []
for camp, d in (("may", MA), ("oct", OA)):
    for r in csv.reader(l for l in open(d + "/time_res.csv") if not l.startswith("#")):
        if r[0] == "cfd50":
            tr.append([camp] + r)
with open(f"{OUT}/time_res_may_vs_oct.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "method", "correction", "min_amp_bin_pe", "n_same", "n_mix_scaled", "n_net", "mu_ns", "mu_err_ns", "sigma_ns", "sigma_err_ns",
                "rms_3sigma_ns", "sigma_strip_ns", "sigma_strip_err_ns", "rms_strip_ns", "chi2", "ndf"])
    w.writerows(tr)
for r in tr:
    print(f"{r[0]} {r[2]:13s} {r[3]:6s} n={r[4]:>5s} sigma={float(r[9]):.3f}±{float(r[10]):.3f} rms3={float(r[11]):.3f}")
