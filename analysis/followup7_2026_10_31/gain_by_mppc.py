#!/usr/bin/python3 -I
# gain_by_mppc.py — seguimiento 2026-10-31 (7), punto 1: razón de ganancias G_mayo/G_oct por MPPC (Y emparejado por barra, como en
# recabling_groups.py) contra la del mismo canal del digitalizador usada en los seguimientos 1 y 6, y los MPPC con |G_mayo/G_oct − 1| ≤ 5 %.
# G: pe_primary_channels.csv del seguimiento 6 (mayo = reajuste de PROC; octubre = c1, mediana para ch22 y ch31).
# Uso: gain_by_mppc.py <dir de recabling_groups.py>
import csv, statistics as st, sys

PE = {int(r["ch"]): r for r in csv.DictReader(open("outputs/hod2018_mip/2026-10-31_pe_primary/pe_primary_channels.csv"))}
rows = list(csv.DictReader(open(sys.argv[1] + "/ratio_by_mppc.csv")))
print("MPPC,ch_may,ch_oct,G_may_mV,G_oct_mV,G_ratio_mppc,G_ratio_same_channel,ratio_prim_mppc,ratio_prim_err,oct_G_is_median")
gm, gy = [], []
for r in rows:
    cm, co = int(r["ch_may"]), int(r["ch_oct"])
    g = float(PE[cm]["G_may_mV"]) / float(PE[co]["G_oct_mV"])
    (gy if r["plane"] == "Y" else gm).append(g)
    print(f"{r['id']},{cm},{co},{float(PE[cm]['G_may_mV']):.4f},{float(PE[co]['G_oct_mV']):.4f},{g:.4f},{float(PE[co]['G_ratio']):.4f},"
          f"{float(r['ratio']):.4f},{float(r['ratio_err']):.4f},{co in (22, 31)}")
print(f"# mediana G_mayo/G_oct por MPPC: X {st.median(gm):.4f} (n {len(gm)}), Y {st.median(gy):.4f} (n {len(gy)})")
print("# MPPC con |G_mayo/G_oct − 1| ≤ 5 %:", [r["id"] for r in rows if abs(float(PE[int(r["ch_may"])]["G_may_mV"]) / float(PE[int(r["ch_oct"])]["G_oct_mV"]) - 1) <= 0.05])
