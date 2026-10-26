#!/usr/bin/python3 -I
# export_inputs.py — XPLUS-4b: entradas comunes para las macros: intervalos UTC excluidos del tiempo bueno de un run de mayo
# (may_common.good_intervals: EXCLUDE y HV_OFF, hora local UTC−4 → UTC) y vectores de G de hod2018 (32 valores, coma):
#   may_refit = G reajustada de work/hd64_may_v8_peaks.csv (canales sin G = mediana, como stage2.py), c1 = pe_calibration_2026_10_c1
#   (mediana para los que faltan), sw0 = G "nuevo" a 0 V de DAQ2-15 (report_es.md; sin G o marcadas "*" = mediana),
#   swm05 = G "nuevo" a −0.5 V (sin G = 0.5758·G(0 V), supuesto).
# Uso: export_inputs.py excluded <run> <salida.txt> | gains <may_refit|c1|sw0|swm05>
import os, re, statistics, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
if sys.argv[1] == "excluded":
    with open(sys.argv[3], "w") as f:
        for _, a, b, _, _ in mc.good_intervals(sys.argv[2]):
            f.write(f"{a + 4 * 3600} {b + 4 * 3600}\n")
    print(open(sys.argv[3]).read().strip())
    sys.exit(0)
kind = sys.argv[2]
if kind == "may_refit":
    g = np.full(32, np.nan)
    for l in open(os.path.join(REPO, "outputs/proceeding_2026_10/work/hd64_may_v8_peaks.csv")):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
elif kind == "c1":
    g = mc.g_c1()
else:
    g0, g05 = np.full(32, np.nan), np.full(32, np.nan)
    for line in open(os.path.join(REPO, "outputs/daq2-15/biasscan_2026_10_23_day/report_es.md")):
        m = re.match(r"\|\s*(\d+)\s*\|\s*([^|]*)\|\s*([^|]*)\|\s*([^|]*)\|\s*([^|]*)\|", line)
        if m and line.count("|") >= 11 and int(m.group(1)) < 32:
            c = int(m.group(1))
            for arr, s in ((g0, m.group(3)), (g05, m.group(5))):
                s = s.strip()
                if re.match(r"^[0-9.]+$", s):
                    arr[c] = float(s)
    g = g0 if kind == "sw0" else np.where(np.isfinite(g05), g05, 0.5758 * fill(g0))
print(",".join(f"{x:.6f}" for x in fill(g)))
