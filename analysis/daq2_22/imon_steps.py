#!/usr/bin/python3 -I
# imon_steps.py — DAQ2-22 etapa 4 (solo lectura de archivos): escalones de IMon de CH00-07 desde el inicio de run 10, del slow control de
# 1 s (/data/daq/slowcontrol/YYYY/hv_YYYYMMDD.csv[.gz]; el log de vigilancia de 15 min no resuelve cambios de menos de 1 min).
# Escalón: |mediana de los 30 s siguientes − mediana de los 30 s anteriores| > 0.3 µA, con las ventanas separadas por ≤ 60 s;
# detecciones del mismo canal a menos de 5 min se agrupan (se informa la de mayor salto). Para cada uno: hora, canal, IMon antes y
# después, VMon antes y después, y si el nivel volvió (mediana de 10-20 min después contra la de antes).
# Además IMon contra aX por bloque (temp_corr_blocks_run010.json, temp_corr_blocks_run008.json): pendiente de ln IMon contra ln aX por canal.
# Uso: imon_steps.py <t0 UTC epoch> <salida dir>
import gzip, json, math, os, statistics as st, sys, time

import numpy as np

T0 = float(sys.argv[1])
OD = sys.argv[2]
T1 = time.time()
rows = {c: [] for c in range(8)}
d = T0 - 86400
days = []
while d <= T1 + 86400:
    days.append(time.strftime("%Y%m%d", time.gmtime(d)))
    d += 86400
for day in sorted(set(days)):
    for fn in (f"/data/daq/slowcontrol/{day[:4]}/hv_{day}.csv.gz", f"/data/daq/slowcontrol/{day[:4]}/hv_{day}.csv"):
        if not os.path.exists(fn):
            continue
        op = gzip.open if fn.endswith(".gz") else open
        with op(fn, "rt") as f:
            h = next(f).strip().split(",")
            it, ic, ii, iv = h.index("t_req_utc_ns"), h.index("ch"), h.index("IMon_A"), h.index("VMon_V")
            for l in f:
                r = l.rstrip().split(",")
                try:
                    c = int(r[ic])
                except (ValueError, IndexError):
                    continue
                if c < 8 and r[ii] and r[iv]:
                    t = int(r[it]) * 1e-9
                    if t >= T0:
                        rows[c].append((t, float(r[ii]) * 1e6, float(r[iv])))
        break
steps, summary = [], {}
for c in range(8):
    a = np.array(sorted(rows[c]))
    if len(a) == 0:
        continue
    t, I, V = a[:, 0], a[:, 1], a[:, 2]
    summary[c] = {"n": int(len(t)), "first_utc": float(t[0]), "last_utc": float(t[-1]), "IMon_first_min_median": float(np.median(I[t < t[0] + 60])),
                  "IMon_last_min_median": float(np.median(I[t > t[-1] - 60])), "IMon_min": float(I.min()), "IMon_max": float(I.max())}
    # medianas de 30 s en una rejilla de 10 s
    grid = np.arange(t[0] + 30, t[-1] - 30, 10.0)
    lo = np.searchsorted(t, grid - 30); mid = np.searchsorted(t, grid); hi = np.searchsorted(t, grid + 30)
    cand = []
    for g, a0, m0, b0 in zip(grid, lo, mid, hi):
        if m0 - a0 < 10 or b0 - m0 < 10:
            continue
        di = np.median(I[m0:b0]) - np.median(I[a0:m0])
        if abs(di) > 0.3:
            cand.append((g, di, a0, m0, b0))
    groups = []
    for x in cand:
        if groups and c == groups[-1][0] and x[0] - groups[-1][1][-1][0] < 300:
            groups[-1][1].append(x)
        else:
            groups.append((c, [x]))
    for _, gl in groups:
        g, di, a0, m0, b0 = max(gl, key=lambda x: abs(x[1]))
        after = (t > g + 600) & (t < g + 1200)
        steps.append({"ch": c, "t_utc": g, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(g)), "chile": time.strftime("%H:%M", time.gmtime(g - 3 * 3600)),
                      "dI_uA": float(di), "I_before": float(np.median(I[a0:m0])), "I_after": float(np.median(I[m0:b0])),
                      "V_before": float(np.median(V[a0:m0])), "V_after": float(np.median(V[m0:b0])),
                      "I_10_20min_after": float(np.median(I[after])) if after.any() else None})
# IMon contra aX por bloque
corr = {}
for run, fn in (("run010", f"{OD}/temp_corr_blocks_run010.json"), ("run008", f"{OD}/temp_corr_blocks_run008.json")):
    bl = sorted(json.load(open(fn)).values(), key=lambda b: b["t0_utc"])
    x = np.log([b["aX"] for b in bl])
    corr[run] = {}
    for c in range(8):
        y = np.array([b["IMon_uA"][c] for b in bl], float)
        ok = np.isfinite(y)
        if ok.sum() < 4:
            continue
        A = np.vstack([np.ones(ok.sum()), x[ok] - x[ok].mean()]).T
        coef, res, *_ = np.linalg.lstsq(A, np.log(y[ok]), rcond=None)
        resid = np.log(y[ok]) - A @ coef
        s = math.sqrt((resid ** 2).sum() / (ok.sum() - 2))
        eB = s / math.sqrt(((x[ok] - x[ok].mean()) ** 2).sum())
        r = float(np.corrcoef(x[ok], np.log(y[ok]))[0, 1])
        corr[run][c] = {"dlnI_dlnaX": float(coef[1]), "err": eB, "pearson": r, "resid_sd": s, "IMon_first": float(y[ok][0]), "IMon_last": float(y[ok][-1])}
json.dump({"t0_utc": T0, "t1_utc": T1, "steps": steps, "summary": summary, "IMon_vs_aX": corr}, open(f"{OD}/imon_steps.json", "w"), indent=1, default=float)
for s in steps:
    print(s)
for run in corr:
    for c, v in corr[run].items():
        print(run, f"CH{c:02d}", {k: round(x, 3) for k, x in v.items()})
