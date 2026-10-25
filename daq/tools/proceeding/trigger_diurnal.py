#!/usr/bin/python3 -I
# trigger_diurnal.py — PROC-2 etapa 5 (solo lectura): tasa del trigger contra la hora del día en v1, v8 y v10.
#  - Tiempos del trigger de work/hd64_may_<run>_hits.bin (UTC = época de mayo + TTT). Tiempo vivo = suma de los ΔT ≤ 20 s (los ΔT
#    mayores son interrupciones y no cuentan); cada intervalo va a la hora (UTC) en que empieza. Tasa por hora = intervalos / tiempo
#    vivo; se usan las horas con ≥ 1800 s vivos.
#  - Plegado en la hora local de Chile (UTC−4 en mayo): media y rms por hora del día (de las horas de todos los días del run).
#  - Ajuste de un coseno a las tasas por hora (sin plegar): r = a + b cos(2π (h − φ)/24 h), mínimos cuadrados pesados; amplitud b ≥ 0
#    y fase φ (hora local del máximo) con errores de la matriz de covarianza.
#  - rms de las tasas por hora relativo a la esperada de Poisson (√N / tiempo vivo).
# Salida: work/trigger_diurnal.json, work/fig_trigger_diurnal.json (proc_plot_diurnal.C). Uso: trigger_diurnal.py <dir del paquete>
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": "PROC-2 etapa 5", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "local_offset_h": -4, "gap_s": 20.0,
     "min_live_s": 1800.0, "runs": {}}
fig = {"runs": []}
for run in ("v1", "v8", "v10"):
    t = np.sort(xy.load_hits(f"{W}/hd64_may_{run}_hits.bin", 64)["t"])
    dt = np.diff(t)
    live = dt <= S["gap_s"]
    hb = np.floor(t[:-1] / 3600).astype(np.int64)
    hours = np.unique(hb)
    rows = []
    for h in hours:
        sel = (hb == h) & live
        lt = float(dt[sel].sum())
        n = int(sel.sum())
        if lt >= S["min_live_s"]:
            rows.append({"hour_utc_s": int(h * 3600), "n": n, "live_s": lt, "rate_Hz": n / lt, "err_Hz": math.sqrt(n) / lt})
    rate = np.array([r["rate_Hz"] for r in rows])
    err = np.array([r["err_Hz"] for r in rows])
    hloc = np.array([((r["hour_utc_s"] / 3600 + S["local_offset_h"]) % 24) + 0.5 for r in rows])  # centro de la hora local
    # ajuste a + c cos(ωh) + s sin(ωh)
    w = 1 / err ** 2
    A = np.vstack([np.ones_like(hloc), np.cos(2 * np.pi * hloc / 24), np.sin(2 * np.pi * hloc / 24)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, c, s = cov @ (A.T @ (w * rate))
    b = math.hypot(c, s)
    phi = (math.degrees(math.atan2(s, c)) / 15.0) % 24  # hora del máximo
    jb = np.array([c / b, s / b])
    eb = math.sqrt(jb @ cov[1:, 1:] @ jb)
    jp = np.array([-s, c]) / b ** 2 * (24 / (2 * math.pi))
    ephi = math.sqrt(jp @ cov[1:, 1:] @ jp)
    chi2_c = float((w * (rate - (A @ np.array([a, c, s]))) ** 2).sum())
    m0 = float((w * rate).sum() / w.sum())
    chi2_0 = float((w * (rate - m0) ** 2).sum())
    fold = []
    for hh in range(24):
        v = rate[(hloc >= hh) & (hloc < hh + 1)]
        fold.append({"hour_local": hh, "n_hours": int(len(v)), "mean_Hz": float(v.mean()) if len(v) else None,
                     "rms_Hz": float(v.std(ddof=1)) if len(v) > 1 else None, "err_mean_Hz": float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None})
    d = {"hours_used": len(rows), "span_h": float((t[-1] - t[0]) / 3600), "rows": rows, "weighted_mean_Hz": m0, "chi2_constant": chi2_0, "ndf_constant": len(rows) - 1,
         "cosine": {"a_Hz": float(a), "amplitude_Hz": b, "amplitude_err_Hz": eb, "rel_amplitude": b / a, "phase_local_h": phi, "phase_err_h": ephi,
                    "chi2": chi2_c, "ndf": len(rows) - 3},
         "rms_over_poisson": float(rate.std(ddof=1) / np.sqrt((err ** 2).mean())), "rms_Hz": float(rate.std(ddof=1)), "poisson_rms_Hz": float(np.sqrt((err ** 2).mean())),
         "folded": fold}
    S["runs"][run] = d
    print(run, "%d horas (%.1f h), media %.4f Hz, χ² constante %.1f/%d; coseno b = %.4f ± %.4f Hz (%.2f %%), máximo a las %.1f ± %.1f h locales, χ² %.1f/%d; rms/Poisson %.2f" % (
        len(rows), d["span_h"], m0, chi2_0, len(rows) - 1, b, eb, 100 * b / a, phi, ephi, chi2_c, len(rows) - 3, d["rms_over_poisson"]))
    t0 = rows[0]["hour_utc_s"]
    fig["runs"].append({"name": run, "t_h": [(r["hour_utc_s"] - t0) / 3600 + 0.5 for r in rows], "rate": rate.tolist(), "err": err.tolist(),
                        "fold_h": [f["hour_local"] + 0.5 for f in fold if f["mean_Hz"] is not None and f["err_mean_Hz"] is not None],
                        "fold_mean": [f["mean_Hz"] for f in fold if f["mean_Hz"] is not None and f["err_mean_Hz"] is not None],
                        "fold_err": [f["err_mean_Hz"] for f in fold if f["mean_Hz"] is not None and f["err_mean_Hz"] is not None],
                        "cos": [float(a), float(c), float(s)], "median_Hz": float(np.median(rate)), "draw_cos": run != "v1"})
json.dump(fig, open(W + "/fig_trigger_diurnal.json", "w"), indent=1)
json.dump(S, open(W + "/trigger_diurnal.json", "w"), indent=1, default=float)
print("hecho")
