#!/usr/bin/python3 -I
# gain_precision.py — seguimiento 2026-11-15 (15, punto 5): precisión estadística de la G horaria de hod2018 (ajustes de followup 12 B,
# fit_gain_hourly.C fuente "spec", p0 fijo al del run: work/{may8,oct8}_G_spec.csv). Solo lectura.
# Bins: los de diurnal.py (B_diurnal/diurnal_bins.csv): mayo run B desde 5.25 h, bins buenos (todos sus eventos en tiempo bueno); octubre run 8, todos.
# Canal-hora aceptado con el criterio de diurnal.py: estado de MIGRAD 0 o 1, error > 0, G en (0.62, 1.38)·semilla, error relativo < 0.10, y ajuste
# del run del canal con estado 0 o 1; ch3 fuera.
#  - 1 h: σ_G/G por canal y hora [M, error del ajuste]; mediana y percentiles 16-84 sobre los canal-hora; por hora, error estadístico de la mediana
#    de la estación ≈ 1.2533·rms(σ_i/G_i)/√n [I]; y el error de la mediana de diurnal.py (1.2533·dispersión entre canales de G/G_ref/√n) [M].
#  - 6 h y 24 h: no hay ajustes directos en esos bins; se combinan los ajustes horarios aceptados del canal dentro del bloque (media ponderada,
#    σ = 1/√Σσ_h⁻²) [I: ajustes horarios independientes; un ajuste del espectro sumado tendría un error parecido]. Bloques consecutivos de 6 h / 24 h UTC
#    desde el primer bin elegido; solo entran las horas elegidas (se da el número de horas por bloque). Octubre cubre 18 h: un solo bloque «24 h».
#  - Control de la escala √N: error del ajuste del run entero (hour = −1, todas las horas del run) contra la combinación de todas las horas del canal.
# Salida: <out>/gain_precision.csv, gain_precision.json. Uso: gain_precision.py <followup12 dir> <out dir>
import csv, json, math, os, sys

import numpy as np

F12, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
H18 = [c for c in range(32) if c != 3]
CAMP = {"may": "may8", "oct": "oct8"}

bins = {}
for r in csv.DictReader(l for l in open(F12 + "/B_diurnal/diurnal_bins.csv") if not l.startswith("#")):
    if r["campaign"] == "oct" or r["good_bin"] == "1":
        bins.setdefault(r["campaign"], {})[int(r["hour_utc_s"])] = float(r["G_rel_median_err"]) / float(r["G_rel_median"])


def q(x):
    x = np.asarray(x, float)
    return float(np.median(x)), float(np.percentile(x, 16)), float(np.percentile(x, 84)), int(len(x))


res, rows = {}, []
for camp, tag in CAMP.items():
    seeds = [float(v) for v in open(f"{F12}/work/{tag}_seeds.txt").read().split(",")]
    run, hr = {}, {}
    for r in csv.DictReader(l for l in open(f"{F12}/work/{tag}_G_spec.csv") if not l.startswith("#")):
        h, c = int(r["hour_utc_s"]), int(r["ch"])
        v = (float(r["G_mV"]), float(r["G_err_mV"]), int(r["fit_status"]), float(r["chi2"]) / max(int(r["ndf"]), 1))
        if h == -1:
            run[c] = v
        else:
            hr[(h, c)] = v
    ok = lambda h, c: ((h, c) in hr and hr[(h, c)][2] in (0, 1) and hr[(h, c)][1] > 0 and 0.62 * seeds[c] < hr[(h, c)][0] < 1.38 * seeds[c]  # noqa: E731
                       and hr[(h, c)][1] / hr[(h, c)][0] < 0.10 and c in run and run[c][2] in (0, 1) and run[c][1] > 0)
    H = sorted(bins[camp])
    out = {"n_hours": len(H)}
    # 1 h
    rel = [hr[(h, c)][1] / hr[(h, c)][0] for h in H for c in H18 if ok(h, c)]
    nacc = [sum(ok(h, c) for c in H18) for h in H]
    st_stat = [1.2533 * math.sqrt(np.mean([(hr[(h, c)][1] / hr[(h, c)][0]) ** 2 for c in H18 if ok(h, c)])) / math.sqrt(n) for h, n in zip(H, nacc) if n]
    out["1h"] = {"channel_rel_err": q(rel), "accepted_channels_per_hour": q(nacc), "accepted_fraction": sum(nacc) / (len(H) * len(H18)),
                 "station_median_stat_rel_err": q(st_stat), "station_median_err_diurnal_rel": q([bins[camp][h] for h in H])}
    # 6 h y 24 h
    for L in (6, 24):
        blocks = {}
        for h in H:
            blocks.setdefault((h - H[0]) // (L * 3600), []).append(h)
        ch_rel, st_rel, nh = [], [], []
        for hs in blocks.values():
            per_ch = []
            for c in H18:
                e = [hr[(h, c)] for h in hs if ok(h, c)]
                if not e:
                    continue
                w = np.array([1 / x[1] ** 2 for x in e])
                g = float((w * np.array([x[0] for x in e])).sum() / w.sum())
                per_ch.append(1 / math.sqrt(w.sum()) / g)
            ch_rel += per_ch
            nh.append(len(hs))
            if per_ch:
                st_rel.append(1.2533 * math.sqrt(np.mean(np.square(per_ch))) / math.sqrt(len(per_ch)))
        out[f"{L}h"] = {"n_blocks": len(blocks), "hours_per_block": q(nh), "channel_rel_err": q(ch_rel), "station_median_stat_rel_err": q(st_rel)}
    # control √N: run entero contra combinación de todas las horas del run (todas las horas del archivo, no solo las elegidas)
    allh = sorted({h for h, _ in hr})
    ctl = []
    for c in H18:
        e = [hr[(h, c)] for h in allh if ok(h, c)]
        if e and c in run and run[c][2] in (0, 1):
            ctl.append((run[c][1] / run[c][0]) / (1 / math.sqrt(sum(1 / x[1] ** 2 for x in e)) / run[c][0]))
    out["check_run_fit_err_over_combined_hourly"] = q(ctl)
    # control empírico: χ²/ndf de los ajustes horarios aceptados, y pull de horas consecutivas elegidas (Δ = G_h+1 − G_h, pull = Δ/√(σ_h² + σ_h+1²));
    # la dispersión 1.4826·MAD del pull > 1 indica errores subestimados (o variación real de G en 1 h)
    out["hourly_chi2_ndf"] = q([hr[(h, c)][3] for h in H for c in H18 if ok(h, c)])
    pulls = [(hr[(b, c)][0] - hr[(a, c)][0]) / math.hypot(hr[(a, c)][1], hr[(b, c)][1]) for a, b in zip(H[:-1], H[1:]) if b - a == 3600
             for c in H18 if ok(a, c) and ok(b, c)]
    pm = float(np.median(pulls))
    out["consecutive_hour_pull"] = {"n": len(pulls), "median": pm, "spread_1p4826MAD": 1.4826 * float(np.median(np.abs(np.array(pulls) - pm))),
                                    "rms": float(np.sqrt(np.mean(np.square(pulls))))}
    res[camp] = out
    for k in ("1h", "6h", "24h"):
        o = out[k]
        rows.append([camp, k, o.get("n_blocks", len(H)), f"{o['hours_per_block'][0]:.0f}" if "hours_per_block" in o else "1",
                     f"{100 * o['channel_rel_err'][0]:.3f}", f"{100 * o['channel_rel_err'][1]:.3f}", f"{100 * o['channel_rel_err'][2]:.3f}", o["channel_rel_err"][3],
                     f"{100 * o['station_median_stat_rel_err'][0]:.3f}", f"{100 * o['station_median_stat_rel_err'][1]:.3f}", f"{100 * o['station_median_stat_rel_err'][2]:.3f}",
                     f"{100 * o['station_median_err_diurnal_rel'][0]:.3f}" if k == "1h" else ""])
with open(f"{OUT}/gain_precision.csv", "w") as f:
    f.write("# gain_precision.py (seguimiento 15, 2026-11-15): relative statistical error of the hod2018 gain G, in percent; 1 h = followup 12 B hourly fits "
            "(fit error); 6 h, 24 h = inverse-variance combination of the accepted hourly fits of the channel in the block [I]; station = 1.2533*rms(sigma_i/G_i)/sqrt(n); "
            "station_diurnal = followup-12 median error from the channel spread of G/G_ref (1 h only); May = run B from 5.25 h, good bins; October run 8 (18 h: one block)\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "bin", "n_bins", "hours_per_bin_median", "channel_rel_err_median_pct", "channel_rel_err_p16_pct", "channel_rel_err_p84_pct",
                "n_channel_bins", "station_stat_rel_err_median_pct", "station_stat_p16_pct", "station_stat_p84_pct", "station_diurnal_rel_err_median_pct"])
    w.writerows(rows)
json.dump(res, open(f"{OUT}/gain_precision.json", "w"), indent=1)
for r in rows:
    print(" ".join(str(x) for x in r))
for camp in res:
    print(camp, "1 h: canales aceptados por hora (mediana, p16, p84, n)", res[camp]["1h"]["accepted_channels_per_hour"], "fracción", f"{res[camp]['1h']['accepted_fraction']:.3f}",
          "| control error del run / combinación horaria", res[camp]["check_run_fit_err_over_combined_hourly"],
          "| χ²/ndf horario", res[camp]["hourly_chi2_ndf"], "| pull de horas consecutivas", res[camp]["consecutive_hour_pull"])
