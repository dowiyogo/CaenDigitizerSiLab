#!/usr/bin/python3 -I
# rm_tables.py — seguimiento 2026-11-19 (matriz de respuesta): tablas Markdown del README a partir de response_matrix*.json y core_tail.json (sin copiar a mano).
# Uso: rm_tables.py <out dir>  > <out>/logs/tables.md
import json, math, sys

import numpy as np

OUT = sys.argv[1]
VARS = (("ls", ""), ("gnorm", "_gnorm"), ("tcut", "_tcut"), ("gnorm_tcut", "_gnorm_tcut"))
R = {v: json.load(open(f"{OUT}/response_matrix{s}.json")) for v, s in VARS}
C = json.load(open(f"{OUT}/core_tail.json"))
print("### Linealidad y migración hacia arriba por capa\n")
print("| Campaña | Capa | Variante | media/n, n = 1 | 5 | 10 | 20 | 40 | RMS(A), n = 1 | 10 | 40 | P(A ≥ n+1), n = 5 | n = 10 |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for camp in ("may", "oct"):
    for L in ("XA", "XB", "YA", "YB"):
        for v, _ in VARS:
            x = R[v][camp]["layers"][L]
            m, r, p = x["mean_over_n"], x["rms"], x["P_ge_n_plus_1"]
            print(f"| {camp} | {L} | {v} | " + " | ".join(f"{m[k]:.4f}" for k in ("1", "5", "10", "20", "40")) + " | "
                  + " | ".join(f"{r[k]:.3f}" for k in ("1", "10", "40")) + f" | {p['5']['p']:.4f} ± {p['5']['err']:.4f} | {p['10']['p']:.4f} ± {p['10']['err']:.4f} |")
print("\n### Inyecciones por n y capa\n")
print("| Campaña | XA | XB | YA | YB |\n|---|---|---|---|---|")
for camp in ("may", "oct"):
    print(f"| {camp} | " + " | ".join(f"{R['ls'][camp]['layers'][L]['injections_per_n_min']:.0f}" for L in ("XA", "XB", "YA", "YB")) + " |")
print("\n### n = 0 (ruido solo): ventana signal contra ventana de control [99,108), mismos segmentos\n")
print("| Campaña | Capa | media(A) | P(A ≥ 0.5) signal / control | P(A ≥ 1.5) | P(A ≥ 3.0) | P(A ≥ 5.0) |\n|---|---|---|---|---|---|---|")
for camp in ("may", "oct"):
    for L in ("XA", "XB", "YA", "YB"):
        x = R["ls"][camp]["layers"][L]
        n0 = x["noise_n0"]
        print(f"| {camp} | {L} | {x['mean_n0']:.4f} | " + " | ".join(f"{n0[k]['P_signal_ge']:.4f} / {n0[k]['P_control_ge']:.4f}" for k in ("0.5", "1.5", "3.0", "5.0")) + " |")
print("\n### Núcleo gaussiano ([0.75 n, 1.25 n]) y cola baja\n")
print("| Campaña | Capa | Variante | núcleo media/n, n = 10 | n = 40 | σ/n, n = 10 | n = 40 | σ/√n, n = 10 | n = 40 | P(A < n/2), n = 10 | n = 40 |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for camp in ("may", "oct"):
    for L in ("XA", "XB", "YA", "YB"):
        for v in ("ls", "tcut", "gnorm", "gnorm_tcut"):
            c = C[camp][v][L]
            print(f"| {camp} | {L} | {v} | {c['10']['core_mean_over_n']:.4f} | {c['40']['core_mean_over_n']:.4f} | {c['10']['core_sigma_over_n']:.4f} | "
                  f"{c['40']['core_sigma_over_n']:.4f} | {c['10']['core_sigma_over_sqrt_n']:.3f} | {c['40']['core_sigma_over_sqrt_n']:.3f} | {c['10']['P_A_lt_half_n']:.4f} | "
                  f"{c['40']['P_A_lt_half_n']:.4f} |")
print("\n### Comprobación de 1 p.e. (medianas sobre canales; tirones = canales con |tirón| < 2 / canales con los dos ajustes)\n")
print("| Variante | Campaña | n = 1 pos. (g+e) | oscuro n = 0 pos. (g+e) | n = 1 σ | oscuro σ | tirones pos. | tirones σ | n = 1 (g) | n = 2 (g) | separación | pendiente 20-60 | tirones pendiente vs G | error rel. de G |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for v, _ in VARS:
    for camp in ("may", "oct"):
        k = R[v][camp]["check_1pe"]
        fin = lambda a: int(np.isfinite(np.array(a, float)).sum())  # noqa: E731
        print(f"| {v} | {camp} | {k['median_inj1_ge_pos']:.4f} | {k['median_dark_ge_pos']:.4f} | {k['median_inj1_ge_sigma']:.4f} | {k['median_dark_ge_sigma']:.4f} | "
              f"{k['n_abs_pull_pos_lt2']}/{fin(k['pull_pos_inj_vs_dark'])} | {k['n_abs_pull_sigma_lt2']}/{fin(k['pull_sigma_inj_vs_dark'])} | {k['median_inj1_g_pos']:.4f} | "
              f"{k['median_inj2_g_pos']:.4f} | {k['median_spacing']:.4f} | {k['median_slope_20_60']:.4f} | {k['n_abs_pull_slope_lt2']}/{fin(k['pull_slope_vs_G'])} | "
              f"{k['median_G_err_rel']:.4f} |")
