#!/usr/bin/python3 -I
# h19_xy.py — PROC-2 etapa 3 (solo lectura): hod2019, qué son los impactos X e Y. Mayo v8 completo (run B), definición única
# (work/hd64_may_v8_hits.bin; kCfg_2026_05: X = ch48-63, Y = ch32-47); octubre run 8 (r4: X = ch32-47, Y = ch48-63) para las
# amplitudes de ch47 y ch39.
#  1. A 3 y 35 mV: rX, rY, rX·rY·N y q·N (mezclados y fuera de tiempo) con xy_estimator.xy_fraction; déficit de X·Y en el mismo evento
#     respecto del producto: N_XY observado contra el esperado de eventos mezclados (X del evento i con Y del i+1, que es el producto
#     de las probabilidades en tiempo), con su significancia de Poisson.
#  2. Respuesta del otro plano: P(algún Y > 3 mV | algún X > 3 mV) en el mismo trigger, en tiempo y en la ventana de control, y al
#     revés; multiplicidad de canales de hod2019 por trigger a 3 mV (1; 2 en el mismo plano; 1 + 1 en los dos planos; ≥ 3), en tiempo
#     menos fuera de tiempo.
#  3. Amplitud de los pulsos en tiempo sobre 3 mV: mediana y percentil 90 por canal (mayo v8) y para ch47 y ch39 en octubre run 8, con
#     error de la mediana por remuestreo (bootstrap, 2000, semilla 42).
# Salida: work/h19_xy.json, tables/h19_xy.tex. Uso: h19_xy.py <dir del paquete>
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": "PROC-2 etapa 3", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
e = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
eo = xy.load_hits(W + "/hd64_oct8_hits.bin", 64)
X, Y = list(range(48, 64)), list(range(32, 48))
N = len(e)
rng = np.random.default_rng(42)

# ---------------------------------------------------------------- 1. déficit X·Y
S["xy"] = {}
for thr in (3, 35):
    xs, ys = (e["aS"][:, X] > thr).any(1), (e["aS"][:, Y] > thr).any(1)
    xc, yc = (e["aC"][:, X] > thr).any(1), (e["aC"][:, Y] > thr).any(1)
    r = xy.xy_fraction(xs, ys, xc, yc)
    exp_ = r["acc_mix"]
    d = {k: r[k] for k in ("N", "N_X", "N_Y", "N_XY", "acc_mix", "acc_ctl", "aX", "aY", "rX", "rY", "f_mix", "f_mix_err", "f_ctl", "f_ctl_err")}
    d.update(rXrY_N=r["rX"] * r["rY"] * N, qN_mix=r["f_mix"] * N, qN_mix_err=r["f_mix_err"] * N, qN_ctl=r["f_ctl"] * N, qN_ctl_err=r["f_ctl_err"] * N,
             product_in_time_N=N * xs.mean() * ys.mean(), deficit=exp_ - r["N_XY"], deficit_z=(exp_ - r["N_XY"]) / math.sqrt(exp_) if exp_ > 0 else None,
             ratio_observed_over_product=r["N_XY"] / exp_ if exp_ > 0 else None)
    # probabilidad de Poisson de observar ≤ N_XY con media exp_
    k, lam = r["N_XY"], exp_
    if k < 50:
        d["p_poisson_le_observed"] = float(sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k + 1)))
    S["xy"][f"{thr}mV"] = d
    print(thr, "mV:", {kk: (round(v, 6) if isinstance(v, float) else v) for kk, v in d.items()})

# ---------------------------------------------------------------- 2. respuesta del otro plano y multiplicidad a 3 mV
thr = 3
hs, hc = e["aS"][:, 32:64] > thr, e["aC"][:, 32:64] > thr  # columnas 0-15 = Y (ch32-47), 16-31 = X (ch48-63)
resp = {}
for win, h in (("in_time", hs), ("control", hc)):
    xh, yh = h[:, 16:].any(1), h[:, :16].any(1)
    resp[win] = {"N_X": int(xh.sum()), "N_Y": int(yh.sum()), "P_Y_given_X": float((xh & yh).sum() / xh.sum()), "P_X_given_Y": float((xh & yh).sum() / yh.sum()),
                 "P_Y_any": float(yh.mean()), "P_X_any": float(xh.mean())}
S["response_3mV"] = resp


def classes(h):
    nx, ny = h[:, 16:].sum(1), h[:, :16].sum(1)
    n = nx + ny
    return {"1": int((n == 1).sum()), "2_same_plane": int(((n == 2) & ((nx == 2) | (ny == 2))).sum()), "1X+1Y": int(((nx == 1) & (ny == 1)).sum()),
            ">=3": int((n >= 3).sum()), ">=3_both_planes": int(((n >= 3) & (nx > 0) & (ny > 0)).sum())}


ci, cc = classes(hs), classes(hc)
S["multiplicity_3mV"] = {"in_time": ci, "control": cc, "in_minus_control": {k: ci[k] - cc[k] for k in ci},
                         "err": {k: math.sqrt(ci[k] + cc[k]) for k in ci}}
print("respuesta:", resp)
print("multiplicidad:", S["multiplicity_3mV"])

# ---------------------------------------------------------------- 3. amplitudes en tiempo sobre 3 mV
def amp_stats(a):
    a = a[a > thr]
    if len(a) < 5:
        return {"n": int(len(a))}
    bs = [np.median(rng.choice(a, len(a))) for _ in range(2000)]
    return {"n": int(len(a)), "median_mV": float(np.median(a)), "median_err_mV": float(np.std(bs)), "p90_mV": float(np.percentile(a, 90))}


S["amplitudes_may8"] = {str(c): amp_stats(e["aS"][:, c]) for c in range(32, 64)}
S["amplitudes_oct8"] = {str(c): amp_stats(eo["aS"][:, c]) for c in (39, 47)}
S["amplitudes_oct8_control"] = {str(c): amp_stats(eo["aC"][:, c]) for c in (39, 47)}
med_may = [v["median_mV"] for v in S["amplitudes_may8"].values() if "median_mV" in v]
p90_may = [v["p90_mV"] for v in S["amplitudes_may8"].values() if "p90_mV" in v]
S["amplitudes_may8_summary"] = {"median_of_medians_mV": float(np.median(med_may)), "range_medians_mV": [float(min(med_may)), float(max(med_may))],
                                "median_of_p90_mV": float(np.median(p90_may)), "range_p90_mV": [float(min(p90_may)), float(max(p90_may))]}
print("mayo:", S["amplitudes_may8_summary"]); print("octubre:", S["amplitudes_oct8"], "control:", S["amplitudes_oct8_control"])

# ---------------------------------------------------------------- tabla
L = ["% PROC-2 etapa 3: hod2019 en mayo v8 (run B completo), definición única; X = ch48-63, Y = ch32-47 (kCfg_2026_05).",
     "% X·Y esperadas del producto = eventos mezclados; significancia del déficit = (esperado − observado)/sqrt(esperado).", f"% Commit {git}.",
     "\\begin{tabular}{lrr}", "\\hline", "Cantidad & 3 mV & 35 mV \\\\", "\\hline"]
a, b = S["xy"]["3mV"], S["xy"]["35mV"]
for name, k, fmt in (("$r_X$", "rX", ".5f"), ("$r_Y$", "rY", ".5f"), ("$r_X r_Y N$", "rXrY_N", ".1f"), ("$N_{XY}$ observadas", "N_XY", "d"),
                     ("$N_{XY}$ del producto (mezclados)", "acc_mix", ".1f"), ("déficit / $\\sqrt{\\mathrm{esperado}}$", "deficit_z", ".1f")):
    L.append(f"{name} & {format(a[k], fmt)} & {format(b[k], fmt)} \\\\")
L.append(f"$qN$ (mezclados) & {a['qN_mix']:.0f} $\\pm$ {a['qN_mix_err']:.0f} & {b['qN_mix']:.1f} $\\pm$ {b['qN_mix_err']:.1f} \\\\")
L.append(f"$qN$ (fuera de tiempo) & {a['qN_ctl']:.0f} $\\pm$ {a['qN_ctl_err']:.0f} & {b['qN_ctl']:.1f} $\\pm$ {b['qN_ctl_err']:.1f} \\\\")
L += ["\\hline", "Multiplicidad a 3 mV (en tiempo $-$ control) & & \\\\"]
for k, v in S["multiplicity_3mV"]["in_minus_control"].items():
    L.append(f"\\quad {k.replace('>=', '$\\geq$ ').replace('_', ' ')} & {v} $\\pm$ {S['multiplicity_3mV']['err'][k]:.0f} & \\\\")
L.append(f"$P(Y \\mid X)$ en tiempo / control (3 mV) & {resp['in_time']['P_Y_given_X']:.4f} / {resp['control']['P_Y_given_X']:.4f} & \\\\")
L.append(f"$P(X \\mid Y)$ en tiempo / control (3 mV) & {resp['in_time']['P_X_given_Y']:.4f} / {resp['control']['P_X_given_Y']:.4f} & \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/h19_xy.tex", "w").write("\n".join(L) + "\n")
json.dump(S, open(W + "/h19_xy.json", "w"), indent=1, default=float)
print("hecho")
