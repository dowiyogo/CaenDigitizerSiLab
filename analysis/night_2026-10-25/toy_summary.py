#!/usr/bin/python3 -I
# toy_summary.py — XPLUS-4 etapa 4: tablas del toy (toy_main.json, toy_delta.json) y el δ = a^oot/a − 1 que reproduce la diferencia
# observada entre estimadores (mezclados − fuera de tiempo), por interpolación lineal en el barrido de δ.
# Uso: toy_summary.py <dir del toy>
import json, os, sys

import numpy as np

D = sys.argv[1]
M = json.load(open(os.path.join(D, "toy_main.json")))
S = json.load(open(os.path.join(D, "toy_delta.json")))
FX = {"v8": (0.027495083052747792, 0.03076897181513901, 0.0022618261790133573), "run8": (0.02386720834654745, 0.026076807182748702, 0.0016840614637531606)}
FG = M["f_2018"] * M["g"]
lines = ["| conjunto | configuración | q_mu | q_T | sesgo mix vs q_T | sesgo ctl vs q_T | sd/σ mix | sd/σ ctl | cobertura 1σ mix | cobertura 1σ ctl | mix − ctl (media ± error; sd) | sesgo mix vs q_mu |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in M["results"]:
    c, p, a, b, d = r["config"], r["params"], r["mix"], r["ctl"], r["mix_minus_ctl"]
    lines.append(f"| {c['ds']} | {c['name']} | {p['q_mu']:.5f} | {p['q_T']:.5f} | {a['bias_vs_qT']:+.2e} ± {a['bias_err']:.1e} | {b['bias_vs_qT']:+.2e} ± {b['bias_err']:.1e} | "
                 f"{a['sd_over_sigma']:.3f} | {b['sd_over_sigma']:.3f} | {a['coverage_1sigma']:.3f} | {b['coverage_1sigma']:.3f} | "
                 f"{d['mean']:+.2e} ± {d['err']:.1e}; {d['sd']:.2e} | {a['bias_vs_qmu']:+.2e} |")
open(os.path.join(D, "toy_table.md"), "w").write("\n".join(lines) + "\n")
out = {"delta": {}}
dl = ["| conjunto | δ | q_T | mix − ctl (media ± error) | sd(mix − ctl) | sesgo mix | sesgo ctl |", "|---|---|---|---|---|---|---|"]
for ds in FX:
    rows = sorted([r for r in S["results"] if r["config"]["ds"] == ds], key=lambda r: r["config"]["delta"])
    x = np.array([r["config"]["delta"] for r in rows])
    y = np.array([r["mix_minus_ctl"]["mean"] for r in rows])
    for r in rows:
        dl.append(f"| {ds} | {r['config']['delta']:+.2f} | {r['params']['q_T']:.5f} | {r['mix_minus_ctl']['mean']:+.3e} ± {r['mix_minus_ctl']['err']:.1e} | "
                  f"{r['mix_minus_ctl']['sd']:.3e} | {r['mix']['bias_vs_qT']:+.2e} | {r['ctl']['bias_vs_qT']:+.2e} |")
    obs = FX[ds][0] - FX[ds][1]
    sd0 = next(r for r in rows if r["config"]["delta"] == 0.0)["mix_minus_ctl"]["sd"]
    slope = np.polyfit(x, y, 1)
    d_obs = (obs - slope[1]) / slope[0]
    out["delta"][ds] = {"observed_mix_minus_ctl": obs, "toy_sd_mix_minus_ctl_at_delta0": sd0, "observed_over_sd": obs / sd0,
                        "linear_fit_slope_per_unit_delta": slope[0], "linear_fit_intercept": slope[1], "delta_reproducing_observed": d_obs,
                        "eps_estimator_term_observed": abs(obs) / FG, "eps_estimator_term_statistical_sd": sd0 / FG}
open(os.path.join(D, "toy_delta_table.md"), "w").write("\n".join(dl) + "\n")
base = {ds: [r for r in M["results"] if r["config"]["ds"] == ds and r["config"]["name"] == "base"] for ds in FX}
out["base_max_abs_bias_vs_qT"] = {ds: {t: max(abs(r[t]["bias_vs_qT"]) for r in base[ds]) for t in ("mix", "ctl")} for ds in FX}
val = [r for r in M["results"] if r["config"]["name"] == "validacion_a0_r0"]
out["validation"] = [{"ds": r["config"]["ds"], "q": r["config"]["q"], "bias_mix_over_err": r["mix"]["bias_vs_qT"] / r["mix"]["bias_err"],
                      "bias_ctl_over_err": r["ctl"]["bias_vs_qT"] / r["ctl"]["bias_err"], "sd_over_sigma_mix": r["mix"]["sd_over_sigma"],
                      "sd_over_sigma_ctl": r["ctl"]["sd_over_sigma"]} for r in val]
json.dump(out, open(os.path.join(D, "toy_summary.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
