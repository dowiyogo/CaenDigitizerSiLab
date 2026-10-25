#!/usr/bin/python3 -I
# headline.py — PROC-2 etapa 1.3: tabla tables/headline.tex con los números principales del artículo para mayo v8 (tiempo bueno) y
# octubre run 8 a 5 y 10 p.e. (R5): q (mezclados), q (fuera de tiempo), ε_XY = q/(f_2018 g) con sus cuatro incertidumbres (stat,
# estimador, aceptancia R1, g R2), ε(X|Y) y ε(Y|X) sin corrección geométrica, y la fracción de doble impacto (pares vecinos) en X e Y.
# Entradas: work/stage2.json, work/stage5.json. Uso: headline.py <dir del paquete>
import json, subprocess, sys

OUT = sys.argv[1]
W = OUT + "/work"
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
s2 = json.load(open(W + "/stage2.json"))
s5 = json.load(open(W + "/stage5.json"))
P = lambda v, e: f"{v:.4f} $\\pm$ {e:.4f}"
L = ["% PROC-2: números principales. Definición única de impacto, ch3 fuera, accidentales por eventos mezclados (central).",
     f"% ε_XY = q_mix/(f_2018 g), f_2018 = {s2['f_2018']:.5f} (CAD, R1), g = {s2['g_used']['value']} (R2); incertidumbres: stat, estimador "
     f"(|q_mix − q_fuera|), aceptancia ({100 * s2['acceptance_rel_syst']:.1f} %), g ({100 * s2['g_used']['syst'] / s2['g_used']['value']:.1f} %).",
     "% ε(X|Y) = q/rY, ε(Y|X) = q/rX: sin corrección geométrica. Doble impacto: pares vecinos verdaderos por trigger con impacto en el plano",
     f"% (etiqueta: el otro plano > 10 p.e.; ± stat ± etiqueta a 15 p.e.). Commit {git}.",
     "\\begin{tabular}{llrrrrrrr}", "\\hline",
     "Conjunto & $T$ (p.e.) & $q_\\mathrm{mix}$ & $q_\\mathrm{fuera}$ & $\\varepsilon_{XY}$ & $\\varepsilon(X|Y)$ & $\\varepsilon(Y|X)$ & "
     "$f_2^X$ & $f_2^Y$ \\\\", "\\hline"]
for key in ("may8_good", "oct8"):
    for T in (5, 10):
        r = s2["datasets"][key]["per_T"][str(T)]
        e = r["eps_XY"]
        f2 = {p: next(x for x in s5["measured"][key][p] if x["T"] == T) for p in "XY"}
        L.append(f"{s2['datasets'][key]['label']} & {T} & {P(r['f_mix'], r['f_mix_err'])} & {P(r['f_ctl'], r['f_ctl_err'])} & "
                 f"{e['value']:.3f} $\\pm$ {e['stat']:.3f} $\\pm$ {e['syst_estimator']:.3f} $\\pm$ {e['syst_acceptance']:.3f} $\\pm$ {e['syst_g']:.3f} & "
                 f"{r['e_X|Y_mix']:.3f} $\\pm$ {r['e_X|Y_mix_err']:.3f} & {r['e_Y|X_mix']:.3f} $\\pm$ {r['e_Y|X_mix_err']:.3f} & "
                 + " & ".join(f"{f2[p]['f2']:.3f} $\\pm$ {f2[p]['f2_stat']:.3f} $\\pm$ {f2[p]['syst_tag']:.3f}" for p in "XY") + " \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/headline.tex", "w").write("\n".join(L) + "\n")
print("\n".join(L[7:-2]))
