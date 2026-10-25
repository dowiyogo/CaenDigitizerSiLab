#!/usr/bin/python3 -I
# stage2.py — PROC-1 etapa 2 (solo lectura): accidentales y eficiencia de hod2018 con la definición única de impacto (D2) y el
# estimador compartido (daq/tools/characterization/xy_estimator.py, D3), ch3 fuera (D6).
#  - Impactos: work/hd64_<dataset>_hits.bin (hitdef_extract.C, nch 64). G: mayo = G reajustada con la definición única
#    (work/hd64_may_v8_peaks.csv, semilla de ch9 = su G de c1; A3) para todos los runs de mayo; octubre = c1 (mediana para los
#    canales sin G). Control: G del archivo de mayo (pe_calibration_v8.txt).
#  - Conjuntos: mayo v8 (tiempo bueno y todo), v7, v9, v10 (tiempo bueno), octubre run 8 y noche del 2026-10-19.
#  - Aceptancia (PROC-2, R1): z del CAD (geometría versión 2), 1e8 rayos, cos³θ: f_2018 de work/geom_sens_CAD_1e8.json
#    (geom_sensitivity.py); syst de aceptancia = error del MC ⊕ 2.7 % (placas z ± 2 mm) ⊕ 2.2 % (cota de posición de 10 mm), un número.
#  - g (PROC-2, R2) = P(X∧Y | cruza el cuadrado) = 0.980 (promedio de los dos órdenes de bloque del MC de barras con el CAD,
#    work/bar_mc_g.json), syst 0.03 (diferencia con Meiga, 0.951). ε_XY = q / (f_2018 · g) con cuatro incertidumbres: stat,
#    estimador (mezclados contra fuera de tiempo), aceptancia (R1) y g. R3: ε_X, ε_Y corregidas por plano ya no se calculan;
#    ε(X|Y) y ε(Y|X) quedan como "sin corrección geométrica".
#  - Convención publicada (p0 + 3.5 G del archivo de mayo, definición de mayo aM), las 8073 del manuscrito (work/msxy_may8.bin,
#    manuscript_xy.C), estabilidad, exclusión de ch2/9/15 y la tabla tab:counts.
# Salida: work/stage2.json, tables/eff_vs_threshold.tex, tables/counts.tex, work/fig_eff.json (para proc_plot_eff.C).
# Uso: stage2.py <dir del paquete>
import json, math, os, statistics, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import bench_mc as bm  # noqa: E402
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
_cad = json.load(open(W + "/geom_sens_CAD_1e8.json"))["cos3"]["nominal"]
F18, F18_ERR, F19, FBOTH, FBOTH_ERR = _cad["f"], _cad["f_err"], _cad["f_2019"], _cad["f_both"], _cad["f_both_err"]
PLATE_Z_SYST, POS_SYST = 0.027, 0.022  # R1 (René): placas z ± 2 mm y cota de posición de 10 mm
ACC_REL = math.sqrt((F18_ERR / F18) ** 2 + PLATE_Z_SYST ** 2 + POS_SYST ** 2)
G_XY, G_SYST = 0.980, 0.03  # R2 (René)
TS = list(range(3, 13))
XCH = [c for c in range(16) if c != 3]
YCH = list(range(16, 32))
NOISY = [2, 9, 15]
S = {"stage": 2, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "f_2018": F18, "f_2018_err": F18_ERR,
     "acceptance_rel_syst": ACC_REL}


# ---------------------------------------------------------------- G
def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


G_MAY, G_MAYFILE, P0_MAYFILE, G_C1 = g_refit(), mc.g_may(), mc.p0_may(), mc.g_c1()
fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
S["G_may_refit"] = G_MAY.tolist()
DS = {"may8_good": ("v8", "mayo v8, tiempo bueno", True), "may8_all": ("v8", "mayo v8, todo", False), "may7_good": ("v7", "mayo v7, tiempo bueno", True),
      "may9_good": ("v9", "mayo v9, tiempo bueno", True), "may10_good": ("v10", "mayo v10, tiempo bueno", True),
      "oct8": ("oct8", "octubre run 8", False), "night": ("night", "octubre, noche 2026-10-19", False)}


def load(key):
    run, lab, good = DS[key]
    f = f"{W}/hd64_may_{run}_hits.bin" if run.startswith("v") else f"{W}/hd64_{run}_hits.bin"
    e = xy.load_hits(f, 64)
    if good:
        e = e[mc.good_time(e["t"], run)]
    G = fill(G_MAY) if run.startswith("v") else fill(G_C1)
    return e, np.concatenate([G, np.full(32, np.nan)])


def frac(e, G, T, xs=XCH, ys=YCH):
    Xs, Xc = xy.plane_hits(e, G, T, xs)
    Ys, Yc = xy.plane_hits(e, G, T, ys)
    return xy.xy_fraction(Xs, Ys, Xc, Yc)


# factores de barra
BG = json.load(open(W + "/bar_mc_g.json"))
old = next(r for r in BG["runs"] if r["geometry"] == "CAD")  # R1: z del CAD
def avg(k):
    a, b = old["X arriba"][k], old["Y arriba"][k]
    return {"value": 0.5 * (a["value"] + b["value"]), "stat": max(a["err"], b["err"]), "syst_top": 0.5 * abs(a["value"] - b["value"]),
            "X_top": a["value"], "Y_top": b["value"]}
GF = {k: avg(k) for k in ("g", "g_XY", "g_YX", "f_XY")}
S["geometry_factors"] = GF
S["geometry_factors_note"] = ("bar_mc.py, 1e8 rayos, semilla 42, placas del CAD; g usado = 0.980 ± 0.03 (R2; MC: %.4f, Meiga 2026-10-24_cosmic_telescope_v2b: 0.951)"
                              % GF["g"]["value"])
S["g_used"] = {"value": G_XY, "syst": G_SYST, "mc_average": GF["g"]["value"], "meiga": 0.951}


def enrich(r):
    """eficiencias: ε_geo = q / f_2018; ε_XY = q / (f_2018 · g) con stat, estimador, aceptancia (R1) y g (R2)."""
    q, qe, qc = r["f_mix"], r["f_mix_err"], r["f_ctl"]
    d = {k: r.get(k) for k in ("N", "N_XY", "acc_mix", "acc_ctl", "aX", "aY", "rX", "rY", "f_mix", "f_mix_err", "f_ctl", "f_ctl_err",
                              "e_X|Y_mix", "e_X|Y_mix_err", "e_Y|X_mix", "e_Y|X_mix_err", "e_X|Y_ctl", "e_Y|X_ctl")}
    d["syst_estimator"] = abs(q - qc)
    eg = q / F18
    d["eps_geo"] = {"value": eg, "stat": qe / F18, "syst_estimator": abs(q - qc) / F18, "syst_acceptance": eg * ACC_REL}
    g = G_XY
    exy = q / (F18 * g)
    d["eps_XY"] = {"value": exy, "stat": qe / (F18 * g), "syst_estimator": abs(q - qc) / (F18 * g), "syst_acceptance": exy * ACC_REL,
                   "syst_g": exy * G_SYST / g}
    return d


# ---------------------------------------------------------------- 2.1-2.3 por conjunto y umbral
DATA = {}
for key in DS:
    e, G = load(key)
    DATA[key] = (e, G)
    S.setdefault("datasets", {})[key] = {"label": DS[key][1], "events": int(len(e)), "per_T": {str(T): enrich(frac(e, G, T)) for T in TS}}
    print(key, len(e), "q(5) %.4f ± %.4f" % (S["datasets"][key]["per_T"]["5"]["f_mix"], S["datasets"][key]["per_T"]["5"]["f_mix_err"]))

# control con la G del archivo de mayo (A3)
e8, G8 = DATA["may8_good"]
Gfile = np.concatenate([fill(G_MAYFILE), np.full(32, np.nan)])
S["G_check"] = {str(T): {"refit": enrich(frac(e8, G8, T)), "may_file": enrich(frac(e8, Gfile, T))} for T in (5, 10)}
for T in (5, 10):  # sin ch9, para ver que la diferencia viene de ch9
    xs9 = [c for c in XCH if c != 9]
    S["G_check"][str(T)]["refit_no_ch9"] = enrich(frac(e8, G8, T, xs9))
    S["G_check"][str(T)]["may_file_no_ch9"] = enrich(frac(e8, Gfile, T, xs9))

# ---------------------------------------------------------------- 2.4 convención publicada
thr_pub = P0_MAYFILE + 3.5 * fill(G_MAYFILE)
ee = e8
Xs, Xc = (ee["aM"][:, XCH] > thr_pub[XCH]).any(1), (ee["aMc"][:, XCH] > thr_pub[XCH]).any(1)
Ys, Yc = (ee["aM"][:, YCH] > thr_pub[YCH]).any(1), (ee["aMc"][:, YCH] > thr_pub[YCH]).any(1)
rp = xy.xy_fraction(Xs, Ys, Xc, Yc)
T_eq = (P0_MAYFILE + 3.5 * fill(G_MAYFILE)) / fill(G_MAY)
S["published_convention"] = {"definition": "definición de mayo (mínimo de la ventana de 10 muestras, pedestal local [k−12, k−4]), umbral p0 + 3.5 G de "
                                           "pe_calibration_v8.txt, v8 tiempo bueno, ch3 fuera",
                             "q_mix": rp["f_mix"], "q_mix_err": rp["f_mix_err"], "q_ctl": rp["f_ctl"], "q_ctl_err": rp["f_ctl_err"],
                             "eps_f2018_mix": rp["f_mix"] / F18, "eps_f2018_ctl": rp["f_ctl"] / F18, "eps_0.0360_mix": rp["f_mix"] / 0.0360,
                             "eps_0.0360_ctl": rp["f_ctl"] / 0.0360, "stat_0.0360": rp["f_mix_err"] / 0.0360, "stat_f2018": rp["f_mix_err"] / F18,
                             "eps_XY_mix": rp["f_mix"] / (F18 * G_XY), "stat_eps_XY": rp["f_mix_err"] / (F18 * G_XY), "f_2018": F18,
                             "T_equivalent_per_channel": {str(c): float(T_eq[c]) for c in XCH + YCH},
                             "T_equivalent_median": float(np.median(T_eq[XCH + YCH])), "handoff": "0.88 ± 0.04 ± 0.11 (v7-v10)"}

# ---------------------------------------------------------------- 2.5 las 8073 del manuscrito
d = np.fromfile(W + "/msxy_may8.bin", dtype=np.dtype([("t", "<f8"), ("b", "u1")]))
bits = lambda b, k: (b >> k) & 1 == 1
MS = {}
for lab, sel in (("run B completo (como el manuscrito)", slice(None)), ("run B, tiempo bueno", mc.good_time(d["t"] * 1e-6 + 1778886208.8, "v8"))):
    b = d["b"][sel]
    X, Y, Xc, Yc, X9, Y9, X9c, Y9c = (bits(b, k) for k in range(8))
    r = xy.xy_fraction(X, Y, Xc, Yc)
    N = len(b)
    one_sided = int((X[:-1] & Y[1:]).sum())
    r19 = xy.xy_fraction(X9, Y9, X9c, Y9c)
    both = int((X & Y & X9 & Y9).sum())
    xy18, xy19 = X & Y, X9 & Y9
    both_mix = 0.5 * (int((xy18[:-1] & xy19[1:]).sum()) + int((xy19[:-1] & xy18[1:]).sum())) * N / (N - 1)
    xy18c, xy19c = (X & Yc) | (Xc & Y) | (Xc & Yc), (X9 & Y9c) | (X9c & Y9) | (X9c & Y9c)
    MS[lab] = {"N": N, "N_X_S": int(X.sum()), "N_X_C": int(Xc.sum()), "N_Y_S": int(Y.sum()), "N_Y_C": int(Yc.sum()), "N_XY": int((X & Y).sum()),
               "mixed_one_sided_as_old_code": one_sided, "acc_mix": r["acc_mix"], "acc_ctl": r["acc_ctl"], "q_mix": r["f_mix"], "q_mix_err": r["f_mix_err"],
               "q_ctl": r["f_ctl"], "q_ctl_err": r["f_ctl_err"], "true_XY_mix": r["f_mix"] * N, "true_XY_ctl": r["f_ctl"] * N,
               "accidental_in_NXY_mix": int((X & Y).sum()) - r["f_mix"] * N, "accidental_in_NXY_ctl": int((X & Y).sum()) - r["f_ctl"] * N,
               "err_true_mix": r["f_mix_err"] * N, "err_true_ctl": r["f_ctl_err"] * N,
               "hod2019": {"N_XY": int(xy19.sum()), "acc_mix": r19["acc_mix"], "acc_ctl": r19["acc_ctl"], "q_mix": r19["f_mix"], "q_mix_err": r19["f_mix_err"],
                           "q_ctl": r19["f_ctl"], "q_ctl_err": r19["f_ctl_err"], "aX": r19.get("aX"), "aY": r19.get("aY")},
               "two_station": {"N": both, "acc_mix_expected": both_mix, "acc_offtime_expected": int((xy18 & xy19c).sum()) + int((xy18c & xy19).sum()),
                               "acc_product_expected": N * xy18.mean() * xy19.mean()}}
S["manuscript_8073"] = MS

# ---------------------------------------------------------------- 2.6 estabilidad
stab = {}
for T in (5, 10):
    pts = []
    for key in ("may7_good", "may8_good", "may9_good", "may10_good"):
        r = S["datasets"][key]["per_T"][str(T)]
        pts.append((key, r["f_mix"], r["f_mix_err"]))
    t8 = e8["t"]
    blocks = np.floor((t8 - t8.min()) / (6 * 3600)).astype(int)
    bl = []
    for k in np.unique(blocks):
        sel = blocks == k
        if sel.sum() < 2000:
            continue
        r = frac(e8[sel], G8, T)
        bl.append((int(k), int(sel.sum()), r["f_mix"], r["f_mix_err"]))
    def chi2(v):
        w = [1 / e ** 2 for _, *_, e in v] if False else [1 / x[-1] ** 2 for x in v]
        m = sum(wi * x[-2] for wi, x in zip(w, v)) / sum(w)
        return m, sum(((x[-2] - m) / x[-1]) ** 2 for x in v), len(v) - 1
    stab[str(T)] = {"runs": pts, "runs_chi2": chi2(pts), "v8_6h_blocks": bl, "blocks_chi2": chi2(bl)}
S["stability"] = stab

# ---------------------------------------------------------------- 2.7 sin ch2, ch9, ch15
S["exclude_noisy"] = {str(T): {"with": enrich(frac(e8, G8, T)), "without_2_9_15": enrich(frac(e8, G8, T, [c for c in XCH if c not in NOISY]))}
                      for T in (5, 8, 10, 12)}

# ---------------------------------------------------------------- 2.8 tab:counts (MC: 1e8 rayos, CAD, cos³θ)
rng = np.random.default_rng(42)
cnt = {"trig": 0, "h18": 0, "h19": 0, "both": 0}
for _ in range(50):
    r = bm.generate(rng, 2_000_000, "cos3")
    t = bm.trigger(r, bm.GEOMS["CAD"])
    a, b = bm.in_square(r, 0.0) & t, bm.in_square(r, 98.4) & t
    cnt["trig"] += int(t.sum()); cnt["h18"] += int(a.sum()); cnt["h19"] += int(b.sum()); cnt["both"] += int((a & b).sum())
fr = {k: cnt[k] / cnt["trig"] for k in ("h18", "h19", "both")}
fe = {k: math.sqrt(fr[k] * (1 - fr[k]) / cnt["trig"]) for k in fr}
S["mc_fractions_1e8"] = {"counts": cnt, "f": fr, "err": fe}
NB = MS["run B completo (como el manuscrito)"]["N"]
S["counts_table"] = {"N_triggers": NB, "expected": {k: fr[k] * NB for k in fr}, "expected_err": {k: fe[k] * NB for k in fe}}

# ---------------------------------------------------------------- tablas y figura
P = lambda v, e: f"{v:.4f} $\\pm$ {e:.4f}"
L = ["% PROC-1/PROC-2 etapa 2: fracción X·Y por trigger con accidentales restados (q) y eficiencia. Definición única de impacto, ch3 fuera.",
     f"% G: mayo = reajustada con la definición única (A3); octubre = c1. f_2018 = {F18:.5f} ± {F18_ERR:.5f} (CAD, R1), g = {G_XY} (R2).",
     f"% ε_XY = q_mix/(f_2018 g) ± stat ± estimador ± aceptancia ({100 * ACC_REL:.1f} %) ± g ({100 * G_SYST / G_XY:.1f} %). "
     f"ε(X|Y), ε(Y|X): sin corrección geométrica. Commit {git}.",
     "\\begin{tabular}{llrrrrr}", "\\hline",
     "Conjunto & $T$ (p.e.) & $q_\\mathrm{mix}$ & $q_\\mathrm{fuera}$ & $\\varepsilon_{XY}$ & $\\varepsilon(X|Y)$ & $\\varepsilon(Y|X)$ \\\\", "\\hline"]
for key in ("may8_good", "oct8"):
    for T in (3, 5, 8, 10, 12):
        r = S["datasets"][key]["per_T"][str(T)]
        e = r["eps_XY"]
        L.append(f"{S['datasets'][key]['label']} & {T} & {P(r['f_mix'], r['f_mix_err'])} & {P(r['f_ctl'], r['f_ctl_err'])} & "
                 f"{e['value']:.3f} $\\pm$ {e['stat']:.3f} $\\pm$ {e['syst_estimator']:.3f} $\\pm$ {e['syst_acceptance']:.3f} $\\pm$ {e['syst_g']:.3f} & "
                 f"{r['e_X|Y_mix']:.3f} $\\pm$ {r['e_X|Y_mix_err']:.3f} & "
                 f"{r['e_Y|X_mix']:.3f} $\\pm$ {r['e_Y|X_mix_err']:.3f} \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/eff_vs_threshold.tex", "w").write("\n".join(L) + "\n")
m = MS["run B completo (como el manuscrito)"]
new5 = S["datasets"]["may8_all"]["per_T"]["5"]
C = ["% PROC-1 etapa 2.8: tab:counts nuevo. Run B = v8 completo (143 814 triggers, se omite el último en el mezclado del código antiguo).",
     "% Esperado: MC geométrico (bench_mc.py, 1e8 rayos, semilla 42, cos^3, placas en las z del CAD, R1), sin eficiencia.",
     "% Observado restado: q·N con el estimador de eventos mezclados (central) y el de la ventana fuera de tiempo.",
     "\\begin{tabular}{lrrrr}", "\\hline", "Condición & Esperado (MC) & Observado & Restado (mezclados) & Restado (fuera de tiempo) \\\\", "\\hline",
     f"X--Y en hod2018, 50 mV, 4000--4400 ns & {S['counts_table']['expected']['h18']:.0f} & {m['N_XY']} & {m['true_XY_mix']:.0f} $\\pm$ {m['err_true_mix']:.0f} & "
     f"{m['true_XY_ctl']:.0f} $\\pm$ {m['err_true_ctl']:.0f} \\\\",
     f"X--Y en hod2018, 5 p.e., ventana de 144 ns & {S['counts_table']['expected']['h18']:.0f} & {new5['N_XY']} & {new5['f_mix'] * new5['N']:.0f} $\\pm$ "
     f"{new5['f_mix_err'] * new5['N']:.0f} & {new5['f_ctl'] * new5['N']:.0f} $\\pm$ {new5['f_ctl_err'] * new5['N']:.0f} \\\\",
     f"X--Y en hod2019, 35 mV & {S['counts_table']['expected']['h19']:.0f} & {m['hod2019']['N_XY']} & azar mezclados {m['hod2019']['acc_mix']:.1f} & azar fuera de tiempo "
     f"{m['hod2019']['acc_ctl']} \\\\",
     f"X--Y en las dos estaciones & {S['counts_table']['expected']['both']:.0f} & {m['two_station']['N']} & azar $N P_{{18}} P_{{19}}$ = "
     f"{m['two_station']['acc_product_expected']:.2f} & azar fuera de tiempo {m['two_station']['acc_offtime_expected']} \\\\",
     "% hod2019: no hay impactos de hod2019 sobre 35 mV en la ventana [6000, 6400] ns (aX = aY = 0); el estimador de mezclados da 26.0 contra 3",
     "% observadas, así que no es una expectativa válida de azar para hod2019. Dos estaciones: el mezclado no sirve con 3 eventos X--Y en hod2019;",
     "% se da el producto N·P(X--Y en hod2018)·P(X--Y en hod2019).",
     "\\hline", "\\end{tabular}"]
open(OUT + "/tables/counts.tex", "w").write("\n".join(C) + "\n")
# figura (PROC-2): eje izquierdo ε_XY = q/(f_2018 g), eje derecho q (escala f_2018·g); banda = syst de aceptancia (R1)
fig = {"T": TS, "f18": F18 * G_XY, "f18_rel_syst": ACC_REL, "series": [], "published": {"T": S["published_convention"]["T_equivalent_median"],
                                                                                         "eps": S["published_convention"]["eps_XY_mix"],
                                                                                         "err": S["published_convention"]["stat_eps_XY"]}}
for key, lab in (("may8_good", "May v8"), ("oct8", "Oct run 8"), ("night", "Oct night 10-19")):
    fig["series"].append({"name": lab, "eps": [S["datasets"][key]["per_T"][str(T)]["eps_XY"]["value"] for T in TS],
                          "stat": [S["datasets"][key]["per_T"][str(T)]["eps_XY"]["stat"] for T in TS],
                          "syst": [S["datasets"][key]["per_T"][str(T)]["eps_XY"]["syst_estimator"] for T in TS]})
json.dump(fig, open(W + "/fig_eff.json", "w"), indent=1)
json.dump(S, open(W + "/stage2.json", "w"), indent=1, default=float)
print("hecho")
