#!/usr/bin/python3 -I
# along_bar.py — PROC-2 etapa 2 (solo lectura): eficiencia a lo largo de las barras de hod2018 con la coordenada ortogonal de la
# misma estación como referencia de posición.
#  - Para cada coordenada Y j (orden de mayo validado; octubre r4): el estimador de xy_estimator.xy_fraction con Y restringido a la
#    tira j (impacto en tiempo y en la ventana de control de esa tira) y X = todo el plano, solo la capa A (coordenadas pares,
#    centros −15, −11, …) o solo la capa B (impares): ε(X | Y_j) = q_j / rY_j (con el término rX·rY_j de los mezclados). La tira j
#    fija la posición a lo largo de las barras X: y_j = −15 + 2 j mm (centros del CAD). Lo mismo con los papeles cambiados:
#    ε(Y | X_i) contra x_i. 5 y 10 p.e., ch3 (coordenada X 12, capa A) fuera.
#  - Geometría: P(capa con cruce | la tira de referencia cruzada) del MC de barras (bench_mc, CAD, cos³θ, detección perfecta,
#    promedio de los dos órdenes de bloque); se da ε/ε_MC para separar la forma geométrica.
#  - Ajustes de ε contra la posición: recta y exponencial (mínimos cuadrados pesados; la exponencial sobre ln ε), pendiente por mm
#    con su error y razón entre los extremos del tramo de 33 mm (±16.5 mm) de cada ajuste.
#  - Fuera del cuadrado: R = (rX − q)/q medido contra R0 = (f_X − f_XY)/f_XY del MC con detección perfecta; se da R/R0.
#  - Lado de los MPPC de cada capa: sin registro (mapeo r4: face "unknown"); las tendencias van como capa A y capa B.
# Salida: work/along_bar.json, work/fig_along_bar.json (proc_plot_along_bar.C), tables/eff_along_bar.tex.
# Uso: along_bar.py <dir del paquete> [N del MC]
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import bench_mc as bm  # noqa: E402
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
NMC = int(sys.argv[2]) if len(sys.argv) > 2 else 30_000_000
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": "PROC-2 etapa 2", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
     "mppc_side": "sin registro (mapeo r4: face unknown en hod2018); René asigna los lados de las capas A y B"}
POS = -15.0 + 2.0 * np.arange(16)  # mm, centro de cada coordenada (CAD)
DEAD = 3
st3 = json.load(open(W + "/stage3.json"))["strip_orders"]
ORDER = {"may": st3["may"]["order"], "oct": st3["oct"]["order"]}


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1())}
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
DATA = {"may8_good": (e8[mc.good_time(e8["t"], "v8")], "may", "mayo v8, tiempo bueno"),
        "oct8": (xy.load_hits(W + "/hd64_oct8_hits.bin", 64), "oct", "octubre run 8")}
LAYERS = {"all": list(range(16)), "A": list(range(0, 16, 2)), "B": list(range(1, 16, 2))}


def hits(e, chans, T, g, field):
    chans = np.asarray(chans)
    h = e[field][:, chans] > T * g[chans]
    h[:, chans == DEAD] = False
    return h


def fits(x, y, ey):
    ok = np.isfinite(y) & np.isfinite(ey) & (ey > 0) & (y > 0)
    x, y, ey = x[ok], y[ok], ey[ok]
    out = {"points": int(ok.sum())}
    if len(x) < 3:
        return out
    w = 1 / ey ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * y))
    chi2 = float((((y - a - b * x) / ey) ** 2).sum())
    lo, hi = a - 16.5 * b, a + 16.5 * b
    elo = math.sqrt(cov[0, 0] + 16.5 ** 2 * cov[1, 1] - 2 * 16.5 * cov[0, 1])
    ehi = math.sqrt(cov[0, 0] + 16.5 ** 2 * cov[1, 1] + 2 * 16.5 * cov[0, 1])
    r = hi / lo
    out["line"] = {"intercept": float(a), "slope_per_mm": float(b), "slope_err": float(math.sqrt(cov[1, 1])), "chi2": chi2, "ndf": int(len(x) - 2),
                   "rel_slope_per_mm": float(b / a), "ratio_plus_over_minus_16p5": float(r),
                   "ratio_err": float(abs(r) * math.hypot(ehi / hi, elo / lo))}  # error sin la correlación de los extremos
    ly, ely = np.log(y), ey / y
    w2 = 1 / ely ** 2
    cov2 = np.linalg.inv(A.T @ (A * w2[:, None]))
    a2, b2 = cov2 @ (A.T @ (w2 * ly))
    out["exp"] = {"amplitude": float(math.exp(a2)), "slope_per_mm": float(b2), "slope_err": float(math.sqrt(cov2[1, 1])),
                  "chi2": float((((ly - a2 - b2 * x) / ely) ** 2).sum()), "ndf": int(len(x) - 2),
                  "ratio_plus_over_minus_16p5": float(math.exp(33 * b2)), "ratio_err": float(math.exp(33 * b2) * 33 * math.sqrt(cov2[1, 1]))}
    return out


# ---------------------------------------------------------------- MC: P(capa cruzada | tira de referencia cruzada)
rng = np.random.default_rng(42)
num = {(p, l): np.zeros(16) for p in "XY" for l in LAYERS}
den = {p: np.zeros(16) for p in "XY"}
for start in range(0, NMC, 2_000_000):
    m = min(2_000_000, NMC - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, bm.GEOMS["CAD"])
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    for top in ("X", "Y"):
        pth = bm.bar_paths(rt, 0.0, top)
        H = {c: pth[c] > 0 for c in "XY"}
        H["X"][:, 12] = False  # ch3 muerto, también en el MC
        for probe, ref in (("X", "Y"), ("Y", "X")):
            den[probe] += H[ref].sum(0)
            for l, cols in LAYERS.items():
                hit = H[probe][:, cols].any(1)
                num[(probe, l)] += (H[ref] & hit[:, None]).sum(0)
MCG = {f"{p}_{l}": (num[(p, l)] / den[p]).tolist() for p in "XY" for l in LAYERS}
S["mc_geometric"] = {"rays": NMC, "geometry": "CAD", "model": "cos3", "seed": 42, "P_layer_given_ref": MCG,
                     "note": "promedio de los dos órdenes de bloque; ch3 (coordenada X 12) apagado"}

# ---------------------------------------------------------------- datos
S["data"] = {}
for key, (e, cal, lab) in DATA.items():
    o, g = ORDER[cal], G[cal]
    d = {"label": lab, "events": int(len(e))}
    for T in (5, 10):
        for probe, ref in (("X", "Y"), ("Y", "X")):
            Rs, Rc = hits(e, o[f"h18_{ref}"], T, g, "aS"), hits(e, o[f"h18_{ref}"], T, g, "aC")
            Ps, Pc = hits(e, o[f"h18_{probe}"], T, g, "aS"), hits(e, o[f"h18_{probe}"], T, g, "aC")
            for l, cols in LAYERS.items():
                xs, xc = Ps[:, cols].any(1), Pc[:, cols].any(1)
                rows = []
                for j in range(16):
                    if ref == "X" and j == 12:  # la tira de referencia muerta
                        rows.append(None)
                        continue
                    rr = xy.xy_fraction(xs, Rs[:, j], xc, Rc[:, j])
                    # xy_fraction(Xs, Ys, …): e_X|Y = q/rY con Ys = la tira de referencia (segundo argumento)
                    eps, err = rr["e_X|Y_mix"], rr["e_X|Y_mix_err"]
                    mcg = MCG[f"{probe}_{l}"][j]
                    rows.append({"coord": j, "pos_mm": float(POS[j]), "q": rr["f_mix"], "q_err": rr["f_mix_err"], "r_ref": rr["rY"], "eps": eps, "eps_err": err,
                                 "eps_ctl": rr["e_X|Y_ctl"], "mc_geometric": mcg, "eps_over_mc": eps / mcg, "eps_over_mc_err": err / mcg})
                ok = [r_ for r_ in rows if r_ is not None]
                x = np.array([r_["pos_mm"] for r_ in ok])
                # control: sin la tira de referencia de ch29 (Y11 de mayo, la barra débil que A4 deja abierta)
                okx = [r_ for r_ in ok if o[f"h18_{ref}"][r_["coord"]] != 29]
                xx = np.array([r_["pos_mm"] for r_ in okx])
                d[f"{probe}_{l}_T{T}_no_ch29_ref"] = fits(xx, np.array([r_["eps_over_mc"] for r_ in okx]), np.array([r_["eps_over_mc_err"] for r_ in okx]))
                d[f"{probe}_{l}_T{T}"] = {"rows": rows,
                                          "fit_eps": fits(x, np.array([r_["eps"] for r_ in ok]), np.array([r_["eps_err"] for r_ in ok])),
                                          "fit_eps_over_mc": fits(x, np.array([r_["eps_over_mc"] for r_ in ok]), np.array([r_["eps_over_mc_err"] for r_ in ok]))}
                if l != "all":
                    f_ = d[f"{probe}_{l}_T{T}"]["fit_eps_over_mc"]
                    print(key, T, f"ε({probe}_{l}|{ref}_j)", "recta %.5f ± %.5f /mm (rel %.5f), razón %.3f ± %.3f; exp %.5f ± %.5f /mm, razón %.3f" % (
                        f_["line"]["slope_per_mm"], f_["line"]["slope_err"], f_["line"]["rel_slope_per_mm"], f_["line"]["ratio_plus_over_minus_16p5"],
                        f_["line"]["ratio_err"], f_["exp"]["slope_per_mm"], f_["exp"]["slope_err"], f_["exp"]["ratio_plus_over_minus_16p5"]))
    S["data"][key] = d

# ---------------------------------------------------------------- fuera del cuadrado: razón de razones
bg = json.load(open(W + "/bar_mc_g.json"))
cad = next(r for r in bg["runs"] if r["geometry"] == "CAD")
R0 = {p: 0.5 * sum((cad[t][f"f_{p}"]["value"] - cad[t]["f_XY"]["value"]) / cad[t]["f_XY"]["value"] for t in ("X arriba", "Y arriba")) for p in "XY"}
s2 = json.load(open(W + "/stage2.json"))
S["outside"] = {"R0_mc": R0, "assumptions": "R = (r − q)/q mide los impactos en un plano sin el otro por cada X·Y; con detección perfecta R0 = (f_P − f_XY)/f_XY "
                                            "(barras fuera del bloque ortogonal más los cruces sin la otra coordenada). R/R0 = 1 si la eficiencia es la misma dentro "
                                            "y fuera del cuadrado y la del otro plano es 1; una eficiencia del otro plano ε_o < 1 sube R en (1 − ε_o)/ε_o/R0.",
                "per_dataset": {}}
for key in ("may8_good", "oct8"):
    for T in ("5", "10"):
        r = s2["datasets"][key]["per_T"][T]
        q, qe = r["f_mix"], r["f_mix_err"]
        out = {}
        for p in "XY":
            rp = r[f"r{p}"]
            R = (rp - q) / q
            out[p] = {"r": rp, "q": q, "R": R, "R_err": rp / q ** 2 * qe, "R0": R0[p], "R_over_R0": R / R0[p], "R_over_R0_err": rp / q ** 2 * qe / R0[p]}
        S["outside"]["per_dataset"][f"{key}_T{T}"] = out
        print("fuera", key, T, {p: (round(v["R"], 3), round(v["R_over_R0"], 3), round(v["R_over_R0_err"], 3)) for p, v in out.items()})

# ---------------------------------------------------------------- tabla y figura
L = ["% PROC-2 etapa 2: eficiencia a lo largo de las barras de hod2018 (mayo v8 tiempo bueno, orden validado de la etapa 1; octubre run 8, r4).",
     "% ε(P_capa | R_j) = q_j / r_j con el estimador de xy_estimator.py restringido a la tira de referencia j, dividido por la probabilidad geométrica",
     "% del MC de barras (CAD, detección perfecta). Pendiente relativa por mm y razón entre los extremos del tramo (+16.5 / −16.5 mm) de la recta",
     f"% y de la exponencial. Lado de los MPPC: sin registro (capa A = coordenadas pares, B = impares). Commit {git}.",
     "\\begin{tabular}{lllrrrr}", "\\hline",
     "Conjunto & Plano y capa & $T$ & pendiente rel. (\\%/mm) & razón (recta) & razón (exp.) & $\\chi^2$/ndf (recta) \\\\", "\\hline"]
for key in ("may8_good", "oct8"):
    for probe in "XY":
        for l in ("A", "B"):
            for T in (5, 10):
                f_ = S["data"][key][f"{probe}_{l}_T{T}"]["fit_eps_over_mc"]
                ln, ex = f_["line"], f_["exp"]
                L.append(f"{S['data'][key]['label']} & {probe} capa {l} & {T} & {100 * ln['rel_slope_per_mm']:+.3f} $\\pm$ {100 * ln['slope_err'] / ln['intercept']:.3f} & "
                         f"{ln['ratio_plus_over_minus_16p5']:.3f} $\\pm$ {ln['ratio_err']:.3f} & {ex['ratio_plus_over_minus_16p5']:.3f} $\\pm$ {ex['ratio_err']:.3f} & "
                         f"{ln['chi2']:.1f}/{ln['ndf']} \\\\")
                if probe == "X":
                    n_ = S["data"][key][f"{probe}_{l}_T{T}_no_ch29_ref"]["line"]
                    L.append(f"\\quad sin la referencia ch29 & {probe} capa {l} & {T} & {100 * n_['rel_slope_per_mm']:+.3f} $\\pm$ "
                             f"{100 * n_['slope_err'] / n_['intercept']:.3f} & {n_['ratio_plus_over_minus_16p5']:.3f} $\\pm$ {n_['ratio_err']:.3f} & -- & "
                             f"{n_['chi2']:.1f}/{n_['ndf']} \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/eff_along_bar.tex", "w").write("\n".join(L) + "\n")
fig = {"panels": []}
for probe, ref in (("X", "Y"), ("Y", "X")):
    for T in (5, 10):
        pan = {"title": f"{probe} bars, reference {ref} strip, {T} p.e.", "xlabel": f"Position along the {probe} bars, {'y' if probe == 'X' else 'x'} (mm)",
               "series": []}
        for l in ("A", "B"):
            d_ = S["data"]["may8_good"][f"{probe}_{l}_T{T}"]
            ok = [r_ for r_ in d_["rows"] if r_ is not None]
            f_ = d_["fit_eps_over_mc"]
            pan["series"].append({"name": f"layer {l}", "x": [r_["pos_mm"] for r_ in ok], "y": [r_["eps_over_mc"] for r_ in ok], "e": [r_["eps_over_mc_err"] for r_ in ok],
                                  "line": [f_["line"]["intercept"], f_["line"]["slope_per_mm"]], "exp": [f_["exp"]["amplitude"], f_["exp"]["slope_per_mm"]]})
        fig["panels"].append(pan)
json.dump(fig, open(W + "/fig_along_bar.json", "w"), indent=1)
json.dump(S, open(W + "/along_bar.json", "w"), indent=1, default=float)
print("hecho")
