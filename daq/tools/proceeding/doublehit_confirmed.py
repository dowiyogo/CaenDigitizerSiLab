#!/usr/bin/python3 -I
# doublehit_confirmed.py — PROC-3 etapa 3 (medido; MC estimado): fracción de doble impacto de hod2018 (definición de stage5.py: pares
# vecinos verdaderos por trigger con impacto verdadero en el plano; etiqueta el otro plano > 10 p.e.) con todos los pares y sin los pares
# cuya unión es derivada. Mayo v8 tiempo bueno: X 3|4 y 11|12 (además 12|13 contiene ch3 muerto), Y 10|11 (etapa 1). Octubre run 8: r4
# marca como derivada solo la posición de ch3 (coordenada X 12), así que salen X 11|12 y 12|13 (ya vacíos por la tira muerta); Y ninguno.
# El MC (bench_mc, CAD, cos³θ, detección perfecta, camino > 0 y > 0.2 mm, los dos órdenes de bloque) excluye los mismos pares.
# T = 3 … 12 p.e. Salida: work/doublehit_confirmed.json, tables/doublehit_confirmed.tex. Uso: doublehit_confirmed.py <dir del paquete> [N MC]
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
TS = list(range(3, 13))
DEAD_X = 12
EXCL = {"may": {"X": [(3, 4), (11, 12), (12, 13)], "Y": [(10, 11)]}, "oct": {"X": [(11, 12), (12, 13)], "Y": []}}
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
DATA = {"may8_good": (e8[mc.good_time(e8["t"], "v8")], "may", "mayo v8, tiempo bueno"), "oct8": (xy.load_hits(W + "/hd64_oct8_hits.bin", 64), "oct", "octubre run 8")}


def hits(e, chans, T, g, field):
    chans = np.asarray(chans)
    h = e[field][:, chans] > T * g[chans]
    h[:, chans == 3] = False
    return h


def f2(hs, hc, keep):
    """pares verdaderos en los pares `keep` (índices k del par (k, k+1)) por trigger con impacto verdadero en la sonda (stage5.py)"""
    N = len(hs)
    a, b = hs[:, :-1][:, keep], hs[:, 1:][:, keep]
    pin = int((a & b).sum())
    pmix = 0.5 * (int((a[:-1] & b[1:]).sum()) + int((a[1:] & b[:-1]).sum())) * N / (N - 1)
    Ps, Pc = hs.any(1).mean(), hc.any(1).mean()
    nt = N * (Ps - Pc) / (1 - Pc)
    return (pin - pmix) / nt, math.sqrt(pin + pmix / 2) / nt


# ---------------------------------------------------------------- MC por par (para excluir los mismos)
rng = np.random.default_rng(42)
mcp = {(cut, p): np.zeros(15) for cut in (0.0, 0.2) for p in "XY"}
mcn = {(cut, p): 0 for cut in (0.0, 0.2) for p in "XY"}
for start in range(0, NMC, 2_000_000):
    m = min(2_000_000, NMC - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, bm.GEOMS["CAD"])
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    for top in ("X", "Y"):
        pth = bm.bar_paths(rt, 0.0, top)
        for cut in (0.0, 0.2):
            H = {c: pth[c] > cut for c in "XY"}
            H["X"][:, DEAD_X] = False
            for probe, tagp in (("X", "Y"), ("Y", "X")):
                hp = H[probe][H[tagp].any(1)]
                mcp[(cut, probe)] += (hp[:, :-1] & hp[:, 1:]).sum(0)
                mcn[(cut, probe)] += int(hp.any(1).sum())


def mc_f2(cut, probe, keep):
    return float(mcp[(cut, probe)][keep].sum() / mcn[(cut, probe)])


S = {"stage": "PROC-3 etapa 3", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "excluded_pairs": EXCL, "mc_rays": NMC,
     "mc": "CAD, cos3, detección perfecta, promedio de los dos órdenes de bloque, misma exclusión", "data": {}}
for key, (e, cal, lab) in DATA.items():
    o, g = ORDER[cal], G[cal]
    d = {"label": lab}
    for probe, tagp in (("X", "Y"), ("Y", "X")):
        tag = hits(e, o[f"h18_{tagp}"], 10, g, "aS").any(1)
        allk = list(range(15))
        conf = [k for k in allk if (k, k + 1) not in EXCL[cal][probe]]
        rows = []
        for T in TS:
            hs, hc = hits(e[tag], o[f"h18_{probe}"], T, g, "aS"), hits(e[tag], o[f"h18_{probe}"], T, g, "aC")
            fa, ea = f2(hs, hc, allk)
            fc, ec = f2(hs, hc, conf)
            row = {"T": T, "all": fa, "all_err": ea, "confirmed": fc, "confirmed_err": ec}
            for cut in (0.0, 0.2):
                ma, mcf = mc_f2(cut, probe, allk), mc_f2(cut, probe, conf)
                row[f"mc_all_cut{cut}"], row[f"mc_conf_cut{cut}"] = ma, mcf
                row[f"ratio_all_cut{cut}"], row[f"ratio_conf_cut{cut}"] = fa / ma, fc / mcf
                row[f"ratio_all_cut{cut}_err"], row[f"ratio_conf_cut{cut}_err"] = ea / ma, ec / mcf
            rows.append(row)
        d[probe] = {"pairs_used": conf, "rows": rows}
        r5 = next(x for x in rows if x["T"] == 5)
        print(key, probe, "5 p.e.: todos %.3f ± %.3f (/MC %.3f), confirmados %.3f ± %.3f (/MC %.3f)" % (r5["all"], r5["all_err"], r5["ratio_all_cut0.0"],
              r5["confirmed"], r5["confirmed_err"], r5["ratio_conf_cut0.0"]))
    S["data"][key] = d
L = ["% PROC-3 etapa 3: fracción de doble impacto de hod2018 (pares vecinos verdaderos por trigger con impacto en el plano; etiqueta > 10 p.e.)",
     "% con todos los pares y sin los pares de unión derivada (mayo: X 3|4, 11|12, 12|13; Y 10|11; octubre r4: X 11|12, 12|13). Razón medido/MC con el",
     f"% MC del CAD (detección perfecta, camino > 0) que excluye los mismos pares. Commit {git}.",
     "\\begin{tabular}{lllrrrr}", "\\hline", "Conjunto & Plano & $T$ & $f_2$ todos & /MC & $f_2$ confirmados & /MC \\\\", "\\hline"]
for key in ("may8_good", "oct8"):
    for probe in "XY":
        for r in S["data"][key][probe]["rows"]:
            if r["T"] in (5, 8, 10, 12):
                L.append(f"{S['data'][key]['label']} & {probe} & {r['T']} & {r['all']:.3f} $\\pm$ {r['all_err']:.3f} & {r['ratio_all_cut0.0']:.3f} $\\pm$ "
                         f"{r['ratio_all_cut0.0_err']:.3f} & {r['confirmed']:.3f} $\\pm$ {r['confirmed_err']:.3f} & {r['ratio_conf_cut0.0']:.3f} $\\pm$ "
                         f"{r['ratio_conf_cut0.0_err']:.3f} \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/doublehit_confirmed.tex", "w").write("\n".join(L) + "\n")
json.dump(S, open(W + "/doublehit_confirmed.json", "w"), indent=1, default=float)
print("hecho")
