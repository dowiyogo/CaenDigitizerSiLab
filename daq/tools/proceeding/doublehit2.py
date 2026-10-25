#!/usr/bin/python3 -I
# doublehit2.py — PROC-2 etapa 4 (solo lectura): asimetría X/Y de la fracción de doble impacto de hod2018 (etapa 5 de PROC-1).
#  1. Por par vecino (k, k+1), mayo v8 tiempo bueno, 5 y 10 p.e., etiqueta = el otro plano > 10 p.e.: f2_k = pares verdaderos (k, k+1)
#     / triggers con impacto verdadero en k o en k+1, con pares verdaderos = mismo evento − mezclados y tiras verdaderas = en tiempo −
#     ventana de control (unión k ∪ k+1). MC de barras (CAD, cos³θ, detección perfecta, camino > 0.2 mm) con la misma definición;
#     razón medido / MC por par. Pares con unión derivada (etapa 1): X 3|4, 11|12, 12|13 (ch3 en la 12); Y 10|11.
#  2. Sensibilidad del MC a la geometría (detección perfecta, camino > 0.2 mm): traslape entre capas de 0.5, 0.75, 1.0, 1.25 y 1.5 mm
#     (desplazamiento entre capas fijo de 2.0 mm; ancho de barra = 2.0 mm + traslape) y separación entre capas de 2.5, 3.0 y 3.5 mm
#     entre centros. Solo para la discusión de tolerancias; no se concluye nada sobre el traslape real.
#  3. Dependencia de la etiqueta: fracción medida (definición de stage5.py) con la etiqueta a 10, 15 y 20 p.e.
# Salida: work/doublehit2.json, tables/doublehit_pairs.tex, tables/doublehit_geometry.tex. Uso: doublehit2.py <dir del paquete> [N MC]
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import bench_mc as bm  # noqa: E402
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
NMC = int(sys.argv[2]) if len(sys.argv) > 2 else 10_000_000
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": "PROC-2 etapa 4", "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
DEAD_X = 12
DERIVED = {"X": [(3, 4), (11, 12), (12, 13)], "Y": [(10, 11)]}
MAY = json.load(open(W + "/stage3.json"))["strip_orders"]["may"]["order"]


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


G = g_refit()
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
e = e8[mc.good_time(e8["t"], "v8")]


def hits(chans, T, field, ev=None):
    ev = e if ev is None else ev
    chans = np.asarray(chans)
    h = ev[field][:, chans] > T * G[chans]
    h[:, chans == 3] = False
    return h


# ---------------------------------------------------------------- 1. por par
def pair_stats(hs, hc):
    N = len(hs)
    out = []
    for k in range(15):
        a, b = hs[:, k], hs[:, k + 1]
        pin = int((a & b).sum())
        pmix = 0.5 * (int((a[:-1] & b[1:]).sum()) + int((a[1:] & b[:-1]).sum())) * N / (N - 1)
        u_s, u_c = int((a | b).sum()), int((hc[:, k] | hc[:, k + 1]).sum())
        den = u_s - u_c
        f = (pin - pmix) / den if den > 0 else float("nan")
        err = math.sqrt(pin + pmix / 2) / den if den > 0 else float("nan")
        out.append({"pair": [k, k + 1], "pairs_same": pin, "pairs_mixed": pmix, "union_in": u_s, "union_ctl": u_c, "f2": f, "f2_err": err})
    return out


S["pairs"] = {}
for probe, tagp in (("X", "Y"), ("Y", "X")):
    tag = hits(MAY[f"h18_{tagp}"], 10, "aS").any(1)
    for T in (5, 10):
        rows = pair_stats(hits(MAY[f"h18_{probe}"], T, "aS")[tag], hits(MAY[f"h18_{probe}"], T, "aC")[tag])
        for r in rows:
            r["derived_junction"] = tuple(r["pair"]) in DERIVED[probe]
            r["dead_strip"] = probe == "X" and DEAD_X in r["pair"]
        S["pairs"][f"{probe}_T{T}"] = rows

# MC por par (CAD, detección perfecta, camino > 0.2 mm, etiqueta = la otra coordenada cruzada)
rng = np.random.default_rng(42)
pp = {p: np.zeros(15) for p in "XY"}
uu = {p: np.zeros(15) for p in "XY"}


def mc_loop(n, cut=0.2, width=None, dz=None, accumulate=None):
    """Fracción de doble impacto (pares por trigger con impacto en la sonda) con barras de ancho `width` y capas a `dz`."""
    w0, d0 = bm.BAR_W, bm.LAYER_DZ
    if width is not None:
        bm.BAR_W = width
    if dz is not None:
        bm.LAYER_DZ = dz
    r_ = np.random.default_rng(42)
    tot = {p: [0, 0] for p in "XY"}
    try:
        for start in range(0, n, 2_000_000):
            m = min(2_000_000, n - start)
            r = bm.generate(r_, m, "cos3")
            t = bm.trigger(r, bm.GEOMS["CAD"])
            rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
            for top in ("X", "Y"):
                pth = bm.bar_paths(rt, 0.0, top)
                H = {c: pth[c] > cut for c in "XY"}
                H["X"][:, DEAD_X] = False
                for probe, tagp in (("X", "Y"), ("Y", "X")):
                    hp = H[probe][H[tagp].any(1)]
                    tot[probe][0] += int((hp[:, :-1] & hp[:, 1:]).sum())
                    tot[probe][1] += int(hp.any(1).sum())
                    if accumulate is not None:
                        accumulate[0][probe] += (hp[:, :-1] & hp[:, 1:]).sum(0)
                        accumulate[1][probe] += (hp[:, :-1] | hp[:, 1:]).sum(0)
    finally:
        bm.BAR_W, bm.LAYER_DZ = w0, d0
    return {p: tot[p][0] / tot[p][1] for p in "XY"}


nominal = mc_loop(NMC, accumulate=(pp, uu))
for probe in "XY":
    mcf = (pp[probe] / np.where(uu[probe] > 0, uu[probe], np.nan)).tolist()
    for T in (5, 10):
        for r, m in zip(S["pairs"][f"{probe}_T{T}"], mcf):
            r["mc_f2"] = m
            r["ratio"] = r["f2"] / m if m and m == m and m > 0 else float("nan")
            r["ratio_err"] = r["f2_err"] / m if m and m == m and m > 0 else float("nan")
for probe in "XY":
    for T in (5, 10):
        rr = [r for r in S["pairs"][f"{probe}_T{T}"] if not r["dead_strip"] and math.isfinite(r["ratio"])]
        w = np.array([1 / r["ratio_err"] ** 2 for r in rr])
        v = np.array([r["ratio"] for r in rr])
        m_ = float((w * v).sum() / w.sum())
        S["pairs"][f"{probe}_T{T}_summary"] = {"weighted_mean_ratio": m_, "chi2_uniform": float((w * (v - m_) ** 2).sum()), "ndf": len(rr) - 1,
                                               "min": float(v.min()), "max": float(v.max())}
        print(probe, T, "razón medido/MC por par: media %.3f, χ² uniforme %.1f/%d, mín %.3f máx %.3f" % (m_, S["pairs"][f"{probe}_T{T}_summary"]["chi2_uniform"],
              len(rr) - 1, v.min(), v.max()))

# ---------------------------------------------------------------- 1b. fase de la alternancia par/impar
# En el MC la fracción por par alterna entre pares (k par) e impares según qué capa está arriba en cada bloque (bench_mc: capa A
# arriba, sin registro) y según la orientación del orden de tiras (fijado salvo una inversión global, etapa 1). Se compara la forma
# medida (5 y 10 p.e.) con las cuatro combinaciones (χ² con un factor de escala libre). Las dos combinaciones {A arriba, orden
# invertido} y {B arriba, orden como está} dan la misma fase, así que los datos solo fijan el producto capa de arriba × orientación.
def layer_z_b_up(z_module, top="X", _orig=bm.layer_z):
    d = _orig(z_module, top)
    return {("X", "A"): d[("X", "B")], ("X", "B"): d[("X", "A")], ("Y", "A"): d[("Y", "B")], ("Y", "B"): d[("Y", "A")]}


def mc_pairs(b_up, reverse, n):
    orig = bm.layer_z
    bm.layer_z = layer_z_b_up if b_up else orig
    r_ = np.random.default_rng(42)
    pq = {p: np.zeros(15) for p in "XY"}
    uq = {p: np.zeros(15) for p in "XY"}
    try:
        for start in range(0, n, 2_000_000):
            m = min(2_000_000, n - start)
            r = bm.generate(r_, m, "cos3")
            t = bm.trigger(r, bm.GEOMS["CAD"])
            rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
            for top in ("X", "Y"):
                pth = bm.bar_paths(rt, 0.0, top)
                H = {c: (pth[c] > 0.2)[:, ::-1] if reverse else pth[c] > 0.2 for c in "XY"}
                H["X"][:, DEAD_X] = False
                for probe, tagp in (("X", "Y"), ("Y", "X")):
                    hp = H[probe][H[tagp].any(1)]
                    pq[probe] += (hp[:, :-1] & hp[:, 1:]).sum(0)
                    uq[probe] += (hp[:, :-1] | hp[:, 1:]).sum(0)
    finally:
        bm.layer_z = orig
    return {p: pq[p] / np.where(uq[p] > 0, uq[p], np.nan) for p in "XY"}


phase = {}
for b_up in (False, True):
    for rev in (False, True):
        lab = f"{'B' if b_up else 'A'} arriba, orden {'invertido' if rev else 'como está'}"
        mcp = mc_pairs(b_up, rev, NMC)
        d = {}
        for probe in "XY":
            for T in (5, 10):
                rr = S["pairs"][f"{probe}_T{T}"]
                idx = [k for k in range(15) if not rr[k]["dead_strip"] and math.isfinite(rr[k]["f2"]) and np.isfinite(mcp[probe][k])]
                y_ = np.array([rr[k]["f2"] for k in idx])
                e_ = np.array([rr[k]["f2_err"] for k in idx])
                m_ = mcp[probe][idx]
                s_ = float((y_ * m_ / e_ ** 2).sum() / (m_ ** 2 / e_ ** 2).sum())
                d[f"{probe}_T{T}"] = {"scale": s_, "chi2": float((((y_ - s_ * m_) / e_) ** 2).sum()), "ndf": len(idx) - 1}
        phase[lab] = {"mc_pairs": {p: mcp[p].tolist() for p in "XY"}, "fit": d}
        print("fase:", lab, {k: "%.1f/%d" % (v["chi2"], v["ndf"]) for k, v in d.items()})
S["alternation_phase"] = {"configs": phase, "note": "los datos fijan el producto (capa de arriba) × (orientación del orden); no cada uno por separado"}

# ---------------------------------------------------------------- 2. sensibilidad geométrica
geo = {"nominal_width_mm": bm.BAR_W, "nominal_layer_dz_mm": bm.LAYER_DZ, "rays": NMC, "path_cut_mm": 0.2, "overlap": [], "layer_spacing": []}
for ov in (0.5, 0.75, 1.0, 1.25, 1.5):
    f = mc_loop(NMC, width=2.0 + ov)
    geo["overlap"].append({"overlap_mm": ov, "bar_width_mm": 2.0 + ov, "f2_X": f["X"], "f2_Y": f["Y"]})
    print("traslape", ov, f)
for dz in (2.5, 3.0, 3.5):
    f = mc_loop(NMC, dz=dz)
    geo["layer_spacing"].append({"layer_dz_mm": dz, "f2_X": f["X"], "f2_Y": f["Y"]})
    print("separación", dz, f)
S["geometry_sensitivity"] = geo

# ---------------------------------------------------------------- 3. dependencia de la etiqueta (definición de stage5.py)
def f2_global(probe, tagp, T, tagT):
    tag = hits(MAY[f"h18_{tagp}"], tagT, "aS").any(1)
    hs, hc = hits(MAY[f"h18_{probe}"], T, "aS")[tag], hits(MAY[f"h18_{probe}"], T, "aC")[tag]
    N = len(hs)
    pin = int((hs[:, :-1] & hs[:, 1:]).sum())
    pmix = 0.5 * (int((hs[:-1, :-1] & hs[1:, 1:]).sum()) + int((hs[1:, :-1] & hs[:-1, 1:]).sum())) * N / (N - 1)
    Ps, Pc = hs.any(1).mean(), hc.any(1).mean()
    nt = N * (Ps - Pc) / (1 - Pc)
    return {"tag_events": N, "f2": (pin - pmix) / nt, "f2_stat": math.sqrt(pin + (pmix * (N - 1) / N) / 2) / nt}


S["tag_dependence"] = {f"{p}_T{T}": {str(tt): f2_global(p, q, T, tt) for tt in (10, 15, 20)} for p, q in (("X", "Y"), ("Y", "X")) for T in (5, 10)}
for k, v in S["tag_dependence"].items():
    print("etiqueta", k, {t: (round(x["f2"], 4), round(x["f2_stat"], 4), x["tag_events"]) for t, x in v.items()})

# ---------------------------------------------------------------- tablas
L = ["% PROC-2 etapa 4: fracción de doble impacto por par vecino (k, k+1) de hod2018, mayo v8 tiempo bueno, etiqueta en el otro plano > 10 p.e.",
     "% f2_k = pares verdaderos / triggers con impacto verdadero en k o k+1; MC: CAD, cos^3, detección perfecta, camino > 0.2 mm.",
     f"% $^d$: unión derivada (etapa 1); $^m$: incluye la tira muerta (ch3). Commit {git}.",
     "\\begin{tabular}{lrrrrrr}", "\\hline", "Par & $X$ 5 p.e. (medido/MC) & $X$ 10 p.e. & $Y$ 5 p.e. & $Y$ 10 p.e. & MC $X$ & MC $Y$ \\\\", "\\hline"]
for k in range(15):
    cells = []
    for probe in "XY":
        for T in (5, 10):
            r = S["pairs"][f"{probe}_T{T}"][k]
            mark = ("$^d$" if r["derived_junction"] else "") + ("$^m$" if r["dead_strip"] else "")
            cells.append("--" if r["dead_strip"] else f"{r['ratio']:.2f} $\\pm$ {r['ratio_err']:.2f}{mark}")
    mx, my = S["pairs"]["X_T5"][k]["mc_f2"], S["pairs"]["Y_T5"][k]["mc_f2"]
    L.append(f"{k}--{k + 1} & " + " & ".join(cells) + f" & {mx:.3f} & {my:.3f} \\\\")
for probe in "XY":
    for T in (5, 10):
        sm = S["pairs"][f"{probe}_T{T}_summary"]
        L.append(f"% {probe} {T} p.e.: media pesada {sm['weighted_mean_ratio']:.3f}, chi2 contra uniforme {sm['chi2_uniform']:.1f}/{sm['ndf']}")
for lab, v in S["alternation_phase"]["configs"].items():
    L.append(f"% fase {lab}: " + ", ".join(f"{k} chi2 {x['chi2']:.1f}/{x['ndf']}" for k, x in v["fit"].items()))
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/doublehit_pairs.tex", "w").write("\n".join(L) + "\n")
G2 = ["% PROC-2 etapa 4: sensibilidad del MC de barras (CAD, cos^3, detección perfecta, camino > 0.2 mm) a la geometría; solo para la discusión",
      f"% de tolerancias (no se concluye el traslape real). Pares vecinos por trigger con impacto en el plano. {NMC:.0e} rayos por punto. Commit {git}.",
      "\\begin{tabular}{lrr}", "\\hline", "Geometría & $f_2^X$ & $f_2^Y$ \\\\", "\\hline"]
for x in geo["overlap"]:
    G2.append(f"traslape {x['overlap_mm']:.2f} mm (ancho {x['bar_width_mm']:.2f} mm) & {x['f2_X']:.4f} & {x['f2_Y']:.4f} \\\\")
for x in geo["layer_spacing"]:
    G2.append(f"separación entre capas {x['layer_dz_mm']:.1f} mm & {x['f2_X']:.4f} & {x['f2_Y']:.4f} \\\\")
G2 += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/doublehit_geometry.tex", "w").write("\n".join(G2) + "\n")
json.dump(S, open(W + "/doublehit2.json", "w"), indent=1, default=float)
print("hecho")
