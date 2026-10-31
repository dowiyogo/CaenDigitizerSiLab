#!/usr/bin/python3 -I
# eff_by_position.py — seguimiento 2026-10-31 (6), punto 2: eficiencia absoluta por capa de hod2018 y por posición a lo largo de la barra,
# octubre run 8 y mayo v8 (tiempo bueno nuevo), sin dividir por ε_MC, y el ε_MC geométrico de cada intervalo.
# Estimador y resta de accidentales: los de daq/tools/proceeding/along_bar.py (xy_estimator.xy_fraction con eventos mezclados):
#   para la capa sonda L (A = coordenadas pares, B = impares) del plano P y un intervalo k de posición a lo largo de sus barras, la referencia es
#   la unión de las tiras de la coordenada ortogonal R cuya posición cae en el intervalo (impacto en la ventana signal / de control de alguna de
#   ellas); ε(L | R_k) = q/r_R = "e_X|Y_mix" de xy_fraction(L_signal, R_k signal, L_control, R_k control), con su error.
#   Impacto = amplitud de pico > T·G (T = 5 y 10 p.e.; G de mayo = reajuste de PROC, de octubre = c1; mediana para los canales sin G);
#   ch3 (coordenada X 12) muerto: fuera de la capa sonda XA y fuera de las referencias, como en along_bar.py. Órdenes de tiras: stage3.json del
#   paquete PROC nuevo (octubre r4, mayo el validado).
# Posición: distancia al extremo de lectura en [21, 29), [29, 37), [37, 45), [45, 54] mm (4 tiras de referencia por intervalo, centros d = 37.5 +
#   (−15 + 2j) mm desde el extremo del índice 0); extremo de lectura de cada capa = el más brillante del seguimiento 3 (layer_from_channels.json,
#   igual en mayo y octubre: XA e YB en el extremo del índice 0, distancia = d; XB e YA en el opuesto, distancia = 75 − d) [supuesto: el extremo
#   más brillante es el del MPPC].
# ε_MC: P(alguna barra de la capa L cruzada | alguna tira de referencia del intervalo cruzada) del MC de barras de along_bar.py (bench_mc, CAD,
#   cos³θ, detección perfecta, los dos órdenes de bloque sumados, semilla 42, 3·10^7 rayos), con ch3 apagado (como along_bar.py: ε/ε_MC es la
#   eficiencia intrínseca de las barras vivas) y, como referencia, con ch3 encendido (geometría pura).
# Variante: sin ch29 (octubre) / ch30 (mayo) como tira de referencia (canales anómalos del seguimiento 1); la referencia se puede quitar sin
#   sesgo porque ε es condicional a ella (el ε_MC se recalcula sin esa tira).
# Salida: <out>/eff_<campaña>_<capa>_T<T>.csv, <out>/eff_by_position.json. Uso: eff_by_position.py <out dir> [N del MC]
import json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import bench_mc as bm  # noqa: E402
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
NMC = int(sys.argv[2]) if len(sys.argv) > 2 else 30_000_000
os.makedirs(OUT, exist_ok=True)
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD = 3
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
LAYERS = {"A": list(range(0, 16, 2)), "B": list(range(1, 16, 2))}
CAMP = {"oct": ("oct2026_run8", "octubre run 8", 29), "may": ("may2026_v8_goodtime", "mayo v8, tiempo bueno nuevo", 30)}


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1())}
e8 = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
DATA = {"may": e8[mc.good_time(e8["t"], "v8")], "oct": xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)}


def dist(probe, layer, j):
    s = SIGN[f"oct2026_run8/{probe}{layer}"]["sign"]
    assert s == SIGN[f"may2026_v8_goodtime/{probe}{layer}"]["sign"]
    return D[j] if s < 0 else 75.0 - D[j]


def bin_js(probe, layer, k):
    lo, hi = BINS[k]
    return [j for j in range(16) if lo <= dist(probe, layer, j) < hi + (1e-9 if hi == 54 else 0)]


# ---------------------------------------------------------------- MC
rng = np.random.default_rng(42)
mcnum, mcden = {}, {}
for start in range(0, NMC, 2_000_000):
    m = min(2_000_000, NMC - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, bm.GEOMS["CAD"])
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    for top in ("X", "Y"):
        pth = bm.bar_paths(rt, 0.0, top)
        for dead in (True, False):
            H = {c: pth[c] > 0 for c in "XY"}
            if dead:
                H["X"][:, 12] = False
            for probe, ref in (("X", "Y"), ("Y", "X")):
                for L, cols in LAYERS.items():
                    hitL = H[probe][:, cols].any(1)
                    for k in range(4):
                        for excl_ch in (None, "oct", "may"):
                            js = bin_js(probe, L, k)
                            if excl_ch:
                                bad = CAMP[excl_ch][2]
                                js = [j for j in js if ORDER[excl_ch][f"h18_{ref}"][j] != bad]
                            if ref == "X":
                                js = [j for j in js if j != 12]
                            key = (dead, probe, L, k, excl_ch)
                            R = H[ref][:, js].any(1)
                            mcden[key] = mcden.get(key, 0) + int(R.sum())
                            mcnum[key] = mcnum.get(key, 0) + int((R & hitL).sum())
    print(f"MC {start + m}/{NMC}", flush=True)


def eps_mc(dead, probe, L, k, excl):
    n, d = mcnum[(dead, probe, L, k, excl)], mcden[(dead, probe, L, k, excl)]
    p = n / d
    return p, math.sqrt(p * (1 - p) / d), n, d


# ---------------------------------------------------------------- datos
def hits(e, chans, T, g, field):
    chans = np.asarray(chans)
    h = e[field][:, chans] > T * g[chans]
    h[:, chans == DEAD] = False
    return h


res = {"bins_mm_from_readout": BINS, "mc": {"rays": NMC, "geometry": "CAD", "model": "cos3", "seed": 42}, "data": {}}
for camp, (tag, label, bad) in CAMP.items():
    e, o, g = DATA[camp], ORDER[camp], G[camp]
    res["data"][camp] = {"label": label, "events": int(len(e))}
    for T in (5, 10):
        for probe, ref in (("X", "Y"), ("Y", "X")):
            Rs, Rc = hits(e, o[f"h18_{ref}"], T, g, "aS"), hits(e, o[f"h18_{ref}"], T, g, "aC")
            Ps, Pc = hits(e, o[f"h18_{probe}"], T, g, "aS"), hits(e, o[f"h18_{probe}"], T, g, "aC")
            for L, cols in LAYERS.items():
                xs, xc = Ps[:, cols].any(1), Pc[:, cols].any(1)
                rows = []
                for variant, excl in (("all_ref", None), ("no_bad_ref", camp)):
                    for k in range(4):
                        js = bin_js(probe, L, k)
                        if excl:
                            js = [j for j in js if o[f"h18_{ref}"][j] != bad]
                        if ref == "X":
                            js = [j for j in js if j != 12]
                        rr = xy.xy_fraction(xs, Rs[:, js].any(1), xc, Rc[:, js].any(1))
                        p_d, ep_d, _, _ = eps_mc(True, probe, L, k, excl)
                        p_a, ep_a, _, _ = eps_mc(False, probe, L, k, excl)
                        eps, err = rr["e_X|Y_mix"], rr["e_X|Y_mix_err"]
                        rows.append({"variant": variant, "bin": k, "dist_lo_mm": BINS[k][0], "dist_hi_mm": BINS[k][1], "ref_coords": " ".join(map(str, js)),
                                     "ref_channels": " ".join(str(o[f"h18_{ref}"][j]) for j in js), "N": rr["N"], "N_XY": rr["N_XY"], "acc_mix": rr["acc_mix"],
                                     "r_ref": rr["rY"], "q": rr["f_mix"], "q_err": rr["f_mix_err"], "eps": eps, "eps_err": err,
                                     "eps_ctl": rr["e_X|Y_ctl"], "eps_ctl_err": rr["e_X|Y_ctl_err"],
                                     "eps_mc_ch3off": p_d, "eps_mc_ch3off_err": ep_d, "eps_mc_ch3on": p_a, "eps_mc_ch3on_err": ep_a,
                                     "eps_over_mc": eps / p_d, "eps_over_mc_err": math.hypot(err / p_d, eps * ep_d / p_d ** 2)})
                fn = f"{OUT}/eff_{camp}_{probe}{L}_T{T}.csv"
                with open(fn, "w") as f:
                    f.write(f"# {label}; capa {probe}{L} (canales {[o[f'h18_{probe}'][i] for i in cols if o[f'h18_{probe}'][i] != DEAD]}), referencia = tiras "
                            f"{ref} del intervalo; umbral {T} p.e.; extremo de lectura = índice {'0' if SIGN[f'oct2026_run8/{probe}{L}']['sign'] < 0 else '15'} "
                            f"de las tiras {ref}; eps = q/r_ref con mezclados (along_bar.py), sin dividir por eps_mc\n")
                    keys = list(rows[0])
                    f.write(",".join(keys) + "\n")
                    for r_ in rows:
                        f.write(",".join(f"{r_[k]:.6f}" if isinstance(r_[k], float) else str(r_[k]) for k in keys) + "\n")
                res["data"][camp][f"{probe}{L}_T{T}"] = rows
                print(f"{camp} T{T} {probe}{L}: " + "; ".join(f"[{r_['dist_lo_mm']},{r_['dist_hi_mm']}) ε {r_['eps']:.4f}±{r_['eps_err']:.4f} "
                                                              f"ε_MC {r_['eps_mc_ch3off']:.4f} ε/ε_MC {r_['eps_over_mc']:.3f}" for r_ in rows if r_["variant"] == "all_ref"))
json.dump(res, open(f"{OUT}/eff_by_position.json", "w"), indent=1, default=float)
print("hecho")
