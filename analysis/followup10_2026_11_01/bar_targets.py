#!/usr/bin/python3 -I
# bar_targets.py — 2026-11-01: objetivos por barra de hod2018 para la comparación con Geant4. Mayo v8 (tiempo bueno, 77 558 triggers) y
# octubre run 8, por separado. Una fila por barra, campaña e intervalo de posición.
# Barra = (coordenada, índice de tira i del orden de stage3.json, capa A = i par, B = i impar); el canal es el de cada campaña (Y se recableó:
# el mismo índice es el mismo MPPC, seguimiento 7). Centro en el banco = −15 + 2i mm (convención de bench_mc: índice creciente con la
# coordenada; el sentido absoluto del orden no está fijado). Extremo de lectura por capa = el más brillante del seguimiento 3
# (layer_from_channels.json, igual en mayo y octubre) [supuesto: el extremo más brillante es el del MPPC].
# Intervalos: distancia al extremo de lectura [21, 29), [29, 37), [37, 45), [45, 54] mm; cada uno = 4 tiras de la coordenada ortogonal
# (centros d = 37.5 − 15 + 2j mm desde el extremo del índice 0), como eff_by_position.py (seguimiento 6).
# Etiqueta: alguna tira ortogonal del intervalo > 10 p.e. en la ventana signal (aS). ch3 (X12) muerto: sin filas útiles como sonda y fuera de
# las etiquetas; ch29/ch30 se usan (son el MPPC Y13, normal; seguimiento 7).
# Amplitud: sonda > 5 p.e. medidos en aS con la etiqueta en el mismo evento, menos los mezclados (etiqueta del evento siguiente) × N/(N − 1);
# histograma en bins de 1 p.e. medido desde 0 (definición del seguimiento 1), bins netos negativos a 0; mediana interpolada linealmente
# dentro del bin (para no cuantizar en 1 p.e. con pocas entradas). p.e. primarios = medidos / F con F de campaña: 1.75 mayo, 1.46 octubre.
# Error: remuestreo de Poisson de los dos histogramas, 200 réplicas. Mediana no fiable si n_net < N_MIN = 15: con la prueba de 2·10^6
# rayos, desde n_net = 15 el error relativo por remuestreo queda bajo ~0.3 (una fila 0.34); debajo llega a 0.9-5.9.
# Eficiencia: ε = q/r_ref de xy_estimator.xy_fraction (mezclados) con la sonda > T·G (T = 5, 10) y la etiqueta > 10 p.e. (signal y control);
# ε_MC = P(barra cruzada | alguna tira de la etiqueta cruzada), MC de barras (bench_mc, CAD, cos³θ, detección perfecta, ch3 apagado, los
# dos órdenes de bloque sumados, semilla 42).
# L_eff: ajuste ponderado de ln(mediana en p.e. medidos) contra la distancia media de las tiras del intervalo, solo con intervalos fiables;
# "not determined" con menos de 3 intervalos fiables o pendiente ≥ 0.
# Salida: <out>/bar_targets.csv, bar_targets.json. Uso: bar_targets.py <out dir> [N del MC]
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
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
TAG, PROBE_T, N_MIN, NBOOT = 10.0, 5.0, 15, 200
FCT = {"may": 1.75, "oct": 1.46}
CAMP = {"may": "mayo v8, tiempo bueno", "oct": "octubre run 8"}
EDG = np.arange(0, 301, 1.0)
RNG = np.random.default_rng(1101)


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1())}
e8 = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
DATA = {"may": e8[mc.good_time(e8["t"], "v8")], "oct": xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)}


def readout_sign(probe, layer):
    s = SIGN[f"oct2026_run8/{probe}{layer}"]["sign"]
    assert s == SIGN[f"may2026_v8_goodtime/{probe}{layer}"]["sign"]
    return s


def dist(probe, layer, j):
    return D[j] if readout_sign(probe, layer) < 0 else 75.0 - D[j]


def bin_js(probe, layer, k):
    lo, hi = BINS[k]
    js = [j for j in range(16) if lo <= dist(probe, layer, j) < hi + (1e-9 if hi == 54 else 0)]
    return [j for j in js if not (probe == "Y" and j == DEAD_X)]  # la referencia X12 (ch3) está muerta


# ---------------------------------------------------------------- MC
rng = np.random.default_rng(42)
mnum = np.zeros((2, 16, 4), np.int64)
mden = np.zeros((2, 16, 4), np.int64)
for start in range(0, NMC, 2_000_000):
    m = min(2_000_000, NMC - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, bm.GEOMS["CAD"])
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    for top in ("X", "Y"):
        pth = bm.bar_paths(rt, 0.0, top)
        H = {c: pth[c] > 0 for c in "XY"}
        H["X"][:, DEAD_X] = False
        for pi, (probe, ref) in enumerate((("X", "Y"), ("Y", "X"))):
            for i in range(16):
                L = "A" if i % 2 == 0 else "B"
                for k in range(4):
                    R = H[ref][:, bin_js(probe, L, k)].any(1)
                    mden[pi, i, k] += int(R.sum())
                    mnum[pi, i, k] += int((R & H[probe][:, i]).sum())
    print(f"MC {start + m}/{NMC}", flush=True)


# ---------------------------------------------------------------- datos
def median_interp(net):
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    before = cs[b - 1] if b > 0 else 0.0
    return float(EDG[b] + (0.5 * n - before) / net[b])


def amp_stats(hs, hm, s):
    net = np.clip(hs - s * hm, 0, None)
    med = median_interp(net)
    boot = []
    for _ in range(NBOOT):
        x = median_interp(np.clip(RNG.poisson(hs) - s * RNG.poisson(hm), 0, None))
        if np.isfinite(x):
            boot.append(x)
    return {"n_same": int(hs.sum()), "n_mix_scaled": float(s * hm.sum()), "n_net": float(hs.sum() - s * hm.sum()), "median_pe_meas": med,
            "median_pe_meas_err": float(np.std(boot)) if len(boot) > 1 else np.nan}


def lfit(d, m, e):
    d, m, e = (np.asarray(v, float) for v in (d, m, e))
    ok = np.isfinite(m) & np.isfinite(e) & (e > 0) & (m > 0)
    d, m, e = d[ok], m[ok], e[ok]
    if len(d) < 3:
        return {"status": f"not determined: {len(d)} reliable bins (< 3)"}
    y, w = np.log(m), (m / e) ** 2
    X = np.vstack([np.ones_like(d), d - 37.5]).T
    C = np.linalg.inv((X * w[:, None]).T @ X)
    b = C @ (X * w[:, None]).T @ y
    k, ek = float(b[1]), math.sqrt(C[1, 1])
    chi = float((((y - X @ b) ** 2) * w).sum())
    out = {"slope_per_mm": k, "slope_err": ek, "chi2": chi, "ndf": int(len(d) - 2), "n_bins": int(len(d))}
    if k >= 0:
        out["status"] = "not determined: slope >= 0"
        return out
    out.update(status="ok", L_eff_mm=-1 / k, L_eff_err_mm=ek / k ** 2)
    return out


rows, fits = [], {}
for camp in ("may", "oct"):
    e, o, g = DATA[camp], ORDER[camp], G[camp]
    N = len(e)
    s = N / (N - 1)
    pe = {f: e[f][:, :32] / g[None, :32] for f in ("aS", "aC")}
    for pi, (probe, ref) in enumerate((("X", "Y"), ("Y", "X"))):
        refch = np.asarray(o[f"h18_{ref}"])
        tagS = pe["aS"][:, refch] > TAG
        tagC = pe["aC"][:, refch] > TAG
        tagS[:, refch == DEAD] = False
        tagC[:, refch == DEAD] = False
        for i, ch in enumerate(o[f"h18_{probe}"]):
            L = "A" if i % 2 == 0 else "B"
            dead = ch == DEAD
            sgn = readout_sign(probe, L)
            base = {"campaign": camp, "coord": probe, "layer": L, "strip_index": i, "channel": int(ch), "bench_centre_mm": -15.0 + 2 * i,
                    "readout_end": f"index {'0' if sgn < 0 else '15'} end of the {ref} strips (bench {ref.lower()} = {'-' if sgn < 0 else '+'}37.5 mm) [S]",
                    "dead": bool(dead)}
            ps = pe["aS"][:, ch]
            dmeds, derrs, dpos = [], [], []
            for k in range(4):
                js = bin_js(probe, L, k)
                r = dict(base, bin=k, dist_lo_mm=BINS[k][0], dist_hi_mm=BINS[k][1], dist_mean_mm=float(np.mean([dist(probe, L, j) for j in js])),
                         tag_strips=" ".join(map(str, js)), tag_channels=" ".join(str(refch[j]) for j in js), F_ct=FCT[camp])
                TS, TC = tagS[:, js].any(1), tagC[:, js].any(1)
                if dead:
                    rows.append(dict(r, reliable=False, note="dead channel (ch3)"))
                    continue
                TSn = np.concatenate([TS[1:], [False]])
                sel = ps > PROBE_T
                hs = np.histogram(ps[sel & TS], bins=EDG)[0].astype(float)
                hm = np.histogram(ps[sel & TSn], bins=EDG)[0].astype(float)
                a = amp_stats(hs, hm, s)
                a["median_pe_prim"] = a["median_pe_meas"] / FCT[camp]
                a["median_pe_prim_err"] = a["median_pe_meas_err"] / FCT[camp]
                rel = a["n_net"] >= N_MIN and np.isfinite(a["median_pe_meas_err"])
                r.update(a, reliable=bool(rel), note="" if rel else f"n_net < {N_MIN}: median unreliable")
                pmc, nmc = mnum[pi, i, k] / mden[pi, i, k], mden[pi, i, k]
                emc = math.sqrt(pmc * (1 - pmc) / nmc)
                r.update(eps_mc=pmc, eps_mc_err=emc)
                for T in (5, 10):
                    hS, hC = e["aS"][:, ch] > T * g[ch], e["aC"][:, ch] > T * g[ch]
                    q = xy.xy_fraction(hS, TS, hC, TC)
                    ep, ee = q["e_X|Y_mix"], q["e_X|Y_mix_err"]
                    r.update({f"eps_T{T}": ep, f"eps_T{T}_err": ee, f"eps_over_mc_T{T}": ep / pmc,
                              f"eps_over_mc_T{T}_err": math.hypot(ee / pmc, ep * emc / pmc ** 2), f"tag_events_T{T}": int(q["N_Y"]),
                              f"r_tag_T{T}": q["rY"]})
                if rel:
                    dmeds.append(a["median_pe_meas"]); derrs.append(a["median_pe_meas_err"]); dpos.append(r["dist_mean_mm"])
                rows.append(r)
            if not dead:
                f = lfit(dpos, dmeds, derrs)
                fits[f"{camp}_{probe}{i:02d}"] = dict(f, channel=int(ch))
                for r in rows[-4:]:
                    r.update(L_eff_status=f["status"], L_eff_mm=f.get("L_eff_mm", np.nan), L_eff_err_mm=f.get("L_eff_err_mm", np.nan),
                             L_eff_chi2=f.get("chi2", np.nan), L_eff_ndf=f.get("ndf", ""))
            print(f"{camp} {probe}{i:02d} ch{ch:02d}: " + ("muerto" if dead else
                  " | ".join(f"n {r_['n_net']:.0f} med {r_['median_pe_meas']:.2f}±{r_['median_pe_meas_err']:.2f} ε/εMC5 {r_['eps_over_mc_T5']:.2f}"
                             for r_ in rows[-4:]) + f" | L {fits[f'{camp}_{probe}{i:02d}']['status']} "
                  + (f"{fits[f'{camp}_{probe}{i:02d}']['L_eff_mm']:.0f}±{fits[f'{camp}_{probe}{i:02d}']['L_eff_err_mm']:.0f}"
                     if fits[f'{camp}_{probe}{i:02d}']['status'] == 'ok' else "")), flush=True)

KEYS = ["campaign", "coord", "layer", "strip_index", "channel", "bench_centre_mm", "readout_end", "dead", "bin", "dist_lo_mm", "dist_hi_mm", "dist_mean_mm",
        "tag_strips", "tag_channels", "n_same", "n_mix_scaled", "n_net", "reliable", "median_pe_meas", "median_pe_meas_err", "F_ct", "median_pe_prim",
        "median_pe_prim_err", "eps_T5", "eps_T5_err", "eps_T10", "eps_T10_err", "eps_mc", "eps_mc_err", "eps_over_mc_T5", "eps_over_mc_T5_err",
        "eps_over_mc_T10", "eps_over_mc_T10_err", "tag_events_T5", "L_eff_status", "L_eff_mm", "L_eff_err_mm", "L_eff_chi2", "L_eff_ndf", "note"]


def fmt(v):
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return "" if not math.isfinite(v) else f"{v:.6g}"
    return "" if v is None else str(v).replace(",", ";")


with open(OUT + "/bar_targets.csv", "w") as f:
    f.write(f"# bar_targets.py (2026-11-01): una fila por barra de hod2018, campaña e intervalo de distancia al extremo de lectura. Etiqueta = alguna tira "
            f"ortogonal del intervalo > {TAG:g} p.e.; amplitud de la sonda > {PROBE_T:g} p.e. medidos, mismo evento − mezclados; mediana interpolada en bins "
            f"de 1 p.e. medido; p.e. primarios = medidos / F_ct; error por remuestreo de Poisson ({NBOOT}); reliable = n_net >= {N_MIN}. "
            f"eps = q/r_ref (mezclados), eps_mc del MC de barras (CAD, cos^3, ch3 apagado, {NMC} rayos). L_eff: ajuste exponencial de la mediana "
            f"(intervalos fiables). Centro en el banco con la convención de bench_mc (sentido absoluto del orden sin fijar).\n")
    f.write(",".join(KEYS) + "\n")
    for r in rows:
        f.write(",".join(fmt(r.get(k)) for k in KEYS) + "\n")
json.dump({"tag_pe": TAG, "probe_threshold_pe": PROBE_T, "n_min_reliable": N_MIN, "bootstrap": NBOOT, "F_ct": FCT, "mc_rays": NMC,
           "events": {k: int(len(v)) for k, v in DATA.items()}, "fits": fits, "rows": rows}, open(OUT + "/bar_targets.json", "w"), indent=1,
          default=lambda x: None if isinstance(x, float) and not math.isfinite(x) else (bool(x) if isinstance(x, np.bool_) else float(x)))
print("hecho")
