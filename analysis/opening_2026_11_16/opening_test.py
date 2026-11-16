#!/usr/bin/python3 -I
# opening_test.py — seguimiento 2026-11-16: ¿la pérdida de luz extra de octubre en hod2018 viene de la apertura de las cajas entre campañas?
# Solo lectura de salidas derivadas: impactos de PROC (hd64_may_v8_hits.bin con el tiempo bueno de v8, hd64_oct8_hits.bin = run 8, los mismos
# eventos que la tarea de la ventana CFD), G del análisis (calibración de outputs/hod2018_ampcal/2026-11-15_may_oct, 447105d: mayo = reajuste de PROC,
# octubre = c1, mediana para los canales sin G), orden de tiras (stage3.json), extremo de lectura (layer_from_channels.json), barras emparejadas por
# MPPC y grupos (seguimiento 7, ratio_by_mppc.csv). Etiquetas de los números en el README: [medido], [calculado], [supuesto].
# A. L_eff: método y cortes de bar_targets.py (seguimiento 10) sin cambios: por barra, intervalos de distancia al extremo de lectura [21,29),
#    [29,37), [37,45), [45,54]; etiqueta = tira ortogonal del intervalo > 10 p.e.; sonda > 5 p.e. medidos; mismo − mezclados; mediana interpolada en
#    bins de 1 p.e. medido, error por remuestreo de Poisson (200, semilla 1101 como bar_targets); intervalos fiables n_net ≥ 15; ajuste ponderado
#    de ln(mediana) contra la distancia; «determinada» con ≥ 3 intervalos fiables y pendiente < 0. Valor común = 1/⟨1/L⟩: media ponderada (1/σ²)
#    de las pendientes de las barras determinadas, L = −1/k̄, intervalo de L desde k̄ ± σ_k̄ (así sale 82.4 +11.1 −8.7 mm de mayo); por plano y
#    combinado. Comprobación: las pendientes por barra deben coincidir con bar_targets.json.
#    I0 (luz en el extremo de lectura, d = 0): con k̄ fijo, por barra ln I0 = media ponderada de ln(m_i) − k̄·d_i sobre sus intervalos fiables;
#    I0 del plano = media ponderada de ln I0 de las barras determinadas; el error suma en cuadratura el estadístico y el de k̄ (σ_k̄ · d̄ ponderada)
#    [calculado]. En p.e. medidos y en primarios con F = 1.75 / 1.46 [supuesto, el de bar_targets].
# B. Razón mayo/octubre contra la posición de cruce: etiqueta = una tira ortogonal j (> 10 p.e.), d_j = su distancia al extremo de lectura de la capa
#    sonda; sonda = toda tira viva del plano (> 3 p.e. medidos, selección de ampcal), mismo − mezclados, sumadas las dos capas del plano; mediana
#    interpolada en bins de 0.05 p.e. medidos; error por remuestreo de Poisson (200). Razón en primarios = razón medida × F_oct/F_mayo
#    (1.46/1.75 [supuesto]). Ajuste ponderado de ln(razón) = a + b·(d − 37.5 mm); χ² de una constante.
# C. Razón por canal (barra = MPPC, emparejada por posición de tira): etiqueta = cualquier tira ortogonal > 10 p.e., sonda > 3 p.e.; mediana,
#    error bootstrap; χ²/ndf de la constante (media ponderada), canales a > 3σ, asimetría (skewness) con error bootstrap de barras, orden físico
#    (índice de tira en el plano; grupo de 8 canales del VX1740 y grupo HV), compatibilidad dentro de cada grupo (χ² de la constante del grupo) y
#    entre grupos (χ² de las medias de grupo).
# D. Alineación de planos: matriz de coincidencias X-barra × Y-barra (las dos > 10 p.e. en la ventana signal), menos mezclados (X del evento i con Y del
#    i + 1, × N/(N − 1)), normalizada al total; por barra X el centroide del índice Y (y al revés), error = desviación/√n_neto; desplazamiento
#    octubre − mayo en unidades de barra; χ² de las dos matrices normalizadas al mismo total (celdas con contenido esperado > 0).
# Salida: <out>/A_leff_bars.csv, A_leff.csv, B_position.csv, B_fit.csv, C_channels.csv, C_groups.csv, D_matrix_{may,oct}.csv, D_centroids.csv,
#         opening_test.json. Uso: opening_test.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
BT = REPO + "/outputs/hod2018_bar_targets/2026-11-01_g4/bar_targets.json"
AMP = REPO + "/outputs/hod2018_ampcal/2026-11-15_may_oct/ampcal_gains.csv"
R7 = REPO + "/outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv"
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
TAG, PROBE_A, PROBE_B, N_MIN, NBOOT = 10.0, 5.0, 3.0, 15, 200
FCT = {"may": 1.75, "oct": 1.46}
EDG1 = np.arange(0, 301, 1.0)
EDG5 = np.arange(0, 300.0001, 0.05)


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))  # noqa: E731
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1())}
# la G debe ser la de ampcal (447105d)
for r in csv.DictReader(open(AMP)):
    if r["gain_set"] in ("refit_PROC (analysis)", "c1 (analysis)"):
        camp = r["campaign"]
        assert abs(G[camp][int(r["ch"])] - float(r["G_mV"])) < 6e-6, ("G distinta de ampcal", r)
e8 = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
DATA = {"may": e8[mc.good_time(e8["t"], "v8")], "oct": xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)}
PE = {c: DATA[c]["aS"][:, :32] / G[c][None, :32] for c in DATA}
SC = {c: len(DATA[c]) / (len(DATA[c]) - 1) for c in DATA}
print("eventos", {c: len(v) for c, v in DATA.items()})


def readout_sign(probe, layer):
    s = SIGN[f"oct2026_run8/{probe}{layer}"]["sign"]
    assert s == SIGN[f"may2026_v8_goodtime/{probe}{layer}"]["sign"]
    return s


def dist(probe, layer, j):
    return D[j] if readout_sign(probe, layer) < 0 else 75.0 - D[j]


def bin_js(probe, layer, k):
    lo, hi = BINS[k]
    js = [j for j in range(16) if lo <= dist(probe, layer, j) < hi + (1e-9 if hi == 54 else 0)]
    return [j for j in js if not (probe == "Y" and j == DEAD_X)]


def median_interp(net, edg):
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    before = cs[b - 1] if b > 0 else 0.0
    return float(edg[b] + (0.5 * n - before) / net[b] * (edg[b + 1] - edg[b]))


def med_boot(hs, hm, s, edg, rng):
    med = median_interp(np.clip(hs - s * hm, 0, None), edg)
    bs = [median_interp(np.clip(rng.poisson(hs) - s * rng.poisson(hm), 0, None), edg) for _ in range(NBOOT)]
    bs = [x for x in bs if np.isfinite(x)]
    return med, (float(np.std(bs)) if len(bs) > 1 else np.nan), float(hs.sum() - s * hm.sum())


def tagmask(camp, ref, js):
    o = ORDER[camp]
    refch = np.asarray(o[f"h18_{ref}"])
    t = PE[camp][:, refch] > TAG
    t[:, refch == DEAD] = False
    T = t[:, js].any(1)
    return T, np.concatenate([T[1:], [False]])


# ---------------------------------------------------------------- A
rngA = np.random.default_rng(1101)
bt = json.load(open(BT))["fits"]
arows, slopes = [], {}
for camp in ("may", "oct"):
    o = ORDER[camp]
    for probe, ref in (("X", "Y"), ("Y", "X")):
        for i, ch in enumerate(o[f"h18_{probe}"]):
            if ch == DEAD:
                arows.append([camp, probe, i, ch, "dead", "", "", "", "", "", ""])
                continue
            L = "A" if i % 2 == 0 else "B"
            ps = PE[camp][:, ch]
            dpos, dm, de = [], [], []
            for k in range(4):
                js = bin_js(probe, L, k)
                T, Tn = tagmask(camp, ref, js)
                sel = ps > PROBE_A
                hs = np.histogram(ps[sel & T], bins=EDG1)[0].astype(float)
                hm = np.histogram(ps[sel & Tn], bins=EDG1)[0].astype(float)
                med, err, nn = med_boot(hs, hm, SC[camp], EDG1, rngA)
                if nn >= N_MIN and np.isfinite(err):
                    dpos.append(float(np.mean([dist(probe, L, j) for j in js])))
                    dm.append(med)
                    de.append(err)
            d_, m_, e_ = (np.asarray(v, float) for v in (dpos, dm, de))
            ok = np.isfinite(m_) & np.isfinite(e_) & (e_ > 0) & (m_ > 0)
            d_, m_, e_ = d_[ok], m_[ok], e_[ok]
            ref_bt = bt[f"{camp}_{probe}{i:02d}"]
            if len(d_) < 3:
                arows.append([camp, probe, i, ch, f"not determined: {len(d_)} reliable bins", "", "", "", "", ref_bt["status"], ""])
                continue
            y, w = np.log(m_), (m_ / e_) ** 2
            X = np.vstack([np.ones_like(d_), d_ - 37.5]).T
            C = np.linalg.inv((X * w[:, None]).T @ X)
            b = C @ (X * w[:, None]).T @ y
            k, ek = float(b[1]), math.sqrt(C[1, 1])
            same = abs(k - ref_bt.get("slope_per_mm", np.nan)) < 1e-9
            status = "ok" if k < 0 else "not determined: slope >= 0"
            arows.append([camp, probe, i, ch, status, f"{k:.6g}", f"{ek:.6g}", f"{-1 / k:.4g}" if k < 0 else "", f"{ek / k ** 2:.4g}" if k < 0 else "",
                          ref_bt["status"], "1" if same else "0"])
            if k < 0:
                slopes[(camp, probe, i)] = (k, ek, d_, y, w, ch)
nmismatch = sum(1 for r in arows if r[10] == "0")
print("A: pendientes distintas de bar_targets.json:", nmismatch)


def combine(keys):
    k = np.array([slopes[x][0] for x in keys])
    e = np.array([slopes[x][1] for x in keys])
    W = 1 / e ** 2
    kb = float((W * k).sum() / W.sum())
    ekb = float(1 / math.sqrt(W.sum()))
    chi = float((W * (k - kb) ** 2).sum())
    L = -1 / kb
    Lhi = -1 / (kb + ekb) if kb + ekb < 0 else float("inf")
    Llo = -1 / (kb - ekb)
    # I0 con k̄ fijo
    lnI, elnI, dbar = [], [], []
    for x in keys:
        _, _, d_, y, w, _ = slopes[x]
        li = y - kb * d_
        lnI.append(float((w * li).sum() / w.sum()))
        elnI.append(float(1 / math.sqrt(w.sum())))
        dbar.append(float((w * d_).sum() / w.sum()))
    lnI, elnI, dbar = map(np.array, (lnI, elnI, dbar))
    Wi = 1 / elnI ** 2
    l0 = float((Wi * lnI).sum() / Wi.sum())
    el0_stat = float(1 / math.sqrt(Wi.sum()))
    el0_k = float(ekb * (Wi * dbar).sum() / Wi.sum())
    el0 = math.hypot(el0_stat, el0_k)
    return {"n": len(keys), "k": kb, "k_err": ekb, "chi2": chi, "ndf": len(keys) - 1, "L": L, "L_plus": Lhi - L, "L_minus": L - Llo,
            "lnI0": l0, "lnI0_err": el0, "lnI0_err_stat": el0_stat, "lnI0_err_k": el0_k, "I0_meas": math.exp(l0), "I0_meas_err": math.exp(l0) * el0}


Ares = {}
for camp in ("may", "oct"):
    for pl in ("X", "Y", "XY"):
        keys = [x for x in slopes if x[0] == camp and x[1] in pl]
        Ares[(camp, pl)] = combine(keys)
with open(OUT + "/A_leff_bars.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "plane", "strip_index", "channel", "status", "slope_per_mm", "slope_err", "L_eff_mm", "L_eff_err_mm", "bar_targets_status",
                "slope_equal_to_bar_targets"])
    w.writerows(arows)
with open(OUT + "/A_leff.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "plane", "n_bars", "k_mean_per_mm", "k_err", "chi2", "ndf", "L_eff_mm", "L_plus_mm", "L_minus_mm", "I0_meas_pe", "I0_meas_err",
                "I0_prim_pe", "I0_prim_err", "lnI0_err_stat", "lnI0_err_from_k"])
    for (camp, pl), a in Ares.items():
        w.writerow([camp, pl, a["n"], f"{a['k']:.6f}", f"{a['k_err']:.6f}", f"{a['chi2']:.2f}", a["ndf"], f"{a['L']:.2f}", f"{a['L_plus']:.2f}",
                    f"{a['L_minus']:.2f}", f"{a['I0_meas']:.3f}", f"{a['I0_meas_err']:.3f}", f"{a['I0_meas'] / FCT[camp]:.3f}",
                    f"{a['I0_meas_err'] / FCT[camp]:.3f}", f"{a['lnI0_err_stat']:.4f}", f"{a['lnI0_err_k']:.4f}"])
Aratio = {}
for pl in ("X", "Y", "XY"):
    m, o_ = Ares[("may", pl)], Ares[("oct", pl)]
    # razón de L por las pendientes: L_o/L_m = k_m/k_o
    rL = m["k"] / o_["k"]
    erL = rL * math.hypot(m["k_err"] / m["k"], o_["k_err"] / o_["k"])
    rI = o_["I0_meas"] / m["I0_meas"]
    erI = rI * math.hypot(m["lnI0_err"], o_["lnI0_err"])
    Aratio[pl] = {"L_oct_over_may": rL, "err": erL, "I0_meas_oct_over_may": rI, "err_I0": erI, "I0_prim_oct_over_may": rI * FCT["may"] / FCT["oct"],
                  "err_I0_prim": erI * FCT["may"] / FCT["oct"]}
    print(f"A {pl}: mayo L {m['L']:.1f} +{m['L_plus']:.1f} −{m['L_minus']:.1f} ({m['n']}), octubre L {o_['L']:.1f} +{o_['L_plus']:.1f} −{o_['L_minus']:.1f} "
          f"({o_['n']}); L_o/L_m {rL:.3f} ± {erL:.3f}; I0 medido {m['I0_meas']:.2f} ± {m['I0_meas_err']:.2f} / {o_['I0_meas']:.2f} ± {o_['I0_meas_err']:.2f}, "
          f"o/m {rI:.3f} ± {erI:.3f}", flush=True)

# ---------------------------------------------------------------- B
rngB = np.random.default_rng(1116)
brows, bfit = [], {}
FR = FCT["oct"] / FCT["may"]
for probe, ref in (("X", "Y"), ("Y", "X")):
    pts = []
    for dval in sorted(set(np.round(D, 3))):
        meds = {}
        for camp in ("may", "oct"):
            o = ORDER[camp]
            hs, hm = np.zeros(len(EDG5) - 1), np.zeros(len(EDG5) - 1)
            nb = 0
            for i, ch in enumerate(o[f"h18_{probe}"]):
                if ch == DEAD:
                    continue
                L = "A" if i % 2 == 0 else "B"
                j = [jj for jj in range(16) if abs(dist(probe, L, jj) - dval) < 1e-6]
                if not j or (probe == "Y" and j[0] == DEAD_X):
                    continue
                T, Tn = tagmask(camp, ref, j)
                ps = PE[camp][:, ch]
                sel = ps > PROBE_B
                hs += np.histogram(ps[sel & T], bins=EDG5)[0]
                hm += np.histogram(ps[sel & Tn], bins=EDG5)[0]
                nb += 1
            meds[camp] = med_boot(hs, hm, SC[camp], EDG5, rngB) + (nb,)
        (mm, em, nm, _), (mo, eo, no, _) = meds["may"], meds["oct"]
        r = mm / mo * FR
        er = r * math.hypot(em / mm, eo / mo)
        pts.append((dval, r, er))
        brows.append([probe, f"{dval:.1f}", meds["may"][3], f"{mm:.4f}", f"{em:.4f}", f"{nm:.1f}", f"{mo:.4f}", f"{eo:.4f}", f"{no:.1f}", f"{r:.4f}", f"{er:.4f}"])
    d_ = np.array([p[0] for p in pts])
    y = np.log([p[1] for p in pts])
    ey = np.array([p[2] / p[1] for p in pts])
    w = 1 / ey ** 2
    X = np.vstack([np.ones_like(d_), d_ - 37.5]).T
    C = np.linalg.inv((X * w[:, None]).T @ X)
    b = C @ (X * w[:, None]).T @ y
    chi = float((w * (y - X @ b) ** 2).sum())
    c0 = float((w * y).sum() / w.sum())
    chi0 = float((w * (y - c0) ** 2).sum())
    bfit[probe] = {"intercept_ln_at_37p5mm": float(b[0]), "intercept_err": math.sqrt(C[0, 0]), "slope_per_mm": float(b[1]), "slope_err": math.sqrt(C[1, 1]),
                   "intercept_ln_at_0mm": float(b[0] - 37.5 * b[1]), "intercept_0_err": math.sqrt(C[0, 0] + 37.5 ** 2 * C[1, 1] - 2 * 37.5 * C[0, 1]),
                   "chi2_line": chi, "ndf_line": len(d_) - 2, "const_ln": c0, "const_ratio": math.exp(c0), "const_err_ln": 1 / math.sqrt(w.sum()),
                   "chi2_const": chi0, "ndf_const": len(d_) - 1}
    print(f"B {probe}: pendiente {b[1]:.5f} ± {math.sqrt(C[1, 1]):.5f} /mm, ordenada (37.5 mm) {b[0]:.4f} ± {math.sqrt(C[0, 0]):.4f}, χ² recta {chi:.2f}/{len(d_) - 2}, "
          f"constante {math.exp(c0):.4f}, χ² {chi0:.2f}/{len(d_) - 1}", flush=True)
with open(OUT + "/B_position.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["plane", "dist_readout_mm", "n_probe_bars", "may_median_meas_pe", "err", "may_n_net", "oct_median_meas_pe", "err", "oct_n_net",
                "ratio_may_over_oct_prim", "err"])
    w.writerows(brows)
with open(OUT + "/B_fit.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    keys = list(next(iter(bfit.values())).keys())
    w.writerow(["plane"] + keys)
    for p, v in bfit.items():
        w.writerow([p] + [f"{v[k]:.6g}" if isinstance(v[k], float) else v[k] for k in keys])

# ---------------------------------------------------------------- C
rngC = np.random.default_rng(1117)
bars = [r for r in csv.DictReader(open(R7))]
crows = []
for b_ in bars:
    if b_["plane"] not in ("X", "Y") or int(b_["ch_may"]) >= 32:
        continue
    cm, co = int(b_["ch_may"]), int(b_["ch_oct"])
    row = {"id": b_["id"], "plane": b_["plane"], "index": int(b_["index"]), "layer": b_["layer"], "ch_may": cm, "ch_oct": co, "hv_group": int(b_["hv_group"]),
           "adc_group_may": int(b_["adc_group_may"]), "adc_group_oct": int(b_["adc_group_oct"])}
    if DEAD in (cm, co):
        row["note"] = "dead (ch3), excluded"
        crows.append(row)
        continue
    ref = "Y" if b_["plane"] == "X" else "X"
    med = {}
    for camp, ch in (("may", cm), ("oct", co)):
        js = [j for j in range(16) if not (ref == "X" and j == DEAD_X)]
        T, Tn = tagmask(camp, ref, js)
        ps = PE[camp][:, ch]
        sel = ps > PROBE_B
        hs = np.histogram(ps[sel & T], bins=EDG5)[0].astype(float)
        hm = np.histogram(ps[sel & Tn], bins=EDG5)[0].astype(float)
        med[camp] = med_boot(hs, hm, SC[camp], EDG5, rngC)
    r = med["may"][0] / med["oct"][0] * FR
    er = r * math.hypot(med["may"][1] / med["may"][0], med["oct"][1] / med["oct"][0])
    row.update(may_med=med["may"][0], may_err=med["may"][1], may_n=med["may"][2], oct_med=med["oct"][0], oct_err=med["oct"][1], oct_n=med["oct"][2],
               ratio=r, ratio_err=er, note="")
    crows.append(row)
live = [r for r in crows if r["note"] == ""]
R = np.array([r["ratio"] for r in live])
E = np.array([r["ratio_err"] for r in live])
W = 1 / E ** 2
mean = float((W * R).sum() / W.sum())
emean = float(1 / math.sqrt(W.sum()))
chi = float((W * (R - mean) ** 2).sum())
for r in live:
    r["pull"] = (r["ratio"] - mean) / r["ratio_err"]


def skew(x):
    x = np.asarray(x, float)
    m = x.mean()
    s = x.std()
    n = len(x)
    g1 = ((x - m) ** 3).mean() / s ** 3
    return float(g1 * math.sqrt(n * (n - 1)) / (n - 2))  # sesgo ajustado (G1)


sk = skew(R)
sk_bs = [skew(R[rngC.integers(0, len(R), len(R))]) for _ in range(1000)]
# sesgo esperado de una muestra normal con estos errores (simulación)
sk_null = [skew(mean + rngC.normal(0, E)) for _ in range(2000)]
grp = []
for kind in ("adc_group_oct", "adc_group_may", "hv_group"):
    gm, ge = [], []
    for gk in sorted({r[kind] for r in live}):
        rr = [r for r in live if r[kind] == gk]
        Rg = np.array([r["ratio"] for r in rr])
        Wg = 1 / np.array([r["ratio_err"] for r in rr]) ** 2
        mg = float((Wg * Rg).sum() / Wg.sum())
        cg = float((Wg * (Rg - mg) ** 2).sum())
        gm.append(mg)
        ge.append(1 / math.sqrt(Wg.sum()))
        grp.append([kind, gk, len(rr), f"{mg:.4f}", f"{1 / math.sqrt(Wg.sum()):.4f}", f"{cg:.2f}", len(rr) - 1, " ".join(r["id"] for r in rr)])
    gm, ge = np.array(gm), np.array(ge)
    wg = 1 / ge ** 2
    mm_ = float((wg * gm).sum() / wg.sum())
    grp.append([kind, "between groups", len(gm), f"{mm_:.4f}", f"{1 / math.sqrt(wg.sum()):.4f}", f"{float((wg * (gm - mm_) ** 2).sum()):.2f}", len(gm) - 1, ""])
CK = ["id", "plane", "index", "layer", "ch_may", "ch_oct", "adc_group_may", "adc_group_oct", "hv_group", "may_med", "may_err", "may_n", "oct_med", "oct_err",
      "oct_n", "ratio", "ratio_err", "pull", "note"]
with open(OUT + "/C_channels.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(CK)
    for r in sorted(crows, key=lambda r: (r["plane"], r["index"])):
        w.writerow([f"{r[k]:.4f}" if isinstance(r.get(k), float) else r.get(k, "") for k in CK])
with open(OUT + "/C_groups.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["grouping", "group", "n_bars", "weighted_mean_ratio", "err", "chi2_const", "ndf", "bars"])
    w.writerows(grp)
Cres = {"n": len(live), "mean": mean, "mean_err": emean, "chi2": chi, "ndf": len(live) - 1, "outliers_3sigma": [r["id"] for r in live if abs(r["pull"]) > 3],
        "skewness_G1": sk, "skewness_bootstrap_err": float(np.std(sk_bs)), "skewness_null_p16_p84": [float(np.percentile(sk_null, 16)), float(np.percentile(sk_null, 84))],
        "skewness_null_p_two_sided": float(np.mean(np.abs(sk_null) >= abs(sk))), "median_ratio": float(np.median(R))}
print("C:", Cres, flush=True)
for g in grp:
    print("C grupo", g[:7], flush=True)

# ---------------------------------------------------------------- D
Dres, drows = {}, []
for camp in ("may", "oct"):
    o = ORDER[camp]
    xc, yc = np.asarray(o["h18_X"]), np.asarray(o["h18_Y"])
    HX = PE[camp][:, xc] > TAG
    HY = PE[camp][:, yc] > TAG
    HX[:, xc == DEAD] = False
    same = HX.T.astype(float) @ HY.astype(float)
    mix = HX[:-1].T.astype(float) @ HY[1:].astype(float)
    s = SC[camp]
    net = same - s * mix
    var = same + s * s * mix
    Dres[camp] = (net, var)
    with open(OUT + f"/D_matrix_{camp}.csv", "w") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["X_index\\Y_index"] + list(range(16)))
        tot = net.sum()
        for i in range(16):
            w.writerow([i] + [f"{net[i, j] / tot:.6f}" for j in range(16)])
cent = {}
for camp in ("may", "oct"):
    net, var = Dres[camp]
    for axis, lab in ((1, "X"), (0, "Y")):
        for i in range(16):
            row = net[i, :] if axis == 1 else net[:, i]
            vr = var[i, :] if axis == 1 else var[:, i]
            idx = np.arange(16)
            n = row.sum()
            if n <= 5:
                cent[(camp, lab, i)] = (np.nan, np.nan, n)
                continue
            c = float((row * idx).sum() / n)
            # error por propagación de Poisson de cada celda
            ec = float(math.sqrt(((idx - c) ** 2 * vr).sum()) / n)
            cent[(camp, lab, i)] = (c, ec, n)
for lab in ("X", "Y"):
    for i in range(16):
        cm_, em_, nm_ = cent[("may", lab, i)]
        co_, eo_, no_ = cent[("oct", lab, i)]
        sh = co_ - cm_
        esh = math.hypot(em_, eo_) if np.isfinite(em_) and np.isfinite(eo_) else np.nan
        drows.append([lab, i, f"{nm_:.1f}", f"{cm_:.4f}", f"{em_:.4f}", f"{no_:.1f}", f"{co_:.4f}", f"{eo_:.4f}", f"{sh:.4f}", f"{esh:.4f}",
                      f"{sh / esh:.2f}" if np.isfinite(esh) and esh > 0 else ""])
(nm, vm), (no, vo) = Dres["may"], Dres["oct"]
k = no.sum() / nm.sum()
ok = (vm + vo) > 0
den = np.where(ok, vo + k * k * vm, 1.0)
chi2D = float((((no - k * nm) ** 2) / den)[ok & ((no > 0) | (nm > 0))].sum())
ndfD = int((ok & ((no > 0) | (nm > 0))).sum()) - 1
sh = np.array([float(r[8]) for r in drows if r[9] not in ("nan", "")])
esh = np.array([float(r[9]) for r in drows if r[9] not in ("nan", "")])
lab = [r[0] for r in drows if r[9] not in ("nan", "")]
shift_mean = {}
for L in ("X", "Y"):
    m_ = np.array([l == L for l in lab])
    w = 1 / esh[m_] ** 2
    mu = float((w * sh[m_]).sum() / w.sum())
    shift_mean[L] = {"mean_shift_bars": mu, "err": float(1 / math.sqrt(w.sum())), "chi2_const": float((w * (sh[m_] - mu) ** 2).sum()), "ndf": int(m_.sum() - 1),
                     "chi2_zero": float((w * sh[m_] ** 2).sum()), "ndf_zero": int(m_.sum())}
with open(OUT + "/D_centroids.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["bar_plane", "bar_index", "may_n_net", "may_centroid_other_index", "err", "oct_n_net", "oct_centroid_other_index", "err", "oct_minus_may_bars", "err",
                "pull"])
    w.writerows(drows)
Dsum = {"chi2_matrices": chi2D, "ndf": ndfD, "n_net_may": float(nm.sum()), "n_net_oct": float(no.sum()), "shift_mean": shift_mean}
print("D:", Dsum, flush=True)


def gammq(a, x):
    """Q(a, x) regularizada: probabilidad de χ² ≥ x con 2a grados de libertad."""
    if x <= 0:
        return 1.0
    if x < a + 1:
        s = t = 1.0 / a
        ap = a
        for _ in range(100000):
            ap += 1
            t *= x / ap
            s += t
            if abs(t) < abs(s) * 1e-14:
                break
        return 1 - s * math.exp(-x + a * math.log(x) - math.lgamma(a))
    b, c, d = x + 1 - a, 1e300, 1 / (x + 1 - a)
    h = d
    for i in range(1, 100000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1 / d if abs(d) > 1e-300 else 1e300
        c = b + an / c
        h *= d * c
        if abs(d * c - 1) < 1e-14:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


PV = {f"B_{p}_const": gammq(bfit[p]["ndf_const"] / 2, bfit[p]["chi2_const"] / 2) for p in bfit}
PV.update({f"B_{p}_line": gammq(bfit[p]["ndf_line"] / 2, bfit[p]["chi2_line"] / 2) for p in bfit})
PV["C_const"] = gammq(Cres["ndf"] / 2, Cres["chi2"] / 2)
for g in grp:
    PV[f"C_{g[0]}_{g[1]}"] = gammq(int(g[6]) / 2, float(g[5]) / 2) if int(g[6]) > 0 else float("nan")
PV["D_matrices"] = gammq(ndfD / 2, chi2D / 2)
for L in ("X", "Y"):
    PV[f"D_{L}_shift_const"] = gammq(shift_mean[L]["ndf"] / 2, shift_mean[L]["chi2_const"] / 2)
    PV[f"D_{L}_shift_zero"] = gammq(shift_mean[L]["ndf_zero"] / 2, shift_mean[L]["chi2_zero"] / 2)
for (camp, pl), a in Ares.items():
    PV[f"A_{camp}_{pl}_common_slope"] = gammq(a["ndf"] / 2, a["chi2"] / 2)
print("p-valores:", {k: round(v, 6) for k, v in PV.items()}, flush=True)
json.dump({"p_values": PV, "A": {f"{k[0]}_{k[1]}": v for k, v in Ares.items()}, "A_ratio": Aratio, "A_slopes_differ_from_bar_targets": nmismatch, "B": bfit, "C": Cres,
           "D": Dsum, "F_assumed": FCT, "events": {c: len(v) for c, v in DATA.items()}}, open(OUT + "/opening_test.json", "w"), indent=1)
