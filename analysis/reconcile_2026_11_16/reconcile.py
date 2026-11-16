#!/usr/bin/python3 -I
# reconcile.py — seguimiento 2026-11-16: reconciliación de las razones mayo/octubre de hod2018, prueba de un intercambio Y07/Y08, residuo de PDE y
# tasas oscuras. Solo lectura de salidas derivadas (y de work/*_dark4.bin de dark4_extract.C).
# Entradas: impactos de PROC (hd64_may_v8_hits.bin con tiempo bueno, hd64_oct8_hits.bin), G del análisis (ampcal 447105d), stage3.json,
# layer_from_channels.json, ratio_by_mppc.csv (barra ↔ canal por campaña, grupo HV), hv_delta.csv (ΔV), muon_refit.json (p de los muones),
# xtalk_ratio.csv (razones primarias por definición de F), seguimiento 12 oct8_G_spec.csv (G horaria).
# 1. B y C del test de apertura con el estimador de 0.6477 (ampcal): mediana interpolada del histograma neto (mismo − mezclados) en bins de 0.05 p.e.
#    medidos, sonda > 3 p.e. medidos, en p.e. MEDIDOS (sin F) y como octubre/mayo. Variantes: sonda > 5 p.e. (corte de bar_targets) y, en C, solo
#    etiquetas del intervalo 29-54 mm (la selección exacta de 0.6477). Comprobación: la mediana de la estación del intervalo lejos debe dar 16.1627 /
#    10.4680.
# 2. Y06-Y09 (y todo par vecino de Y con desvíos opuestos > 2σ en C): G de cada canal, tasa > 0.5 p.e., columna de la matriz X × Y (distribución sobre
#    las barras X de las coincidencias de la barra Y, ambas > 10 p.e., menos mezclados), comparada entre mayo y octubre con χ² de formas; C con los
#    canales de octubre del par intercambiados.
# 4. Residuo = razón primaria oct/mayo por canal (C de 1, con F = 1/(1 − p) de los muones, muon_refit.json ajuste R) / razón de PDE;
#    PDE(OV) de la hoja de S12572-100P a 408 nm [hoja de datos], interpolación lineal [S], sin extrapolar bajo 0.6 V; OV_mayo = 1.0, 1.2, 1.4 V [S],
#    OV_oct = OV_mayo · G_oct/G_mayo [C]; MPPC con G de octubre rellenada (sin c1) fuera. Por canal, grupo HV y plano: media ponderada, χ² de la
#    constante, pendiente contra ΔV (grupos). R_prim/R_G por definición de F: cota superior del factor de luz extra para toda PDE(OV) cóncava con
#    PDE(0) = 0 (entonces PDE(x·OV)/PDE(OV) ≥ x para x < 1).
# 5. Tasa > 0.5 y > 1.5 p.e. en [16,37) por canal (mayo tiempo bueno, octubre run 8), octubre/mayo por MPPC contra G_oct/G_mayo; octubre por hora:
#    ln R_h regresado sobre ln G_h (G horaria del seguimiento 12, fuente spec, criterio de diurnal.py) por canal; residuo por hora (mediana de la
#    estación) contra la hora local (UTC−3), día 07-19 h contra noche [S: luces de la sala encendidas de día].
# Salida en <out>: r1_B.csv, r1_C.csv, r1_summary.json, r2_pairs.csv, r2_swap.json, r4_residual_channels.csv, r4_residual_groups.csv, r4_summary.json,
#   r4_bound.csv, r5_dark_channels.csv, r5_hourly.csv, r5_summary.json. Uso: reconcile.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
WK = OUT + "/work"
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
R7 = REPO + "/outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv"
HVD = REPO + "/outputs/hod2018_hv_may_oct/2026-11-15/hv_delta.csv"
XT = REPO + "/outputs/hod2018_xtalk/2026-11-15_borel"
F12 = REPO + "/outputs/hod_stability/2026-11-03_followup12/work"
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
NBOOT = 200
EDG5 = np.arange(0, 300.0001, 0.05)


def gammq(a, x):
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


def pchi(c, n):
    return gammq(n / 2, c / 2) if n > 0 else float("nan")


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))  # noqa: E731
C1 = mc.g_c1()[:32]
G = {"may": fill(g_refit()), "oct": fill(C1)}
OCT_FILLED = [c for c in range(32) if not np.isfinite(C1[c]) and c != DEAD]
e8 = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
GOOD = mc.good_time(e8["t"], "v8")
DATA = {"may": e8[GOOD], "oct": xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)}
PE = {c: DATA[c]["aS"][:, :32] / G[c][None, :32] for c in DATA}
SC = {c: len(DATA[c]) / (len(DATA[c]) - 1) for c in DATA}


def readout_sign(probe, layer):
    return SIGN[f"oct2026_run8/{probe}{layer}"]["sign"]


def dist(probe, layer, j):
    return D[j] if readout_sign(probe, layer) < 0 else 75.0 - D[j]


def median_interp(net, edg=EDG5):
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    before = cs[b - 1] if b > 0 else 0.0
    return float(edg[b] + (0.5 * n - before) / net[b] * (edg[b + 1] - edg[b]))


def med_boot(hs, hm, s, rng):
    med = median_interp(np.clip(hs - s * hm, 0, None))
    bs = [median_interp(np.clip(rng.poisson(hs) - s * rng.poisson(hm), 0, None)) for _ in range(NBOOT)]
    bs = [x for x in bs if np.isfinite(x)]
    return med, float(np.std(bs)) if len(bs) > 1 else np.nan, float(hs.sum() - s * hm.sum())


def tagmask(camp, ref, js):
    refch = np.asarray(ORDER[camp][f"h18_{ref}"])
    t = PE[camp][:, refch] > 10.0
    t[:, refch == DEAD] = False
    T = t[:, js].any(1)
    return T, np.concatenate([T[1:], [False]])


def hist_probe(camp, ch, T, Tn, thr):
    ps = PE[camp][:, ch]
    sel = ps > thr
    return np.histogram(ps[sel & T], bins=EDG5)[0].astype(float), np.histogram(ps[sel & Tn], bins=EDG5)[0].astype(float)


def far_js(probe, L):
    js = [j for j in range(16) if 29 <= dist(probe, L, j) < 54 + 1e-9]
    return [j for j in js if not (probe == "Y" and j == DEAD_X)]


# ---------------------------------------------------------------- comprobación de 0.6477
chk = {}
for camp in ("may", "oct"):
    o = ORDER[camp]
    hs, hm = np.zeros(len(EDG5) - 1), np.zeros(len(EDG5) - 1)
    for probe, ref in (("X", "Y"), ("Y", "X")):
        for i, ch in enumerate(o[f"h18_{probe}"]):
            if ch == DEAD:
                continue
            L = "A" if i % 2 == 0 else "B"
            refch = np.asarray(o[f"h18_{ref}"])
            t = PE[camp][:, refch] > 10.0
            t[:, refch == DEAD] = False
            js = [j for j in range(16) if 29 <= dist(probe, L, j) < 54 + 1e-9]  # como ampcal: el corte X12 lo hace la máscara muerta
            T = t[:, js].any(1)
            a, b = hist_probe(camp, ch, T, np.concatenate([T[1:], [False]]), 3.0)
            hs += a
            hm += b
    chk[camp] = median_interp(np.clip(hs - SC[camp] * hm, 0, None))
print(f"comprobación de 0.6477: mediana lejos mayo {chk['may']:.4f} (ampcal 16.1627), octubre {chk['oct']:.4f} (ampcal 10.4680), razón {chk['oct'] / chk['may']:.4f}",
      flush=True)

# ---------------------------------------------------------------- 1. B y C en p.e. medidos, octubre/mayo
rng = np.random.default_rng(1601)
bars = [r for r in csv.DictReader(open(R7)) if r["plane"] in ("X", "Y") and int(r["ch_may"]) < 32]
BY = {r["id"]: r for r in bars}


def wmean(R, E):
    R, E = np.asarray(R, float), np.asarray(E, float)
    W = 1 / E ** 2
    m = float((W * R).sum() / W.sum())
    return m, float(1 / math.sqrt(W.sum())), float((W * (R - m) ** 2).sum()), len(R) - 1


def C_ratios(thr, far_only, oct_map=None):
    """Razón oct/mayo en p.e. medidos por barra. oct_map: id -> canal de octubre (para intercambios)."""
    out = {}
    for b in bars:
        cm = int(b["ch_may"])
        co = int(oct_map.get(b["id"], b["ch_oct"])) if oct_map else int(b["ch_oct"])
        if DEAD in (cm, co):
            continue
        probe = b["plane"]
        ref = "Y" if probe == "X" else "X"
        L = b["layer"]
        meds = {}
        for camp, ch in (("may", cm), ("oct", co)):
            js = far_js(probe, L) if far_only else [j for j in range(16) if not (ref == "X" and j == DEAD_X)]
            T, Tn = tagmask(camp, ref, js)
            hs, hm = hist_probe(camp, ch, T, Tn, thr)
            meds[camp] = med_boot(hs, hm, SC[camp], rng)
        r = meds["oct"][0] / meds["may"][0]
        er = r * math.hypot(meds["oct"][1] / meds["oct"][0], meds["may"][1] / meds["may"][0])
        out[b["id"]] = {"may": meds["may"], "oct": meds["oct"], "ratio": r, "err": er}
    return out


def C_summary(cr):
    ids = list(cr)
    m, e, chi, ndf = wmean([cr[i]["ratio"] for i in ids], [cr[i]["err"] for i in ids])
    pulls = {i: (cr[i]["ratio"] - m) / cr[i]["err"] for i in ids}
    return {"mean": m, "err": e, "chi2": chi, "ndf": ndf, "p": pchi(chi, ndf), "median": float(np.median([cr[i]["ratio"] for i in ids])),
            "outliers_3sigma": {i: round(p, 3) for i, p in pulls.items() if abs(p) > 3}, "pulls": pulls}


CV = {"thr3_alltags": C_ratios(3.0, False), "thr5_alltags": C_ratios(5.0, False), "thr3_far": C_ratios(3.0, True)}
CS = {k: C_summary(v) for k, v in CV.items()}
for k, s in CS.items():
    print(f"1 C {k}: media {s['mean']:.4f} ± {s['err']:.4f}, mediana {s['median']:.4f}, χ² {s['chi2']:.2f}/{s['ndf']} (p {s['p']:.3g}), >3σ {s['outliers_3sigma']}",
          flush=True)
with open(OUT + "/r1_C.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["variant", "id", "plane", "index", "ch_may", "ch_oct", "may_median_meas", "err", "may_n_net", "oct_median_meas", "err", "oct_n_net", "ratio_oct_over_may",
                "err", "pull"])
    for k, v in CV.items():
        for i, r in sorted(v.items(), key=lambda x: (BY[x[0]]["plane"], int(BY[x[0]]["index"]))):
            b = BY[i]
            w.writerow([k, i, b["plane"], b["index"], b["ch_may"], b["ch_oct"], f"{r['may'][0]:.4f}", f"{r['may'][1]:.4f}", f"{r['may'][2]:.1f}", f"{r['oct'][0]:.4f}",
                        f"{r['oct'][1]:.4f}", f"{r['oct'][2]:.1f}", f"{r['ratio']:.4f}", f"{r['err']:.4f}", f"{CS[k]['pulls'][i]:.3f}"])

# B en medidos, octubre/mayo
brows, bfit = [], {}
for thr in (3.0, 5.0):
    for probe, ref in (("X", "Y"), ("Y", "X")):
        pts = []
        for dval in sorted(set(np.round(D, 3))):
            meds = {}
            for camp in ("may", "oct"):
                o = ORDER[camp]
                hs, hm = np.zeros(len(EDG5) - 1), np.zeros(len(EDG5) - 1)
                for i, ch in enumerate(o[f"h18_{probe}"]):
                    if ch == DEAD:
                        continue
                    L = "A" if i % 2 == 0 else "B"
                    j = [jj for jj in range(16) if abs(dist(probe, L, jj) - dval) < 1e-6]
                    if not j or (probe == "Y" and j[0] == DEAD_X):
                        continue
                    T, Tn = tagmask(camp, ref, j)
                    a, b_ = hist_probe(camp, ch, T, Tn, thr)
                    hs += a
                    hm += b_
                meds[camp] = med_boot(hs, hm, SC[camp], rng)
            r = meds["oct"][0] / meds["may"][0]
            er = r * math.hypot(meds["oct"][1] / meds["oct"][0], meds["may"][1] / meds["may"][0])
            pts.append((dval, r, er))
            brows.append([f"thr{thr:g}", probe, f"{dval:.1f}", f"{meds['may'][0]:.4f}", f"{meds['may'][1]:.4f}", f"{meds['oct'][0]:.4f}", f"{meds['oct'][1]:.4f}", f"{r:.4f}",
                          f"{er:.4f}"])
        d_ = np.array([p[0] for p in pts])
        y = np.log([p[1] for p in pts])
        ey = np.array([p[2] / p[1] for p in pts])
        w_ = 1 / ey ** 2
        X = np.vstack([np.ones_like(d_), d_ - 37.5]).T
        Cm = np.linalg.inv((X * w_[:, None]).T @ X)
        bb = Cm @ (X * w_[:, None]).T @ y
        chi = float((w_ * (y - X @ bb) ** 2).sum())
        c0 = float((w_ * y).sum() / w_.sum())
        chi0 = float((w_ * (y - c0) ** 2).sum())
        bfit[f"thr{thr:g}_{probe}"] = {"slope_per_mm": float(bb[1]), "slope_err": math.sqrt(Cm[1, 1]), "ln_at_37p5": float(bb[0]), "ln_err": math.sqrt(Cm[0, 0]),
                                       "chi2_line": chi, "ndf_line": len(d_) - 2, "const_ratio": math.exp(c0), "const_err_rel": 1 / math.sqrt(w_.sum()),
                                       "chi2_const": chi0, "ndf_const": len(d_) - 1, "p_const": pchi(chi0, len(d_) - 1)}
        print(f"1 B thr{thr:g} {probe}: pendiente {bb[1]:.5f} ± {math.sqrt(Cm[1, 1]):.5f} /mm, constante {math.exp(c0):.4f}, χ² {chi0:.2f}/{len(d_) - 1} "
              f"(p {pchi(chi0, len(d_) - 1):.3g})", flush=True)
with open(OUT + "/r1_B.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["threshold", "plane", "dist_readout_mm", "may_median_meas", "err", "oct_median_meas", "err", "ratio_oct_over_may", "err"])
    w.writerows(brows)
json.dump({"check_0p6477": chk, "C": {k: {kk: vv for kk, vv in s.items() if kk != "pulls"} for k, s in CS.items()}, "B": bfit},
          open(OUT + "/r1_summary.json", "w"), indent=1)

# ---------------------------------------------------------------- 5 (primero: las tasas se usan en 2)
DT4 = np.dtype([("t", "<f8"), ("n05", "u1", (32,)), ("n15", "u1", (32,))])
DK = {}
for camp, pre in (("may", "may_v8"), ("oct", "oct_run8")):
    d = np.fromfile(f"{WK}/{pre}_dark4.bin", dtype=DT4)
    if camp == "may":
        d = d[mc.good_time(d["t"], "v8")]
    DK[camp] = d
EXPO = 21 * 16e-9  # s por evento y canal
RATE = {camp: {k: DK[camp][k].sum(0) / (len(DK[camp]) * EXPO) for k in ("n05", "n15")} for camp in DK}
RERR = {camp: {k: np.sqrt(DK[camp][k].sum(0)) / (len(DK[camp]) * EXPO) for k in ("n05", "n15")} for camp in DK}
drows = []
for b in sorted(bars, key=lambda b: (b["plane"], int(b["index"]))):
    cm, co = int(b["ch_may"]), int(b["ch_oct"])
    if DEAD in (cm, co):
        drows.append([b["id"], cm, co] + [""] * 11 + ["dead"])
        continue
    gr = G["oct"][co] / G["may"][cm]
    r05 = RATE["oct"]["n05"][co] / RATE["may"]["n05"][cm]
    r15 = RATE["oct"]["n15"][co] / RATE["may"]["n15"][cm]
    e05 = r05 * math.hypot(RERR["oct"]["n05"][co] / RATE["oct"]["n05"][co], RERR["may"]["n05"][cm] / RATE["may"]["n05"][cm])
    e15 = r15 * math.hypot(RERR["oct"]["n15"][co] / RATE["oct"]["n15"][co], RERR["may"]["n15"][cm] / RATE["may"]["n15"][cm])
    drows.append([b["id"], cm, co, f"{RATE['may']['n05'][cm] / 1e3:.2f}", f"{RATE['oct']['n05'][co] / 1e3:.2f}", f"{r05:.4f}", f"{e05:.4f}",
                  f"{RATE['may']['n15'][cm] / 1e3:.2f}", f"{RATE['oct']['n15'][co] / 1e3:.2f}", f"{r15:.4f}", f"{e15:.4f}", f"{G['may'][cm]:.5f}", f"{G['oct'][co]:.5f}",
                  f"{gr:.4f}", "oct G filled (no c1)" if co in OCT_FILLED else ""])
with open(OUT + "/r5_dark_channels.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "ch_may", "ch_oct", "rate05_may_kHz", "rate05_oct_kHz", "rate05_oct_over_may", "err", "rate15_may_kHz", "rate15_oct_kHz", "rate15_oct_over_may", "err",
                "G_may_mV", "G_oct_mV", "G_oct_over_may", "note"])
    w.writerows(drows)
ok = [r for r in drows if r[-1] == ""]
gq = np.array([float(r[13]) for r in ok])
r05 = np.array([float(r[5]) for r in ok])
r15 = np.array([float(r[9]) for r in ok])


def ols(x, y):
    X = np.vstack([np.ones_like(x), x]).T
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    s2 = float((res ** 2).sum() / (len(y) - 2))
    Cm = s2 * np.linalg.inv(X.T @ X)
    return float(b[0]), float(b[1]), math.sqrt(Cm[0, 0]), math.sqrt(Cm[1, 1]), float(np.corrcoef(x, y)[0, 1])


r5 = {"median_rate05_oct_over_may": float(np.median(r05)), "median_rate15_oct_over_may": float(np.median(r15)), "median_G_oct_over_may": float(np.median(gq)),
      "station_rate05_kHz": {c: float(np.median(np.delete(RATE[c]["n05"], DEAD))) / 1e3 for c in RATE},
      "station_rate15_kHz": {c: float(np.median(np.delete(RATE[c]["n15"], DEAD))) / 1e3 for c in RATE},
      "fit_ln_rate05_vs_ln_G": ols(np.log(gq), np.log(r05)), "fit_ln_rate15_vs_ln_G": ols(np.log(gq), np.log(r15)), "n_mppc": len(ok)}
# horario de octubre
hG = {}
seeds = [float(v) for v in open(F12 + "/oct8_seeds.txt").read().split(",")]
run_ok = {}
for r in csv.DictReader(l for l in open(F12 + "/oct8_G_spec.csv") if not l.startswith("#")):
    h, c = int(r["hour_utc_s"]), int(r["ch"])
    v = (float(r["G_mV"]), float(r["G_err_mV"]), int(r["fit_status"]))
    if h == -1:
        run_ok[c] = v[2] in (0, 1) and v[1] > 0
    elif v[2] in (0, 1) and v[1] > 0 and 0.62 * seeds[c] < v[0] < 1.38 * seeds[c] and v[1] / v[0] < 0.10:
        hG[(h, c)] = v[0]
d = DK["oct"]
hr = (np.floor(d["t"] / 3600) * 3600).astype(np.int64)
hours = sorted(set(hr.tolist()))
hrows, resid = [], {}
for c in range(32):
    if c == DEAD or not run_ok.get(c):
        continue
    xs, ys, hs_ = [], [], []
    for h in hours:
        m = hr == h
        n = int(m.sum())
        if n < 500 or (h, c) not in hG:
            continue
        R = d["n05"][m, c].sum() / (n * EXPO)
        xs.append(math.log(hG[(h, c)]))
        ys.append(math.log(R))
        hs_.append(h)
    if len(xs) < 6:
        continue
    xs, ys = np.array(xs), np.array(ys)
    a0, a1, _, _, _ = ols(xs, ys)
    for h, x_, y_ in zip(hs_, xs, ys):
        resid.setdefault(h, []).append(y_ - (a0 + a1 * x_))
    resid.setdefault("slopes", []).append(a1)
day = {"day": [], "night": []}
for h in hours:
    if h not in resid:
        continue
    v = np.array(resid[h])
    loc = ((h - 3 * 3600) % 86400) / 3600.0
    k = "day" if 7 <= loc < 19 else "night"
    med = float(np.median(v))
    err = float(1.2533 * np.std(v) / math.sqrt(len(v)))
    day[k].append((med, err))
    n = int((hr == h).sum())
    hrows.append([h, f"{loc:.0f}", k, n, len(v), f"{med:.5f}", f"{err:.5f}", f"{float(np.median([d['n05'][hr == h, c].sum() / (n * EXPO) for c in range(32) if c != DEAD])) / 1e3:.2f}"])
with open(OUT + "/r5_hourly.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hour_utc_s", "local_hour_utc-3", "day_night", "events", "n_channels", "median_resid_ln_rate05_fixed_G", "err", "median_rate05_kHz_raw"])
    w.writerows(hrows)
dn = {}
for k, v in day.items():
    if v:
        m_ = np.array([x[0] for x in v])
        e_ = np.array([x[1] for x in v])
        W = 1 / e_ ** 2
        mm = float((W * m_).sum() / W.sum())
        dn[k] = {"n_hours": len(v), "mean_resid": mm, "err": float(1 / math.sqrt(W.sum())), "chi2": float((W * (m_ - mm) ** 2).sum())}
if "day" in dn and "night" in dn:
    dn["day_minus_night"] = dn["day"]["mean_resid"] - dn["night"]["mean_resid"]
    dn["day_minus_night_err"] = math.hypot(dn["day"]["err"], dn["night"]["err"])
allr = np.array([x[0] for v in day.values() for x in v])
alle = np.array([x[1] for v in day.values() for x in v])
W = 1 / alle ** 2
m0 = float((W * allr).sum() / W.sum())
dn["chi2_const_all_hours"] = float((W * (allr - m0) ** 2).sum())
dn["ndf"] = len(allr) - 1
dn["p_const"] = pchi(dn["chi2_const_all_hours"], dn["ndf"])
dn["median_dlnR_dlnG"] = float(np.median(resid["slopes"]))
dn["n_channels"] = len(resid["slopes"])
r5["october_hourly"] = dn
json.dump(r5, open(OUT + "/r5_summary.json", "w"), indent=1)
print("5:", json.dumps(r5, indent=None), flush=True)

# ---------------------------------------------------------------- 2. Y06-Y09 y pares con desvíos opuestos
def coinc(camp):
    o = ORDER[camp]
    xc, yc = np.asarray(o["h18_X"]), np.asarray(o["h18_Y"])
    HX = PE[camp][:, xc] > 10.0
    HY = PE[camp][:, yc] > 10.0
    HX[:, xc == DEAD] = False
    same = HX.T.astype(float) @ HY.astype(float)
    mix = HX[:-1].T.astype(float) @ HY[1:].astype(float)
    return same - SC[camp] * mix, same + SC[camp] ** 2 * mix


CM = {c: coinc(c) for c in ("may", "oct")}


def shape_chi2(jm, jo):
    (nm, vm), (no, vo) = CM["may"], CM["oct"]
    a, va, b, vb = nm[:, jm], vm[:, jm], no[:, jo], vo[:, jo]
    k = b.sum() / a.sum()
    m = (va + vb) > 0
    m[DEAD_X] = False
    return float(((b - k * a) ** 2 / np.where(m, vb + k * k * va, 1))[m].sum()), int(m.sum() - 1)


ypulls = CS["thr3_alltags"]["pulls"]
pairs = [("Y06", "Y07"), ("Y07", "Y08"), ("Y08", "Y09")]
for i in range(15):
    a, b = f"Y{i:02d}", f"Y{i + 1:02d}"
    if a in ypulls and b in ypulls and ypulls[a] * ypulls[b] < 0 and abs(ypulls[a]) > 2 and abs(ypulls[b]) > 2 and (a, b) not in pairs:
        pairs.append((a, b))
pr, sw = [], {}
for bid in sorted({x for p in pairs for x in p}):
    b = BY[bid]
    j = int(b["index"])
    cm, co = int(b["ch_may"]), int(b["ch_oct"])
    row = [bid, j, cm, co, f"{G['may'][cm]:.5f}", f"{G['oct'][co]:.5f}", f"{G['oct'][co] / G['may'][cm]:.4f}", f"{RATE['may']['n05'][cm] / 1e3:.2f}",
           f"{RATE['oct']['n05'][co] / 1e3:.2f}", f"{RATE['oct']['n05'][co] / RATE['may']['n05'][cm]:.4f}"]
    for jj in (j - 1, j, j + 1):
        if 0 <= jj < 16:
            c2, n2 = shape_chi2(j, jj)
            row.append(f"Oct Y{jj:02d}: {c2:.2f}/{n2}")
        else:
            row.append("")
    pr.append(row)
with open(OUT + "/r2_pairs.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["bar", "index", "ch_may", "ch_oct", "G_may_mV", "G_oct_mV", "G_oct_over_may", "rate05_may_kHz", "rate05_oct_kHz", "rate05_oct_over_may",
                "matrix_column_chi2_May_bar_vs_Oct_index_minus1", "vs_Oct_same_index", "vs_Oct_index_plus1"])
    w.writerows(pr)
base = CS["thr3_alltags"]
for a, b in pairs[1:] if False else pairs:
    if a not in BY or b not in BY:
        continue
    m = {a: BY[b]["ch_oct"], b: BY[a]["ch_oct"]}
    cr = C_ratios(3.0, False, oct_map=m)
    s = C_summary(cr)
    ga, gb = BY[a], BY[b]
    gsw = {a: G["oct"][int(gb["ch_oct"])] / G["may"][int(ga["ch_may"])], b: G["oct"][int(ga["ch_oct"])] / G["may"][int(gb["ch_may"])]}
    rsw = {a: RATE["oct"]["n05"][int(gb["ch_oct"])] / RATE["may"]["n05"][int(ga["ch_may"])],
           b: RATE["oct"]["n05"][int(ga["ch_oct"])] / RATE["may"]["n05"][int(gb["ch_may"])]}
    sw[f"{a}<->{b}"] = {"chi2_base": base["chi2"], "chi2_swapped": s["chi2"], "ndf": s["ndf"], "p_swapped": s["p"], "mean_swapped": s["mean"],
                        "ratios_base": {a: CV["thr3_alltags"][a]["ratio"], b: CV["thr3_alltags"][b]["ratio"]}, "ratios_swapped": {a: cr[a]["ratio"], b: cr[b]["ratio"]},
                        "errs_swapped": {a: cr[a]["err"], b: cr[b]["err"]}, "pulls_swapped": {a: s["pulls"][a], b: s["pulls"][b]},
                        "outliers_swapped": s["outliers_3sigma"], "G_oct_over_may_swapped": gsw, "rate05_oct_over_may_swapped": rsw,
                        "pulls_base": {a: base["pulls"][a], b: base["pulls"][b]}}
    print(f"2 intercambio {a}<->{b}: χ² {base['chi2']:.2f} -> {s['chi2']:.2f} /{s['ndf']}; razones {cr[a]['ratio']:.4f}, {cr[b]['ratio']:.4f}", flush=True)
json.dump({"pairs_tested": pairs, "swaps": sw}, open(OUT + "/r2_swap.json", "w"), indent=1)
for r in pr:
    print("2", r, flush=True)

# ---------------------------------------------------------------- 4. residuo de PDE
mr = json.load(open(OUT + "/muon_refit.json"))["R"]["x"]
Fm, Fo = 1 / (1 - mr["p_may"]), 1 / (1 - mr["p_oct"])
PDE_OV = [0.6, 0.8, 1.0, 1.2, 1.4]
PDE_V = [18.0, 24.0, 29.0, 33.0, 36.0]


def pde(ov):
    if ov < 0.6 or ov > 1.4:
        return float("nan")
    return float(np.interp(ov, PDE_OV, PDE_V))


dv = {int(r["hv_group"]): float(r["dVMon_V"]) for r in csv.DictReader(l for l in open(HVD) if not l.startswith("#"))}
cr = CV["thr3_alltags"]
r4rows, r4 = [], {"F_may": Fm, "F_oct": Fo, "p_may": mr["p_may"], "p_oct": mr["p_oct"], "excluded_oct_G_filled": [], "OV": {}}
for ovm in (1.0, 1.2, 1.4):
    items = []
    for bid, r in cr.items():
        b = BY[bid]
        cm, co = int(b["ch_may"]), int(b["ch_oct"])
        if co in OCT_FILLED:
            if ovm == 1.2:
                r4["excluded_oct_G_filled"].append(bid)
            continue
        rp = r["ratio"] * Fm / Fo
        erp = r["err"] * Fm / Fo
        gq_ = G["oct"][co] / G["may"][cm]
        ovo = ovm * gq_
        pr_ = pde(ovo) / pde(ovm)
        if not np.isfinite(pr_):
            r4rows.append([ovm, bid, b["plane"], b["hv_group"], f"{rp:.4f}", f"{erp:.4f}", f"{gq_:.4f}", f"{ovo:.4f}", "", "", "", "OV_oct outside the datasheet range"])
            continue
        res_ = rp / pr_
        items.append({"id": bid, "plane": b["plane"], "group": int(b["hv_group"]), "res": res_, "err": erp / pr_})
        r4rows.append([ovm, bid, b["plane"], b["hv_group"], f"{rp:.4f}", f"{erp:.4f}", f"{gq_:.4f}", f"{ovo:.4f}", f"{pr_:.4f}", f"{res_:.4f}", f"{erp / pr_:.4f}", ""])
    m, e, chi, ndf = wmean([x["res"] for x in items], [x["err"] for x in items])
    out = {"n": len(items), "mean": m, "err": e, "chi2": chi, "ndf": ndf, "p": pchi(chi, ndf), "groups": {}, "planes": {}}
    for gk in range(1, 9):
        it = [x for x in items if x["group"] == gk]
        if it:
            gm, ge, gc, gn = wmean([x["res"] for x in it], [x["err"] for x in it])
            out["groups"][gk] = {"n": len(it), "mean": gm, "err": ge, "chi2": gc, "ndf": gn, "dV": dv[gk]}
    for pl in ("X", "Y"):
        it = [x for x in items if x["plane"] == pl]
        pm, pe_, pc, pn = wmean([x["res"] for x in it], [x["err"] for x in it])
        out["planes"][pl] = {"n": len(it), "mean": pm, "err": pe_, "chi2": pc, "ndf": pn}
    gs = sorted(out["groups"])
    gx = np.array([dv[g] for g in gs])
    gy = np.array([out["groups"][g]["mean"] for g in gs])
    ge_ = np.array([out["groups"][g]["err"] for g in gs])
    Wg = 1 / ge_ ** 2
    X = np.vstack([np.ones_like(gx), gx]).T
    Cm = np.linalg.inv((X * Wg[:, None]).T @ X)
    bb = Cm @ (X * Wg[:, None]).T @ gy
    g0 = float((Wg * gy).sum() / Wg.sum())
    out["group_fit"] = {"intercept": float(bb[0]), "intercept_err": math.sqrt(Cm[0, 0]), "slope_per_V": float(bb[1]), "slope_err": math.sqrt(Cm[1, 1]),
                        "chi2_line": float((Wg * (gy - X @ bb) ** 2).sum()), "ndf_line": len(gs) - 2, "const": g0,
                        "chi2_const": float((Wg * (gy - g0) ** 2).sum()), "ndf_const": len(gs) - 1}
    out["group_fit"]["p_const"] = pchi(out["group_fit"]["chi2_const"], out["group_fit"]["ndf_const"])
    out["X_minus_Y"] = out["planes"]["X"]["mean"] - out["planes"]["Y"]["mean"]
    out["X_minus_Y_err"] = math.hypot(out["planes"]["X"]["err"], out["planes"]["Y"]["err"])
    r4["OV"][str(ovm)] = out
    print(f"4 OV_mayo {ovm}: residuo {m:.4f} ± {e:.4f}, χ² {chi:.2f}/{ndf} (p {pchi(chi, ndf):.3g}); grupos: pendiente {bb[1]:.4f} ± {math.sqrt(Cm[1, 1]):.4f} /V, "
          f"χ² const {out['group_fit']['chi2_const']:.2f}/{len(gs) - 1}; X − Y {out['X_minus_Y']:.4f} ± {out['X_minus_Y_err']:.4f}", flush=True)
with open(OUT + "/r4_residual_channels.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["OV_may_V", "id", "plane", "hv_group", "R_prim_oct_over_may", "err", "G_oct_over_may", "OV_oct_V", "PDE_ratio", "residual", "err", "note"])
    w.writerows(r4rows)
with open(OUT + "/r4_residual_groups.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["OV_may_V", "unit", "key", "n", "mean_residual", "err", "chi2", "ndf", "dV_V"])
    for ovm, out in r4["OV"].items():
        for gk, g in out["groups"].items():
            w.writerow([ovm, "group", gk, g["n"], f"{g['mean']:.4f}", f"{g['err']:.4f}", f"{g['chi2']:.2f}", g["ndf"], f"{g['dV']:.3f}"])
        for pl, g in out["planes"].items():
            w.writerow([ovm, "plane", pl, g["n"], f"{g['mean']:.4f}", f"{g['err']:.4f}", f"{g['chi2']:.2f}", g["ndf"], ""])
# cota superior R_prim/R_G por definición de F (razones de la estación de xtalk_ratio.csv)
RG = float(np.median([G["oct"][int(b["ch_oct"])] / G["may"][int(b["ch_may"])] for b in bars if DEAD not in (int(b["ch_may"]), int(b["ch_oct"]))
                      and int(b["ch_oct"]) not in OCT_FILLED]))
brow = []
for r in csv.DictReader(open(XT + "/xtalk_ratio.csv")):
    rp, erp = float(r["prim_ratio_oct_over_may"]), float(r["err"])
    brow.append([r["F_definition"], f"{rp:.4f}", f"{erp:.4f}", f"{RG:.4f}", f"{rp / RG:.4f}", f"{erp / RG:.4f}"])
rpm = chk["oct"] / chk["may"] * Fm / Fo
brow.append([f"espectro de los muones, ajuste R (p {mr['p_may']:.4f} / {mr['p_oct']:.4f})", f"{rpm:.4f}", "", f"{RG:.4f}", f"{rpm / RG:.4f}", ""])
with open(OUT + "/r4_bound.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["F_definition", "R_prim_oct_over_may", "err", "R_G_oct_over_may_median", "R_prim_over_R_G_upper_bound", "err"])
    w.writerows(brow)
r4["R_G_median"] = RG
json.dump(r4, open(OUT + "/r4_summary.json", "w"), indent=1)
for r in brow:
    print("4 cota", r, flush=True)
