#!/usr/bin/python3 -I
# xtalk.py — seguimiento 2026-11-15: diafonía de hod2018 con una sola definición para mayo y octubre, y su efecto en la razón de p.e. primarios.
# Solo lectura de salidas derivadas. Entradas: <out>/work/{may_v8,oct_run8}_borel.csv (fit_borel.C), xtalk_prep.csv; seguimiento 10
# (charge.bin: aS por evento; _dark.csv: G del análisis); ampcal (medianas de la estación en p.e. medidos); seguimiento 7 (ratio_by_mppc.csv:
# barra, canal de mayo y de octubre, grupo HV); seguimiento 15.4 (hv_delta.csv: ΔV por grupo).
# 1. Por canal y variante del ajuste oscuro: p, F = 1/(1 − p) (media de Borel), P(N ≥ 2) = 1 − e^{−p}; mediana de la estación y error
#    1.2533·σ/√n entre canales; ajustes en el límite de p (≥ 0.94 o ≤ 0.002) se cuentan y se dejan dentro.
# 2. Comparación con F de los seguimientos 5-6 (1.75 / 1.46) y con N(>1.5)/N(>0.5) (0.326 / 0.259; corregido 0.264 / 0.199), leído como
#    P(N ≥ 2) de Borel: p = −ln(1 − P).
# 3. Razón de p.e. primarios octubre/mayo de las medianas de los muones etiquetados (29-54 mm) = razón medida × F_mayo/F_oct, para cada F.
#    Razón medida = 10.4680/16.1627 (G del análisis; error de las dos medianas, bootstrap de ampcal) [M]; intervalo 0.648-0.681 por juego de G.
# 4. Espectro de los muones etiquetados del intervalo lejos (selección de mip_tagged.py: etiqueta > 10 G en la coordenada ortogonal, sonda > 3 p.e.,
#    mismo − mezclados), en unidades de p.e. de la cascada: N = (aS − off)/Gf con Gf y off del ajuste oscuro E del canal (si el ajuste no tiene estado
#    0: G del análisis × mediana(Gf/G) y la mediana de off). Histograma conjunto de las 31 tiras, bins de 1 en [6, 80).
#    Modelo: n primarios ~ Poisson(μ·x), x ~ Landau(MPV 1, anchura w) en [0.02, 12] (camino y pérdida de energía; forma común a las dos campañas [S]),
#    N | n ~ Borel–Tanner(p), resolución gaussiana σ_N² = s0² + N·s1² (medianas del ajuste E, en unidades de Gf), normalizado a los datos en el rango.
#    χ² con error² = mismo + escala²·mezclados. Hipótesis: p fijo por campaña (cada definición) con μ_mayo, μ_oct, w libres; y p_mayo, p_oct libres
#    (perfil Δχ² = 1 en una rejilla de paso 0.04 refinada a 0.01 alrededor del mínimo). Mínimo con Nelder–Mead.
# 5. Anexo: razón mayo/octubre de p.e. primarios de la mediana de los muones etiquetados por grupo HV (barras emparejadas, grupos del seguimiento 7)
#    y por plano, con el p de cada canal (variantes D y E y p común de los muones); error bootstrap de impactos (200); regresión ponderada contra ΔV.
# Salida: <out>/xtalk_channels.csv, xtalk_station.csv, xtalk_ratio.csv, muon_fit.csv, muon_spectra.csv, group_ratio.csv, group_fit.csv, xtalk.json.
# Uso: xtalk.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

OUT = sys.argv[1]
WK = OUT + "/work"
F10 = REPO + "/outputs/hod2018_mip_tagged/2026-11-01_tagged/work"
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
R7 = REPO + "/outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv"
HVD = REPO + "/outputs/hod2018_hv_may_oct/2026-11-15/hv_delta.csv"
DEAD = 3
LIVE = [c for c in range(32) if c != DEAD]
PRE = {"may": "may_v8", "oct": "oct_run8"}
VARS = "ABCDEFGH"
VAR_TXT = {"A": "Borel", "B": "Borel + apilamiento R·τ", "C": "Borel + apilamiento libre", "D": "Borel + apilamiento R·τ + off1 (definición de la tarea)",
           "E": "Borel + off1 + fondo de centelleo", "F": "E + apilamiento libre", "G": "E + apilamiento R·τ", "H": "picos 1-2 + apilamiento R·τ"}
MEAS = {"may": (16.1627, 0.4458), "oct": (10.4680, 0.2954)}  # ampcal, G del análisis [M]
MEAS_RANGE = (0.648, 0.681)


def med_err(v):
    v = np.asarray(v, float)
    return float(np.median(v)), float(1.2533 * np.std(v) / math.sqrt(len(v)))


# ---------------------------------------------------------------- 1. ajustes oscuros
FITS = {}
for camp, pre in PRE.items():
    for r in csv.DictReader(open(f"{WK}/{pre}_borel.csv")):
        if r["p"]:
            FITS[(camp, r["variant"], int(r["ch"]))] = r
ch_rows, st = [], {}
for camp in PRE:
    for v in VARS:
        ps = []
        for c in LIVE:
            r = FITS[(camp, v, c)]
            p = float(r["p"])
            ps.append(p)
            ch_rows.append([camp, v, c, r["p"], r["p_err"], f"{1 / (1 - p):.5f}", f"{1 - math.exp(-p):.5f}", r["q"], r["b"], r["r"], r["off1_mV"], r["Gf_mV"],
                            r["off_mV"], r["chi2_BC"], r["ndf"], r["status"]])
        pm, pe = med_err(ps)
        F = [1 / (1 - p) for p in ps]
        P2 = [1 - math.exp(-p) for p in ps]
        chi = [float(FITS[(camp, v, c)]["chi2_BC"]) / int(FITS[(camp, v, c)]["ndf"]) for c in LIVE]
        st[(camp, v)] = {"p": (pm, pe), "F": med_err(F), "F_of_median_p": 1 / (1 - pm), "P2": med_err(P2), "chi2ndf": float(np.median(chi)),
                         "n_at_upper": int(sum(p >= 0.94 for p in ps)), "n_at_lower": int(sum(p <= 0.002 for p in ps)),
                         "n_status_ne0": int(sum(FITS[(camp, v, c)]["status"] != "0" for c in LIVE))}
with open(OUT + "/xtalk_channels.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "variant", "ch", "p", "p_err", "F", "P_N_ge_2", "q_pileup", "b_background", "r_background", "off1_mV", "Gf_mV", "off_mV", "chi2_BC",
                "ndf", "status"])
    w.writerows(ch_rows)
with open(OUT + "/xtalk_station.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "variant", "model", "p_median", "p_err", "F_median", "F_err", "P_N_ge_2_median", "P_err", "chi2_ndf_median", "n_p_at_upper_limit",
                "n_p_at_lower_limit", "n_status_ne0"])
    for (camp, v), s in st.items():
        w.writerow([camp, v, VAR_TXT[v], f"{s['p'][0]:.5f}", f"{s['p'][1]:.5f}", f"{s['F'][0]:.5f}", f"{s['F'][1]:.5f}", f"{s['P2'][0]:.5f}", f"{s['P2'][1]:.5f}",
                    f"{s['chi2ndf']:.3f}", s["n_at_upper"], s["n_at_lower"], s["n_status_ne0"]])

# ---------------------------------------------------------------- 2-3. definiciones de F y razón primaria
HYP = {"F_amp_f5f6": ((1.75, 0.0), (1.46, 0.0)),
       "Nratio_raw": ((1 / math.exp(-0.0), 0), (1, 0))}  # se rellena abajo
def F_from_P2(P):
    p = -math.log(1 - P)
    return 1 / (1 - p)
HYP = {"F_amp_f5f6 (1.75/1.46)": ((1.75, 0.0), (1.46, 0.0)),
       "N(>1.5)/N(>0.5) raw (0.326/0.259)": ((F_from_P2(0.326), 0.0), (F_from_P2(0.259), 0.0)),
       "N(>1.5)/N(>0.5) corr (0.264/0.199)": ((F_from_P2(0.264), 0.0), (F_from_P2(0.199), 0.0)),
       "F = 1 (no crosstalk)": ((1.0, 0.0), (1.0, 0.0))}
for v in VARS:
    HYP[f"Borel {v}: {VAR_TXT[v]}"] = (st[("may", v)]["F"], st[("oct", v)]["F"])
rm = MEAS["oct"][0] / MEAS["may"][0]
rm_err = rm * math.hypot(MEAS["oct"][1] / MEAS["oct"][0], MEAS["may"][1] / MEAS["may"][0])
ratio_rows = []
for name, ((Fm, eFm), (Fo, eFo)) in HYP.items():
    fr = Fo / Fm
    efr = fr * math.hypot(eFm / Fm, eFo / Fo)
    prim = rm / fr
    eprim = prim * math.hypot(rm_err / rm, efr / fr)
    ratio_rows.append([name, f"{Fm:.4f}", f"{eFm:.4f}", f"{Fo:.4f}", f"{eFo:.4f}", f"{fr:.4f}", f"{efr:.4f}", f"{prim:.4f}", f"{eprim:.4f}",
                       f"{MEAS_RANGE[0] / fr:.4f}", f"{MEAS_RANGE[1] / fr:.4f}", f"{1 / prim:.4f}", f"{eprim / prim ** 2:.4f}"])
with open(OUT + "/xtalk_ratio.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["F_definition", "F_may", "F_may_err", "F_oct", "F_oct_err", "F_oct_over_F_may", "err", "prim_ratio_oct_over_may", "err",
                "prim_ratio_lo_Gset", "prim_ratio_hi_Gset", "prim_ratio_may_over_oct", "err"])
    w.writerows(ratio_rows)

# ---------------------------------------------------------------- 4. espectro de los muones
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("aC", "<f4", 32), ("QS", "<f4", 32), ("QC", "<f4", 32)])
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
D = 37.5 - 15.0 + 2.0 * np.arange(16)
GU = {}
for camp, pre in PRE.items():
    GU[camp] = np.array([float(r["G_mV"]) for r in csv.DictReader(open(f"{F10}/{pre}_dark.csv"))])
EV = {}
for camp, pre in PRE.items():
    e = np.fromfile(f"{F10}/{pre}_charge.bin", dtype=DT)
    e = e[np.argsort(e["t"], kind="stable")]
    if camp == "may":
        e = e[mc.good_time(e["t"], "v8")]
    EV[camp] = e


def tags(camp):
    e, o, g = EV[camp], ORDER[camp], GU[camp]
    out = {}
    for probe, ref in (("X", "Y"), ("Y", "X")):
        refch = np.asarray(o[f"h18_{ref}"])
        tagS = e["aS"][:, refch] > 10.0 * g[refch]
        tagS[:, refch == DEAD] = False
        for L in ("A", "B"):
            s = SIGN[f"oct2026_run8/{probe}{L}"]["sign"]
            d = D if s < 0 else 75.0 - D
            js = [j for j in range(16) if 29 <= d[j] < 54 + 1e-9]
            T = tagS[:, js].any(1)
            for i, c in enumerate(o[f"h18_{probe}"]):
                if i % 2 == (0 if L == "A" else 1) and c != DEAD:
                    out[c] = (T, np.concatenate([T[1:], [False]]))
    return out


TAGS = {camp: tags(camp) for camp in PRE}


def conv(camp, v="E"):
    """Gf y off por canal del ajuste oscuro v (relleno si el estado no es 0)."""
    gf, of = np.full(32, np.nan), np.full(32, np.nan)
    for c in LIVE:
        r = FITS[(camp, v, c)]
        if r["status"] == "0":
            gf[c], of[c] = float(r["Gf_mV"]), float(r["off_mV"])
    ok = np.isfinite(gf)
    rat = np.median(gf[ok] / GU[camp][ok])
    return np.where(ok, gf, GU[camp] * rat), np.where(ok, of, np.median(of[ok])), [c for c in LIVE if not ok[c]]


EDG = np.arange(-0.5, 300.5, 1.0)
SPEC, CONV = {}, {}
for camp in PRE:
    gf, of, miss = conv(camp)
    CONV[camp] = {"filled": miss, "Gf_over_G_median": float(np.median(gf[LIVE] / GU[camp][LIVE])), "off_median_mV": float(np.median(of[LIVE]))}
    e = EV[camp]
    N = len(e)
    s = N / (N - 1)
    hs, hm = np.zeros(len(EDG) - 1), np.zeros(len(EDG) - 1)
    for c in LIVE:
        T, Tn = TAGS[camp][c]
        a = e["aS"][:, c]
        sel = a > 3.0 * GU[camp][c]
        x = (a - of[c]) / gf[c]
        hs += np.histogram(x[sel & T], bins=EDG)[0]
        hm += np.histogram(x[sel & Tn], bins=EDG)[0]
    SPEC[camp] = (hs - s * hm, np.sqrt(hs + s * s * hm))
    sres = [float(FITS[(camp, "E", c)][k]) / float(FITS[(camp, "E", c)]["Gf_mV"]) for c in LIVE for k in ("s0_mV",)]
    s1 = [float(FITS[(camp, "E", c)]["s1_mV"]) / float(FITS[(camp, "E", c)]["Gf_mV"]) for c in LIVE]
    CONV[camp]["s0_rel"] = float(np.median(sres))
    CONV[camp]["s1_rel"] = float(np.median(s1))
with open(OUT + "/muon_spectra.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["N_center", "may_net", "may_err", "oct_net", "oct_err"])
    for k in range(len(EDG) - 1):
        w.writerow([k, f"{SPEC['may'][0][k]:.2f}", f"{SPEC['may'][1][k]:.3f}", f"{SPEC['oct'][0][k]:.2f}", f"{SPEC['oct'][1][k]:.3f}"])

# Landau (convención de TMath::Landau): f(λ) = (1/π) ∫0^∞ exp(−t ln t − λ t) sin(π t) dt
_t = np.linspace(1e-6, 60, 60001)
_dt = _t[1] - _t[0]


def landau_std(lam):
    lam = np.atleast_1d(lam)
    out = np.empty(len(lam))
    for i0 in range(0, len(lam), 50):
        L = lam[i0:i0 + 50, None]
        out[i0:i0 + 50] = (np.exp(-_t * np.log(_t) - L * _t) * np.sin(np.pi * _t)).sum(1) * _dt / np.pi
    return out


XG = np.linspace(0.02, 12.0, 600)
LAM = np.linspace(-4.0, 120.0, 6201)
LANDAU_TAB = landau_std(LAM)


def landau(x, w):
    return np.interp((x - 1.0) / w - 0.22278, LAM, LANDAU_TAB, left=0.0, right=0.0) / w


NMAXP, NMAXN = 260, 300
_n = np.arange(NMAXP + 1)
_lgn = np.array([math.lgamma(k + 1) for k in range(NMAXN + 1)])


def bt_matrix(p):
    """M[N, n] = P(N | n primarios), N = 0..NMAXN, n = 0..NMAXP."""
    N = np.arange(NMAXN + 1)[:, None]
    n = _n[None, :]
    M = np.zeros((NMAXN + 1, NMAXP + 1))
    if p <= 1e-9:
        for k in range(min(NMAXN, NMAXP) + 1):
            M[k, k] = 1.0
        return M
    k = N - n
    ok = (k >= 0) & (n >= 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log(np.where(ok, n, 1)) - np.log(np.where(N > 0, N, 1)) - p * N + k * np.log(p * np.where(N > 0, N, 1)) - _lgn[np.clip(k, 0, NMAXN)]
    M[ok] = np.exp(lg[ok])
    M[0, 0] = 1.0
    return M


def smear_matrix(s0, s1):
    """S[bin, N]: probabilidad de que N caiga en el bin k (centro k, ancho 1)."""
    N = np.arange(NMAXN + 1, dtype=float)
    sig = np.sqrt(s0 * s0 + N * s1 * s1) + 1e-6
    k = np.arange(len(EDG) - 1, dtype=float)
    erf = np.vectorize(math.erf)
    lo = erf((k[:, None] - 0.5 - N[None, :]) / (math.sqrt(2) * sig[None, :]))
    hi = erf((k[:, None] + 0.5 - N[None, :]) / (math.sqrt(2) * sig[None, :]))
    return 0.5 * (hi - lo)


SM = {camp: smear_matrix(CONV[camp]["s0_rel"], CONV[camp]["s1_rel"]) for camp in PRE}
RNG = (6, 80)
SEL = slice(RNG[0], RNG[1])
_lgp = np.array([math.lgamma(k + 1) for k in range(NMAXP + 1)])
_btcache = {}


def model(camp, mu, w, p):
    key = round(p, 6)
    if key not in _btcache:
        if len(_btcache) > 400:
            _btcache.clear()
        _btcache[key] = bt_matrix(p)
    wx = landau(XG, w)
    wx = wx / wx.sum()
    lam = mu * XG
    pois = np.exp(_n[:, None] * np.log(lam[None, :]) - lam[None, :] - _lgp[:, None])
    pn = pois @ wx
    PN = _btcache[key] @ pn
    return SM[camp] @ PN


def chi2_camp(camp, mu, w, p):
    if mu <= 0.2 or w <= 0.02 or not (0 <= p < 0.9):
        return 1e12
    m = model(camp, mu, w, p)[SEL]
    d, e = SPEC[camp][0][SEL], SPEC[camp][1][SEL]
    e = np.where(e > 0, e, 1.0)
    m = m * d.sum() / max(m.sum(), 1e-30)
    return float((((d - m) / e) ** 2).sum())


def nelder_mead(f, x0, step, it=1500, tol=1e-6):
    n = len(x0)
    pts = [np.array(x0, float)] + [np.array(x0, float) + np.eye(n)[i] * step[i] for i in range(n)]
    val = [f(p) for p in pts]
    for _ in range(it):
        o = np.argsort(val)
        pts, val = [pts[i] for i in o], [val[i] for i in o]
        if abs(val[-1] - val[0]) < tol * (1 + abs(val[0])):
            break
        c = np.mean(pts[:-1], axis=0)
        xr = c + (c - pts[-1])
        fr = f(xr)
        if fr < val[0]:
            xe = c + 2 * (c - pts[-1])
            fe = f(xe)
            pts[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            pts[-1], val[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c)
            fc = f(xc)
            if fc < val[-1]:
                pts[-1], val[-1] = xc, fc
            else:
                pts = [pts[0]] + [pts[0] + 0.5 * (p - pts[0]) for p in pts[1:]]
                val = [val[0]] + [f(p) for p in pts[1:]]
    i = int(np.argmin(val))
    return pts[i], val[i]


def fit_fixed(pm, po, x0=(12.0, 9.0, 0.25), quick=False):
    f = lambda x: chi2_camp("may", x[0], x[2], pm) + chi2_camp("oct", x[1], x[2], po)  # noqa: E731
    best = None
    starts = (x0,) if quick else (x0, (x0[0] * 0.6, x0[1] * 0.6, x0[2] * 1.5))
    for start in starts:
        x, v = nelder_mead(f, start, (2.0, 2.0, 0.05))
        x, v = nelder_mead(f, x, (0.5, 0.5, 0.01))
        if best is None or v < best[1]:
            best = (x, v)
    return best


ndf_base = 2 * (RNG[1] - RNG[0])
mu_rows = []
hyp_p = {}
for name, ((Fm, _), (Fo, _)) in HYP.items():
    if name.startswith("Borel"):
        v = name.split()[1].rstrip(":")
        pm, po = st[("may", v)]["p"][0], st[("oct", v)]["p"][0]
    else:
        pm, po = 1 - 1 / Fm, 1 - 1 / Fo
    if pm >= 0.9 or po >= 0.9:
        mu_rows.append([name, f"{pm:.4f}", f"{po:.4f}", "", "", "", "n.d. (p ≥ 0.9)", ndf_base - 3])
        continue
    x, v = fit_fixed(pm, po)
    hyp_p[name] = (pm, po, x, v)
    mu_rows.append([name, f"{pm:.4f}", f"{po:.4f}", f"{x[0]:.3f}", f"{x[1]:.3f}", f"{x[2]:.4f}", f"{v:.2f}", ndf_base - 3])
    print(f"muones, {name}: p {pm:.4f}/{po:.4f} μ {x[0]:.3f}/{x[1]:.3f} w {x[2]:.4f} χ² {v:.2f}/{ndf_base - 3}", flush=True)
# p libres: perfil en una rejilla de p_mayo × p_oct
grid = np.round(np.arange(0.0, 0.61, 0.04), 3)
prof = {}
xbest = (12.0, 9.0, 0.25)
for pm in grid:
    for po in grid:
        x, v = fit_fixed(pm, po, xbest, quick=True)
        prof[(pm, po)] = v
best = min(prof, key=prof.get)
x, v = fit_fixed(best[0], best[1])
# refinamiento local
pmb, pob = best
fine = np.round(np.arange(-0.04, 0.0401, 0.01), 4)
for dm in fine:
    for do in fine:
        a, b = pmb + dm, pob + do
        if a < 0 or b < 0:
            continue
        xx, vv = fit_fixed(a, b, x, quick=True)
        prof[(round(a, 4), round(b, 4))] = vv
best = min(prof, key=prof.get)
x, v = fit_fixed(best[0], best[1])
cmin = prof[best]


def prof_1d(i):
    keys = sorted({k[i] for k in prof})
    vals = [min(val for k, val in prof.items() if k[i] == kk) for kk in keys]
    inside = [k for k, val in zip(keys, vals) if val <= cmin + 1.0]
    return min(inside), max(inside)


lim_m, lim_o = prof_1d(0), prof_1d(1)
mu_rows.append(["p libres (muones)", f"{best[0]:.4f}", f"{best[1]:.4f}", f"{x[0]:.3f}", f"{x[1]:.3f}", f"{x[2]:.4f}", f"{cmin:.2f}", ndf_base - 5])
print(f"muones, p libres: p {best[0]:.4f} [{lim_m[0]:.3f}, {lim_m[1]:.3f}] / {best[1]:.4f} [{lim_o[0]:.3f}, {lim_o[1]:.3f}] χ² {cmin:.2f}/{ndf_base - 5}", flush=True)
# perfil de la diferencia: p común (p_mayo = p_oct)
diag = {k: val for k, val in prof.items() if abs(k[0] - k[1]) < 1e-9}
bd = min(diag, key=diag.get)
with open(OUT + "/muon_fit.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hypothesis", "p_may", "p_oct", "mu_may", "mu_oct", "landau_w", "chi2", "ndf"])
    w.writerows(mu_rows)
    w.writerow(["p común (muones, p_mayo = p_oct)", f"{bd[0]:.4f}", f"{bd[1]:.4f}", "", "", "", f"{diag[bd]:.2f}", ndf_base - 4])
with open(OUT + "/muon_profile.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["p_may", "p_oct", "chi2"])
    for k in sorted(prof):
        w.writerow([k[0], k[1], f"{prof[k]:.3f}"])

# ---------------------------------------------------------------- 5. anexo: grupos HV y planos
bars = [r for r in csv.DictReader(open(R7))]
dv = {int(r["hv_group"]): float(r["dVMon_V"]) for r in csv.DictReader(l for l in open(HVD) if not l.startswith("#"))}


def prim_values(camp, c, Fc):
    e = EV[camp]
    T, Tn = TAGS[camp][c]
    a = e["aS"][:, c] / GU[camp][c]
    sel = a > 3.0
    return a[sel & T] / Fc, a[sel & Tn] / Fc, len(e) / (len(e) - 1)


def wmedian(same, mix, s):
    edg = np.arange(0, 200.0001, 0.05)
    net = np.clip(np.histogram(same, bins=edg)[0] - s * np.histogram(mix, bins=edg)[0], 0, None)
    cs = np.cumsum(net)
    if cs[-1] <= 0:
        return float("nan")
    k = int(np.argmax(cs >= 0.5 * cs[-1]))
    prev = cs[k - 1] if k else 0.0
    return float(edg[k] + (0.5 * cs[-1] - prev) / max(net[k], 1e-12) * 0.05)


def group_median(camp, chans, Fch, rng=None):
    S, M, s = [], [], 1.0
    for c in chans:
        a, b, s = prim_values(camp, c, Fch[c])
        S.append(a)
        M.append(b)
    S, M = np.concatenate(S), np.concatenate(M)
    med = wmedian(S, M, s)
    bs = []
    if rng is not None:
        for _ in range(200):
            bs.append(wmedian(S[rng.integers(0, len(S), len(S))], M[rng.integers(0, len(M), len(M))] if len(M) else M, s))
    return med, (float(np.std(bs)) if bs else float("nan"))


PSETS = {"D": {camp: np.array([float(FITS[(camp, "D", c)]["p"]) if c != DEAD else 0.0 for c in range(32)]) for camp in PRE},
         "E": {camp: np.array([float(FITS[(camp, "E", c)]["p"]) if c != DEAD else 0.0 for c in range(32)]) for camp in PRE},
         "muon_common_p": {"may": np.full(32, best[0]), "oct": np.full(32, best[1])},
         "none (F = 1)": {camp: np.zeros(32) for camp in PRE}}
grows, gfit = [], []
rng = np.random.default_rng(2026111601)
for pname, pset in PSETS.items():
    Fch = {camp: 1.0 / (1.0 - np.clip(pset[camp], 0, 0.95)) for camp in PRE}
    units = [("group", g, [b for b in bars if int(b["hv_group"]) == g]) for g in range(1, 9)] + \
            [("plane", pl, [b for b in bars if b["plane"] == pl]) for pl in ("X", "Y")]
    for kind, key, bl in units:
        bl = [b for b in bl if int(b["ch_may"]) != DEAD and int(b["ch_oct"]) != DEAD and int(b["ch_may"]) < 32]
        if not bl:
            continue
        mm, em = group_median("may", [int(b["ch_may"]) for b in bl], Fch["may"], rng)
        mo, eo = group_median("oct", [int(b["ch_oct"]) for b in bl], Fch["oct"], rng)
        r = mm / mo
        er = r * math.hypot(em / mm, eo / mo)
        Fm = float(np.median([Fch["may"][int(b["ch_may"])] for b in bl]))
        Fo = float(np.median([Fch["oct"][int(b["ch_oct"])] for b in bl]))
        grows.append([pname, kind, key, len(bl), f"{mm:.4f}", f"{em:.4f}", f"{mo:.4f}", f"{eo:.4f}", f"{r:.4f}", f"{er:.4f}", f"{Fm:.4f}", f"{Fo:.4f}",
                      f"{dv[key]:.3f}" if kind == "group" else ""])
    g = [r for r in grows if r[0] == pname and r[1] == "group"]
    X = np.array([float(r[12]) for r in g])
    Y = np.array([float(r[8]) for r in g])
    E = np.array([float(r[9]) for r in g])
    W = 1 / E ** 2
    A = np.vstack([np.ones_like(X), X]).T
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    beta = cov @ (A.T @ (W * Y))
    res = Y - A @ beta
    chi = float((res ** 2 * W).sum())
    w0 = float((W * Y).sum() / W.sum())
    chi0 = float(((Y - w0) ** 2 * W).sum())
    pl = [r for r in grows if r[0] == pname and r[1] == "plane"]
    dxy = float(pl[0][8]) - float(pl[1][8])
    edxy = math.hypot(float(pl[0][9]), float(pl[1][9]))
    gfit.append([pname, f"{beta[0]:.4f}", f"{math.sqrt(cov[0, 0]):.4f}", f"{beta[1]:.4f}", f"{math.sqrt(cov[1, 1]):.4f}", f"{chi:.2f}", len(g) - 2, f"{w0:.4f}",
                 f"{math.sqrt(1 / W.sum()):.4f}", f"{chi0:.2f}", len(g) - 1, f"{dxy:.4f}", f"{edxy:.4f}"])
    print(f"grupos, p {pname}: ordenada {beta[0]:.4f} ± {math.sqrt(cov[0, 0]):.4f}, pendiente {beta[1]:.4f} ± {math.sqrt(cov[1, 1]):.4f} /V, χ² {chi:.2f}/{len(g) - 2}; "
          f"constante {w0:.4f}, χ² {chi0:.2f}/{len(g) - 1}; X − Y {dxy:.4f} ± {edxy:.4f}", flush=True)
with open(OUT + "/group_ratio.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["p_set", "unit", "key", "n_bars", "may_prim_median", "err", "oct_prim_median", "err", "ratio_may_over_oct", "err", "F_may_median", "F_oct_median",
                "dVMon_V"])
    w.writerows(grows)
with open(OUT + "/group_fit.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["p_set", "intercept", "err", "slope_per_V", "err", "chi2", "ndf", "constant", "err", "chi2_constant", "ndf_constant", "X_minus_Y", "err"])
    w.writerows(gfit)

json.dump({"station": {f"{k[0]}_{k[1]}": v for k, v in st.items()}, "conversion": CONV, "meas_ratio_oct_over_may": [rm, rm_err],
           "muon_free_p": {"p_may": best[0], "p_may_1sigma": lim_m, "p_oct": best[1], "p_oct_1sigma": lim_o, "chi2": cmin, "ndf": ndf_base - 5},
           "muon_common_p": {"p": bd[0], "chi2": diag[bd]}}, open(OUT + "/xtalk.json", "w"), indent=1)
print("== estación")
for (camp, v), s in st.items():
    print(f"  {camp} {v} p {s['p'][0]:.4f} ± {s['p'][1]:.4f}  F {s['F'][0]:.4f} ± {s['F'][1]:.4f}  P(N≥2) {s['P2'][0]:.4f}  χ²/ndf {s['chi2ndf']:.2f}  "
          f"límite sup {s['n_at_upper']} inf {s['n_at_lower']} estado≠0 {s['n_status_ne0']}")
print("== razón primaria octubre/mayo, medida", f"{rm:.4f} ± {rm_err:.4f}")
for r in ratio_rows:
    print("  " + " | ".join(r))
print("== conversión", CONV)
