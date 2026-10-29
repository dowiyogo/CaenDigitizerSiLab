#!/usr/bin/python3 -I
# temp_corr.py — DAQ2-22 etapa 2 (y datos de la etapa 4): fracción de accidentales del trigger por correlación con la temperatura.
# Idea: el ruido de los SiPM de las tejas (y las accidentales T1·T2) depende de la temperatura; la tasa de muones no. La fracción de
# pulsos oscuros fuera de tiempo de hod2018 (aX, aY a 5 p.e.) hace de termómetro.
# Bloques: run 10 = un segmento cerrado (~1 h; incremental: solo procesa segmentos nuevos), run 8 (control, placas viejas) = horas UTC.
# Por bloque: R = (n − 1)/(t_último − t_primero); G por canal = la del run entero (gain_blocks.C de DAQ2-21: run 10 segmentos 1-14,
# run 8 completo; canales sin ajuste: c1 × mediana de G/c1 de los ajustados), fija en todos los bloques porque con ~1000 eventos por
# segmento el ajuste por bloque no converge (probado: 31 de 32 canales "insufficient" en el segmento 1). Con G fija, aX sigue al DCR
# y a la ganancia, que dependen los dos de la temperatura. aX, aY a 5 p.e.; rX, rY a 10 p.e.; IMon mediano de CH00-07 (slow control, 1 s).
# Modelo [S]: R = R_μ + R_acc·(aX/⟨aX⟩)^p (el ruido de las tejas escala como el DCR de hod2018 con exponente p: p = 1 si domina
# ruido × muón, p = 2 si domina ruido × ruido). Entonces d ln R / d ln aX = a·p y d ln rX / d ln aX = −a·p, con a = R_acc/R.
# Ajustes ponderados (varianza efectiva) de ln R y ln rX contra ln(aX/⟨aX⟩), y lo mismo con aY.
# Uso: temp_corr.py [--since YYYY-MM-DDTHH:MM] [--out DIR]   (salidas por defecto en /data/daq/derived/daq2_22)
import argparse, calendar, csv, gzip, json, math, os, statistics as st, subprocess, sys, time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--since", default=None, help="UTC; solo bloques que empiezan después (para la regresión; la tabla se guarda entera)")
ap.add_argument("--out", default="/data/daq/derived/daq2_22")
A = ap.parse_args()
OUT = A.out
NTD, BLD, OD = f"{OUT}/ntuples", f"{OUT}/build", f"{OUT}/out"
for d in (NTD, BLD, OD):
    os.makedirs(d, exist_ok=True)
RAW = "/data/daq/raw/hodoscope_2026_10"
B10 = "hodoscope_2026_10_20261028T212101Z_run010"
RUN8_NT = "/data/daq/derived/night_2026-10-25/ntuples/oct_run8"
XCH, YCH = [c for c in range(16) if c != 3], list(range(16, 32))
H18 = XCH + YCH
C1 = mc.g_c1()
C1 = np.where(np.isfinite(C1), C1, np.nanmedian(C1))
ROOTENV = dict(os.environ, PATH=os.path.expanduser("~/root6/bin") + ":" + os.environ["PATH"],
               LD_LIBRARY_PATH=os.path.expanduser("~/root6/lib") + ":" + os.environ.get("LD_LIBRARY_PATH", ""))


def root(call):
    for m in ("ntuple_extract.C", "gain_blocks.C"):
        src = f"{REPO}/analysis/night_2026-10-25/{m}"
        dst = f"{BLD}/{m}"
        if not os.path.exists(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
            subprocess.run(["cp", "-p", src, dst], check=True)
    cmd = ["bash", "-c", f"ulimit -v 4194304; cd {BLD} && exec nice -n 10 ionice -c3 timeout 3600 root -l -b -q -e "
           f"'gSystem->AddIncludePath(\"-I{REPO}/daq/core/include\")' '{call}'"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=ROOTENV)
    if r.returncode != 0:
        sys.exit(r.stdout[-2000:] + r.stderr[-2000:])


def gains(of):
    """G por canal del ajuste de todo el run (CSV de gain_blocks.C)."""
    g = {int(r["ch"]): (float(r["G_mV"]), r["status"] in ("ok", "changed")) for r in csv.DictReader(open(of))}
    ok = [c for c in H18 if g[c][1]]
    m = st.median([g[c][0] / C1[c] for c in ok])
    G = np.array([g[c][0] if g[c][1] else C1[c] * m for c in range(32)])
    return np.concatenate([G, np.full(32, np.nan)]), len(ok), m


def fractions(e, G):
    o = {}
    for T, keys in ((5, ("aX", "aY")), (10, ("rX", "rY", "aX10"))):
        Xs, Xc = xy.plane_hits(e, G, T, XCH)
        Ys, Yc = xy.plane_hits(e, G, T, YCH)
        r = xy.xy_fraction(Xs, Ys, Xc, Yc)
        N = r["N"]
        if T == 5:
            for k in ("X", "Y"):
                o[f"a{k}"] = r[f"a{k}"]
                o[f"a{k}_err"] = math.sqrt(r[f"a{k}"] * (1 - r[f"a{k}"]) / N)
        else:
            for k in ("X", "Y"):
                aa, pp = r[f"a{k}"], r[f"N_{k}"] / N
                o[f"r{k}"] = r[f"r{k}"]
                o[f"r{k}_err"] = math.hypot(math.sqrt(pp * (1 - pp) / N) / (1 - aa), (1 - pp) / (1 - aa) ** 2 * math.sqrt(aa * (1 - aa) / N))
    o["N"] = int(len(e))
    return o


def imon(t0, t1):
    """IMon mediano (µA) de CH00-07 en [t0, t1] del slow control de 1 s (hv_YYYYMMDD.csv[.gz])."""
    vals = {c: [] for c in range(8)}
    days = sorted({time.strftime("%Y%m%d", time.gmtime(x)) for x in np.arange(t0, t1 + 1, 3600.0)} | {time.strftime("%Y%m%d", time.gmtime(t1))})
    for d in days:
        for fn in (f"/data/daq/slowcontrol/{d[:4]}/hv_{d}.csv.gz", f"/data/daq/slowcontrol/{d[:4]}/hv_{d}.csv"):
            if not os.path.exists(fn):
                continue
            op = gzip.open if fn.endswith(".gz") else open
            with op(fn, "rt") as f:
                h = next(f).strip().split(",")
                it, ic, ii = h.index("t_req_utc_ns"), h.index("ch"), h.index("IMon_A")
                for l in f:
                    r = l.split(",")
                    c = int(r[ic])
                    if c < 8:
                        t = int(r[it]) * 1e-9
                        if t0 <= t <= t1 and r[ii]:
                            vals[c].append(float(r[ii]) * 1e6)
            break
    return [st.median(vals[c]) if vals[c] else None for c in range(8)]


G10 = gains("/data/daq/derived/daq2_21/out/run10_gain.csv")
G8 = gains("/data/daq/derived/daq2_21/out/run8_gain.csv")
# ---------------------------------------------------------------- bloques de run 10 (incremental)
cache_f = f"{OD}/temp_corr_blocks_run010.json"
cache = json.load(open(cache_f)) if os.path.exists(cache_f) else {}
ix = json.load(open(f"{RAW}/{B10}.index.json"))
for s in ix["segments"]:
    key = f"seg{s['seg']:03d}"
    if key in cache or not s.get("closed_cleanly"):
        continue
    fn = f"{RAW}/{s['file']}"
    pre = f"{NTD}/run010_{key}"
    if not os.path.exists(pre + "_hits.bin"):
        root(f'ntuple_extract.C+("{fn}", "daqv2", "{pre}")')
    e = xy.load_hits(pre + "_hits.bin", 64)
    G, nok, gm = G10
    t0, t1 = s["first_utc_ns"] * 1e-9, s["last_utc_ns"] * 1e-9
    b = {"seg": s["seg"], "t0_utc": t0, "t1_utc": t1, "n_events": s["n_events"], "R_Hz": (s["n_events"] - 1) / (t1 - t0),
         "R_err_Hz": math.sqrt(s["n_events"] - 1) / (t1 - t0), "G_over_c1_median": gm, "n_ch_gain_fit": nok, **fractions(e, G), "IMon_uA": imon(t0, t1)}
    cache[key] = b
    json.dump(cache, open(cache_f, "w"), indent=1)
    print(f"run 10 {key}: R {b['R_Hz']:.4f} aX {b['aX']:.4f} rX {b['rX']:.4f} G/c1 {gm:.4f}", flush=True)
# ---------------------------------------------------------------- control: run 8 por hora (ntupla entera ya validada)
c8f = f"{OD}/temp_corr_blocks_run008.json"
if not os.path.exists(c8f):
    e8 = xy.load_hits(RUN8_NT + "_hits.bin", 64)
    hrs = np.floor(e8["t"] / 3600).astype(np.int64) * 3600
    out8 = {}
    for h in np.unique(hrs):
        m = hrs == h
        if m.sum() < 1500:
            continue
        G, nok, gm = G8
        t = e8["t"][m]
        out8[str(int(h))] = {"t0_utc": float(t[0]), "t1_utc": float(t[-1]), "n_events": int(m.sum()), "R_Hz": (m.sum() - 1) / (t[-1] - t[0]),
                             "R_err_Hz": math.sqrt(m.sum() - 1) / (t[-1] - t[0]), "G_over_c1_median": gm, "n_ch_gain_fit": nok, **fractions(e8[m], G),
                             "IMon_uA": imon(float(t[0]), float(t[-1]))}
        print(f"run 8 {time.strftime('%m-%dT%H', time.gmtime(h))}: R {out8[str(int(h))]['R_Hz']:.4f} aX {out8[str(int(h))]['aX']:.4f}", flush=True)
    json.dump(out8, open(c8f, "w"), indent=1)
out8 = json.load(open(c8f))


# ---------------------------------------------------------------- regresiones
def wfit(x, y, sy, sx):
    """y = A + B x con varianza efectiva σ² = σy² + (B σx)²; devuelve B, σB, A, χ², ndf."""
    x, y, sy, sx = map(np.asarray, (x, y, sy, sx))
    B = 0.0
    for _ in range(10):
        w = 1 / (sy ** 2 + (B * sx) ** 2)
        X = np.vstack([np.ones_like(x), x]).T
        C = np.linalg.inv((X * w[:, None]).T @ X)
        A_, B = C @ (X * w[:, None]).T @ y
    chi = float((((y - A_ - B * x) ** 2) * w).sum())
    return float(B), float(math.sqrt(C[1, 1])), float(A_), chi, len(x) - 2


def regress(blocks, since=None):
    bl = [b for b in blocks if since is None or b["t0_utc"] >= since]
    res = {"n_blocks": len(bl)}
    if len(bl) < 4:
        return res
    for th in ("aX", "aY"):
        a = np.array([b[th] for b in bl]); ea = np.array([b[th + "_err"] for b in bl])
        x = np.log(a / a.mean()); sx = ea / a
        res[f"{th}_range"] = [float(a.min()), float(a.max())]
        res[f"ln_{th}_sd"] = float(x.std())
        R = np.array([b["R_Hz"] for b in bl]); eR = np.array([b["R_err_Hz"] for b in bl])
        B, eB, _, chi, ndf = wfit(x, np.log(R), eR / R, sx)
        res[f"lnR_vs_ln{th}"] = {"slope": B, "slope_err": eB, "chi2": chi, "ndf": ndf}
        r = np.array([b["rX"] for b in bl]); er = np.array([b["rX_err"] for b in bl])
        B2, eB2, _, chi2, ndf2 = wfit(x, np.log(r), er / r, sx)
        res[f"lnrX_vs_ln{th}"] = {"slope": B2, "slope_err": eB2, "chi2": chi2, "ndf": ndf2}
        # a implicada (p = 1 y p = 2) y cota a 95 % si es compatible con 0
        for nm, s, es, sg in (("R", B, eB, 1), ("rX", B2, eB2, -1)):
            v = sg * s
            res[f"a_from_{nm}_{th}"] = {f"p{p}": {"a": v / p, "a_err": es / p, "upper95": max(v, 0) / p + 1.645 * es / p} for p in (1, 2)}
    # rX contra R: con ruido, rX ∝ 1/R (pendiente −1 en log-log)
    R = np.array([b["R_Hz"] for b in bl]); eR = np.array([b["R_err_Hz"] for b in bl])
    r = np.array([b["rX"] for b in bl]); er = np.array([b["rX_err"] for b in bl])
    B3, eB3, _, chi3, ndf3 = wfit(np.log(R / R.mean()), np.log(r), er / r, eR / R)
    res["lnrX_vs_lnR"] = {"slope": B3, "slope_err": eB3, "chi2": chi3, "ndf": ndf3}
    return res


since = calendar.timegm(time.strptime(A.since, "%Y-%m-%dT%H:%M")) if A.since else None
b10 = sorted(cache.values(), key=lambda b: b["t0_utc"])
b8 = sorted(out8.values(), key=lambda b: b["t0_utc"])
REG = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "since": A.since, "run010": regress(b10, since), "run008_control": regress(b8)}
# horas necesarias para distinguir a = 0.1 de a = 0.4 (p = 1: Δpendiente = 0.3) a 3σ (σ ≤ 0.1) con el error actual de ln R vs ln aX
r10 = REG["run010"]
if "lnR_vs_lnaX" in r10:
    s_now, n_now = r10["lnR_vs_lnaX"]["slope_err"], r10["n_blocks"]
    REG["hours_needed"] = {f"sigma_target_{tgt}": {"same_aX_spread_h": n_now * (s_now / tgt) ** 2, "double_aX_spread_h": n_now * (s_now / tgt) ** 2 / 4}
                           for tgt in (0.15, 0.10)}
json.dump(REG, open(f"{OD}/temp_corr_regression.json", "w"), indent=1)
print(json.dumps(REG, indent=1))
