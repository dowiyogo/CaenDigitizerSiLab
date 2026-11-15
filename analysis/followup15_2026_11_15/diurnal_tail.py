#!/usr/bin/python3 -I
# diurnal_tail.py — seguimiento 2026-11-15 (15, punto 3): prueba diurna de la cola del ruido de hod2018 (¿entra luz de día?). Solo lectura de las
# salidas del seguimiento 12 (outputs/hod_stability/2026-11-03_followup12/work/{may8,oct8}_dark.bin de ev_dump.C: por hora UTC y canal, espectro de
# los pulsos oscuros de la región pre-trigger en bins de 0.25 mV de 0 a 100 mV + desborde, y exposición en muestras descontando el tiempo muerto;
# mayo [20, 258), octubre [20, 49)) y de B_diurnal/diurnal_bins.csv (G/G_ref horaria de la estación, followup 12 B).
# Bins y exclusiones como diurnal.py: bins de 1 h UTC con ≥ 1000 eventos; mayo desde t0 + 5.25 h; canal-hora fuera si el canal está en un HV_OFF que
# toca la hora; bin «bueno» = todos sus eventos en tiempo bueno (mayo), todos (octubre). ch3 fuera. Hora local: mayo UTC−4, octubre UTC−3.
# Tasa por canal y hora sobre T p.e. medidos (T = 0.5, 3, 5, 10, 15): pulsos con amplitud > T·G_PROC (G fija del canal, escala de p.e. medidos del
# proceeding: G_fit de hd64_*_peaks.csv o la semilla), / exposición. Variante «G seguida»: umbral T·G_PROC·G_h/G_run como diurnal.py (quita el efecto
# de la ganancia horaria sobre una cola empinada); se toma de la columna G_rel de la estación (misma corrección para todos los canales).
# Por hora: mediana de la estación de la tasa absoluta y de la tasa relativa (tasa/mediana del canal en los bins buenos), dispersión 1.4826·MAD,
# error de la mediana 1.2533·dispersión/√n; tasa agregada de la estación = Σ cuentas / Σ exposición (error de Poisson); exposición por canal (mediana)
# y total. Día = hora local del centro del bin en [6, 18), noche = el resto. Cociente día/noche de la media de la tasa relativa por bin (error de la
# dispersión entre bins: std/√n de cada grupo) y de la tasa agregada (Σ/Σ, Poisson). Exceso de luz: (día/noche a T)/(día/noche a 0.5 p.e.).
# Correlación por bin: Pearson y Spearman de la tasa relativa con cos(2π(h − 14 h)/24) (máximo a media tarde) y con G/G_ref; p por permutación (20000;
# optimista porque bins vecinos no son independientes); ajuste y = c + a cos ωh + b sin ωh (pico a pico 2√(a²+b²)/c, hora del máximo).
# Salida: <out>/diurnal_tail_bins.csv, diurnal_tail.json. Uso: diurnal_tail.py <followup12 dir> <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

F12, OUT = sys.argv[1], sys.argv[2]
W = F12 + "/work"
os.makedirs(OUT, exist_ok=True)
PEAKS = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work/hd64_%s_peaks.csv"
TSMP = 16e-9
H18 = [c for c in range(32) if c != 3]
TPE = [0.5, 3, 5, 10, 15]
OMEGA = 2 * math.pi / 24
RNG = np.random.default_rng(1115)
DARK_DT = np.dtype([("hour", "<i8"), ("ch", "<i4"), ("pad", "<i4"), ("events", "<i8"), ("expo", "<f8"), ("counts", "<u4", 401)])
EV_DT = np.dtype([("t", "<f8"), ("aS", "<f4", 64), ("kS", "<u2", 64), ("aC", "<f4", 64)])
CAMP = {"may": {"tag": "may8", "peaks": "may_v8", "utc_off": -4, "t_start_h": 5.25}, "oct": {"tag": "oct8", "peaks": "oct8", "utc_off": -3, "t_start_h": None}}


def g_proc(name, seeds):
    g = np.full(32, np.nan)
    for l in open(PEAKS % name):
        f = l.strip().split(",")
        if f[0] != "ch" and int(f[0]) < 32 and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, seeds)


def med(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0
    m = float(np.median(x))
    s = 1.4826 * float(np.median(np.abs(x - m)))
    return m, s, 1.2533 * s / math.sqrt(len(x)), len(x)


def rank(x):
    o = np.argsort(x, kind="stable")
    r = np.empty(len(x))
    r[o] = np.arange(len(x))
    return r


def corr(x, y, nperm=20000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.isfinite(x) & np.isfinite(y)
    x, y = x[k], y[k]
    if len(x) < 4:
        return {"n": int(len(x))}
    pr = float(np.corrcoef(x, y)[0, 1])
    perm = np.array([np.corrcoef(x, RNG.permutation(y))[0, 1] for _ in range(nperm)])
    return {"n": int(len(x)), "pearson": pr, "spearman": float(np.corrcoef(rank(x), rank(y))[0, 1]),
            "p_perm_two_sided": float((np.abs(perm) >= abs(pr)).mean())}


def harm(h, y, e):
    X = np.vstack([np.ones_like(h), np.cos(OMEGA * h), np.sin(OMEGA * h)]).T
    w = 1 / e ** 2
    C = np.linalg.inv((X * w[:, None]).T @ X)
    p = C @ (X * w[:, None]).T @ y
    chi = float((((y - X @ p) ** 2) * w).sum())
    ndf = len(y) - 3
    C = C * max(1.0, chi / ndf)
    c, a, b = p
    A = math.hypot(a, b)
    J = np.array([-2 * A / c ** 2, 2 * a / (A * c), 2 * b / (A * c)])
    return {"n_bins": int(len(y)), "p2p_frac": 2 * A / c, "p2p_frac_err": math.sqrt(float(J @ C @ J)), "max_local_h": (math.atan2(b, a) / OMEGA) % 24,
            "chi2": chi, "ndf": int(ndf)}


gb = {}
for r in csv.DictReader(l for l in open(F12 + "/B_diurnal/diurnal_bins.csv") if not l.startswith("#")):
    gb[(r["campaign"], int(r["hour_utc_s"]))] = float(r["G_rel_median"])

res, rows = {"thresholds_pe": TPE, "campaigns": {}}, []
for camp, cfg in CAMP.items():
    seeds = np.array([float(v) for v in open(f"{W}/{cfg['tag']}_seeds.txt").read().split(",")])
    gp = g_proc(cfg["peaks"], seeds)
    dark = np.fromfile(f"{W}/{cfg['tag']}_dark.bin", dtype=DARK_DT)
    ev = np.fromfile(f"{W}/{cfg['tag']}_ev.bin", dtype=EV_DT)
    t = np.sort(ev["t"])
    good_ev = mc.good_time(t, "v8") if camp == "may" else np.ones(len(t), bool)
    hv_off = [(a + 4 * 3600, b + 4 * 3600, [int(c) for c in chs.split()]) for st_, a, b, chs, _ in mc.good_intervals("v8") if st_ == "HV_OFF"] \
        if camp == "may" else []
    t_lo = t[0] + cfg["t_start_h"] * 3600 if cfg["t_start_h"] is not None else -np.inf
    D = {(int(r["hour"]), int(r["ch"])): r for r in dark}
    bins = []
    for h in sorted(set(int(x) for x in dark["hour"])):
        inb = (t >= h) & (t < h + 3600)
        n = int(inb.sum())
        if n < 1000 or h < t_lo:
            continue
        bins.append({"hour": h, "events": n, "good": bool(good_ev[inb].all()), "t_mid": float(np.median(t[inb])),
                     "local_h": (float(np.median(t[inb])) / 3600 + cfg["utc_off"]) % 24, "G_rel": gb.get((camp, h), np.nan)})
    val = {}  # (variante, h, c) -> (cuentas[nT], exposición s)
    for var in ("fixedG", "trackedG"):
        for b in bins:
            h = b["hour"]
            gscale = b["G_rel"] if var == "trackedG" and np.isfinite(b["G_rel"]) else 1.0
            for c in H18:
                if any(c in chs and a < h + 3600 and bb > h for a, bb, chs in hv_off):
                    continue
                r = D[(h, c)]
                cnt = r["counts"].astype(float)
                n = np.array([cnt[min(400, int(T * gp[c] * gscale / 0.25)):].sum() for T in TPE])
                val[(var, h, c)] = (n, r["expo"] * TSMP)
    goodh = [b["hour"] for b in bins if b["good"]]
    out = {"utc_offset_h": cfg["utc_off"], "n_bins": len(bins), "n_good_bins": len(goodh), "coverage_h": (bins[-1]["hour"] + 3600 - bins[0]["hour"]) / 3600,
           "G_PROC_mV": gp.tolist(), "variants": {}}
    for var in ("fixedG", "trackedG"):
        ref = {c: [np.median([val[(var, h, c)][0][k] / val[(var, h, c)][1] for h in goodh if (var, h, c) in val]) for k in range(len(TPE))] for c in H18}
        per = []
        for b in bins:
            h = b["hour"]
            chs = [c for c in H18 if (var, h, c) in val]
            d = {"hour": h}
            ex = np.array([val[(var, h, c)][1] for c in chs])
            for k, T in enumerate(TPE):
                rate = [val[(var, h, c)][0][k] / val[(var, h, c)][1] for c in chs]
                rel = [val[(var, h, c)][0][k] / val[(var, h, c)][1] / ref[c][k] for c in chs if ref[c][k] > 0]
                N = sum(val[(var, h, c)][0][k] for c in chs)
                d[T] = {"med": med(rate), "rel": med(rel), "pool_counts": N, "pool_expo": float(ex.sum()), "pool_rate": N / ex.sum(),
                        "pool_err": math.sqrt(N) / ex.sum(), "nch": len(chs), "expo_med": float(np.median(ex)),
                        "counts_med_per_ch": float(np.median([val[(var, h, c)][0][k] for c in chs]))}
            per.append(d)
            if var == "fixedG" or True:
                for k, T in enumerate(TPE):
                    x = d[T]
                    rows.append([camp, var, h, b["events"], int(b["good"]), f"{(b['t_mid'] - t[0]) / 3600:.3f}", f"{b['local_h']:.3f}",
                                 int((6 <= b["local_h"] < 18)), f"{b['G_rel']:.5f}", T, x["nch"], f"{x['expo_med']:.6g}", f"{x['pool_expo']:.6g}",
                                 f"{x['counts_med_per_ch']:.1f}", f"{x['med'][0]:.6g}", f"{x['med'][1]:.4g}", f"{x['med'][2]:.4g}", f"{x['rel'][0]:.5f}",
                                 f"{x['rel'][1]:.4f}", f"{x['rel'][2]:.4f}", int(x["pool_counts"]), f"{x['pool_rate']:.6g}", f"{x['pool_err']:.4g}"])
        vres = {}
        for sel in ("good", "all"):
            idx = [i for i, b in enumerate(bins) if b["good"] or sel == "all"]
            day = [i for i in idx if 6 <= bins[i]["local_h"] < 18]
            night = [i for i in idx if not 6 <= bins[i]["local_h"] < 18]
            sres = {"n_day": len(day), "n_night": len(night)}
            dn = {}
            for T in TPE:
                yd = np.array([per[i][T]["rel"][0] for i in day])
                yn = np.array([per[i][T]["rel"][0] for i in night])
                rd, rn = yd.mean() if len(yd) else np.nan, yn.mean() if len(yn) else np.nan
                ed = yd.std(ddof=1) / math.sqrt(len(yd)) if len(yd) > 1 else np.nan
                en = yn.std(ddof=1) / math.sqrt(len(yn)) if len(yn) > 1 else np.nan
                Nd, Ed = sum(per[i][T]["pool_counts"] for i in day), sum(per[i][T]["pool_expo"] for i in day)
                Nn, En = sum(per[i][T]["pool_counts"] for i in night), sum(per[i][T]["pool_expo"] for i in night)
                pr = (Nd / Ed) / (Nn / En) if Nd > 0 and Nn > 0 and Ed > 0 and En > 0 else np.nan
                pre = pr * math.sqrt(1 / Nd + 1 / Nn) if np.isfinite(pr) else np.nan
                h_ = np.array([bins[i]["local_h"] for i in idx])
                y_ = np.array([per[i][T]["rel"][0] for i in idx])
                e_ = np.array([per[i][T]["rel"][2] for i in idx])
                e_ = np.where(e_ > 0, e_, np.nanmedian(e_[e_ > 0]) if (e_ > 0).any() else 1.0)
                dn[str(T)] = {"rel_day_mean": rd, "rel_day_err": ed, "rel_night_mean": rn, "rel_night_err": en,
                              "ratio_rel": rd / rn, "ratio_rel_err": (rd / rn) * math.sqrt((ed / rd) ** 2 + (en / rn) ** 2),
                              "pooled_day_counts": int(Nd), "pooled_day_live_s": Ed, "pooled_night_counts": int(Nn), "pooled_night_live_s": En,
                              "ratio_pooled": pr, "ratio_pooled_err": pre,
                              "corr_vs_cos14h": corr(y_, np.cos(OMEGA * (h_ - 14.0))),
                              "corr_vs_G_rel": corr(y_, [bins[i]["G_rel"] for i in idx]),
                              "harmonic": harm(h_, y_, e_) if len(idx) > 4 else None}
            for T in TPE[1:]:
                a, b0 = dn[str(T)], dn["0.5"]
                a["excess_over_0p5_rel"] = a["ratio_rel"] / b0["ratio_rel"]
                a["excess_over_0p5_rel_err"] = a["excess_over_0p5_rel"] * math.hypot(a["ratio_rel_err"] / a["ratio_rel"], b0["ratio_rel_err"] / b0["ratio_rel"])
                a["excess_over_0p5_pooled"] = a["ratio_pooled"] / b0["ratio_pooled"]
            sres["per_T"] = dn
            vres[sel] = sres
        out["variants"][var] = vres
    res["campaigns"][camp] = out

hdr = ["campaign", "variant", "hour_utc_s", "events", "good_bin", "t_since_first_trigger_h", "local_h_mid", "day_06_18", "G_rel_station", "T_pe", "n_channels",
       "live_s_median_per_channel", "live_s_station_total", "counts_median_per_channel", "rate_median_cps", "rate_spread_cps", "rate_median_err_cps",
       "rel_median", "rel_spread", "rel_median_err", "pooled_counts", "pooled_rate_cps", "pooled_rate_err_cps"]
with open(f"{OUT}/diurnal_tail_bins.csv", "w") as f:
    f.write("# diurnal_tail.py (seguimiento 15, 2026-11-15): hod2018 pre-trigger pulse rate above T measured p.e. per 1 h UTC bin; station median over live "
            "channels (ch3 and HV_OFF channel-hours out), spread 1.4826 MAD, median error 1.2533 spread/sqrt(n); rel = rate / channel median over good "
            "bins; pooled = sum counts / sum live time (Poisson); variant fixedG = threshold T*G_PROC, trackedG = T*G_PROC*G_rel(hour); local time May UTC-4, Oct UTC-3\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(hdr)
    w.writerows(rows)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if math.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


json.dump(clean(res), open(f"{OUT}/diurnal_tail.json", "w"), indent=1)
for camp, o in res["campaigns"].items():
    print(f"== {camp}: {o['n_bins']} bins ({o['n_good_bins']} buenos), cobertura {o['coverage_h']:.1f} h")
    for var, v in o["variants"].items():
        for sel, s in v.items():
            print(f" [{var}, {sel}] día {s['n_day']} bins, noche {s['n_night']} bins")
            for T, d in s["per_T"].items():
                ex = f"; exceso sobre 0.5 p.e. {d['excess_over_0p5_rel']:.4f} ± {d['excess_over_0p5_rel_err']:.4f} (agregado {d['excess_over_0p5_pooled']:.4f})" \
                    if "excess_over_0p5_rel" in d else ""
                cc, cg, hm = d["corr_vs_cos14h"], d["corr_vs_G_rel"], d["harmonic"]
                print(f"   T {T:>4}: día/noche rel {d['ratio_rel']:.4f} ± {d['ratio_rel_err']:.4f}; agregado {d['ratio_pooled']:.4f} ± {d['ratio_pooled_err']:.4f} "
                      f"(cuentas día {d['pooled_day_counts']} en {d['pooled_day_live_s']:.4g} s, noche {d['pooled_night_counts']} en {d['pooled_night_live_s']:.4g} s)"
                      f"{ex}; r(cos 14 h) {cc.get('pearson', float('nan')):.3f} p {cc.get('p_perm_two_sided', float('nan')):.4f}; "
                      f"r(G_rel) {cg.get('pearson', float('nan')):.3f} p {cg.get('p_perm_two_sided', float('nan')):.4f}; "
                      + (f"pico a pico {hm['p2p_frac']:.4f} ± {hm['p2p_frac_err']:.4f}, máx {hm['max_local_h']:.1f} h" if hm else ""))
