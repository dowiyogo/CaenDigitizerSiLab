#!/usr/bin/python3 -I
# diurnal.py — seguimiento 2026-11-03 (12), punto B: estabilidad diurna de hod2018 en mayo (run B = v8, desde t = 5.25 h) y en octubre (run 8),
# en bins de 1 h UTC (los del árbol "spec" de la ntupla XPLUS-4b), solo lectura de salidas derivadas.
# Por canal (ch3 fuera) y bin:
#  - G: ajuste de picos de fit_gain_hourly.C (fuente central "spec" = método del reajuste de PROC, p0 fijo al del run); se acepta con estado de MIGRAD
#    0 o 1, error > 0 (y el ajuste del run del canal con estado 0 o 1), G dentro de (0.62, 1.38)·semilla (lejos de los límites) y error relativo < 0.10;
#  - DCR = N(> 0.5 p.e.) / exposición y p_ct = N(> 1.5 p.e.) / N(> 0.5 p.e.) de los pulsos de la región pre-trigger (ev_dump.C, como dark_stage2.C /
#    seguimiento 4), con el umbral en p.e. de la G de la hora: G_umbral = G_PROC · G_h / G_run (G_PROC = G_fit de hd64_*_peaks.csv, la escala del
#    seguimiento 4; G_run = el ajuste del run entero con el mismo método que G_h); sin G aceptada en la hora, G_PROC. Corrección de apilamiento del
#    seguimiento 4: p_pile = 1 − exp(−DCR·τ), τ = 32 ns [S], p_ct,corr = (p_ct − p_pile)/(1 − p_pile). (La corrección por colas partidas del
#    seguimiento 5 es un factor por canal y no se aplica: no cambia la variación relativa entre horas.)
# Bins: ≥ 1000 eventos; mayo además inicio ≥ t0 + 5.25 h (t0 = primer trigger de v8). Canal-hora fuera si el canal está en un intervalo HV_OFF de
# runs/good_intervals_2026.csv que toca la hora. Bin "bueno" (mayo): todos sus eventos en tiempo bueno (may_common.good_time); octubre: todos.
# Por bin: mediana de la estación (canales vivos) de G (mV), G/G_ref, DCR (kHz), DCR/DCR_ref, p_ct, p_ct,corr; G_ref, DCR_ref = mediana del canal
# sobre los bins buenos. Dispersión = 1.4826·MAD entre canales; error de la mediana = 1.2533·dispersión/√n.
# Modulación diurna de G/G_ref (bins buenos, central; todos los bins, variante): mínimos cuadrados ponderados con
#   y = c + d·(t − t̄) + a·cos(ω·h) + b·sin(ω·h), ω = 2π/24 h, h = hora local del centro del bin (mayo UTC−4, octubre UTC−3);
# pico a pico = 2·√(a² + b²)/c, fase = hora local del máximo = atan2(b, a)/ω; errores de la covarianza, escalados por √(χ²/ndf) si > 1. Variante sin
# deriva (d = 0). Sin modelo: perfil plegado por hora local (media de los bins buenos de cada hora), máximo − mínimo y hora del máximo.
# Octubre cubre < 24 h: el ajuste con deriva queda degenerado y se da solo sin deriva, con la advertencia.
# Correlación entre bins: Pearson y Spearman de las medianas de la estación (G/G_ref contra DCR/DCR_ref y contra p_ct), p por permutación (20000);
# los bins vecinos no son independientes (deriva lenta): el p es optimista. Pendiente d ln y / d ln G (mínimos cuadrados sin pesos, error de los
# residuos). También la mediana de la correlación de Pearson por canal.
# Temperatura equivalente [inferida bajo supuesto]: ΔT_pp = (ΔG/G)_pp · OV / (dV_bd/dT) con OV = 1.2 V y dV_bd/dT = 50 mV/°C (supuestos de la tarea):
# supone G ∝ OV y que toda la modulación de G es de V_bd(T) con V_bias constante.
# Uso: diurnal.py <work dir> <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

W, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
PEAKS = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work/hd64_%s_peaks.csv"
TS, TAU, OV, DVDT = 16e-9, 32e-9, 1.2, 0.050
H18 = [c for c in range(32) if c != 3]
OMEGA = 2 * math.pi / 24.0
RNG = np.random.default_rng(1103)
DARK_DT = np.dtype([("hour", "<i8"), ("ch", "<i4"), ("pad", "<i4"), ("events", "<i8"), ("expo", "<f8"), ("counts", "<u4", 401)])
EV_DT = np.dtype([("t", "<f8"), ("aS", "<f4", 64), ("kS", "<u2", 64), ("aC", "<f4", 64)])
CAMP = {"may": {"tag": "may8", "peaks": "may_v8", "utc_off": -4, "t_start_h": 5.25, "label": "May 2026, run B (v8)"},
        "oct": {"tag": "oct8", "peaks": "oct8", "utc_off": -3, "t_start_h": None, "label": "October 2026, run 8"}}


def g_proc(name):
    g = np.full(32, np.nan)
    for l in open(PEAKS % name):
        f = l.strip().split(",")
        if f[0] != "ch" and int(f[0]) < 32 and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


def read_g(path):
    run, hr = {}, {}
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        h, c = int(r["hour_utc_s"]), int(r["ch"])
        v = {k: float(r[k]) for k in ("G_mV", "G_err_mV", "chi2")}
        v.update(st=int(r["fit_status"]), ndf=int(r["ndf"]), events=int(r["events"]))
        (run if h == -1 else hr)[(h, c) if h != -1 else c] = v
    return run, hr


def g_ok(v, seed):
    return (v["st"] in (0, 1) and v["G_err_mV"] > 0 and 0.62 * seed < v["G_mV"] < 1.38 * seed and v["G_err_mV"] / v["G_mV"] < 0.10)


def med(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0
    m = float(np.median(x))
    s = 1.4826 * float(np.median(np.abs(x - m)))
    return m, s, 1.2533 * s / math.sqrt(len(x)), len(x)


def wls(t, h, y, e, drift):
    X = [np.ones_like(t)] + ([t - t.mean()] if drift else []) + [np.cos(OMEGA * h), np.sin(OMEGA * h)]
    X = np.vstack(X).T
    w = 1 / e ** 2
    C = np.linalg.inv((X * w[:, None]).T @ X)
    p = C @ (X * w[:, None]).T @ y
    chi = float((((y - X @ p) ** 2) * w).sum())
    ndf = len(y) - X.shape[1]
    sc = max(1.0, chi / ndf) if ndf > 0 else 1.0
    C = C * sc
    ia, ib = X.shape[1] - 2, X.shape[1] - 1
    a, b, c = p[ia], p[ib], p[0]
    A = math.hypot(a, b)
    J = np.zeros(X.shape[1])
    J[ia], J[ib], J[0] = 2 * a / (A * c), 2 * b / (A * c), -2 * A / c ** 2
    pp, pp_e = 2 * A / c, math.sqrt(float(J @ C @ J))
    Jp = np.zeros(X.shape[1])
    Jp[ia], Jp[ib] = -b / A ** 2, a / A ** 2
    ph = (math.atan2(b, a) / OMEGA) % 24.0
    ph_e = math.sqrt(float(Jp @ C @ Jp)) / OMEGA
    out = {"n_bins": int(len(y)), "drift": drift, "c": float(c), "p2p_frac": pp, "p2p_frac_err": pp_e, "phase_local_h": ph, "phase_err_h": ph_e,
           "chi2": chi, "ndf": int(ndf), "err_scale": math.sqrt(sc)}
    if drift:
        out.update(drift_per_day=float(p[1] * 24 / c), drift_per_day_err=float(math.sqrt(C[1, 1]) * 24 / c))
    return out


def rank(x):
    o = np.argsort(x, kind="stable")
    r = np.empty(len(x))
    r[o] = np.arange(len(x))
    return r


def corr(x, y, nperm=20000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.isfinite(x) & np.isfinite(y)
    x, y = x[k], y[k]
    pr = float(np.corrcoef(x, y)[0, 1])
    lx, ly = np.log(x), np.log(y)
    X = np.vstack([np.ones_like(lx), lx - lx.mean()]).T
    C = np.linalg.inv(X.T @ X)
    b = C @ X.T @ ly
    s2 = float(((ly - X @ b) ** 2).sum()) / (len(lx) - 2)
    sr = float(np.corrcoef(rank(x), rank(y))[0, 1])
    perm = np.array([np.corrcoef(x, RNG.permutation(y))[0, 1] for _ in range(nperm)])
    return {"n": int(len(x)), "pearson": pr, "spearman": sr, "p_perm_two_sided": float((np.abs(perm) >= abs(pr)).mean()),
            "dlny_dlnG": float(b[1]), "dlny_dlnG_err": float(math.sqrt(s2 * C[1, 1]))}


summary = {}
rows_bin, rows_ch = [], []
for camp, cfg in CAMP.items():
    gp = g_proc(cfg["peaks"])
    seeds = np.array([float(v) for v in open(f"{W}/{cfg['tag']}_seeds.txt").read().split(",")])
    gp = np.where(np.isfinite(gp), gp, seeds)  # sin G_fit aceptada en PROC: la semilla (G_seed)
    grun, ghr = read_g(f"{W}/{cfg['tag']}_G_spec.csv")
    _, ghr_dark = read_g(f"{W}/{cfg['tag']}_G_dark.csv")
    _, ghr_free = read_g(f"{W}/{cfg['tag']}_G_spec_p0free.csv") if os.path.exists(f"{W}/{cfg['tag']}_G_spec_p0free.csv") else ({}, {})
    grun_dark, _ = read_g(f"{W}/{cfg['tag']}_G_dark.csv")
    dark = np.fromfile(f"{W}/{cfg['tag']}_dark.bin", dtype=DARK_DT)
    ev = np.fromfile(f"{W}/{cfg['tag']}_ev.bin", dtype=EV_DT)
    t = np.sort(ev["t"])
    good_ev = mc.good_time(t, "v8") if camp == "may" else np.ones(len(t), bool)
    hours = sorted(set(int(h) for h in dark["hour"]))
    hv_off = [(a + 4 * 3600, b + 4 * 3600, [int(c) for c in chs.split()]) for st_, a, b, chs, _ in mc.good_intervals("v8") if st_ == "HV_OFF"] \
        if camp == "may" else []
    t_lo = t[0] + cfg["t_start_h"] * 3600 if cfg["t_start_h"] is not None else -np.inf
    D = {(int(r["hour"]), int(r["ch"])): r for r in dark}
    bins = []
    for h in hours:
        inb = (t >= h) & (t < h + 3600)
        n = int(inb.sum())
        if n < 1000 or h < t_lo:
            continue
        bins.append({"hour": h, "events": n, "good": bool(good_ev[inb].all()), "good_frac": float(good_ev[inb].mean()),
                     "t_mid": float(np.median(t[inb]))})
    # canal-hora
    val = {}
    for b in bins:
        h = b["hour"]
        for c in H18:
            if any(c in chs and a < h + 3600 and bb > h for a, bb, chs in hv_off):
                continue
            v = ghr.get((h, c))
            gok = v is not None and g_ok(v, seeds[c]) and c in grun and grun[c]["st"] in (0, 1) and grun[c]["G_err_mV"] > 0
            G = v["G_mV"] if gok else np.nan
            gthr = gp[c] * G / grun[c]["G_mV"] if gok else gp[c]
            r = D[(h, c)]
            cnt = r["counts"].astype(float)
            n05 = cnt[int(0.5 * gthr / 0.25):].sum()
            n15 = cnt[int(1.5 * gthr / 0.25):].sum()
            expo = r["expo"] * TS
            dcr = n05 / expo
            pct = n15 / n05 if n05 > 0 else np.nan
            pile = 1 - math.exp(-dcr * TAU)
            vd = ghr_dark.get((h, c))
            gd = vd["G_mV"] if vd is not None and g_ok(vd, seeds[c]) and grun_dark.get(c, {}).get("st") in (0, 1) else np.nan
            vf = ghr_free.get((h, c))
            gf = vf["G_mV"] if vf is not None and g_ok(vf, seeds[c]) else np.nan
            val[(h, c)] = {"G": G, "G_err": v["G_err_mV"] if gok else np.nan, "G_dark": gd, "G_p0free": gf, "G_thr": gthr, "g_accepted": gok,
                           "DCR": dcr, "N05": n05, "N15": n15, "expo_s": expo, "pct": pct, "pct_corr": (pct - pile) / (1 - pile)}
    goodh = [b["hour"] for b in bins if b["good"]]
    ref = {}
    for c in H18:
        rr = {}
        for k in ("G", "DCR", "pct", "G_dark", "G_p0free"):
            x = [val[(h, c)][k] for h in goodh if (h, c) in val and np.isfinite(val[(h, c)][k])]
            rr[k] = float(np.median(x)) if x else np.nan
        ref[c] = rr
    for b in bins:
        h = b["hour"]
        rel = lambda k: [val[(h, c)][k] / ref[c][k] for c in H18 if (h, c) in val and ref[c][k] > 0]  # noqa: E731
        q = {}
        for name, x in (("G_mV", [val[(h, c)]["G"] for c in H18 if (h, c) in val]), ("G_rel", rel("G")), ("G_dark_rel", rel("G_dark")),
                        ("G_p0free_rel", rel("G_p0free")), ("DCR_kHz", [val[(h, c)]["DCR"] / 1e3 for c in H18 if (h, c) in val]),
                        ("DCR_rel", rel("DCR")), ("pct", [val[(h, c)]["pct"] for c in H18 if (h, c) in val]), ("pct_rel", rel("pct")),
                        ("pct_corr", [val[(h, c)]["pct_corr"] for c in H18 if (h, c) in val])):
            m, s, e, n = med(x)
            q[name] = (m, s, e, n)
        loc = (b["t_mid"] / 3600 + cfg["utc_off"]) % 24
        b.update(local_h=loc, q=q)
        rows_bin.append([camp, h, b["events"], int(b["good"]), f"{b['good_frac']:.4f}",
                         f"{(b['t_mid'] - t[0]) / 3600:.3f}", f"{loc:.3f}"] +
                        sum(([f"{q[k][0]:.5f}", f"{q[k][1]:.5f}", f"{q[k][2]:.5f}", q[k][3]] for k in q), []))
    for (h, c), v in sorted(val.items()):
        rows_ch.append([camp, h, c, int(v["g_accepted"]), f"{v['G']:.5f}", f"{v['G_err']:.5f}", f"{v['G_thr']:.5f}", f"{v['DCR'] / 1e3:.3f}",
                        int(v["N05"]), int(v["N15"]), f"{v['expo_s']:.6f}", f"{v['pct']:.5f}", f"{v['pct_corr']:.5f}", f"{v['G_dark']:.5f}",
                        f"{v['G_p0free']:.5f}"])
    S = {"label": cfg["label"], "utc_offset_h": cfg["utc_off"], "n_bins": len(bins), "n_good_bins": len(goodh),
         "first_bin_utc_s": bins[0]["hour"], "last_bin_utc_s": bins[-1]["hour"], "t0_first_trigger_utc_s": float(t[0]),
         "coverage_h": (bins[-1]["hour"] + 3600 - bins[0]["hour"]) / 3600}
    S["channel_ref"] = {str(c): ref[c] for c in H18}
    S["g_accepted_fraction"] = float(np.mean([v["g_accepted"] for v in val.values()]))
    for sel, bb in (("good", [b for b in bins if b["good"]]), ("all", bins)):
        out = {}
        for key in ("G_rel", "G_dark_rel", "G_p0free_rel", "DCR_rel", "pct_rel"):
            k = [b for b in bb if np.isfinite(b["q"][key][0]) and b["q"][key][2] > 0]
            if len(k) < 6:
                continue
            tt = np.array([b["t_mid"] for b in k]) / 3600
            hh = np.array([b["local_h"] for b in k])
            y = np.array([b["q"][key][0] for b in k])
            e = np.array([b["q"][key][2] for b in k])
            fits = {"no_drift": wls(tt, hh, y, e, False)}
            if S["coverage_h"] > 36:
                fits["drift"] = wls(tt, hh, y, e, True)
            hb = np.floor(hh).astype(int)
            prof = {int(x): float(y[hb == x].mean()) for x in sorted(set(hb))}
            pv = np.array(list(prof.values()))
            imax = max(prof, key=prof.get)
            fits["folded"] = {"n_hours": len(prof), "max_minus_min_frac": float((pv.max() - pv.min()) / np.median(y)), "hour_of_max_local": imax,
                              "hour_of_min_local": min(prof, key=prof.get), "profile": prof}
            fits["range"] = {"max_minus_min_bins_frac": float((y.max() - y.min()) / np.median(y)),
                             "local_h_of_max_bin": float(hh[np.argmax(y)]), "local_h_of_min_bin": float(hh[np.argmin(y)])}
            out[key] = fits
        gq = [b["q"]["G_rel"][0] for b in bb]
        out["corr_G_DCR"] = corr(gq, [b["q"]["DCR_rel"][0] for b in bb])
        out["corr_G_pct"] = corr(gq, [b["q"]["pct_rel"][0] for b in bb])
        per = []
        hs = [b["hour"] for b in bb]
        for c in H18:
            x = np.array([val[(h, c)]["G"] for h in hs if (h, c) in val])
            y = np.array([val[(h, c)]["DCR"] for h in hs if (h, c) in val])
            k = np.isfinite(x) & np.isfinite(y)
            if k.sum() >= 6:
                per.append(float(np.corrcoef(x[k], y[k])[0, 1]))
        out["per_channel_pearson_G_DCR"] = {"n_channels": len(per), "median": float(np.median(per)), "n_positive": int(sum(r > 0 for r in per))}
        main = "drift" if "drift" in out["G_rel"] else "no_drift"
        f = out["G_rel"][main]
        out["temperature_equivalent"] = {
            "status": "derived under assumption: G proportional to OV, V_bias constant, OV = 1.2 V, dV_bd/dT = 50 mV/degC",
            "fit": main, "dT_pp_degC": f["p2p_frac"] * OV / DVDT, "dT_pp_err_degC": f["p2p_frac_err"] * OV / DVDT,
            "dT_pp_degC_folded": out["G_rel"]["folded"]["max_minus_min_frac"] * OV / DVDT}
        S[sel] = out
    summary[camp] = S

hdr = ["campaign", "hour_utc_s", "events", "good_bin", "good_frac", "t_since_first_trigger_h", "local_h_mid"]
for k in ("G_mV", "G_rel", "G_dark_rel", "G_p0free_rel", "DCR_kHz", "DCR_rel", "pct", "pct_rel", "pct_corr"):
    hdr += [f"{k}_median", f"{k}_spread", f"{k}_median_err", f"{k}_n"]
with open(f"{OUT}/diurnal_bins.csv", "w") as f:
    f.write("# diurnal.py: medianas de la estación hod2018 por bin de 1 h UTC (ch3 fuera, canal-hora en HV_OFF fuera); spread = 1.4826 MAD entre canales,"
            " median_err = 1.2533 spread/sqrt(n); *_rel = valor/mediana del canal en los bins buenos; G del ajuste de picos (spec, p0 fijo); DCR > 0.5 p.e.,"
            " pct = N(>1.5)/N(>0.5) de la región pre-trigger; pct_corr con apilamiento tau = 32 ns [S]; local_h: mayo UTC-4, octubre UTC-3\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(hdr)
    w.writerows(rows_bin)
with open(f"{OUT}/diurnal_channel_hours.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "hour_utc_s", "ch", "g_accepted", "G_mV", "G_err_mV", "G_threshold_mV", "DCR_kHz", "N_gt0p5", "N_gt1p5", "exposure_s",
                "pct", "pct_corr", "G_dark_mV", "G_p0free_mV"])
    w.writerows(rows_ch)
json.dump(summary, open(f"{OUT}/diurnal_summary.json", "w"), indent=1)

for camp, S in summary.items():
    print(f"== {S['label']}: {S['n_bins']} bins ({S['n_good_bins']} buenos), cobertura {S['coverage_h']:.0f} h, G aceptada en "
          f"{S['g_accepted_fraction']:.3f} de los canal-hora")
    for sel in ("good", "all"):
        o = S[sel]
        print(f" [{sel}]")
        for key in ("G_rel", "G_dark_rel", "G_p0free_rel", "DCR_rel", "pct_rel"):
            if key not in o:
                continue
            for fk in ("drift", "no_drift"):
                if fk in o[key]:
                    x = o[key][fk]
                    extra = f" deriva {x['drift_per_day'] * 100:+.2f} ± {x['drift_per_day_err'] * 100:.2f} %/día" if fk == "drift" else ""
                    print(f"  {key:13s} {fk:8s}: pico a pico {x['p2p_frac'] * 100:.2f} ± {x['p2p_frac_err'] * 100:.2f} %, máximo a las "
                          f"{x['phase_local_h']:.1f} ± {x['phase_err_h']:.1f} h locales, χ² {x['chi2']:.1f}/{x['ndf']}{extra}")
            fo = o[key]["folded"]
            print(f"  {key:13s} plegado : máx − mín {fo['max_minus_min_frac'] * 100:.2f} % ({fo['n_hours']} horas), máximo {fo['hour_of_max_local']} h, "
                  f"mínimo {fo['hour_of_min_local']} h; bins: máx − mín {o[key]['range']['max_minus_min_bins_frac'] * 100:.2f} %")
        for k in ("corr_G_DCR", "corr_G_pct"):
            c = o[k]
            print(f"  {k}: n {c['n']} Pearson {c['pearson']:+.3f} Spearman {c['spearman']:+.3f} p(perm) {c['p_perm_two_sided']:.4f}; "
                  f"d ln y / d ln G = {c['dlny_dlnG']:+.2f} ± {c['dlny_dlnG_err']:.2f}")
        pc = o["per_channel_pearson_G_DCR"]
        print(f"  Pearson G-DCR por canal: mediana {pc['median']:+.3f}, positivos {pc['n_positive']}/{pc['n_channels']}")
        te = o["temperature_equivalent"]
        print(f"  ΔT equivalente [I, supuesto]: {te['dT_pp_degC']:.2f} ± {te['dT_pp_err_degC']:.2f} °C ({te['fit']}); plegado {te['dT_pp_degC_folded']:.2f} °C")
