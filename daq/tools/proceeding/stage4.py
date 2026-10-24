#!/usr/bin/python3 -I
# stage4.py — PROC-1 etapa 4 (solo lectura): espectros y tiempos.
#  - Espectros de amplitud en p.e. (0-30 p.e., 0.5 p.e.) de 4 canales de hod2018 elegidos por regla (G mediana, G mínima,
#    más ruidoso a 5 p.e. = mayor P(impacto fuera de tiempo), canal Y de ocupación mediana a 5 p.e.); mayo v8 tiempo bueno,
#    definición única (aS/G en tiempo, aC/G fuera de tiempo), G reajustada (A3).
#  - Tiempos de pulso en todo el registro (work/tm_may_v8_*.csv de timing_extract.C): pico de hod2018 y hod2019, componente
#    plana de ch0 (comparada con "0.3"), tasa fuera de tiempo de ch0 sobre 50 mV, los "53 σ" de hod2019 y los pulsos de 3510 ns.
#  - Trigger: salida de analysis/bin/poisson_analysis (work/poisson/run.log, se copia el número tal cual), tasa sin
#    interrupciones, índice de dispersión con intervalos de 1-60 s dentro de los tramos vivos, estabilidad por hora,
#    intervalos más cortos (tiempo muerto) y el ancho del barrido (60 ms o 60 s).
# Salida: work/stage4.json, work/fig_stage4.json (para proc_plot_stage4.C), tables/trigger_stats.tex.
# Uso: stage4.py <dir del paquete>
import json, math, re, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": 4, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
DEAD = 3
LIVE = [c for c in range(32) if c != DEAD]


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


G0 = g_refit()
G = np.where(np.isfinite(G0), G0, np.nanmedian(G0))
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
eg = e8[mc.good_time(e8["t"], "v8")]

# ---------------------------------------------------------------- 4.1 espectros
fit_ok = [c for c in LIVE if np.isfinite(G0[c])]
g_med = float(np.median(G0[fit_ok]))
ch_med = min(fit_ok, key=lambda c: (abs(G0[c] - g_med), c))
ch_low = min(fit_ok, key=lambda c: G0[c])
p_oot = {c: float((eg["aC"][:, c] > 5 * G[c]).mean()) for c in LIVE}
ch_noisy = max(LIVE, key=lambda c: p_oot[c])
net = {c: float((eg["aS"][:, c] > 5 * G[c]).sum() - (eg["aC"][:, c] > 5 * G[c]).sum()) for c in range(16, 32)}
y_med = float(np.median(list(net.values())))
ch_ymed = min(net, key=lambda c: (abs(net[c] - y_med), c))
chosen = [("G mediana", ch_med), ("G mínima", ch_low), ("más ruidoso a 5 p.e.", ch_noisy), ("Y de ocupación mediana", ch_ymed)]
edges = np.arange(0, 30.0001, 0.5)
spec = []
for rule, c in chosen:
    hi, _ = np.histogram(eg["aS"][:, c] / G[c], edges)
    ho, _ = np.histogram(eg["aC"][:, c] / G[c], edges)
    over = int((eg["aS"][:, c] / G[c] >= 30).sum())
    spec.append({"rule": rule, "ch": int(c), "G_mV": float(G[c]), "G_fit": bool(np.isfinite(G0[c])), "in_time": hi.tolist(), "out_of_time": ho.tolist(),
                 "overflow_in_time": over, "P_oot_5pe": p_oot[c], "net_5pe": net.get(c)})
S["spectra"] = {"edges_pe": edges.tolist(), "events": int(len(eg)), "channels": spec, "G_median_fit": g_med, "Y_net_median": y_med,
                "P_oot_5pe": p_oot}
print("canales:", [(r, c) for r, c in chosen])

# ---------------------------------------------------------------- 4.2 tiempos de pulso
PT, NEV = {}, None
for l in open(W + "/tm_may_v8_pulsetimes.csv"):
    if l.startswith("#"):
        NEV = int(re.search(r"eventos (\d+)", l).group(1))
        continue
    if l.startswith("ch,"):
        continue
    f = l.split(",")
    PT[(int(f[0]), float(f[1]) >= 34.9)] = np.array(f[2:], float)  # clave: (canal, umbral alto)
MS = []
for l in open(W + "/tm_may_v8_manuscript.csv"):
    if l[0] in "#c":
        continue
    f = l.strip().split(",")
    MS.append({"ch": int(f[0]), "thr_mV": float(f[1]), "lo_ns": float(f[2]), "hi_ns": float(f[3]), "h": np.array(f[4:], float)})
SW = slice(250, 275)               # 4000-4400 ns
SIDE = np.r_[150:250, 275:624]     # resto del registro donde se buscan pulsos
side_ns = len(SIDE) * 16.0


def tstats(h):
    tot = h[150:].sum()
    win, side = h[SW].sum(), h[SIDE].sum()
    bkg = side * 400.0 / side_ns
    return {"peak_ns": int(np.argmax(h)) * 16, "total": int(tot), "in_4000_4400": int(win), "outside": int(side), "bkg_expected_in_window": float(bkg),
            "excess": float(win - bkg), "excess_z": float((win - bkg) / math.sqrt(bkg)) if bkg > 0 else None,
            "rate_outside_per_event_per_us": float(side / NEV / (side_ns * 1e-3))}


T = {"events": NEV, "per_channel": {str(c): {"low": tstats(PT[(c, False)]), "high": tstats(PT[(c, True)])} for c in range(64)}}
hi18 = [c for c in LIVE if PT[(c, True)][SW].sum() >= 100]
hi19 = [c for c in range(32, 64) if PT[(c, True)][SW].sum() >= 50]
T["peak_hod2018_high_ns"] = {str(c): T["per_channel"][str(c)]["high"]["peak_ns"] for c in hi18}
T["peak_hod2019_high_ns"] = {str(c): T["per_channel"][str(c)]["high"]["peak_ns"] for c in hi19}
s18 = sum(PT[(c, True)] for c in hi18)
s19 = sum(PT[(c, True)] for c in hi19)


def centroid(h):
    k = int(np.argmax(h))
    w = h[k - 3:k + 4]
    return float(((np.arange(k - 3, k + 4) * w).sum() / w.sum()) * 16)


T["summed_peak_ns"] = {"hod2018 (>50 mV, canales con >= 100 en la ventana)": int(np.argmax(s18)) * 16, "hod2019 (>35 mV, >= 50)": int(np.argmax(s19)) * 16,
                       "hod2018_centroid_ns": centroid(s18), "hod2019_centroid_ns": centroid(s19), "channels_18": hi18, "channels_19": hi19}
T["hod2019_outside_35mV"] = {str(c): int(PT[(c, True)][SIDE].sum()) for c in range(32, 64)}
# ch0: componente plana. Definición del manuscrito (80 ns, normalizada al máximo); la única, agrupada a 80 ns.
m0 = next(m for m in MS if m["ch"] == 0 and m["thr_mV"] == 50)
h0 = PT[(0, True)][150:620].reshape(-1, 5).sum(1)  # 2400-9920 ns en 80 ns
T["ch0_flat"] = {"manuscript_def_50mV": {"flat_over_max": float(np.median(m0["h"][63:87]) / m0["h"].max()), "range_ns": "5040-6960, mediana de bins de 80 ns",
                                         "max_bin_ns": int(np.argmax(m0["h"])) * 80},
                 "single_def_50mV_80ns": {"flat_over_max": float(np.median(h0[33:94]) / h0.max()), "range_ns": "5040-9920"},
                 "manuscript_value": 0.3}
# tasa fuera de tiempo de ch0 sobre 50 mV
T["ch0_oot_50mV"] = {"pulses_per_event_per_us": T["per_channel"]["0"]["high"]["rate_outside_per_event_per_us"],
                     "P_window144_oot_hits_file_good_time": float((eg["aC"][:, 0] > 50).mean()),
                     "P_window144_in_time_hits_file_good_time": float((eg["aS"][:, 0] > 50).mean())}
# los 53 σ: fórmula de outputs/timesync_test_20260531_174521/summary_timesync.txt
fifty = []
for m in MS:
    if m["ch"] < 32:
        continue
    n = m["h"].sum()
    nwin = m["h"][51:54].sum()  # 4080-4320 ns: los bins de 80 ns que cubren [4100, 4300]
    flat = n * 200.0 / 2000.0
    fifty.append({"ch": m["ch"], "N_ph_gt_35": int(n), "N_in_4100_4300_bins": int(nwin), "flat_expectation_200_of_2000": flat,
                  "z_as_published": float((nwin - flat) / math.sqrt(flat)),
                  "single_def_outside_window": int(PT[(m["ch"], True)][SIDE].sum()), "single_def_in_window": int(PT[(m["ch"], True)][SW].sum())})
T["fifty_three_sigma"] = {"channels": fifty, "definition": "(N_ventana − N·200/2000)/sqrt(N·200/2000), N = pulsos con ph > 35 mV y mínimo global en "
                                                              "[3000, 5000] ns; la expectativa 'plana' usa el mismo N, no una medida del fondo"}
# pulsos de ~3510 ns: la definición del manuscrito acumula el mínimo global en el primer bin de su rango de búsqueda
m5 = next(m for m in MS if m["ch"] == 0 and m["thr_mV"] == 5)
T["edge_3510"] = {"manuscript_def_ch0_5mV_first_bin_2960_3040": int(m5["h"][37]), "next_10_bins_median": float(np.median(m5["h"][38:48])),
                  "single_def_ch0_low_samples_212_225": PT[(0, False)][212:226].astype(int).tolist(),
                  "single_def_ch0_low_median_150_250": float(np.median(PT[(0, False)][150:250]))}
S["timing"] = T
print("pico hod2018 %d ns, hod2019 %d ns; centroides %.1f / %.1f" % (T["summed_peak_ns"]["hod2018 (>50 mV, canales con >= 100 en la ventana)"],
      T["summed_peak_ns"]["hod2019 (>35 mV, >= 50)"], T["summed_peak_ns"]["hod2018_centroid_ns"], T["summed_peak_ns"]["hod2019_centroid_ns"]))
print("ch0 plana:", T["ch0_flat"])
print("53σ:", [(f["ch"], f["N_ph_gt_35"], round(f["z_as_published"], 1), f["single_def_outside_window"]) for f in fifty])

# ---------------------------------------------------------------- 4.3 trigger
log = open(W + "/poisson/run.log").read()
tab = {}
for l in log.splitlines():
    m = re.match(r"(.+?)\s{2,}(\S+)\s+(\S+)\s*$", l)
    if m and not l.startswith(("  Gap", "[", "Info", "===", "Output", "v7:", "v8:", "Out:")):
        tab[m.group(1).strip()] = (m.group(2), m.group(3))
gaps = {"v7": [], "v8": []}
cur = None
for l in log.splitlines():
    if l.startswith("[v7]"):
        cur = "v7"
    elif l.startswith("[v8]"):
        cur = "v8"
    m = re.match(r"\s+Gap @event (\d+)\s+dur=([\d.]+)s", l)
    if m and cur:
        gaps[cur].append(float(m.group(2)))
TR = {"poisson_analysis_table": tab, "poisson_analysis_gaps_s": gaps, "poisson_analysis_method": {
    "fit": "histograma de 200 bins en [0, 10·T_media], ajuste [0]·exp(−[1]·x) en [max(1.05·τ, 0.02·T_media), 5·T_media], sin los ΔT > 50·T_media",
    "KS": "ΔT sin interrupciones contra Exp(λ_fit), λ del mismo ajuste", "tau_dead": "ΔT mínimo observado",
    "dispersion_scan_us": [100, 500, 1000, 5000, 10000, 60000], "rate_window_us": 60000,
    "source": "analysis/src/poisson_analysis.cxx:117-150, 306, 327 (binario analysis/bin/poisson_analysis)"}}
inv = {json.loads(l)["label"]: json.loads(l) for l in open(W + "/may_inventory.jsonl")}
for run in ("v7", "v8"):
    e = xy.load_hits(W + f"/hd64_may_{run}_hits.bin", 64)
    t = np.sort(e["t"])
    dt = np.diff(t)
    tm = (t[-1] - t[0]) / (len(t) - 1)
    big = dt > 50 * tm
    d = {"N": int(len(t)), "T_total_s": float(t[-1] - t[0]), "lambda_hat_Hz": float((len(t) - 1) / (t[-1] - t[0])),
         "gaps_gt_50Tmean": int(big.sum()), "gap_sum_50Tmean_s": float(dt[big].sum()),
         "lambda_without_gaps_Hz": float((~big).sum() / dt[~big].sum()), "lambda_without_gaps_err_Hz": float(math.sqrt((~big).sum()) / dt[~big].sum()),
         "gaps_gt_20s": int((dt > 20).sum()), "gap_sum_20s_s": float(dt[dt > 20].sum()),
         "lambda_without_gt20s_Hz": float((dt <= 20).sum() / dt[dt <= 20].sum()), "min_dt_ms": float(dt.min() * 1e3),
         "inventory_min_dt_s": inv[run]["min_dt_s"], "inventory_gaps_gt20s": inv[run]["gaps_gt20s"], "inventory_gap_sum_s": inv[run]["gap_sum_s"]}
    # intervalos cortos: con un tiempo muerto τ, la distribución es plana sobre τ
    lam = d["lambda_without_gaps_Hz"]
    for lo, hi in ((0, 1.122e-3), (1.122e-3, 2e-3), (2e-3, 3e-3), (3e-3, 5e-3)):
        n_obs = int(((dt >= lo) & (dt < hi)).sum())
        n_exp = (len(dt) - big.sum()) * (math.exp(-lam * max(lo, 0)) - math.exp(-lam * hi)) if hi > 1.122e-3 else None
        d.setdefault("short_intervals", []).append({"range_ms": [lo * 1e3, hi * 1e3], "observed": n_obs, "expected_poisson_no_deadtime": n_exp})
    # índice de dispersión con intervalos largos, solo con bins enteros dentro de tramos vivos (cortes en ΔT > 20 s)
    cuts = np.where(dt > 20)[0]
    segs = np.split(t, cuts + 1)
    for w in (1, 10, 30, 60, 300):
        cnt = []
        for sg in segs:
            if len(sg) < 2 or sg[-1] - sg[0] < w:
                continue
            nb = int((sg[-1] - sg[0]) // w)
            c, _ = np.histogram(sg, bins=nb, range=(sg[0], sg[0] + nb * w))
            cnt += c.tolist()
        cnt = np.array(cnt, float)
        d.setdefault("dispersion", []).append({"window_s": w, "bins": int(len(cnt)), "mean": float(cnt.mean()), "D": float(cnt.var(ddof=1) / cnt.mean()),
                                               "D_err": float(math.sqrt(2.0 / (len(cnt) - 1)))})
    # estabilidad: tasa por hora de tiempo vivo
    rates, errs = [], []
    for sg in segs:
        if len(sg) < 2:
            continue
        nb = int((sg[-1] - sg[0]) // 3600)
        if nb < 1:
            continue
        c, _ = np.histogram(sg, bins=nb, range=(sg[0], sg[0] + nb * 3600))
        rates += (c / 3600).tolist()
        errs += (np.sqrt(c) / 3600).tolist()
    rates, errs = np.array(rates), np.array(errs)
    wm = (rates / errs ** 2).sum() / (1 / errs ** 2).sum()
    d["hourly"] = {"blocks": int(len(rates)), "weighted_mean_Hz": float(wm), "chi2": float((((rates - wm) / errs) ** 2).sum()), "ndf": int(len(rates) - 1),
                   "min_Hz": float(rates.min()), "max_Hz": float(rates.max()), "rel_sd": float(rates.std() / rates.mean())}
    TR[run] = d
    print(run, {k: d[k] for k in ("lambda_hat_Hz", "lambda_without_gaps_Hz", "lambda_without_gaps_err_Hz", "gaps_gt_50Tmean", "gap_sum_50Tmean_s", "min_dt_ms")})
    print("  dispersión", [(x["window_s"], round(x["D"], 3), round(x["D_err"], 3)) for x in d["dispersion"]])
    print("  hora", d["hourly"], "cortos", d["short_intervals"])
S["trigger"] = TR

# ---------------------------------------------------------------- tabla y figura
P = tab
L = ["% PROC-1 etapa 4: estadística del trigger. Columnas 'poisson_analysis': copiadas de analysis/bin/poisson_analysis (ajuste binned del",
     "% ΔT sin interrupciones > 50 T_media; KS contra el mismo λ). Resto: stage4.py sobre los tiempos del trigger (TTT, µs).", f"% Commit {git}.",
     "\\begin{tabular}{lrr}", "\\hline", "Cantidad & Run A (v7) & Run B (v8) \\\\", "\\hline"]
row = lambda name, a, b: L.append(f"{name} & {a} & {b} \\\\")
row("Triggers", TR["v7"]["N"], TR["v8"]["N"])
row("Duración (h)", f"{TR['v7']['T_total_s'] / 3600:.2f}", f"{TR['v8']['T_total_s'] / 3600:.2f}")
row("$N/T$ (Hz)", f"{TR['v7']['lambda_hat_Hz']:.3f}", f"{TR['v8']['lambda_hat_Hz']:.3f}")
row("Interrupciones $> 50\\,\\bar{T}$ (número, s)", f"{TR['v7']['gaps_gt_50Tmean']}, {TR['v7']['gap_sum_50Tmean_s']:.0f}",
    f"{TR['v8']['gaps_gt_50Tmean']}, {TR['v8']['gap_sum_50Tmean_s']:.0f}")
row("Interrupciones $> 20$ s (número, s)", f"{TR['v7']['gaps_gt_20s']}, {TR['v7']['gap_sum_20s_s']:.0f}", f"{TR['v8']['gaps_gt_20s']}, {TR['v8']['gap_sum_20s_s']:.0f}")
row("$N/T$ sin interrupciones $> 50\\,\\bar{T}$ (Hz)", f"{TR['v7']['lambda_without_gaps_Hz']:.4f} $\\pm$ {TR['v7']['lambda_without_gaps_err_Hz']:.4f}",
    f"{TR['v8']['lambda_without_gaps_Hz']:.4f} $\\pm$ {TR['v8']['lambda_without_gaps_err_Hz']:.4f}")
fit = P.get("λ_fit (exp) [Hz]", ("–", "–"))
row("$\\lambda$ ajuste exponencial (Hz, poisson\\_analysis)", fit[0].replace("±", " $\\pm$ "), fit[1].replace("±", " $\\pm$ "))
ks = P.get("p-value KS", ("–", "–"))
row("$p$ KS (poisson\\_analysis)", ks[0], ks[1])
row("$\\Delta T$ mínimo (ms)", f"{TR['v7']['min_dt_ms']:.3f}", f"{TR['v8']['min_dt_ms']:.3f}")
for w in (10, 60):
    a = next(x for x in TR["v7"]["dispersion"] if x["window_s"] == w)
    b = next(x for x in TR["v8"]["dispersion"] if x["window_s"] == w)
    row(f"$\\sigma^2/\\mu$ en {w} s", f"{a['D']:.3f} $\\pm$ {a['D_err']:.3f}", f"{b['D']:.3f} $\\pm$ {b['D_err']:.3f}")
row("Tasa por hora: $\\chi^2$/ndf", f"{TR['v7']['hourly']['chi2']:.1f}/{TR['v7']['hourly']['ndf']}", f"{TR['v8']['hourly']['chi2']:.1f}/{TR['v8']['hourly']['ndf']}")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/trigger_stats.tex", "w").write("\n".join(L) + "\n")
json.dump({"edges": edges.tolist(), "channels": spec}, open(W + "/fig_stage4.json", "w"), indent=1)
json.dump(S, open(W + "/stage4.json", "w"), indent=1, default=float)
print("hecho")
