#!/usr/bin/python3 -I
# noise_tail.py — seguimiento 2026-11-15 (15, punto 2): cola del ruido de hod2018 para el estudio de accidentales. Solo lectura de salidas existentes
# de dark_stage2.C (XPLUS-4b; buscador de pulsos del DAQ, separación 3 muestras, piso 1 mV, pedestal local; exposición por canal descontando
# 3 muestras de tiempo muerto por pulso hallado, así la tasa es directa):
#  - mayo, tiempo bueno de v8 (HV-DIAG-5 excluido, 77 558 disparos), región pre-trigger [20, 258) (3.81 µs por registro), G del reajuste de PROC:
#    outputs/hod_dark/2026-10-31_dcr_ct/v8goodnew_pre_{stair,big}.csv;
#  - octubre run 8 (44 176 disparos), región pre-trigger [20, 49) (464 ns por registro), G de c1:
#    outputs/night_2026-10-25/stage2/run8_pre_{stair.csv,big.csv.gz}.
# Tasa sobre T p.e. medidos (amplitud > T·G): T = 3, 5 de la escalera (stair); T = 7, 10, 12, 15 de la lista de pulsos > 6 p.e. (big), con
# amp_mV > T·G y G = mediana de amp_mV/amp_pe del canal (la G con que se escribió); exposición = la de la escalera (misma pasada). Comprobación:
# la cuenta de big sobre 6 p.e. contra la de la escalera a 6 p.e. (difieren por el redondeo de amp_mV a 3 decimales, columna check_6pe).
# Error: Poisson, sqrt(n)/t; con n = 0, límite superior al 95 % = 3.0/t. ch3 (sin señal) fuera de la estación.
# Modelo de Borel: R0·P(N ≥ T | p), R0 = R(> 0.5 p.e.) del canal (escalera), P(N = n) = e^{−pn}(pn)^{n−1}/n!, p = 0.429 (mayo), 0.315 (octubre).
# Estación (por campaña y T): mediana, percentiles 16 y 84 y dispersión = (p84 − p16)/2 sobre los 31 canales vivos; idem para el cociente
# medido/modelo.
# Salida: <out>/noise_tail.csv (una fila por campaña, canal, T; filas "station" con ch = median), noise_tail_summary.json.
# Uso: noise_tail.py <out dir>
import csv, gzip, hashlib, json, math, os, statistics, sys
from collections import defaultdict

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
CAMP = {
    "may_v8_goodtime": {"stair": "outputs/hod_dark/2026-10-31_dcr_ct/v8goodnew_pre_stair.csv", "big": "outputs/hod_dark/2026-10-31_dcr_ct/v8goodnew_pre_big.csv",
                        "expo": "outputs/hod_dark/2026-10-31_dcr_ct/v8goodnew_pre_expo.csv", "p_borel": 0.429, "region": "pre-trigger samples 20-257"},
    "oct_run8": {"stair": "outputs/night_2026-10-25/stage2/run8_pre_stair.csv", "big": "outputs/night_2026-10-25/stage2/run8_pre_big.csv.gz",
                 "expo": "outputs/night_2026-10-25/stage2/run8_pre_expo.csv", "p_borel": 0.315, "region": "pre-trigger samples 20-48"},
}
TS = [3, 5, 7, 10, 12, 15]
DEAD = 3


def borel_ge(T, p):
    s = sum(math.exp(-p * n + (n - 1) * math.log(p * n) - math.lgamma(n + 1)) for n in range(1, T))
    return 1.0 - s


def sha(p):
    return hashlib.sha256(open(os.path.join(REPO, p), "rb").read()).hexdigest()


rows, summ = [], {"inputs": {}, "thresholds_pe": TS, "campaigns": {}}
for camp, c in CAMP.items():
    for k in ("stair", "big", "expo"):
        summ["inputs"][c[k]] = sha(c[k])
    ev = list(csv.DictReader(open(os.path.join(REPO, c["expo"]))))[0]
    st = {(int(r["ch"]), float(r["threshold"]), r["unit"]): r for r in csv.DictReader(open(os.path.join(REPO, c["stair"])))}
    op = gzip.open if c["big"].endswith(".gz") else open
    amps, ratio = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(op(os.path.join(REPO, c["big"]), "rt")):
        ch = int(r["ch"])
        if ch < 32:
            a, ape = float(r["amp_mV"]), float(r["amp_pe"])
            amps[ch].append(a)
            if ape > 0:
                ratio[ch].append(a / ape)
    p = c["p_borel"]
    per = defaultdict(dict)
    for ch in range(32):
        if ch == DEAD:
            continue
        r05 = st[(ch, 0.5, "pe")]
        t = float(r05["exposure_s"])
        R0, R0e = float(r05["rate_cps"]), float(r05["rate_err_cps"])
        G = statistics.median(ratio[ch]) if ratio[ch] else float("nan")
        a = np.asarray(amps[ch])
        n6_big, n6_st = int((a > 6.0 * G).sum()), int(st[(ch, 6.0, "pe")]["counts"])
        for T in TS:
            if (ch, float(T), "pe") in st:
                n, src = int(st[(ch, float(T), "pe")]["counts"]), "stair"
            else:
                n, src = int((a > T * G).sum()), "big"
            rate, err = n / t, math.sqrt(n) / t
            mod = R0 * borel_ge(T, p)
            mod_e = R0e * borel_ge(T, p)
            ratio_m = rate / mod if mod > 0 else float("nan")
            per[T][ch] = (rate, ratio_m)
            rows.append([camp, c["region"], ev["events"], ch, "X" if ch < 16 else "Y", f"{G:.4f}", T, src, n, f"{t:.6g}", f"{rate:.6g}", f"{err:.4g}",
                         f"{3.0 / t:.4g}" if n == 0 else "", f"{R0:.6g}", f"{R0e:.4g}", p, f"{borel_ge(T, p):.6g}", f"{mod:.6g}", f"{mod_e:.4g}",
                         f"{ratio_m:.4g}", f"{(err / mod if mod > 0 else float('nan')):.3g}", f"{n6_big}/{n6_st}"])
    cs = {}
    for T in TS:
        v = np.array([per[T][ch][0] for ch in sorted(per[T])])
        q = np.array([per[T][ch][1] for ch in sorted(per[T])])
        p16, med, p84 = np.percentile(v, [16, 50, 84])
        q16, qmed, q84 = np.percentile(q, [16, 50, 84])
        mod_med = float(np.median([float(st[(ch, 0.5, "pe")]["rate_cps"]) for ch in per[T]])) * borel_ge(T, p)
        cs[str(T)] = {"n_channels": len(v), "rate_median_cps": med, "rate_p16_cps": p16, "rate_p84_cps": p84, "rate_spread_cps": (p84 - p16) / 2,
                      "rate_min_cps": float(v.min()), "rate_max_cps": float(v.max()), "model_with_median_R0_cps": mod_med,
                      "ratio_median": qmed, "ratio_p16": q16, "ratio_p84": q84}
        rows.append([camp, c["region"], ev["events"], "median", "station", "", T, "", "", "", f"{med:.6g}", f"{(p84 - p16) / 2:.4g}", "", "", "", p,
                     f"{borel_ge(T, p):.6g}", f"{mod_med:.6g}", "", f"{qmed:.4g}", f"{(q84 - q16) / 2:.3g}", ""])
    summ["campaigns"][camp] = {"events": int(ev["events"]), "region": c["region"], "p_borel": p, "per_T": cs}

hdr = ["campaign", "region", "triggers", "ch", "plane", "G_mV", "T_pe", "source", "counts", "live_time_s", "rate_cps", "rate_err_cps", "rate_ul95_cps_if_zero",
       "R0_gt0p5pe_cps", "R0_err_cps", "p_borel", "P_borel_N_ge_T", "model_rate_cps", "model_err_cps", "ratio_meas_over_model", "ratio_err",
       "check_6pe_big_over_stair"]
with open(os.path.join(OUT, "noise_tail.csv"), "w") as f:
    f.write("# noise_tail.py (seguimiento 15, 2026-11-15): hod2018 pulse rate above T measured p.e. (amplitude > T*G) in the pre-trigger region, per channel "
            "and campaign; live time = exposure with 3-sample dead time per found pulse removed; error Poisson; model = R0*P_Borel(N>=T|p). "
            "Station rows (ch = median, plane = station): rate_err_cps = (p84-p16)/2 over the 31 live channels, ratio_err = idem for the ratio; ch3 excluded.\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(hdr)
    w.writerows(rows)
json.dump(summ, open(os.path.join(OUT, "noise_tail_summary.json"), "w"), indent=1)
for camp, s in summ["campaigns"].items():
    print(f"{camp}: {s['events']} disparos, {s['region']}, p = {s['p_borel']}")
    for T, x in s["per_T"].items():
        print(f"  T {T:>2} p.e.: mediana {x['rate_median_cps']:.4g} cps (p16 {x['rate_p16_cps']:.4g}, p84 {x['rate_p84_cps']:.4g}, min {x['rate_min_cps']:.4g}, "
              f"max {x['rate_max_cps']:.4g}); Borel con R0 mediana {x['model_with_median_R0_cps']:.4g} cps; medido/modelo mediana {x['ratio_median']:.3g} "
              f"(p16 {x['ratio_p16']:.3g}, p84 {x['ratio_p84']:.3g})")
