#!/usr/bin/python3 -I
# stage7.py — PROC-1 etapa 7 (solo lectura): respuesta de hod2019 y el barrido de mayo.
#  1. Mayo v1-v9 (y v10 como referencia): tasa fuera de tiempo por canal (ventana de control de 144 ns de
#     hitdef_extract.C, tasa = −ln(1 − P)/144 ns) sobre 1 y 3 mV en hod2019 y sobre 0.5 y 5 p.e. en hod2018 (G de v8
#     reajustada); fecha de cada run (may_inventory). Sin registro de HV (etapa 1): contra el índice del run. Tiras de hod2019
#     con exceso en tiempo > 3σ a 3 mV por run.
#  2. Barrido de día de octubre (outputs/daq2-15/biasscan_2026_10_23_day/noise_step{0,1}.csv: espectro de pulsos fuera de
#     tiempo por canal y exposición, runs con trigger de software a 0 y −0.5 V): tasa > 1 mV por canal de hod2019, razón
#     −0.5 V / 0 V con error de Poisson; canales compatibles con 1 a 2σ; comparación con las marcas de DAQ2-15.
#  3. Cambio de nivel de hod2019 entre 2026-10-20 13Z y 2026-10-21 10Z: tasa > 1 mV por hora, suma de los 32 canales
#     (outputs/characterization_2026_10_23/work/<run>_hourly.csv de charact_extract.C). Causa: no registrada.
#  4. Exceso en tiempo por canal de hod2019 en mayo v8 a 3 y 35 mV, con y sin eventos coherentes (≥ 8 de los 32 canales sobre
#     3 mV en la misma ventana, la definición de charact_analyse.py para octubre).
# Salida: work/stage7.json, work/fig_stage7.json (proc_plot_stage7.C), tables/hod2019_bias.tex.
# Uso: stage7.py <dir del paquete>
import csv, json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": 7, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
WIN = 144e-9
H19 = list(range(32, 64))
H18 = [c for c in range(32) if c != 3]


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


G8 = g_refit()
inv = {json.loads(l)["label"]: json.loads(l) for l in open(W + "/may_inventory.jsonl")}
rate = lambda k, n: -math.log(1 - k / n) / WIN if k < n else float("inf")

# ---------------------------------------------------------------- 1. mayo por run
runs = {}
for i in range(1, 11):
    run = f"v{i}"
    e = xy.load_hits(f"{W}/hd64_may_{run}_hits.bin", 64)
    n = len(e)
    d = {"index": i, "events": n, "utc_first": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime(inv[run]["utc_first_s"])),
         "utc_last": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime(inv[run]["utc_last_s"])), "per_channel": {}}
    for c in H19:
        k1, k3 = int((e["aC"][:, c] > 1).sum()), int((e["aC"][:, c] > 3).sum())
        s3 = int((e["aS"][:, c] > 3).sum())
        d["per_channel"][str(c)] = {"oot_rate_1mV_Hz": rate(k1, n), "oot_rate_3mV_Hz": rate(k3, n), "n_oot_1mV": k1, "n_oot_3mV": k3,
                                    "in_3mV": s3, "excess_3mV": s3 - k3, "z_3mV": (s3 - k3) / math.sqrt(s3 + k3) if s3 + k3 else 0.0}
    for c in H18:
        k05, k5 = int((e["aC"][:, c] > 0.5 * G8[c]).sum()), int((e["aC"][:, c] > 5 * G8[c]).sum())
        d["per_channel"][str(c)] = {"oot_rate_0p5pe_Hz": rate(k05, n), "oot_rate_5pe_Hz": rate(k5, n), "n_oot_0p5pe": k05, "n_oot_5pe": k5}
    pc = d["per_channel"]
    d["h19_strips_z3_3mV"] = [c for c in H19 if pc[str(c)]["z_3mV"] > 3]
    d["summary"] = {"h19_median_1mV_Hz": float(np.median([pc[str(c)]["oot_rate_1mV_Hz"] for c in H19])),
                    "h19_median_3mV_Hz": float(np.median([pc[str(c)]["oot_rate_3mV_Hz"] for c in H19])),
                    "h18_median_0p5pe_Hz": float(np.median([pc[str(c)]["oot_rate_0p5pe_Hz"] for c in H18])),
                    "h18_median_5pe_Hz": float(np.median([pc[str(c)]["oot_rate_5pe_Hz"] for c in H18])),
                    "h19_n_strips_z3_3mV": len(d["h19_strips_z3_3mV"])}
    gfit = []
    for l in open(f"{W}/hd64_may_{run}_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed") and int(f[0]) != 3:
            gfit.append(float(f[2]))
    d["summary"]["h18_G_median_mV"] = float(np.median(gfit)) if gfit else None
    d["summary"]["h18_G_channels"] = len(gfit)
    runs[run] = d
    print(run, d["utc_first"], n, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in d["summary"].items()})
sets = [set(runs[f"v{i}"]["h19_strips_z3_3mV"]) for i in range(1, 10)]
S["may_runs"] = runs
S["may_hv_record"] = "no hay registro de HV para v1-v10 (etapa 1); el sesgo es desconocido y se grafica contra el índice del run"
S["may_strip_set_changes"] = {"union_v1_v9": sorted(set().union(*sets)), "intersection_v1_v9": sorted(set.intersection(*sets)),
                              "per_run": {f"v{i}": sorted(sets[i - 1]) for i in range(1, 10)}}

# ---------------------------------------------------------------- 2. barrido de día de octubre
BS = REPO + "/outputs/daq2-15/biasscan_2026_10_23_day"


def noise(step):
    exp, cnt = {}, {}
    for l in open(f"{BS}/noise_step{step}.csv"):
        if l.startswith("#exposure"):
            _, c, x = l.split()
            exp[int(c)] = float(x)
            continue
        if l[0] in "#c":
            continue
        c, b, k = l.strip().split(",")
        if float(b) >= 1.0 - 1e-9:
            cnt[int(c)] = cnt.get(int(c), 0) + int(k)
    return exp, cnt


e0, n0 = noise(0)
e1, n1 = noise(1)
FLAGS = {45: "no_hv_dependence", 50: "no_hv_dependence", 52: "no_hv_dependence", 41: "ok (grupo HV sin medir)", 42: "ok (grupo HV sin medir)",
         51: "ok (grupo HV sin medir)"}
r4 = json.load(open(REPO + "/daq/config/mapping/hodo_mapping_2026_10_r4.json"))
grp = {c["dig_ch"]: c.get("hv_group") for c in r4["channels"]}
scan = {}
for c in H19:
    a, b = n0.get(c, 0), n1.get(c, 0)
    r0, r1 = a / e0[c], b / e1[c]
    d = {"count_0V": a, "count_m0p5V": b, "exposure_s": e0[c], "rate_0V_Hz": r0, "rate_m0p5V_Hz": r1, "hv_group_r4": grp.get(c), "flag_daq2_15": FLAGS.get(c, "")}
    if a > 0 and b > 0:
        rr = r1 / r0
        err = rr * math.sqrt(1 / a + 1 / b)
        d.update(ratio=rr, ratio_err=err, consistent_with_1_2sigma=abs(rr - 1) < 2 * err)
    elif a == 0 and b == 0:
        d.update(ratio=None, ratio_err=None, consistent_with_1_2sigma=None, note="0 cuentas en los dos escalones")
    else:
        d.update(ratio=(r1 / r0) if a else float("inf"), ratio_err=None, consistent_with_1_2sigma=False if max(a, b) >= 4 else None,
                 note="un escalón con 0 cuentas")
    scan[str(c)] = d
S["oct_day_scan"] = {"source": BS + "/noise_step{0,1}.csv", "steps": {"0V": "2026-10-23T13:36:35Z", "-0.5V": "2026-10-23T13:58:27Z"}, "per_channel": scan,
                     "median_0V_Hz": float(np.median([scan[str(c)]["rate_0V_Hz"] for c in H19])),
                     "median_m0p5V_Hz": float(np.median([scan[str(c)]["rate_m0p5V_Hz"] for c in H19])),
                     "unchanged_2sigma": [c for c in H19 if scan[str(c)]["consistent_with_1_2sigma"]],
                     "no_counts_both": [c for c in H19 if scan[str(c)].get("note") == "0 cuentas en los dos escalones"]}
print("barrido: sin cambio a 2σ", S["oct_day_scan"]["unchanged_2sigma"], "sin cuentas", S["oct_day_scan"]["no_counts_both"])
for c in sorted(FLAGS):
    d = scan[str(c)]
    print("  ch%d %s: %d / %d cuentas, razón %s" % (c, FLAGS[c], d["count_0V"], d["count_m0p5V"], d.get("ratio")))

# ---------------------------------------------------------------- 3. cambio de nivel (por hora)
CH = REPO + "/outputs/characterization_2026_10_23/work"
hourly = []
for k, lab in (("night", "night_2026_10_19 run 1"), ("run004", "hodoscope_2026_10 run 4"), ("run007", "hodoscope_2026_10 run 7"),
               ("run008", "hodoscope_2026_10 run 8")):
    by = {}
    for r in csv.DictReader(open(f"{CH}/{k}_hourly.csv")):
        c = int(r["ch"])
        if c < 32:
            continue
        assert abs(float(r["thr0_mV"]) - 1.0) < 1e-6
        h = int(r["hour_start_utc_s"])
        x = by.setdefault(h, [0.0, 0, 0.0])
        x[0] += int(r["cnt0"]) / float(r["exposure_s"])
        x[1] += int(r["cnt0"])
        x[2] = float(r["exposure_s"])
    for h in sorted(by):
        s_, n_, ex = by[h]
        hourly.append({"run": lab, "hour_utc": time.strftime("%Y-%m-%dT%HZ", time.gmtime(h)), "t_s": h, "sum_rate_1mV_Hz": s_, "counts": n_,
                       "err_Hz": s_ / math.sqrt(n_) if n_ else None})
S["oct_hourly"] = {"source": CH + "/<run>_hourly.csv (charact_extract.C; thr0 = 1 mV)", "rows": hourly,
                   "cause": "no registrada (René responde qué se hizo en hod2019 esos días)",
                   "before_max_Hz": max(r["sum_rate_1mV_Hz"] for r in hourly if r["t_s"] < 1792504800),
                   "before_min_Hz": min(r["sum_rate_1mV_Hz"] for r in hourly if r["t_s"] < 1792504800),
                   "after_run4_range_Hz": [min(r["sum_rate_1mV_Hz"] for r in hourly if r["run"].endswith("run 4")),
                                           max(r["sum_rate_1mV_Hz"] for r in hourly if r["run"].endswith("run 4"))]}
S["oct_hourly"]["ratio_night_median_over_run4_median"] = float(np.median([r["sum_rate_1mV_Hz"] for r in hourly if r["run"].startswith("night")]) /
                                                               np.median([r["sum_rate_1mV_Hz"] for r in hourly if r["run"].endswith("run 4")]))
print("por hora: noche %.0f-%.0f Hz, run 4 %s, razón de medianas %.1f" % (S["oct_hourly"]["before_min_Hz"], S["oct_hourly"]["before_max_Hz"],
      S["oct_hourly"]["after_run4_range_Hz"], S["oct_hourly"]["ratio_night_median_over_run4_median"]))

# ---------------------------------------------------------------- 4. mayo v8: exceso en tiempo con y sin coherentes
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
coh_in = (e8["aS"][:, 32:64] > 3).sum(1) >= 8
coh_ctl = (e8["aC"][:, 32:64] > 3).sum(1) >= 8
coh = coh_in | coh_ctl
exc = {}
for c in H19:
    d = {}
    for thr in (3, 35):
        for lab, sel in (("all", np.ones(len(e8), bool)), ("no_coherent", ~coh)):
            a, b = int((e8["aS"][sel, c] > thr).sum()), int((e8["aC"][sel, c] > thr).sum())
            d[f"{thr}mV_{lab}"] = {"in": a, "ctl": b, "excess": a - b, "z": (a - b) / math.sqrt(a + b) if a + b else 0.0}
    exc[str(c)] = d
S["may8_h19_excess"] = {"coherent_in": int(coh_in.sum()), "coherent_ctl": int(coh_ctl.sum()), "events": int(len(e8)), "per_channel": exc,
                        "n_z3": {f"{thr}mV_{lab}": sum(1 for c in H19 if exc[str(c)][f"{thr}mV_{lab}"]["z"] > 3) for thr in (3, 35) for lab in ("all", "no_coherent")},
                        "sum_excess": {f"{thr}mV_{lab}": sum(exc[str(c)][f"{thr}mV_{lab}"]["excess"] for c in H19) for thr in (3, 35) for lab in ("all", "no_coherent")}}
print("mayo v8 coherentes: en tiempo", int(coh_in.sum()), "control", int(coh_ctl.sum()), S["may8_h19_excess"]["n_z3"], S["may8_h19_excess"]["sum_excess"])

# ---------------------------------------------------------------- tabla y figura
L = ["% PROC-1 etapa 7: barrido de día de octubre (biasscan_2026_10_23_day, runs con trigger de software): tasa fuera de tiempo > 1 mV de",
     "% hod2019 por canal a 0 V y −0.5 V (exposición 0.1198 s por canal y escalón: 1 cuenta = 8.3 Hz), razón con error de Poisson, y la",
     f"% marca de DAQ2-15. Grupo HV: mapeo r4. Commit {git}.", "\\begin{tabular}{rrrrrl}", "\\hline",
     "Canal & Grupo HV & 0 V (Hz) & $-0.5$ V (Hz) & Razón & Marca DAQ2-15 \\\\", "\\hline"]
for c in H19:
    d = scan[str(c)]
    rt = "--" if d.get("ratio") is None else ("$\\infty$" if d["ratio"] == float("inf") else (f"{d['ratio']:.2f}" + (f" $\\pm$ {d['ratio_err']:.2f}" if d.get("ratio_err") else "")))
    L.append(f"{c} & {d['hv_group_r4'] if d['hv_group_r4'] is not None else '--'} & {d['rate_0V_Hz']:.1f} & {d['rate_m0p5V_Hz']:.1f} & {rt} & "
             f"{d['flag_daq2_15'].replace('_', '\\_')} \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/hod2019_bias.tex", "w").write("\n".join(L) + "\n")
fig = {"may": {"index": list(range(1, 11)), "labels": [f"v{i}" for i in range(1, 11)],
               "h19_1mV": [runs[f"v{i}"]["summary"]["h19_median_1mV_Hz"] for i in range(1, 11)],
               "h19_3mV": [runs[f"v{i}"]["summary"]["h19_median_3mV_Hz"] for i in range(1, 11)],
               "h18_0p5pe": [runs[f"v{i}"]["summary"]["h18_median_0p5pe_Hz"] for i in range(1, 11)],
               "h18_5pe": [runs[f"v{i}"]["summary"]["h18_median_5pe_Hz"] for i in range(1, 11)]},
       "oct": {"t_s": [r["t_s"] for r in hourly], "rate": [r["sum_rate_1mV_Hz"] for r in hourly], "err": [r["err_Hz"] or 0.0 for r in hourly],
               "run": [r["run"] for r in hourly], "change_lo_s": 1792501200, "change_hi_s": 1792576800}}
json.dump(fig, open(W + "/fig_stage7.json", "w"), indent=1)
json.dump(S, open(W + "/stage7.json", "w"), indent=1, default=float)
print("hecho")
