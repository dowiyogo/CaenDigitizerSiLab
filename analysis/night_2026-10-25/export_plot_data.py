#!/usr/bin/python3 -I
# export_plot_data.py — XPLUS-4b etapa 6: CSV pequeños para plot_day.C desde las salidas de las etapas 2, 3, 4b y 5.
# Uso: export_plot_data.py <dir outputs/night_2026-10-25>
import csv, json, os, sys

B = sys.argv[1]
P = os.path.join(B, "plots")
os.makedirs(P, exist_ok=True)
# etapa 2: escalera por canal (hod2018) por condición; intervalos sumados con la predicción
for tag, lab in (("sw0V_all", "SW 0 V"), ("swm0p5V_all", "SW -0.5 V"), ("v8good_pre", "mayo v8"), ("run8_post", "run 8 (post)")):
    with open(f"{P}/stair_{tag}.csv", "w") as f:
        f.write("# ch T_pe R_cps\n")
        for r in csv.DictReader(open(f"{B}/stage2/{tag}_stair.csv")):
            c = int(r["ch"])
            if c < 32 and c != 3 and float(r["rate_cps"]) > 0:
                f.write(f"{c} {r['threshold']} {r['rate_cps']}\n")
S = json.load(open(f"{B}/stage2/stage2_summary.json"))
for tag in ("sw0V_all", "v8good_pre"):
    i = S[tag]["intervals"]
    with open(f"{P}/dt_{tag}.csv", "w") as f:
        f.write("# dt_ns observado predicho\n")
        for k, v in i["n_by_dt"].items():
            f.write(f"{int(k) * 16} {v} {i['pred_by_dt'][k]}\n")
# DCR contra ΔV: mediana de R(0.5 p.e.) en mV fijos de G(0 V) y en p.e. del escalón
with open(f"{P}/dcr_vs_dV.csv", "w") as f:
    f.write("# dV R_fixed_mV_0p5 R_pe_step_0p5\n")
    for dv, tags in ((0.0, ("sw0V_all", "sw0V_all")), (-0.5, ("swm0p5V_fixmV", "swm0p5V_all")), (-1.0, ("swm1p0V_fixmV", None)), (-1.5, ("swm1p5V_fixmV", None))):
        vals = []
        for t in tags:
            if t is None:
                vals.append(0.0)
                continue
            rs = [float(r["rate_cps"]) for r in csv.DictReader(open(f"{B}/stage2/{t}_stair.csv")) if int(r["ch"]) < 32 and int(r["ch"]) != 3 and float(r["threshold"]) == 0.5]
            rs.sort()
            vals.append(rs[len(rs) // 2])
        f.write(f"{dv} {vals[0]} {vals[1]}\n")
# etapa 4b: perfiles
for tag in ("v8_T5", "run8_T5", "sw0V_T5"):
    with open(f"{P}/profile_{tag}.csv", "w") as f:
        f.write("# coord(0=X,1=Y) w_ns frac\n")
        for r in csv.DictReader(open(f"{B}/stage4b/{tag}_profile.csv")):
            if r["overlaps_signal"] == "0":
                f.write(f"{0 if r['coord'] == 'X' else 1} {r['w_ns']} {r['frac']}\n")
# etapa 3: q contra g_b
for run in ("v8", "v10", "run8"):
    d = json.load(open(f"{B}/stage3/{run}_stage3.json"))
    with open(f"{P}/stage3_{run}.csv", "w") as f:
        f.write("# g gerr qA qAerr qB qBerr\n")
        for b in d["blocks"]:
            f.write(f"{b['g_b']} {b['g_b_err']} {b['q_A_T5']} {b['q_A_T5_err']} {b['q_B_T5']} {b['q_B_T5_err']}\n")
    with open(f"{P}/stage3_{run}_fit.csv", "w") as f:
        for m in ("A", "B"):
            ft = d["fits"][f"q_{m}_T5"]
            f.write(f"{m} {ft['A']} {ft['s']} {ft['s_err']}\n")
# etapa 5
d = json.load(open(f"{B}/stage5/stage5.json"))
with open(f"{P}/stage5_eps.csv", "w") as f:
    f.write("# i eps5 err valid t_start G_median\n")
    for i, r in enumerate(d["runs"]):
        f.write(f"{i + 1} {r['T5']['eps_XY']} {r['T5']['eps_XY_err']} {1 if r['strip_order_valid'] else 0} {r['start_utc']} {r['G_median_mV']} {1 if r['G_source'] == 'run' else 0}\n")
c = d["combined"]["T5_validated"]
open(f"{P}/stage5_comb.csv", "w").write(f"{c['mean']} {c['err']}\n")
print("ok", P)
