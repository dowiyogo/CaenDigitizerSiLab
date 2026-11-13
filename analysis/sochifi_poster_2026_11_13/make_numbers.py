#!/usr/bin/python3 -I
# make_numbers.py — SOCHIFI 2026 (póster), 2026-11-13: poster_numbers.json con todos los números del póster, copiados de salidas existentes (sin
# selecciones, ajustes ni estimadores nuevos). Cada entrada: value, error (o errors), unit, sample, source, commit, tag (measured | calculated |
# assumed) y note. "calculated" = cuenta aritmética sobre números existentes (por ejemplo q·N), indicada en note.
# Uso: make_numbers.py <out dir>
import csv, json, math, os, re, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
O = os.path.join(REPO, "outputs")
PK = "outputs/proceeding_hvdiag5/2026-10-31_full"
GOOD = "run B (mayo v8), tiempo bueno, 77 558 triggers"
STAT = "run B, segmento estacionario t ≥ 5.25 h, 142 544 triggers"
nums = {}


def add(key, value, error, unit, sample, source, commit, tag, note=""):
    nums[key] = {"value": value, "error": error, "unit": unit, "sample": sample, "source": source, "commit": commit, "tag": tag, "note": note}


# trigger (log de interarrival_stationary.C, revisión 2026-11-02)
log = open(os.path.join(O, "proceeding_revision/2026-11-02_revision/logs/interarrival_stationary.log")).read()
m = re.search(r"lambda_fit = ([\d.]+) \+- ([\d.]+) Hz", log)
add("trigger_rate", float(m.group(1)), float(m.group(2)), "Hz", STAT, "outputs/proceeding_revision/2026-11-02_revision/logs/interarrival_stationary.log",
    "e0ae42a (interarrival_stationary.C)", "measured", "ajuste exponencial de los intervalos; KS contra el mismo λ")
m = re.search(r"KS vs exp\(lambda_fit\): D = ([\d.]+), p = ([\d.]+)", log)
add("trigger_ks_p", float(m.group(2)), None, "", STAT, "idem", "e0ae42a", "measured", f"D = {m.group(1)}")
m = re.search(r"duration ([\d.]+) h", log)
add("stationary_duration", float(m.group(1)), None, "h", STAT, "idem", "e0ae42a", "measured")
add("stationary_triggers", 142544, None, "triggers", STAT, "idem", "e0ae42a", "measured")
m = re.search(r"dead time (\d+) us", log)
add("dead_time_min_interval", float(m.group(1)), None, "µs", STAT, "idem", "e0ae42a", "measured", "intervalo mínimo entre triggers (tiempo muerto)")
s4 = json.load(open(os.path.join(REPO, PK, "work/stage4.json")))["trigger"]
pt = s4["poisson_analysis_table"]
add("dead_time_fraction", float(pt["Fracción tiempo muerto [%]"][1]), None, "%", "run B completo (v8, 143 814 triggers)", f"{PK}/work/stage4.json (poisson_analysis_table)",
    "e654684 (stage4.py, analysis/bin/poisson_analysis)", "measured")
disp = {d["window_s"]: d for d in s4["v8"]["dispersion"]}
for w in (10, 60):
    add(f"dispersion_index_{w}s", disp[w]["D"], disp[w]["D_err"], "", "run B completo dentro de los tramos vivos (v8)", f"{PK}/work/stage4.json (v8.dispersion); tables/trigger_stats.tex",
        "e654684 (stage4.py)", "measured", f"σ²/μ de las cuentas en ventanas de {w} s")

# eficiencia
fe = json.load(open(os.path.join(O, "proceeding_revision/2026-11-03_band_legend/fig_eff_g09551_bandrunB.json")))
s = fe["series"][0]
for T in (5, 10):
    i = fe["T"].index(T)
    e = s["eps"][i]
    add(f"eps_XY_{T}pe", e, {"stat": s["stat"][i], "estimator": s["syst"][i], "acceptance": e * fe["f18_rel_syst"]}, "", GOOD,
        "outputs/proceeding_revision/2026-11-03_band_legend/fig_eff_g09551_bandrunB.json", "5862882 (g_dead_ch3.py) / 4debf58 (proc_plot_eff.C)", "measured",
        "ε_XY = q/(f_2018·g); stat = σ_q/(f·g); estimator = |q_mix − q_ctl|/(f·g); acceptance = ε × f18_rel_syst (calculated)")
add("acceptance_rel_uncertainty", fe["f18_rel_syst"], None, "", "MC del banco", "idem (f18_rel_syst)", "5862882", "calculated", "f_2018 (3.50 %) ⊕ g (3.06 %)")
cad = json.load(open(os.path.join(REPO, PK, "work/geom_sens_CAD_1e8.json")))["cos3"]["nominal"]
s2 = json.load(open(os.path.join(REPO, PK, "work/stage2.json")))
add("f_2018", cad["f"], {"mc_stat": cad["f_err"], "relative_total": s2["acceptance_rel_syst"]}, "", "MC del banco (CAD, cos³θ, 1e8 rayos)", f"{PK}/work/geom_sens_CAD_1e8.json; stage2.json",
    "a5ba581 (geom_sensitivity.py) / dae37a1 (stage2.py)", "calculated", "fracción de triggers que cruza el cuadrado de 33 mm de hod2018; syst = MC ⊕ placas z ± 2 mm ⊕ 10 mm")
gd = json.load(open(os.path.join(O, "proceeding_revision/2026-11-02_revision/g_dead_ch3/g_dead_ch3.json")))["stage5_sample"]
add("g_dead_ch3", gd["average"]["g_dead"], {"mc_stat": max(gd["X_top"]["g_dead"]["err"], gd["Y_top"]["g_dead"]["err"]), "relative_syst": 0.03 / 0.98}, "",
    "MC del banco, muestra de stage5", "outputs/proceeding_revision/2026-11-02_revision/g_dead_ch3/g_dead_ch3.json", "5862882 (g_dead_ch3.py)", "calculated",
    "P(X ∧ Y | cruza el cuadrado) con la tira X12 (ch3) apagada; promedio de los dos órdenes de bloque")
add("g_all_bars", gd["average"]["g16"], {"mc_stat": max(gd["X_top"]["g16"]["err"], gd["Y_top"]["g16"]["err"])}, "", "MC del banco, muestra de stage5", "idem", "5862882",
    "calculated", "16 barras X; el proceeding de octubre usaba 0.980 (bar_mc_g.json, 1e8 rayos)")

# muones etiquetados después de restar accidentales (q·N)
b0 = open(os.path.join(O, "hod2018_impact_map/2026-11-12_impact_map_runB/logs/b0_checks.log")).read()
for T in (5, 10):
    m = re.search(rf"tiempo bueno: N 77558 T {T}: N_XY (\d+) acc_mix ([\d.]+) acc_ctl (\d+) q ([\d.]+) ± ([\d.]+)", b0)
    q, qe = float(m.group(4)), float(m.group(5))
    add(f"tagged_muons_after_subtraction_{T}pe", q * 77558, qe * 77558, "muones", GOOD, "outputs/hod2018_impact_map/2026-11-12_impact_map_runB/logs/b0_checks.log", "6e25075 (b0_checks.py)",
        "calculated", f"q·N con q = {q} ± {qe} (xy_estimator, mezclados); X·Y en tiempo N_XY = {m.group(1)}, accidentales mezclados {m.group(2)}")
im = json.load(open(os.path.join(O, "hod2018_impact_map/2026-11-12_impact_map_runB/stationary_goodtime/summary.json")))
for T in (5, 10):
    add(f"tagged_muons_impact_sample_{T}pe", im[f"T{T}"]["qN"], im[f"T{T}"]["qN_err"], "muones", "run B, estacionario ∩ tiempo bueno, 76 784 triggers",
        "outputs/hod2018_impact_map/2026-11-12_impact_map_runB/stationary_goodtime/summary.json", "6e25075 (impact_map.py)", "calculated", "q·N de la muestra del perfil de impacto")

# fracción de dos barras
dh = json.load(open(os.path.join(REPO, PK, "work/doublehit_confirmed.json")))["data"]["may8_good"]
for pl in "XY":
    for r in dh[pl]["rows"]:
        if r["T"] in (5, 10):
            add(f"two_bar_fraction_{pl}_{r['T']}pe", r["all"], r["all_err"], "", GOOD + ", etiqueta > 10 p.e.", f"{PK}/work/doublehit_confirmed.json", "b503286 (doublehit_confirmed.py)",
                "measured", "todos los pares vecinos (sin los excluidos: X 3|4, 11|12, 12|13; Y 10|11)")
    add(f"two_bar_fraction_{pl}_MC", dh[pl]["rows"][0]["mc_all_cut0.0"], None, "", "MC del banco (CAD, cos³θ, 3e7 rayos), detección perfecta", f"{PK}/work/doublehit_confirmed.json",
        "b503286", "calculated", "promedio de los dos órdenes de bloque, misma exclusión de pares")

# resolución temporal
for r in csv.DictReader(l for l in open(os.path.join(O, "hod2018_timing/2026-11-03_followup13/A_time/time_res.csv")) if not l.startswith(('#', '"#'))):
    if r["method"] == "cfd50" and r["correction"] == "offsets+prop":
        b = r["min_amp_bin_pe"].replace(">", "gt")
        add(f"sigma_tX_minus_tY_{b}", float(r["sigma_ns"]), float(r["sigma_err_ns"]), "ns", GOOD + ", muones etiquetados X·Y > 10 p.e.",
            "outputs/hod2018_timing/2026-11-03_followup13/A_time/time_res.csv", "4debf58 (time_res.py)", "measured",
            f"CFD 0.5 con desplazamientos por canal y propagación; bin de la menor amplitud {r['min_amp_bin_pe']} p.e.; n = {r['n_same']} (mismo evento)")

# estabilidad diurna de G (el JSON tiene NaN: se lee con el módulo json de Python)
ds = json.load(open(os.path.join(O, "hod_stability/2026-11-03_followup12/B_diurnal/diurnal_summary.json")))
for camp, var, lab in (("may", "drift", "mayo, run B, bins de 1 h en tiempo bueno (34)"), ("oct", "no_drift", "octubre run 8, 18 bins de 1 h (menos de 24 h)")):
    f = ds[camp]["good"]["G_rel"][var]
    add(f"G_diurnal_peak_to_peak_{camp}", f["p2p_frac"] * 100, f["p2p_frac_err"] * 100, "%", lab, "outputs/hod_stability/2026-11-03_followup12/B_diurnal/diurnal_summary.json",
        "4debf58 (diurnal.py)", "measured", f"ajuste c + a·cos(ωh) + b·sin(ωh){' + deriva' if var == 'drift' else ''}; fase (hora local del máximo) "
        f"{f['phase_local_h']:.2f} ± {f['phase_err_h']:.2f} h; errores escalados por √(χ²/ndf)")

# moda del espectro XB
for r in csv.DictReader(open(os.path.join(O, "proceeding_revision/2026-11-02_revision/paper_table.csv"))):
    if r["layer"] == "XB" and r["pos_bin"] == "far":
        add("mode_XB_29_54mm", float(r["mode_pe_meas"]), None, "p.e.", GOOD + ", etiqueta > 10 p.e.", "outputs/proceeding_revision/2026-11-02_revision/paper_table.csv",
            "e0ae42a (mip_paper.C)", "measured", f"máximo de la Landau ⊗ Gauss (4-105 p.e.); n en el ajuste = {r['n_net_in_fit_range']}; χ² {r['chi2']}/{r['ndf']}")
bt = [r for r in csv.DictReader(open(os.path.join(O, "hod2018_may_targets/2026-11-04_bootstrap/boot_fit.csv"))) if r["layer"] == "XB" and r["pos_bin"] == "29-54" and r["rep"] != "0" and r["status"] == "0"]
md = [float(r["mode"]) for r in bt]
mm = sum(md) / len(md)
nums["mode_XB_29_54mm"]["error"] = math.sqrt(sum((x - mm) ** 2 for x in md) / (len(md) - 1))
nums["mode_XB_29_54mm"]["note"] += f"; error = SD de la moda sobre {len(md)} réplicas del bootstrap de eventos (outputs/hod2018_may_targets/2026-11-04_bootstrap/boot_fit.csv, commit 141eeeb)"

os.makedirs(OUT, exist_ok=True)
json.dump({"generated_by": "analysis/sochifi_poster_2026_11_13/make_numbers.py", "numbers": nums}, open(os.path.join(OUT, "poster_numbers.json"), "w"), indent=1, ensure_ascii=False)
for k, v in nums.items():
    print(k, v["value"], v["error"], v["unit"], "|", v["tag"])
