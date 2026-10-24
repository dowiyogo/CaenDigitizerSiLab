#!/usr/bin/python3 -I
# package.py — PROC-1: arma outputs/proceeding_2026_10/results.json (esquema de DAQ2-11/DAQ2-15: provenance, results, series,
# maps, tables) desde work/stage<N>.json de las etapas hechas. Se vuelve a correr al cerrar cada etapa.
# Uso: package.py <dir del paquete>
import json, os, subprocess, sys, time

OUT = sys.argv[1]
W = OUT + "/work"
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
R = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "git_head": git, "provenance": [], "results": [], "series": [], "maps": {},
     "tables": {}}


def res(id_, title, value, stat=None, syst=None, unit="", status="measured", dataset="", hitdef="", threshold="", code="", commit=git, note=""):
    unc = None if stat is None and syst is None else {"stat": stat, "syst": syst}
    R["results"].append({"id": id_, "title_es": title, "value": value, "uncertainty": unc, "unit": unit, "status": status, "dataset": dataset,
                         "hit_definition": hitdef, "threshold": threshold, "code": code, "commit": commit, "note_es": note})


def ser(id_, title, xl, yl, data, status="measured", source="", note=""):
    for d in data:
        assert len(d["x"]) <= 200
    R["series"].append({"id": id_, "title_es": title, "x_label": xl, "y_label": yl, "status": status, "source": source, "note_es": note, "data": data})


SINGLE = "definición única (DAQ2-15): pulso del buscador del DAQ con su mínimo en la ventana signal de 144 ns, pedestal = media recortada de las 40 muestras previas, amplitud > T·G"
# ---------------------------------------------------------------- etapa 1
if os.path.exists(W + "/stage1.json"):
    s1 = json.load(open(W + "/stage1.json"))
    c1 = s1["commit"]
    R["provenance"] += [
        {"stage": 1, "claim": "inventario de los runs de mayo v1-v10", "source_file": "/data/daq/raw/hodoscope_extnim_2026/*.root", "commit": c1,
         "code": "daq/tools/characterization/may_inventory.C; daq/tools/proceeding/stage1.py",
         "method": "árbol data_dig, branch timestamp; UTC con la época de hod2018_good_time_step2.C, validada con la fecha de los archivos", "status": "measured"},
        {"stage": 1, "claim": "calibración de mayo (archivo y reajuste con la definición única), G_mayo/G_c1", "commit": c1, "code": "daq/tools/proceeding/stage1.py",
         "source_file": "outputs/pe_calibration_v8.txt; work/hd64_may_v8_peaks.csv; pe_calibration_2026_10_c1.json", "method": "gain_fit.h sobre pulsos oscuros",
         "status": "measured"},
        {"stage": 1, "claim": "orden de tiras de mayo v8 por coincidencias entre vecinas", "commit": c1,
         "code": "daq/tools/characterization/hitdef_extract.C (nch 64); daq/tools/proceeding/stage1.py",
         "source_file": "/data/daq/raw/hodoscope_extnim_2026/lectura_vx1740_ext_nim_1M_v8.root",
         "method": "exceso en tiempo − mezclados, z > 5, top-2 mutuo; 5 p.e. (G de mayo) en hod2018, 3 mV en hod2019; tiempo bueno",
         "status": "measured (cadenas); derivado (uniones, ch3); orden fijado salvo una inversión global (iluminación uniforme)"}]
    res("v8_triggers", "v8 (run B): triggers", 143814, dataset="v8", code="may_inventory.C", commit=c1)
    res("v8_duration_h", "v8 (run B): duración", 71.83, unit="h", dataset="v8", code="may_inventory.C", commit=c1)
    res("v7_triggers", "v7 (run A): triggers", 17005, dataset="v7", code="may_inventory.C", commit=c1)
    res("v7_duration_h", "v7 (run A): duración", 8.29, unit="h", dataset="v7", code="may_inventory.C", commit=c1)
    res("v8_good_time_events", "v8: triggers en tiempo bueno", s1["v8_good_time_events"], dataset="v8", code="may_common.good_time", commit=c1)
    cs = s1["calibration_summary"]
    res("may_G_refit_over_file", "G reajustada / G del archivo de mayo, mediana", cs["ratio_fit_over_file_median"], dataset="v8", code="stage1.py", commit=c1)
    res("may_p0_refit_median", "p0 reajustado, mediana", cs["p0_fit_median_mV"], unit="mV", dataset="v8", code="stage1.py", commit=c1)
    res("G_may_over_c1_median", "G mayo (archivo) / G c1, mediana", cs["ratio_may_c1_median"], dataset="v8, c1", code="stage1.py", commit=c1)
    R["tables"]["may_runs"] = s1["inventory"]
    R["tables"]["calibration"] = s1["calibration"]
    R["tables"]["strip_order_may_v8"] = s1["strip_order"]
    R["maps"]["neighbour_coincidences_may_v8"] = s1["strip_order_matrices"]
# ---------------------------------------------------------------- etapa 2
if os.path.exists(W + "/stage2.json"):
    s2 = json.load(open(W + "/stage2.json"))
    c2 = s2["commit"]
    R["provenance"].append({"stage": 2, "claim": "accidentales y eficiencia de hod2018", "commit": c2,
                            "code": "daq/tools/proceeding/stage2.py; daq/tools/characterization/xy_estimator.py; daq/tools/proceeding/bar_mc.py",
                            "source_file": "work/hd64_*_hits.bin (hitdef_extract.C, nch 64), work/msxy_may8.bin (manuscript_xy.C), work/bar_mc_g.json",
                            "method": SINGLE + "; accidentales: eventos mezclados (central, D3) y ventana fuera de tiempo (syst); ch3 fuera",
                            "status": "measured (q); estimated (÷ aceptancia del MC, ÷ g)"})
    for key, d in s2["datasets"].items():
        for T in ("5", "10"):
            r = d["per_T"][T]
            res(f"q_{key}_T{T}", f"{d['label']}: fracción X·Y por trigger (q, mezclados), T = {T} p.e.", r["f_mix"], r["f_mix_err"], r["syst_estimator"],
                dataset=d["label"], hitdef=SINGLE, threshold=f"{T} p.e.", code="stage2.py", commit=c2)
            e = r["eps_XY"]
            res(f"epsXY_{key}_T{T}", f"{d['label']}: ε_XY = q / (f_2018 · g), T = {T} p.e.", e["value"], e["stat"],
                (e["syst_estimator"] ** 2 + e["syst_acceptance"] ** 2) ** 0.5, status="estimated", dataset=d["label"], hitdef=SINGLE,
                threshold=f"{T} p.e.", code="stage2.py; bar_mc.py", commit=c2,
                note="f_2018 = 0.03646 ± 0.00026 (cos³θ), g = P(X∧Y | cruza el cuadrado) del MC de barras; syst = estimador ⊕ aceptancia ⊕ g")
        ser(f"q_vs_T_{key}", f"{d['label']}: q contra T", "T (p.e.)", "q", [
            {"name": "mezclados", "x": [int(T) for T in d["per_T"]], "y": [d["per_T"][T]["f_mix"] for T in d["per_T"]],
             "yerr": [d["per_T"][T]["f_mix_err"] for T in d["per_T"]]},
            {"name": "fuera de tiempo", "x": [int(T) for T in d["per_T"]], "y": [d["per_T"][T]["f_ctl"] for T in d["per_T"]],
             "yerr": [d["per_T"][T]["f_ctl_err"] for T in d["per_T"]]}], source="stage2.py")
        ser(f"aX_aY_rX_rY_{key}", f"{d['label']}: aX, aY (fuera de tiempo) y rX, rY (verdaderas) contra T", "T (p.e.)", "probabilidad",
            [{"name": n, "x": [int(T) for T in d["per_T"]], "y": [d["per_T"][T][n] for T in d["per_T"]], "yerr": None} for n in ("aX", "aY", "rX", "rY")],
            source="stage2.py")
    p = s2["published_convention"]
    res("published_convention_v8", "v8, convención del taller (p0 + 3.5 G, definición de mayo): q / 0.0360", p["eps_0.0360_mix"], p["stat_0.0360"],
        abs(p["eps_0.0360_mix"] - p["eps_0.0360_ctl"]), status="estimated", dataset="v8 tiempo bueno", hitdef="definición de mayo", threshold="p0 + 3.5 G",
        code="stage2.py", commit=c2, note="handoff: 0.88 ± 0.04 ± 0.11 (v7-v10); con 0.03646: %.4f" % p["eps_0.03646_mix"])
    m = s2["manuscript_8073"]["run B completo (como el manuscrito)"]
    res("ms8073_true_mix", "Run B, definición del manuscrito (50 mV, 4000-4400 ns): coincidencias X·Y verdaderas (q·N, mezclados)", m["true_XY_mix"],
        m["err_true_mix"], abs(m["true_XY_mix"] - m["true_XY_ctl"]), dataset="v8 completo", hitdef="manuscrito", threshold="50 mV", code="manuscript_xy.C; stage2.py",
        commit=c2, note=f"N_XY = {m['N_XY']}; accidentales {m['accidental_in_NXY_mix']:.0f} (mezclados) y {m['accidental_in_NXY_ctl']:.0f} (fuera de tiempo)")
    R["tables"]["stage2"] = {k: s2[k] for k in ("geometry_factors", "published_convention", "manuscript_8073", "stability", "exclude_noisy", "G_check",
                                               "mc_fractions_1e8", "counts_table")}
# ---------------------------------------------------------------- etapa 3
if os.path.exists(W + "/stage3.json"):
    s3 = json.load(open(W + "/stage3.json"))
    c3 = s3["commit"]
    R["provenance"].append({"stage": 3, "claim": "mapas X × Y, ocupación por tira y multiplicidad de hod2018; ocupación de hod2019", "commit": c3,
                            "code": "daq/tools/proceeding/stage3.py; proc_plot_stage3.C; xy_estimator.xy_map; bar_mc.py (cruces por tira)",
                            "source_file": "work/hd64_may_v8_hits.bin, work/hd64_oct8_hits.bin, work/bar_mc_g.json",
                            "method": SINGLE + "; por tira en tiempo − fuera de tiempo; mapa en tiempo − mezclados + término rX·rY; orden de mayo de la etapa 1, "
                                               "octubre r4", "status": "measured (cuentas); estimated (banda del MC)"})
    for key in ("may8_good", "oct8"):
        d = s3["occupancy_2018"][key]
        for pl in "XY":
            for T in (5, 10):
                r = d[f"{pl}_T{T}"]["ratio"]
                res(f"occ_maxmin_{key}_{pl}_T{T}", f"{d['label']}: ocupación por tira de hod2018 {pl}, máx/mín (en tiempo − fuera de tiempo)", r["max_over_min"],
                    r["err"], dataset=d["label"], hitdef=SINGLE, threshold=f"{T} p.e.", code="stage3.py", commit=c3,
                    note=f"coordenadas {r['max_coord']} / {r['min_coord']}; MC {d[pl + '_mc']['max_over_min']:.3f}; chi2 contra MC "
                         f"{d[f'{pl}_T{T}']['chi2_vs_mc']:.1f}/{d[f'{pl}_T{T}']['ndf_vs_mc']}")
        m = s3["maps"][key]
        res(f"map_qN_{key}", f"{m['label']}: suma del mapa X × Y restado (= q_mix·N)", m["sum_sub_corrected"], dataset=m["label"], hitdef=SINGLE,
            threshold="5 p.e.", code="stage3.py", commit=c3, note=f"término rX·rY devuelto {m['rXrY_term_total']:.1f}; celdas negativas {m['sub_corrected_cells_negative']}")
    for k, v in s3["mc_strip_spread"].items():
        res(f"mc_strip_maxmin_{k}", f"MC de barras: cruces por tira, máx/mín ({k})", v["max_over_min"], status="estimated", dataset="MC 1e8 rayos, cos³θ, OLD",
            code="bar_mc.py", commit=c3, note=f"desviación relativa {v['rel_sd']:.4f}")
    hc = s3["hod2019_claims"]
    res("h19_top4_kCfg", "mayo v8 (run B): las 4 tiras X de hod2019 con más impactos a 35 mV (coordenada kCfg_2026_05)",
        hc["reading_k_eq_kCfg_coordinate"]["top4"], dataset="v8 completo", hitdef=SINGLE, threshold="35 mV", code="stage3.py", commit=c3,
        note="cuentas " + str([int(x) for x in hc["reading_k_eq_kCfg_coordinate"]["counts"]]) + "; orden de hod2019 no validado")
    z = s3["occupancy_2019"]["may8_all"]
    res("h19_strips_3sigma_35mV", "mayo v8 (run B): tiras de hod2019 con exceso en tiempo > 3σ a 35 mV (X, Y)",
        [len(z["X_35mV"]["strips_above_3sigma"]), len(z["Y_35mV"]["strips_above_3sigma"])], dataset="v8 completo", hitdef=SINGLE, threshold="35 mV",
        code="stage3.py", commit=c3, note="fuera de tiempo: 0 cuentas en todas las tiras")
    for ch, v in s3["oct_ch39_ch47"].items():
        res(f"oct8_ch{ch}_3mV_net", f"octubre run 8: exceso en tiempo de ch{ch} a 3 mV", v["3mV"]["net"], v["3mV"]["net_err"], dataset="octubre run 8",
            hitdef=SINGLE, threshold="3 mV", code="stage3.py", commit=c3, note=f"z = {v['3mV']['z']:.2f}; 35 mV: {v['35mV']['in_time']:.0f} en tiempo")
    for key in ("may8_good", "oct8"):
        d = s3["occupancy_2018"][key]
        for pl in "XY":
            ser(f"occ_{key}_{pl}", f"{d['label']}: impactos por 1000 triggers por tira de hod2018 {pl}", "coordenada", "impactos / 1000 triggers",
                [{"name": f"{T} p.e.", "x": list(range(16)), "y": d[f"{pl}_T{T}"]["net_per_1000"], "yerr": d[f"{pl}_T{T}"]["net_per_1000_err"]} for T in (5, 10)],
                source="stage3.py")
    for pl in "XY":
        ser(f"occ_h19_may8_{pl}", f"mayo v8 (run B): impactos por 1000 triggers por tira de hod2019 {pl}", "coordenada (kCfg_2026_05)", "impactos / 1000 triggers",
            [{"name": f"{T} mV", "x": list(range(16)), "y": z[f"{pl}_{T}mV"]["net_per_1000"], "yerr": z[f"{pl}_{T}mV"]["net_per_1000_err"]} for T in (35, 3)],
            source="stage3.py")
    R["maps"]["hod2018_xy_stage3"] = {k: {a: s3["maps"][k][a] for a in ("label", "sub_corrected", "mixed", "in_time")} for k in s3["maps"]}
    R["tables"]["stage3"] = {k: s3[k] for k in ("strip_orders", "multiplicity", "hod2019_claims", "oct_ch39_ch47", "occupancy_2018_50mV_proxy", "mc_strip_spread")}
json.dump(R, open(OUT + "/results.json", "w"), indent=1, ensure_ascii=False, default=float)
print("results.json", os.path.getsize(OUT + "/results.json"), "bytes;", len(R["results"]), "resultados;", len(R["series"]), "series")
