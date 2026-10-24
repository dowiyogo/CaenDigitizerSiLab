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
json.dump(R, open(OUT + "/results.json", "w"), indent=1, ensure_ascii=False, default=float)
print("results.json", os.path.getsize(OUT + "/results.json"), "bytes;", len(R["results"]), "resultados;", len(R["series"]), "series")
