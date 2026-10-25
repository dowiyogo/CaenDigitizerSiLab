#!/usr/bin/python3 -I
# stage8.py — PROC-1 etapa 8: cierra el paquete.
#  1. Último cálculo de datos (línea 326 del manuscrito): fracción por trigger con impacto en un plano de hod2019 (35 y 3 mV,
#     run B completo, definición única, estimador de xy_estimator.py) contra f_2019 del MC.
#  2. results.json (package.py), comprobación de tamaño (< 2 MB) y de series (< 200 puntos).
#  3. handover/: figs/, tables/, results.json, figures_manifest.json, deltas.md, report_es.md, MANIFEST.sha256 y README.md
#     con los comandos, el commit y los archivos crudos usados con su SHA256 del índice crudo (SHA256SUMS de cada campaña).
#  4. Convención del monorepo pedida por René para la etapa 8: manifest.json, SHA256SUMS de los archivos pesados de work/ y una
#     línea en outputs/INDEX.md.
# Uso: stage8.py <dir del paquete>
import glob, hashlib, json, os, platform, shutil, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
git_full = subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
branch = subprocess.run(["git", "-C", REPO, "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": 8, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

# ---------------------------------------------------------------- 1. hod2019 por coordenada
e = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
X19, Y19 = list(range(48, 64)), list(range(32, 48))  # kCfg_2026_05: X = ch48-63, Y = ch32-47
f19 = json.load(open(W + "/stage2.json"))["mc_fractions_1e8"]["f"]["h19"]
S["h19_per_coordinate"] = {}
for thr in (35, 3):
    r = xy.xy_fraction((e["aS"][:, X19] > thr).any(1), (e["aS"][:, Y19] > thr).any(1), (e["aC"][:, X19] > thr).any(1), (e["aC"][:, Y19] > thr).any(1))
    S["h19_per_coordinate"][f"{thr}mV"] = {k: r[k] for k in ("N", "N_X", "N_Y", "N_XY", "aX", "aY", "rX", "rY", "f_mix", "f_mix_err", "f_ctl")} | {
        "f_2019_mc": f19, "rX_over_f2019": r["rX"] / f19, "rY_over_f2019": r["rY"] / f19, "rXrY_times_N": r["rX"] * r["rY"] * r["N"]}
    print(thr, "mV:", {k: round(v, 5) if isinstance(v, float) else v for k, v in S["h19_per_coordinate"][f"{thr}mV"].items()})
json.dump(S, open(W + "/stage8.json", "w"), indent=1, default=float)

# ---------------------------------------------------------------- 2. results.json
subprocess.run([sys.executable, "-I", "-B", REPO + "/daq/tools/proceeding/package.py", OUT], check=True)
R = json.load(open(OUT + "/results.json"))
size = os.path.getsize(OUT + "/results.json")
assert size < 2_000_000, size
assert all(len(d["x"]) <= 200 for s in R["series"] for d in s["data"])
stages = sorted({p["stage"] for p in R["provenance"]})
print("results.json", size, "bytes; etapas con provenance:", stages)


# ---------------------------------------------------------------- 3. handover
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


HO = OUT + "/handover"
if os.path.exists(HO):  # no se borra: la versión anterior va al attic
    shutil.move(HO, OUT + "/attic/handover_" + time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()))
    open(OUT + "/attic/README.md", "a").write(f"- handover_*: versión anterior de handover/, reemplazada por stage8.py ({S['generated_utc']}).\n")
os.makedirs(HO)
shutil.copytree(OUT + "/figs", HO + "/figs")
shutil.copytree(OUT + "/tables", HO + "/tables")
for f in ("results.json", "figures_manifest.json", "deltas.md", "report_es.md"):
    shutil.copy2(f"{OUT}/{f}", f"{HO}/{f}")
# archivos crudos usados y su SHA256 del índice de cada campaña
used = {"hodoscope_extnim_2026": [f"lectura_vx1740_ext_nim_1M{s}.root" for s in [""] + [f"_v{i}" for i in range(2, 11)]],
        "hodoscope_2026_10": ["hodoscope_2026_10_20261022T202019Z_run008_seg", "hodoscope_2026_10_20261021T171007Z_run004_seg",
                              "hodoscope_2026_10_20261022T135553Z_run007_seg"],
        "night_2026_10_19": ["night_2026_10_19_20261019T211349Z_run001_seg"],
        "biasscan_2026_10_23_day": ["biasscan_2026_10_23_day_20261023T133835Z_run001_seg", "biasscan_2026_10_23_day_20261023T140147Z_run004_seg"]}
raw_lines, missing = [], []
for camp, pats in used.items():
    idx = {}
    for l in open(f"/data/daq/raw/{camp}/SHA256SUMS"):
        if l.strip() and not l.startswith("#"):
            h, n = l.split(maxsplit=1)
            idx[n.strip().lstrip("*")] = h
    for p in pats:
        hits = [n for n in idx if n == p or (n.startswith(p) and n.endswith(".root"))]
        if not hits:
            missing.append(f"{camp}/{p}")
        raw_lines += [f"{idx[n]}  /data/daq/raw/{camp}/{n}" for n in sorted(hits)]
README = f"""# PROC-1: paquete para el proceeding de Muographers 2026 (JAIS)

Generado por `daq/tools/proceeding/stage8.py` el {S['generated_utc']}. Repositorio `CaenDigitizerSiLab-jai`, rama `{branch}`, commit `{git_full}`.
René copia esta carpeta a `results/proceeding_jais_2026/` del monorepo y corre ahí la sincronización con Overleaf.

## Contenido

- `figs/`: figuras en PDF y PNG (300 ppp), estilo `daq/tools/proceeding/proc_style.C`.
- `tables/`: fragmentos LaTeX.
- `results.json`: todos los números, con su procedencia (provenance por etapa).
- `figures_manifest.json`: pie de figura propuesto, datos y script de cada figura.
- `deltas.md`: cambios línea por línea del manuscrito y respuesta a las marcas `\\verif`.
- `report_es.md`: informe.
- `MANIFEST.sha256`: SHA256 de cada archivo de esta carpeta.

## Comandos (desde la raíz del repositorio; `nice -n 19`, `ulimit -v 8388608`)

```
P=outputs/proceeding_2026_10
INC='gSystem->AddIncludePath("-I{REPO}/daq/core/include")'
# extracción (una vez por archivo; 64 canales)
root -l -b -q -e "$INC" 'daq/tools/characterization/hitdef_extract.C+("<archivo crudo>", "may"|"daqv2", "<G semilla>", "$P/work/hd64_<conjunto>", 64)'
root -l -b -q 'daq/tools/characterization/may_inventory.C(...)'           # inventario de mayo -> work/may_inventory.jsonl
root -l -b -q 'daq/tools/proceeding/manuscript_xy.C(...)'                 # definición del manuscrito -> work/msxy_may8.bin
root -l -b -q -e "$INC" 'daq/tools/proceeding/timing_extract.C+("/data/daq/raw/hodoscope_extnim_2026/lectura_vx1740_ext_nim_1M_v8.root", "<G v8>", "$P/work/tm_may_v8")'
analysis/bin/poisson_analysis --v7 <v7.root> --v8 <v8.root>               # desde $P/work/poisson
python3 -I -B daq/tools/proceeding/bar_mc.py $P/work/bar_mc_g.json 100000000
python3 -I -B daq/tools/characterization/geom_sensitivity.py $P/work/geom_sens_OLD_1e8.json 100000000 OLD   # y CAD
# etapas
python3 -I -B daq/tools/proceeding/stage1.py $P      # … hasta stage7.py; stage5.py $P 30000000; stage6.py $P 100000000
# PROC-2 (versión 2: R1 aceptancia del CAD, R2 g = 0.980)
python3 -I -B daq/tools/characterization/geom_sensitivity.py $P/work/geom_sens_CAD_1e8.json 100000000 CAD
python3 -I -B daq/tools/proceeding/headline.py $P
python3 -I -B daq/tools/proceeding/along_bar.py $P 30000000
python3 -I -B daq/tools/proceeding/h19_xy.py $P
python3 -I -B daq/tools/proceeding/doublehit2.py $P 10000000
python3 -I -B daq/tools/proceeding/trigger_diurnal.py $P
python3 -I -B daq/tools/proceeding/stage8.py $P
# figuras
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_eff.C("$P/work/fig_eff.json","$P/figs/eff_vs_threshold")'
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage3.C("$P/work/fig_stage3.json","$P/figs")'
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage4.C("$P/work/fig_stage4.json","$P/figs/amp_spectra_hod2018_pe")'
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage5.C("$P/work/fig_stage5.json","$P/figs/double_hit_fraction")'
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_stage7.C("$P/work/fig_stage7.json","$P/figs/hod2019_rates")'
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_along_bar.C("$P/work/fig_along_bar.json","$P/figs/eff_along_bar")'
root -l -b -q daq/tools/proceeding/proc_style.C 'daq/tools/proceeding/proc_plot_diurnal.C("$P/work/fig_trigger_diurnal.json","$P/figs/trigger_rate_vs_time")'
```

Entradas derivadas usadas además: `outputs/daq2-15/biasscan_2026_10_23_day/noise_step{{0,1}}.csv` (barrido de día) y
`outputs/characterization_2026_10_23/work/<run>_hourly.csv` (tasas por hora de octubre).

## Archivos crudos usados (SHA256 copiado del `SHA256SUMS` de cada campaña; no recalculado)

```
{chr(10).join(raw_lines)}
```
{"Sin entrada en el índice: " + ", ".join(missing) if missing else "Todos los archivos usados tienen entrada en su índice."}
"""
open(HO + "/README.md", "w").write(README)
files = sorted(p for p in glob.glob(HO + "/**", recursive=True) if os.path.isfile(p) and not p.endswith("MANIFEST.sha256"))
open(HO + "/MANIFEST.sha256", "w").write("".join(f"{sha(p)}  {os.path.relpath(p, HO)}\n" for p in files))
print("handover:", len(files), "archivos; crudos", len(raw_lines), "; sin índice:", missing)

# ---------------------------------------------------------------- 4. convención del monorepo (René, para la etapa 8)
heavy = sorted(p for p in glob.glob(W + "/**", recursive=True) if os.path.isfile(p) and os.path.getsize(p) > 1_000_000)
open(OUT + "/SHA256SUMS", "w").write("".join(f"{sha(p)}  {os.path.relpath(p, OUT)}\n" for p in heavy))
commits = subprocess.run(["git", "-C", REPO, "log", "--format=%h %s", "8027437^..HEAD"], capture_output=True, text=True).stdout.strip().splitlines()
manifest = {"module": "proceeding", "tag": "proceeding_2026_10", "task": "PROC-1 + PROC-2", "package_version": 2, "generated_utc": S["generated_utc"],
            "repo": {"name": "CaenDigitizerSiLab-jai", "branch": branch, "commit": git_full, "proc1_commits": commits},
            "dependencies": {"ROOT": "6.32.10 (~/root6)", "python": platform.python_version(), "numpy": np.__version__,
                             "nlohmann_json": "sistema (/usr/include/nlohmann)", "analysis/bin/poisson_analysis": "binario del repositorio, sin cambios"},
            "parameters": {"hit_definition": "única (DAQ2-15)", "thresholds_pe": list(range(3, 13)), "main_threshold_pe": 5, "hod2019_thresholds_mV": [1, 3, 35],
                           "accidentals": "eventos mezclados (central), ventana fuera de tiempo (syst)", "dead_channel": 3,
                           "G_may": "reajustada con la definición única (A3)", "G_oct": "pe_calibration_2026_10_c1.json"},
            "geometry_version": {"note": "1 = z anteriores (geom_acceptance.py, placas OLD), 2 = CAD 2026-10-24", "used": 2,
                                 "stage6_previous_values": 1, "bars": "CAD (René, 2026-10-24) en bar_mc.py y bench_mc.py",
                                 "acceptance": "R1: f_2018, f_2019, f_both del CAD con 1e8 rayos, cos3", "g": "R2: 0.980 ± 0.03"},
            "seeds": {"bench_mc": 42, "bar_mc": 42, "stage5_mc": 42, "geom_sensitivity": 42, "along_bar_mc": 42, "doublehit2_mc": 42, "h19_bootstrap": 42},
            "duration": "no medida (las etapas se corrieron por separado entre el 2026-10-24 y la fecha de generación)",
            "machine": {"host": platform.node(), "os": platform.platform()},
            "heavy_files": "SHA256SUMS (archivos de work/ de más de 1 MB)", "handover": "handover/ (MANIFEST.sha256)"}
json.dump(manifest, open(OUT + "/manifest.json", "w"), indent=1, ensure_ascii=False)
idx = REPO + "/outputs/INDEX.md"
line = (f"- `proceeding_2026_10/` — PROC-1 + PROC-2 (versión 2: aceptancia del CAD, g = 0.980), paquete del proceeding de Muographers 2026 (JAIS): "
        f"mayo v8 y octubre, figuras, tablas, results.json; rama {branch}, commit {git}, {S['generated_utc'][:10]}\n")
old = open(idx).read() if os.path.exists(idx) else "# Índice de outputs/\n\n"
keep = [l for l in old.splitlines(keepends=True) if not l.startswith("- `proceeding_2026_10/`")]  # la línea de este paquete se reemplaza
open(idx, "w").write("".join(keep) + line)
print("manifest.json, SHA256SUMS (", len(heavy), "archivos ), INDEX.md")
