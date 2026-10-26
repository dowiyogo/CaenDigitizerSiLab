#!/usr/bin/python3 -I
# write_manifests.py — XPLUS-4b cierre: manifest.json y SHA256SUMS por carpeta de salida, y el manifiesto de las ntuplas (entradas crudas
# por tamaño + mtime; SHA256 de mayo del MANIFEST.sha256 de XPLUS-1 D). Uso: write_manifests.py
import glob, hashlib, json, os, subprocess, time

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
B = os.path.join(REPO, "outputs/night_2026-10-25")
ND = "/data/daq/derived/night_2026-10-25/ntuples"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sums(d):
    fs = sorted(os.path.relpath(os.path.join(r, f), d) for r, _, ff in os.walk(d) for f in ff if f != "SHA256SUMS")
    open(os.path.join(d, "SHA256SUMS"), "w").write("".join(f"{sha(os.path.join(d, f))}  {f}\n" for f in fs))
    return len(fs)


C = "analysis/night_2026-10-25/"
STAGES = {
    "stage2": (["dark_stage2.C", "stage2_summary.py", "export_inputs.py"], "ntuplas sw_0V, sw_m0p5V, sw_m1p0V, sw_m1p5V, may_v8 (tiempo bueno, región pre), oct_run8 (pre/post)",
               {"dead_time_samples": 3, "floor_mV": 1.0, "thresholds_pe": [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6], "interval_fit_samples": [40, 150], "short_samples": [8, 31]}),
    "stage3": (["stage3.py", "gain_blocks.C"], "ntuplas may_v8, may_v10, oct_run8 (_hits.bin y árbol spec); slow control hv_202610{22,23}.csv.gz",
               {"min_events_block": 8000, "T_pe": [5, 10], "dlnG_dV_per_V": 1.104, "f_2018": 0.03522, "g": 0.980}),
    "stage4b": (["profile_4b.C", "profile_summary.py", "toy_4b.py", "toy_estimator.py"], "ntuplas may_v8 (tiempo bueno), oct_run8, sw_0V",
                {"T_pe": 5, "window_samples": 9, "toy_replicas": 1000, "toy_seeds": [20261028, 20261029]}),
    "stage5": (["stage5.py", "stage5_table.py", "gain_blocks.C"], "ntuplas may_v1..v10 (_hits.bin y spec); work/stage3.json (orden de tiras de mayo)",
               {"T_pe": [5, 10], "strip_order_min": {"X": 12, "Y": 14}, "f_2018": 0.03522, "g": 0.980}),
    "hod2019_hv": (["hod2019_hv_check.py"], "/data/daq/slowcontrol/2026/hv_2026101{9,20}, hv_202610{22,23}.csv.gz; outputs/characterization_2026_10_23/work/{night,run008}_hourly.csv", {}),
    "plots": (["export_plot_data.py"], "salidas de stage2, stage3, stage4b, stage5", {}),
    "figs": (["plot_day.C", "export_plot_data.py"], "plots/*.csv; daq/tools/proceeding/proc_style.C", {}),
}
for d, (code, inputs, par) in STAGES.items():
    p = os.path.join(B, d)
    if not os.path.isdir(p):
        continue
    json.dump({"stage": f"XPLUS-4b {d}", "code": [C + c for c in code], "code_commit": git, "inputs": inputs, "parameters": par, "generated_utc": now},
              open(os.path.join(p, "manifest.json"), "w"), indent=1, ensure_ascii=False)
    print(d, sums(p), "archivos")
# ntuplas
may = {}
mf = os.path.join(REPO, "outputs/may_raw_manifest_2026-10-25/MANIFEST.sha256")
for l in open(mf):
    if not l.startswith("#"):
        h, size, mt, name = l.split()
        may[name] = {"sha256": h, "size_bytes": int(size), "mtime_utc": mt}
srcs = {"may_v1": "hodoscope_extnim_2026/lectura_vx1740_ext_nim_1M.root"}
for k in range(2, 11):
    srcs[f"may_v{k}"] = f"hodoscope_extnim_2026/lectura_vx1740_ext_nim_1M_v{k}.root"
srcs.update({"oct_run8": "hodoscope_2026_10/hodoscope_2026_10_20261022T202019Z_run008_seg*.root",
             "oct_night": "night_2026_10_19/night_2026_10_19_20261019T211349Z_run001_seg*.root",
             "sw_0V": "biasscan_2026_10_23_day/biasscan_2026_10_23_day_20261023T133835Z_run001_seg*.root",
             "sw_m0p5V": "biasscan_2026_10_23_day/biasscan_2026_10_23_day_20261023T140147Z_run004_seg*.root",
             "sw_m1p0V": "biasscan_2026_10_23_day/biasscan_2026_10_23_day_20261023T150620Z_run011_seg*.root",
             "sw_m1p5V": "biasscan_2026_10_23_day/biasscan_2026_10_23_day_20261023T161052Z_run018_seg*.root"})
ent = {}
for name, pat in srcs.items():
    files = sorted(glob.glob("/data/daq/raw/" + pat))
    ent[name] = {"inputs": [{"path": f, "size_bytes": os.path.getsize(f), "mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(os.path.getmtime(f))),
                             **({"sha256 (XPLUS-1 D)": may[os.path.basename(f)]["sha256"]} if os.path.basename(f) in may else {})} for f in files],
                 "outputs": {os.path.basename(o): {"size_bytes": os.path.getsize(o)} for o in sorted(glob.glob(f"{ND}/{name}.root") + glob.glob(f"{ND}/{name}_hits.bin"))},
                 "log": f"/home/xplus/night_2026-10-25/logs/ntuple_{name}.log" if os.path.exists(f"/home/xplus/night_2026-10-25/logs/ntuple_{name}.log") else "primer plano (prueba de v5/v6)"}
tot = sum(v["size_bytes"] for e in ent.values() for v in e["outputs"].values())
json.dump({"stage": "XPLUS-4b etapa 1 (ntuplas)", "code": [C + "ntuple_extract.C", C + "validate_ntuple.py"], "code_commit": git, "generated_utc": now,
           "validation": "v5, v6, v8 y run 8: _hits.bin idéntico byte a byte a work/hd64_*_hits.bin de PROC; v8 y run 8 reproducen el fixture exactamente",
           "total_bytes": tot, "ntuples": ent}, open(os.path.join(ND, "manifest.json"), "w"), indent=1)
print("ntuplas", len(ent), "total", tot, "B")
