#!/usr/bin/python3 -I
# hod2019_hv_check.py — XPLUS-4b: V0Set, VMon, IMon y Pw de CH00-07 (los ocho grupos HV; hod2019 tiene canales en todos) durante la
# noche 2026-10-19 y el run 8 de octubre, desde la slow control de hvmon (/data/daq/slowcontrol/2026/hv_*.csv.gz, solo lectura).
# Uso: hod2019_hv_check.py <salida.json>
import gzip, json, statistics, sys

RUNS = {"night_2026_10_19": (1792444429.595, 1792501184.771, ["20261019", "20261020"]),
        "run8_hodoscope_2026_10": (1792700419.379, 1792762587.820, ["20261022", "20261023"])}
out = {}
for run, (t0, t1, days) in RUNS.items():
    vals = {c: {"V0Set_V": [], "VMon_V": [], "IMon_A": [], "Pw": []} for c in range(8)}
    for d in days:
        with gzip.open(f"/data/daq/slowcontrol/2026/hv_{d}.csv.gz", "rt") as f:
            next(f)
            for line in f:
                r = line.split(",")
                t = int(r[0]) * 1e-9
                c = int(r[4])
                if t0 <= t <= t1 and c < 8:
                    vals[c]["V0Set_V"].append(float(r[5]))
                    vals[c]["VMon_V"].append(float(r[7]))
                    vals[c]["IMon_A"].append(float(r[8]))
                    vals[c]["Pw"].append(int(r[10]))
    out[run] = {c: {"n": len(v["VMon_V"]), "V0Set_V_set": sorted(set(v["V0Set_V"])), "VMon_V_median": statistics.median(v["VMon_V"]) if v["VMon_V"] else None,
                    "VMon_V_min": min(v["VMon_V"]) if v["VMon_V"] else None, "VMon_V_max": max(v["VMon_V"]) if v["VMon_V"] else None,
                    "IMon_uA_median": statistics.median(v["IMon_A"]) * 1e6 if v["IMon_A"] else None, "Pw_set": sorted(set(v["Pw"]))} for c, v in vals.items()}
json.dump(out, open(sys.argv[1], "w"), indent=1)
for c in range(8):
    a, b = out["night_2026_10_19"][c], out["run8_hodoscope_2026_10"][c]
    print(f"CH{c:02d} noche: V0Set {a['V0Set_V_set']} VMon {a['VMon_V_median']} IMon {a['IMon_uA_median']:.2f} µA Pw {a['Pw_set']} n {a['n']} | "
          f"run 8: V0Set {b['V0Set_V_set']} VMon {b['VMon_V_median']} IMon {b['IMon_uA_median']:.2f} µA Pw {b['Pw_set']} n {b['n']}")
