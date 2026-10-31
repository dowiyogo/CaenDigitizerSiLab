#!/usr/bin/python3 -I
# runs_hv.py — seguimiento 2026-10-31 (8), punto 1: V0Set, VMon e IMon de CH00-07 (grupos HV 1-8) y temperatura de la placa A7030DP durante
# cada toma usada para el voltaje de ruptura: los runs SW de los barridos que terminaron bien (2026-10-22 run 1; 2026-10-23 runs 1, 4, 11,
# 18), la toma SW de la calibración c1 (commissioning2_2026_10_21 run 2) y los runs de física de octubre con calibración c1 (runs 4, 7, 8, 9,
# 10 de hodoscope_2026_10).
# Ventana de cada toma: [primer evento, fin] del index.json del run. Slow control de 1 s (hv_YYYYMMDD.csv[.gz]): media y desviación de V0Set,
# VMon e IMon por canal en la ventana; temperatura: bd_Temp_degC de hv_full (60 s), media, mínimo y máximo en la ventana.
# Salida: <out>/runs_hv.json y runs_hv.csv (una fila por toma y canal). Uso: runs_hv.py <out dir>
import csv, glob, gzip, io, json, os, statistics as st, sys, time

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
RAW = "/data/daq/raw"
SC = "/data/daq/slowcontrol/2026"
RUNS = [("scan22_dV0", "biasscan_2026_10_22", 1, "SW", 0.0), ("scan23_dV0", "biasscan_2026_10_23_day", 1, "SW", 0.0),
        ("scan23_dVm0p5", "biasscan_2026_10_23_day", 4, "SW", -0.5), ("scan23_dVm1p0", "biasscan_2026_10_23_day", 11, "SW", -1.0),
        ("scan23_dVm1p5", "biasscan_2026_10_23_day", 18, "SW", -1.5), ("c1_source", "commissioning2_2026_10_21", 2, "SW", 0.0)] + \
       [(f"phys_run{r}", "hodoscope_2026_10", r, "ext", 0.0) for r in (4, 7, 8, 9, 10)]


def index(camp, run):
    f = glob.glob(f"{RAW}/{camp}/{camp}_*_run{run:03d}.index.json")
    assert len(f) == 1, (camp, run, f)
    d = json.load(open(f[0]))
    return d, sorted(glob.glob(f[0].replace(".index.json", "_seg*.root")))


def days(t0, t1):
    out, t = [], t0 - (t0 % 86400)
    while t <= t1:
        out.append(time.strftime("%Y%m%d", time.gmtime(t)))
        t += 86400
    return out


def rows(kind, day):
    for p in (f"{SC}/{kind}_{day}.csv", f"{SC}/{kind}_{day}.csv.gz"):
        if os.path.exists(p):
            f = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
            yield from csv.DictReader(f)
            return


res = []
for label, camp, run, trig, dv in RUNS:
    d, segs = index(camp, run)
    t0, t1 = d["first_event_utc_ns"] / 1e9, d["end_utc_ns"] / 1e9
    acc = {c: {"V0Set_V": [], "VMon_V": [], "IMon_A": []} for c in range(8)}
    temp = []
    for day in days(t0, t1):
        for r in rows("hv", day):
            t = int(r["t_resp_utc_ns"]) / 1e9
            c = int(r["ch"])
            if c < 8 and t0 <= t <= t1:
                for k in acc[c]:
                    acc[c][k].append(float(r[k]))
        for r in rows("hv_full", day):
            t = int(r["t_resp_utc_ns"]) / 1e9
            if t0 - 60 <= t <= t1 + 60 and r["ch"] == "0" and r.get("bd_Temp_degC", "") != "":
                temp.append(float(r["bd_Temp_degC"]))
    e = {"label": label, "campaign": camp, "run": run, "trigger": trig, "delta_V_scan": dv, "n_events": d["n_events"],
         "t0_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)), "t1_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t1)),
         "segments": segs, "board_temp_degC": {"mean": st.mean(temp) if temp else None, "min": min(temp) if temp else None, "max": max(temp) if temp else None, "n": len(temp)},
         "channels": {}}
    for c in range(8):
        a = acc[c]
        e["channels"][c] = {k: {"mean": st.mean(v) if v else None, "sd": st.pstdev(v) if len(v) > 1 else None, "n": len(v)} for k, v in a.items()}
    res.append(e)
    print(f"{label:14s} {camp} run {run:2d} {trig:3s} ΔV {dv:+.1f} n {d['n_events']:6d} {e['t0_utc']}–{e['t1_utc']} T {e['board_temp_degC']['mean']} "
          f"({e['board_temp_degC']['min']}–{e['board_temp_degC']['max']}, n {len(temp)}) | V0Set " +
          " ".join(f"{e['channels'][c]['V0Set_V']['mean']:.2f}" for c in range(8)) + " | VMon " +
          " ".join(f"{e['channels'][c]['VMon_V']['mean']:.3f}" for c in range(8)) + f" | n_sc {e['channels'][0]['VMon_V']['n']}")
json.dump(res, open(OUT + "/runs_hv.json", "w"), indent=1)
with open(OUT + "/runs_hv.csv", "w") as f:
    f.write("label,campaign,run,trigger,delta_V_scan,n_events,t0_utc,t1_utc,board_temp_mean_degC,board_temp_min_degC,board_temp_max_degC,hv_ch,hv_group,"
            "V0Set_mean_V,VMon_mean_V,VMon_sd_V,IMon_mean_A,IMon_sd_A,n_sc\n")
    for e in res:
        for c in range(8):
            ch = e["channels"][c]
            fm = lambda x: "" if x is None else f"{x:.6g}"
            f.write(f"{e['label']},{e['campaign']},{e['run']},{e['trigger']},{e['delta_V_scan']},{e['n_events']},{e['t0_utc']},{e['t1_utc']},"
                    f"{fm(e['board_temp_degC']['mean'])},{fm(e['board_temp_degC']['min'])},{fm(e['board_temp_degC']['max'])},{c},{c + 1},"
                    f"{fm(ch['V0Set_V']['mean'])},{fm(ch['VMon_V']['mean'])},{fm(ch['VMon_V']['sd'])},{fm(ch['IMon_A']['mean'])},{fm(ch['IMon_A']['sd'])},{ch['VMon_V']['n']}\n")
