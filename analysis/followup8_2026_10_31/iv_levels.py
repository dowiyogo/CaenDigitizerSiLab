#!/usr/bin/python3 -I
# iv_levels.py — seguimiento 2026-10-31 (8), punto 2: IMon contra VMon de CH00-07 en régimen estable, del slow control de 1 s (hv_YYYYMMDD.csv.gz),
# en las ventanas de los barridos: 2026-10-22 19:16:56–20:17:49 (DAQ2-14, ΔV = 0), 2026-10-23 13:36:35–17:12:04 (ΔV = 0, −0.5, −1.0, −1.5) y
# 2026-10-30 02:00:00–02:12:00 (+0.25 V en CH00-02, trip de CH03).
# Estable: Pw = 1, Status = 0x1 (encendido, sin rampa ni otros bits), al menos 10 s después del último cambio de V0Set del canal (pedido) y
# 300 s después de un cambio de Pw (variante: 60 s después del cambio de V0Set, para ver la cola de la descarga del PCB, τ ≈ 115 s).
# Por canal y nivel de V0Set (y barrido): mediana y desviación de VMon e IMon y número de muestras. Además, el registro de 1 s de CH00-03 de
# 02:07:40 a 02:09:00 del 2026-10-30 (trip de CH03).
# Salida: <out>/iv_levels.csv, iv_levels.json, trip_ch03_1s.csv. Uso: iv_levels.py <out dir>
import csv, gzip, json, statistics as st, sys, time, calendar

OUT = sys.argv[1]
SC = "/data/daq/slowcontrol/2026"
WIN = [("scan22", "20261022", "2026-10-22T19:16:56", "2026-10-22T20:17:49"), ("scan23", "20261023", "2026-10-23T13:36:35", "2026-10-23T17:12:04"),
       ("scan30", "20261030", "2026-10-30T02:00:00", "2026-10-30T02:12:00")]
ts = lambda s: calendar.timegm(time.strptime(s, "%Y-%m-%dT%H:%M:%S"))
levels, trip = {}, []
for tag, day, a, b in WIN:
    t0, t1 = ts(a), ts(b)
    last = {c: {"v0": None, "tv0": -1e18, "pw": None, "tpw": -1e18} for c in range(8)}
    for r in csv.DictReader(gzip.open(f"{SC}/hv_{day}.csv.gz", "rt")):
        c = int(r["ch"])
        if c > 7:
            continue
        t = int(r["t_resp_utc_ns"]) / 1e9
        v0, pw, stt = float(r["V0Set_V"]), int(r["Pw"]), int(r["Status"])
        L = last[c]
        if L["v0"] is not None and v0 != L["v0"]:
            L["tv0"] = t
        if L["pw"] is not None and pw != L["pw"]:
            L["tpw"] = t
        L["v0"], L["pw"] = v0, pw
        if tag == "scan30" and c <= 3 and ts("2026-10-30T02:07:40") <= t <= ts("2026-10-30T02:09:00"):
            trip.append({"t_utc": time.strftime("%H:%M:%S", time.gmtime(t)), "t_s": t, "ch": c, "V0Set_V": v0, "VMon_V": float(r["VMon_V"]),
                         "IMon_uA": float(r["IMon_A"]) * 1e6, "Status": stt, "Pw": pw})
        if not (t0 <= t <= t1) or pw != 1 or stt != 1 or t - L["tpw"] < 300:
            continue
        for var, wait in (("10s", 10), ("60s", 60)):
            if t - L["tv0"] >= wait:
                k = (var, tag, c, v0)
                levels.setdefault(k, {"V": [], "I": []})
                levels[k]["V"].append(float(r["VMon_V"]))
                levels[k]["I"].append(float(r["IMon_A"]))
rows = []
for (var, tag, c, v0), d in sorted(levels.items()):
    rows.append({"variant": var, "scan": tag, "hv_ch": c, "hv_group": c + 1, "V0Set_V": v0, "n": len(d["V"]), "VMon_median_V": st.median(d["V"]),
                 "VMon_sd_V": st.pstdev(d["V"]), "IMon_median_A": st.median(d["I"]), "IMon_sd_A": st.pstdev(d["I"]),
                 "IMon_first10_median_A": st.median(d["I"][:10]), "IMon_last10_median_A": st.median(d["I"][-10:])})
with open(OUT + "/iv_levels.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
json.dump(rows, open(OUT + "/iv_levels.json", "w"), indent=1)
with open(OUT + "/trip_ch03_1s.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(trip[0]))
    w.writeheader()
    for r in trip:
        w.writerow(r)
for r in rows:
    if r["variant"] == "10s":
        print(f"{r['scan']} CH{r['hv_ch']:02d} V0Set {r['V0Set_V']:.2f} n {r['n']:5d} VMon {r['VMon_median_V']:.3f} IMon {r['IMon_median_A'] * 1e6:.3f} µA "
              f"(sd {r['IMon_sd_A'] * 1e6:.3f}; primeras 10 {r['IMon_first10_median_A'] * 1e6:.3f}, últimas 10 {r['IMon_last10_median_A'] * 1e6:.3f})")
