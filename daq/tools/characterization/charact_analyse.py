#!/usr/bin/python3 -I
# charact_analyse.py — DAQ2-15, parte B (solo lectura): primera caracterización del detector desde datos ya tomados.
# Entradas (work/, de charact_extract.C y charact_gain.C): <run>_events.bin, <run>_hourly.csv, <run>_ootspec.csv,
# <run>_q1.csv, gainwin_<run>.csv; G por hora de daqoffline (/data/daq/sim/offline/charact_<run>/gains.json);
# control lento hv_full_*.csv(.gz) (VMon, IMon, temperatura de la placa A7030DP); mapeo r4; calibración c1.
# Salida: <out>/results.json (esquema de DAQ2-11: provenance, results, series, episodes), <out>/figs.json (lo
# dibuja charact_plots.C) y tablas para report_es.md (<out>/tables.json).
# Convenciones: hora de Chile = UTC − 3 (octubre); día = 12-17 h Chile, noche = 00-06 h Chile. Errores: Poisson
# para conteos, binomial para fracciones; los de las razones por propagación. Estado: measured | estimated.
# Uso: charact_analyse.py <dir de salida>
import calendar, csv, glob, gzip, json, math, os, statistics, subprocess, sys, time

import numpy as np

OUT = sys.argv[1]
W = OUT + "/work"
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
MAP = os.path.expanduser("~/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r4.json")
C1 = os.path.expanduser("~/.local/opt/daqv2/share/calibration/pe_calibration_2026_10_c1.json")
RAW = {"run008": "/data/daq/raw/hodoscope_2026_10/hodoscope_2026_10_20261022T202019Z_run008",
       "night": "/data/daq/raw/night_2026_10_19/night_2026_10_19_20261019T211349Z_run001",
       "run004": "/data/daq/raw/hodoscope_2026_10/hodoscope_2026_10_20261021T171007Z_run004",
       "run007": "/data/daq/raw/hodoscope_2026_10/hodoscope_2026_10_20261022T135553Z_run007"}
NAME = {"run008": "hodoscope_2026_10 run 8", "night": "night_2026_10_19 run 1", "run004": "hodoscope_2026_10 run 4",
        "run007": "hodoscope_2026_10 run 7"}
ACC_GEOM = 0.0360   # fracción de triggers externos que cruza el cuadrado activo de hod2018 (MC; handoff línea 125)
ACC_GEOM_SRC = "docs/handoff_analysis_2026-10-18.md:125 (MC geométrico: 3.60 % de los triggers cruza el área activa)"
MAY = {"value": 0.88, "stat": 0.04, "syst": 0.11, "src": "docs/handoff_analysis_2026-10-18.md:125 (v7-v10, 5 p.e., tiempo bueno)"}
NOISY18 = [2, 9, 15]  # fuera del auto-trigger en DAQ2-13 (ruido: a 10 p.e. 49565, 28980 y 19072 Hz contra mediana 2858 Hz)
CHILE = -3 * 3600
DT = np.dtype([("t", "<f8"), ("s", "<f4", 64), ("c", "<f4", 64), ("q", "<f4", 64), ("k", "u1", 64)])
SIG = (49, 58)

prov, results, series, episodes, figs, tables = [], [], [], [], [], {}


def utc(t):
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t))


def chile(t):
    return time.strftime("%Y-%m-%d %H:%M", time.gmtime(t + CHILE))


def r6(x):
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else round(float(x), 6)


def res(id_, title, value, unc, unit, status, source, note=""):
    results.append({"id": id_, "title_es": title, "value": r6(value), "uncertainty": unc, "unit": unit, "status": status, "source": source,
                    "note_es": note})


def ser(id_, title, xl, yl, xu, yu, kind, status, source, data, note=""):
    for d in data:  # como mucho ~200 puntos por curva
        assert len(d["x"]) <= 210, (id_, d["name"], len(d["x"]))
        d["x"] = [r6(v) if isinstance(v, float) else v for v in d["x"]]
        d["y"] = [r6(v) for v in d["y"]]
        if d.get("yerr") is not None:
            d["yerr"] = [r6(v) for v in d["yerr"]]
    series.append({"id": id_, "title_es": title, "x_label": xl, "y_label": yl, "x_unit": xu, "y_unit": yu, "kind": kind, "status": status,
                   "source": source, "note_es": note, "data": data})


# ---------------------------------------------------------------- entradas
M = json.load(open(MAP))
CH = {c["dig_ch"]: c for c in M["channels"]}
cal = json.load(open(C1))
G = np.full(64, np.nan)
for e in cal["channels"]:
    if e["ch"] < 32 and isinstance(e.get("G_mV"), (int, float)):
        G[e["ch"]] = e["G_mV"]
G_FALLBACK = [c for c in range(32) if not np.isfinite(G[c])]
G[:32][~np.isfinite(G[:32])] = np.median(G[:32][np.isfinite(G[:32])])
PLANES = {}
for st in ("hod2018", "hod2019"):
    for pl in ("X", "Y"):
        chs = sorted((c for c in CH.values() if c["station"] == st and c["plane"] == pl), key=lambda c: c["strip"])
        PLANES[(st, pl)] = [c["dig_ch"] for c in chs]  # índice = coordenada (strip) del mapeo
STRIP_CONF = {f"{st} {pl}": dict(__import__("collections").Counter(CH[c]["confirmed"]["strip"] for c in chs)) for (st, pl), chs in PLANES.items()}
ev = {k: np.fromfile(f"{W}/{k}_events.bin", dtype=DT) for k in RAW}
for k in ev:
    ev[k] = ev[k][np.argsort(ev[k]["t"], kind="stable")]
idx = {k: json.load(open(RAW[k] + ".index.json")) for k in RAW}
git_head = subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()


def hourly(k):
    rows = list(csv.DictReader(open(f"{W}/{k}_hourly.csv")))
    return [{kk: (float(v) if kk not in ("ch",) else int(v)) for kk, v in r.items()} for r in rows]


HR = {k: hourly(k) for k in RAW}


def pe(e):
    return e["s"][:, :32] / G[:32], e["c"][:, :32] / G[:32]


# ---------------------------------------------------------------- 1. resumen de datos
def gaps(t, thr=20.0):
    d = np.diff(t)
    i = np.where(d > thr)[0]
    return [(float(t[j]), float(t[j + 1]), float(d[j])) for j in i]


over = {}
for k, e in ev.items():
    t = e["t"]
    g = gaps(t)
    span = float(t[-1] - t[0])
    live = span - sum(x[2] for x in g)
    over[k] = {"run": NAME[k], "events": int(len(t)), "first_utc": utc(t[0]), "last_utc": utc(t[-1]), "span_s": span, "live_s": live,
               "gaps_gt20s": [{"start_utc": utc(a), "end_utc": utc(b), "dur_s": round(d, 1)} for a, b, d in g],
               "rate_Hz": len(t) / live if live > 0 else float("nan"), "rate_err_Hz": math.sqrt(len(t)) / live if live > 0 else float("nan")}
    for a, b, d in g:
        episodes.append({"id": f"gap_{k}_{utc(a)}", "start": utc(a), "end": utc(b), "duration": f"{d:.1f} s", "cause": "sin triggers externos (causa no registrada)",
                         "status_times": "measured", "status_cause": "unknown", "source": RAW[k] + "_seg*.root (t_utc_ns)"})
tables["overview"] = over
prov.append({"claim": "eventos, tiempo vivo, huecos > 20 s y tasa de trigger por run", "source_file": "; ".join(RAW[k] + "_seg*.root" for k in RAW),
             "commit": git_head, "code": "daq/tools/characterization/charact_extract.C + charact_analyse.py",
             "method": "measured: t_utc_ns (TTT ajustado por run); tiempo vivo = primer a último evento menos los huecos > 20 s",
             "inputs": "árbol events", "plots": ["rate_vs_time.png"], "status": "measured",
             "assumptions": ["huecos > 20 s = sin trigger (a 0.7 Hz, P(Δt > 20 s) = e^-14)"]})
# tasa por hora (run 8 y noche): conteo / parte de la hora cubierta por el run (sin huecos)
rate_series = []
for k in ("run008", "night"):
    t = ev[k]["t"]
    h0 = math.floor(t[0] / 3600) * 3600
    xs, ys, es = [], [], []
    g = gaps(t)
    for h in np.arange(h0, t[-1], 3600):
        a, b = max(h, t[0]), min(h + 3600, t[-1])
        cov = b - a - sum(max(0, min(b, y) - max(a, x)) for x, y, _ in g)
        n = int(((t >= h) & (t < h + 3600)).sum())
        if cov > 300:
            xs.append(chile(h)[5:])
            ys.append(n / cov)
            es.append(math.sqrt(n) / cov)
    rate_series.append({"name": NAME[k], "x": xs, "y": ys, "yerr": es})
    tables.setdefault("rate_per_hour", {})[k] = list(zip(xs, [round(y, 4) for y in ys], [round(x, 4) for x in es]))
ser("trigger_rate_hour", "Tasa de trigger externo por hora", "hora de Chile (inicio)", "tasa", "", "Hz", "points", "measured",
    "run 8 y night_2026_10_19 (t_utc_ns)", rate_series, "conteo / tiempo cubierto de la hora, sin huecos > 20 s")
figs.append({"file": "rate_vs_time.png", "title": "Tasa de trigger externo por hora", "xlabel": "hora de Chile", "ylabel": "Hz", "kind": "points",
             "series": rate_series})


# control lento: hv_full (10 s) por bins de 10 min, CH00-07 y temperatura de la placa
def read_full(t0, t1):
    rows = []
    days = sorted({time.strftime("%Y%m%d", time.gmtime(x)) for x in np.arange(t0 - 86400, t1 + 86400, 3600)})
    for d in days:
        for f in sorted(glob.glob(f"/data/daq/slowcontrol/{d[:4]}/hv_full_{d}*.csv*")):
            if f.endswith(".schema.json"):
                continue
            op = gzip.open if f.endswith(".gz") else open
            with op(f, "rt") as fh:
                rd = csv.reader(l for l in fh if not l.startswith("#"))
                hd = next(rd)
                ix = {h: i for i, h in enumerate(hd)}
                for r in rd:
                    try:
                        ts = int(r[ix["t_req_utc_ns"]]) / 1e9
                    except (ValueError, IndexError):
                        continue
                    if t0 <= ts < t1 and int(r[ix["ch"]]) < 8:
                        rows.append((ts, int(r[ix["ch"]]), float(r[ix["VMon_V"]]), float(r[ix["IMon_A"]]), float(r[ix["V0Set_V"]]), int(r[ix["Pw"]]),
                                     float(r[ix["bd_Temp_degC"]]) if r[ix["bd_Temp_degC"]] not in ("", "nan") else float("nan")))
    return rows


def sc_bins(t0, t1, w=600):
    rows = read_full(t0, t1)
    b = {}
    for ts, ch, vm, im, v0, pw, tp in rows:
        k = int(ts // w) * w
        d = b.setdefault(k, {"T": [], "ch": {}})
        d["T"].append(tp)
        d["ch"].setdefault(ch, []).append((vm, im, v0, pw))
    out = []
    for k in sorted(b):
        d = b[k]
        T = [x for x in d["T"] if x == x]
        out.append({"t": k, "T": statistics.mean(T) if T else float("nan"),
                    "ch": {ch: (statistics.mean(x[0] for x in v), statistics.mean(x[1] for x in v), v[-1][2], min(x[3] for x in v)) for ch, v in d["ch"].items()}})
    return out


SC = {k: sc_bins(ev[k]["t"][0], ev[k]["t"][-1]) for k in RAW}
sc8 = SC["run008"]
step = max(1, len(sc8) // 200 + 1)
sel = sc8[::step]
ser("hv_imon_run8", "IMon por grupo HV durante el run 8 (medias de 10 min)", "hora de Chile", "IMon", "", "µA", "lines", "measured",
    "/data/daq/slowcontrol/2026/hv_full_2026102[23].csv*",
    [{"name": f"CH{ch:02d} (grupo {ch + 1})", "x": [chile(b["t"])[5:] for b in sel], "y": [b["ch"].get(ch, (0, float("nan")))[1] * 1e6 for b in sel], "yerr": None}
     for ch in range(8)])
ser("hv_vmon_run8", "VMon por grupo HV durante el run 8 (medias de 10 min)", "hora de Chile", "VMon", "", "V", "lines", "measured",
    "/data/daq/slowcontrol/2026/hv_full_2026102[23].csv*",
    [{"name": f"CH{ch:02d}", "x": [chile(b["t"])[5:] for b in sel], "y": [b["ch"].get(ch, (float("nan"),))[0] for b in sel], "yerr": None} for ch in range(8)])
ser("board_temp_run8", "Temperatura de la placa A7030DP durante el run 8 (única referencia de temperatura; sigue la del rack)", "hora de Chile",
    "temperatura", "", "°C", "lines", "measured", "/data/daq/slowcontrol/2026/hv_full_2026102[23].csv* (bd_Temp_degC, resolución 1 °C)",
    [{"name": "A7030DP", "x": [chile(b["t"])[5:] for b in sel], "y": [b["T"] for b in sel], "yerr": None}])
figs.append({"file": "hv_imon_run8.png", "title": "IMon por grupo HV, run 8 (10 min)", "xlabel": "hora de Chile", "ylabel": "#muA", "kind": "lines",
             "series": series[-3]["data"]})
figs.append({"file": "hv_vmon_run8.png", "title": "VMon por grupo HV, run 8 (10 min)", "xlabel": "hora de Chile", "ylabel": "V", "kind": "lines",
             "series": series[-2]["data"]})
figs.append({"file": "board_temp_run8.png", "title": "Temperatura de la placa A7030DP, run 8", "xlabel": "hora de Chile", "ylabel": "#circC", "kind": "lines",
             "series": series[-1]["data"]})
hv_sum = {}
for ch in range(8):
    im = [b["ch"][ch][1] * 1e6 for b in sc8 if ch in b["ch"]]
    vm = [b["ch"][ch][0] for b in sc8 if ch in b["ch"]]
    v0 = sorted({b["ch"][ch][2] for b in sc8 if ch in b["ch"]})
    pw = min(b["ch"][ch][3] for b in sc8 if ch in b["ch"])
    hv_sum[ch] = {"IMon_uA_min": min(im), "IMon_uA_median": statistics.median(im), "IMon_uA_max": max(im), "VMon_V_min": min(vm), "VMon_V_max": max(vm),
                  "V0Set_V": v0, "Pw_min": pw}
TT = [b["T"] for b in sc8 if b["T"] == b["T"]]
tables["hv_run8"] = {"per_channel": hv_sum, "board_T_min": min(TT), "board_T_max": max(TT), "n_bins": len(sc8)}

# ---------------------------------------------------------------- 2. por canal
def chan_table(k, hours=None):
    rows = [r for r in HR[k] if hours is None or hours(r["hour_start_utc_s"])]
    out = {}
    for c in range(64):
        rr = [r for r in rows if r["ch"] == c]
        n = sum(r["events"] for r in rr)
        ex = sum(r["exposure_s"] for r in rr)
        cnt = [sum(r[f"cnt{i}"] for r in rr) for i in range(4)]
        ps, pq, rs = (sum(r[x] for r in rr) for x in ("ped_sum_mV", "ped_sumsq_mV2", "rms_sum_mV"))
        m = ps / n if n else float("nan")
        out[c] = {"events": n, "exposure_s": ex, "counts": cnt, "thr_mV": [rr[0][f"thr{i}_mV"] for i in range(4)] if rr else [],
                  "rate_Hz": [x / ex if ex else float("nan") for x in cnt], "rate_err_Hz": [math.sqrt(x) / ex if ex else float("nan") for x in cnt],
                  "ped_mV": m, "ped_sd_mV": math.sqrt(max(0.0, pq / n - m * m)) if n else float("nan"), "ped_rms_mV": rs / n if n else float("nan")}
    return out


CT = {k: chan_table(k) for k in RAW}
# G por hora (fase 5) de daqoffline
gains = {}
for k in RAW:
    p = f"/data/daq/sim/offline/charact_{k}/gains.json"
    if os.path.exists(p):
        gains[k] = [g for g in json.load(open(p)) if g["ch"] < 32]
ghour = {}
for c in range(32):
    v = [g["value"] for g in gains["run008"] if g["ch"] == c and g["status"] in ("ok", "changed")]
    st = __import__("collections").Counter(g["status"] for g in gains["run008"] if g["ch"] == c)
    ghour[c] = {"n_ok": len(v), "median": statistics.median(v) if v else float("nan"), "min": min(v) if v else float("nan"),
                "max": max(v) if v else float("nan"), "sd": statistics.stdev(v) if len(v) > 1 else float("nan"), "status_counts": dict(st)}
# pulsos de 1 p.e.: integral media (para pasar la integral en tiempo a p.e.)
q1 = {int(r["ch"]): float(r["q1_mean_mVns"]) for r in csv.DictReader(open(f"{W}/run008_q1.csv"))}
e8 = ev["run008"]
pe8, cpe8 = pe(e8)


def hist(x, lo, hi, nb):
    h, _ = np.histogram(x, bins=nb, range=(lo, hi))
    return h.tolist()


AMP18 = (0, 40, 80)       # p.e., 0.5 p.e.
AMP19 = (0, 100, 100)     # mV, 1 mV
INT18 = (0, 80, 80)       # p.e. (integral / integral de 1 p.e.)
INT19 = (0, 4000, 80)     # mV·ns
spec_amp, spec_int, flags = {}, {}, {}
st_med = {}
for st, r in (("hod2018", range(32)), ("hod2019", range(32, 64))):
    i = 3 if st == "hod2018" else 1  # 5 p.e. o 3 mV
    st_med[st] = statistics.median(CT["run008"][c]["rate_Hz"][i] for c in r if CT["run008"][c]["counts"][i] > 0)
    st_med[st + "_low"] = statistics.median(CT["run008"][c]["rate_Hz"][0] for c in r)
for c in range(64):
    t = CT["run008"][c]
    if c < 32:
        a, b = pe8[:, c], cpe8[:, c]
        spec_amp[c] = {"in": hist(a, *AMP18), "ctl": hist(b, *AMP18)}
        qq = e8["q"][:, c] / q1[c] if q1.get(c) and q1[c] == q1[c] and q1[c] > 0 else np.full(len(e8), np.nan)
        spec_int[c] = hist(qq[np.isfinite(qq)], *INT18)
        i = 3
    else:
        spec_amp[c] = {"in": hist(e8["s"][:, c], *AMP19), "ctl": hist(e8["c"][:, c], *AMP19)}
        spec_int[c] = hist(e8["q"][:, c], *INT19)
        i = 1
    test = CH[c].get("hv_group_test", {})
    reason = test.get("reason", "")
    if t["rate_Hz"][0] < 0.01 * st_med[("hod2018" if c < 32 else "hod2019") + "_low"]:
        f = "dead"
    elif "no dependen de la HV" in reason or "no depende de la HV" in reason:
        f = "no_hv_dependence"
    elif t["rate_Hz"][i] > 3 * st_med["hod2018" if c < 32 else "hod2019"]:
        f = "noisy"
    elif test.get("result") != "medido":
        f = "ok (grupo HV sin medir)"
    else:
        f = "ok"
    flags[c] = f
per_channel = {}
for c in range(64):
    t = CT["run008"][c]
    d = {"station": CH[c]["station"], "plane": CH[c]["plane"], "strip": CH[c]["strip"], "hv_group": CH[c]["hv_group"],
         "pedestal_mV": t["ped_mV"], "pedestal_sd_mV": t["ped_sd_mV"], "pedestal_rms_mV": t["ped_rms_mV"], "flag": flags[c],
         "dark": [{"thr": lab, "thr_mV": t["thr_mV"][i], "rate_Hz": t["rate_Hz"][i], "err_Hz": t["rate_err_Hz"][i], "counts": t["counts"][i]}
                  for i, lab in enumerate(["0.5 p.e.", "1.5 p.e.", "3 p.e.", "5 p.e."] if c < 32 else ["1 mV", "3 mV"])],
         "exposure_s": t["exposure_s"]}
    if c < 32:
        d["G_c1_mV"] = float(G[c])
        d["G_c1_fallback"] = c in G_FALLBACK
        d["G_hourly"] = ghour[c]
        d["q1_mVns"] = q1.get(c)
    per_channel[c] = d
tables["per_channel"] = per_channel
for st, r, labs in (("hod2018", range(32), ["0.5 p.e.", "3 p.e.", "5 p.e."]), ("hod2019", range(32, 64), ["1 mV", "3 mV"])):
    data = []
    for lab in labs:
        i = ["0.5 p.e.", "1.5 p.e.", "3 p.e.", "5 p.e."].index(lab) if st == "hod2018" else ["1 mV", "3 mV"].index(lab)
        data.append({"name": f"> {lab}", "x": list(r), "y": [CT["run008"][c]["rate_Hz"][i] / 1e3 for c in r],
                     "yerr": [CT["run008"][c]["rate_err_Hz"][i] / 1e3 for c in r]})
    ser(f"dark_rate_{st}", f"{st}: tasa fuera de tiempo por canal (run 8, por canal, no sumada)", "canal del digitizador", "tasa", "", "kHz", "scatter",
        "measured", "run 8; charact_extract.C (pedestal + control + post, 1.456 µs por evento)", data,
        "hod2018: umbral T·G (G de c1; ch3 y ch22 con la mediana); hod2019 en mV")
    figs.append({"file": f"dark_rate_{st}.png", "title": f"{st}: tasa fuera de tiempo por canal, run 8", "xlabel": "canal", "ylabel": "kHz", "kind": "points",
                 "logy": True, "series": data})
ser("pedestal", "Pedestal por canal (run 8): nivel medio y RMS recortado", "canal del digitizador", "mV", "", "mV", "scatter", "measured",
    "run 8; charact_extract.C (media recortada ±2 mV de [0,40); RMS de las muestras a ±2 mV de la mediana)",
    [{"name": "nivel", "x": list(range(64)), "y": [CT["run008"][c]["ped_mV"] for c in range(64)], "yerr": [CT["run008"][c]["ped_sd_mV"] for c in range(64)]},
     {"name": "RMS", "x": list(range(64)), "y": [CT["run008"][c]["ped_rms_mV"] for c in range(64)], "yerr": None}],
    "barra de error del nivel = desviación estándar evento a evento")
ser("gain_hourly_summary", "hod2018: G por hora (reglas de la fase 5, daqoffline) en el run 8: mediana, mínimo y máximo por canal",
    "canal del digitizador", "G", "", "mV/p.e.", "scatter", "measured", "/data/daq/sim/offline/charact_run008/gains.json",
    [{"name": n, "x": list(range(32)), "y": [ghour[c][k] for c in range(32)], "yerr": None} for n, k in (("mediana", "median"), ("mínimo", "min"), ("máximo", "max"))])
figs.append({"file": "gain_hourly_run8.png", "title": "hod2018: G por hora (fase 5), run 8", "xlabel": "canal", "ylabel": "mV/p.e.", "kind": "points",
             "series": series[-1]["data"]})
# series de G por hora (por canal, 17 puntos)
xs = sorted({g["hour_start_utc"] for g in gains["run008"]})
ser("gain_hourly_series", "hod2018: G por hora y canal (run 8; sin valor si la hora no es ok/changed)", "hora de Chile (inicio)", "G", "", "mV/p.e.",
    "lines", "measured", "/data/daq/sim/offline/charact_run008/gains.json",
    [{"name": f"ch{c}", "x": [chile(calendar.timegm(time.strptime(h[:19], "%Y-%m-%dT%H:%M:%S")))[5:] for h in xs],
      "y": [next((g["value"] for g in gains["run008"] if g["ch"] == c and g["hour_start_utc"] == h and g["status"] in ("ok", "changed")), float("nan")) for h in xs],
      "yerr": None} for c in range(32)])
# espectros en tiempo por canal (en results.json: amplitud en tiempo y de control, integral)
for st, r, ax, ai, unit, iunit in (("hod2018", range(32), AMP18, INT18, "p.e.", "p.e. (integral / integral de 1 p.e.)"),
                                    ("hod2019", range(32, 64), AMP19, INT19, "mV", "mV·ns")):
    lo, hi, nb = ax
    xa = [lo + (hi - lo) / nb * (i + 0.5) for i in range(nb)]
    ser(f"amp_intime_{st}", f"{st}: amplitud en tiempo [49,58) por canal, run 8 (y ventana de control [40,49), mismo ancho)", "amplitud", "eventos",
        unit, "", "hist", "measured", "run 8; charact_extract.C", [{"name": f"ch{c} en tiempo", "x": xa, "y": spec_amp[c]["in"], "yerr": None} for c in r]
        + [{"name": f"ch{c} control", "x": xa, "y": spec_amp[c]["ctl"], "yerr": None} for c in r])
    lo, hi, nb = ai
    xi = [lo + (hi - lo) / nb * (i + 0.5) for i in range(nb)]
    ser(f"int_intime_{st}", f"{st}: integral en tiempo por canal, run 8", "integral", "eventos", iunit, "", "hist",
        "estimated" if st == "hod2018" else "measured", "run 8; charact_extract.C",
        [{"name": f"ch{c}", "x": xi, "y": spec_int[c], "yerr": None} for c in r],
        "hod2018: integral / integral media de los pulsos fuera de tiempo de 1 p.e. (0.5-1.5 G, 9 muestras): estimado" if st == "hod2018" else "")
    for part in (0, 1):
        chs = list(r)[16 * part:16 * part + 16]
        figs.append({"file": f"amp_intime_{st}_{part}.png", "title": f"{st}: amplitud en tiempo (negro) y control (rojo), run 8", "xlabel": unit,
                     "ylabel": "eventos", "kind": "grid", "logy": True,
                     "pads": [{"title": f"ch{c} {CH[c]['plane']}{CH[c]['strip']} [{flags[c]}]", "x": xa,
                               "ys": [spec_amp[c]["in"], spec_amp[c]["ctl"]]} for c in chs]})
        figs.append({"file": f"int_intime_{st}_{part}.png", "title": f"{st}: integral en tiempo, run 8", "xlabel": iunit, "ylabel": "eventos", "kind": "grid",
                     "logy": True, "pads": [{"title": f"ch{c} {CH[c]['plane']}{CH[c]['strip']}", "x": xi, "ys": [spec_int[c]]} for c in chs]})
prov.append({"claim": "pedestal, tasas fuera de tiempo, espectros en tiempo y estado por canal (64 canales)", "source_file": RAW["run008"] + "_seg*.root",
             "commit": git_head, "code": "daq/tools/characterization/charact_extract.C, charact_analyse.py", "inputs": "árbol events; calibración c1; mapeo r4",
             "method": "measured: buscador de pulsos del DAQ (pulse_finder.h, separación 4) fuera de tiempo; amplitud = línea base − mínimo",
             "plots": ["dark_rate_hod2018.png", "dark_rate_hod2019.png", "amp_intime_*.png", "int_intime_*.png"], "status": "measured",
             "assumptions": ["fuera de tiempo = pedestal [0,40) + control [40,49) + post [102,144); la cola [58,102) no se usa",
                             "estado: dead = tasa sobre el umbral menor (0.5 p.e. o 1 mV) < 1 % de la mediana de la estación; noisy = tasa > 3 × la mediana de la estación "
                             "(5 p.e. o 3 mV); no_hv_dependence = prueba de grupos HV (mapeo r4)",
                             f"canales sin G en c1: {G_FALLBACK} (se usa la mediana)"]})

# ---------------------------------------------------------------- 3. eficiencia de hod2018
X18, Y18 = PLANES[("hod2018", "X")], PLANES[("hod2018", "Y")]


def eff(e, T, excl=()):
    p, cp = pe(e)
    xs = [c for c in X18 if c not in excl]
    ys = [c for c in Y18 if c not in excl]
    X = (p[:, xs] > T).any(1); Y = (p[:, ys] > T).any(1)
    Xc = (cp[:, xs] > T).any(1); Yc = (cp[:, ys] > T).any(1)
    N = len(e)
    nxy = int((X & Y).sum())
    a1, a2, a3 = int((X & Yc).sum()), int((Xc & Y).sum()), int((Xc & Yc).sum())
    acc_ctl = a1 + a2 - a3
    m1, m2 = int((X[:-1] & Y[1:]).sum()), int((Y[:-1] & X[1:]).sum())
    acc_mix = 0.5 * (m1 + m2) * N / (N - 1)
    out = {"T": T, "N": N, "N_XY": nxy, "acc_ctl": acc_ctl, "acc_mix": acc_mix, "N_X": int(X.sum()), "N_Y": int(Y.sum()),
           "N_Xc": int(Xc.sum()), "N_Yc": int(Yc.sum())}
    for tag, acc, vacc in (("mix", acc_mix, (m1 + m2) / 4), ("ctl", acc_ctl, a1 + a2 + a3)):
        f = (nxy - acc) / N
        ef = math.sqrt(nxy + vacc) / N
        out[f"frac_{tag}"] = f
        out[f"frac_{tag}_err"] = ef
        out[f"eps_geo_{tag}"] = f / ACC_GEOM
        out[f"eps_geo_{tag}_err"] = ef / ACC_GEOM
        # condicionadas: impactos verdaderos del plano de referencia = sus impactos menos los de ruido (ventana de control)
        for a, b, nb_, nbc in (("X|Y", "Y", out["N_Y"], out["N_Yc"]), ("Y|X", "X", out["N_X"], out["N_Xc"])):
            den = nb_ - nbc
            num = nxy - acc
            v = num / den if den > 0 else float("nan")
            out[f"cond_{a}_{tag}"] = v
            out[f"cond_{a}_{tag}_err"] = math.sqrt(max(v * (1 - v), 0) / den + (nxy + vacc) / den ** 2) if den > 0 else float("nan")
    out["frac_raw"] = nxy / N
    return out


TS = list(range(3, 13))
EF = {k: {tag: [eff(ev[k], T, ex) for T in TS] for tag, ex in (("all", ()), ("sin_2_9_15", tuple(NOISY18)))} for k in ("run008", "night")}
tables["eff"] = EF
for k in ("run008", "night"):
    for tag in ("all", "sin_2_9_15"):
        lab = f"{NAME[k]}, {'todos los canales' if tag == 'all' else 'sin ch2, ch9, ch15'}"
        ser(f"eff_xy_{k}_{tag}", f"hod2018: fracción X·Y con accidentales restados por trigger contra umbral ({lab})", "umbral", "fracción", "p.e.", "",
            "points", "measured", RAW[k] + "_seg*.root",
            [{"name": "eventos mezclados", "x": TS, "y": [d["frac_mix"] for d in EF[k][tag]], "yerr": [d["frac_mix_err"] for d in EF[k][tag]]},
             {"name": "ventana fuera de tiempo", "x": TS, "y": [d["frac_ctl"] for d in EF[k][tag]], "yerr": [d["frac_ctl_err"] for d in EF[k][tag]]},
             {"name": "sin restar", "x": TS, "y": [d["frac_raw"] for d in EF[k][tag]], "yerr": None}])
        ser(f"eff_cond_{k}_{tag}", f"hod2018: eficiencias condicionadas ε(X|Y) y ε(Y|X) contra umbral ({lab})", "umbral", "eficiencia", "p.e.", "", "points",
            "estimated", RAW[k] + "_seg*.root",
            [{"name": f"ε({a}) {tg}", "x": TS, "y": [d[f"cond_{a}_{tg}"] for d in EF[k][tag]], "yerr": [d[f"cond_{a}_{tg}_err"] for d in EF[k][tag]]}
             for a in ("X|Y", "Y|X") for tg in ("mix",)],
            "sin corrección geométrica: un muón que cruza una tira Y puede no cruzar el área de X (supuesto no corregido)")
figs.append({"file": "eff_xy_vs_threshold.png", "title": "hod2018: fracción X·Y con accidentales restados (mezclados)", "xlabel": "umbral (p.e.)",
             "ylabel": "fracción por trigger", "kind": "points",
             "series": [{"name": f"{NAME[k]} {tag}", "x": TS, "y": [d["frac_mix"] for d in EF[k][tag]], "yerr": [d["frac_mix_err"] for d in EF[k][tag]]}
                        for k in ("run008", "night") for tag in ("all", "sin_2_9_15")]
             + [{"name": "aceptancia MC 0.0360", "x": [3, 12], "y": [ACC_GEOM, ACC_GEOM], "yerr": None}]})
figs.append({"file": "eff_cond_vs_threshold.png", "title": "hod2018: #varepsilon(X|Y) y #varepsilon(Y|X), run 8", "xlabel": "umbral (p.e.)", "ylabel": "eficiencia",
             "kind": "points", "series": [{"name": f"{a} {tag}", "x": TS, "y": [d[f"cond_{a}_mix"] for d in EF["run008"][tag]],
                                           "yerr": [d[f"cond_{a}_mix_err"] for d in EF["run008"][tag]]} for a in ("X|Y", "Y|X") for tag in ("all", "sin_2_9_15")]})
d5 = EF["run008"]["all"][2]
res("hod2018_frac_xy_5pe_run8", "hod2018: fracción X·Y por trigger con accidentales restados (mezclados), 5 p.e., run 8", d5["frac_mix"],
    {"stat": r6(d5["frac_mix_err"]), "syst": r6(abs(d5["frac_mix"] - d5["frac_ctl"]))}, "", "measured", RAW["run008"] + "_seg*.root",
    "syst = diferencia entre los dos estimadores de accidentales (mezclados y fuera de tiempo)")
res("hod2018_eps_geo_5pe_run8", "hod2018: fracción X·Y / aceptancia geométrica del MC, 5 p.e., run 8", d5["eps_geo_mix"],
    {"stat": r6(d5["eps_geo_mix_err"]), "syst": r6(abs(d5["eps_geo_mix"] - d5["eps_geo_ctl"]))}, "", "estimated", ACC_GEOM_SRC,
    "aceptancia 0.0360 del MC del handoff, sin error propio publicado; mayo: 0.88 ± 0.04 ± 0.11")
prov.append({"claim": "eficiencia X·Y de hod2018 y eficiencias condicionadas contra umbral", "source_file": RAW["run008"] + "_seg*.root; " + RAW["night"] + "_seg*.root",
             "commit": git_head, "code": "charact_analyse.py eff()", "inputs": "amplitud en tiempo [49,58) y de control [40,49) por canal; G de c1",
             "method": "measured: X = alguna tira X sobre T en tiempo, Y igual; accidentales: (a) X del evento i con Y del i+1 y al revés (promedio), "
                       "(b) X en tiempo·Y en control + X en control·Y en tiempo − ambos en control",
             "plots": ["eff_xy_vs_threshold.png", "eff_cond_vs_threshold.png"], "status": "measured (fracción); estimated (÷ aceptancia del MC, condicionadas)",
             "assumptions": [ACC_GEOM_SRC, "ε(X|Y) = (N_XY − acc) / (N_Y − N_Y en control): sin corrección por geometría (muones que cruzan Y fuera del área de X)",
                             "ch3 (X, coordenada 12) muerto en todo el análisis; ch2, ch9, ch15 ruidosos: se dan los dos casos"]})

# ---------------------------------------------------------------- 4. mapas, impactos por tira, vecinas, multiplicidad
def planes_hits(e, st, thr, window="s"):
    xs, ys = PLANES[(st, "X")], PLANES[(st, "Y")]
    a = e[window]
    if st == "hod2018":
        Xa = a[:, xs] / G[xs]; Ya = a[:, ys] / G[ys]
    else:
        Xa = a[:, xs]; Ya = a[:, ys]
    return Xa > thr, Ya > thr


def xy_map(X, Y):
    nx, ny = X.sum(1), Y.sum(1)
    w = np.where((nx > 0) & (ny > 0), 1.0 / np.maximum(nx * ny, 1), 0.0)
    m_in = (X * w[:, None]).T.astype(float) @ Y.astype(float)
    nx2, ny2 = nx[:-1], ny[1:]
    w2 = np.where((nx2 > 0) & (ny2 > 0), 1.0 / np.maximum(nx2 * ny2, 1), 0.0)
    m_mix = (X[:-1] * w2[:, None]).T.astype(float) @ Y[1:].astype(float) * len(X) / (len(X) - 1)
    return m_in, m_mix


MAPS = {}
MAPTHR = {"hod2018": [5, 10], "hod2019": [3, 5]}
for st in ("hod2018", "hod2019"):
    u = "p.e." if st == "hod2018" else "mV"
    for thr in MAPTHR[st]:
        X, Y = planes_hits(e8, st, thr)
        Xc, Yc = planes_hits(e8, st, thr, "c")
        m_in, m_mix = xy_map(X, Y)
        sub = m_in - m_mix
        hits = {pl: {"in": H.sum(0).tolist(), "ctl": Hc.sum(0).tolist(), "sub": (H.sum(0) - Hc.sum(0)).tolist()} for pl, H, Hc in (("X", X, Xc), ("Y", Y, Yc))}
        neigh = {}
        for pl, H in (("X", X), ("Y", Y)):
            Hf = H.astype(float)
            co = Hf.T @ Hf
            mix = Hf[:-1].T @ Hf[1:] * len(H) / (len(H) - 1)
            mix = 0.5 * (mix + mix.T)
            np.fill_diagonal(co, 0); np.fill_diagonal(mix, 0)
            adj = [(co[i, i + 1], mix[i, i + 1]) for i in range(15)]
            far = [(co[i, j], mix[i, j]) for i in range(16) for j in range(i + 3, 16)]
            neigh[pl] = {"coinc": co.tolist(), "mixed": mix.tolist(), "adjacent_obs": sum(a for a, _ in adj), "adjacent_mixed": sum(b for _, b in adj),
                         "far_obs": sum(a for a, _ in far), "far_mixed": sum(b for _, b in far)}
        mult = {}
        for pl, H, Hc in (("X", X, Xc), ("Y", Y, Yc)):
            m_in_ = np.bincount(H.sum(1), minlength=17)[:17]
            m_c = np.bincount(Hc.sum(1), minlength=17)[:17]
            mu = float(Hc.mean(0).sum())
            pois = [len(H) * math.exp(-mu) * mu ** k / math.factorial(k) for k in range(17)]
            mult[pl] = {"in": m_in_.tolist(), "ctl": m_c.tolist(), "poisson": pois, "mu_noise": mu, "frac_ge2_in": float((H.sum(1) >= 2).mean()),
                        "frac_ge2_poisson": 1 - math.exp(-mu) * (1 + mu)}
        MAPS[(st, thr)] = {"in": m_in.tolist(), "mixed": m_mix.tolist(), "sub": sub.tolist(), "hits": hits, "neigh": neigh, "mult": mult,
                           "sum_in": float(m_in.sum()), "sum_mixed": float(m_mix.sum())}
        labx = [f"{c}" for c in PLANES[(st, "X")]]
        laby = [f"{c}" for c in PLANES[(st, "Y")]]
        figs.append({"file": f"map_xy_{st}_{thr}.png", "title": f"{st} X #times Y, run 8, > {thr} {u}: en tiempo − mezclados (peso 1/(n_{{X}}n_{{Y}}))",
                     "kind": "map", "xlabel": "X: coordenada del mapeo (canal)", "ylabel": "Y: coordenada (canal)", "z": sub.tolist(), "xlabels": labx, "ylabels": laby})
        figs.append({"file": f"hits_{st}_{thr}.png", "title": f"{st}: impactos por tira en tiempo − control, run 8, > {thr} {u}", "xlabel": "coordenada",
                     "ylabel": "eventos", "kind": "points",
                     "series": [{"name": pl, "x": list(range(16)), "y": hits[pl]["sub"], "yerr": [math.sqrt(a + b) for a, b in zip(hits[pl]["in"], hits[pl]["ctl"])]}
                                for pl in ("X", "Y")]})
        for pl in ("X", "Y"):
            co, mix = np.array(neigh[pl]["coinc"]), np.array(neigh[pl]["mixed"])
            figs.append({"file": f"neigh_{st}_{pl}_{thr}.png", "title": f"{st} {pl}: coincidencias entre tiras − mezcladas, run 8, > {thr} {u}", "kind": "map",
                         "xlabel": "coordenada (canal)", "ylabel": "coordenada (canal)", "z": (co - mix).tolist(),
                         "xlabels": labx if pl == "X" else laby, "ylabels": labx if pl == "X" else laby})
        figs.append({"file": f"mult_{st}_{thr}.png", "title": f"{st}: multiplicidad por plano en tiempo, run 8, > {thr} {u}, y Poisson del ruido",
                     "xlabel": "tiras sobre el umbral", "ylabel": "eventos", "kind": "points", "logy": True,
                     "series": [{"name": f"{pl} {w}", "x": list(range(17)), "y": mult[pl][k], "yerr": None}
                                for pl in ("X", "Y") for w, k in (("en tiempo", "in"), ("control", "ctl"), ("Poisson", "poisson"))]})
        ser(f"hits_strip_{st}_{thr}", f"{st}: impactos por tira (en tiempo − control), run 8, > {thr} {u}", "coordenada del mapeo", "eventos", "", "", "scatter",
            "measured", RAW["run008"] + "_seg*.root",
            [{"name": pl, "x": list(range(16)), "y": hits[pl]["sub"], "yerr": [math.sqrt(a + b) for a, b in zip(hits[pl]["in"], hits[pl]["ctl"])]} for pl in ("X", "Y")],
            "canales por coordenada: X " + str(PLANES[(st, "X")]) + "; Y " + str(PLANES[(st, "Y")]))
        ser(f"mult_{st}_{thr}", f"{st}: multiplicidad por plano, run 8, > {thr} {u}", "tiras sobre el umbral", "eventos", "", "", "hist", "measured",
            RAW["run008"] + "_seg*.root", [{"name": f"{pl} {w}", "x": list(range(17)), "y": mult[pl][k], "yerr": None}
                                            for pl in ("X", "Y") for w, k in (("en tiempo", "in"), ("control", "ctl"), ("Poisson del ruido", "poisson"))],
            "Poisson: media = suma de las probabilidades de impacto en la ventana de control (mismo ancho); supone tiras independientes")
tables["maps"] = {f"{st}_{thr}": {k: v for k, v in d.items() if k in ("sum_in", "sum_mixed", "mult", "hits")} | {"neigh": {pl: {k: v for k, v in d["neigh"][pl].items()
                                                                                                                                   if k not in ("coinc", "mixed")} for pl in ("X", "Y")}}
                  for (st, thr), d in MAPS.items()}
mapres = {f"{st}_{thr}": {"sub": [[r6(v) for v in row] for row in d["sub"]], "mixed": [[r6(v) for v in row] for row in d["mixed"]],
                          "neigh_X_minus_mixed": [[r6(a - b) for a, b in zip(r1, r2)] for r1, r2 in zip(d["neigh"]["X"]["coinc"], d["neigh"]["X"]["mixed"])],
                          "neigh_Y_minus_mixed": [[r6(a - b) for a, b in zip(r1, r2)] for r1, r2 in zip(d["neigh"]["Y"]["coinc"], d["neigh"]["Y"]["mixed"])]}
          for (st, thr), d in MAPS.items() if thr in (5, 3)}
prov.append({"claim": "mapas X × Y, impactos por tira, matrices de vecinas y multiplicidad por plano", "source_file": RAW["run008"] + "_seg*.root", "commit": git_head,
             "code": "charact_analyse.py", "inputs": "amplitud en tiempo y de control; mapeo r4 (coordenada strip)",
             "method": "measured: mapa con peso 1/(nX·nY) por evento; accidentales con X del evento i e Y del i+1; impactos en tiempo − ventana de control",
             "plots": ["map_xy_*.png", "hits_*.png", "neigh_*.png", "mult_*.png"], "status": "measured",
             "assumptions": [f"orden de tiras = coordenada del mapeo r4; confirmación por plano: {STRIP_CONF}",
                             "hod2018 en p.e. (G de c1), hod2019 en mV"]})

# ---------------------------------------------------------------- 5. hod2019
tag = (pe8[:, X18] > 5).any(1) & (pe8[:, Y18] > 5).any(1)  # eventos con X·Y de hod2018 en tiempo (enriquecidos en muones)
# Eventos coherentes: ≥ 8 de los 32 canales de hod2019 sobre 3 mV en la misma ventana (todos los canales a la vez:
# no es un muón en una tira). Se cuentan aparte y se excluyen del exceso por canal.
coh_in = (e8["s"][:, 32:64] > 3).sum(1) >= 8
coh_ctl = (e8["c"][:, 32:64] > 3).sum(1) >= 8
coh = coh_in | coh_ctl
coh_info = {"n_in": int(coh_in.sum()), "n_ctl": int(coh_ctl.sum()), "n_in_with_h18_xy5": int((coh_in & tag).sum()),
            "n_in_h18_any_gt50pe": int((coh_in & (pe8 > 50).any(1)).sum()),
            "times_chile": [chile(x) for x in e8["t"][coh_in][:40]],
            "h18_max_pe_median": float(np.median(pe8[coh_in].max(1))) if coh_in.any() else float("nan"),
            "h19_min_ch_over3": int((e8["s"][coh_in][:, 32:64] > 3).sum(1).min()) if coh_in.any() else 0}
h19 = {}
for c in range(32, 64):
    s, cc, k = e8["s"][:, c], e8["c"][:, c], e8["k"][:, c]
    d = {}
    for thr in (3, 5, 10):
        a, b = int((s[~coh] > thr).sum()), int((cc[~coh] > thr).sum())
        at, bt = int((s[tag & ~coh] > thr).sum()), int((cc[tag & ~coh] > thr).sum())
        d[f"{thr}mV"] = {"in": a, "ctl": b, "excess": a - b, "z": (a - b) / math.sqrt(a + b) if a + b else 0.0,
                         "in_tag": at, "ctl_tag": bt, "excess_tag": at - bt, "z_tag": (at - bt) / math.sqrt(at + bt) if at + bt else 0.0}
    big = (s > 3) & ~coh
    d["frac_min_in_signal"] = float(((k[big] >= SIG[0]) & (k[big] < SIG[1])).mean()) if big.sum() else float("nan")
    d["n_big"] = int(big.sum())
    h19[c] = d
tables["hod2019_intime"] = {"tag_events": int(tag.sum()), "coherent": coh_info, "per_channel": h19}
ser("hod2019_intime_excess", "hod2019: exceso en tiempo sobre la ventana de control (> 3 mV) por canal, sin los eventos coherentes (≥ 8 canales), todos y con X·Y de hod2018 > 5 p.e.",
    "canal del digitizador", "exceso", "", "eventos", "scatter", "measured", RAW["run008"] + "_seg*.root",
    [{"name": "todos", "x": list(range(32, 64)), "y": [h19[c]["3mV"]["excess"] for c in range(32, 64)],
      "yerr": [math.sqrt(h19[c]["3mV"]["in"] + h19[c]["3mV"]["ctl"]) for c in range(32, 64)]},
     {"name": "con X·Y de hod2018", "x": list(range(32, 64)), "y": [h19[c]["3mV"]["excess_tag"] for c in range(32, 64)],
      "yerr": [math.sqrt(h19[c]["3mV"]["in_tag"] + h19[c]["3mV"]["ctl_tag"]) for c in range(32, 64)]}])
figs.append({"file": "hod2019_intime_excess.png", "title": "hod2019: en tiempo − control (> 3 mV) por canal, run 8", "xlabel": "canal", "ylabel": "eventos",
             "kind": "points", "series": series[-1]["data"]})


# tasa oscura por hora contra temperatura de la placa y hora del día (todos los runs)
def temp_by_hour(k):
    b = {}
    for x in SC[k]:
        b.setdefault(int(x["t"] // 3600) * 3600, []).append(x["T"])
    return {h: statistics.mean(v) for h, v in b.items() if v}


dn = []
for k in RAW:
    th = temp_by_hour(k)
    hrs = sorted({r["hour_start_utc_s"] for r in HR[k]})
    for h in hrs:
        rows = [r for r in HR[k] if r["hour_start_utc_s"] == h]
        n = rows[0]["events"]
        if n < 200:
            continue
        r19 = [r["cnt0"] / r["exposure_s"] for r in rows if r["ch"] >= 32]
        r19_3 = [r["cnt1"] / r["exposure_s"] for r in rows if r["ch"] >= 32]
        s19 = sum(r["cnt0"] for r in rows if r["ch"] >= 32)
        ex = rows[0]["exposure_s"]
        r18 = [r["cnt2"] / r["exposure_s"] for r in rows if r["ch"] < 32 and r["ch"] != 3]
        dn.append({"run": k, "hour_utc": utc(h), "hour_chile": chile(h), "chile_h": int(((h + CHILE) % 86400) // 3600), "events": n,
                   "board_T": th.get(int(h), float("nan")), "h19_1mV_median_Hz": statistics.median(r19), "h19_3mV_median_Hz": statistics.median(r19_3),
                   "h19_1mV_sum_Hz": s19 / ex, "h19_1mV_sum_err_Hz": math.sqrt(s19) / ex, "h18_3pe_median_Hz": statistics.median(r18),
                   "h19_1mV_per_ch": [r["cnt0"] / r["exposure_s"] for r in sorted(rows, key=lambda r: r["ch"]) if r["ch"] >= 32]})
tables["dark_vs_time"] = [{k: v for k, v in d.items() if k != "h19_1mV_per_ch"} for d in dn]
ok = [d for d in dn if d["board_T"] == d["board_T"]]


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


corr = {"n_hours": len(ok), "r_T_h19_sum_log": pearson([d["board_T"] for d in ok], [math.log(max(d["h19_1mV_sum_Hz"], 1e-3)) for d in ok]),
        "r_T_h18_median_log": pearson([d["board_T"] for d in ok], [math.log(d["h18_3pe_median_Hz"]) for d in ok]),
        "T_range": [min(d["board_T"] for d in ok), max(d["board_T"] for d in ok)] if ok else []}
if len(ok) >= 3:
    A = np.vstack([[d["board_T"] for d in ok], np.ones(len(ok))]).T
    yv = np.log([max(d["h19_1mV_sum_Hz"], 1e-3) for d in ok])
    coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
    corr["h19_sum_dlog_dT_per_degC"] = float(coef[0])
tables["dark_vs_temp"] = corr
ser("h19_dark_vs_time", "hod2019: tasa oscura > 1 mV por hora (suma de los 32 canales y mediana por canal) y temperatura de la placa", "hora de Chile",
    "tasa", "", "kHz", "points", "measured", "; ".join(RAW[k] + "_seg*.root" for k in RAW) + "; hv_full_*.csv",
    [{"name": f"{NAME[k]} suma", "x": [d["hour_chile"][5:] for d in dn if d["run"] == k], "y": [d["h19_1mV_sum_Hz"] / 1e3 for d in dn if d["run"] == k],
      "yerr": [d["h19_1mV_sum_err_Hz"] / 1e3 for d in dn if d["run"] == k]} for k in RAW]
    + [{"name": f"{NAME[k]} mediana por canal", "x": [d["hour_chile"][5:] for d in dn if d["run"] == k], "y": [d["h19_1mV_median_Hz"] / 1e3 for d in dn if d["run"] == k],
        "yerr": None} for k in RAW]
    + [{"name": f"{NAME[k]} placa (°C)", "x": [d["hour_chile"][5:] for d in dn if d["run"] == k], "y": [d["board_T"] for d in dn if d["run"] == k], "yerr": None}
       for k in RAW],
    "horas con ≥ 200 eventos; las sumas sobre canales anteriores (p. ej. 9942.3 kHz) eran sumas, aquí también se da la mediana por canal")
ser("h19_dark_per_channel_hour", "hod2019: tasa oscura > 1 mV por canal y hora (run 8)", "hora de Chile", "tasa", "", "Hz", "lines", "measured",
    RAW["run008"] + "_seg*.root", [{"name": f"ch{c}", "x": [d["hour_chile"][5:] for d in dn if d["run"] == "run008"],
                                    "y": [d["h19_1mV_per_ch"][c - 32] for d in dn if d["run"] == "run008"], "yerr": None} for c in range(32, 64)])
figs.append({"file": "h19_dark_vs_temp.png", "title": "hod2019: tasa oscura > 1 mV (suma) contra temperatura de la placa, por hora", "xlabel": "placa A7030DP (#circC)",
             "ylabel": "kHz", "kind": "points", "logy": True,
             "series": [{"name": NAME[k], "x": [d["board_T"] for d in ok if d["run"] == k], "y": [d["h19_1mV_sum_Hz"] / 1e3 for d in ok if d["run"] == k],
                         "yerr": [d["h19_1mV_sum_err_Hz"] / 1e3 for d in ok if d["run"] == k]} for k in RAW]})
figs.append({"file": "h19_dark_vs_hour.png", "title": "hod2019: tasa oscura > 1 mV (suma) por hora de Chile", "xlabel": "hora de Chile (h)", "ylabel": "kHz",
             "kind": "points", "logy": True,
             "series": [{"name": NAME[k], "x": [d["chile_h"] for d in dn if d["run"] == k], "y": [d["h19_1mV_sum_Hz"] / 1e3 for d in dn if d["run"] == k],
                         "yerr": [d["h19_1mV_sum_err_Hz"] / 1e3 for d in dn if d["run"] == k]} for k in RAW]})

# ---------------------------------------------------------------- 6. día / noche
gw = {}
for k in RAW:
    p = f"{W}/gainwin_{k}.csv"
    for r in csv.DictReader(open(p)):
        gw.setdefault(r["label"], {})[int(r["ch"])] = (float(r["G_mV"]), r["status"], int(r["events"]))


def window(k, h0, h1):
    """Eventos del run k entre las horas de Chile h0 y h1 (todas las fechas del run)."""
    e = ev[k]
    ch_h = ((e["t"] + CHILE) % 86400) / 3600
    return e[(ch_h >= h0) & (ch_h < h1)]


def window_stats(k, h0, h1, label):
    e = window(k, h0, h1)
    if len(e) < 50:
        return {"label": label, "events": int(len(e))}
    t = e["t"]
    g = gaps(t)
    # tiempo cubierto: por día, primer a último evento menos huecos
    days = np.floor((t + CHILE) / 86400)
    cov = 0.0
    for d in np.unique(days):
        tt = t[days == d]
        cov += tt[-1] - tt[0] - sum(x[2] for x in gaps(tt))
    hrs = [r for r in HR[k] if h0 <= ((r["hour_start_utc_s"] + CHILE) % 86400) / 3600 < h1]
    dark = {}
    for st, chs, i in (("hod2018_3pe", [c for c in range(32) if c != 3], 2), ("hod2019_1mV", list(range(32, 64)), 0)):
        cnt = sum(r[f"cnt{i}"] for r in hrs if r["ch"] in chs)
        ex = sum(r["exposure_s"] for r in hrs if r["ch"] in chs) / len(chs)
        dark[st] = {"sum_Hz": cnt / ex if ex else float("nan"), "err_Hz": math.sqrt(cnt) / ex if ex else float("nan")}
    e5 = eff(e, 5)
    T = [x["T"] for x in SC[k] if h0 <= ((x["t"] + CHILE) % 86400) / 3600 < h1 and x["T"] == x["T"]]
    return {"label": label, "events": int(len(e)), "covered_s": cov, "rate_Hz": len(e) / cov, "rate_err_Hz": math.sqrt(len(e)) / cov, "dark": dark,
            "frac_xy_5pe_mix": e5["frac_mix"], "frac_xy_5pe_mix_err": e5["frac_mix_err"], "board_T_mean": statistics.mean(T) if T else float("nan"),
            "first_chile": chile(t[0]), "last_chile": chile(t[-1])}


DN = {"noche_run8": window_stats("run008", 0, 6, "run 8, 00-06 h Chile (2026-10-23)"),
      "dia_run8": window_stats("run008", 12, 17, "run 8, 12-17 h Chile"),
      "dia_run7": window_stats("run007", 12, 17, "run 7, 12-17 h Chile (2026-10-22; datos solo desde 15:45)"),
      "dia_run4": window_stats("run004", 12, 17, "run 4, 12-17 h Chile (2026-10-21)"),
      "noche_night": window_stats("night", 0, 6, "night_2026_10_19, 00-06 h Chile (2026-10-20)")}


def gain_ratio(lab_d, lab_n):
    a, b = gw.get(lab_d, {}), gw.get(lab_n, {})
    r = [a[c][0] / b[c][0] for c in range(32) if c in a and c in b and a[c][1] in ("ok", "changed") and b[c][1] in ("ok", "changed")]
    if len(r) < 2:
        return {"n": len(r)}
    return {"n": len(r), "median": statistics.median(r), "mean": statistics.mean(r), "err_mean": statistics.stdev(r) / math.sqrt(len(r))}


def ratio(a, ea, b, eb):
    if not (a == a and b == b) or b == 0:
        return float("nan"), float("nan")
    r = a / b
    return r, abs(r) * math.sqrt((ea / a) ** 2 + (eb / b) ** 2) if a else float("nan")


dnr = {"windows": DN, "gain_ratio_run7day_over_run8night": gain_ratio("dia_run7_12-17", "noche_run8_00-06"),
       "gain_ratio_night191020_over_run8night": gain_ratio("noche_191020_00-06", "noche_run8_00-06")}
d, n = DN["dia_run7"], DN["noche_run8"]
if "rate_Hz" in d:
    dnr["trigger_rate_ratio"] = ratio(d["rate_Hz"], d["rate_err_Hz"], n["rate_Hz"], n["rate_err_Hz"])
    dnr["dark18_ratio"] = ratio(d["dark"]["hod2018_3pe"]["sum_Hz"], d["dark"]["hod2018_3pe"]["err_Hz"], n["dark"]["hod2018_3pe"]["sum_Hz"], n["dark"]["hod2018_3pe"]["err_Hz"])
    dnr["dark19_ratio"] = ratio(d["dark"]["hod2019_1mV"]["sum_Hz"], d["dark"]["hod2019_1mV"]["err_Hz"], n["dark"]["hod2019_1mV"]["sum_Hz"], n["dark"]["hod2019_1mV"]["err_Hz"])
    dnr["eff_ratio"] = ratio(d["frac_xy_5pe_mix"], d["frac_xy_5pe_mix_err"], n["frac_xy_5pe_mix"], n["frac_xy_5pe_mix_err"])
tables["day_night"] = dnr
prov.append({"claim": "razones día (12-17 h Chile) / noche (00-06 h Chile)", "source_file": "; ".join(RAW[k] + "_seg*.root" for k in ("run008", "run007", "run004")),
             "commit": git_head, "code": "charact_analyse.py window_stats(); charact_gain.C (G con gain_fit.h sobre el espectro de la ventana)",
             "method": "measured; G por ventana con las reglas de la fase 5 sobre el espectro fuera de tiempo acumulado en la ventana",
             "plots": [], "status": "measured (con muy poca estadística de día)",
             "assumptions": ["el run 8 no tiene horas de 12-17 h Chile (17:20 del 22 a 10:36 del 23): el día sale del run 7 (2026-10-22, 15:45-16:16 h, 1274 eventos)",
                             "el run 4 (2026-10-21) tiene 50 eventos de 12-17 h (sin triggers 17:2x-20:00Z): no sirve",
                             "tasas oscuras: suma sobre canales (hod2018 sin ch3) por canal-tiempo expuesto"]})

# ---------------------------------------------------------------- resultados principales (results)
o8 = over["run008"]
res("run8_trigger_rate", "Run 8: tasa media de trigger externo", o8["rate_Hz"], o8["rate_err_Hz"], "Hz", "measured", RAW["run008"] + "_seg*.root",
    f"{o8['events']} eventos en {o8['live_s']:.0f} s; huecos > 20 s: {len(o8['gaps_gt20s'])}")
res("run8_board_T_range", "Run 8: temperatura de la placa A7030DP (mín.-máx. de medias de 10 min)", tables["hv_run8"]["board_T_max"] - tables["hv_run8"]["board_T_min"],
    None, "°C", "measured", "hv_full_2026102[23].csv*", f"{tables['hv_run8']['board_T_min']:.1f}-{tables['hv_run8']['board_T_max']:.1f} °C; resolución de lectura 1 °C")
for ch in range(8):
    h = tables["hv_run8"]["per_channel"][ch]
    res(f"run8_imon_ch{ch:02d}", f"Run 8: IMon mediana de CH{ch:02d} (grupo {ch + 1})", h["IMon_uA_median"], {"min": r6(h["IMon_uA_min"]), "max": r6(h["IMon_uA_max"])},
        "µA", "measured", "hv_full_2026102[23].csv*", f"VMon {h['VMon_V_min']:.2f}-{h['VMon_V_max']:.2f} V; V0Set {h['V0Set_V']}; Pw mínimo {h['Pw_min']}")
res("hod2018_dead_noisy", "hod2018: canales muertos y ruidosos (run 8)", None, None, "", "measured", RAW["run008"] + "_seg*.root",
    "muertos: " + str([c for c in range(32) if flags[c] == "dead"]) + "; ruidosos (> 3 × mediana a 5 p.e.): " + str([c for c in range(32) if flags[c] == "noisy"]))
res("hod2019_flags", "hod2019: estados por canal (run 8)", None, None, "", "measured", RAW["run008"] + "_seg*.root; mapeo r4",
    "; ".join(f"{f}: {[c for c in range(32, 64) if flags[c] == f]}" for f in sorted(set(flags[c] for c in range(32, 64)))))
d5n = EF["night"]["all"][2]
res("hod2018_frac_xy_5pe_night", "hod2018: fracción X·Y por trigger con accidentales restados (mezclados), 5 p.e., noche 2026-10-19", d5n["frac_mix"],
    {"stat": r6(d5n["frac_mix_err"]), "syst": r6(abs(d5n["frac_mix"] - d5n["frac_ctl"]))}, "", "measured", RAW["night"] + "_seg*.root",
    "DAQ2-13 dio 0.0167 con otra definición de impacto (pico de pulso dentro de la ventana, selftrig_check.C)")
for a_ in ("X|Y", "Y|X"):
    res(f"hod2018_cond_{a_.replace('|', '_given_')}_5pe_run8", f"hod2018: ε({a_}) a 5 p.e., run 8 (sin corrección geométrica)", d5[f"cond_{a_}_mix"],
        d5[f"cond_{a_}_mix_err"], "", "estimated", RAW["run008"] + "_seg*.root", "accidentales por eventos mezclados")
dnn = tables["day_night"]
for k_, lab in (("trigger_rate_ratio", "tasa de trigger"), ("dark18_ratio", "tasa oscura hod2018 > 3 p.e. (suma)"), ("dark19_ratio", "tasa oscura hod2019 > 1 mV (suma)"),
                ("eff_ratio", "fracción X·Y a 5 p.e.")):
    if k_ in dnn:
        res(f"daynight_{k_}", f"Día (run 7, 2026-10-22 15:45-16:16 h) / noche (run 8, 2026-10-23 00-06 h): {lab}", dnn[k_][0], dnn[k_][1], "", "measured",
            RAW["run007"] + "; " + RAW["run008"], "el run 8 no tiene horas de día; un solo tramo de 31 min de día")
g_ = dnn["gain_ratio_run7day_over_run8night"]
if g_.get("n", 0) >= 2:
    res("daynight_gain_ratio", "Día (run 7) / noche (run 8): G de hod2018 (mediana sobre canales con G aceptada en ambas ventanas)", g_["median"], g_["err_mean"], "",
        "measured", "charact_gain.C (gain_fit.h)", f"{g_['n']} canales; error = desviación estándar de la media de las razones")


# ---------------------------------------------------------------- salida
def root_text(x):
    for a_, b_ in (("−", "-"), ("×", "#times"), ("·", "#upoint"), ("≥", "#geq"), ("ε", "#varepsilon"), ("µ", "#mu")):
        x = x.replace(a_, b_)
    return x


for f_ in figs:
    for k_ in ("title", "xlabel", "ylabel"):
        if k_ in f_:
            f_[k_] = root_text(f_[k_])
    for s_ in f_.get("series", []):
        s_["name"] = root_text(s_["name"])
out = {"generated_utc": utc(time.time()), "git_head": git_head, "provenance": prov, "results": results, "series": series, "episodes": episodes,
       "maps": mapres, "per_channel": {str(c): {k: (r6(v) if isinstance(v, float) else v) for k, v in d.items()} for c, d in per_channel.items()}}
json.dump(out, open(f"{OUT}/results.json", "w"), ensure_ascii=False, default=lambda x: r6(x) if isinstance(x, (float, np.floating)) else int(x) if isinstance(x, np.integer) else str(x))
json.dump(figs, open(f"{W}/figs.json", "w"), ensure_ascii=False, default=float)
json.dump(tables, open(f"{W}/tables.json", "w"), ensure_ascii=False, indent=1, default=lambda x: float(x) if isinstance(x, np.floating) else int(x) if isinstance(x, np.integer) else str(x))
print("results.json", os.path.getsize(f"{OUT}/results.json"), "bytes;", len(series), "series;", len(figs), "figuras")
