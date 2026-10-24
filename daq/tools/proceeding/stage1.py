#!/usr/bin/python3 -I
# stage1.py — PROC-1, etapa 1 (solo lectura): inventario de los runs de mayo, calibraciones y validación del orden de las
# tiras de mayo v8 por coincidencias entre vecinas.
#  - Inventario: work/may_inventory.jsonl (may_inventory.C) + tiempo bueno de v8 (may_common.good_time) -> tables/may_runs.tex.
#  - Calibración: G y p0 de outputs/pe_calibration_v8.txt (convención de mayo), G y p0 reajustados con la definición única
#    (outputs/characterization_2026_10_23/work/hd_may8_peaks.csv, hitdef_extract.C: posiciones de picos con el pedestal
#    común; p0 = media de pos_i − (n0 + i)·G, n0 = round(pos_0 / G)), G de c1 y G_mayo / G_c1 -> tables/calibration.tex.
#    A3 (René): el análisis de mayo usa la G reajustada en todos los canales (work/hd64_may_v8_peaks.csv, semilla de ch9 = su G
#    de c1, porque la G del archivo de mayo de ch9 es probablemente la mitad de la verdadera).
#  - Orden de tiras: impactos con la definición única (work/hd64_may8_hits.bin, 64 canales, tiempo bueno), umbral 5 p.e.
#    (G de mayo) en hod2018 y 3 mV en hod2019; por grupo de 16 canales: coincidencias en tiempo co[i,j], esperadas por
#    azar con eventos mezclados (i del evento k con j del k+1 y al revés, promedio, × N/(N−1)); exceso E = co − mix,
#    z = E / sqrt(co + mix). Arista aceptada: z > Z_MIN y cada canal está entre los dos de mayor exceso del otro. Las
#    aristas aceptadas forman cadenas (grado ≤ 2, sin ciclos): el orden medido es el de cada cadena (orientación absoluta
#    no determinada). Se compara con kCfg_2026_05 (analysis/src/hodo_mapping.h), la hoja de mayo (canal k = coordenada k)
#    y el mapeo r4 de octubre.
# Uso: stage1.py <dir de salida del paquete>
import json, math, os, re, statistics, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
Z_MIN = 5.0
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
res = {"stage": 1, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

# ---------------------------------------------------------------- 1. inventario
inv = [json.loads(l) for l in open(W + "/may_inventory.jsonl")]
e8 = xy.load_hits(W + "/hd64_may8_hits.bin", 64)
good8 = mc.good_time(e8["t"], "v8")
res["inventory"] = inv
res["v8_good_time_events"] = int(good8.sum())
res["v8_good_intervals"] = mc.good_intervals("v8")
L = ["% PROC-1 etapa 1: runs de mayo 2026 (archivos de /data/daq/raw/hodoscope_extnim_2026). Generado por daq/tools/proceeding/stage1.py, "
     f"commit {git}.",
     "% UTC: timestamp del contador + época de referencia de hod2018_good_time_step2.C (validada con la fecha de modificación de cada archivo).",
     "\\begin{tabular}{lrrrllr}", "\\hline",
     "Run & Triggers & Duración (h) & Huecos $>$20 s (s) & Inicio (UTC) & Fin (UTC) & HV \\\\", "\\hline"]
for d in inv:
    L.append(f"{d['label']} & {d['entries']} & {d['duration_s'] / 3600:.2f} & {d['gap_sum_s']:.0f} & "
             f"{time.strftime('%Y-%m-%d %H:%M', time.gmtime(d['utc_first_s']))} & {time.strftime('%Y-%m-%d %H:%M', time.gmtime(d['utc_last_s']))} & no \\\\")
L += ["\\hline", "\\end{tabular}",
      f"% v8 en tiempo bueno: {int(good8.sum())} de {len(e8)} triggers (runs/good_intervals_2026.csv, fuera EXCLUDE y HV_OFF).",
      "% HV: no hay registro de V0Set por run para v1-v10 (runs/runs_2026.csv: bias_voltage_V = not recorded). Solo la captura de "
      "GECO2020 citada en docs/handoff_analysis_2026-10-18.md:21, sin run asociado."]
os.makedirs(OUT + "/tables", exist_ok=True)
open(OUT + "/tables/may_runs.tex", "w").write("\n".join(L) + "\n")

# ---------------------------------------------------------------- 2. calibración
Gm, p0m, Gc1 = mc.g_may(), mc.p0_may(), mc.g_c1()
peaks = {}
for l in open(W + "/hd64_may_v8_peaks.csv"):  # A3: reajuste de v8 con la semilla de ch9 = su G de c1
    if l.startswith("ch,"):
        continue
    f = l.strip().split(",")
    ch, st = int(f[0]), f[3]
    pos = [float(x) for x in f[5].split()] if len(f) > 5 and f[5] else []
    G = float(f[2])
    p0 = None
    if st in ("ok", "changed") and pos and G > 0:
        n0 = round(pos[0] / G)
        p0 = statistics.mean(p - (n0 + i) * G for i, p in enumerate(pos))
    peaks[ch] = {"G_fit_mV": G if st in ("ok", "changed") else None, "status": st, "p0_fit_mV": p0, "positions_mV": pos}
cal = []
for ch in range(32):
    pk = peaks.get(ch, {})
    cal.append({"ch": ch, "G_may_mV": None if not np.isfinite(Gm[ch]) else float(Gm[ch]), "p0_may_mV": float(p0m[ch]) if np.isfinite(Gm[ch]) else None,
                "G_fit_mV": pk.get("G_fit_mV"), "p0_fit_mV": pk.get("p0_fit_mV"), "fit_status": pk.get("status"),
                "G_c1_mV": None if not np.isfinite(Gc1[ch]) else float(Gc1[ch]),
                "ratio_may_c1": float(Gm[ch] / Gc1[ch]) if np.isfinite(Gm[ch]) and np.isfinite(Gc1[ch]) else None})
res["calibration"] = cal
fits = [c for c in cal if c["G_fit_mV"] and c["G_may_mV"]]
res["calibration_summary"] = {"n_fit": len(fits), "ratio_fit_over_file_median": statistics.median(c["G_fit_mV"] / c["G_may_mV"] for c in fits),
                              "p0_fit_median_mV": statistics.median(c["p0_fit_mV"] for c in fits if c["p0_fit_mV"] is not None),
                              "ratio_may_c1_median": statistics.median(c["ratio_may_c1"] for c in cal if c["ratio_may_c1"])}
T = ["% PROC-1 etapa 1: calibración de hod2018. G_mayo, p0_mayo: outputs/pe_calibration_v8.txt (convención de mayo: pedestal local, primer pico "
     "visible numerado 1 p.e.). G_aj, p0_aj: picos de pulsos oscuros de v8 con el pedestal común de la definición única (hitdef_extract.C, "
     "gain_fit.h; semilla de ch9 = su G de c1). G_c1: calibración de octubre 2026_10-c1. Se usa G_aj en mayo (A3) y G_c1 en octubre; "
     "G_mayo y p0_mayo solo para la convención publicada (p0 + 3.5 G). Umbral de 50 mV en p.e.: 50 / G_aj.",
     "\\begin{tabular}{rrrrrrr}", "\\hline", "Canal & $G_\\mathrm{mayo}$ (mV) & $p_0$ mayo (mV) & $G_\\mathrm{aj}$ (mV) & $p_0$ aj (mV) & "
     "$G_\\mathrm{c1}$ (mV) & $G_\\mathrm{mayo}/G_\\mathrm{c1}$ & 50 mV en p.e. \\\\", "\\hline"]
fmt = lambda v, d=3: "--" if v is None else f"{v:.{d}f}"
for c in cal:
    note = "$^a$" if c["ch"] == 9 else ""
    T.append(f"{c['ch']} & {fmt(c['G_may_mV'])}{note} & {fmt(c['p0_may_mV'])} & {fmt(c['G_fit_mV'])} & {fmt(c['p0_fit_mV'])} & {fmt(c['G_c1_mV'])} & "
             f"{fmt(c['ratio_may_c1'])} & {fmt(50.0 / c['G_fit_mV'] if c['G_fit_mV'] else None, 1)} \\\\")
T += ["\\hline", "\\end{tabular}", "% $^a$ ch9: la G del archivo de mayo es probablemente la mitad de la ganancia verdadera (G_aj = 5.437 mV, G_c1 = 5.050 mV)."]
open(OUT + "/tables/calibration.tex", "w").write("\n".join(T) + "\n")

# ---------------------------------------------------------------- 3. orden de tiras
src = open(REPO + "/analysis/src/hodo_mapping.h").read()
blk = src[src.index("kCfg_2026_05 = {"):]
arrs = [[int(x) for x in a.split(",")] for a in re.findall(r"\{([0-9 ,]+)\}", blk)[:4]]
CFG05 = {"hod2018 X": arrs[0], "hod2018 Y": arrs[1], "hod2019 (ch48-63)": arrs[2], "hod2019 (ch32-47)": arrs[3]}
r4 = json.load(open(os.path.expanduser("~/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r4.json")))
R4 = {}
for c in r4["channels"]:
    R4.setdefault((c["station"], c["plane"]), []).append((c["strip"], c["dig_ch"]))
R4 = {k: [ch for _, ch in sorted(v)] for k, v in R4.items()}

e = e8[good8]
Gfill = np.where(np.isfinite(Gm), Gm, np.nanmedian(Gm))
N = len(e)


def neighbours(chans, thr):
    H = (e["aS"][:, chans] > thr).astype(float)
    co = H.T @ H
    mix = (H[:-1].T @ H[1:]) * N / (N - 1)
    mix = 0.5 * (mix + mix.T)
    np.fill_diagonal(co, 0)
    np.fill_diagonal(mix, 0)
    E = co - mix
    z = np.where(co + mix > 0, E / np.sqrt(np.maximum(co + mix, 1)), 0)
    n = len(chans)
    top2 = {i: set(sorted([j for j in range(n) if j != i and z[i, j] > Z_MIN], key=lambda j: -E[i, j])[:2]) for i in range(n)}
    edges = sorted({(min(i, j), max(i, j)) for i in range(n) for j in top2[i] if i in top2[j]})
    adj = {i: [] for i in range(n)}
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    chains, seen = [], set()
    cyc = False
    for s0 in [i for i in range(n) if len(adj[i]) == 1] + [i for i in range(n) if len(adj[i]) == 2]:
        if s0 in seen:
            continue
        path, prev, cur = [s0], None, s0
        seen.add(s0)
        while True:
            nxt = [k for k in adj[cur] if k != prev and k not in seen]
            if not nxt:
                if len(adj[cur]) == 2 and adj[s0] and len(path) > 2 and s0 in adj[cur]:
                    cyc = True
                break
            prev, cur = cur, nxt[0]
            path.append(cur)
            seen.add(cur)
        chains.append([chans[k] for k in path])
    isolated = [chans[i] for i in range(n) if not adj[i]]
    chains = [c for c in chains if len(c) > 1]
    return {"threshold_mV_or_pe": None, "events": N, "edges": [[chans[i], chans[j], float(E[i, j]), float(z[i, j])] for i, j in edges],
            "chains": chains, "isolated": isolated, "cycle": cyc, "co": co.tolist(), "mixed": mix.tolist(), "channels": chans,
            "max_degree": max(len(v) for v in adj.values())}


def adjacency_check(order, edges):
    es = {(min(a, b), max(a, b)) for a, b, _, _ in edges}
    pairs = [(min(order[k], order[k + 1]), max(order[k], order[k + 1])) for k in range(len(order) - 1)]
    return {"pairs": len(pairs), "confirmed": sum(p in es for p in pairs), "missing": [list(p) for p in pairs if p not in es],
            "extra_edges": [list(p) for p in sorted(es) if p not in set(pairs)]}


ORD = {}
for name, chans, thr in (("hod2018 X", list(range(16)), 5 * Gfill[:16]), ("hod2018 Y", list(range(16, 32)), 5 * Gfill[16:32]),
                         ("hod2019 ch32-47", list(range(32, 48)), 3.0), ("hod2019 ch48-63", list(range(48, 64)), 3.0)):
    r = neighbours(chans, thr)
    r["threshold"] = "5 p.e. (G de mayo)" if name.startswith("hod2018") else "3 mV"
    cmp = {}
    if name == "hod2018 X":
        cmp = {"kCfg_2026_05": adjacency_check(CFG05["hod2018 X"], r["edges"]), "hoja_mayo (canal k = coordenada k)": adjacency_check(chans, r["edges"]),
               "r4_octubre": adjacency_check(R4[("hod2018", "X")], r["edges"])}
    elif name == "hod2018 Y":
        cmp = {"kCfg_2026_05": adjacency_check(CFG05["hod2018 Y"], r["edges"]), "hoja_mayo (canal k = coordenada k)": adjacency_check(chans, r["edges"]),
               "r4_octubre": adjacency_check(R4[("hod2018", "Y")], r["edges"])}
    elif name == "hod2019 ch48-63":
        cmp = {"kCfg_2026_05 (hoja)": adjacency_check(CFG05["hod2019 (ch48-63)"], r["edges"])}
    else:
        cmp = {"kCfg_2026_05 (hoja)": adjacency_check(CFG05["hod2019 (ch32-47)"], r["edges"])}
    r["comparison"] = cmp
    ORD[name] = r
res["strip_order"] = {k: {kk: vv for kk, vv in v.items() if kk not in ("co", "mixed")} for k, v in ORD.items()}
res["strip_order_matrices"] = {k: {"channels": v["channels"], "co": v["co"], "mixed": v["mixed"]} for k, v in ORD.items()}
json.dump(res, open(W + "/stage1.json", "w"), indent=1, default=float)
for k, v in ORD.items():
    print(k, v["threshold"], "N", v["events"], "aristas", len(v["edges"]), "grado máx", v["max_degree"], "ciclo", v["cycle"])
    print("  cadenas:", v["chains"], "aisladas:", v["isolated"])
    for m, c in v["comparison"].items():
        print(f"  {m}: {c['confirmed']}/{c['pairs']} pares vecinos confirmados; faltan {c['missing']}; aristas extra {c['extra_edges']}")
print("calibración:", res["calibration_summary"], "v8 tiempo bueno:", res["v8_good_time_events"])
