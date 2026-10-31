#!/usr/bin/python3 -I
# vbd_analysis.py — seguimiento 2026-10-31 (9), punto 1: análisis post-barrido del voltaje de ruptura de hod2018 (solo lectura de crudos y del
# slow control; no toca hardware).
# Por toma SW (un escalón del barrido): espectro de pulsos oscuros por canal con vbd_gain.C (reglas de la fase 5, las de c1; bins de 1 LSB) y dos
# medidas de G:
#   "c1":     gain_fit.h con las opciones de c1 (solo ajustes aceptados, ok/changed); error de G = error estándar de la pendiente de las posiciones
#             de los picos contra su índice (≥ 3 picos; piso 0.5 % [supuesto]);
#   "period": período del peine de picos (period_gain.py), error por réplicas de Poisson (piso 0.5 %). Converge a −0.5 V donde c1 no (seguimiento 9:
#             31/31 contra 5/31) y lee ~3.5-3.9 % más que c1 a nominal: un factor común cambia k pero no V_bd.
# Semilla / G esperada de cada toma: G de c1 × (OV0 + ΔV)/OV0 con OV0 = --ov0 (1.24 V, seguimiento 8 [supuesto]) o el seed_scale de la toma.
# Voltaje de cada canal: VMon medio de su canal SY5527 (grupo HV del mapeo, por omisión r4) en la ventana de la toma (index.json: primer evento a fin),
# del slow control de 1 s (/data/daq/slowcontrol/YYYY/hv_YYYYMMDD.csv[.gz]).
# Ajuste por canal y método: G = k (V − V_bd), mínimos cuadrados ponderados (1/σ_G²; --unweighted: sin pesos) con ≥ 2 escalones distintos; con
# ndf > 0 los errores se escalan por √(χ²/ndf) si > 1; V_bd = −b/k con la covarianza; OV nominal = VMon de la toma ΔV = 0 (media) − V_bd.
# Resumen por grupo HV y por estación: mediana y media ponderada de V_bd, OV y k.
# Tomas: --biasscan-results <biasscan_results.json> (steps[].sw_run.segments y delta_V; los escalones omitidos con sw_run también) y/o
# --take etiqueta:ΔV:patrón_de_segmentos[:seed_scale] (varias veces).
# Salida: <out>/vbd_channels.csv, vbd_summary.json, takes.json, gain/S_<toma>_{spec,fits}.csv.
# Uso: vbd_analysis.py --out DIR [--biasscan-results F] [--take ...] [--repo R] [--ov0 1.24] [--unweighted] [--reuse-spectra]
import argparse, csv, glob, gzip, json, math, os, statistics as st, subprocess, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from period_gain import load_spec, period_gain  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--biasscan-results", action="append", default=[])
ap.add_argument("--take", action="append", default=[])
ap.add_argument("--repo", default=os.path.abspath(os.path.join(HERE, "..", "..", "..", "..")))
ap.add_argument("--c1", default=None)
ap.add_argument("--mapping", default=None)
ap.add_argument("--ov0", type=float, default=1.24)
ap.add_argument("--unweighted", action="store_true")
ap.add_argument("--reuse-spectra", action="store_true")
ap.add_argument("--slowcontrol", default="/data/daq/slowcontrol")
a = ap.parse_args()
REPO = a.repo
C1 = a.c1 or f"{REPO}/daq/config/calibration/pe_calibration_2026_10_c1.json"
MAPF = a.mapping or f"{REPO}/daq/config/mapping/hodo_mapping_2026_10_r4.json"
os.makedirs(f"{a.out}/gain", exist_ok=True)
c1 = {e["ch"]: e.get("G_mV") for e in json.load(open(C1))["channels"]}
c1med = float(np.median([v for c, v in c1.items() if c < 32 and isinstance(v, (int, float))]))
G1 = {c: (c1.get(c) if isinstance(c1.get(c), (int, float)) else c1med) for c in range(32)}
HVG = {c["dig_ch"]: c["hv_group"] for c in json.load(open(MAPF))["channels"] if c["dig_ch"] < 32}
CH = [c for c in range(32) if c != 3 and HVG.get(c)]

# ---------------------------------------------------------------- tomas
takes = []
for f in a.biasscan_results:
    d = json.load(open(f))
    camp = os.path.basename(os.path.dirname(os.path.abspath(f)))
    for s in list(d.get("steps", [])) + [x.get("partial", {}) for x in d.get("skipped_steps", [])]:
        sw = s.get("sw_run") or {}
        if sw.get("segments"):
            takes.append({"label": f"{camp}_dV{s['delta_V']:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p"), "delta_V": float(s["delta_V"]),
                          "segments": sw["segments"], "seed_scale": None})
for t in a.take:
    p = t.split(":")
    takes.append({"label": p[0], "delta_V": float(p[1]), "segments": sorted(glob.glob(p[2])), "seed_scale": float(p[3]) if len(p) > 3 else None})
assert takes, "sin tomas"
for t in takes:
    if t["seed_scale"] is None:
        t["seed_scale"] = max(0.15, (a.ov0 + t["delta_V"]) / a.ov0)
    idx = t["segments"][0].rsplit("_seg", 1)[0] + ".index.json"
    j = json.load(open(idx))
    t["t0"], t["t1"], t["n_events"] = j["first_event_utc_ns"] / 1e9, j["end_utc_ns"] / 1e9, j["n_events"]


def sc_rows(day):
    y = day[:4]
    for p in (f"{a.slowcontrol}/{y}/hv_{day}.csv", f"{a.slowcontrol}/{y}/hv_{day}.csv.gz"):
        if os.path.exists(p):
            yield from csv.DictReader(gzip.open(p, "rt") if p.endswith(".gz") else open(p))
            return


for t in takes:
    acc = {c: [] for c in range(8)}
    day, end = t["t0"] - t["t0"] % 86400, t["t1"]
    while day <= end:
        for r in sc_rows(time.strftime("%Y%m%d", time.gmtime(day))):
            tt = int(r["t_resp_utc_ns"]) / 1e9
            c = int(r["ch"])
            if c < 8 and t["t0"] <= tt <= t["t1"]:
                acc[c].append((float(r["V0Set_V"]), float(r["VMon_V"])))
        day += 86400
    t["V0Set"] = {c: (st.mean(x[0] for x in v) if v else None) for c, v in acc.items()}
    t["VMon"] = {c: (st.mean(x[1] for x in v) if v else None) for c, v in acc.items()}
    t["n_sc"] = {c: len(v) for c, v in acc.items()}

# ---------------------------------------------------------------- G por toma
env = dict(os.environ)
env["LD_LIBRARY_PATH"] = os.path.expanduser("~/root6/lib") + ":" + env.get("LD_LIBRARY_PATH", "")
root = os.path.expanduser("~/root6/bin/root")
for t in takes:
    pre = f"{a.out}/gain/S_{t['label']}"
    if not (a.reuse_spectra and os.path.exists(pre + "_spec.csv")):
        seeds = f"{a.out}/gain/seeds_{t['label']}.json"
        json.dump({"channels": [{"ch": c, "G_mV": G1[c] * t["seed_scale"]} for c in range(32)]}, open(seeds, "w"), indent=1)
        pattern = t["segments"][0].replace("_seg001.root", "_seg*.root") if len(t["segments"]) > 1 else t["segments"][0]
        cmd = [root, "-l", "-b", "-q", "-e", f'gSystem->AddIncludePath("-I{REPO}/daq/core/include")',
               f'{HERE}/vbd_gain.C+("{pattern}", "{seeds}", "{pre}")']
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=f"{a.out}/gain", timeout=3600)
        open(pre + ".log", "w").write(r.stdout + r.stderr)
        if r.returncode != 0 or not os.path.exists(pre + "_spec.csv"):
            sys.exit(f"vbd_gain.C falló en {t['label']}: ver {pre}.log")
    sp, b = load_spec(pre + "_spec.csv")
    t["G"] = {"c1": {}, "period": {}}
    for r in csv.DictReader(open(pre + "_fits.csv")):
        c = int(r["ch"])
        if r["options"] != "c1" or r["status"] not in ("ok", "changed") or c not in CH:
            continue
        pos = [float(x) for x in r["positions_mV"].split()]
        g = float(r["G_mV"])
        if len(pos) >= 3:
            i = np.arange(len(pos))
            res = np.array(pos) - np.polyval(np.polyfit(i, pos, 1), i)
            se = math.sqrt((res ** 2).sum() / (len(pos) - 2) / ((i - i.mean()) ** 2).sum())
        else:
            se = 0.0
        t["G"]["c1"][c] = (g, max(se, 0.005 * g))
    for c in CH:
        g, sd, n = period_gain(sp[c], b, G1[c] * t["seed_scale"])
        if g is not None and sd is not None:
            t["G"]["period"][c] = (g, max(sd, 0.005 * g))
    print(f"toma {t['label']}: ΔV {t['delta_V']:+.2f} V, {t['n_events']} eventos, semilla {t['seed_scale']:.3f} × c1; G aceptadas: c1 {len(t['G']['c1'])}, "
          f"período {len(t['G']['period'])}; VMon grupos 1-8 " + " ".join(f"{t['VMon'][c]:.3f}" if t["VMon"][c] is not None else "-" for c in range(8)))


# ---------------------------------------------------------------- ajuste por canal
def linfit(v, g, s):
    v, g, s = map(np.asarray, (v, g, s))
    w = np.ones_like(s) if a.unweighted else 1 / s ** 2
    A = np.vstack([v, np.ones_like(v)]).T
    C = np.linalg.inv(A.T @ (A * w[:, None]))
    k, b0 = C @ (A.T @ (w * g))
    ndf = len(v) - 2
    chi = float((((g - k * v - b0) ** 2) / s ** 2).sum())
    if a.unweighted:  # covarianza de mínimos cuadrados ordinarios con σ_G propagados
        C = C @ (A.T @ (A * (s ** 2)[:, None])) @ C
    if ndf > 0 and chi / ndf > 1:
        C = C * chi / ndf
    vbd = -b0 / k
    J = np.array([b0 / k ** 2, -1 / k])  # dV_bd/dk, dV_bd/db
    evbd = float(math.sqrt(J @ C @ J))
    return {"k_mV_per_V": float(k), "k_err": float(math.sqrt(C[0, 0])), "V_bd_V": float(vbd), "V_bd_err_V": evbd, "chi2": chi, "ndf": ndf}


rows = []
nom_takes = [t for t in takes if abs(t["delta_V"]) < 1e-9]
for c in CH:
    hv = HVG[c] - 1
    r = {"ch": c, "hv_group": HVG[c]}
    vnom = [t["VMon"][hv] for t in nom_takes if t["VMon"][hv] is not None]
    r["V_nom_VMon_V"] = st.mean(vnom) if vnom else None
    for m in ("c1", "period"):
        pts = [(t["VMon"][hv], *t["G"][m][c], t["label"], t["delta_V"]) for t in takes if c in t["G"][m] and t["VMon"][hv] is not None]
        r[f"{m}_points"] = [{"take": p[3], "VMon_V": p[0], "G_mV": p[1], "G_err_mV": p[2]} for p in pts]
        nv = len({round(p[4], 3) for p in pts})  # escalones distintos (ΔV de la toma), no VMon redondeado: las tomas nominales difieren ~0.02 V
        r[f"{m}_n_voltages"] = nv
        if nv >= 2:
            f = linfit([p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts])
            if r["V_nom_VMon_V"] is not None:
                f["OV_nominal_V"] = r["V_nom_VMon_V"] - f["V_bd_V"]
            r[m] = f
    rows.append(r)


def summarize(sel, m, key):
    v = [(r[m][key], r[m]["V_bd_err_V"] if key != "k_mV_per_V" else r[m]["k_err"]) for r in sel if m in r and key in r[m]]
    if not v:
        return None
    x, e = np.array([p[0] for p in v]), np.array([p[1] for p in v])
    w = 1 / np.maximum(e, 1e-6) ** 2
    return {"n": len(v), "median": float(np.median(x)), "wmean": float((w * x).sum() / w.sum()), "wmean_err": float(1 / math.sqrt(w.sum())), "sd": float(np.std(x))}


summary = {"takes": [{k: t[k] for k in ("label", "delta_V", "n_events", "seed_scale", "VMon", "V0Set")} | {"n_G": {m: len(t["G"][m]) for m in t["G"]}} for t in takes],
           "ov0_assumed_V": a.ov0, "weighted": not a.unweighted, "by_group": {}, "station": {}}
for m in ("c1", "period"):
    summary["station"][m] = {key: summarize(rows, m, key) for key in ("V_bd_V", "OV_nominal_V", "k_mV_per_V")}
    summary["station"][m]["n_channels_fitted"] = sum(m in r for r in rows)
    for g in sorted({r["hv_group"] for r in rows}):
        summary["by_group"].setdefault(str(g), {})[m] = {key: summarize([r for r in rows if r["hv_group"] == g], m, key) for key in ("V_bd_V", "OV_nominal_V", "k_mV_per_V")}
json.dump(summary, open(f"{a.out}/vbd_summary.json", "w"), indent=1, default=float)
json.dump([{k: v for k, v in t.items() if k != "G"} | {"G": {m: {str(c): list(x) for c, x in t["G"][m].items()} for m in t["G"]}} for t in takes],
          open(f"{a.out}/takes.json", "w"), indent=1, default=float)
with open(f"{a.out}/vbd_channels.csv", "w") as f:
    keys = ["ch", "hv_group", "V_nom_VMon_V"]
    for m in ("c1", "period"):
        keys += [f"{m}_n_voltages", f"{m}_k_mV_per_V", f"{m}_k_err", f"{m}_V_bd_V", f"{m}_V_bd_err_V", f"{m}_OV_nominal_V", f"{m}_chi2", f"{m}_ndf"]
    f.write(",".join(keys) + "\n")
    for r in rows:
        vals = [r["ch"], r["hv_group"], r["V_nom_VMon_V"]]
        for m in ("c1", "period"):
            fm = r.get(m, {})
            vals += [r[f"{m}_n_voltages"], fm.get("k_mV_per_V"), fm.get("k_err"), fm.get("V_bd_V"), fm.get("V_bd_err_V"), fm.get("OV_nominal_V"), fm.get("chi2"), fm.get("ndf")]
        f.write(",".join("" if v is None else (f"{v:.5f}" if isinstance(v, float) else str(v)) for v in vals) + "\n")
for m in ("c1", "period"):
    s = summary["station"][m]
    print(f"método {m}: {s['n_channels_fitted']} canales con ajuste; estación: V_bd mediana "
          + (f"{s['V_bd_V']['median']:.3f} V, OV nominal mediana {s['OV_nominal_V']['median']:.3f} V (media ponderada {s['OV_nominal_V']['wmean']:.3f} ± "
             f"{s['OV_nominal_V']['wmean_err']:.3f}), k mediana {s['k_mV_per_V']['median']:.3f} mV/V" if s["V_bd_V"] else "sin ajustes"))
    for r in rows:
        if m in r:
            f_ = r[m]
            print(f"   {m} ch{r['ch']:02d} grupo {r['hv_group']}: {r[f'{m}_n_voltages']} voltajes, k {f_['k_mV_per_V']:.3f} ± {f_['k_err']:.3f} mV/V, V_bd {f_['V_bd_V']:.3f} ± "
                  f"{f_['V_bd_err_V']:.3f} V, OV {f_.get('OV_nominal_V', float('nan')):.3f} V, χ² {f_['chi2']:.2f}/{f_['ndf']}")
