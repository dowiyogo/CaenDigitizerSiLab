#!/usr/bin/python3 -I
# pe_primary.py — seguimiento 2026-10-31 (6), punto 1: espectro MIP de hod2018 en p.e. primarios (corregido por crosstalk), mayo contra octubre.
# Mismas entradas y definiciones que mip_gain_vs_mv.py (seguimiento 1): archivos de impactos de PROC, amplitud de pico aS (signal) y aC (fuera
# de tiempo) en la ventana de 144 ns, espectro sustraído = signal − fuera de tiempo (bins negativos a 0), umbral 5 p.e. con la G de cada canal,
# mayo v8 con el tiempo bueno nuevo (may_common.good_time, HV-DIAG-5 excluido), octubre run 8 completo, hod2018 sin ch3.
# G: mayo = reajuste de PROC (hd64_may_v8_peaks.csv), octubre = c1 (mediana para ch22 y ch31).
# Factor de crosstalk F = 1/(1 − p_ct) por canal y campaña, del seguimiento 5 (outputs/hod_dark/2026-10-31_pairs/pct_corrected_<muestra>.csv):
#   central = corregido por colas partidas, E < 0.5 G, sin apilamiento ("factor_corr"; f medido en ch4, 13, 19, 28 y mediana de f en el
#   resto), región pre-trigger (mayo v8 tiempo bueno nuevo, octubre run 8); variantes: sin corregir ("factor_orig"), con apilamiento
#   ("factor_corr_pile") y E < 0.75 G.
# Amplitud en p.e. primarios = aS / (G · F).
# Espectro sustraído con la definición del seguimiento 1 (la que dio 19.5 / 13.5 = 1.44): por canal, bins de 1 p.e. MEDIDO desde 0, umbral
# 5 p.e. medidos, signal − fuera de tiempo, bins negativos a 0. En p.e. primarios cada bin medido [n, n+1) se lleva a su centro (n + 0.5)/F
# del canal y la mediana (por canal o de la suma) es la mediana ponderada de esos centros: así la población y el recorte son exactamente los
# del seguimiento 1 y la razón en p.e. medidos reproduce 1.44. Variante "umbral primario": solo los bins con centro/F > 5.
# Sensibilidad: el mismo cálculo SIN recortar los bins negativos (mediana = primer centro donde la acumulada neta llega a la mitad). En los
# bins de 5-10 p.e. la ventana signal tiene menos pulsos que la de control (ch00 de mayo: 1371 contra 1517), así que sin recorte la parte baja
# queda negativa y la mediana sube; con recorte baja: las dos versiones acotan el efecto de la resta.
# Error de la mediana: remuestreo de Poisson de signal y control por canal (200 réplicas).
# Salida: <out>/pe_primary.json y pe_primary_channels.csv. Uso: pe_primary.py <out dir>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
PAIRS = os.path.join(REPO, "outputs/hod_dark/2026-10-31_pairs")
CH = [c for c in range(32) if c != 3]
gm = {int(r["ch"]): float(r["G_fit"]) for r in csv.DictReader(open(W + "/hd64_may_v8_peaks.csv")) if r["status"] in ("ok", "changed")}
c1 = json.load(open(os.path.join(REPO, "daq/config/calibration/pe_calibration_2026_10_c1.json")))
go = {e["ch"]: e["G_mV"] for e in c1["channels"] if e["ch"] < 32 and isinstance(e.get("G_mV"), (int, float))}
med_o = float(np.median([go[c] for c in CH if c in go]))
G = {"may": {c: gm[c] for c in CH}, "oct": {c: go.get(c, med_o) for c in CH}}


def factors(sample):
    F = {}
    for r in csv.DictReader(open(f"{PAIRS}/pct_corrected_{sample}.csv")):
        c, e = int(r["ch"]), float(r["Ecut_pe"])
        d = F.setdefault(c, {"f_source": r["f_source"], "none": 1.0})
        if abs(e - 0.5) < 1e-6:
            d.update(central=float(r["factor_corr"]), orig=float(r["factor_orig"]), pile=float(r["factor_corr_pile"]), f_spurious=float(r["f_spurious"]),
                     p_ct_corr=float(r["p_ct_corr"]))
        else:
            d["ecut075"] = float(r["factor_corr"])
    return F


F = {"may": factors("may_v8good_pre"), "oct": factors("oct_run8_pre")}
E = {"may": xy.load_hits(W + "/hd64_may_v8_hits.bin", 64), "oct": xy.load_hits(W + "/hd64_oct8_hits.bin", 64)}
E["may"] = E["may"][mc.good_time(E["may"]["t"], "v8")]
EDG1 = np.arange(0, 301, 1.0)
CEN1 = EDG1[:-1] + 0.5
RNG = np.random.default_rng(6)
NB = 200
# histogramas signal y control por campaña y canal en p.e. medidos, umbral 5 p.e. medidos
H = {}
for k in ("may", "oct"):
    for c in CH:
        e, g = E[k], G[k][c]
        s, o = e["aS"][:, c] / g, e["aC"][:, c] / g
        H[(k, c)] = (np.histogram(s[s > 5], bins=EDG1)[0].astype(float), np.histogram(o[o > 5], bins=EDG1)[0].astype(float))


def wmedian(vals, w, clip):
    """mediana ponderada; clip=False: primer valor donde la acumulada (con pesos negativos) llega a la mitad del total."""
    o = np.argsort(vals)
    v, ww = vals[o], w[o]
    cs = np.cumsum(ww)
    if cs[-1] <= 0:
        return None
    return float(v[int(np.argmax(cs >= 0.5 * cs[-1]))])


def med_set(k, chans, var, clip=True, thr_prim=False, hs_override=None):
    vals, ws = [], []
    for c in chans:
        hs, ho = H[(k, c)] if hs_override is None else hs_override[c]
        net = hs - ho
        if clip:
            net = np.clip(net, 0, None)
        x = CEN1 / F[k][c][var]
        m = x > 5 if thr_prim else np.ones_like(x, bool)
        vals.append(x[m]); ws.append(net[m])
    return wmedian(np.concatenate(vals), np.concatenate(ws), clip)


def med_err(k, chans, var, clip=True, thr_prim=False):
    v = []
    for _ in range(NB):
        ov = {c: (RNG.poisson(H[(k, c)][0]).astype(float), RNG.poisson(H[(k, c)][1]).astype(float)) for c in chans}
        x = med_set(k, chans, var, clip, thr_prim, ov)
        if x is not None:
            v.append(x)
    return float(np.std(v))


def netsum(k, chans):
    return float(sum(np.clip(H[(k, c)][0] - H[(k, c)][1], 0, None).sum() for c in chans)), float(sum((H[(k, c)][0] - H[(k, c)][1]).sum() for c in chans))


VAR = ("none", "central", "orig", "pile", "ecut075")
rows = []
for c in CH:
    r = {"ch": c, "G_may_mV": G["may"][c], "G_oct_mV": G["oct"][c], "G_ratio": G["may"][c] / G["oct"][c],
         "f_source_may": F["may"][c]["f_source"], "f_source_oct": F["oct"][c]["f_source"],
         "F_may_central": F["may"][c]["central"], "F_oct_central": F["oct"][c]["central"]}
    for k in ("may", "oct"):
        r[f"{k}_median_pe_meas"] = med_set(k, [c], "none")
        r[f"{k}_median_pe_prim"] = med_set(k, [c], "central")
        r[f"{k}_median_pe_prim_err"] = med_err(k, [c], "central")
        r[f"{k}_median_pe_meas_noclip"] = med_set(k, [c], "none", clip=False)
        r[f"{k}_net_hits"] = netsum(k, [c])[1]
    r["ratio_meas"] = r["may_median_pe_meas"] / r["oct_median_pe_meas"]
    r["ratio_prim"] = r["may_median_pe_prim"] / r["oct_median_pe_prim"]
    r["F_ratio"] = r["F_may_central"] / r["F_oct_central"]
    rows.append(r)

res = {"definition": "seguimiento 1 (bins de 1 p.e. medido, umbral 5 p.e. medidos, signal − fuera de tiempo, negativos a 0); p.e. primarios = "
                     "centro del bin / F, F = 1/(1 − p_ct) del seguimiento 5", "F_source": {"may": "pct_corrected_may_v8good_pre.csv",
                                                                                         "oct": "pct_corrected_oct_run8_pre.csv"},
       "channels": rows, "sums": {}}
SETS = {"all31": CH, "no_ch29_ch30": [c for c in CH if c not in (29, 30)], "ch11_ch26": [11, 26], "ch11": [11], "ch26": [26]}
for name, chans in SETS.items():
    d = {"channels": chans}
    for k in ("may", "oct"):
        d[f"{k}_net_hits_clipped"], d[f"{k}_net_hits"] = netsum(k, chans)
        for v in VAR:
            d[f"{k}_median_{v}"] = med_set(k, chans, v)
        d[f"{k}_median_none_err"] = med_err(k, chans, "none")
        d[f"{k}_median_central_err"] = med_err(k, chans, "central")
        d[f"{k}_median_central_thrprim"] = med_set(k, chans, "central", thr_prim=True)
        d[f"{k}_median_none_noclip"] = med_set(k, chans, "none", clip=False)
        d[f"{k}_median_central_noclip"] = med_set(k, chans, "central", clip=False)
        d[f"{k}_F_median_central"] = float(np.median([F[k][c]["central"] for c in chans]))
    for v in VAR:
        d[f"ratio_{v}"] = d[f"may_median_{v}"] / d[f"oct_median_{v}"]
    d["ratio_central_err"] = d["ratio_central"] * float(np.hypot(d["may_median_central_err"] / d["may_median_central"],
                                                                 d["oct_median_central_err"] / d["oct_median_central"]))
    d["ratio_none_err"] = d["ratio_none"] * float(np.hypot(d["may_median_none_err"] / d["may_median_none"], d["oct_median_none_err"] / d["oct_median_none"]))
    d["ratio_central_thrprim"] = d["may_median_central_thrprim"] / d["oct_median_central_thrprim"]
    d["ratio_none_noclip"] = d["may_median_none_noclip"] / d["oct_median_none_noclip"]
    d["ratio_central_noclip"] = d["may_median_central_noclip"] / d["oct_median_central_noclip"]
    d["ratio_F_median_central"] = d["may_F_median_central"] / d["oct_F_median_central"]
    d["ratio_F_effective"] = d["ratio_none"] / d["ratio_central"]
    d["G_ratio_median"] = float(np.median([G["may"][c] / G["oct"][c] for c in chans]))
    res["sums"][name] = d
json.dump(res, open(OUT + "/pe_primary.json", "w"), indent=1)
with open(OUT + "/pe_primary_channels.csv", "w") as f:
    keys = list(rows[0])
    f.write(",".join(keys) + "\n")
    for r in rows:
        f.write(",".join("" if r[k] is None else (f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])) for k in keys) + "\n")
for r in rows:
    print(f"ch{r['ch']:02d} G {r['G_may_mV']:.3f}/{r['G_oct_mV']:.3f} F mayo {r['F_may_central']:.4f} ({r['f_source_may']}) oct {r['F_oct_central']:.4f} "
          f"({r['f_source_oct']}) | mediana p.e. medidos {r['may_median_pe_meas']:.2f}/{r['oct_median_pe_meas']:.2f} | primarios "
          f"{r['may_median_pe_prim']:.2f}±{r['may_median_pe_prim_err']:.2f} / {r['oct_median_pe_prim']:.2f}±{r['oct_median_pe_prim_err']:.2f} "
          f"razón {r['ratio_prim']:.3f} | sin recorte (medidos) {r['may_median_pe_meas_noclip']}/{r['oct_median_pe_meas_noclip']}")
for name, d in res["sums"].items():
    print(name, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items()}))
