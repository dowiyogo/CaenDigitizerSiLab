#!/usr/bin/python3 -I
# ampcal.py — seguimiento 2026-11-15: consistencia de la calibración de amplitud de hod2018 entre mayo (v8, tiempo bueno) y octubre (run 8, DAQ v2).
# Solo lectura de salidas derivadas. Entradas:
#  - seguimiento 10 (outputs/hod2018_mip_tagged/2026-11-01_tagged/work): {may_v8,oct_run8}_charge.bin (aS, QS por evento) y _dark.csv (G usada,
#    Q1 = mediana de la carga de pulsos oscuros aislados de 1 p.e., regiones [16,37) ∪ [102,132));
#  - este seguimiento (<out>/work): {pre}_dark2.bin (pulsos oscuros aislados con a y q, dark2_extract.C) y {pre}_darkfit_r0n{14,25}.csv
#    (fit_dark2.C: picos de amplitud de los pulsos aislados de la región pre-trigger [16,37), n = 1-4 o 2-5);
#  - G: PROC (hd64_{may_v8,oct8}_peaks.csv), c1 (may_common.g_c1), seguimiento 12 (work/{may8,oct8}_G_{dark,spec}.csv, run entero, hour = −1).
# (a) Por canal: R = (amplitud/G)/(carga/Q1) para
#     - muones etiquetados del intervalo lejos 29-54 mm (selección de mip_tagged.py: etiqueta = tira de la coordenada ortogonal del intervalo con
#       aS > 10 G; sonda = la tira con aS > 3 G; mezclados = etiqueta del evento siguiente, × N/(N − 1)). R_mu = Σ_neto(aS/G)/Σ_neto(QS/Q1),
#       neto = mismo − mezclados; también la mediana por impacto (mismo evento) y ρ = Σ_neto aS / Σ_neto QS (sin calibración, 1/ns);
#     - pulsos oscuros aislados de la región pre-trigger [16,37) con a en [(n − ½)G, (n + ½)G], n = 1, 2: R_n = mediana(a/G)/mediana(q/Q1),
#       ρ_n = mediana(a)/mediana(q). Variante: Q1 solo de la región pre-trigger (Q1_r0).
# (b) G de octubre: c1 (la del análisis) contra los reajustes de run 8: seguimiento 12 «dark» (pulsos de la región pre-trigger [20,49), picos 2-5,
#     p0 libre), seguimiento 12 «spec» y PROC (método de PROC, tramo [0,49) ∪ [102,144)), y los de este seguimiento (pulsos aislados de [16,37),
#     picos 2-5 = iso25, 1-4 = iso14). Lo mismo en mayo con la G del análisis (reajuste de PROC) como referencia. Ajuste aceptado: estado 0 o 1
#     (PROC: ok/changed) y error > 0; los canales sin G aceptada se llenan con G_análisis × mediana(G_juego/G_análisis).
# (c) Mediana de la estación de la amplitud de los muones etiquetados del intervalo lejos en p.e. medidos, aS/G con cada juego de G (sonda > 3 p.e.
#     del juego; etiqueta con la G del análisis para que la muestra de eventos no cambie): histograma neto de todas las tiras vivas (bins 0.05 p.e.).
# Salida: <out>/ampcal_channels.csv, ampcal_gains.csv, ampcal_station.csv, ampcal.json. Uso: ampcal.py <out dir>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

OUT = sys.argv[1]
WK = OUT + "/work"
F10 = os.path.join(REPO, "outputs/hod2018_mip_tagged/2026-11-01_tagged/work")
F12 = os.path.join(REPO, "outputs/hod_stability/2026-11-03_followup12/work")
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD = 3
LIVE = [c for c in range(32) if c != DEAD]
D = 37.5 - 15.0 + 2.0 * np.arange(16)
FAR = (29, 54)
TAG, PROBE_T = 10.0, 3.0
PRE = {"may": "may_v8", "oct": "oct_run8"}
DT = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("aC", "<f4", 32), ("QS", "<f4", 32), ("QC", "<f4", 32)])
DT2 = np.dtype([("ev", "<u4"), ("ch", "u1"), ("reg", "u1"), ("k", "<u2"), ("a", "<f4"), ("q", "<f4")])


def dark_csv(pre):
    g, q1 = np.full(32, np.nan), np.full(32, np.nan)
    for r in csv.DictReader(open(f"{F10}/{pre}_dark.csv")):
        c = int(r["ch"])
        g[c] = float(r["G_mV"])
        if r["Q1_median_mVns"] and c != DEAD:
            q1[c] = float(r["Q1_median_mVns"])
    return g, q1


def proc_peaks(name):
    g = np.full(32, np.nan)
    for r in csv.DictReader(open(f"{PK}/hd64_{name}_peaks.csv")):
        if int(r["ch"]) < 32 and r["status"] in ("ok", "changed"):
            g[int(r["ch"])] = float(r["G_fit"])
    return g


def f12(tag, src):
    g = np.full(32, np.nan)
    for r in csv.DictReader(l for l in open(f"{F12}/{tag}_G_{src}.csv") if not l.startswith("#")):
        if int(r["hour_utc_s"]) == -1 and int(r["fit_status"]) in (0, 1) and float(r["G_err_mV"]) > 0:
            g[int(r["ch"])] = float(r["G_mV"])
    return g


def iso(pre, v):
    g = np.full(32, np.nan)
    for r in csv.DictReader(open(f"{WK}/{pre}_darkfit_{v}.csv")):
        if r["Gf_mV"] and int(r["status"]) in (0, 1) and float(r["Gf_err_mV"]) > 0:
            g[int(r["ch"])] = float(r["Gf_mV"])
    return g


GU, Q1 = {}, {}
for camp in PRE:
    GU[camp], Q1[camp] = dark_csv(PRE[camp])
c1 = mc.g_c1()[:32]
SETS = {"may": {"refit_PROC (analysis)": GU["may"], "f12_dark": f12("may8", "dark"), "f12_spec": f12("may8", "spec"),
                "iso25": iso("may_v8", "r0n25"), "iso14": iso("may_v8", "r0n14")},
        "oct": {"c1 (analysis)": GU["oct"], "c1_raw": c1, "refit_PROC": proc_peaks("oct8"), "f12_dark": f12("oct8", "dark"), "f12_spec": f12("oct8", "spec"),
                "iso25": iso("oct_run8", "r0n25"), "iso14": iso("oct_run8", "r0n14")}}


def filled(camp, g):
    ok = np.isfinite(g) & (np.arange(32) != DEAD)
    r = np.median(g[ok] / GU[camp][ok])
    return np.where(ok, g, GU[camp] * r), float(r), [c for c in LIVE if not ok[c]]


# --- datos de eventos
EV = {}
for camp, pre in PRE.items():
    e = np.fromfile(f"{F10}/{pre}_charge.bin", dtype=DT)
    e = e[np.argsort(e["t"], kind="stable")]
    if camp == "may":
        e = e[mc.good_time(e["t"], "v8")]
    EV[camp] = e
assert np.allclose(GU["oct"][np.isfinite(c1)], c1[np.isfinite(c1)], atol=6e-6, rtol=0), "G de octubre del seguimiento 10 ≠ c1"


def readout_sign(probe, layer):
    s = SIGN[f"oct2026_run8/{probe}{layer}"]["sign"]
    assert s == SIGN[f"may2026_v8_goodtime/{probe}{layer}"]["sign"]
    return s


def far_tags(camp):
    """Por canal sonda: (máscara de etiqueta del mismo evento, del evento siguiente)."""
    e, o, g = EV[camp], ORDER[camp], GU[camp]
    out = {}
    for probe, ref in (("X", "Y"), ("Y", "X")):
        refch = np.asarray(o[f"h18_{ref}"])
        tagS = e["aS"][:, refch] > TAG * g[refch]
        tagS[:, refch == DEAD] = False
        for L in ("A", "B"):
            d = D if readout_sign(probe, L) < 0 else 75.0 - D
            js = [j for j in range(16) if FAR[0] <= d[j] < FAR[1] + 1e-9]
            T = tagS[:, js].any(1)
            Tn = np.concatenate([T[1:], [False]])
            for i, c in enumerate(o[f"h18_{probe}"]):
                if i % 2 == (0 if L == "A" else 1) and c != DEAD:
                    out[c] = (T, Tn, f"{probe}{L}", i)
    return out


TAGS = {camp: far_tags(camp) for camp in PRE}
DARK = {camp: np.fromfile(f"{WK}/{pre}_dark2.bin", dtype=DT2) for camp, pre in PRE.items()}


def muon_ch(camp, c, g, q1):
    e = EV[camp]
    T, Tn, lay, pos = TAGS[camp][c]
    N = len(e)
    s = N / (N - 1)
    a, q = e["aS"][:, c], e["QS"][:, c]
    sel = a > PROBE_T * g
    ms, mm = sel & T, sel & Tn
    na = (a[ms] / g).sum() - s * (a[mm] / g).sum()
    nq = (q[ms] / q1).sum() - s * (q[mm] / q1).sum()
    ra = a[ms].sum() - s * a[mm].sum()
    rq = q[ms].sum() - s * q[mm].sum()
    r_hit = (a[ms] / g) / (q[ms] / q1)
    # error de R_mu: bootstrap de los impactos del mismo evento (los mezclados son pocos)
    rng = np.random.default_rng(2026111503 + c + (0 if camp == "may" else 100))
    xa, xq = a[ms] / g, q[ms] / q1
    bs = [xa[i].sum() / xq[i].sum() for i in (rng.integers(0, len(xa), len(xa)) for _ in range(200))]
    return {"layer": lay, "strip_position": pos, "n_same": int(ms.sum()), "n_mix_scaled": s * int(mm.sum()), "R_mu": na / nq, "R_mu_err": float(np.std(bs)),
            "R_mu_median_hit": float(np.median(r_hit)), "rho_mu_per_ns": ra / rq, "amp_pe_mean_net": na / (ms.sum() - s * mm.sum()),
            "charge_pe_mean_net": nq / (ms.sum() - s * mm.sum())}


def dark_ch(camp, c, g, q1):
    d = DARK[camp]
    d = d[(d["ch"] == c) & (d["reg"] == 0)]
    out = {}
    s1 = d[(d["a"] >= 0.5 * g) & (d["a"] <= 1.5 * g)]
    q1r0 = float(np.median(s1["q"]))
    out["Q1_r0_mVns"] = q1r0
    for n in (1, 2):
        s = d[(d["a"] >= (n - 0.5) * g) & (d["a"] <= (n + 0.5) * g)]
        ma, mq = float(np.median(s["a"])), float(np.median(s["q"]))
        out[f"n_{n}pe"] = int(len(s))
        out[f"amp_{n}pe_over_G"] = ma / g
        out[f"charge_{n}pe_over_Q1"] = mq / q1
        out[f"R_{n}pe"] = (ma / g) / (mq / q1)
        out[f"rho_{n}pe_per_ns"] = ma / mq
        out[f"R_{n}pe_err"] = 1.2533 * float(np.std(s["a"] / g)) / np.sqrt(len(s)) / (mq / q1)
    return out


res = {"definitions": "ver cabecera de ampcal.py", "events": {k: int(len(v)) for k, v in EV.items()}, "gain_sets": {}, "station": {}}
rows = []
for camp in PRE:
    g, q1 = GU[camp], Q1[camp]
    for c in LIVE:
        m = muon_ch(camp, c, g[c], q1[c])
        dk = dark_ch(camp, c, g[c], q1[c])
        rows.append({"campaign": camp, "ch": c, "G_used_mV": g[c], "Q1_mVns": q1[c], **m, **dk})
cols = list(rows[0].keys())
with open(OUT + "/ampcal_channels.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(cols)
    for r in rows:
        w.writerow([f"{r[k]:.5g}" if isinstance(r[k], float) else r[k] for k in cols])

# resumen por campaña y por barra (misma MPPC = misma posición de tira en el plano)
by = {(r["campaign"], r["ch"]): r for r in rows}
summ = {}
for camp in PRE:
    rr = [r for r in rows if r["campaign"] == camp]
    summ[camp] = {k: [float(np.median([r[k] for r in rr])), float(np.percentile([r[k] for r in rr], 16)), float(np.percentile([r[k] for r in rr], 84))]
                  for k in ("R_mu", "R_mu_median_hit", "rho_mu_per_ns", "R_1pe", "R_2pe", "rho_1pe_per_ns", "rho_2pe_per_ns", "amp_1pe_over_G", "amp_2pe_over_G",
                            "charge_1pe_over_Q1", "charge_2pe_over_Q1", "amp_pe_mean_net", "charge_pe_mean_net")}
    summ[camp]["Q1_r0_over_Q1"] = float(np.median([r["Q1_r0_mVns"] / r["Q1_mVns"] for r in rr]))
pairs = []
for coord in ("X", "Y"):
    for j in range(16):
        cm, co = ORDER["may"][f"h18_{coord}"][j], ORDER["oct"][f"h18_{coord}"][j]
        if DEAD in (cm, co):
            continue
        pairs.append((coord, j, cm, co))
ratio = {k: [by[("oct", co)][k] / by[("may", cm)][k] for _, _, cm, co in pairs]
         for k in ("R_mu", "rho_mu_per_ns", "R_1pe", "R_2pe", "rho_1pe_per_ns", "amp_pe_mean_net", "charge_pe_mean_net", "G_used_mV", "Q1_mVns")}
summ["oct_over_may_by_bar"] = {k: [float(np.median(v)), float(np.percentile(v, 16)), float(np.percentile(v, 84)), len(v)] for k, v in ratio.items()}
res["per_channel_summary_median_p16_p84"] = summ
with open(OUT + "/ampcal_bars.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["coord", "strip_position", "ch_may", "ch_oct"] + [f"{k}_{c}" for k in ("R_mu", "R_1pe", "R_2pe", "rho_mu_per_ns") for c in ("may", "oct")])
    for coord, j, cm, co in pairs:
        w.writerow([coord, j, cm, co] + [f"{by[(c, ch)][k]:.5g}" for k in ("R_mu", "R_1pe", "R_2pe", "rho_mu_per_ns") for c, ch in (("may", cm), ("oct", co))])


# (b) juegos de G y (c) medianas de la estación
def station_median(camp, g, var="amp"):
    e = EV[camp]
    N = len(e)
    s = N / (N - 1)
    edg = np.arange(0, 300.0001, 0.05)
    hs, hm = np.zeros(len(edg) - 1), np.zeros(len(edg) - 1)
    for c in LIVE:
        T, Tn, _, _ = TAGS[camp][c]
        x = e["aS"][:, c] / g[c] if var == "amp" else e["QS"][:, c] / Q1[camp][c]
        sel = x > PROBE_T
        hs += np.histogram(x[sel & T], bins=edg)[0]
        hm += np.histogram(x[sel & Tn], bins=edg)[0]
    net = np.clip(hs - s * hm, 0, None)
    cs = np.cumsum(net)
    k = int(np.argmax(cs >= 0.5 * cs[-1]))
    # interpolación lineal dentro del bin
    prev = cs[k - 1] if k else 0.0
    med = edg[k] + (0.5 * cs[-1] - prev) / max(net[k], 1e-12) * (edg[k + 1] - edg[k])
    # error: bootstrap de Poisson de los histogramas
    rng = np.random.default_rng(2026111504)
    bs = []
    for _ in range(200):
        n2 = np.clip(rng.poisson(hs) - s * rng.poisson(hm), 0, None)
        c2 = np.cumsum(n2)
        k2 = int(np.argmax(c2 >= 0.5 * c2[-1]))
        bs.append(edg[k2] + 0.025)
    return float(med), float(np.std(bs)), float(hs.sum()), float(s * hm.sum())


grows, srows = [], []
for camp, sets in SETS.items():
    res["gain_sets"][camp] = {}
    for name, g in sets.items():
        gf, rmed, miss = filled(camp, g)
        rel = [g[c] / GU[camp][c] for c in LIVE if np.isfinite(g[c])]
        res["gain_sets"][camp][name] = {"median_ratio_to_analysis_G": rmed, "p16_p84": [float(np.percentile(rel, 16)), float(np.percentile(rel, 84))],
                                        "filled_channels": miss, "ch22": None if not np.isfinite(g[22]) else float(g[22]),
                                        "ch31": None if not np.isfinite(g[31]) else float(g[31])}
        med, err, ns, nm = station_median(camp, gf)
        srows.append([camp, name, f"{rmed:.5f}", " ".join(map(str, miss)), f"{med:.4f}", f"{err:.4f}", f"{ns:.0f}", f"{nm:.1f}"])
        res["station"].setdefault(camp, {})[name] = {"median_amp_pe": med, "err": err, "n_same": ns, "n_mix_scaled": nm}
        for c in range(32):
            grows.append([camp, name, c, "" if not np.isfinite(g[c]) else f"{g[c]:.5f}", f"{gf[c]:.5f}", f"{GU[camp][c]:.5f}"])
    med, err, ns, nm = station_median(camp, None, "charge")
    srows.append([camp, "charge QS/Q1", "", "", f"{med:.4f}", f"{err:.4f}", f"{ns:.0f}", f"{nm:.1f}"])
    res["station"][camp]["charge QS/Q1"] = {"median_charge_pe": med, "err": err}
with open(OUT + "/ampcal_gains.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "gain_set", "ch", "G_mV", "G_filled_mV", "G_analysis_mV"])
    w.writerows(grows)
with open(OUT + "/ampcal_station.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "gain_set", "median_G_over_G_analysis", "filled_channels", "station_median_far_pe", "err", "n_same", "n_mix_scaled"])
    w.writerows(srows)
json.dump(res, open(OUT + "/ampcal.json", "w"), indent=1)

# impresión
for camp in PRE:
    print(f"== {camp}: eventos {len(EV[camp])}")
    for k, v in summ[camp].items():
        print(f"  {k:22s} " + (" ".join(f"{x:.4f}" for x in v) if isinstance(v, list) else f"{v:.4f}"))
print("== octubre / mayo por barra (mediana, p16, p84, n)")
for k, v in summ["oct_over_may_by_bar"].items():
    print(f"  {k:22s} {v[0]:.4f} {v[1]:.4f} {v[2]:.4f} {v[3]}")
print("== juegos de G y mediana de la estación (intervalo lejos, p.e. medidos)")
for r in srows:
    print("  " + " | ".join(r))
for camp in PRE:
    for name, v in res["gain_sets"][camp].items():
        print(f"  {camp} {name:24s} G/G_an {v['median_ratio_to_analysis_G']:.4f} (p16 {v['p16_p84'][0]:.4f}, p84 {v['p16_p84'][1]:.4f}) ch22 {v['ch22']} ch31 {v['ch31']} "
              f"llenados {v['filled_channels']}")
