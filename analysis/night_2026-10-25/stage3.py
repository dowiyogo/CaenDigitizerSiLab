#!/usr/bin/python3 -I
# stage3.py — XPLUS-4b etapa 3: eficiencia contra ganancia dentro de un run (barrido sustituto).
#  - Bloques: horas UTC enteras (las de tiempo bueno en mayo) acumuladas hasta ≥ min_events eventos; un resto menor se une al bloque
#    anterior. G_b por canal y bloque con gain_blocks.C (gain_fit.h sobre el espectro "spec" de la ntupla, método del reajuste de PROC),
#    semilla = G del run. g_b = mediana_c(G_b) / mediana_c(G_run), con G_run del mismo método sobre todo el run (hod2018 sin ch3).
#  - q_b con xy_estimator.xy_fraction sobre el _hits.bin de la ntupla (idéntico al de PROC): modo A umbral 5·G_A (G del fixture: reajuste
#    de mayo para v8, c1 para run 8, G_run para v10), modo B umbral 5·G_b (canales sin ajuste aceptado: G_run·g_b). También T = 10.
#  - Ajuste ponderado q_b = A + B·(g_b − 1) → s = B/A, σ_s por covarianza (sin error de g_b), χ²/ndf, Pearson. |s|/σ_s < 2 → cota a 95 %:
#    |s| < |s| + 1.96 σ_s. dε/dV = s · dlnG/dV · ε, dlnG/dV = 1.104 V⁻¹ (XPLUS-4), ε = A/(f_2018·g).
#  - Controles por bloque: tasa de trigger (eventos / tiempo vivo, sin huecos > 20 s), DCR (mediana de pulsos > 0.5·G por evento de
#    gain_blocks.C / 91 muestras de 16 ns), VMon mediano de CH00-07 (slow control; octubre), ε(X|Y), ε(Y|X).
# Uso: stage3.py <run: v8|v10|run8> <ntupla sin extensión> <salida dir> [min_events]
import csv, gzip, json, math, os, statistics as st, subprocess, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

RUN, NT, OUT = sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])
MINEV = int(sys.argv[4]) if len(sys.argv) > 4 else 8000
os.makedirs(OUT, exist_ok=True)
H18 = [c for c in range(32) if c != 3]
XCH, YCH = [c for c in range(16) if c != 3], list(range(16, 32))
F18, GG, DLNG = 0.03522, 0.980, 1.104
fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))


def macro(call):
    r = subprocess.run(["/home/xplus/night_2026-10-25/run_macro.sh", call], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stdout + r.stderr)


def gains(blocks, seeds, tag):
    bf, of = f"{OUT}/{RUN}_{tag}_blocks.txt", f"{OUT}/{RUN}_{tag}_gain.csv"
    with open(bf, "w") as f:
        for lab, hrs in blocks:
            for h in hrs:
                f.write(f"{lab} {h} {h + 3600}\n")
    macro(f'gain_blocks.C+("{NT}.root", "{bf}", "{",".join(f"{x:.6f}" for x in seeds)}", "{of}")')
    res = {}
    for r in csv.DictReader(open(of)):
        res.setdefault(r["block"], {})[int(r["ch"])] = (float(r["G_mV"]), r["status"] in ("ok", "changed"), float(r["dark_gt0p5G_per_event"]), int(r["events"]))
    return res


def g_refit():
    g = np.full(32, np.nan)
    for l in open(os.path.join(REPO, "outputs/proceeding_2026_10/work/hd64_may_v8_peaks.csv")):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


e = xy.load_hits(NT + "_hits.bin", 64)
if RUN in ("v8", "v10"):
    e = e[mc.good_time(e["t"], RUN)]
t = e["t"]
hours = np.floor(t / 3600).astype(np.int64) * 3600
uh, cnt = np.unique(hours, return_counts=True)
# G del run (mismo método, todo el run / tiempo bueno) y G del modo A
seed0 = fill(g_refit()) if RUN != "run8" else fill(mc.g_c1())
whole = gains([("run", [int(h) for h in uh])], seed0, "run")["run"]
G_run = np.array([whole[c][0] if whole.get(c, (0, False))[1] else np.nan for c in range(32)])
G_run = fill(G_run)
G_A = fill(g_refit()) if RUN == "v8" else (fill(mc.g_c1()) if RUN == "run8" else G_run)
# bloques
blocks, cur, n = [], [], 0
for h, k in zip(uh, cnt):
    cur.append(int(h))
    n += int(k)
    if n >= MINEV:
        blocks.append(cur)
        cur, n = [], 0
if cur:
    if blocks:
        blocks[-1] += cur
    else:
        blocks.append(cur)
blk = [(f"b{i:02d}", b) for i, b in enumerate(blocks)]
bg = gains(blk, G_run, "blocks")
# slow control para octubre
vmon = {}
if RUN == "run8":
    rows = []
    for d in ("20261022", "20261023"):
        with gzip.open(f"/data/daq/slowcontrol/2026/hv_{d}.csv.gz", "rt") as f:
            next(f)
            for line in f:
                r = line.split(",")
                if int(r[4]) < 8:
                    rows.append((int(r[0]) * 1e-9, float(r[7])))
    vt = np.array([x[0] for x in rows]); vv = np.array([x[1] for x in rows])
med_run = st.median([G_run[c] for c in H18])
out = []
for lab, hrs in blk:
    a, b = hrs[0], hrs[-1] + 3600
    m = np.isin(hours, np.array(hrs))
    eb = e[m]
    G_b = np.array([bg[lab][c][0] if bg[lab].get(c, (0, False))[1] else np.nan for c in range(32)])
    ok = [c for c in H18 if np.isfinite(G_b[c])]
    rat = [G_b[c] / G_run[c] for c in ok]
    g_b = st.median(rat)  # mediana de las razones por canal (robusta a qué canales ajustaron)
    g_b_err = 1.2533 * st.pstdev(rat) / math.sqrt(len(rat)) if len(rat) > 1 else float("nan")
    G_B = np.where(np.isfinite(G_b), G_b, G_run * g_b)
    tb = eb["t"]
    dtb = np.diff(tb)
    live = float(dtb[dtb <= 20].sum())
    row = {"block": lab, "t0_utc": a, "t1_utc": b, "hours": len(hrs), "events": int(m.sum()), "g_b": g_b, "g_b_err": g_b_err, "n_ch_fit": len(ok), "live_s": live,
           "trigger_rate_Hz": (len(tb) - 1) / live if live > 0 else None,
           "DCR_proxy_cps": st.median([bg[lab][c][2] for c in H18]) / (91 * 16e-9)}
    if RUN == "run8":
        sel = np.isin(np.floor(vt / 3600).astype(np.int64) * 3600, np.array(hrs))
        row["VMon_median_V"] = float(np.median(vv[sel])) if sel.any() else None
    for T in (5, 10):
        for mode, Gm in (("A", G_A), ("B", G_B)):
            G64 = np.concatenate([Gm, np.full(32, np.nan)])
            Xs, Xc = xy.plane_hits(eb, G64, T, XCH)
            Ys, Yc = xy.plane_hits(eb, G64, T, YCH)
            r = xy.xy_fraction(Xs, Ys, Xc, Yc)
            row[f"q_{mode}_T{T}"] = r["f_mix"]
            row[f"q_{mode}_T{T}_err"] = r["f_mix_err"]
            row[f"qctl_{mode}_T{T}"] = r["f_ctl"]
            row[f"eXgY_{mode}_T{T}"], row[f"eXgY_{mode}_T{T}_err"] = r["e_X|Y_mix"], r["e_X|Y_mix_err"]
            row[f"eYgX_{mode}_T{T}"], row[f"eYgX_{mode}_T{T}_err"] = r["e_Y|X_mix"], r["e_Y|X_mix_err"]
    out.append(row)


def fit(x, y, s, sx=None, is_q=True):
    x, y, s = map(np.asarray, (x, y, s))
    sx = np.zeros_like(x) if sx is None else np.asarray(sx)
    X = np.vstack([np.ones_like(x), x - 1]).T
    B = 0.0
    for _ in range(5):  # varianza efectiva: σ² = σ_q² + (B·σ_g)²
        w = 1 / (s ** 2 + (B * sx) ** 2)
        C = np.linalg.inv(X.T @ (X * w[:, None]))
        A, B = C @ (X.T @ (w * y))
    chi2 = float(np.sum(w * (y - A - B * (x - 1)) ** 2))
    sv = B / A
    var = C[1, 1] / A ** 2 + B ** 2 * C[0, 0] / A ** 4 - 2 * B * C[0, 1] / A ** 3
    ss = math.sqrt(var)
    r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else None
    o = {"A": float(A), "A_err": math.sqrt(C[0, 0]), "s": float(sv), "s_err": ss, "chi2": chi2, "ndf": len(x) - 2, "pearson": r,
         "significance": abs(sv) / ss}
    if abs(sv) / ss < 2:
        o["upper_95_abs_s"] = abs(sv) + 1.96 * ss
    if not is_q:
        return o
    eps = A / (F18 * GG)
    o["eps0"] = float(eps)
    o["deps_dV_per_V"] = float(sv * DLNG * eps)
    o["deps_dV_err"] = float(ss * DLNG * eps)
    for dv in (0.5, 1.0):
        o[f"eps_at_+{dv}V_linear"] = float(eps * (1 + sv * DLNG * dv))
        o[f"eps_at_+{dv}V_err"] = float(eps * DLNG * dv * ss)
    return o


g = [r["g_b"] for r in out]
gerr = [r["g_b_err"] for r in out]
fits = {}
for T in (5, 10):
    for mode in ("A", "B"):
        fits[f"q_{mode}_T{T}"] = fit(g, [r[f"q_{mode}_T{T}"] for r in out], [r[f"q_{mode}_T{T}_err"] for r in out], gerr)
        for k in ("eXgY", "eYgX"):
            fits[f"{k}_{mode}_T{T}"] = fit(g, [r[f"{k}_{mode}_T{T}"] for r in out], [r[f"{k}_{mode}_T{T}_err"] for r in out], gerr, is_q=False)
ctrl = {}
for var in ("trigger_rate_Hz", "DCR_proxy_cps") + (("VMon_median_V",) if RUN == "run8" else ()):
    xs = [r[var] for r in out]
    ctrl[var] = {"pearson_qB5": float(np.corrcoef(xs, [r["q_B_T5"] for r in out])[0, 1]), "pearson_qA5": float(np.corrcoef(xs, [r["q_A_T5"] for r in out])[0, 1]),
                 "pearson_with_g": float(np.corrcoef(xs, g)[0, 1]), "range": [min(xs), max(xs)]}
summary = {"run": RUN, "min_events": MINEV, "n_blocks": len(out), "g_range": [min(g), max(g)], "g_span": max(g) - min(g),
           "G_run_median_mV": med_run, "blocks": out, "fits": fits, "controls": ctrl,
           "notes": "error de g_b no incluido en el ajuste; temperatura de placa: no registrada en results.json ni en los archivos de run 8"}
json.dump(summary, open(f"{OUT}/{RUN}_stage3.json", "w"), indent=1, default=float)
print(f"{RUN}: {len(out)} bloques, g_b {min(g):.4f}-{max(g):.4f} (rango {max(g) - min(g):.4f}), σ(g_b) mediana {st.median(gerr):.4f}, canales con G por bloque {[r['n_ch_fit'] for r in out]}")
for k in ("q_A_T5", "q_B_T5", "q_A_T10", "q_B_T10", "eXgY_B_T5", "eYgX_B_T5"):
    f = fits[k]
    print(f"  {k}: s = {f['s']:+.3f} ± {f['s_err']:.3f} ({f['significance']:.1f}σ) χ²/ndf {f['chi2']:.1f}/{f['ndf']} r {f['pearson']:+.2f}"
          + (f"  cota95 |s| < {f['upper_95_abs_s']:.2f}" if "upper_95_abs_s" in f else "") + (f"  dε/dV {f['deps_dV_per_V']:+.3f} ± {f['deps_dV_err']:.3f} /V" if "deps_dV_per_V" in f else ""))
for k, v in ctrl.items():
    print(f"  control {k}: r(q_B5) {v['pearson_qB5']:+.2f}  r(q_A5) {v['pearson_qA5']:+.2f}  r(g) {v['pearson_with_g']:+.2f}  rango {v['range']}")
