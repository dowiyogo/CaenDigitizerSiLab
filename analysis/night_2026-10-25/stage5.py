#!/usr/bin/python3 -I
# stage5.py — XPLUS-4b etapa 5: mayo v1-v10 run por run, sobre los _hits.bin de las ntuplas (idénticos a los de PROC).
#  - Tiempo bueno: may_common.good_time(t, run) (runs/good_intervals_2026.csv, EXCLUDE y HV_OFF, hora local UTC−4), el mismo criterio
#    que definió el tiempo bueno de v8, aplicado a todos los runs.
#  - Orden de tiras: el de mayo de PROC (work/stage3.json, strip_orders.may.order). Validación por run con neighbours() y
#    adjacency_check() de stage1.py (copia literal; z > 5 y top-2 mutuo, 5 p.e. con la G del run): pares confirmados de los 15 por
#    coordenada; acepta X ≥ 12/15 e Y ≥ 14/15 (los valores de v8).
#  - G por run: gain_blocks.C sobre las horas de tiempo bueno (gain_fit.h, método del reajuste de PROC), semilla = G reajustada de v8;
#    canales sin ajuste aceptado = mediana del run; si el run no da ≥ 16 canales aceptados, G de v8 (se informa).
#  - q a 5 y 10 p.e. con los dos estimadores, aX, aY, ε_XY = q/(f_2018·g); regresión: q con la G de v8 en todos los runs (como stage2.py)
#    para v7-v10 contra results.json. Ocupación máx/mín (en tiempo − fuera de tiempo, stage3.py) y f_2 (stage5.py, etiqueta 10 p.e.) por run.
#  - χ² de compatibilidad entre runs validados (error estadístico), media ponderada, error escalado por √(χ²/ndf) si p < 0.05, runs a > 2σ.
# Uso: stage5.py <dir de ntuplas> <salida dir>
import csv, json, math, os, statistics as st, subprocess, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

ND, OUT = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
os.makedirs(OUT, exist_ok=True)
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
F18, GG = 0.03522, 0.980
XCH, YCH = [c for c in range(16) if c != 3], list(range(16, 32))
Z_MIN = 5.0
RUNS = ["v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"]
NT = {r: os.path.join(ND, f"may_{r}") for r in RUNS}
fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
ORDER = json.load(open(W + "/stage3.json"))["strip_orders"]["may"]["order"]
R3 = json.load(open(os.path.join(REPO, "outputs/proceeding_2026_10/results.json")))


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


G8 = fill(g_refit())


def run_gain(run, hours):
    bf, of = f"{OUT}/{run}_hours.txt", f"{OUT}/{run}_gain.csv"
    with open(bf, "w") as f:
        for h in hours:
            f.write(f"run {h} {h + 3600}\n")
    call = f'gain_blocks.C+("{NT[run]}.root", "{bf}", "{",".join(f"{x:.6f}" for x in G8)}", "{of}")'
    r = subprocess.run(["/home/xplus/night_2026-10-25/run_macro.sh", call], capture_output=True, text=True)
    if r.returncode:
        sys.exit(r.stdout + r.stderr)
    g = np.full(32, np.nan)
    for row in csv.DictReader(open(of)):
        if row["status"] in ("ok", "changed"):
            g[int(row["ch"])] = float(row["G_mV"])
    return g


# --- copia literal de stage1.py (neighbours, adjacency_check), con e y N como argumentos
def neighbours(e, chans, thr):
    N = len(e)
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
    return [[chans[i], chans[j], float(E[i, j]), float(z[i, j])] for i, j in edges]


def adjacency_check(order, edges):
    es = {(min(a, b), max(a, b)) for a, b, _, _ in edges}
    pairs = [(min(order[k], order[k + 1]), max(order[k], order[k + 1])) for k in range(len(order) - 1)]
    return {"pairs": len(pairs), "confirmed": sum(p in es for p in pairs), "missing": [list(p) for p in pairs if p not in es]}


# --- copia de stage5.py (hits, f2_at) y de stage3.py (per_strip, ratio)
def hits(e, chans, T, g, field):
    chans = np.asarray(chans)
    h = e[field][:, chans] > T * g[chans]
    h[:, chans == 3] = False
    return h


def f2_at(e, g, probe, tagp, T, tagT):
    tag = hits(e, ORDER[f"h18_{tagp}"], tagT, g, "aS").any(1)
    hs = hits(e[tag], ORDER[f"h18_{probe}"], T, g, "aS")
    hc = hits(e[tag], ORDER[f"h18_{probe}"], T, g, "aC")
    N = len(hs)
    pin = int((hs[:, :-1] & hs[:, 1:]).sum())
    pmix = 0.5 * (int((hs[:-1, :-1] & hs[1:, 1:]).sum()) + int((hs[1:, :-1] & hs[:-1, 1:]).sum())) * N / (N - 1)
    Ps, Pc = hs.any(1).mean(), hc.any(1).mean()
    nt = N * (Ps - Pc) / (1 - Pc)
    return (pin - pmix) / nt, math.sqrt(pin + (pmix * (N - 1) / N) / 2) / nt


def occupancy(e, g, coord, T):
    chans = ORDER[f"h18_{coord}"]
    hin = hits(e, chans, T, g, "aS").sum(0).astype(float)
    hoot = hits(e, chans, T, g, "aC").sum(0).astype(float)
    net, err = hin - hoot, np.sqrt(hin + hoot)
    idx = np.array([k for k in range(16) if chans[k] != 3])
    v, ev = net[idx], err[idx]
    i, j = int(np.argmax(v)), int(np.argmin(v))
    r = v[i] / v[j]
    return float(r), float(r * math.hypot(ev[i] / v[i], ev[j] / v[j]))


def q_of(e, G, T):
    G64 = np.concatenate([G, np.full(32, np.nan)])
    Xs, Xc = xy.plane_hits(e, G64, T, XCH)
    Ys, Yc = xy.plane_hits(e, G64, T, YCH)
    return xy.xy_fraction(Xs, Ys, Xc, Yc)


def chi2_sf(x, k):
    """P(χ²_k > x) por la función gamma incompleta regularizada (serie; sin scipy)."""
    a, xx = k / 2.0, x / 2.0
    if xx <= 0:
        return 1.0
    s, term, n = 1.0 / a, 1.0 / a, 1
    while abs(term) > 1e-14 * abs(s) and n < 10000:
        term *= xx / (a + n)
        s += term
        n += 1
    return 1.0 - math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * s


rows = []
for run in RUNS:
    e_all = xy.load_hits(NT[run] + "_hits.bin", 64)
    gt = mc.good_time(e_all["t"], run)
    e = e_all[gt]
    t = e["t"]
    hours = sorted(set((np.floor(t / 3600) * 3600).astype(np.int64).tolist()))
    Gr = run_gain(run, hours)
    nacc = int(np.isfinite(Gr[[c for c in range(32) if c != 3]]).sum())
    G = fill(Gr) if nacc >= 16 else G8
    dur = float(t[-1] - t[0]) if len(t) > 1 else 0.0
    r = {"run": run, "events_all": int(len(e_all)), "events_good": int(len(e)), "good_fraction": float(len(e) / len(e_all)),
         "start_utc": float(e_all["t"][0]), "end_utc": float(e_all["t"][-1]), "duration_good_h": dur / 3600,
         "G_channels_accepted": nacc, "G_source": "run" if nacc >= 16 else "v8 (el run no da 16 canales)", "G_median_mV": float(np.median(G[[c for c in range(32) if c != 3]])),
         "G_spread_mV": float(np.std(G[[c for c in range(32) if c != 3]]))}
    val = {}
    for coord, chans in (("X", list(range(16))), ("Y", list(range(16, 32)))):
        edges = neighbours(e, chans, 5 * G[chans])
        val[coord] = adjacency_check(ORDER[f"h18_{coord}"], edges)
    r["strip_order"] = {k: f"{v['confirmed']}/{v['pairs']}" for k, v in val.items()}
    r["strip_order_valid"] = val["X"]["confirmed"] >= 12 and val["Y"]["confirmed"] >= 14
    for T in (5, 10):
        q = q_of(e, G, T)
        q8 = q_of(e, G8, T)
        r[f"T{T}"] = {"N": q["N"], "N_XY": q["N_XY"], "aX": q.get("aX"), "aY": q.get("aY"), "rX": q.get("rX"), "rY": q.get("rY"),
                      "q_mix": q.get("f_mix"), "q_mix_err": q.get("f_mix_err"), "q_ctl": q.get("f_ctl"), "q_ctl_err": q.get("f_ctl_err"),
                      "eps_XY": q.get("f_mix", float("nan")) / (F18 * GG), "eps_XY_err": q.get("f_mix_err", float("nan")) / (F18 * GG),
                      "q_mix_with_G_v8": q8.get("f_mix"), "q_mix_err_with_G_v8": q8.get("f_mix_err")}
        try:
            r[f"T{T}"]["occ_X"], r[f"T{T}"]["occ_X_err"] = occupancy(e, G, "X", T)
            r[f"T{T}"]["occ_Y"], r[f"T{T}"]["occ_Y_err"] = occupancy(e, G, "Y", T)
        except (ZeroDivisionError, ValueError, FloatingPointError):
            pass
        if T == 5:
            for probe, tagp in (("X", "Y"), ("Y", "X")):
                try:
                    r[f"T{T}"][f"f2_{probe}"], r[f"T{T}"][f"f2_{probe}_err"] = f2_at(e, G, probe, tagp, T, 10)
                except ZeroDivisionError:
                    r[f"T{T}"][f"f2_{probe}"] = None
    rows.append(r)
    t5 = r["T5"]
    print(f"{run}: ev {r['events_good']}/{r['events_all']} G {r['G_source']} ({nacc} canales, mediana {r['G_median_mV']:.3f} mV) orden X {r['strip_order']['X']} Y {r['strip_order']['Y']} "
          f"q5 {t5['q_mix']:.4f} ± {t5['q_mix_err']:.4f} (G v8: {t5['q_mix_with_G_v8']:.6f}) ε {t5['eps_XY']:.3f}", flush=True)

# regresión contra results.json (G de v8 en todos los runs, como stage2.py)
reg = {}
for run in ("v7", "v8", "v9", "v10"):
    for T in (5, 10):
        ref = next(x for x in R3["results"] if x["id"] == f"q_may{run[1:]}_good_T{T}")
        got = next(r for r in rows if r["run"] == run)[f"T{T}"]["q_mix_with_G_v8"]
        reg[f"{run}_T{T}"] = {"results_json": ref["value"], "here": got, "equal": got == ref["value"]}


def combine(T, which):
    sel = [r for r in rows if r["strip_order_valid"] and r[f"T{T}"]["q_mix"] is not None and r[f"T{T}"]["q_mix_err"] and r[f"T{T}"]["q_mix_err"] > 0]
    if which == "v7-v10":
        sel = [r for r in rows if r["run"] in ("v7", "v8", "v9", "v10")]
    x = np.array([r[f"T{T}"]["eps_XY"] for r in sel])
    s = np.array([r[f"T{T}"]["eps_XY_err"] for r in sel])
    w = 1 / s ** 2
    m = float((w * x).sum() / w.sum())
    em = float(math.sqrt(1 / w.sum()))
    chi2 = float((w * (x - m) ** 2).sum())
    ndf = len(x) - 1
    p = chi2_sf(chi2, ndf) if ndf > 0 else None
    scale = math.sqrt(chi2 / ndf) if (p is not None and p < 0.05) else 1.0
    pulls = {r["run"]: float((xi - m) / si) for r, xi, si in zip(sel, x, s)}
    return {"runs": [r["run"] for r in sel], "mean": m, "err": em * scale, "err_unscaled": em, "chi2": chi2, "ndf": ndf, "p": p, "scale": scale,
            "beyond_2sigma": [k for k, v in pulls.items() if abs(v) > 2], "pulls": pulls}


comb = {f"T{T}_{w}": combine(T, w) for T in (5, 10) for w in ("validated", "v7-v10")}
# compatibilidad de q (no de ε) v7-v10 con la G de v8, para comparar con stage2.py (χ² = 1.75/3)
qq = [(r[f"T5"]["q_mix_with_G_v8"], r["T5"]["q_mix_err_with_G_v8"]) for r in rows if r["run"] in ("v7", "v8", "v9", "v10")]
x, s = np.array([a for a, _ in qq]), np.array([b for _, b in qq])
w = 1 / s ** 2
m = (w * x).sum() / w.sum()
reg["chi2_q_v7_v10_T5_G_v8"] = {"chi2": float((w * (x - m) ** 2).sum()), "ndf": 3, "stage2_results_json": 1.746021745096426}
json.dump({"runs": rows, "regression": reg, "combined": comb}, open(f"{OUT}/stage5.json", "w"), indent=1, default=float)
print(json.dumps(reg, indent=0, default=float))
for k, v in comb.items():
    print(f"{k}: ε = {v['mean']:.4f} ± {v['err']:.4f} (sin escalar {v['err_unscaled']:.4f}) χ²/ndf {v['chi2']:.2f}/{v['ndf']} p {v['p']} escala {v['scale']:.2f} "
          f"runs {v['runs']} > 2σ {v['beyond_2sigma']}")
