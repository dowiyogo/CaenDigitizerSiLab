#!/usr/bin/python3 -I
# stage5.py — PROC-1 etapa 5 (solo lectura; A8): fracción de doble impacto de hod2018 contra el umbral, medida y del MC.
#  - Definición: f2 = pares de tiras vecinas (coordenadas k, k+1, que se traslapan 1 mm) con impacto verdadero por trigger
#    con impacto verdadero en el plano sonda. Pares verdaderos = pares en el mismo trigger − pares de eventos mezclados (tira k
#    del trigger i con k+1 del i+1 y al revés, × N/(N−1)): restan ruido·ruido, muón·ruido y ruido·muón, como el estimador X·Y
#    (xy_estimator.py). Triggers con impacto verdadero en la sonda = N (P_s − P_c)/(1 − P_c) (P_s, P_c: algún impacto en la
#    ventana signal y fuera de tiempo). Etiqueta: el otro plano con impacto en tiempo sobre 10 p.e. (misma muestra para todos
#    los T). Syst: diferencia al subir la etiqueta a 15 p.e.; la versión por clases (1 tira o 2 vecinas, restando la ventana fuera
#    de tiempo) se guarda solo como referencia (sobrecorrige a T bajo).
#  - Datos: mayo v8 tiempo bueno (orden de la etapa 1, G reajustada) y octubre run 8 (r4, c1); ch3 (coordenada X 12) muerto.
#  - MC: bench_mc.py, cos³θ, placas OLD, triggers T1 ∧ T2; tira con impacto = camino > corte (0, 0.2, 0.5 mm); la coordenada
#    X 12 se apaga como en los datos; etiqueta = la otra coordenada cruzada (camino > corte); bloque X arriba y bloque Y
#    arriba. Comprobación a incidencia normal sobre |x|, |y| ≤ 16.5 mm: 15 traslapes de 1 mm / 33 mm = 0.455.
# Salida: work/stage5.json, work/fig_stage5.json (proc_plot_stage5.C).
# Uso: stage5.py <dir del paquete> [N del MC]
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import bench_mc as bm  # noqa: E402
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
NMC = int(sys.argv[2]) if len(sys.argv) > 2 else 30_000_000
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": 5, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
TS = list(range(3, 13))
DEAD_X = 12  # coordenada de ch3 en el orden de mayo y en r4
st3 = json.load(open(W + "/stage3.json"))["strip_orders"]
ORDER = {"may": st3["may"]["order"], "oct": st3["oct"]["order"]}


def classes(H):
    """(N con exactamente 1 tira, N con exactamente 2 tiras vecinas) por evento de H (n, 16)"""
    n = H.sum(1)
    adj = (H[:, :-1] & H[:, 1:]).any(1)
    return int((n == 1).sum()), int(((n == 2) & adj).sum())


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
G = {"may": fill(g_refit()), "oct": fill(mc.g_c1())}
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
DATA = {"may8_good": (e8[mc.good_time(e8["t"], "v8")], "may", "mayo v8, tiempo bueno"),
        "oct8": (xy.load_hits(W + "/hd64_oct8_hits.bin", 64), "oct", "octubre run 8")}


def hits(e, chans, T, g, field):
    chans = np.asarray(chans)
    h = e[field][:, chans] > T * g[chans]
    h[:, chans == 3] = False
    return h


S["measured"] = {}
for key, (e, cal, lab) in DATA.items():
    o, g = ORDER[cal], G[cal]
    d = {"label": lab, "events": int(len(e))}
    for probe, tagp in (("X", "Y"), ("Y", "X")):
        def f2_at(T, tagT):
            tag = hits(e, o[f"h18_{tagp}"], tagT, g, "aS").any(1)
            hs = hits(e[tag], o[f"h18_{probe}"], T, g, "aS")
            hc = hits(e[tag], o[f"h18_{probe}"], T, g, "aC")
            N = len(hs)
            pin = int((hs[:, :-1] & hs[:, 1:]).sum())
            pmix = 0.5 * (int((hs[:-1, :-1] & hs[1:, 1:]).sum()) + int((hs[1:, :-1] & hs[:-1, 1:]).sum())) * N / (N - 1)
            Ps, Pc = hs.any(1).mean(), hc.any(1).mean()
            nt = N * (Ps - Pc) / (1 - Pc)
            return {"tag_events": N, "pairs_same": pin, "pairs_mixed": pmix, "probe_true": nt, "P_s": float(Ps), "P_c": float(Pc),
                    "f2": (pin - pmix) / nt, "f2_stat": math.sqrt(pin + (pmix * (N - 1) / N) / 2) / nt}, hs, hc
        rows = []
        for T in TS:
            r, hs, hc = f2_at(T, 10)
            r15, _, _ = f2_at(T, 15)
            s1, s2 = classes(hs)
            c1, c2 = classes(hc)
            n1, n2 = s1 - c1, s2 - c2
            # versión por clases (1 tira o 2 vecinas, restando la ventana fuera de tiempo): sobrecorrige a T bajo (patrones de ruido
            # y de muón en el mismo trigger no son excluyentes); solo como referencia
            r.update(T=T, f2_tag15=r15["f2"], syst_tag=abs(r15["f2"] - r["f2"]), single_in=s1, double_in=s2, single_oot=c1, double_oot=c2,
                     f2_classes=n2 / (n1 + n2) if n1 + n2 > 0 else None, f2_raw_classes=s2 / (s1 + s2))
            rows.append(r)
        d[probe] = rows
        print(key, probe, " ".join("%d:%.3f±%.3f±%.3f" % (r["T"], r["f2"], r["f2_stat"], r["syst_tag"]) for r in rows))
    S["measured"][key] = d

# ---------------------------------------------------------------- MC
rng = np.random.default_rng(42)
CUTS = (0.0, 0.2, 0.5)
cnt = {(top, cut, p): [0, 0, 0, 0] for top in ("X", "Y") for cut in CUTS for p in "XY"}
trig = 0
for start in range(0, NMC, 2_000_000):
    m = min(2_000_000, NMC - start)
    r = bm.generate(rng, m, "cos3")
    t = bm.trigger(r, bm.GEOMS["OLD"])
    rt = bm.Rays(r.x0[t], r.y0[t], r.dx[t], r.dy[t], r.dz[t])
    trig += int(t.sum())
    for top in ("X", "Y"):
        p = bm.bar_paths(rt, 0.0, top)
        for cut in CUTS:
            H = {c: p[c] > cut for c in "XY"}
            H["X"][:, DEAD_X] = False
            for probe, tagp in (("X", "Y"), ("Y", "X")):
                sel = H[tagp].any(1)
                hp = H[probe][sel]
                a, b = classes(hp)
                c = cnt[(top, cut, probe)]
                c[0] += a
                c[1] += b
                c[2] += int((hp[:, :-1] & hp[:, 1:]).sum())
                c[3] += int(hp.any(1).sum())
MCR = {"rays": NMC, "triggers": trig, "model": "cos3", "geometry": "OLD", "seed": 42, "values": []}
for (top, cut, probe), (a, b, pr, nh) in cnt.items():
    f = pr / nh
    MCR["values"].append({"top_block": top, "path_cut_mm": cut, "probe": probe, "single": a, "double": b, "pairs": pr, "probe_events": nh, "f2": f,
                          "f2_err": math.sqrt(f * (1 - f) / nh), "f2_classes": b / (a + b)})
    print("MC", top, "arriba, corte", cut, probe, "f2 %.4f clases %.4f" % (f, b / (a + b)))
# incidencia normal, sin tira muerta y sin corte
n = 4_000_000
u = rng.uniform(-16.5, 16.5, n)
rv = bm.Rays(u, u * 0.0, np.zeros(n), np.zeros(n), -np.ones(n))
pv = bm.bar_paths(rv, 0.0, "X")["X"] > 0
kv = pv.sum(1)
MCR["normal_incidence"] = {"rays": n, "f_two_bars": float((kv == 2).mean()), "f_one_bar": float((kv == 1).mean()), "expected": 15 / 33,
                           "f2_definition": float((kv == 2).sum() / ((kv == 1) | (kv == 2)).sum())}
print("incidencia normal:", MCR["normal_incidence"])
S["mc"] = MCR

# ---------------------------------------------------------------- figura
fig = {"T": TS, "series": [], "mc": []}
for key, lab in (("may8_good", "May v8"), ("oct8", "Oct run 8")):
    for probe in "XY":
        rr = S["measured"][key][probe]
        fig["series"].append({"name": f"{lab}, {probe}", "f2": [r["f2"] for r in rr], "stat": [r["f2_stat"] for r in rr],
                              "syst": [r["syst_tag"] for r in rr]})
for cut in CUTS:
    v = [x for x in MCR["values"] if x["path_cut_mm"] == cut]
    fig["mc"].append({"cut_mm": cut, "min": min(x["f2"] for x in v), "max": max(x["f2"] for x in v)})
fig["normal"] = MCR["normal_incidence"]["f2_definition"]
json.dump(fig, open(W + "/fig_stage5.json", "w"), indent=1)
json.dump(S, open(W + "/stage5.json", "w"), indent=1, default=float)
print("hecho")
