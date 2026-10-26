#!/usr/bin/python3 -I
# profile_summary.py — XPLUS-4b etapa 4b.2/4b.3: resumen del perfil de profile_4b.C por coordenada.
#  - regiones (inicio de la ventana de 9 muestras): pre lejano = [20, 200]; fuera de tiempo del estimador = w_oot; post lejano = las
#    últimas 150 muestras del registro (w desde n − 150 hasta n − 9); en SW trigger, "todo" = [20, n − 9].
#  - media de una región con ventanas que no se solapan (cada 9 muestras); error = máx(dispersión entre ventanas/√k, binomial/√k).
#  - a^oot (vector) en w_oot con error binomial; a^oot (estimador, aC) de _windows.csv; razones con propagación de errores.
#  - estructura: mínimo y máximo de la fracción en ±60 muestras alrededor de la señal (fuera de ella) contra el pre lejano.
# Uso: profile_summary.py <prefijo> <w_oot o -1 (SW)> <n del registro> <salida.json> [pre_lo,pre_hi post_lo,post_hi s0]
#      (octubre, 144 muestras: pre lejano [20,40], post lejano [110,135], s0 = 49)
import csv, json, math, sys

pre, w_oot, n, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
PRE = tuple(map(int, sys.argv[5].split(","))) if len(sys.argv) > 5 else (20, 200)
POST = tuple(map(int, sys.argv[6].split(","))) if len(sys.argv) > 6 else (n - 150, n - 9)
S0 = int(sys.argv[7]) if len(sys.argv) > 7 else 258
P = {}
for r in csv.DictReader(open(pre + "_profile.csv")):
    P.setdefault(r["coord"], {})[int(r["w_sample"])] = (int(r["n_events"]), int(r["n_hit"]), float(r["frac"]), r["overlaps_signal"] == "1")
W = {r["coord"]: r for r in csv.DictReader(open(pre + "_windows.csv"))}


def region(d, lo, hi):
    ws = [w for w in range(lo, hi + 1, 9) if w in d and not d[w][3]]
    f = [d[w][2] for w in ws]
    k = len(f)
    m = sum(f) / k
    sd = math.sqrt(sum((x - m) ** 2 for x in f) / (k - 1)) if k > 1 else 0.0
    nE = d[ws[0]][0]
    binom = math.sqrt(m * (1 - m) / nE)
    return {"mean": m, "err": max(sd, binom) / math.sqrt(k), "spread_sd": sd, "k_windows": k, "w_range": [lo, hi]}


def ratio(a, ea, b, eb):
    r = a / b
    return r, r * math.sqrt((ea / a) ** 2 + (eb / b) ** 2)


res = {}
for co, d in P.items():
    o = {}
    if w_oot < 0:
        o["all"] = region(d, 20, n - 9)
    else:
        o["pre_far"] = region(d, *PRE)
        o["post_far"] = region(d, *POST)
        nE, nh, f, _ = d[w_oot]
        o["oot_vector"] = {"w": w_oot, "frac": f, "err": math.sqrt(f * (1 - f) / nE), "n_events": nE}
        ac = float(W[co]["ctl_frac"])
        o["oot_estimator_aC"] = {"frac": ac, "err": math.sqrt(ac * (1 - ac) / int(W[co]["n_events"]))}
        o["signal_aS"] = float(W[co]["signal_frac"])
        for key in ("pre_far", "post_far"):
            r, er = ratio(f, o["oot_vector"]["err"], o[key]["mean"], o[key]["err"])
            o[f"ratio_oot_over_{key}"] = {"value": r, "err": er}
        s0 = S0
        near = {w: d[w][2] for w in range(max(20, s0 - 60), min(n - 9, s0 + 9 + 60)) if w in d and not d[w][3]}
        o["near_signal"] = {"min": min(near.values()), "w_min": min(near, key=near.get), "max": max(near.values()), "w_max": max(near, key=near.get),
                            "pre_far_mean": o["pre_far"]["mean"]}
    res[co] = o
json.dump(res, open(out, "w"), indent=1)
print(json.dumps(res, indent=1))
