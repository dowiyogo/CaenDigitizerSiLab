#!/usr/bin/python3 -I
# recabling_groups.py — seguimiento 2026-10-31 (7), punto 1: ¿la pérdida de luz de mayo a octubre y el cambio de forma del pulso de 1 p.e.
# siguen al recableado del plano Y de hod2018?
# Emparejamiento mayo ↔ octubre:
#   - por MPPC (barra): X igual canal (no se recableó); Y la tira i del orden de mayo (stage3.json, validado con datos de mayo) con la tira i del
#     orden de octubre (r4). Comprobación: el grupo HV medido en mayo (hv_group_previous del mapeo r4, apagones de mayo) del canal de mayo debe ser
#     el grupo HV medido en octubre (r4) del canal de octubre de la misma barra (la línea de HV va con el MPPC, no con el cable de señal);
#     se prueban las dos orientaciones relativas (i ↔ i, i ↔ 15 − i).
#   - por canal del digitalizador (como pe_primary_channels.csv): mismo canal en las dos campañas.
# Razón mayo/octubre en p.e. primarios (pe_primary_channels.csv del seguimiento 6) por MPPC, con el error de las dos medianas.
# Grupos (atributos del MPPC; los del cable se dan para el canal de octubre y, en Y, también para el de mayo): plano, capa (A = índice par de
# la coordenada, B = impar), grupo HV (r4), grupo de 8 entradas del VX1740 (canal // 8) y posición en el grupo (canal % 8) [supuesto: la
# "posición en el conector" es la posición dentro del grupo de 8 entradas].
# Por grupo: media ponderada (1/σ²), su error, desviación estándar ponderada, χ²/ndf y error escalado por √(χ²/ndf) si > 1; diferencia X − Y con
# su significancia (error escalado).
# Exclusión principal: los pares con ch29 o ch30 en cualquiera de las dos campañas (por MPPC: barras Y11, Y13, Y14); variante: solo la barra
# con la anomalía en las dos campañas (Y13 = ch30 de mayo = ch29 de octubre).
# Forma: T(d), d = 3..8, de pair_shapes.C (pulsos oscuros aislados de 0.5-1.5 p.e., sin otro pulso en ±12 muestras; mayo v8 pre-trigger [20,258),
# octubre run 8 post-señal [58,144) y pre-trigger [20,49)). ΔT(d) = T_oct − T_mayo por MPPC y por canal; agrupado igual. Error de T: s_d/√n
# con s_d (dispersión de un pulso) estimado de la diferencia entre las regiones pre y post de octubre del mismo canal [inferido: supone la
# misma forma en las dos regiones].
# Salida: <out>/recabling_groups.json, ratio_by_mppc.csv, shape_by_mppc.csv. Uso: recabling_groups.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
ST3 = json.load(open(REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work/stage3.json"))["strip_orders"]
MAP = {c["dig_ch"]: c for c in json.load(open(REPO + "/daq/config/mapping/hodo_mapping_2026_10_r4.json"))["channels"]}
PE = {int(r["ch"]): r for r in csv.DictReader(open(REPO + "/outputs/hod2018_mip/2026-10-31_pe_primary/pe_primary_channels.csv"))}
SH = REPO + "/outputs/hod_dark/2026-10-31_shapes31"
DS = range(3, 9)
res = {}


def tmpl(f):
    T, N = {}, {}
    for r in csv.DictReader(open(f"{SH}/{f}_template.csv")):
        T.setdefault(int(r["ch"]), {})[int(r["d"])] = float(r["T"])
        N[int(r["ch"])] = int(r["n_isolated"])
    return T, N


Tm, Nm = tmpl("may_v8good_pre")
To, No = tmpl("oct_run8_post")
Tp, Np = tmpl("oct_run8_pre")
# s_d: Var(T_pre − T_post) = s_d² (1/n_pre + 1/n_post), agrupado sobre canales
S = {}
for d in DS:
    num = [(Tp[c][d] - To[c][d]) ** 2 for c in To]
    den = [1 / Np[c] + 1 / No[c] for c in To]
    S[d] = math.sqrt(sum(num) / sum(den))
res["s_single_pulse"] = S
res["oct_pre_minus_post_mean"] = {d: float(np.mean([Tp[c][d] - To[c][d] for c in To])) for d in DS}

# ---------------------------------------------------------------- emparejamiento
MX, MY = ST3["may"]["order"]["h18_X"], ST3["may"]["order"]["h18_Y"]
OX, OY = ST3["oct"]["order"]["h18_X"], ST3["oct"]["order"]["h18_Y"]
assert MX == OX
check = {}
for name, f in (("same", lambda i: i), ("reversed", lambda i: 15 - i)):
    rows = []
    for i in range(16):
        cm, co = MY[i], OY[f(i)]
        gprev = (MAP[cm].get("hv_group_previous") or {}).get("hv_group")
        if gprev is not None:
            rows.append({"bar_may": i, "ch_may": cm, "hv_group_may": gprev, "ch_oct": co, "hv_group_oct": MAP[co]["hv_group"], "match": gprev == MAP[co]["hv_group"]})
    check[name] = {"n": len(rows), "n_match": sum(r["match"] for r in rows), "rows": rows}
res["y_pairing_hv_check"] = check
# X: el mismo control con el canal
rows = []
for c in MX:
    gprev = (MAP[c].get("hv_group_previous") or {}).get("hv_group")
    if gprev is not None and MAP[c]["hv_group"] is not None:
        rows.append({"ch": c, "hv_group_may": gprev, "hv_group_oct": MAP[c]["hv_group"], "match": gprev == MAP[c]["hv_group"]})
res["x_same_channel_hv_check"] = {"n": len(rows), "n_match": sum(r["match"] for r in rows), "rows": rows}

mppc = []
for i, c in enumerate(OX):
    if c == 3:
        continue
    mppc.append({"id": f"X{i:02d}", "plane": "X", "index": i, "layer": "A" if i % 2 == 0 else "B", "ch_may": c, "ch_oct": c})
for i in range(16):
    mppc.append({"id": f"Y{i:02d}", "plane": "Y", "index": i, "layer": "A" if i % 2 == 0 else "B", "ch_may": MY[i], "ch_oct": OY[i]})
for m in mppc:
    co, cm = m["ch_oct"], m["ch_may"]
    m.update(hv_group=MAP[co]["hv_group"], adc_group_oct=co // 8, pos_oct=co % 8, adc_group_may=cm // 8, pos_may=cm % 8)
    a, ea = float(PE[cm]["may_median_pe_prim"]), float(PE[cm]["may_median_pe_prim_err"])
    b, eb = float(PE[co]["oct_median_pe_prim"]), float(PE[co]["oct_median_pe_prim_err"])
    m.update(may_prim=a, may_prim_err=ea, oct_prim=b, oct_prim_err=eb, ratio=a / b, ratio_err=a / b * math.hypot(ea / a, eb / b))
    for d in DS:
        m[f"dT{d}"] = To[co][d] - Tm[cm][d]
        m[f"dT{d}_err"] = math.hypot(S[d] / math.sqrt(No[co]), S[d] / math.sqrt(Nm[cm]))
        m[f"Tmay{d}"], m[f"Toct{d}"] = Tm[cm][d], To[co][d]
        # por canal del digitalizador (mismo canal, en Y otro MPPC)
        m[f"dT{d}_samech"] = To[co][d] - Tm[co][d]
    m["ratio_samech"] = float(PE[co]["ratio_prim"])


def wstats(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float)
    w = 1 / e ** 2
    mu = float((w * v).sum() / w.sum())
    err = float(1 / math.sqrt(w.sum()))
    chi2 = float((w * (v - mu) ** 2).sum())
    ndf = len(v) - 1
    sd = float(math.sqrt((w * (v - mu) ** 2).sum() / w.sum())) if len(v) > 1 else 0.0
    sc = err * math.sqrt(max(1.0, chi2 / ndf)) if ndf > 0 else err
    return {"n": len(v), "mean": mu, "err": err, "err_scaled": sc, "wsd": sd, "chi2": chi2, "ndf": ndf, "unweighted_median": float(np.median(v))}


EXCL = {"main": lambda m: m["ch_may"] in (29, 30) or m["ch_oct"] in (29, 30), "only_Y13": lambda m: m["id"] == "Y13"}
KEYS = ["plane", "layer", "hv_group", "adc_group_oct", "pos_oct", "adc_group_may", "pos_may"]


def grouped(sel, val, err):
    out = {}
    for k in KEYS:
        g = {}
        for m in sel:
            g.setdefault(str(m[k]), []).append(m)
        out[k] = {kk: wstats([val(x) for x in v], [err(x) for x in v]) for kk, v in sorted(g.items())}
    X = [m for m in sel if m["plane"] == "X"]
    Y = [m for m in sel if m["plane"] == "Y"]
    sx, sy = wstats([val(x) for x in X], [err(x) for x in X]), wstats([val(x) for x in Y], [err(x) for x in Y])
    dif = sx["mean"] - sy["mean"]
    e = math.hypot(sx["err_scaled"], sy["err_scaled"])
    e0 = math.hypot(sx["err"], sy["err"])
    out["X_minus_Y"] = {"diff": dif, "err_scaled": e, "z_scaled": dif / e, "err_stat": e0, "z_stat": dif / e0}
    out["all"] = wstats([val(x) for x in sel], [err(x) for x in sel])
    return out


res["ratio"], res["shape"] = {}, {}
for ex, f in EXCL.items():
    sel = [m for m in mppc if not f(m)]
    res["ratio"][ex] = {"excluded": [m["id"] for m in mppc if f(m)], "by_mppc": grouped(sel, lambda m: m["ratio"], lambda m: m["ratio_err"])}
    res["shape"][ex] = {f"dT{d}": grouped(sel, lambda m, d=d: m[f"dT{d}"], lambda m, d=d: m[f"dT{d}_err"]) for d in DS}
# por canal del digitalizador (ratio_prim del seguimiento 6; excluye ch29 y ch30)
selc = [m for m in mppc if m["ch_oct"] not in (29, 30)]
res["ratio"]["same_channel_excl_ch29_ch30"] = {
    "X_mean": wstats([float(PE[m["ch_oct"]]["ratio_prim"]) for m in selc if m["plane"] == "X"],
                     [m["ratio_err"] for m in selc if m["plane"] == "X"]),
    "Y_mean_same_dig_channel": wstats([m["ratio_samech"] for m in selc if m["plane"] == "Y" and m["ch_oct"] not in (29, 30)],
                                      [float(PE[m["ch_oct"]]["ratio_prim"]) * math.hypot(float(PE[m["ch_oct"]]["may_median_pe_prim_err"]) / float(PE[m["ch_oct"]]["may_median_pe_prim"]),
                                                                                         float(PE[m["ch_oct"]]["oct_median_pe_prim_err"]) / float(PE[m["ch_oct"]]["oct_median_pe_prim"]))
                                       for m in selc if m["plane"] == "Y"])}
# ¿la forma sigue al MPPC o al canal? dispersión de ΔT(d) en Y con los dos emparejamientos
res["shape_pairing_test_Y"] = {}
for d in DS:
    Y = [m for m in mppc if m["plane"] == "Y" and not EXCL["main"](m)]
    a = [m[f"dT{d}"] for m in Y]
    b = [m[f"dT{d}_samech"] for m in Y]
    res["shape_pairing_test_Y"][d] = {"sd_by_mppc": float(np.std(a)), "sd_by_channel": float(np.std(b)), "mean_by_mppc": float(np.mean(a)), "mean_by_channel": float(np.mean(b))}
json.dump(res, open(OUT + "/recabling_groups.json", "w"), indent=1, default=float)
keys = list(mppc[0])
for fn, ks in (("ratio_by_mppc.csv", [k for k in keys if not k.startswith(("dT", "Tmay", "Toct"))]), ("shape_by_mppc.csv", [k for k in keys if not k.startswith(("may_prim", "oct_prim", "ratio"))])):
    with open(f"{OUT}/{fn}", "w") as f:
        f.write(",".join(ks) + "\n")
        for m in mppc:
            f.write(",".join(f"{m[k]:.5f}" if isinstance(m[k], float) else str(m[k]) for k in ks) + "\n")

# ---------------------------------------------------------------- salida textual
print("comprobación del emparejamiento Y por grupo HV:", {k: f"{v['n_match']}/{v['n']}" for k, v in check.items()},
      "| X mismo canal:", f"{res['x_same_channel_hv_check']['n_match']}/{res['x_same_channel_hv_check']['n']}")
for r in check["same"]["rows"]:
    print("   ", r)
print("s_d (un pulso):", {d: round(v, 3) for d, v in S.items()}, "| octubre pre − post medio:", {d: round(v, 4) for d, v in res["oct_pre_minus_post_mean"].items()})
for m in mppc:
    print(f"{m['id']} ch mayo {m['ch_may']:2d} oct {m['ch_oct']:2d} HV {m['hv_group']} capa {m['layer']} | razón {m['ratio']:.3f} ± {m['ratio_err']:.3f} "
          f"(mismo canal {m['ratio_samech']:.3f}) | ΔT(3..8) " + " ".join(f"{m[f'dT{d}']:+.3f}" for d in DS))
for ex in EXCL:
    print(f"== razón, exclusión {ex}: {res['ratio'][ex]['excluded']}")
    g = res["ratio"][ex]["by_mppc"]
    for k in KEYS:
        print(f"  {k}: " + "; ".join(f"{kk}: {v['mean']:.3f} ± {v['err_scaled']:.3f} (n {v['n']}, sd {v['wsd']:.3f}, χ²/ndf {v['chi2']:.1f}/{v['ndf']})" for kk, v in g[k].items()))
    print(f"  todos: {g['all']['mean']:.3f} ± {g['all']['err_scaled']:.3f} (χ²/ndf {g['all']['chi2']:.1f}/{g['all']['ndf']}); X − Y = {g['X_minus_Y']['diff']:+.3f} ± "
          f"{g['X_minus_Y']['err_scaled']:.3f} (z = {g['X_minus_Y']['z_scaled']:.2f}; solo estadístico z = {g['X_minus_Y']['z_stat']:.2f})")
    for d in DS:
        g = res["shape"][ex][f"dT{d}"]
        print(f"  ΔT({d}): " + " | ".join(f"{k}: " + ", ".join(f"{kk} {v['mean']:+.4f}±{v['err_scaled']:.4f}" for kk, v in g[k].items()) for k in ("plane", "layer", "hv_group", "adc_group_oct", "pos_oct"))
              + f" | todos {g['all']['mean']:+.4f}±{g['all']['err_scaled']:.4f} (χ²/ndf {g['all']['chi2']:.0f}/{g['all']['ndf']}) | X−Y {g['X_minus_Y']['diff']:+.4f}±{g['X_minus_Y']['err_scaled']:.4f} "
                f"(z {g['X_minus_Y']['z_scaled']:.2f})")
print("mismo canal (excl. ch29, ch30):", json.dumps({k: {kk: round(vv, 4) if isinstance(vv, float) else vv for kk, vv in v.items()} for k, v in res["ratio"]["same_channel_excl_ch29_ch30"].items()}))
print("forma en Y, ¿sigue al MPPC o al canal? (sd de ΔT):", {d: (round(v["sd_by_mppc"], 4), round(v["sd_by_channel"], 4)) for d, v in res["shape_pairing_test_Y"].items()})
