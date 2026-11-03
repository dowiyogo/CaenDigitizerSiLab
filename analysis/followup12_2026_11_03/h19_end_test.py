#!/usr/bin/python3 -I
# h19_end_test.py — seguimiento 2026-11-03 (12b): ¿la hipótesis «hod2019 solo detecta a menos de z del extremo con MPPC» reproduce r_X, r_Y, el
# déficit de coincidencias X·Y y la cuenta de dos estaciones de run B (mayo v8 completo, 143 814 triggers, como h19_xy.py)?
# MC: h19_end_mc.npz de h19_end_mc.py (1e8 rayos, 2 639 979 triggers). Promedio de los dos órdenes de bloque (X arriba, Y arriba), como
# sensitive_length.py, y las asignaciones de extremo A+/A− de cada coordenada (rango).
#  1. P_X(z), P_Y(z) para z = 2, 3, 5, 8 mm; z que reproduce r_X = 0.0144 y r_Y = 0.0093 (35 mV) y 0.088 / 0.076 (3 mV) (interpolación lineal de la CDF
#     en la malla de 0.25 mm).
#  2. Para esas trayectorias: P(Y a menos de z_Y | X a menos de z_X) y P(cuadrado de hod2018 | ·); predichos N·P_XY(z_X, z_Y) contra 3 X·Y observadas
#     (35 mV) y N·P_XY∧18 contra 1 de dos estaciones (cota geométrica: eficiencia de hod2018 = 1; con ε_XY de hod2018 se reduce).
#  3. Datos (work/may8_ev.bin = hits de PROC, v8 completo): impactos en tiempo > 35 mV por canal de hod2019 (kCfg_2026_05: X = ch48-63, Y = ch32-47;
#     tira k por centro; k par = capa A = Delantera (X) / Opuesta (Y), impar = Trasera / Frontal); fuera de tiempo (aC > 35 mV) como control.
#     Comparación con la distribución por barra del MC a la z ajustada, para las 4 combinaciones de orientación (sentido del orden de tiras s = ±1,
#     cara Delantera/Opuesta en + o −) — ninguna está registrada. Las X·Y > 35 mV y el evento de dos estaciones (msxy_may8.bin, bits del manuscrito):
#     canales, amplitudes, kS; para cada X·Y, el rango (1 = el más cercano) del canal Y encendido entre las 16 barras Y por distancia entre el extremo
#     con MPPC de la barra X y el de la barra Y, en las 16 combinaciones de orientación.
# Salida: <out>/h19_end_test.json, h19_end_P.csv, h19_end_bars.csv, h19_end_events.csv. Uso: h19_end_test.py <mc work dir> <ev.bin> <msxy.bin> <out>
import csv, itertools, json, math, os, sys

import numpy as np

MCW, EVB, MSXY, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
HL = 37.5
EV_DT = np.dtype([("t", "<f8"), ("aS", "<f4", 64), ("kS", "<u2", 64), ("aC", "<f4", 64)])
H19X = [53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56]  # kCfg_2026_05, tira k -> canal
H19Y = [37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40]
MEAS = {"35mV": {"rX": 0.0144, "rY": 0.0093, "N_XY": 3, "N_two_station": 1}, "3mV": {"rX": 0.088, "rY": 0.076}}
ZS = (2.0, 3.0, 5.0, 8.0)

z = np.load(f"{MCW}/h19_end_mc.npz")
meta = json.load(open(f"{MCW}/h19_end_mc_meta.json"))
TRIG = meta["triggers"]
E = z["edges"]
UP = E[1:-1]  # bordes superiores de los bins regulares (0.25 … 75)
ZB = z["zb"]


def cdf(top, cx, cy, s="all"):
    h = z[f"{top}|{cx}|{cy}|{s}"][:-1, :-1]  # sin el bin «sin cruce»
    full = z[f"{top}|{cx}|{cy}|{s}"]
    PX = np.cumsum(full[:-1, :].sum(1)) / TRIG
    PY = np.cumsum(full[:, :-1].sum(0)) / TRIG
    C = np.cumsum(np.cumsum(h, 0), 1) / TRIG
    return PX, PY, C


def at(P, v):
    return float(np.interp(v, UP, P))


def z_for(P, r):
    return float(np.interp(r, P, UP)) if r <= P[-1] else None


def c2(C, zx, zy):  # bilineal en la malla de bordes superiores
    i = np.interp(zx, UP, np.arange(len(UP)))
    j = np.interp(zy, UP, np.arange(len(UP)))
    i0, j0 = int(i), int(j)
    i1, j1 = min(i0 + 1, len(UP) - 1), min(j0 + 1, len(UP) - 1)
    fi, fj = i - i0, j - j0
    return float((1 - fi) * (1 - fj) * C[i0, j0] + fi * (1 - fj) * C[i1, j0] + (1 - fi) * fj * C[i0, j1] + fi * fj * C[i1, j1])


CONV = [(cx, cy) for cx in ("A+", "A-") for cy in ("A+", "A-")]
res = {"mc": meta, "measured": MEAS, "per_z": {}, "matched": {}}
rows_p = []
for zz in ZS:
    d = {}
    for cx, cy in CONV:
        v = {}
        for top in "XY":
            PX, PY, C = cdf(top, cx, cy)
            PX18, PY18, C18 = cdf(top, cx, cy, "h18")
            pxy = c2(C, zz, zz)
            v[top] = {"PX": at(PX, zz), "PY": at(PY, zz), "PXY": pxy, "P_Y_given_X": pxy / at(PX, zz) if at(PX, zz) > 0 else None,
                      "P_X_given_Y": pxy / at(PY, zz) if at(PY, zz) > 0 else None, "P_h18_given_X": at(PX18, zz) / at(PX, zz),
                      "P_h18_given_Y": at(PY18, zz) / at(PY, zz), "PXY_h18": c2(C18, zz, zz)}
        avg = {k: 0.5 * (v["X"][k] + v["Y"][k]) if v["X"][k] is not None and v["Y"][k] is not None else None for k in v["X"]}
        d[f"{cx}/{cy}"] = {"Xtop": v["X"], "Ytop": v["Y"], "avg": avg}
        rows_p.append([zz, cx, cy] + [f"{avg[k]:.6g}" if avg[k] is not None else "" for k in avg])
    res["per_z"][str(zz)] = d
with open(f"{OUT}/h19_end_P.csv", "w") as f:
    f.write(f"# h19_end_test.py: bench MC 1e8 rays ({TRIG} triggers), hod2019 at z = 98.4 mm, efficiency 1 within z of the MPPC end, average of the two "
            "block orders; conv = MPPC end of layer A (even strips; X Delantera, Y Opuesta) at +37.5 (A+) or -37.5 (A-), layer B opposite\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["z_mm", "conv_X", "conv_Y", "PX", "PY", "PXY", "P_Y_given_X", "P_X_given_Y", "P_h18_given_X", "P_h18_given_Y", "PXY_h18"])
    w.writerows(rows_p)

N = 143814
for thr, m in MEAS.items():
    out = {}
    for cx, cy in CONV:
        v = {}
        for top in "XY":
            PX, PY, C = cdf(top, cx, cy)
            _, _, C18 = cdf(top, cx, cy, "h18")
            zx, zy = z_for(PX, m["rX"]), z_for(PY, m["rY"])
            pxy = c2(C, zx, zy)
            v[top] = {"z_X": zx, "z_Y": zy, "PXY": pxy, "N_XY_pred": N * pxy, "P_Y_given_X": pxy / m["rX"], "PXY_h18": c2(C18, zx, zy),
                      "N_two_station_pred_geom": N * c2(C18, zx, zy)}
        avg = {k: 0.5 * (v["X"][k] + v["Y"][k]) for k in v["X"]}
        out[f"{cx}/{cy}"] = {"Xtop": v["X"], "Ytop": v["Y"], "avg": avg}
    res["matched"][thr] = out

# ---------------------------------------------------------------- datos
e = np.fromfile(EVB, dtype=EV_DT)
e = e[np.argsort(e["t"], kind="stable")]
assert len(e) == N, len(e)
T = 35.0
hitS, hitC = e["aS"][:, 32:] > T, e["aC"][:, 32:] > T
cnt = {c: int(hitS[:, c - 32].sum()) for c in range(32, 64)}
cntC = {c: int(hitC[:, c - 32].sum()) for c in range(32, 64)}
NX = int(hitS[:, [c - 32 for c in H19X]].any(1).sum())
NY = int(hitS[:, [c - 32 for c in H19Y]].any(1).sum())
res["data"] = {"N": N, "N_X_35mV": NX, "N_Y_35mV": NY, "rX": NX / N, "rY": NY / N, "N_out_of_time_any": int(hitC.any(1).sum()),
               "per_channel": {str(c): {"in_time": cnt[c], "out_of_time": cntC[c]} for c in range(32, 64)}}

# distribución por barra: datos contra MC a la z ajustada (35 mV)
bars = []
for coord, chans in (("X", H19X), ("Y", H19Y)):
    dk = np.array([cnt[chans[k]] for k in range(16)], float)
    for s in (+1, -1):
        for face_plus in (True, False):  # Delantera (X) / Opuesta (Y) con MPPC en +37.5
            pred = np.zeros(16)
            for top in "XY":
                key = "rX" if coord == "X" else "rY"
                # z ajustada con la convención de MC correspondiente (índice de MC j = k si s = +1, 15 − k si s = −1)
                for k in range(16):
                    j = k if s > 0 else 15 - k
                    data_A = k % 2 == 0
                    mc_A = j % 2 == 0
                    plus = face_plus if data_A else not face_plus  # extremo de la barra k
                    cv = "A+" if (plus == mc_A) else "A-"
                    PX, PY, _ = cdf(top, cv, "A+") if coord == "X" else cdf(top, "A+", cv)
                    zc = z_for(PX if coord == "X" else PY, MEAS["35mV"][key])
                    pb = z[f"bars|{top}|{coord}|{cv}"]
                    pred[k] += 0.5 * np.array([np.interp(zc, ZB, pb[:, jj]) for jj in range(16)])[j]
            pf = pred / pred.sum()
            expn = pf * dk.sum()
            chi2 = float((((dk - expn) ** 2) / np.maximum(expn, 1e-9)).sum())
            bars.append({"coord": coord, "strip_order_sign": s, "face_A_at_plus": face_plus, "chi2": chi2, "ndf": 15,
                         "data": dk.tolist(), "pred_frac": pf.tolist()})
res["bars"] = bars
with open(f"{OUT}/h19_end_bars.csv", "w") as f:
    f.write("# in-time hits > 35 mV per strip (May v8 complete) and MC fraction per bar at the z that reproduces r_X / r_Y, for each orientation; "
            "strip k by centre, kCfg_2026_05 channel; face: X even k = Delantera, odd = Trasera; Y even = Opuesta, odd = Frontal\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["coord", "strip_k", "channel", "face", "data_in_time", "data_out_of_time"] +
               [f"pred_frac_s{b['strip_order_sign']:+d}_{'Aplus' if b['face_A_at_plus'] else 'Aminus'}" for b in bars if b["coord"] == "X"])
    for coord, chans in (("X", H19X), ("Y", H19Y)):
        bb = [b for b in bars if b["coord"] == coord]
        for k in range(16):
            face = (("Delantera" if k % 2 == 0 else "Trasera") if coord == "X" else ("Opuesta" if k % 2 == 0 else "Frontal"))
            w.writerow([coord, k, chans[k], face, cnt[chans[k]], cntC[chans[k]]] + [f"{b['pred_frac'][k]:.5f}" for b in bb])

# eventos X·Y y dos estaciones
xev = hitS[:, [c - 32 for c in H19X]].any(1) & hitS[:, [c - 32 for c in H19Y]].any(1)
ms = np.fromfile(MSXY, dtype=np.dtype([("t", "<f8"), ("b", "u1")]))
bit = lambda k: (ms["b"] >> k) & 1 == 1  # noqa: E731
two = bit(0) & bit(1) & bit(4) & bit(5)
t_two = ms["t"][two] * 1e-6 + 1778886208.8
ev_rows = []
idx_two = [int(np.argmin(np.abs(e["t"] - t))) for t in t_two]
res["two_station_msxy"] = {"n": int(two.sum()), "t_utc_s": t_two.tolist(), "dt_to_nearest_ev_s": [float(e["t"][i] - t) for i, t in zip(idx_two, t_two)]}


def rank_table(kx, ky_list):
    out = []
    for sx, fx, sy, fy in itertools.product((1, -1), (True, False), (1, -1), (True, False)):
        xc = sx * (-15 + 2 * kx)
        yX = HL if ((kx % 2 == 0) == fx) else -HL
        dist = []
        for k2 in range(16):
            yc = sy * (-15 + 2 * k2)
            xY = HL if ((k2 % 2 == 0) == fy) else -HL
            dist.append(math.hypot(xc - xY, yX - yc))
        order = np.argsort(dist, kind="stable")
        rk = {int(k2): int(np.where(order == k2)[0][0]) + 1 for k2 in range(16)}
        out.append({"sx": sx, "X_face_A_plus": fx, "sy": sy, "Y_face_A_plus": fy, "ranks": [rk[k2] for k2 in ky_list],
                    "nearest_Y_strip": int(order[0]), "dist_hit_mm": [round(dist[k2], 2) for k2 in ky_list], "dist_min_mm": round(min(dist), 2)})
    return out


res["xy_events"] = []
for i in list(np.where(xev)[0]) + [j for j in idx_two if not xev[j]]:
    xs = [k for k in range(16) if hitS[i, H19X[k] - 32]]
    ys = [k for k in range(16) if hitS[i, H19Y[k] - 32]]
    d = {"index": int(i), "t_utc_s": float(e["t"][i]), "two_station": i in idx_two, "is_h19_xy_35mV": bool(xev[i]),
         "X": [{"k": k, "ch": H19X[k], "a_mV": float(e["aS"][i, H19X[k]]), "kS": int(e["kS"][i, H19X[k]])} for k in xs],
         "Y": [{"k": k, "ch": H19Y[k], "a_mV": float(e["aS"][i, H19Y[k]]), "kS": int(e["kS"][i, H19Y[k]])} for k in ys],
         "h19_all_gt3mV": [int(c) for c in range(32, 64) if e["aS"][i, c] > 3.0]}
    if xs and ys:
        d["rank"] = {f"X{kx}": rank_table(kx, ys) for kx in xs}
    res["xy_events"].append(d)
    for x in d["X"] + d["Y"]:
        ev_rows.append([i, f"{d['t_utc_s']:.3f}", int(d["two_station"]), "X" if x in d["X"] else "Y", x["k"], x["ch"], f"{x['a_mV']:.2f}", x["kS"]])
with open(f"{OUT}/h19_end_events.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["event_index", "t_utc_s", "two_station", "coord", "strip_k", "channel", "aS_mV", "kS"])
    w.writerows(ev_rows)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, np.generic):
        return o.item()
    return o


json.dump(clean(res), open(f"{OUT}/h19_end_test.json", "w"), indent=1)

# ---------------------------------------------------------------- resumen
print(f"MC: {TRIG} triggers; datos: N {N}, N_X {NX} (rX {NX / N:.5f}), N_Y {NY} (rY {NY / N:.5f}), triggers con aC > 35 mV: {res['data']['N_out_of_time_any']}")
for zz in ZS:
    for c, v in res["per_z"][str(zz)].items():
        a = v["avg"]
        print(f"z {zz:>3g} mm {c}: P_X {a['PX']:.5f} P_Y {a['PY']:.5f} P_XY {a['PXY']:.2e} P(Y|X) {a['P_Y_given_X']:.2e} P(18|X) {a['P_h18_given_X']:.4f} "
              f"P(18|Y) {a['P_h18_given_Y']:.4f} P_XY∧18 {a['PXY_h18']:.2e}")
for thr, o in res["matched"].items():
    for c, v in o.items():
        a, x, y = v["avg"], v["Xtop"], v["Ytop"]
        print(f"{thr} {c}: z_X {a['z_X']:.2f} (X/Y arriba {x['z_X']:.2f}/{y['z_X']:.2f}) z_Y {a['z_Y']:.2f} ({x['z_Y']:.2f}/{y['z_Y']:.2f}) mm; "
              f"N_XY pred {a['N_XY_pred']:.3f} ({x['N_XY_pred']:.3f}/{y['N_XY_pred']:.3f}); P(Y|X) {a['P_Y_given_X']:.2e}; "
              f"dos estaciones pred (geom.) {a['N_two_station_pred_geom']:.4f}")
print("por canal > 35 mV (X k: canal n):", " ".join(f"{k}:{H19X[k]}:{cnt[H19X[k]]}" for k in range(16)))
print("por canal > 35 mV (Y k: canal n):", " ".join(f"{k}:{H19Y[k]}:{cnt[H19Y[k]]}" for k in range(16)))
for b in bars:
    print(f"barras {b['coord']} s {b['strip_order_sign']:+d} cara A en {'+' if b['face_A_at_plus'] else '-'}: χ² {b['chi2']:.1f}/15")
print("dos estaciones (msxy):", res["two_station_msxy"])
for d in res["xy_events"]:
    print(f"evento {d['index']} t {d['t_utc_s']:.3f} dos estaciones {d['two_station']} X·Y35 {d['is_h19_xy_35mV']}: X {[(x['k'], x['ch'], round(x['a_mV'], 1), x['kS']) for x in d['X']]} "
          f"Y {[(y['k'], y['ch'], round(y['a_mV'], 1), y['kS']) for y in d['Y']]}; hod2019 > 3 mV: {d['h19_all_gt3mV']}")
    for kx, tab in d.get("rank", {}).items():
        rks = [t["ranks"] for t in tab]
        n1 = sum(1 for r in rks if 1 in r)
        print(f"   {kx}: rango del Y encendido en las 16 orientaciones: {rks}; rango 1 en {n1}/16")
