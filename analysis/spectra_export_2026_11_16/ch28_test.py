#!/usr/bin/python3 -I
# ch28_test.py — seguimiento 2026-11-16 (spectra-export), punto 3: qué barra y capa leía ch28 (Y09 en el orden de mayo) en mayo. No cambia el mapeo.
# Entradas: impactos de PROC de mayo v8 (tiempo bueno) y octubre run 8, G del análisis, stage3.json, layer_from_channels.json.
# Pruebas [calculado]:
#  a) perfil a lo largo de la barra de todos los canales Y de cada campaña (como r2_profile.py: etiqueta = tiras X por grupos de 4 índices, sonda
#     > 3 p.e., mediana neta, pendiente de ln(mediana) por índice X); el signo debe ser el de la capa (YA +, YB −);
#  b) coincidencias Y–Y: las dos > 10 p.e. en la ventana signal, mismo − mezclados (canal b del evento siguiente) × N/(N − 1); las capas A y B están
#     desplazadas media barra, así que una barra coincide sobre todo con sus dos vecinas de la otra capa;
#  c) fila de la matriz X × Y del canal (distribución sobre las barras X de sus coincidencias con X > 10 p.e., menos mezclados): χ² de forma contra la de
#     cada barra Y de octubre y de mayo.
# Objetivos por capa (método de layer_targets.py, modo spectra, sin bootstrap; mediana: sonda > 5 p.e. en bins de 1 p.e., error = SD de 200 réplicas de
# Poisson con la semilla de mayo 1102; espectro de ajuste: sonda > 3 p.e.) para el orden Y de mayo:
#   nominal; ch28_to_YA = ch28 sondea en YA (YB lo pierde), en la posición de Y08; swap_20_28 = ch28 en Y08 (YA) y ch20 en Y09 (YB).
#   Cada variante escribe <out>/ch28/<variante>/spectra_layers.csv y medians.json; la moda sale de layer_fit.C sobre esos espectros.
# Salida: <out>/ch28/profiles.csv, yy_coinc.csv, xy_rows.csv, <variante>/...; ch28.json. Uso: ch28_test.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1] + "/ch28"
os.makedirs(OUT, exist_ok=True)
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
ST3 = json.load(open(PK + "/stage3.json"))["strip_orders"]
ORDER = {"may": ST3["may"]["order"], "oct": ST3["oct"]["order"]}
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
SKEY = {"may": "may2026_v8_goodtime", "oct": "oct2026_run8"}
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
CHECK = (29, 54)
EDG5 = np.arange(0, 300.0001, 0.05)
EDG1 = np.arange(0, 301, 1.0)


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


G = {"may": g_refit(), "oct": mc.g_c1()}
G = {k: np.where(np.isfinite(v), v, np.nanmedian(v)) for k, v in G.items()}
e8 = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
e8 = e8[np.argsort(e8["t"], kind="stable")]
DATA = {"may": e8[mc.good_time(e8["t"], "v8")], "oct": xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)}
DATA["oct"] = DATA["oct"][np.argsort(DATA["oct"]["t"], kind="stable")]
PE = {c: DATA[c]["aS"][:, :32] / G[c][None, :32] for c in DATA}
SC = {c: len(DATA[c]) / (len(DATA[c]) - 1) for c in DATA}
rng = np.random.default_rng(1803)


def med(net, edg):
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    before = cs[b - 1] if b > 0 else 0.0
    return float(edg[b] + (0.5 * n - before) / net[b] * (edg[b + 1] - edg[b]))


res = {}
# a) perfiles
prow = []
for camp in ("may", "oct"):
    o = ORDER[camp]
    xc = np.asarray(o["h18_X"])
    t = PE[camp][:, xc] > 10.0
    t[:, xc == DEAD] = False
    for j, ch in enumerate(o["h18_Y"]):
        L = "A" if j % 2 == 0 else "B"
        pts = []
        for k in range(4):
            js = [x for x in range(4 * k, 4 * k + 4) if x != DEAD_X]
            T = t[:, js].any(1)
            Tn = np.concatenate([T[1:], [False]])
            ps = PE[camp][:, ch]
            sel = ps > 3.0
            hs = np.histogram(ps[sel & T], bins=EDG5)[0].astype(float)
            hm = np.histogram(ps[sel & Tn], bins=EDG5)[0].astype(float)
            m = med(np.clip(hs - SC[camp] * hm, 0, None), EDG5)
            bs = [med(np.clip(rng.poisson(hs) - SC[camp] * rng.poisson(hm), 0, None), EDG5) for _ in range(200)]
            pts.append((float(np.mean(js)), m, float(np.nanstd(bs))))
        x = np.array([p[0] for p in pts])
        y = np.log([p[1] for p in pts])
        w = (np.array([p[1] for p in pts]) / np.array([p[2] for p in pts])) ** 2
        X = np.vstack([np.ones_like(x), x - 7.5]).T
        C = np.linalg.inv((X * w[:, None]).T @ X)
        b = C @ (X * w[:, None]).T @ y
        sgn = SIGN[f"{SKEY[camp]}/Y{L}"]["sign"]
        z = b[1] / math.sqrt(C[1, 1])
        prow.append([camp, f"Y{j:02d}", L, sgn, ch, f"{b[1]:.4f}", f"{math.sqrt(C[1, 1]):.4f}", f"{z:.2f}", "yes" if np.sign(b[1]) == sgn else "no"])
with open(OUT + "/profiles.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "bar", "layer", "layer_sign", "channel", "slope_ln_per_X_index", "err", "z", "sign_matches_layer"])
    w.writerows(prow)
for r in prow:
    print("perfil", " | ".join(map(str, r)), flush=True)


# b) coincidencias Y–Y
def yy(camp):
    yc = np.asarray(ORDER[camp]["h18_Y"])
    H = PE[camp][:, yc] > 10.0
    same = H.T.astype(float) @ H.astype(float)
    mix = H[:-1].T.astype(float) @ H[1:].astype(float)
    return same - SC[camp] * mix, same, mix, yc


crow = []
for camp in ("may", "oct"):
    net, same, mix, yc = yy(camp)
    for j in range(16):
        tot = sum(net[j, k] for k in range(16) if k != j)
        top = sorted([(net[j, k], k) for k in range(16) if k != j], reverse=True)[:3]
        crow.append([camp, f"Y{j:02d}", int(yc[j]), f"{tot:.1f}"] + [f"Y{k:02d}(ch{yc[k]}): {v:.1f} ({v / tot:.3f})" for v, k in top])
    res[f"yy_{camp}"] = net.tolist()
with open(OUT + "/yy_coinc.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "bar", "channel", "net_coinc_with_other_Y", "first_partner", "second_partner", "third_partner"])
    w.writerows(crow)
for r in crow:
    if r[1] in ("Y07", "Y08", "Y09", "Y10"):
        print("YY", " | ".join(map(str, r)), flush=True)


# c) filas X × Y
def xyrows(camp):
    o = ORDER[camp]
    xc, yc = np.asarray(o["h18_X"]), np.asarray(o["h18_Y"])
    HX = PE[camp][:, xc] > 10.0
    HY = PE[camp][:, yc] > 10.0
    HX[:, xc == DEAD] = False
    same = HY.T.astype(float) @ HX.astype(float)
    mix = HY[:-1].T.astype(float) @ HX[1:].astype(float)
    return same - SC[camp] * mix, same + SC[camp] ** 2 * mix


XR = {c: xyrows(c) for c in ("may", "oct")}
nm, vm = XR["may"]
xrow = []
j28 = ORDER["may"]["h18_Y"].index(28)
for camp in ("may", "oct"):
    no, vo = XR[camp]
    for k in range(16):
        if camp == "may" and k == j28:
            continue
        a, va, b, vb = nm[j28], vm[j28], no[k], vo[k]
        kk = b.sum() / a.sum()
        m = ((va + vb) > 0)
        m[DEAD_X] = False
        chi = float(((b - kk * a) ** 2 / np.where(m, vb + kk * kk * va, 1))[m].sum())
        xrow.append([camp, f"Y{k:02d}", int(ORDER[camp]["h18_Y"][k]), f"{chi:.2f}", int(m.sum() - 1)])
with open(OUT + "/xy_rows.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["compared_campaign", "bar", "channel", "chi2_shape_vs_may_ch28_row", "ndf"])
    w.writerows(xrow)
best = sorted(xrow, key=lambda r: float(r[3]))[:5]
print("XY fila de ch28 (mayo), las 5 más parecidas:", best, flush=True)

# objetivos por capa con variantes del orden Y de mayo
def dist(L, j):
    return D[j] if SIGN[f"{SKEY['may']}/Y{L}"]["sign"] < 0 else 75.0 - D[j]


def layer_spectra(yorder, name):
    e, g = DATA["may"], G["may"]
    N = len(e)
    s = N / (N - 1)
    xo = np.asarray(ORDER["may"]["h18_X"])
    rows, meds = [], {}
    rngm = np.random.default_rng(1102)
    tagS_X = e["aS"][:, xo] > 10.0 * g[xo]
    tagS_X[:, xo == DEAD] = False
    for L in ("A", "B"):
        chans = [c for i, c in enumerate(yorder) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
        for lo, hi in BINS + [CHECK]:
            js = [j for j in range(16) if lo <= dist(L, j) < hi + (1e-9 if hi == 54 else 0)]
            js = [j for j in js if j != DEAD_X]
            T = tagS_X[:, js].any(1)
            Tn = np.concatenate([T[1:], [False]])
            hf_s, hf_m, hd_s, hd_m = (np.zeros(len(EDG1) - 1) for _ in range(4))
            for c in chans:
                meas = e["aS"][:, c] / g[c]
                for thr, hs_, hm_ in ((3.0, hf_s, hf_m), (5.0, hd_s, hd_m)):
                    sel = meas > thr
                    hs_ += np.histogram(meas[sel & T], bins=EDG1)[0]
                    hm_ += np.histogram(meas[sel & Tn], bins=EDG1)[0]
            net, err = hf_s - s * hf_m, np.sqrt(hf_s + s * s * hf_m)
            for k in range(len(net)):
                rows.append([f"Y{L}", f"{lo}-{hi}", EDG1[k], EDG1[k + 1], int(hf_s[k]), s * hf_m[k], net[k], err[k]])
            m = med(np.clip(hd_s - s * hd_m, 0, None), EDG1)
            bm = [med(np.clip(rngm.poisson(hd_s) - s * rngm.poisson(hd_m), 0, None), EDG1) for _ in range(200)]
            bm = [b for b in bm if np.isfinite(b)]
            meds[f"Y{L}/{lo}-{hi}"] = {"probe_channels": [int(c) for c in chans], "median_pe": m, "median_pe_err": float(np.std(bm)),
                                       "n_net_5": float(hd_s.sum() - s * hd_m.sum()), "n_net_3": float(net.sum())}
    d = OUT + "/" + name
    os.makedirs(d, exist_ok=True)
    with open(d + "/spectra_layers.csv", "w") as f:
        f.write(f"# ch28_test.py: mayo, orden Y {name}: {list(map(int, yorder))}; amplitud aS/G en p.e. medidos, etiqueta > 10 p.e. en el intervalo, sonda > 3 p.e.\n")
        f.write("layer,pos_bin,pe_lo,pe_hi,n_same,n_mix_scaled,n_net,err\n")
        for r in rows:
            f.write(",".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    json.dump({"y_order": list(map(int, yorder)), "rows": meds}, open(d + "/medians.json", "w"), indent=1)
    return meds


nom = list(ORDER["may"]["h18_Y"])
i20, i28 = nom.index(20), nom.index(28)
to_ya = nom.copy()
to_ya[i28] = -1  # sale de YB
swap = nom.copy()
swap[i20], swap[i28] = 28, 20
V = {}
V["nominal"] = layer_spectra(nom, "nominal")
# ch28 en YA además de ch20 (ch28 fuera de YB): se agrega una sonda a YA
ya_extra = [c for c in nom]
ya_extra[i28] = -1


def layer_spectra_extra():
    # YA con ch28 agregado y YB sin ch28: se construye un orden de 17 entradas en el que ch28 aparece en una posición par
    yo = [c for c in nom]
    yo[i28] = DEAD  # DEAD se salta en las sondas
    yo = yo + [28]  # índice 16: par → capa A
    return yo


V["ch28_to_YA"] = layer_spectra(layer_spectra_extra(), "ch28_to_YA")
V["swap_20_28"] = layer_spectra(swap, "swap_20_28")
ref = json.load(open(REPO + "/outputs/hod2018_may_targets/2026-11-04_bootstrap/medians.json"))["rows"]
chk = max(abs(V["nominal"][k]["median_pe"] - ref[k]["median_pe"]) for k in V["nominal"])
print("comprobación nominal contra may_targets medians.json: máx |Δ mediana| =", chk, flush=True)
res["check_nominal_max_abs_dmedian"] = chk
res["medians"] = V
res["profiles"] = prow
res["xy_rows_best5"] = best
json.dump(res, open(OUT + "/ch28.json", "w"), indent=1, default=float)
for k in V["nominal"]:
    print(k, " | ".join(f"{n} {V[n][k]['median_pe']:.3f} ± {V[n][k]['median_pe_err']:.3f}" for n in V), flush=True)
