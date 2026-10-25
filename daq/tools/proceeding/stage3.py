#!/usr/bin/python3 -I
# stage3.py — PROC-1 etapa 3 (solo lectura): mapas y ocupación.
#  - Impactos: work/hd64_<dataset>_hits.bin (hitdef_extract.C, nch 64), definición única; hod2018 con T·G (mayo: G reajustada,
#    A3; octubre: c1), hod2019 con umbral en mV (35 y 3 mV). ch3 fuera (D6).
#  - Orden de tiras (D4, A4): mayo = orden medido de la etapa 1 (kCfg_2026_05 en hod2018; hod2019 kCfg_2026_05 sin validar);
#    octubre = r4 (daq/config/mapping/hodo_mapping_2026_10_r4.json). Tiras cuya posición no fija una cadena confirmada: "hatched".
#  - Mapa X × Y de hod2018 (5 p.e., peso 1/(nX·nY), en tiempo − mezclados; xy_estimator.xy_map), impactos por tira en tiempo −
#    fuera de tiempo (5 y 10 p.e.) con la expectativa del MC de barras (work/bar_mc_g.json, strip_fraction_per_trigger),
#    hod2019 por tira (35 y 3 mV, mayo v8 y octubre run 8), multiplicidad por plano contra Poisson.
#    PROC-2 (R1): el MC de barras usa las placas en las z del CAD.
# Salida: work/stage3.json, work/fig_stage3.json (para proc_plot_stage3.C), tables/occupancy_ratio.tex, tables/multiplicity.tex.
# Uso: stage3.py <dir del paquete>
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": 3, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
DEAD = 3

# ---------------------------------------------------------------- órdenes (coordenada 0 … 15 = centro −15 … 15 mm)
MAY = {"h18_X": [15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14],
       "h18_Y": [16, 24, 17, 25, 18, 26, 19, 27, 20, 28, 21, 29, 22, 30, 23, 31],
       "h19_X": [53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56],
       "h19_Y": [37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40]}
# validada = la tira está en la cadena confirmada más larga de la etapa 1 (X: 10-5-8-7-6-9-4-11, coordenadas 4-11;
# Y: 16 … 21, coordenadas 0-10); el resto se coloca con uniones derivadas o sin confirmar; hod2019 sin validar.
MAY_VALID = {"h18_X": [4 <= k <= 11 for k in range(16)], "h18_Y": [k <= 10 for k in range(16)],
             "h19_X": [False] * 16, "h19_Y": [False] * 16}
r4 = json.load(open(REPO + "/daq/config/mapping/hodo_mapping_2026_10_r4.json"))
OCT, OCT_VALID = {}, {}
for st, key in (("hod2018", "h18"), ("hod2019", "h19")):
    for pl in "XY":
        cs = sorted([c for c in r4["channels"] if c["station"] == st and c["plane"] == pl], key=lambda c: c["strip"])
        assert [c["strip"] for c in cs] == list(range(16))
        OCT[f"{key}_{pl}"] = [c["dig_ch"] for c in cs]
        OCT_VALID[f"{key}_{pl}"] = [c["confirmed"]["strip"] == "measurement" for c in cs]
S["strip_orders"] = {"may": {"order": MAY, "validated": MAY_VALID, "source": "etapa 1 (coincidencias entre vecinas, v8) / kCfg_2026_05"},
                     "oct": {"order": OCT, "validated": OCT_VALID, "source": r4["mapping_version"]}}


# ---------------------------------------------------------------- datos
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
G = {"may": np.concatenate([fill(g_refit()), np.full(32, np.nan)]), "oct": np.concatenate([fill(mc.g_c1()), np.full(32, np.nan)])}
e8 = xy.load_hits(W + "/hd64_may_v8_hits.bin", 64)
DATA = {"may8_good": (e8[mc.good_time(e8["t"], "v8")], "may", "mayo v8, tiempo bueno"), "may8_all": (e8, "may", "mayo v8, todo (run B)"),
        "oct8": (xy.load_hits(W + "/hd64_oct8_hits.bin", 64), "oct", "octubre run 8")}
ORD = {"may": (MAY, MAY_VALID), "oct": (OCT, OCT_VALID)}


def thr(ch, T, cal):
    """umbral en mV: T·G en hod2018 (T en p.e.), T en mV en hod2019"""
    ch = np.asarray(ch)
    return np.where(ch < 32, T * G[cal][np.minimum(ch, 63)], T)


def strip_hits(e, chans, T, cal, field="aS"):
    """booleanos (n, 16) por coordenada; ch3 nunca da impacto"""
    h = e[field][:, chans] > thr(chans, T, cal)
    h[:, np.asarray(chans) == DEAD] = False
    return h


# ---------------------------------------------------------------- 3.1 mapa X × Y de hod2018
BM = json.load(open(W + "/bar_mc_g.json"))
MC = next(r for r in BM["runs"] if r["geometry"] == "CAD")["strip_fraction_per_trigger"]
S["maps"] = {}
for key in ("may8_good", "oct8"):
    e, cal, lab = DATA[key]
    o, _ = ORD[cal]
    X, Y = strip_hits(e, o["h18_X"], 5, cal), strip_hits(e, o["h18_Y"], 5, cal)
    m_in, m_mix = xy.xy_map(X, Y)
    sub = m_in - m_mix
    # los mezclados combinan también la X verdadera de un muón con la Y verdadera del siguiente: el estimador lo devuelve como
    # rX·rY (xy_estimator.py). En el mapa se reparte ese término con el producto de los perfiles por tira en tiempo − fuera
    # de tiempo (aproximación: ignora el peso 1/(nX·nY) en esa parte) para que la suma sea q_mix·N.
    pX = np.clip((X.sum(0) - strip_hits(e, o["h18_X"], 5, cal, "aC").sum(0)).astype(float), 0, None)
    pY = np.clip((Y.sum(0) - strip_hits(e, o["h18_Y"], 5, cal, "aC").sum(0)).astype(float), 0, None)
    Xs, Xc = xy.plane_hits(e, G[cal], 5, [c for c in range(16) if c != DEAD])
    Ys, Yc = xy.plane_hits(e, G[cal], 5, range(16, 32))
    r = xy.xy_fraction(Xs, Ys, Xc, Yc)
    live = [k for k in range(16) if o["h18_X"][k] != DEAD]
    corr_total = r["f_mix"] * r["N"] - sub.sum()
    corr = np.outer(pX, pY) / (pX.sum() * pY.sum()) * corr_total
    subc = sub + corr
    rows, cols = subc[live].sum(1), subc[live].sum(0)
    S["maps"][key] = {"label": lab, "events": int(len(e)), "threshold": "5 p.e.", "in_time": m_in.tolist(), "mixed": m_mix.tolist(), "sub": sub.tolist(), "sub_corrected": subc.tolist(), "rXrY_term_total": float(corr_total),
                      "q_mix_times_N": float(r["f_mix"] * r["N"]), "sum_sub_corrected": float(subc.sum()), "sub_corrected_cells_negative": int((subc < 0).sum()),
                      "sum_in": float(m_in.sum()), "sum_mixed": float(m_mix.sum()), "sum_sub": float(sub.sum()), "N_XY": r["N_XY"], "acc_mix": r["acc_mix"],
                      "sub_projection_X": subc.sum(1).tolist(), "sub_projection_Y": subc.sum(0).tolist(),
                      "sub_cells_negative": int((sub < 0).sum()), "sub_min": float(sub.min()), "sub_max": float(sub.max()),
                      "sub_row_max_over_min_live": float(rows.max() / rows.min()), "sub_col_max_over_min": float(cols.max() / cols.min()),
                      "note": "filas = coordenada X (0 … 15), columnas = coordenada Y; la fila de ch3 queda en cero; con iluminación uniforme el mapa no "
                              "puede validar el orden de tiras"}
    print(key, "mapa: suma en tiempo %.1f (N_XY %d), mezclados %.1f (acc_mix %.1f)" % (m_in.sum(), r["N_XY"], m_mix.sum(), r["acc_mix"]))


# ---------------------------------------------------------------- 3.2 impactos por tira de hod2018
def ratio(v, ev, idx):
    v, ev = np.asarray(v)[idx], np.asarray(ev)[idx]
    i, j = int(np.argmax(v)), int(np.argmin(v))
    r = v[i] / v[j]
    return {"max_over_min": float(r), "err": float(r * math.hypot(ev[i] / v[i], ev[j] / v[j])), "max_coord": int(idx[i]), "min_coord": int(idx[j]),
            "rel_sd": float(v.std() / v.mean())}


def per_strip(e, chans, T, cal):
    n_in = strip_hits(e, chans, T, cal, "aS").sum(0).astype(float)
    n_oot = strip_hits(e, chans, T, cal, "aC").sum(0).astype(float)
    net = n_in - n_oot
    return {"channels": list(chans), "in_time": n_in.tolist(), "out_of_time": n_oot.tolist(), "net": net.tolist(), "net_err": np.sqrt(n_in + n_oot).tolist(),
            "net_per_1000": (net / len(e) * 1e3).tolist(), "net_per_1000_err": (np.sqrt(n_in + n_oot) / len(e) * 1e3).tolist(),
            "z": np.where(n_in + n_oot > 0, net / np.sqrt(np.maximum(n_in + n_oot, 1)), 0).tolist()}


S["occupancy_2018"] = {}
for key in ("may8_good", "oct8"):
    e, cal, lab = DATA[key]
    o, val = ORD[cal]
    d = {"label": lab, "events": int(len(e))}
    for pl in "XY":
        k = f"h18_{pl}"
        idx = np.array([i for i in range(16) if o[k][i] != DEAD])
        mcv = np.asarray(MC[k])
        for T in (5, 10):
            p = per_strip(e, o[k], T, cal)
            net, ne = np.asarray(p["net"]), np.asarray(p["net_err"])
            p["ratio"] = ratio(net, ne, idx)
            # expectativa del MC con detección perfecta, normalizada a la suma medida de las tiras vivas (la eficiencia no se conoce)
            scale = net[idx].sum() / mcv[idx].sum()
            p["mc_expected"] = (mcv * scale).tolist()
            p["mc_scale_eff_equiv"] = float(scale / len(e))
            p["chi2_vs_mc"] = float((((net - mcv * scale)[idx] / ne[idx]) ** 2).sum())
            p["ndf_vs_mc"] = int(len(idx) - 1)
            p["strips_above_3sigma"] = int((np.asarray(p["z"])[idx] > 3).sum())
            d[f"{pl}_T{T}"] = p
            print(key, pl, T, "max/min %.2f ± %.2f (coord %d / %d)" % (p["ratio"]["max_over_min"], p["ratio"]["err"], p["ratio"]["max_coord"],
                                                                      p["ratio"]["min_coord"]), "chi2/ndf vs MC %.1f/%d" % (p["chi2_vs_mc"], p["ndf_vs_mc"]))
        d[f"{pl}_mc"] = ratio(mcv, np.zeros(16), idx) | {"fraction_per_trigger": mcv.tolist()}
    S["occupancy_2018"][key] = d
# aproximación a la figura del manuscrito (50 mV): mismo buscador y ventana, umbral fijo de 50 mV, run B completo, sin restar
S["occupancy_2018_50mV_proxy"] = {}
for pl in "XY":
    k = f"h18_{pl}"
    hi = (e8["aS"][:, MAY[k]] > 50).sum(0).astype(float)
    ho = (e8["aC"][:, MAY[k]] > 50).sum(0).astype(float)
    idx = np.array([i for i in range(16) if MAY[k][i] != DEAD])
    S["occupancy_2018_50mV_proxy"][pl] = {"channels": MAY[k], "in_time": hi.tolist(), "out_of_time": ho.tolist(),
                                          "ratio_in_time": ratio(hi, np.sqrt(hi), idx), "ratio_net": ratio(hi - ho, np.sqrt(hi + ho), idx)}
    print("50 mV", pl, "en tiempo", hi.astype(int).tolist(), "máx/mín %.2f" % S["occupancy_2018_50mV_proxy"][pl]["ratio_in_time"]["max_over_min"])
S["mc_strip_spread"] = next(r for r in BM["runs"] if r["geometry"] == "CAD")["strip_spread"]

# ---------------------------------------------------------------- 3.3 hod2019 por tira
S["occupancy_2019"] = {}
for key in ("may8_all", "may8_good", "oct8"):
    e, cal, lab = DATA[key]
    o, _ = ORD[cal]
    d = {"label": lab, "events": int(len(e))}
    for pl in "XY":
        k = f"h19_{pl}"
        for T in (35, 3):
            p = per_strip(e, o[k], T, cal)
            z = np.asarray(p["z"])
            p["strips_above_3sigma"] = [int(i) for i in np.where(z > 3)[0]]
            p["channels_above_3sigma"] = [int(o[k][i]) for i in np.where(z > 3)[0]]
            d[f"{pl}_{T}mV"] = p
            print(key, k, T, "mV: z>3 en coordenadas", p["strips_above_3sigma"], "canales", p["channels_above_3sigma"])
    S["occupancy_2019"][key] = d
# afirmaciones del manuscrito (índice k de las X de hod2019): se prueban dos lecturas de k
mx = S["occupancy_2019"]["may8_all"]["X_35mV"]
n_in = np.asarray(mx["in_time"])
by_ch = {ch: n_in[i] for i, ch in enumerate(MAY["h19_X"])}
k_ch = np.array([by_ch[48 + k] for k in range(16)])          # lectura 1: k = canal − 48
k_pos = n_in                                                  # lectura 2: k = coordenada de kCfg_2026_05
S["hod2019_claims"] = {
    "manuscript": "4 tiras activas con índices 0, 2, 10, 14; las otras 12 con 30-155 cuentas frente a 280-360; pares > impares vecinas",
    "definition_used": "definición única, ventana signal de 144 ns, amplitud > 35 mV, run B completo, cuentas en tiempo sin restar",
    "reading_k_eq_ch_minus_48": {"counts": k_ch.tolist(), "top4": sorted(np.argsort(k_ch)[-4:].tolist()),
                                 "even_minus_odd_mean": float(k_ch[0::2].mean() - k_ch[1::2].mean())},
    "reading_k_eq_kCfg_coordinate": {"counts": k_pos.tolist(), "top4": sorted(np.argsort(k_pos)[-4:].tolist()),
                                     "even_minus_odd_mean": float(k_pos[0::2].mean() - k_pos[1::2].mean())},
    "may_oct_relation": "no relacionado: en r4 la tira de hod2019 viene de la hoja de octubre (cabling_record) y la regla cara/conector de mayo no "
                        "se cumple (r2); no hay registro que una barra de mayo con un canal de octubre"}
print("hod2019 X 35 mV, k = canal − 48:", k_ch.astype(int).tolist())
print("hod2019 X 35 mV, k = coordenada:", k_pos.astype(int).tolist())
oc = S["occupancy_2019"]["oct8"]
S["oct_ch39_ch47"] = {}
for ch in (39, 47):
    pl = "X" if ch in OCT["h19_X"] else "Y"
    i = OCT[f"h19_{pl}"].index(ch)
    S["oct_ch39_ch47"][str(ch)] = {"plane": pl, "coordinate": i, **{f"{T}mV": {f: oc[f"{pl}_{T}mV"][f][i] for f in ("in_time", "out_of_time", "net", "net_err", "z")}
                                                                    for T in (35, 3)}}
print("octubre ch39/ch47:", json.dumps(S["oct_ch39_ch47"]))


# ---------------------------------------------------------------- 3.4 multiplicidad
def mult(h, kmax=5):
    m = h.sum(1)
    return np.array([(m == k).sum() for k in range(kmax)] + [(m >= kmax).sum()], float), float(m.mean())


S["multiplicity"] = {}
for key in ("may8_good", "oct8"):
    e, cal, lab = DATA[key]
    o, _ = ORD[cal]
    d = {"label": lab, "events": int(len(e))}
    for pl in "XY":
        for T in (5, 10):
            hi, ho = strip_hits(e, o[f"h18_{pl}"], T, cal, "aS"), strip_hits(e, o[f"h18_{pl}"], T, cal, "aC")
            ci, mi = mult(hi)
            co, mo = mult(ho)
            pois = np.array([math.exp(-mi) * mi ** k / math.factorial(k) for k in range(5)])
            pois = np.append(pois, 1 - pois.sum()) * len(e)
            # condicionado a ≥ 1 impacto: multiplicidad del muón (más vecinas) frente a Poisson truncado
            hit = ci[1:].sum()
            d[f"{pl}_T{T}"] = {"in_time": ci.tolist(), "out_of_time": co.tolist(), "poisson_same_mean": pois.tolist(), "mean_in": mi, "mean_oot": mo,
                               "frac_ge2_given_ge1_in": float(ci[2:].sum() / hit), "frac_ge2_given_ge1_oot": float(co[2:].sum() / co[1:].sum()),
                               "frac_ge2_given_ge1_poisson": float(pois[2:].sum() / pois[1:].sum())}
    S["multiplicity"][key] = d

# ---------------------------------------------------------------- tablas
L = ["% PROC-1 etapa 3: impactos por tira de hod2018, en tiempo − fuera de tiempo, definición única; razón máx/mín sobre las 15 tiras X vivas",
     "% (ch3 fuera) y las 16 Y. MC: bar_mc.py (1e8 rayos, cos^3, placas en las z del CAD, R1), cruces por tira con detección perfecta.",
     f"% Commit {git}.", "\\begin{tabular}{llrrrr}", "\\hline",
     "Conjunto & Plano & $\\max/\\min$ (5 p.e.) & $\\max/\\min$ (10 p.e.) & $\\chi^2/\\mathrm{ndf}$ vs MC (5 p.e.) & MC $\\max/\\min$ \\\\", "\\hline"]
for key in ("may8_good", "oct8"):
    d = S["occupancy_2018"][key]
    for pl in "XY":
        a, b = d[f"{pl}_T5"], d[f"{pl}_T10"]
        L.append(f"{d['label']} & {pl} & {a['ratio']['max_over_min']:.2f} $\\pm$ {a['ratio']['err']:.2f} & {b['ratio']['max_over_min']:.2f} $\\pm$ "
                 f"{b['ratio']['err']:.2f} & {a['chi2_vs_mc']:.0f}/{a['ndf_vs_mc']} & {d[pl + '_mc']['max_over_min']:.3f} \\\\")
ms = S["mc_strip_spread"]
L += ["\\hline", f"% MC hod2019: máx/mín X {ms['h19_X']['max_over_min']:.3f}, Y {ms['h19_Y']['max_over_min']:.3f}; desviación relativa hod2018 X "
      f"{ms['h18_X']['rel_sd']:.3f}, Y {ms['h18_Y']['rel_sd']:.3f}, hod2019 X {ms['h19_X']['rel_sd']:.3f}, Y {ms['h19_Y']['rel_sd']:.3f}.", "\\end{tabular}"]
open(OUT + "/tables/occupancy_ratio.tex", "w").write("\n".join(L) + "\n")
M = ["% PROC-1 etapa 3: multiplicidad por plano de hod2018 (número de tiras con impacto), en tiempo, fuera de tiempo y Poisson con la media",
     f"% en tiempo; definición única, ch3 fuera. Commit {git}.", "\\begin{tabular}{lllrrrrrr}", "\\hline",
     "Conjunto & Plano & $T$ & Serie & 0 & 1 & 2 & 3 & 4 & $\\geq 5$ \\\\", "\\hline"]
M[2] = "\\begin{tabular}{llllrrrrrr}"
for key in ("may8_good", "oct8"):
    d = S["multiplicity"][key]
    for pl in "XY":
        for T in (5, 10):
            r = d[f"{pl}_T{T}"]
            for nm, v in (("en tiempo", r["in_time"]), ("fuera de tiempo", r["out_of_time"]), ("Poisson", r["poisson_same_mean"])):
                M.append(f"{d['label']} & {pl} & {T} p.e. & {nm} & " + " & ".join(f"{x:.0f}" for x in v) + " \\\\")
M += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/multiplicity.tex", "w").write("\n".join(M) + "\n")

# ---------------------------------------------------------------- figura
fig = {"maps": {k: {"sub": S["maps"][k]["sub_corrected"], "mixed": S["maps"][k]["mixed"], "label": S["maps"][k]["label"], "dead_X_coord": ORD[c][0]["h18_X"].index(DEAD),
                    "valid_X": ORD[c][1]["h18_X"], "valid_Y": ORD[c][1]["h18_Y"]} for k, c in (("may8_good", "may"), ("oct8", "oct"))},
       "occ": []}
d = S["occupancy_2018"]["may8_good"]
for pl in "XY":
    fig["occ"].append({"panel": f"hod2018 {pl}", "ylabel": "Hits per 1000 triggers", "valid": MAY_VALID[f"h18_{pl}"],
                       "dead": [ch == DEAD for ch in MAY[f"h18_{pl}"]],
                       "series": [{"name": f"{T} p.e.", "y": d[f"{pl}_T{T}"]["net_per_1000"], "e": d[f"{pl}_T{T}"]["net_per_1000_err"]} for T in (5, 10)],
                       "mc": (np.asarray(d[f"{pl}_T5"]["mc_expected"]) / d["events"] * 1e3).tolist(),
                       "mc_err": (np.asarray(d[f"{pl}_T5"]["mc_expected"]) / d["events"] * 1e3 * np.sqrt(1 / np.maximum(np.asarray(
                           next(r for r in BM["runs"] if r["geometry"] == "CAD")["strip_crossings"][f"h18_{pl}"]), 1))).tolist()})
d = S["occupancy_2019"]["may8_all"]
for pl in "XY":
    fig["occ"].append({"panel": f"hod2019 {pl}", "ylabel": "Hits per 1000 triggers", "valid": MAY_VALID[f"h19_{pl}"], "dead": [False] * 16,
                       "series": [{"name": f"{T} mV", "y": d[f"{pl}_{T}mV"]["net_per_1000"], "e": d[f"{pl}_{T}mV"]["net_per_1000_err"]} for T in (35, 3)],
                       "mc": None, "mc_err": None})
json.dump(fig, open(W + "/fig_stage3.json", "w"), indent=1)
json.dump(S, open(W + "/stage3.json", "w"), indent=1, default=float)
print("hecho")
