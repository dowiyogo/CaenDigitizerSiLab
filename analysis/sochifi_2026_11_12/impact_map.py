#!/usr/bin/python3 -I
# impact_map.py — SOCHIFI 2026 (póster), 2026-11-12: proyecciones medidas de la posición de impacto de los muones en hod2018, run B (mayo v8).
# Reutiliza, sin reimplementarlas, las definiciones del análisis del proceeding (revisión be9e8de):
#  - impactos: work/hd64_may_v8_hits.bin de hitdef_extract.C (aS = pulso más alto de la ventana signal, aC = de la ventana fuera de tiempo);
#  - G por canal: reajuste de PROC (work/hd64_may_v8_peaks.csv, estado ok/changed; mediana para el canal sin G);
#  - impacto de canal: aS > T·G (xy_estimator.plane_hits, comprobado con assert); ch3 (X, posición 12) muerto y fuera;
#  - muón etiquetado: X·Y en tiempo (algún canal X y algún canal Y sobre el umbral); q, rX, rY, aX, aY de xy_estimator.xy_fraction (mezclados);
#  - orden de tiras de mayo: stage3.json strip_orders["may"] (cadenas medidas en v8 + kCfg_2026_05); centro de la tira i = −15 + 2i mm
#    (convención de bench_mc, sentido absoluto sin fijar).
# Celdas (31 por eje, centradas): celda 2i = tira i sola, celda 2i+1 = solapamiento de las tiras i e i+1; bordes −16.5, −14.5, −13.5, …,
# 14.5, 16.5. Topología por eje: una tira sobre el umbral -> celda 2i; dos tiras vecinas (i, i+1) -> celda 2i+1; el resto -> "otra".
# Accidentales por celda (no existe un estimador por celda en el código del proceeding; este es nuevo y usa solo sus ingredientes): con ruido
# independiente del trigger, entre los X·Y en tiempo el patrón de X de un evento sin X verdadero tiene la distribución n_c del patrón en la
# ventana fuera de tiempo, y esos eventos tienen Y con probabilidad P(Ys ∧ sin X verdadero) = P(Ys) − P_XYs, P_XYs = q + (rX − q)·aY. Los
# eventos con X verdadero e Y de ruido (X sin Y verdadero, (rX − q)·aY) también son accidentales y tienen el patrón de los verdaderos:
#   raw_c = (q / P_XYs) · (on_c − N·n_c·(P(Ys) − P_XYs)),   accidental_c = on_c − raw_c;
# la suma sobre las celdas y "otra" es q·N. Con un toy (ruido independiente por tira, dos tiras ruidosas, aX = 0.04-0.35) reproduce el número
# de X·Y verdaderos por patrón observado dentro de la estadística; la resta de mezclados celda a celda (columna accidental_mixed, solo
# diagnóstico) no, porque cuenta como accidentales en la celda c los X de ruido que en tiempo migran a "otra" al sumarse a un muón verdadero.
# Error estadístico: on_c (Poisson) y n_c (binomial de la ventana fuera de tiempo); la escala q/P_XYs es común a todas las celdas y su
# variación (q_mix ± error, q_ctl) se da en summary.json. Lo mismo por canal para la ocupación por tira (n_c -> ocupación de ruido del canal).
# Para el eje y se intercambian X e Y. No se corrige por eficiencia ni aceptancia; los muones con un X de ruido extra quedan en "otra".
# Segmentos: "stationary" = t ≥ 5.25 h desde el primer trigger (142 544 triggers, pedido de la tarea); "stationary_goodtime" = además el
# tiempo bueno de may_common.good_time (la muestra de ε_XY = 0.81).
# Salida (<out>/<segmento>/): impact_map_runB_T5.csv, impact_map_runB_T10.csv, per_bar_occupancy.csv, topology.csv, asymmetry.csv, summary.json.
# Uso: impact_map.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ORD = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]
XO, YO = np.asarray(ORD["h18_X"]), np.asarray(ORD["h18_Y"])
DEAD, T_STAT_H = 3, 5.25
EDGES = np.array([-16.5] + [-14.5 + k for k in range(30)] + [16.5])
assert len(EDGES) == 32 and EDGES[-2] == 14.5
NC = 31
# celdas afectadas por ch3 (posición 12 de X): 23 = solap. 11|12 (no se puede formar), 24 = tira 12 sola (sin cuentas), 25 = solap. 12|13;
# 22 y 26 = tiras 11 y 13 solas, reciben los cruces de los solapamientos con la 12
DEAD_FLAG = {"x": {22: "neighbour_of_dead_ch3", 23: "dead_ch3", 24: "dead_ch3", 25: "dead_ch3", 26: "neighbour_of_dead_ch3"}, "y": {}}
# uniones del orden de mayo que no salen de las cadenas medidas en v8 (stage1): X 3|4 (ch2–ch11), X 5|6 (ch4–ch9), Y 10|11 (ch21–ch29)
JUNCTION_FLAG = {"x": {7: "junction_from_cabling", 11: "junction_from_cabling"}, "y": {21: "junction_from_cabling"}}
# canales con intervalos fuera del tiempo bueno dentro del segmento estacionario (good_intervals_2026.csv): grupo 2 apagado (ch6, ch20) y
# HV-DIAG-5 (ch0, ch6, ch20, ch22, ch31; todos en el arranque inestable)
HV_CH = {0, 6, 20, 22, 31}
NOISY = {2, 9, 15}


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


def cells(H):
    """H: (n, 16) impactos en orden de posición -> celda (0..30), 31 = otra topología, −1 = ninguna tira."""
    n = H.sum(1)
    first = np.argmax(H, 1)
    c = np.where(n > 0, NC, -1)
    c[n == 1] = 2 * first[n == 1]
    adj = (n == 2) & (first < 15) & H[np.arange(len(H)), np.minimum(first + 1, 15)]
    c[adj] = 2 * first[adj] + 1
    return c


def hist(c):
    """31 celdas + "otra" (índice 31); los eventos sin tiras (−1) no cuentan."""
    return np.bincount(c[c >= 0], minlength=NC + 1).astype(float)


def cell_flags(axis, k, seg):
    out = []
    o = XO if axis == "x" else YO
    bars = [k // 2] if k % 2 == 0 else [k // 2, k // 2 + 1]
    if k in DEAD_FLAG[axis]:
        out.append(DEAD_FLAG[axis][k])
    if k in JUNCTION_FLAG[axis]:
        out.append(JUNCTION_FLAG[axis][k])
    if any(int(o[b]) in NOISY for b in bars):
        out.append("noisy_channel")
    if seg == "stationary" and any(int(o[b]) in HV_CH for b in bars):
        out.append("hv_interval_in_segment")
    return ";".join(out)


def cond_estimate(on, nz, N, q, PXYs, PYs, nz_count):
    """raw = (q/P_XYs)(on − N·n·(P(Ys) − P_XYs)); n = nz_count/N. Devuelve raw, error estadístico."""
    k = q / PXYs
    pn = PYs - PXYs
    raw = k * (on - nz_count * pn)
    n = nz_count / N
    err = k * np.sqrt(on + (N * pn) ** 2 * n * (1 - n) / N)
    return raw, err


G = g_refit()
e_all = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
t0 = float(e_all["t"].min())
stat = e_all["t"] >= t0 + T_STAT_H * 3600
SEGS = {"stationary": stat, "stationary_goodtime": stat & mc.good_time(e_all["t"], "v8")}
XCH = [c for c in range(16) if c != DEAD]
YCH = list(range(16, 32))
GG = np.concatenate([G, np.full(32, np.nan)])
summary_all = {"t_first_utc_s": t0, "stationary_cut_h": T_STAT_H, "edges_mm": EDGES.tolist(), "strip_order_may": {"X": XO.tolist(), "Y": YO.tolist()}}
for seg, msk in SEGS.items():
    e = e_all[msk]
    N = len(e)
    s = N / (N - 1)
    od = os.path.join(OUT, seg)
    os.makedirs(od, exist_ok=True)
    summ = {"segment": seg, "triggers": int(N), "duration_h": float((e["t"].max() - e["t"].min()) / 3600)}
    bars_rows, topo_rows, asym_rows = [], [], []
    for T in (5, 10):
        thr = T * G
        HX = e["aS"][:, XO] > thr[XO]
        HY = e["aS"][:, YO] > thr[YO]
        HXc = e["aC"][:, XO] > thr[XO]
        HYc = e["aC"][:, YO] > thr[YO]
        HX[:, XO == DEAD] = False
        HXc[:, XO == DEAD] = False
        Xs, Ys, Xc, Yc = HX.any(1), HY.any(1), HXc.any(1), HYc.any(1)
        Xs2, Xc2 = xy.plane_hits(e, GG, T, XCH)
        Ys2, Yc2 = xy.plane_hits(e, GG, T, YCH)
        assert (Xs == Xs2).all() and (Ys == Ys2).all() and (Xc == Xc2).all() and (Yc == Yc2).all()
        r = xy.xy_fraction(Xs, Ys, Xc, Yc)
        q, rX, rY, aX, aY = r["f_mix"], r["rX"], r["rY"], r["aX"], r["aY"]
        on = Xs & Ys
        m1 = np.zeros(N, bool)
        m1[:-1] = Xs[:-1] & Ys[1:]
        m2 = np.zeros(N, bool)
        m2[:-1] = Ys[:-1] & Xs[1:]
        res_T = {"N_XY": r["N_XY"], "acc_mix": r["acc_mix"], "acc_ctl": r["acc_ctl"], "rX": rX, "rY": rY, "aX": aX, "aY": aY, "q_mix": q,
                 "q_mix_err": r["f_mix_err"], "q_ctl": r["f_ctl"], "q_ctl_err": r["f_ctl_err"], "qN": q * N, "qN_err": r["f_mix_err"] * N}
        rows = []
        for axis, H, Hc, Ps, Po, ro, ao in (("x", HX, HXc, Ys.mean(), Xs.mean(), rX, aY), ("y", HY, HYc, Xs.mean(), Ys.mean(), rY, aX)):
            # P_XYs para el eje: P(este plano verdadero ∧ el otro con impacto) = q + (r_este − q)·a_otro; P(otro con impacto) = Ps
            PXYs = q + (ro - q) * ao
            c, cc = cells(H), cells(Hc)
            h_on = hist(c[on])
            h_nz = hist(cc)  # patrón de ruido solo, todos los triggers (ventana fuera de tiempo)
            raw, err = cond_estimate(h_on, None, N, q, PXYs, Ps, h_nz)
            # diagnóstico: resta de mezclados celda a celda (sesgada, ver cabecera)
            if axis == "x":
                hm = 0.5 * (hist(c[:-1][m1[:-1]]) + hist(c[1:][m2[:-1]])) * s
            else:
                hm = 0.5 * (hist(c[1:][m1[:-1]]) + hist(c[:-1][m2[:-1]])) * s
            for qq, tag in ((q + r["f_mix_err"], "q_mix_plus_err"), (q - r["f_mix_err"], "q_mix_minus_err"), (r["f_ctl"], "q_ctl")):
                rr_, _ = cond_estimate(h_on, None, N, qq, qq + (ro - qq) * ao, Ps, h_nz)
                res_T[f"{axis}_cells_total_with_{tag}"] = float(rr_[:NC].sum())
            res_T[f"{axis}_P_XYs"] = PXYs
            res_T[f"{axis}_scale_q_over_PXYs"] = q / PXYs
            res_T[f"{axis}_other_on_time"] = float(h_on[NC])
            res_T[f"{axis}_other_raw"] = float(raw[NC])
            res_T[f"{axis}_other_err"] = float(err[NC])
            res_T[f"{axis}_cells_raw_total"] = float(raw[:NC].sum())
            res_T[f"{axis}_cells_raw_total_err"] = float(math.sqrt((err[:NC] ** 2).sum()))
            res_T[f"{axis}_all_raw_total"] = float(raw.sum())
            for k in range(NC):
                w = EDGES[k + 1] - EDGES[k]
                o = XO if axis == "x" else YO
                rows.append({"axis": axis, "cell_index": k, "edge_low_mm": EDGES[k], "edge_high_mm": EDGES[k + 1],
                             "edge_low_from_corner_mm": EDGES[k] + 16.5, "edge_high_from_corner_mm": EDGES[k + 1] + 16.5,
                             "bars": f"{k // 2}" if k % 2 == 0 else f"{k // 2}|{k // 2 + 1}",
                             "channels": f"{o[k // 2]}" if k % 2 == 0 else f"{o[k // 2]}|{o[k // 2 + 1]}",
                             "raw_count": raw[k], "triggers_denominator": N, "counts_per_trigger": raw[k] / N, "stat_error_per_trigger": err[k] / N,
                             "on_time_count": int(h_on[k]), "accidental_estimate": h_on[k] - raw[k], "stat_error_count": err[k],
                             "counts_per_mm": raw[k] / w, "counts_per_mm_err": err[k] / w, "noise_pattern_count_offtime": int(h_nz[k]),
                             "accidental_mixed_diagnostic": hm[k], "dead_channel_flag": int(k in DEAD_FLAG[axis]), "flags": cell_flags(axis, k, seg)})
            # asimetría: celdas 0-14 (< −0.5 mm) contra 16-30 (> 0.5 mm); la celda 15 [−0.5, 0.5] aparte
            for variant, excl in (("all_cells", set()), ("without_dead_ch3_cells_and_mirror", {22, 23, 24, 25, 26, 4, 5, 6, 7, 8})):
                neg = [k for k in range(15) if k not in excl]
                pos = [k for k in range(16, 31) if k not in excl]
                Nn, Np = raw[neg].sum(), raw[pos].sum()
                en, ep = math.sqrt((err[neg] ** 2).sum()), math.sqrt((err[pos] ** 2).sum())
                A = (Np - Nn) / (Np + Nn)
                eA = 2 * math.sqrt((Nn * ep) ** 2 + (Np * en) ** 2) / (Np + Nn) ** 2
                asym_rows.append({"segment": seg, "T_pe": T, "axis": axis, "variant": variant, "N_negative": Nn, "N_negative_err": en, "N_positive": Np,
                                  "N_positive_err": ep, "N_centre_cell15": raw[15], "N_centre_cell15_err": err[15], "asymmetry_pos_minus_neg": A,
                                  "asymmetry_err": eA})
            topo_rows.append({"segment": seg, "T_pe": T, "axis": axis, "tagged_on_time": int(on.sum()), "tagged_after_subtraction_qN": q * N,
                              "tagged_after_subtraction_err": r["f_mix_err"] * N, "other_on_time": int(h_on[NC]),
                              "other_fraction_on_time": h_on[NC] / on.sum(), "other_raw": float(raw[NC]), "other_raw_err": float(err[NC]),
                              "other_fraction_raw": float(raw[NC] / raw.sum()), "n_strips_ge3_on_time": int(((H.sum(1) >= 3) & on).sum()),
                              "two_non_adjacent_on_time": int(((H.sum(1) == 2) & on & (c == NC)).sum())})
            # ocupación por tira de los mismos eventos: mismo estimador con la ocupación de ruido del canal (ventana fuera de tiempo)
            occ_on = H[on].sum(0).astype(float)
            occ_nz = Hc.sum(0).astype(float)
            ro_, eo_ = cond_estimate(occ_on, None, N, q, PXYs, Ps, occ_nz)
            o = XO if axis == "x" else YO
            for i in range(16):
                bars_rows.append({"segment": seg, "T_pe": T, "plane": axis.upper(), "position_index": i, "sublayer": "A" if i % 2 == 0 else "B",
                                  "bench_centre_mm": -15.0 + 2 * i, "channel": int(o[i]), "dead": int(o[i] == DEAD), "on_time": int(occ_on[i]),
                                  "accidental_estimate": occ_on[i] - ro_[i], "raw_count": ro_[i], "stat_error": eo_[i], "per_trigger": ro_[i] / N,
                                  "noise_hits_offtime": int(occ_nz[i]), "noisy_channel": int(int(o[i]) in NOISY),
                                  "hv_interval_in_segment": int(seg == "stationary" and int(o[i]) in HV_CH)})
        summ[f"T{T}"] = res_T
        keys = list(rows[0])
        with open(os.path.join(od, f"impact_map_runB_T{T}.csv"), "w", newline="") as f:
            f.write(f"# impact_map.py: run B (mayo v8), segmento {seg}, {N} triggers, umbral {T} p.e. por canal (aS > T·G, G reajuste de PROC); muones "
                    "etiquetados X·Y en tiempo; raw_count = (q/P_XYs)(on − N·n_c·(P(Ys) − P_XYs)), accidental_estimate = on_time_count − raw_count (ver "
                    "impact_map.py); sin corrección de eficiencia ni aceptancia; coordenadas centradas en mm (centro de la tira i = −15 + 2i, sentido "
                    "absoluto del orden sin fijar); edge_*_from_corner_mm: origen en la esquina inferior izquierda\n")
            w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
            w.writeheader()
            for rr in rows:
                w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) else v) for k, v in rr.items()})
        print(f"{seg} T{T}: N {N} N_XY {r['N_XY']} acc_mix {r['acc_mix']:.1f} q·N {q * N:.1f} ± {r['f_mix_err'] * N:.1f}; celdas x {res_T['x_cells_raw_total']:.1f} "
              f"otra x {res_T['x_other_raw']:.1f} (en tiempo {res_T['x_other_on_time']:.0f}); celdas y {res_T['y_cells_raw_total']:.1f} otra y {res_T['y_other_raw']:.1f} "
              f"(en tiempo {res_T['y_other_on_time']:.0f}); escala x {res_T['x_scale_q_over_PXYs']:.3f} y {res_T['y_scale_q_over_PXYs']:.3f}", flush=True)
    for name, rr in (("per_bar_occupancy.csv", bars_rows), ("topology.csv", topo_rows), ("asymmetry.csv", asym_rows)):
        with open(os.path.join(od, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rr[0]), lineterminator="\n")
            w.writeheader()
            for x in rr:
                w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) else v) for k, v in x.items()})
    json.dump(summ, open(os.path.join(od, "summary.json"), "w"), indent=1, default=float)
    summary_all[seg] = summ
json.dump(summary_all, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
