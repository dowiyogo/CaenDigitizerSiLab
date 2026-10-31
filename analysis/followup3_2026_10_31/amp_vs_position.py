#!/usr/bin/python3 -I
# amp_vs_position.py — seguimiento 2026-10-31 (3), punto 2: amplitud de hod2018 contra la posición a lo largo de la barra.
# Entradas: los archivos de impactos de PROC (outputs/proceeding_2026_10/work/hd64_<dataset>_hits.bin): aS = amplitud de pico (mV) del pulso
# más alto con su mínimo en la ventana signal de 144 ns. G: octubre c1 (mediana para los canales sin G), mayo reajuste de PROC.
# Mayo v8 con el tiempo bueno de runs/good_intervals_2026.csv (HV-DIAG-5 excluido, commit e654684); octubre run 8 completo.
# Orden de tiras: el de along_bar.py (stage3.json del paquete nuevo): octubre = r4; mayo = el validado en los datos de mayo (X igual a r4).
# Para cada canal sonda P (X o Y) y cada tira de referencia R de la coordenada ortogonal (16 posiciones, centros −15 … +15 mm del área activa,
# es decir d = 21 … 54 mm desde el extremo del índice 0: d = 37.5 + u):
#   mismo evento: P > T·G_P y alguna tira de la otra coordenada > T·G; la referencia es la tira de mayor amplitud (en p.e.) de esa coordenada;
#   mezclados: P del evento k con la otra coordenada del evento k + 1 (escala N/(N − 1)), la misma regla;
#   se histograma la amplitud de P en p.e. (bins de 0.25 p.e.); resta = mismo − mezclados; media y mediana de la resta con errores por
#   remuestreo de Poisson (200 réplicas) de los dos histogramas. T = 5 p.e. (la condición P > T trunca el espectro por abajo: es la misma en
#   todas las posiciones).
# Por capa (A = índices pares, B = impares de la coordenada de la sonda): suma de los histogramas de los canales de la capa.
# Ajustes (mínimos cuadrados ponderados de la media de la resta contra d): exponencial m(d) = A exp(d/λ) (sobre ln m) → L_eff = |λ| y el signo
# de λ; lineal m(d) = a + b (d − 37.5) → L_lin = a/|b|; razón de los extremos m(54)/m(21) de cada ajuste.
# Espectro de coincidencias X–Y verdaderas: eventos con algún impacto X y alguno Y (> T); amplitud de cada canal de la capa con impacto;
# mezclados con la coordenada ortogonal del evento siguiente; resta en bins de 1 p.e. hasta el recorte del digitalizador de cada canal
# (máximo de aS observado / G).
# Excluidos como sonda y como referencia: ch3 (muerto), ch29 en octubre, ch30 en mayo.
# Salida: <out>/<run>/profile_<layer o chNN>.csv, fits.json, spectrum_<layer>.csv, clip.csv. Uso: amp_vs_position.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
ST3 = json.load(open(os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work/stage3.json")))["strip_orders"]
T = 5.0
U = -15.0 + 2.0 * np.arange(16)
D = 37.5 + U  # mm desde el extremo del índice 0
BW, PMAX = 0.25, 300.0
EDGES = np.arange(0, PMAX + BW, BW)
CEN = 0.5 * (EDGES[1:] + EDGES[:-1])
RNG = np.random.default_rng(31)
NBOOT = 200


def g_may():
    g = np.full(32, np.nan)
    for r in csv.DictReader(open(W + "/hd64_may_v8_peaks.csv")):
        if r["status"] in ("ok", "changed"):
            g[int(r["ch"])] = float(r["G_fit"])
    return g


def fill(g):
    return np.where(np.isfinite(g), g, np.nanmedian(g[np.isfinite(g)]))


RUNS = {"oct2026_run8": ("hd64_oct8_hits.bin", None, "oct", fill(mc.g_c1()), {3, 29}),
        "may2026_v8_goodtime": ("hd64_may_v8_hits.bin", "v8", "may", fill(g_may()), {3, 30})}


def stats(hs, hm, s):
    """media y mediana de la resta (hs − s·hm), con errores por remuestreo de Poisson."""
    def mm(a, b):
        net = np.clip(a - s * b, 0, None)
        n = net.sum()
        if n <= 0:
            return np.nan, np.nan, 0.0
        mean = float((net * CEN).sum() / n)
        cs = np.cumsum(net)
        med = float(CEN[np.searchsorted(cs, 0.5 * n)])
        return mean, med, n
    mean, med, n = mm(hs, hm)
    if not np.isfinite(mean):
        return {"n_same": int(hs.sum()), "n_mix_scaled": float(s * hm.sum()), "n_net": 0.0}
    bm, bd = [], []
    for _ in range(NBOOT):
        m_, d_, _ = mm(RNG.poisson(hs), RNG.poisson(hm))
        if np.isfinite(m_):
            bm.append(m_); bd.append(d_)
    def m_only(h):
        n_ = h.sum()
        if n_ <= 0:
            return None, None
        cs = np.cumsum(h)
        return float((h * CEN).sum() / n_), float(CEN[np.searchsorted(cs, 0.5 * n_)])
    ms, ds = m_only(hs); mx, dx = m_only(hm)
    return {"n_same": int(hs.sum()), "n_mix_scaled": float(s * hm.sum()), "n_net": float(n), "mean_same": ms, "median_same": ds,
            "mean_mix": mx, "median_mix": dx, "mean_net": mean, "mean_net_err": float(np.std(bm)), "median_net": med, "median_net_err": float(np.std(bd))}


def fit(d, m, e):
    ok = np.isfinite(m) & np.isfinite(e) & (e > 0) & (m > 0)
    d, m, e = d[ok], m[ok], e[ok]
    if len(d) < 4:
        return None
    out = {"n_points": int(len(d))}
    # exponencial sobre ln m
    y, sy = np.log(m), e / m
    X = np.vstack([np.ones_like(d), d - 37.5]).T
    w = 1 / sy ** 2
    C = np.linalg.inv((X * w[:, None]).T @ X); b = C @ (X * w[:, None]).T @ y
    chi = float((((y - X @ b) ** 2) * w).sum())
    k, ek = b[1], math.sqrt(C[1, 1])
    out["exp"] = {"slope_per_mm": float(k), "slope_err": ek, "lambda_mm": float(1 / k) if k != 0 else None, "L_eff_mm": float(abs(1 / k)) if k != 0 else None,
                  "L_eff_err_mm": float(ek / k ** 2) if k != 0 else None, "sign": int(np.sign(k)), "ratio_54_over_21": float(math.exp(33 * k)),
                  "ratio_err": float(math.exp(33 * k) * 33 * ek), "chi2": chi, "ndf": int(len(d) - 2), "amp_at_37p5_pe": float(math.exp(b[0]))}
    # lineal
    w = 1 / e ** 2
    C = np.linalg.inv((X * w[:, None]).T @ X); b = C @ (X * w[:, None]).T @ m
    chi = float((((m - X @ b) ** 2) * w).sum())
    a, sl = b
    ea, es = math.sqrt(C[0, 0]), math.sqrt(C[1, 1])
    L = abs(a / sl) if sl != 0 else None
    eL = L * math.sqrt((ea / a) ** 2 + (es / sl) ** 2 - 2 * C[0, 1] / (a * sl)) if L else None
    r = (a + 16.5 * sl) / (a - 16.5 * sl)
    dr_da = ((a - 16.5 * sl) - (a + 16.5 * sl)) / (a - 16.5 * sl) ** 2
    dr_ds = (16.5 * (a - 16.5 * sl) + 16.5 * (a + 16.5 * sl)) / (a - 16.5 * sl) ** 2
    er = math.sqrt(dr_da ** 2 * C[0, 0] + dr_ds ** 2 * C[1, 1] + 2 * dr_da * dr_ds * C[0, 1])
    out["lin"] = {"a_pe": float(a), "a_err": ea, "slope_pe_per_mm": float(sl), "slope_err": es, "L_lin_mm": L, "L_lin_err_mm": eL, "sign": int(np.sign(sl)),
                  "ratio_54_over_21": float(r), "ratio_err": float(er), "chi2": chi, "ndf": int(len(d) - 2)}
    return out


summary = {"T_pe": T, "bins_pe": BW, "positions_mm_from_index0_end": D.tolist(), "runs": {}}
for run, (fn, good, ok, G, excl) in RUNS.items():
    e = xy.load_hits(os.path.join(W, fn), 64)
    if good:
        e = e[mc.good_time(e["t"], good)]
    N = len(e)
    s = N / (N - 1)
    order = {"X": [c for c in ST3[ok]["order"]["h18_X"]], "Y": [c for c in ST3[ok]["order"]["h18_Y"]]}
    pe = e["aS"][:, :32] / G[None, :32]
    hit = pe > T
    for c in excl:
        hit[:, c] = False
    od = os.path.join(OUT, run)
    os.makedirs(od, exist_ok=True)
    clip = {c: float(e["aS"][:, c].max() / G[c]) for c in range(32) if c not in excl}
    with open(od + "/clip.csv", "w") as f:
        f.write("ch,G_mV,max_aS_mV,clip_pe\n")
        for c in sorted(clip):
            f.write(f"{c},{G[c]:.4f},{e['aS'][:, c].max():.1f},{clip[c]:.1f}\n")
    R = {"N_triggers": int(N), "excluded": sorted(excl), "order": order, "channels": {}, "layers": {}}
    for coord, other in (("X", "Y"), ("Y", "X")):
        ref = order[other]  # índice j -> canal
        refpe = np.where(hit[:, ref], pe[:, ref], -1.0)  # (N,16)
        anyref = (refpe > 0).any(1)
        jbest = np.argmax(refpe, axis=1)
        # mezclados: referencia del evento siguiente
        anyref_n = np.concatenate([anyref[1:], [False]])
        jbest_n = np.concatenate([jbest[1:], [0]])
        layer_h = {L: (np.zeros((16, len(CEN))), np.zeros((16, len(CEN)))) for L in ("A", "B")}
        for i, c in enumerate(order[coord]):
            if c in excl:
                continue
            L = "A" if i % 2 == 0 else "B"
            hs, hm = np.zeros((16, len(CEN))), np.zeros((16, len(CEN)))
            m_s = hit[:, c] & anyref
            m_m = hit[:, c] & anyref_n
            for j in range(16):
                if ref[j] in excl:
                    continue
                hs[j] = np.histogram(pe[m_s & (jbest == j), c], bins=EDGES)[0]
                hm[j] = np.histogram(pe[m_m & (jbest_n == j), c], bins=EDGES)[0]
            layer_h[L][0][:] += hs; layer_h[L][1][:] += hm
            rows = [dict(j=j, d_mm=float(D[j]), ref_ch=int(ref[j]), **stats(hs[j], hm[j], s)) if ref[j] not in excl else dict(j=j, d_mm=float(D[j]), ref_ch=int(ref[j]), excluded=True)
                    for j in range(16)]
            dd = np.array([r["d_mm"] for r in rows]); mm_ = np.array([r.get("mean_net", np.nan) for r in rows]); ee = np.array([r.get("mean_net_err", np.nan) for r in rows])
            R["channels"][f"{coord}{i:02d}_ch{c:02d}"] = {"coord": coord, "index": i, "ch": c, "layer": L, "rows": rows, "fit_mean": fit(dd, mm_, ee)}
            with open(f"{od}/profile_{coord}_ch{c:02d}.csv", "w") as f:
                f.write("j,d_mm,ref_ch,n_same,n_mix_scaled,n_net,mean_same,median_same,mean_mix,median_mix,mean_net,mean_net_err,median_net,median_net_err\n")
                for r in rows:
                    f.write(",".join(str(r.get(k, "")) for k in ("j", "d_mm", "ref_ch", "n_same", "n_mix_scaled", "n_net", "mean_same", "median_same", "mean_mix", "median_mix",
                                                                   "mean_net", "mean_net_err", "median_net", "median_net_err")) + "\n")
        for L in ("A", "B"):
            hs, hm = layer_h[L]
            rows = [dict(j=j, d_mm=float(D[j]), ref_ch=int(ref[j]), **stats(hs[j], hm[j], s)) if ref[j] not in excl else dict(j=j, d_mm=float(D[j]), ref_ch=int(ref[j]), excluded=True)
                    for j in range(16)]
            dd = np.array([r["d_mm"] for r in rows])
            fm = fit(dd, np.array([r.get("mean_net", np.nan) for r in rows]), np.array([r.get("mean_net_err", np.nan) for r in rows]))
            fd = fit(dd, np.array([r.get("median_net", np.nan) for r in rows]), np.array([r.get("median_net_err", np.nan) for r in rows]))
            R["layers"][f"{coord}{L}"] = {"rows": rows, "fit_mean": fm, "fit_median": fd,
                                          "channels": [c for i, c in enumerate(order[coord]) if (i % 2 == (0 if L == "A" else 1)) and c not in excl]}
            with open(f"{od}/profile_{coord}{L}.csv", "w") as f:
                f.write("j,d_mm,ref_ch,n_same,n_mix_scaled,n_net,mean_same,median_same,mean_mix,median_mix,mean_net,mean_net_err,median_net,median_net_err\n")
                for r in rows:
                    f.write(",".join(str(r.get(k, "")) for k in ("j", "d_mm", "ref_ch", "n_same", "n_mix_scaled", "n_net", "mean_same", "median_same", "mean_mix", "median_mix",
                                                                   "mean_net", "mean_net_err", "median_net", "median_net_err")) + "\n")
    # espectro de coincidencias X–Y verdaderas por capa (1 p.e.)
    anyX = hit[:, order["X"]].any(1); anyY = hit[:, order["Y"]].any(1)
    anyX_n = np.concatenate([anyX[1:], [False]]); anyY_n = np.concatenate([anyY[1:], [False]])
    e1 = np.arange(0, 301, 1.0)
    R["spectra"] = {}
    for coord, oth_s, oth_m in (("X", anyY, anyY_n), ("Y", anyX, anyX_n)):
        for L in ("A", "B"):
            chans = [c for i, c in enumerate(order[coord]) if (i % 2 == (0 if L == "A" else 1)) and c not in excl]
            hs = sum(np.histogram(pe[hit[:, c] & oth_s, c], bins=e1)[0] for c in chans)
            hm = sum(np.histogram(pe[hit[:, c] & oth_m, c], bins=e1)[0] for c in chans)
            net = hs - s * hm
            clip_min = min(clip[c] for c in chans)
            with open(f"{od}/spectrum_{coord}{L}.csv", "w") as f:
                f.write(f"# {run} capa {coord}{L} (canales {chans}); coincidencia X–Y, mismo evento − mezclados × {s:.6f}; recorte mínimo de la capa {clip_min:.1f} p.e.\n")
                f.write("pe_lo,pe_hi,n_same,n_mix_scaled,n_net\n")
                for k in range(len(net)):
                    f.write(f"{k},{k + 1},{hs[k]},{s * hm[k]:.3f},{net[k]:.3f}\n")
            n = net[5:].clip(0)
            cs = np.cumsum(n)
            R["spectra"][f"{coord}{L}"] = {"channels": chans, "n_net_total": float(net.sum()), "median_net_pe": float(5 + np.searchsorted(cs, 0.5 * cs[-1]) + 0.5),
                                           "clip_pe_min_of_layer": clip_min, "clip_pe_by_channel": {c: clip[c] for c in chans}}
    summary["runs"][run] = R
    print(run, "N", N)
    for k, v in R["layers"].items():
        fm = v["fit_mean"]
        if fm:
            ex, li = fm["exp"], fm["lin"]
            print(f"  capa {k}: exp L_eff {ex['L_eff_mm']:.1f} ± {ex['L_eff_err_mm']:.1f} mm, signo {ex['sign']:+d}, razón 54/21 {ex['ratio_54_over_21']:.3f} ± {ex['ratio_err']:.3f}, "
                  f"χ² {ex['chi2']:.1f}/{ex['ndf']} | lin L {li['L_lin_mm']:.1f} ± {li['L_lin_err_mm']:.1f} mm, razón {li['ratio_54_over_21']:.3f} ± {li['ratio_err']:.3f} | "
                  f"amp(37.5) {ex['amp_at_37p5_pe']:.1f} p.e.")
    for k, v in R["spectra"].items():
        print(f"  espectro {k}: netos {v['n_net_total']:.0f}, mediana {v['median_net_pe']}, recorte mínimo {v['clip_pe_min_of_layer']:.0f} p.e.")
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=lambda o: None if (isinstance(o, float) and not math.isfinite(o)) else float(o))
