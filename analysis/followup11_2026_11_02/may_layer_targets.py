#!/usr/bin/python3 -I
# may_layer_targets.py — 2026-11-02: objetivos por capa de hod2018 en mayo (v8, tiempo bueno, 77 558 triggers) para Geant4, en p.e. MEDIDOS
# (aS/G, celdas disparadas, sin dividir por F). Filas = capa (XA, XB, YA, YB) × intervalo de distancia al extremo de lectura (21-29, 29-37,
# 37-45, 45-54 mm, los del seguimiento 6).
#  - Selección etiquetada: etiqueta = alguna tira de la coordenada ortogonal del intervalo con aS > 10 G (ch3 fuera; tiras de
#    eff_by_position.py / bar_targets.py); sonda = cada tira viva de la capa (ch3 fuera de XA), una entrada por tira sobre el umbral; mezclados =
#    etiqueta del evento siguiente, × N/(N − 1); neto = mismo − mezclados. Igual que mip_paper.py, con los 4 intervalos.
#  - Mediana: sonda > 5 p.e., bins de 1 p.e. medido, netos negativos a 0, mediana interpolada en el bin; error = desviación estándar de 200
#    réplicas de Poisson de los dos histogramas (bar_targets.py).
#  - Ajuste: espectro con sonda > 3 p.e. (spectra_layers.csv) ajustado por layer_fit.C (Landau ⊗ Gauss de 4 a 105 p.e., como mip_paper.C).
#  - ε/ε_MC a 5 y 10 p.e.: definición del seguimiento 6, tomada de outputs/hod2018_eff/2026-10-31_by_position/eff_by_position.json
#    (mayo, variante all_ref; mismo tiempo bueno, 77 558 eventos).
#  - L_eff por capa: ajuste ponderado de ln(mediana) contra la distancia media de las 4 tiras de etiqueta del intervalo (lfit de bar_targets.py).
# Uso: may_layer_targets.py spectra <out>   (escribe spectra_layers.csv y medians.json)
#      may_layer_targets.py merge <out>     (lee además layer_fit.csv y escribe may_layer_targets.csv)
import json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

MODE, OUT = sys.argv[1], sys.argv[2]
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
EFF = os.path.join(REPO, "outputs/hod2018_eff/2026-10-31_by_position/eff_by_position.json")
ORDER = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
CHECK = (29, 54)  # comprobación contra paper_table.csv de mip_paper.C
TAG, PROBE_FIT, PROBE_MED, NBOOT, F, P_CT = 10.0, 3.0, 5.0, 200, 1.75, 0.429
EDG = np.arange(0, 301, 1.0)
RNG = np.random.default_rng(1102)
LAYERS = ("XA", "XB", "YA", "YB")


def dist(probe, L, j):
    return D[j] if SIGN[f"may2026_v8_goodtime/{probe}{L}"]["sign"] < 0 else 75.0 - D[j]


def bin_js(probe, L, lo, hi):
    js = [j for j in range(16) if lo <= dist(probe, L, j) < hi + (1e-9 if hi == 54 else 0)]
    return [j for j in js if not (probe == "Y" and j == DEAD_X)]


def median_interp(net):
    n = net.sum()
    if n <= 0:
        return np.nan
    cs = np.cumsum(net)
    b = int(np.searchsorted(cs, 0.5 * n))
    before = cs[b - 1] if b > 0 else 0.0
    return float(EDG[b] + (0.5 * n - before) / net[b])


def lfit(d, m, e):
    d, m, e = (np.asarray(v, float) for v in (d, m, e))
    y, w = np.log(m), (m / e) ** 2
    X = np.vstack([np.ones_like(d), d - 37.5]).T
    C = np.linalg.inv((X * w[:, None]).T @ X)
    b = C @ (X * w[:, None]).T @ y
    k, ek = float(b[1]), math.sqrt(C[1, 1])
    chi = float((((y - X @ b) ** 2) * w).sum())
    out = {"slope_per_mm": k, "slope_err": ek, "chi2": chi, "ndf": int(len(d) - 2)}
    if k < 0:
        out.update(L_eff_mm=-1 / k, L_eff_err_mm=ek / k ** 2)
    return out


if MODE == "spectra":
    PKG_G = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            PKG_G[int(f[0])] = float(f[2])
    g = np.where(np.isfinite(PKG_G), PKG_G, np.nanmedian(PKG_G))
    e = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
    e = e[np.argsort(e["t"], kind="stable")]
    e = e[mc.good_time(e["t"], "v8")]
    N = len(e)
    s = N / (N - 1)
    print("eventos de mayo v8 en tiempo bueno:", N)
    rows, meds = [], {}
    for probe, refc in (("X", "Y"), ("Y", "X")):
        refch = np.asarray(ORDER[f"h18_{refc}"])
        tagS = e["aS"][:, refch] > TAG * g[refch]
        tagS[:, refch == DEAD] = False
        for L in ("A", "B"):
            chans = [c for i, c in enumerate(ORDER[f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
            for lo, hi in BINS + [CHECK]:
                js = bin_js(probe, L, lo, hi)
                T = tagS[:, js].any(1)
                Tn = np.concatenate([T[1:], [False]])
                hf_s, hf_m, hd_s, hd_m = (np.zeros(len(EDG) - 1) for _ in range(4))
                for c in chans:
                    meas = e["aS"][:, c] / g[c]
                    for thr, hs_, hm_ in ((PROBE_FIT, hf_s, hf_m), (PROBE_MED, hd_s, hd_m)):
                        sel = meas > thr
                        hs_ += np.histogram(meas[sel & T], bins=EDG)[0]
                        hm_ += np.histogram(meas[sel & Tn], bins=EDG)[0]
                key = f"{probe}{L}", f"{lo}-{hi}"
                net, err = hf_s - s * hf_m, np.sqrt(hf_s + s * s * hf_m)
                for k in range(len(net)):
                    rows.append([key[0], key[1], EDG[k], EDG[k + 1], int(hf_s[k]), s * hf_m[k], net[k], err[k]])
                med = median_interp(np.clip(hd_s - s * hd_m, 0, None))
                boot = [median_interp(np.clip(RNG.poisson(hd_s) - s * RNG.poisson(hd_m), 0, None)) for _ in range(NBOOT)]
                boot = [b for b in boot if np.isfinite(b)]
                meds["/".join(key)] = {"probe_channels": [int(c) for c in chans], "tag_strips": js, "tag_channels": [int(refch[j]) for j in js],
                                       "dist_mean_mm": float(np.mean([dist(probe, L, j) for j in js])), "n_same_5": int(hd_s.sum()),
                                       "n_mix_5": float(s * hd_m.sum()), "n_net_5": float(hd_s.sum() - s * hd_m.sum()), "median_pe": med,
                                       "median_pe_err": float(np.std(boot)), "n_net_3": float(net.sum())}
                m = meds["/".join(key)]
                print(f"{key[0]} {key[1]}: etiqueta {m['tag_channels']} n_net(>5) {m['n_net_5']:.1f} mediana {med:.3f} ± {m['median_pe_err']:.3f}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(OUT + "/spectra_layers.csv", "w") as f:
        f.write(f"# may_layer_targets.py: mayo v8 tiempo bueno, amplitud aS/G en p.e. medidos, etiqueta > {TAG:g} p.e. en el intervalo, "
                f"sonda > {PROBE_FIT:g} p.e.; neto = mismo − mezclados\n")
        f.write("layer,pos_bin,pe_lo,pe_hi,n_same,n_mix_scaled,n_net,err\n")
        for r in rows:
            f.write(",".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    json.dump({"events": N, "bootstrap": NBOOT, "seed": 1102, "rows": meds}, open(OUT + "/medians.json", "w"), indent=1)

elif MODE == "merge":
    meds = json.load(open(OUT + "/medians.json"))["rows"]
    fit = {}
    lines = [l for l in open(OUT + "/layer_fit.csv") if not l.startswith("#")]
    nsamp = int(open(OUT + "/layer_fit.csv").readline().split(" con ")[1].split()[0])
    hdr = lines[0].strip().split(",")
    for l in lines[1:]:
        v = dict(zip(hdr, l.strip().split(",")))
        fit[v["layer"] + "/" + v["pos_bin"]] = v
    eff = json.load(open(EFF))["data"]["may"]
    assert json.load(open(EFF))["data"]["may"]["events"] == json.load(open(OUT + "/medians.json"))["events"]
    # comprobación: el ajuste de 29-54 contra paper_table.csv (mip_paper.C, misma selección y mismo ajuste)
    pt = {}
    for l in open(os.path.join(REPO, "outputs/proceeding_revision/2026-11-02_revision/paper_table.csv")):
        f = l.strip().split(",")
        if f[0] != "layer" and f[1] == "far":
            pt[f[0]] = (float(f[3]), float(f[2]))
    for L in LAYERS:
        x = fit[f"{L}/29-54"]
        print(f"comprobación {L} 29-54: moda {float(x['mode_pe']):.3f} (paper_table {pt[L][0]:.3f}), n_fit {float(x['n_net_fit']):.1f} ({pt[L][1]:.1f})")
    out = []
    for L in LAYERS:
        ks = [f"{L}/{lo}-{hi}" for lo, hi in BINS]
        lf = lfit([meds[k]["dist_mean_mm"] for k in ks], [meds[k]["median_pe"] for k in ks], [meds[k]["median_pe_err"] for k in ks])
        print(f"{L}: L_eff {lf.get('L_eff_mm', float('nan')):.1f} ± {lf.get('L_eff_err_mm', float('nan')):.1f} mm, χ² {lf['chi2']:.2f}/{lf['ndf']}")
        for b, k in enumerate(ks):
            m, x = meds[k], fit[k]
            r = {"layer": L, "dist_lo_mm": BINS[b][0], "dist_hi_mm": BINS[b][1], "dist_mean_mm": m["dist_mean_mm"],
                 "tag_channels": " ".join(map(str, m["tag_channels"])), "probe_channels": " ".join(map(str, m["probe_channels"])),
                 "n_net": m["n_net_5"], "n_same": m["n_same_5"], "n_mix_scaled": m["n_mix_5"], "median_pe": m["median_pe"], "median_pe_err": m["median_pe_err"],
                 "n_net_fit_4_105": float(x["n_net_fit"]), "mode_pe": float(x["mode_pe"]), "mode_pe_err": float(x["mode_pe_err"]),
                 "fit_chi2": float(x["chi2"]), "fit_ndf": int(x["ndf"]), "fit_status": int(x["status"]), "fit_samples_used": int(x["n_samples_used"]),
                 "f_above_5": float(x["f5"]), "f_above_5_err": float(x["f5_err"]), "f_above_10": float(x["f10"]), "f_above_10_err": float(x["f10_err"])}
            # moda en el borde inferior del ajuste (4 p.e.): la función ajustada decrece en todo el rango, la moda no está determinada
            r["mode_note"] = "mode at the lower fit edge (4 p.e.): not determined" if r["mode_pe"] < 4.0 + 1e-3 else ""
            for T in (5, 10):
                ee = next(z for z in eff[f"{L}_T{T}"] if z["variant"] == "all_ref" and z["bin"] == b)
                assert (ee["dist_lo_mm"], ee["dist_hi_mm"]) == BINS[b]
                r[f"eps_over_mc_T{T}"], r[f"eps_over_mc_T{T}_err"] = ee["eps_over_mc"], ee["eps_over_mc_err"]
            r.update(L_eff_mm=lf.get("L_eff_mm", float("nan")), L_eff_err_mm=lf.get("L_eff_err_mm", float("nan")), L_eff_chi2=lf["chi2"], L_eff_ndf=lf["ndf"])
            out.append(r)
            print(f"  {k}: n_net {r['n_net']:.1f} mediana {r['median_pe']:.2f} ± {r['median_pe_err']:.2f} moda {r['mode_pe']:.2f} ± {r['mode_pe_err']:.2f} "
                  f"f5 {r['f_above_5']:.3f} f10 {r['f_above_10']:.3f} ε/εMC(5) {r['eps_over_mc_T5']:.3f} ± {r['eps_over_mc_T5_err']:.3f} "
                  f"ε/εMC(10) {r['eps_over_mc_T10']:.3f} ± {r['eps_over_mc_T10_err']:.3f}")
    fmt = lambda v: "" if isinstance(v, float) and not math.isfinite(v) else (f"{v:.6g}" if isinstance(v, float) else str(v))
    with open(OUT + "/may_layer_targets.csv", "w") as f:
        f.write(f"# crosstalk: F = 1/(1 - p) = {F:.2f} (p = {P_CT:.3f}); all amplitudes in MEASURED p.e. (fired cells, amplitude / G, not divided by F); "
                f"the simulation applies the same F\n")
        f.write("# hod2018, May 2026 run B (v8), good time, 77558 triggers; G = PROC refit per channel; ch3 (X strip 12) dead. Tag = any strip of the "
                "other coordinate in the position bin with amplitude > 10 p.e.; probe = every live strip of the layer (one entry per strip above threshold); "
                "net = same event - mixed (tag of the next event) x N/(N-1)\n")
        f.write("# n_net, median: probe > 5 p.e., median interpolated in 1 p.e. bins, error = SD of 200 Poisson replicas. mode, f_above_5/10: "
                "Landau(x)Gauss fitted from 4 to 105 p.e. on the probe > 3 p.e. spectrum (n_net_fit_4_105 = net entries in the fit range); mode = maximum of "
                f"the fitted function; f_above_T = fraction of the fitted function above T over [0, 262.5]; errors = SD over {nsamp} parameter samples of the fit "
                f"covariance, samples outside the fit limits dropped (fit_samples_used). eps_over_mc: followup 6 (eff_by_position.py, mixed-event estimator, probe > T*G, reference = bin strips > T*G, eps_MC from the "
                "bench MC, CAD, cos^3, ch3 off). L_eff: weighted fit of ln(median) = a - d/L_eff vs dist_mean_mm (mean distance of the 4 tag strips to the "
                "readout end), per layer, repeated on the 4 rows\n")
        f.write("# position = distance to the readout end of the layer, taken as its brighter end (followup 3) [assumption]; the absolute sense of the "
                "strip order of hod2018 is not fixed\n")
        keys = list(out[0])
        f.write(",".join(keys) + "\n")
        for r in out:
            f.write(",".join(fmt(r[k]) for k in keys) + "\n")
    print("hecho")
