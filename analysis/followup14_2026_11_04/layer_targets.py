#!/usr/bin/python3 -I
# layer_targets.py — 2026-11-04: objetivos por capa de hod2018 para Geant4, mayo (v8, tiempo bueno) u octubre (run 8 completo); misma selección, mismo ajuste y
# mismo formato que analysis/followup11_2026_11_02/may_layer_targets.py (con "may" el CSV de covarianza reproduce el de 2026-11-02, comprobado con cmp).
#  - G: mayo = reajuste de PROC; octubre = c1, mediana para los canales sin G (ch22, ch31; ch3 muerto). Orden de tiras de cada campaña (stage3.json).
#  - Modos: spectra <camp> <out>: spectra_layers.csv, medians.json y boot_spectra.bin + boot_rows.csv (bootstrap de eventos);
#           merge <camp> <out>: <camp>_layer_targets.csv (errores de f(>T) por covarianza, layer_fit.C);
#           merge_boot <camp> <out>: <camp>_layer_targets_v2_bootstrap.csv (errores de f(>T) del bootstrap, boot_fit.C).
#  - Bootstrap: para cada capa × intervalo, los eventos con la etiqueta del intervalo en el mismo evento o en el siguiente (los que entran en mismo o en
#    mezclados, incluidos los que no tienen sonda sobre el umbral) se remuestrean con reemplazo, NBOOT veces; cada evento lleva todas sus sondas (tiras de la
#    capa > 3 p.e.). Cada réplica da mismo, mezclados, neto y error por bin con la fórmula del espectro nominal, y boot_fit.C repite el ajuste completo.
#    La réplica 0 de boot_spectra.bin es el nominal. Error de f(>T) = SD sobre las réplicas con ajuste convergido (no se descarta por límites de parámetros).
import json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

MODE, CAMP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
EFF = os.path.join(REPO, "outputs/hod2018_eff/2026-10-31_by_position/eff_by_position.json")
ORDER = json.load(open(PK + "/stage3.json"))["strip_orders"][CAMP]["order"]
SIGN = json.load(open(os.path.join(REPO, "outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json")))
SKEY = {"may": "may2026_v8_goodtime", "oct": "oct2026_run8"}[CAMP]
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
CHECK = (29, 54)
TAG, PROBE_FIT, PROBE_MED, NBOOT = 10.0, 3.0, 5.0, 200
F, P_CT = {"may": (1.75, 0.429), "oct": (1.46, 0.315)}[CAMP]
EDG = np.arange(0, 301, 1.0)
RNG = np.random.default_rng({"may": 1102, "oct": 1104}[CAMP])
RNG_BOOT = np.random.default_rng({"may": 2102, "oct": 2104}[CAMP])
LAYERS = ("XA", "XB", "YA", "YB")


def dist(probe, L, j):
    return D[j] if SIGN[f"{SKEY}/{probe}{L}"]["sign"] < 0 else 75.0 - D[j]


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


def load_events():
    if CAMP == "may":
        PKG_G = np.full(32, np.nan)
        for l in open(PK + "/hd64_may_v8_peaks.csv"):
            f = l.strip().split(",")
            if f[0] != "ch" and f[3] in ("ok", "changed"):
                PKG_G[int(f[0])] = float(f[2])
        g = np.where(np.isfinite(PKG_G), PKG_G, np.nanmedian(PKG_G))
        e = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
        e = e[np.argsort(e["t"], kind="stable")]
        return e[mc.good_time(e["t"], "v8")], g
    g = mc.g_c1()
    g = np.where(np.isfinite(g), g, np.nanmedian(g))
    e = xy.load_hits(PK + "/hd64_oct8_hits.bin", 64)
    return e[np.argsort(e["t"], kind="stable")], g


if MODE == "spectra":
    e, g = load_events()
    N = len(e)
    s = N / (N - 1)
    print(f"eventos {CAMP}:", N)
    rows, meds, keys, boot = [], {}, [], []
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
                bm = [median_interp(np.clip(RNG.poisson(hd_s) - s * RNG.poisson(hd_m), 0, None)) for _ in range(NBOOT)]
                bm = [b for b in bm if np.isfinite(b)]
                meds["/".join(key)] = {"probe_channels": [int(c) for c in chans], "tag_strips": js, "tag_channels": [int(refch[j]) for j in js],
                                       "dist_mean_mm": float(np.mean([dist(probe, L, j) for j in js])), "n_same_5": int(hd_s.sum()),
                                       "n_mix_5": float(s * hd_m.sum()), "n_net_5": float(hd_s.sum() - s * hd_m.sum()), "median_pe": med,
                                       "median_pe_err": float(np.std(bm)), "n_net_3": float(net.sum())}
                m = meds["/".join(key)]
                print(f"{key[0]} {key[1]}: etiqueta {m['tag_channels']} n_net(>5) {m['n_net_5']:.1f} mediana {med:.3f} ± {m['median_pe_err']:.3f}", flush=True)
                # bootstrap de eventos: entradas (evento, bin) de todas las sondas > 3 p.e. de los eventos de mismo o mezclados
                U = np.nonzero(T | Tn)[0]
                pos = np.full(N, -1)
                pos[U] = np.arange(len(U))
                evs, bns = [], []
                for c in chans:
                    meas = e["aS"][U, c] / g[c]
                    m3 = (meas > PROBE_FIT) & (meas < EDG[-1])
                    evs.append(np.nonzero(m3)[0])
                    bns.append(meas[m3].astype(int))
                ev, bn = np.concatenate(evs), np.concatenate(bns)
                isS, isM = T[U][ev], Tn[U][ev]
                A = np.zeros((NBOOT + 1, 2, len(EDG) - 1))
                A[0, 0], A[0, 1] = net, err
                for b in range(1, NBOOT + 1):
                    w = RNG_BOOT.multinomial(len(U), np.full(len(U), 1.0 / len(U))).astype(float)
                    hs = np.bincount(bn[isS], w[ev[isS]], len(EDG) - 1)
                    hm = np.bincount(bn[isM], w[ev[isM]], len(EDG) - 1)
                    A[b, 0], A[b, 1] = hs - s * hm, np.sqrt(hs + s * s * hm)
                assert np.allclose(np.bincount(bn[isS], minlength=300), hf_s) and np.allclose(np.bincount(bn[isM], minlength=300), hf_m)
                keys.append(f"{key[0]},{key[1]},{len(U)}")
                boot.append(A)
    os.makedirs(OUT, exist_ok=True)
    with open(OUT + "/spectra_layers.csv", "w") as f:
        f.write(f"# layer_targets.py: {CAMP}, amplitud aS/G en p.e. medidos, etiqueta > {TAG:g} p.e. en el intervalo, "
                f"sonda > {PROBE_FIT:g} p.e.; neto = mismo − mezclados\n")
        f.write("layer,pos_bin,pe_lo,pe_hi,n_same,n_mix_scaled,n_net,err\n")
        for r in rows:
            f.write(",".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
    json.dump({"events": N, "bootstrap": NBOOT, "seed": {"may": 1102, "oct": 1104}[CAMP], "rows": meds}, open(OUT + "/medians.json", "w"), indent=1)
    np.stack(boot).astype("<f8").tofile(OUT + "/boot_spectra.bin")  # (fila, réplica, {neto, error}, bin)
    with open(OUT + "/boot_rows.csv", "w") as f:
        f.write("layer,pos_bin,n_tagged_events_resampled\n" + "\n".join(keys) + "\n")

elif MODE in ("merge", "merge_boot"):
    meds = json.load(open(OUT + "/medians.json"))["rows"]
    fit = {}
    lines = [l for l in open(OUT + "/layer_fit.csv") if not l.startswith("#")]
    nsamp = int(open(OUT + "/layer_fit.csv").readline().split(" con ")[1].split()[0])
    hdr = lines[0].strip().split(",")
    for l in lines[1:]:
        v = dict(zip(hdr, l.strip().split(",")))
        fit[v["layer"] + "/" + v["pos_bin"]] = v
    effall = json.load(open(EFF))["data"][CAMP]
    assert effall["events"] == json.load(open(OUT + "/medians.json"))["events"]
    NEV = effall["events"]
    if CAMP == "may":
        pt = {}
        for l in open(os.path.join(REPO, "outputs/proceeding_revision/2026-11-02_revision/paper_table.csv")):
            f = l.strip().split(",")
            if f[0] != "layer" and f[1] == "far":
                pt[f[0]] = (float(f[3]), float(f[2]))
        for L in LAYERS:
            x = fit[f"{L}/29-54"]
            print(f"comprobación {L} 29-54: moda {float(x['mode_pe']):.3f} (paper_table {pt[L][0]:.3f}), n_fit {float(x['n_net_fit']):.1f} ({pt[L][1]:.1f})")
    boot = {}
    if MODE == "merge_boot":
        import csv
        recs = {}
        for r in csv.DictReader(open(OUT + "/boot_fit.csv")):
            recs.setdefault(r["layer"] + "/" + r["pos_bin"], []).append(r)
        for k, rr in recs.items():
            nom = next(r for r in rr if int(r["rep"]) == 0)
            rep = [r for r in rr if int(r["rep"]) > 0]
            ok = [r for r in rep if int(r["status"]) == 0]
            f5, f10 = np.array([float(r["f5"]) for r in ok]), np.array([float(r["f10"]) for r in ok])
            f5a, f10a = np.array([float(r["f5"]) for r in rep]), np.array([float(r["f10"]) for r in rep])
            boot[k] = dict(nom5=float(nom["f5"]), nom10=float(nom["f10"]), n=len(ok), nall=len(rep), e5=float(f5.std(ddof=1)), e10=float(f10.std(ddof=1)),
                           e5a=float(f5a.std(ddof=1)), e10a=float(f10a.std(ddof=1)),
                           l5=float(np.percentile(f5, 16)), h5=float(np.percentile(f5, 84)), l10=float(np.percentile(f10, 16)), h10=float(np.percentile(f10, 84)))
    out = []
    for L in LAYERS:
        ks = [f"{L}/{lo}-{hi}" for lo, hi in BINS]
        lf = lfit([meds[k]["dist_mean_mm"] for k in ks], [meds[k]["median_pe"] for k in ks], [meds[k]["median_pe_err"] for k in ks])
        print(f"{L}: L_eff {lf.get('L_eff_mm', float('nan')):.1f} ± {lf.get('L_eff_err_mm', float('nan')):.1f} mm, χ² {lf['chi2']:.2f}/{lf['ndf']}, pendiente {lf['slope_per_mm']:.5f} ± {lf['slope_err']:.5f}")
        for b, k in enumerate(ks):
            m, x = meds[k], fit[k]
            r = {"layer": L, "dist_lo_mm": BINS[b][0], "dist_hi_mm": BINS[b][1], "dist_mean_mm": m["dist_mean_mm"],
                 "tag_channels": " ".join(map(str, m["tag_channels"])), "probe_channels": " ".join(map(str, m["probe_channels"])),
                 "n_net": m["n_net_5"], "n_same": m["n_same_5"], "n_mix_scaled": m["n_mix_5"], "median_pe": m["median_pe"], "median_pe_err": m["median_pe_err"],
                 "n_net_fit_4_105": float(x["n_net_fit"]), "mode_pe": float(x["mode_pe"]), "mode_pe_err": float(x["mode_pe_err"]),
                 "fit_chi2": float(x["chi2"]), "fit_ndf": int(x["ndf"]), "fit_status": int(x["status"]), "fit_samples_used": int(x["n_samples_used"]),
                 "f_above_5": float(x["f5"]), "f_above_5_err": float(x["f5_err"]), "f_above_10": float(x["f10"]), "f_above_10_err": float(x["f10_err"])}
            if MODE == "merge_boot":
                z = boot[k]
                r.update(f_above_5_err=z["e5"], f_above_10_err=z["e10"], f_above_5_p16=z["l5"], f_above_5_p84=z["h5"], f_above_10_p16=z["l10"], f_above_10_p84=z["h10"],
                         boot_fits_used=z["n"], f_above_5_err_cov_old=float(x["f5_err"]), f_above_10_err_cov_old=float(x["f10_err"]))
            r["mode_note"] = "mode at the lower fit edge (4 p.e.): not determined" if r["mode_pe"] < 4.0 + 1e-3 else ""
            for T in (5, 10):
                ee = next(z for z in effall[f"{L}_T{T}"] if z["variant"] == "all_ref" and z["bin"] == b)
                assert (ee["dist_lo_mm"], ee["dist_hi_mm"]) == BINS[b]
                r[f"eps_over_mc_T{T}"], r[f"eps_over_mc_T{T}_err"] = ee["eps_over_mc"], ee["eps_over_mc_err"]
            r.update(L_eff_mm=lf.get("L_eff_mm", float("nan")), L_eff_err_mm=lf.get("L_eff_err_mm", float("nan")), L_eff_chi2=lf["chi2"], L_eff_ndf=lf["ndf"])
            out.append(r)
            line = (f"  {k}: n_net {r['n_net']:.1f} mediana {r['median_pe']:.2f} ± {r['median_pe_err']:.2f} moda {r['mode_pe']:.2f} ± {r['mode_pe_err']:.2f} "
                    f"f5 {r['f_above_5']:.4f} ± {r['f_above_5_err']:.4f} f10 {r['f_above_10']:.4f} ± {r['f_above_10_err']:.4f} χ² {r['fit_chi2']:.1f}/{r['fit_ndf']} st {r['fit_status']}")
            if MODE == "merge_boot":
                z = boot[k]
                line += (f" | nominal réplica 0: f5 {z['nom5']:.4f} f10 {z['nom10']:.4f}; convergidas {z['n']}/{z['nall']}; SD todas {z['e5a']:.4f} {z['e10a']:.4f}; "
                         f"p16-84 f5 {z['l5']:.4f}-{z['h5']:.4f} f10 {z['l10']:.4f}-{z['h10']:.4f}; covarianza {x['f5_err']} {x['f10_err']}")
            print(line)
    fmt = lambda v: "" if isinstance(v, float) and not math.isfinite(v) else (f"{v:.6g}" if isinstance(v, float) else str(v))
    desc = {"may": f"hod2018, May 2026 run B (v8), good time, {NEV} triggers; G = PROC refit per channel; ch3 (X strip 12) dead.",
            "oct": f"hod2018, October 2026 run 8, whole run, {NEV} triggers; G = calibration c1, station median for ch22 and ch31 (no c1 G); ch3 (X strip 12) dead."}[CAMP]
    ferr = (f"errors = SD over {nsamp} parameter samples of the fit covariance, samples outside the fit limits dropped (fit_samples_used)." if MODE == "merge" else
            f"mode_pe_err = SD over {nsamp} parameter samples of the fit covariance, samples outside the fit limits dropped (fit_samples_used). f_above_T_err = SD of f(>T) over an event "
            f"bootstrap: {NBOOT} resamples with replacement of the tagged events (same or mixed, each with all its probe strips), full refit of each, converged fits only (boot_fits_used; none dropped for "
            "parameter limits); f_above_T_p16/p84 = 16th and 84th percentiles of the same fits (use these near 1, f cannot exceed 1); f_above_T_err_cov_old = previous covariance-sampling error.")
    name = f"{CAMP}_layer_targets.csv" if MODE == "merge" else f"{CAMP}_layer_targets_v2_bootstrap.csv"
    with open(f"{OUT}/{name}", "w") as f:
        f.write(f"# crosstalk: F = 1/(1 - p) = {F:.2f} (p = {P_CT:.3f}); all amplitudes in MEASURED p.e. (fired cells, amplitude / G, not divided by F); "
                f"the simulation applies the same F\n")
        f.write(f"# {desc} Tag = any strip of the "
                "other coordinate in the position bin with amplitude > 10 p.e.; probe = every live strip of the layer (one entry per strip above threshold); "
                "net = same event - mixed (tag of the next event) x N/(N-1)\n")
        f.write("# n_net, median: probe > 5 p.e., median interpolated in 1 p.e. bins, error = SD of 200 Poisson replicas. mode, f_above_5/10: "
                "Landau(x)Gauss fitted from 4 to 105 p.e. on the probe > 3 p.e. spectrum (n_net_fit_4_105 = net entries in the fit range); mode = maximum of "
                f"the fitted function; f_above_T = fraction of the fitted function above T over [0, 262.5]; {ferr} eps_over_mc: followup 6 (eff_by_position.py, mixed-event estimator, probe > T*G, reference = bin strips > T*G, eps_MC from the "
                "bench MC, CAD, cos^3, ch3 off). L_eff: weighted fit of ln(median) = a - d/L_eff vs dist_mean_mm (mean distance of the 4 tag strips to the "
                "readout end), per layer, repeated on the 4 rows\n")
        f.write("# position = distance to the readout end of the layer, taken as its brighter end (followup 3) [assumption]; the absolute sense of the "
                "strip order of hod2018 is not fixed\n")
        keys = list(out[0])
        f.write(",".join(keys) + "\n")
        for r in out:
            f.write(",".join(fmt(r[k]) for k in keys) + "\n")
    print("hecho", name)
