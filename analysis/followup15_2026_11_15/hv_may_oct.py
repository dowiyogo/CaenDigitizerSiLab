#!/usr/bin/python3 -I
# hv_may_oct.py — seguimiento 2026-11-15 (15, punto 4): HV de mayo (capturas del panel GECO, docs/hardware/hv_may2026_geco_screenshots.csv) contra el
# registro de hvmon en octubre run 8, y la ganancia y la luz de mayo/octubre por grupo HV. Solo lectura (no toca el SY5527).
#  1. Referencia de octubre: mediana de VMon, IMon y V0Set por canal (CH00-07, Pw = 1) en /data/daq/slowcontrol/2026/hv_2026102{2,3}.csv.gz durante
#     run 8 (primer a último disparo de hd64_oct8_hits.bin). ΔV = VMon(mayo 2026-05-28 01:10 local) − VMon(octubre run 8) [M], con el supuesto [S] de
#     V0Set de mayo constante durante el tiempo bueno.
#  2. Por grupo HV (grupo del MPPC de seguimiento 7, ratio_by_mppc.csv: mayo y octubre emparejados por MPPC, grupo comprobado 5 de 5; solo hod2018,
#     hod2019 no tiene G): G_mayo/G_oct por MPPC (gain_by_mppc.csv; fuera ch3 y los MPPC sin G de octubre medida) y razón de p.e. primarios mayo/octubre
#     (seguimiento 7; exclusión principal de seguimiento 7: barras Y11, Y13, Y14; variante con todas). Mediana del grupo, dispersión 1.4826·MAD.
#     Regresión de ln(razón) contra ΔV por grupo: error de cada mediana = 1.2533·s/√n con s = dispersión conjunta dentro de los grupos (MAD de los
#     residuos respecto de la mediana de su grupo, todos los MPPC), mínimos cuadrados ponderados; variante sin pesos. Interpretación [I]: con
#     G ∝ (V − V_bd), ln(G_M/G_O) ≈ (ΔV − ΔV_bd)/OV, pendiente = 1/OV, ΔV_bd(mayo − octubre) = −ordenada/pendiente; ΔT = ΔV_bd / 60 mV/°C [S].
#     Lo mismo para la luz [S: luz ∝ OV].
#  3. V_bd de mayo por grupo: los tres puntos I-V de mayo de cada canal con el modelo (b) de seguimiento 8, I = I_0 + a (V − V_bd)^n, I_0 y n fijos al
#     ajuste de octubre del grupo (b_free si no está en el borde de la rejilla; si no, b_fixed; si tampoco, n = mediana de los grupos determinados [S]);
#     a y V_bd libres (1 grado de libertad); σ_I = √(0.002² + (0.05·I)²) µA como seguimiento 8 (sin sd: lectura de pantalla). Rejilla en V_bd
#     (paso 1 mV), a por mínimos cuadrados, error del perfil Δχ² = 1. Variante con n libre (0 grados de libertad, solo si tiene solución).
# Salida: <out>/hv_oct_run8.csv, hv_delta.csv, ratio_vs_dv.csv, vbd_may.csv, hv_may_oct.json. Uso: hv_may_oct.py <out dir>
import csv, gzip, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
MAYCSV = REPO + "/docs/hardware/hv_may2026_geco_screenshots.csv"
R7 = REPO + "/outputs/hod2018_recabling/2026-10-31_groups/"
IV = REPO + "/outputs/hod2018_iv/2026-10-31_iv/iv_fit.json"
SC = "/data/daq/slowcontrol/2026/hv_%s.csv.gz"
DVBD_DT = 0.060  # V/°C [S]
EXCL_Y = {"Y11", "Y13", "Y14"}

# ---------------------------------------------------------------- mayo
may = [r for r in csv.DictReader(l for l in open(MAYCSV) if not l.startswith("#"))]
MAYPTS = {g: [] for g in range(1, 9)}
for r in may:
    MAYPTS[int(r["hv_group"])].append((r["local_time_utc-4"], float(r["V0Set_V"]), float(r["VMon_V"]), float(r["IMon_A"]) * 1e6))
VM_MAY = {g: [p for p in MAYPTS[g] if p[0].startswith("2026-05-28")][0] for g in MAYPTS}

# ---------------------------------------------------------------- 1. octubre run 8
e = np.fromfile(REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work/hd64_oct8_hits.bin", dtype=xy.hits_dtype(64))
T0, T1 = float(e["t"].min()), float(e["t"].max())
acc = {c: {"V0Set": [], "VMon": [], "IMon": [], "status_not1": 0, "pw0": 0, "n": 0} for c in range(8)}
for day in ("20261022", "20261023"):
    for r in csv.DictReader(gzip.open(SC % day, "rt")):
        t = int(r["t_req_utc_ns"]) * 1e-9
        c = int(r["ch"])
        if not (T0 <= t <= T1) or int(r["slot"]) != 3 or c > 7:
            continue
        a = acc[c]
        a["n"] += 1
        if r["Pw"] != "1":
            a["pw0"] += 1
            continue
        a["status_not1"] += r["Status"] != "1"
        a["V0Set"].append(float(r["V0Set_V"]))
        a["VMon"].append(float(r["VMon_V"]))
        a["IMon"].append(float(r["IMon_A"]) * 1e6)
OCT = {}
rows1 = []
for c in range(8):
    a = acc[c]
    v, i = np.array(a["VMon"]), np.array(a["IMon"])
    OCT[c + 1] = {"V0Set_V": float(np.median(a["V0Set"])), "V0Set_distinct": sorted(set(a["V0Set"])), "VMon_V": float(np.median(v)),
                  "VMon_p16_p84": [float(np.percentile(v, 16)), float(np.percentile(v, 84))], "IMon_uA": float(np.median(i)),
                  "IMon_p16_p84": [float(np.percentile(i, 16)), float(np.percentile(i, 84))], "n_reads": a["n"], "n_pw0": a["pw0"], "n_status_not1": a["status_not1"]}
    o = OCT[c + 1]
    rows1.append([c, c + 1, o["n_reads"], o["n_pw0"], o["n_status_not1"], " ".join(f"{x:.2f}" for x in o["V0Set_distinct"]), f"{o['VMon_V']:.3f}",
                  f"{o['VMon_p16_p84'][0]:.3f}", f"{o['VMon_p16_p84'][1]:.3f}", f"{o['IMon_uA']:.4f}", f"{o['IMon_p16_p84'][0]:.4f}", f"{o['IMon_p16_p84'][1]:.4f}"])
with open(f"{OUT}/hv_oct_run8.csv", "w") as f:
    f.write(f"# hvmon log /data/daq/slowcontrol/2026/hv_2026102[23].csv.gz, slot 3, between the first and last trigger of October run 8 (UTC {T0:.1f} - {T1:.1f} s); "
            "medians over reads with Pw = 1 [M]\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hv_ch", "hv_group", "n_reads", "n_pw0", "n_status_not1", "V0Set_V_distinct", "VMon_median_V", "VMon_p16_V", "VMon_p84_V", "IMon_median_uA",
                "IMon_p16_uA", "IMon_p84_uA"])
    w.writerows(rows1)
DV = {g: VM_MAY[g][2] - OCT[g]["VMon_V"] for g in range(1, 9)}
with open(f"{OUT}/hv_delta.csv", "w") as f:
    f.write("# Delta V = VMon(May 2026-05-28 01:10 local, GECO screenshot) - VMon(October run 8, hvmon median) [M]; V0Set May assumed constant in the "
            "May good time [S]\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hv_group", "hv_ch", "V0Set_may_V", "VMon_may_V", "IMon_may_uA", "V0Set_oct_V", "VMon_oct_V", "IMon_oct_uA", "dV0Set_V", "dVMon_V", "dIMon_uA"])
    for g in range(1, 9):
        m, o = VM_MAY[g], OCT[g]
        w.writerow([g, g - 1, f"{m[1]:.2f}", f"{m[2]:.2f}", f"{m[3]:.3f}", f"{o['V0Set_V']:.2f}", f"{o['VMon_V']:.3f}", f"{o['IMon_uA']:.4f}",
                    f"{m[1] - o['V0Set_V']:.2f}", f"{DV[g]:.3f}", f"{m[3] - o['IMon_uA']:.4f}"])

# ---------------------------------------------------------------- 2. razones por grupo
rb = {r["id"]: r for r in csv.DictReader(open(R7 + "ratio_by_mppc.csv"))}
gb = {r["MPPC"]: r for r in csv.DictReader(open(R7 + "gain_by_mppc.csv"))}
items = []
for mid, r in rb.items():
    g = gb.get(mid)
    gr = float(g["G_ratio_mppc"]) if g and g["oct_G_is_median"] == "False" and int(r["ch_may"]) != 3 and int(r["ch_oct"]) != 3 else np.nan
    lr = float(r["ratio"]) if r["ratio"] not in ("", "nan") else np.nan
    items.append({"id": mid, "group": int(r["hv_group"]), "G_ratio": gr, "light_ratio": lr, "excl": mid in EXCL_Y})


def group_stats(key, use_excl):
    vals = {g: [math.log(it[key]) for it in items if it["group"] == g and np.isfinite(it[key]) and it[key] > 0 and not (use_excl and it["excl"])]
            for g in range(1, 9)}
    res_all = [x - np.median(v) for v in vals.values() if len(v) for x in v]
    s_pool = 1.4826 * float(np.median(np.abs(res_all)))
    out = {}
    for g, v in vals.items():
        if not v:
            continue
        m = float(np.median(v))
        sp = 1.4826 * float(np.median(np.abs(np.array(v) - m))) if len(v) > 1 else float("nan")
        out[g] = {"n": len(v), "ln_median": m, "ratio_median": math.exp(m), "ln_spread_group": sp, "ln_err": 1.2533 * s_pool / math.sqrt(len(v)),
                  "ratio_min": math.exp(min(v)), "ratio_max": math.exp(max(v))}
    return out, s_pool


def wls(x, y, s):
    x, y, s = map(np.asarray, (x, y, s))
    w = 1 / s ** 2
    X = np.vstack([np.ones_like(x), x]).T
    C = np.linalg.inv((X * w[:, None]).T @ X)
    p = C @ (X * w[:, None]).T @ y
    chi = float((((y - X @ p) ** 2) * w).sum())
    ndf = len(y) - 2
    return {"intercept": float(p[0]), "slope": float(p[1]), "intercept_err": float(math.sqrt(C[0, 0])), "slope_err": float(math.sqrt(C[1, 1])),
            "cov": float(C[0, 1]), "chi2": chi, "ndf": ndf}


def interp(f):
    """ΔV_bd = −b/k, ΔT = ΔV_bd/0.060; errores por propagación con la covarianza."""
    b, k = f["intercept"], f["slope"]
    if k == 0:
        return {}
    dvb = -b / k
    J = np.array([-1 / k, b / k ** 2])
    C = np.array([[f["intercept_err"] ** 2, f["cov"]], [f["cov"], f["slope_err"] ** 2]])
    e = math.sqrt(float(J @ C @ J))
    return {"OV_from_slope_V": 1 / k if k > 0 else None, "dVbd_may_minus_oct_V": dvb, "dVbd_err_V": e, "dT_may_minus_oct_C": dvb / DVBD_DT,
            "dT_err_C": e / DVBD_DT}


REG = {}
rows2 = []
for key in ("G_ratio", "light_ratio"):
    for use_excl in (True, False):
        st, sp = group_stats(key, use_excl)
        gs = sorted(st)
        x = [DV[g] for g in gs]
        y = [st[g]["ln_median"] for g in gs]
        s = [st[g]["ln_err"] for g in gs]
        fw = wls(x, y, s)
        fu = wls(x, y, [1.0] * len(gs))
        # sin pesos: error de los residuos
        sc = math.sqrt(fu["chi2"] / fu["ndf"]) if fu["ndf"] > 0 else float("nan")
        for k in ("intercept_err", "slope_err"):
            fu[k] *= sc
        fu["cov"] *= sc ** 2
        tag = f"{key}|{'excl_Y11_Y13_Y14' if use_excl else 'all'}"
        REG[tag] = {"pooled_ln_spread": sp, "groups": st, "weighted": fw, "weighted_interp": interp(fw), "unweighted": fu, "unweighted_interp": interp(fu)}
        for g in gs:
            rows2.append([key, "excl" if use_excl else "all", g, f"{DV[g]:.3f}", st[g]["n"], f"{st[g]['ratio_median']:.4f}", f"{st[g]['ln_spread_group']:.4f}",
                          f"{st[g]['ln_err']:.4f}", f"{st[g]['ratio_min']:.4f}", f"{st[g]['ratio_max']:.4f}"])
with open(f"{OUT}/ratio_vs_dv.csv", "w") as f:
    f.write("# per HV group (hod2018 MPPCs matched by MPPC, followup 7): median of the May/October gain ratio (gain_by_mppc.csv) and primary-p.e. light ratio "
            "(ratio_by_mppc.csv); ln_spread = 1.4826 MAD of ln ratio within the group; ln_err = 1.2533 * pooled within-group spread / sqrt(n); "
            "excl = followup-7 main exclusion (bars Y11, Y13, Y14); all = no exclusion; dVMon from hv_delta.csv\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["quantity", "selection", "hv_group", "dVMon_V", "n_mppc", "ratio_median", "ln_spread_group", "ln_err", "ratio_min", "ratio_max"])
    w.writerows(rows2)

# ---------------------------------------------------------------- 3. V_bd de mayo
ivf = json.load(open(IV))["groups"]
good_n = [ivf[str(g)]["b_free"]["n"] for g in range(1, 9) if not ivf[str(g)]["b_free"]["at_grid_edge"]]
n_med = float(np.median(good_n))


def fit_may(pts, i0, n, vb_grid, n_free=False):
    v = np.array([p[2] for p in pts])
    i = np.array([p[3] for p in pts])
    s = np.sqrt(0.002 ** 2 + (0.05 * i) ** 2)
    w = 1 / s ** 2
    y = i - i0
    best, prof = None, []
    ngrid = np.arange(1.0, 10.0001, 0.01) if n_free else np.array([n])
    for vb in vb_grid:
        x = np.clip(v - vb, 0, None)
        F = np.where(x[None, :] > 0, np.clip(x, 1e-12, None)[None, :] ** ngrid[:, None], 0.0)  # (n, puntos), vectorizado en n
        den = (w * F * F).sum(1)
        a = np.where(den > 0, (w * F * y).sum(1) / np.where(den > 0, den, 1), -1.0)
        chi = np.where(a >= 0, (w * (y[None, :] - a[:, None] * F) ** 2).sum(1), np.inf)
        k = int(np.argmin(chi))
        prof.append(float(chi[k]))
        if np.isfinite(chi[k]) and (best is None or chi[k] < best[0]):
            best = (float(chi[k]), float(a[k]), float(vb), float(ngrid[k]))
    if best is None:
        return None
    prof = np.array(prof)
    ok = vb_grid[prof <= best[0] + 1.0]
    edge = abs(best[2] - vb_grid.min()) < 2e-3 or abs(best[2] - vb_grid.max()) < 2e-3
    return {"chi2": best[0], "ndf": len(v) - (3 if n_free else 2), "a_uA": best[1], "V_bd_V": best[2], "n": best[3],
            "V_bd_err_V": [float(ok.min() - best[2]), float(ok.max() - best[2])], "at_grid_edge": bool(edge),
            "OV_at_may28_V": float(max(p[2] for p in pts) - best[2])}


VB = {}
rows3 = []
grid = np.arange(60.0, 65.4, 0.001)
grid_free = np.arange(60.0, 65.4, 0.005)  # n libre: paso 5 mm V
for g in range(1, 9):
    G_ = ivf[str(g)]
    i0 = G_["I0_floor_uA"]
    if not G_["b_free"]["at_grid_edge"]:
        ref, src = G_["b_free"], "b_free"
    elif "b_fixed" in G_:
        ref, src = G_["b_fixed"], "b_fixed"
    else:
        ref, src = None, "none"
    n = ref["n"] if ref else n_med
    pts = MAYPTS[g]
    above = [p for p in pts if p[3] > i0]
    vmin_ok = min(p[2] for p in pts)
    gg = grid[grid < max(p[2] for p in pts) - 0.01]
    fx = fit_may(pts, i0, n, gg)
    ff = fit_may(pts, i0, n, grid_free[grid_free < max(p[2] for p in pts) - 0.01], n_free=True)
    VB[g] = {"I0_oct_uA": i0, "n_source": src if ref else f"median of determined groups ({n_med:.3f}) [S]", "n_fixed": n,
             "V_bd_oct_V": ref["V_bd_V"] if ref else None, "V_bd_oct_err_V": ref["V_bd_err_V"] if ref else None,
             "V_bd_oct_b_free_V": G_["b_free"]["V_bd_V"], "b_free_oct_at_grid_edge": G_["b_free"]["at_grid_edge"],
             "may_points_above_I0": len(above), "fixed_n": fx, "free_n": ff}
    d = (fx["V_bd_V"] - ref["V_bd_V"]) if (fx and ref) else None
    VB[g]["dVbd_may_minus_oct_V"] = d
    rows3.append([g, g - 1, f"{i0:.3f}", src, f"{n:.3f}", len(above),
                  f"{fx['V_bd_V']:.3f}" if fx else "", f"{fx['V_bd_err_V'][0]:+.3f}/{fx['V_bd_err_V'][1]:+.3f}" if fx else "",
                  f"{fx['chi2']:.2f}/{fx['ndf']}" if fx else "", int(fx["at_grid_edge"]) if fx else "", f"{fx['OV_at_may28_V']:.3f}" if fx else "",
                  f"{ref['V_bd_V']:.3f}" if ref else "", f"{ref['V_bd_err_V'][0]:+.3f}/{ref['V_bd_err_V'][1]:+.3f}" if ref else "",
                  f"{d:+.3f}" if d is not None else "",
                  f"{ff['V_bd_V']:.3f}" if ff else "", f"{ff['n']:.3f}" if ff else "", f"{ff['chi2']:.3g}" if ff else ""])
with open(f"{OUT}/vbd_may.csv", "w") as f:
    f.write("# May breakdown voltage per HV group: fit of the three May GECO I-V points with the followup-8 model (b) I = I0 + a (V - V_bd)^n, I0 and n "
            "fixed to the October fit of the group (b_free, or b_fixed if b_free is at the grid edge, else the median n of the determined groups [S]); "
            "a, V_bd free (1 dof); sigma_I = sqrt(0.002^2 + (0.05 I)^2) uA; V_bd error = profile dchi2 = 1; variant with n free (0 dof)\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hv_group", "hv_ch", "I0_oct_uA", "oct_ref_fit", "n_fixed", "may_points_above_I0", "V_bd_may_V", "V_bd_may_err_V", "chi2_ndf", "at_grid_edge",
                "OV_may28_V", "V_bd_oct_V", "V_bd_oct_err_V", "dVbd_may_minus_oct_V", "V_bd_may_free_n_V", "n_free", "chi2_free_n"])
    w.writerows(rows3)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(x) for x in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if math.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


json.dump(clean({"oct_run8_window_utc_s": [T0, T1], "oct": OCT, "may_28": {g: VM_MAY[g] for g in VM_MAY}, "dV": DV, "regressions": REG, "vbd": VB,
                 "dVbd_dT_V_per_C": DVBD_DT}), open(f"{OUT}/hv_may_oct.json", "w"), indent=1)
print(f"octubre run 8: UTC {T0:.1f}-{T1:.1f}")
for r in rows1:
    print("oct", r)
for g in range(1, 9):
    print(f"grupo {g}: VMon mayo {VM_MAY[g][2]:.2f} oct {OCT[g]['VMon_V']:.3f} -> ΔV {DV[g]:+.3f} V; IMon mayo {VM_MAY[g][3]:.3f} oct {OCT[g]['IMon_uA']:.3f} uA")
for tag, R in REG.items():
    print(f"== {tag} (dispersión conjunta ln {R['pooled_ln_spread']:.4f})")
    for g, s in R["groups"].items():
        print(f"   grupo {g}: n {s['n']} mediana {s['ratio_median']:.4f} (ln disp. {s['ln_spread_group']:.4f}, err ln {s['ln_err']:.4f}) min {s['ratio_min']:.4f} max {s['ratio_max']:.4f}")
    for k in ("weighted", "unweighted"):
        f_, i_ = R[k], R[k + "_interp"]
        print(f"   {k}: pendiente {f_['slope']:.4f} ± {f_['slope_err']:.4f} /V, ordenada {f_['intercept']:.4f} ± {f_['intercept_err']:.4f}, χ² {f_['chi2']:.2f}/{f_['ndf']}; "
              f"ΔV_bd(mayo−oct) {i_.get('dVbd_may_minus_oct_V', float('nan')):+.3f} ± {i_.get('dVbd_err_V', float('nan')):.3f} V, "
              f"ΔT {i_.get('dT_may_minus_oct_C', float('nan')):+.2f} ± {i_.get('dT_err_C', float('nan')):.2f} °C")
for r in rows3:
    print("vbd", r)
