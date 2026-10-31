#!/usr/bin/python3 -I
# iv_fit.py — seguimiento 2026-10-31 (8), punto 2: ajuste de la corriente contra el voltaje por grupo HV (CH00-07) con los niveles estables de
# iv_levels.py (variante 10 s).
# Modelos (VMon en V, IMon en µA):
#   (a) pedido:            I = a (V − V_bd)^n, sobre los niveles por encima del piso (−0.5 V, nominal, +0.25 V): 3 parámetros;
#   (b) con piso de fuga:  I = I_0 + a (V − V_bd)^n para V > V_bd e I_0 debajo, I_0 = IMon a −1.5 V (por debajo de la ruptura) [supuesto: la
#       corriente a −1.5 V es fuga sin multiplicación], sobre todos los niveles (−1.0, −0.5, nominal, +0.25);
#   en los dos, V_bd libre y fijado al V_bd de dos puntos del punto 1 (mediana del grupo, solo grupos 2, 6 y 8).
# χ² con σ_I = √(sd² + 0.002² + (0.05·I)²) µA (sd de la muestra, resolución de IMon y 5 % de deriva entre días: el nivel nominal de un mismo
# canal cambia 4-10 % entre el 22, el 23 y el 30 de octubre); mínimo por rejilla
# en V_bd (paso 1 mV) y n (paso 0.005) con a por mínimos cuadrados; error de V_bd y n del perfil Δχ² = 1.
# Un ajuste con V_bd o n en el borde de la rejilla se marca "at_grid_edge" (no determinado).
# Corriente esperada a nominal + 0.25, 0.5, 0.75 y 1.0 V (VMon nominal del 2026-10-23 + ΔV) frente a I0Set = 15 µA.
# Salida: <out>/iv_fit.json, iv_fit.csv. Uso: iv_fit.py <dir iv> <dir vbd>
import csv, json, math, statistics as st, sys

import numpy as np

D, V = sys.argv[1], sys.argv[2]
L = [r for r in json.load(open(D + "/iv_levels.json")) if r["variant"] == "10s"]
VBD = json.load(open(V + "/vbd.json"))["summary"]["two_point"]
vbd_fixed = {}
for t in VBD:
    vbd_fixed.setdefault(t["hv_group"], []).append(t["V_bd_V"])
vbd_fixed = {g: st.median(v) for g, v in vbd_fixed.items()}
DV_OUT = (0.25, 0.5, 0.75, 1.0)


def model(v, i0, a, vb, n):
    x = np.clip(v - vb, 0, None)
    return i0 + a * x ** n


def fit(v, i, s, i0, vb_grid, n_grid):
    """mínimo de χ² en la rejilla (V_bd, n), a por mínimos cuadrados; vectorizado en n."""
    y = i - (i0 or 0.0)
    w = 1 / s ** 2
    best = None
    prof_vb, prof_n = {}, np.full(len(n_grid), np.inf)
    for vb in vb_grid:
        x = v - vb
        if i0 is None and np.any(x <= 0):
            continue
        F = np.where(x[None, :] > 0, np.clip(x, 1e-12, None)[None, :] ** n_grid[:, None], 0.0)  # (n, puntos)
        den = (w * F * F).sum(1)
        a = np.where(den > 0, (w * F * y).sum(1) / np.where(den > 0, den, 1), -1)
        chi = np.where(a > 0, (w * (y[None, :] - a[:, None] * F) ** 2).sum(1), np.inf)
        k = int(np.argmin(chi))
        if not np.isfinite(chi[k]):
            continue
        prof_vb[round(float(vb), 4)] = float(chi[k])
        prof_n = np.minimum(prof_n, chi)
        if best is None or chi[k] < best[0]:
            best = (float(chi[k]), float(a[k]), float(vb), float(n_grid[k]))
    if best is None:
        return None
    chi, a, vb, n = best

    def err(keys, vals, x0):
        ok = [k for k, c in zip(keys, vals) if c <= chi + 1]
        return (min(ok) - x0, max(ok) - x0) if ok else (None, None)
    return {"chi2": chi, "a_uA": a, "V_bd_V": vb, "n": n,
            "V_bd_err_V": err(list(prof_vb), list(prof_vb.values()), vb) if len(vb_grid) > 1 else (0.0, 0.0),
            "n_err": err(list(n_grid), list(prof_n), n)}


out = {"groups": {}, "I0Set_uA": 15.0}
rows = []
for g in range(1, 9):
    ch = g - 1
    P = [r for r in L if r["hv_ch"] == ch]
    v = np.array([r["VMon_median_V"] for r in P])
    i = np.array([r["IMon_median_A"] * 1e6 for r in P])
    sd = np.array([r["IMon_sd_A"] * 1e6 for r in P])
    s = np.sqrt(sd ** 2 + 0.002 ** 2 + (0.05 * i) ** 2)
    floor_pts = [r for r in P if r["scan"] == "scan23" and abs(r["V0Set_V"] - max(x["V0Set_V"] for x in P if x["scan"] == "scan23") + 1.5) < 1e-6]
    i0 = floor_pts[0]["IMon_median_A"] * 1e6
    vfloor = floor_pts[0]["VMon_median_V"]
    v_nom23 = [r for r in P if r["scan"] == "scan23"]
    v_nom = max(r["VMon_median_V"] for r in v_nom23)
    above = v > vfloor + 0.6  # −0.5 V, nominal, +0.25 V (sin −1.0 ni −1.5)
    allp = v > vfloor + 0.01  # sin el piso mismo
    res = {"hv_ch": ch, "points": [{"scan": r["scan"], "V0Set_V": r["V0Set_V"], "VMon_V": r["VMon_median_V"], "IMon_uA": r["IMon_median_A"] * 1e6, "n": r["n"]} for r in P],
           "I0_floor_uA": i0, "V_floor_VMon_V": vfloor, "V_nom_VMon_V": v_nom}
    ng = np.arange(1.0, 8.0001, 0.005)
    vbg = np.arange(vfloor - 1.0, v.max() - 0.3, 0.001)
    vbg_a = vbg[vbg < v[above].min() - 1e-4]
    res["a_free"] = fit(v[above], i[above], s[above], None, vbg_a, ng) if above.sum() >= 3 else None
    res["b_free"] = fit(v[allp], i[allp], s[allp], i0, vbg, ng)
    if g in vbd_fixed:
        res["V_bd_fixed_V"] = vbd_fixed[g]
        res["a_fixed"] = fit(v[above], i[above], s[above], None, np.array([vbd_fixed[g]]), ng)
        res["b_fixed"] = fit(v[allp], i[allp], s[allp], i0, np.array([vbd_fixed[g]]), ng)
    for k in ("a_free", "b_free", "a_fixed", "b_fixed"):
        f = res.get(k)
        if f:
            f["ndf"] = int((above if k.startswith("a") else allp).sum()) - (3 if k.endswith("free") else 2)
            f["OV_nominal_V"] = v_nom - f["V_bd_V"]
            grid = (vbg_a if k == "a_free" else vbg) if k.endswith("free") else np.array([f["V_bd_V"]])
            f["at_grid_edge"] = bool(k.endswith("free") and (abs(f["V_bd_V"] - grid.min()) < 2e-3 or abs(f["V_bd_V"] - grid.max()) < 2e-3)) or \
                f["n"] <= 1.0 + 1e-6 or f["n"] >= 8.0 - 1e-6
            f["n_distinct_V0Set"] = len({round(p_["V0Set_V"], 2) for p_, m_ in zip(P, (above if k.startswith("a") else allp)) if m_})
            i0_ = 0.0 if k.startswith("a") else i0
            f["I_expected_uA"] = {str(dv): float(model(np.array([v_nom + dv]), i0_, f["a_uA"], f["V_bd_V"], f["n"])[0]) for dv in (0.0,) + DV_OUT}
            f["first_dV_over_15uA"] = next((dv for dv in np.arange(0, 1.0001, 0.01) if model(np.array([v_nom + dv]), i0_, f["a_uA"], f["V_bd_V"], f["n"])[0] > 15.0), None)
    out["groups"][g] = res
    for k in ("a_free", "b_free", "a_fixed", "b_fixed"):
        f = res.get(k)
        if f:
            e = f["V_bd_err_V"]
            print(f"grupo {g} (CH{ch:02d}) {k:7s}: V_bd {f['V_bd_V']:.3f} V ({e[0] if e[0] is None else f'{e[0]:+.3f}'}/{e[1] if e[1] is None else f'{e[1]:+.3f}'}) "
                  f"OV_nom {f['OV_nominal_V']:.3f} V, n {f['n']:.2f} ({f['n_err'][0]:+.2f}/{f['n_err'][1]:+.2f}), χ² {f['chi2']:.2f}/{f['ndf']} | I(nom) "
                  f"{f['I_expected_uA']['0.0']:.2f}, +0.25 {f['I_expected_uA']['0.25']:.2f}, +0.5 {f['I_expected_uA']['0.5']:.2f}, +0.75 {f['I_expected_uA']['0.75']:.2f}, "
                  f"+1.0 {f['I_expected_uA']['1.0']:.2f} µA; > 15 µA desde +{f['first_dV_over_15uA'] if f['first_dV_over_15uA'] is None else round(f['first_dV_over_15uA'], 2)} V"
                  f"{' | BORDE DE LA REJILLA' if f['at_grid_edge'] else ''} | {f['n_distinct_V0Set']} V0Set distintos")
            rows.append({"hv_group": g, "model": k, "V_bd_V": f["V_bd_V"], "V_bd_err_lo": e[0], "V_bd_err_hi": e[1], "OV_nominal_V": f["OV_nominal_V"],
                         "n": f["n"], "chi2": f["chi2"], "ndf": f["ndf"], **{f"I_{dv}V_uA": f["I_expected_uA"][str(dv)] for dv in DV_OUT},
                         "first_dV_over_15uA": f["first_dV_over_15uA"], "at_grid_edge": f["at_grid_edge"], "n_distinct_V0Set": f["n_distinct_V0Set"]})
json.dump(out, open(D + "/iv_fit.json", "w"), indent=1, default=float)
with open(D + "/iv_fit.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
