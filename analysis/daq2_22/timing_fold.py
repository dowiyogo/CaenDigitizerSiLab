#!/usr/bin/python3 -I
# timing_fold.py — DAQ2-22 etapa 3: tiempo del muón (impactos > 20 p.e. de hod2018, timing_interp.C) plegado módulo 32 ns.
# En run 8 el tiempo en el registro toma tres valores separados 32 ns con poblaciones parecidas, sin relación con ttt_raw mod 4/8
# (el bit 0 de ttt_raw vale siempre 1): es la colocación del registro respecto del trigger, de a 32 ns. Plegando módulo 32 ns se
# quita y queda la respuesta temporal trigger − muón. Ajuste de máxima verosimilitud en el círculo (gaussianas con imágenes ±32 ns):
#   modelo 1: una gaussiana (μ, σ); modelo 2: dos gaussianas con σ común (μ1, μ2 = μ1 + Δ, fracción w) + fondo plano b
#   (impactos de ruido > 20 p.e. y pulsos mal medidos). Errores por la curvatura (diferencias finitas).
# Además, por tira (X e Y, orden r4) el μ plegado con el modelo 1 sobre una ventana ±8 ns del pico global.
# Uso: timing_fold.py <dir con timing_run008.csv, timing_run010.csv> <mapping r4 json>
import json, math, sys
import numpy as np

D, MAP = sys.argv[1], sys.argv[2]
P = 32.0
CH = {c["dig_ch"]: c for c in json.load(open(MAP))["channels"]}


def wrap(x):
    return (x + P / 2) % P - P / 2


def pdf(t, mus, sig, ws, b):
    v = np.full_like(t, b / P)
    for m, w in zip(mus, ws):
        g = sum(np.exp(-0.5 * ((t - m + k * P) / sig) ** 2) for k in (-2, -1, 0, 1, 2)) / (sig * math.sqrt(2 * math.pi))
        v += (1 - b) * w * g
    return v


def nll1(t, p):
    m, s, b = p
    if s <= 0.3 or not 0 <= b < 1:
        return 1e30
    return -np.log(pdf(t, [m], s, [1.0], b)).sum()


def nll2(t, p):
    m, d, s, w, b = p
    if s <= 0.3 or not 0 < w < 1 or not 0 <= b < 1 or not 2 < d < 30:
        return 1e30
    return -np.log(pdf(t, [m, m + d], s, [w, 1 - w], b)).sum()


def box(t, a, W, sig):
    """Caja [a, a + W] convolucionada con una gaussiana σ, en el círculo de 32 ns (normalizada)."""
    from math import erf
    v = np.zeros_like(t)
    for k in (-2, -1, 0, 1, 2):
        u = t + k * P
        v += 0.5 * (np.vectorize(erf)((u - a) / (sig * math.sqrt(2))) - np.vectorize(erf)((u - a - W) / (sig * math.sqrt(2))))
    return v / W


def nll3(t, p):
    # modelo 3: pico gaussiano (μ, σ, fracción w) + meseta (caja desde μ + a0, ancho W, suavizada con σ) + fondo plano b
    m, s, w, a0, W, b = p
    if s <= 0.3 or not 0 < w < 1 or not 0 <= b < 1 or w + b >= 1 or not 1 < W < 26 or not -20 < a0 < 20:
        return 1e30
    g = sum(np.exp(-0.5 * ((t - m + k * P) / s) ** 2) for k in (-2, -1, 0, 1, 2)) / (s * math.sqrt(2 * math.pi))
    v = w * g + (1 - w - b) * box(t, m + a0, W, s) + b / P
    return -np.log(v).sum()


def minimize(f, p0, steps):
    p = np.array(p0, float)
    st = np.array(steps, float)
    best = f(p)
    for _ in range(400):
        imp = False
        for i in range(len(p)):
            for sg in (1, -1):
                q = p.copy(); q[i] += sg * st[i]
                v = f(q)
                if v < best:
                    best, p, imp = v, q, True
        if not imp:
            st /= 2
            if st.max() < 1e-4:
                break
    # errores por la hessiana numérica
    n = len(p); H = np.zeros((n, n)); h = np.maximum(np.abs(p) * 1e-3, 1e-3)
    for i in range(n):
        for j in range(n):
            pp = p.copy(); pp[i] += h[i]; pp[j] += h[j]; fpp = f(pp)
            pm = p.copy(); pm[i] += h[i]; pm[j] -= h[j]; fpm = f(pm)
            mp = p.copy(); mp[i] -= h[i]; mp[j] += h[j]; fmp = f(mp)
            mm = p.copy(); mm[i] -= h[i]; mm[j] -= h[j]; fmm = f(mm)
            H[i, j] = (fpp - fpm - fmp + fmm) / (4 * h[i] * h[j])
    try:
        err = np.sqrt(np.diag(np.linalg.inv(H)))
    except np.linalg.LinAlgError:
        err = np.full(n, np.nan)
    return p, err, best


R = {}
for run in ("run008", "run010"):
    a = np.genfromtxt(f"{D}/timing_{run}.csv", delimiter=",", names=True)
    out = {"n": int(len(a))}
    for tk in ("t_par", "t_cfd"):
        t = a[tk] * 16.0
        # centro: máximo del histograma plegado (1 ns)
        h, e = np.histogram(t % P, bins=32, range=(0, P))
        c0 = e[np.argmax(h)] + 0.5
        tw = wrap(t - c0)
        p1, e1, l1 = minimize(lambda p: nll1(tw, p), [0.0, 3.0, 0.1], [0.5, 0.5, 0.05])
        best2 = None
        for d0 in (8, 12, 16, 20, 24):
            p2, e2, l2 = minimize(lambda p: nll2(tw, p), [0.0, d0, 2.5, 0.5, 0.05], [0.5, 1.0, 0.3, 0.05, 0.02])
            if best2 is None or l2 < best2[2]:
                best2 = (p2, e2, l2)
        p2, e2, l2 = best2
        best3 = None
        for a00 in (-12, -6, 2, 6):
            for W0 in (6, 10, 14):
                p3, e3, l3 = minimize(lambda p: nll3(tw, p), [0.0, 2.0, 0.4, a00, W0, 0.05], [0.5, 0.3, 0.05, 1.0, 1.0, 0.02])
                if best3 is None or l3 < best3[2]:
                    best3 = (p3, e3, l3)
        p3, e3, l3 = best3
        hf, ef = np.histogram(tw, bins=32, range=(-P / 2, P / 2))
        out[tk] = {"centre_ns_mod32": float(c0), "model1": {"mu": p1[0], "mu_err": e1[0], "sigma": p1[1], "sigma_err": e1[1], "flat": p1[2], "nll": l1},
                   "model2": {"mu1": p2[0], "delta": p2[1], "delta_err": e2[1], "sigma": p2[2], "sigma_err": e2[2], "w1": p2[3], "w1_err": e2[3], "flat": p2[4], "nll": l2},
                   "model3": {"mu": p3[0], "sigma": p3[1], "sigma_err": e3[1], "w_peak": p3[2], "w_peak_err": e3[2], "box_start_rel": p3[3], "box_start_err": e3[3],
                              "box_width": p3[4], "box_width_err": e3[4], "flat": p3[5], "nll": l3},
                   "two_minus_one_dnll": l1 - l2, "three_minus_one_dnll": l1 - l3, "hist_folded_1ns": hf.tolist()}
        # por tira (modelo 1, ventana ±8 ns alrededor del pico del modelo 1 global): media y error
        if tk == "t_par":
            strips = {}
            for c in sorted(set(a["ch"].astype(int))):
                m = a["ch"].astype(int) == c
                x = wrap(tw[m] - p1[0])
                x = x[np.abs(x) < 8]
                if len(x) >= 8:
                    strips[int(c)] = {"plane": CH[c]["plane"], "strip": CH[c]["strip"], "n": int(len(x)), "mean_ns": float(x.mean()), "err_ns": float(x.std() / math.sqrt(len(x)))}
            out["strips"] = strips
            for pl in ("X", "Y"):
                v = [(s["strip"], s["mean_ns"], s["err_ns"]) for s in strips.values() if s["plane"] == pl]
                if len(v) >= 3:
                    u = np.array([x[0] for x in v], float); y = np.array([x[1] for x in v]); w = 1 / np.array([x[2] for x in v]) ** 2
                    A = np.vstack([np.ones_like(u), u - 7.5]).T
                    C = np.linalg.inv((A * w[:, None]).T @ A); coef = C @ (A * w[:, None]).T @ y
                    chi = float((((y - A @ coef) ** 2) * w).sum())
                    m0 = float((y * w).sum() / w.sum()); chi0 = float((((y - m0) ** 2) * w).sum())
                    out[f"strip_slope_{pl}"] = {"slope_ns_per_strip": float(coef[1]), "slope_err": float(math.sqrt(C[1, 1])), "chi2_const": chi0, "ndf_const": len(v) - 1,
                                                "chi2_lin": chi, "ndf_lin": len(v) - 2}
    R[run] = out
    for tk in ("t_par", "t_cfd"):
        o = out[tk]
        print(f"{run} {tk}: 1 gauss σ {o['model1']['sigma']:.2f}±{o['model1']['sigma_err']:.2f} ns (plano {o['model1']['flat']:.2f}); "
              f"2 gauss Δ {o['model2']['delta']:.2f}±{o['model2']['delta_err']:.2f} ns, σ {o['model2']['sigma']:.2f}±{o['model2']['sigma_err']:.2f}, "
              f"w1 {o['model2']['w1']:.2f}±{o['model2']['w1_err']:.2f}, plano {o['model2']['flat']:.2f}; Δ(−lnL) = {o['two_minus_one_dnll']:.1f}")
        m3 = o["model3"]
        print(f"   pico + meseta: σ {m3['sigma']:.2f}±{m3['sigma_err']:.2f}, w_pico {m3['w_peak']:.2f}±{m3['w_peak_err']:.2f}, meseta desde {m3['box_start_rel']:+.1f}±{m3['box_start_err']:.1f} ns, "
              f"ancho {m3['box_width']:.1f}±{m3['box_width_err']:.1f} ns, plano {m3['flat']:.2f}; Δ(−lnL) vs 1 gauss = {o['three_minus_one_dnll']:.1f}")
    for pl in ("X", "Y"):
        if f"strip_slope_{pl}" in out:
            print(f"   tiras {pl}: {out[f'strip_slope_{pl}']}")
json.dump(R, open(f"{D}/timing_fold.json", "w"), indent=1, default=float)
