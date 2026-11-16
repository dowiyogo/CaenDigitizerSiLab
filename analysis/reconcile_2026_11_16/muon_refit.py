#!/usr/bin/python3 -I
# muon_refit.py — seguimiento 2026-11-16, punto 3: reajuste del espectro de los muones etiquetados (29-54 mm) con Landau ⊗ Poisson ⊗ Borel–Tanner
# ⊗ resolución, el modelo de xtalk.py (analysis/xtalk_2026_11_15, 85e8498) sin cambios, sobre el mismo espectro (muon_spectra.csv de
# outputs/hod2018_xtalk/2026-11-15_borel, unidades N = (aS − off)/Gf, bins de 1 en [6, 80)) y la misma resolución (xtalk.json, conversion).
# Ajustes (χ² con error² = mismo + escala²·mezclados; Nelder–Mead con varios arranques):
#   R  referencia: μ_m, μ_o, w común, p_m, p_o libres (el mínimo de xtalk.py se buscó en rejilla; aquí con los 5 libres);
#   a  p_oct = 0.74·p_mayo [hoja de datos: 27 % a 1.2 V, 20 % a 1.0 V]: μ_m, μ_o, w común, p_m;
#   b  anchura de Landau libre por campaña: μ_m, μ_o, w_m, w_o, p_m, p_o;
#   b0 = b con p_m = p_o.
# Covarianza: 2·H⁻¹ con H la hessiana numérica de χ² en el mínimo (diferencias centrales); correlación entre p y w.
# Salida: <out>/muon_refit.csv, muon_refit.json. Uso: muon_refit.py <xtalk dir> <out dir>
import csv, json, math, sys

import numpy as np

XT, OUT = sys.argv[1], sys.argv[2]
conv = json.load(open(XT + "/xtalk.json"))["conversion"]
spec = {"may": [], "oct": []}
for r in csv.DictReader(open(XT + "/muon_spectra.csv")):
    spec["may"].append((float(r["may_net"]), float(r["may_err"])))
    spec["oct"].append((float(r["oct_net"]), float(r["oct_err"])))
SPEC = {k: (np.array([a for a, _ in v]), np.array([b for _, b in v])) for k, v in spec.items()}
EDG = np.arange(-0.5, 300.5, 1.0)
assert len(SPEC["may"][0]) == len(EDG) - 1

# ---- modelo (copiado de xtalk.py sin cambios)
_t = np.linspace(1e-6, 60, 60001)
_dt = _t[1] - _t[0]


def landau_std(lam):
    lam = np.atleast_1d(lam)
    out = np.empty(len(lam))
    for i0 in range(0, len(lam), 50):
        L = lam[i0:i0 + 50, None]
        out[i0:i0 + 50] = (np.exp(-_t * np.log(_t) - L * _t) * np.sin(np.pi * _t)).sum(1) * _dt / np.pi
    return out


XG = np.linspace(0.02, 12.0, 600)
LAM = np.linspace(-4.0, 120.0, 6201)
LANDAU_TAB = landau_std(LAM)


def landau(x, w):
    return np.interp((x - 1.0) / w - 0.22278, LAM, LANDAU_TAB, left=0.0, right=0.0) / w


NMAXP, NMAXN = 260, 300
_n = np.arange(NMAXP + 1)
_lgn = np.array([math.lgamma(k + 1) for k in range(NMAXN + 1)])
_lgp = np.array([math.lgamma(k + 1) for k in range(NMAXP + 1)])


def bt_matrix(p):
    N = np.arange(NMAXN + 1)[:, None]
    n = _n[None, :]
    M = np.zeros((NMAXN + 1, NMAXP + 1))
    if p <= 1e-9:
        for k in range(min(NMAXN, NMAXP) + 1):
            M[k, k] = 1.0
        return M
    k = N - n
    ok = (k >= 0) & (n >= 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log(np.where(ok, n, 1)) - np.log(np.where(N > 0, N, 1)) - p * N + k * np.log(p * np.where(N > 0, N, 1)) - _lgn[np.clip(k, 0, NMAXN)]
    M[ok] = np.exp(lg[ok])
    M[0, 0] = 1.0
    return M


def smear_matrix(s0, s1):
    N = np.arange(NMAXN + 1, dtype=float)
    sig = np.sqrt(s0 * s0 + N * s1 * s1) + 1e-6
    k = np.arange(len(EDG) - 1, dtype=float)
    erf = np.vectorize(math.erf)
    lo = erf((k[:, None] - 0.5 - N[None, :]) / (math.sqrt(2) * sig[None, :]))
    hi = erf((k[:, None] + 0.5 - N[None, :]) / (math.sqrt(2) * sig[None, :]))
    return 0.5 * (hi - lo)


SM = {c: smear_matrix(conv[c]["s0_rel"], conv[c]["s1_rel"]) for c in ("may", "oct")}
SEL = slice(6, 80)
_bt = {}


def model(camp, mu, w, p):
    key = round(p, 7)
    if key not in _bt:
        if len(_bt) > 300:
            _bt.clear()
        _bt[key] = bt_matrix(p)
    wx = landau(XG, w)
    wx = wx / wx.sum()
    lam = mu * XG
    pois = np.exp(_n[:, None] * np.log(lam[None, :]) - lam[None, :] - _lgp[:, None])
    return SM[camp] @ (_bt[key] @ (pois @ wx))


def chi2c(camp, mu, w, p):
    if mu <= 0.2 or w <= 0.02 or not (0 <= p < 0.9):
        return 1e12
    m = model(camp, mu, w, p)[SEL]
    d, e = SPEC[camp][0][SEL], SPEC[camp][1][SEL]
    e = np.where(e > 0, e, 1.0)
    m = m * d.sum() / max(m.sum(), 1e-30)
    return float((((d - m) / e) ** 2).sum())


def nelder_mead(f, x0, step, it=3000, tol=1e-8):
    n = len(x0)
    pts = [np.array(x0, float)] + [np.array(x0, float) + np.eye(n)[i] * step[i] for i in range(n)]
    val = [f(p) for p in pts]
    for _ in range(it):
        o = np.argsort(val)
        pts, val = [pts[i] for i in o], [val[i] for i in o]
        if abs(val[-1] - val[0]) < tol * (1 + abs(val[0])):
            break
        c = np.mean(pts[:-1], axis=0)
        xr = c + (c - pts[-1])
        fr = f(xr)
        if fr < val[0]:
            xe = c + 2 * (c - pts[-1])
            fe = f(xe)
            pts[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            pts[-1], val[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c)
            fc = f(xc)
            if fc < val[-1]:
                pts[-1], val[-1] = xc, fc
            else:
                pts = [pts[0]] + [pts[0] + 0.5 * (p - pts[0]) for p in pts[1:]]
                val = [val[0]] + [f(p) for p in pts[1:]]
    i = int(np.argmin(val))
    return pts[i], val[i]


def fit(f, starts, step):
    best = None
    for s in starts:
        x, v = nelder_mead(f, s, step)
        x, v = nelder_mead(f, x, [0.2 * a for a in step])
        x, v = nelder_mead(f, x, [0.05 * a for a in step])
        if best is None or v < best[1]:
            best = (x, v)
    return best


def hess(f, x, h):
    n = len(x)
    H = np.zeros((n, n))
    f0 = f(x)
    for i in range(n):
        for j in range(i, n):
            ei, ej = np.eye(n)[i] * h[i], np.eye(n)[j] * h[j]
            if i == j:
                H[i, i] = (f(x + ei) - 2 * f0 + f(x - ei)) / h[i] ** 2
            else:
                H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h[i] * h[j])
    return H


def cov_corr(f, x, h, names):
    H = hess(f, np.asarray(x, float), h)
    try:
        C = 2 * np.linalg.inv(H)
        d = np.sqrt(np.clip(np.diag(C), 0, None))
        R = C / np.outer(d, d)
        return {"err": dict(zip(names, map(float, d))), "corr": {f"{names[i]}|{names[j]}": float(R[i, j]) for i in range(len(names)) for j in range(i + 1, len(names))},
                "hessian_posdef": bool(np.all(np.linalg.eigvalsh(H) > 0))}
    except np.linalg.LinAlgError:
        return {"err": None, "corr": None, "hessian_posdef": False}


NBIN = 2 * (SEL.stop - SEL.start)
res, rows = {}, []
# R: referencia, 5 libres
fR = lambda x: chi2c("may", x[0], x[2], x[3]) + chi2c("oct", x[1], x[2], x[4])  # noqa: E731
xR, vR = fit(fR, [(6.275, 3.965, 0.0858, 0.55, 0.52), (7.25, 5.5, 0.15, 0.43, 0.32), (10.0, 7.0, 0.2, 0.2, 0.2)], (1.0, 1.0, 0.03, 0.05, 0.05))
res["R"] = {"x": dict(zip(["mu_may", "mu_oct", "w", "p_may", "p_oct"], map(float, xR))), "chi2": vR, "ndf": NBIN - 5 - 2,
            **cov_corr(fR, xR, [0.05, 0.05, 0.003, 0.005, 0.005], ["mu_may", "mu_oct", "w", "p_may", "p_oct"])}
print("R", res["R"], flush=True)
# a: p_oct = 0.74 p_may
fA = lambda x: chi2c("may", x[0], x[2], x[3]) + chi2c("oct", x[1], x[2], 0.74 * x[3])  # noqa: E731
xA, vA = fit(fA, [(xR[0], xR[1], xR[2], xR[3]), (7.25, 5.5, 0.15, 0.43), (10.0, 7.0, 0.2, 0.2)], (1.0, 1.0, 0.03, 0.05))
res["a"] = {"x": dict(zip(["mu_may", "mu_oct", "w", "p_may"], map(float, xA))), "p_oct": 0.74 * float(xA[3]), "chi2": vA, "ndf": NBIN - 4 - 2,
            "dchi2_vs_R": vA - vR, **cov_corr(fA, xA, [0.05, 0.05, 0.003, 0.005], ["mu_may", "mu_oct", "w", "p_may"])}
print("a", res["a"], flush=True)
# b: w libre por campaña
fB = lambda x: chi2c("may", x[0], x[2], x[4]) + chi2c("oct", x[1], x[3], x[5])  # noqa: E731
xB, vB = fit(fB, [(xR[0], xR[1], xR[2], xR[2], xR[3], xR[4]), (7.25, 5.5, 0.15, 0.15, 0.43, 0.32), (10.0, 7.0, 0.2, 0.2, 0.2, 0.2),
                  (12.0, 7.6, 0.22, 0.22, 0.01, 0.01)], (1.0, 1.0, 0.03, 0.03, 0.05, 0.05))
res["b"] = {"x": dict(zip(["mu_may", "mu_oct", "w_may", "w_oct", "p_may", "p_oct"], map(float, xB))), "chi2": vB, "ndf": NBIN - 6 - 2, "dchi2_vs_R": vB - vR,
            **cov_corr(fB, xB, [0.05, 0.05, 0.003, 0.003, 0.005, 0.005], ["mu_may", "mu_oct", "w_may", "w_oct", "p_may", "p_oct"])}
print("b", res["b"], flush=True)
# b0: w libre por campaña, p común
fB0 = lambda x: chi2c("may", x[0], x[2], x[4]) + chi2c("oct", x[1], x[3], x[4])  # noqa: E731
xB0, vB0 = fit(fB0, [(xB[0], xB[1], xB[2], xB[3], 0.5 * (xB[4] + xB[5])), (7.25, 5.5, 0.15, 0.15, 0.4), (12.0, 7.6, 0.22, 0.22, 0.01)],
               (1.0, 1.0, 0.03, 0.03, 0.05))
res["b0"] = {"x": dict(zip(["mu_may", "mu_oct", "w_may", "w_oct", "p"], map(float, xB0))), "chi2": vB0, "ndf": NBIN - 5 - 2, "dchi2_vs_b": vB0 - vB}
print("b0", res["b0"], flush=True)
# perfiles 1D de p en el ajuste b (p de una campaña fijo, todo lo demás libre), para ver la degeneración p–w
prof = []
for camp, ip in (("may", 4), ("oct", 5)):
    for pv in np.round(np.arange(0.0, 0.81, 0.1), 2):
        def fP(x, pv=pv, ip=ip):
            y = list(x[:4]) + ([pv, x[4]] if ip == 4 else [x[4], pv])
            return fB(y)
        x0 = list(xB[:4]) + [xB[5] if ip == 4 else xB[4]]
        xp, vp = fit(fP, [x0], (1.0, 1.0, 0.03, 0.03, 0.05))
        prof.append([camp, pv, f"{vp:.3f}", f"{vp - vB:.3f}", f"{xp[2]:.4f}", f"{xp[3]:.4f}"])
        print("perfil b", camp, pv, round(vp, 3), flush=True)
res["profile_b"] = prof
with open(OUT + "/muon_refit.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["fit", "parameters", "chi2", "ndf", "dchi2"])
    w.writerow(["R (5 libres, w común)", json.dumps(res["R"]["x"]), f"{vR:.3f}", res["R"]["ndf"], "0"])
    w.writerow(["a (p_oct = 0.74 p_may)", json.dumps(res["a"]["x"]), f"{vA:.3f}", res["a"]["ndf"], f"{vA - vR:.3f}"])
    w.writerow(["b (w por campaña)", json.dumps(res["b"]["x"]), f"{vB:.3f}", res["b"]["ndf"], f"{vB - vR:.3f}"])
    w.writerow(["b0 (w por campaña, p común)", json.dumps(res["b0"]["x"]), f"{vB0:.3f}", res["b0"]["ndf"], f"{vB0 - vR:.3f}"])
    w.writerow([])
    w.writerow(["profile_b_campaign", "p_fixed", "chi2", "dchi2_vs_b", "w_may", "w_oct"])
    w.writerows(prof)
json.dump(res, open(OUT + "/muon_refit.json", "w"), indent=1)
