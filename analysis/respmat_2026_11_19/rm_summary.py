#!/usr/bin/python3 -I
# rm_summary.py — seguimiento 2026-11-19 (matriz de respuesta), resumen de las inyecciones de rm_inject.C para mayo y octubre.
#  - Matriz de respuesta por capa (XA, XB, YA, YB; canales de la capa en el orden de tiras de cada campaña, ch3 fuera, suma de canales con peso igual por
#    inyección [supuesto]): P(A_meas bin | n), bins de 0.25 p.e. medidos de 0 a 60 + desborde, n = 0 … 60 (n = 0 = ruido solo).
#  - Linealidad: media(A)/n y RMS(A) contra n; P(A_meas ≥ n + 1 | n) (migración hacia arriba de al menos 1 p.e. por el ruido; sin diafonía inyectada).
#  - Distribución de n = 0 en la ventana signal y en la de control [99,108) (mismos segmentos).
#  - Comprobación de 1 p.e.: pico de n = 1 (gaussiana en ±0.2 p.e. alrededor del máximo en [0.5, 1.3], bins de 0.01) por canal, contra
#    (a) la calibración de los espectros (calibration_<camp>.csv: G = separación de picos; 1 p.e. = 1 p.e. medido por definición; el archivo no da anchura ni
#        error; error de G estimado aquí de la recta de las posiciones de picos de su fuente: hd64_may_v8_peaks.csv, pe_calibration_2026_10_c1.json);
#    (b) el pico de 1 p.e. de los pulsos oscuros reales en la ventana signal de los mismos segmentos (n = 0, mismo ajuste);
#    (c) el pico de 1 p.e. del estimador de pico en los pulsos aislados pre-trigger (prep_<camp>.json);
#    y la separación de picos inyectada (n = 2 − n = 1) contra 1 p.e. medido (= G).
#    Modelos de pico: "ge" = gaussiana + exponencial en [0.6, 1.5] p.e. (el mismo para el pico oscuro de n = 0 y para n = 1); "g" = gaussiana sola en ±0.2
#    alrededor del máximo (n = 1 en [0.5, 1.3], n = 2 en [1.2, 2.4]); ganancia inyectada = pendiente de media(A) entre n = 20 y 60.
# Salida (sufijo _<variante> si se da una): <out>/response_<camp>.csv, linearity_<camp>.csv, linearity_channels_<camp>.csv, noise_n0_<camp>.csv, check_1pe.csv,
#         response_matrix.json.
# Uso: rm_summary.py <out dir> [variante]   (sin variante: plantilla "ls", work/<camp>_resp.bin; "gnorm": work/<camp>_gnorm_resp.bin)
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(os.environ.get("ROOTSYS", os.path.expanduser("~/root6")), "lib"))  # python3 -I ignora PYTHONPATH
import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)
OUT = sys.argv[1]
VAR = sys.argv[2] if len(sys.argv) > 2 else ""
SUF = f"_{VAR}" if VAR else ""
W = OUT + "/work"
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
ORDER = json.load(open(PK + "/stage3.json"))["strip_orders"]
C1 = os.path.expanduser("~/.local/opt/daqv2/share/calibration/pe_calibration_2026_10_c1.json")
NN, NB, NF = 61, 241, 600
EDG = np.round(np.arange(0, 60.0001, 0.25), 4)
LAYERS = ("XA", "XB", "YA", "YB")
DEAD = 3


def layers(camp):
    out = {}
    for probe in ("X", "Y"):
        for L in ("A", "B"):
            out[probe + L] = [int(c) for i, c in enumerate(ORDER[camp]["order"][f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
    return out


def load(camp):
    v = np.fromfile(f"{W}/{camp}{SUF}_resp.bin", dtype="<f8")
    o = 0
    out = {}
    for name, shp in (("H", (32, NN, NB)), ("F", (32, 4, NF)), ("S", (32, NN, 3)), ("Hc", (32, NB)), ("Fc", (32, NF))):
        k = int(np.prod(shp))
        out[name] = v[o:o + k].reshape(shp)
        o += k
    assert o == len(v)
    return out


def gfit(h, lo=0.5, hi=1.3, half=0.2):
    """Gaussiana en ±half alrededor del máximo de h (bins de 0.01 desde 0) dentro de [lo, hi)."""
    th = ROOT.TH1D("hf", "", NF, 0.0, 6.0)
    for i, v in enumerate(h):
        th.SetBinContent(i + 1, float(v))
        th.SetBinError(i + 1, math.sqrt(max(float(v), 1.0)))
    a, b = int(round(lo / 0.01)), int(round(hi / 0.01))
    mx = (a + int(np.argmax(h[a:b])) + 0.5) * 0.01
    r = th.Fit("gaus", "QSN0L", "", mx - half, mx + half)
    out = None
    if int(r) == 0:
        out = (r.Parameter(1), r.ParError(1), abs(r.Parameter(2)), r.ParError(2))
    th.Delete()
    return out


def gefit(h, lo=0.6, hi=1.5, mu=1.0):
    """Gaussiana + exponencial en [lo, hi] (bins de 0.01 desde 0), máxima verosimilitud de Poisson."""
    th = ROOT.TH1D("hf", "", NF, 0.0, 6.0)
    for i, v in enumerate(h):
        th.SetBinContent(i + 1, float(v))
    f = ROOT.TF1("fge", "gaus(0)+expo(3)", lo, hi)
    f.SetParameters(max(float(h[int(mu / 0.01)]), 1.0), mu, 0.18, 5.0, -3.0)
    f.SetParLimits(1, lo, hi)
    f.SetParLimits(2, 0.03, 0.6)
    r = th.Fit(f, "QSN0L", "", lo, hi)
    out = (r.Parameter(1), r.ParError(1), abs(r.Parameter(2)), r.ParError(2)) if int(r) == 0 else None
    th.Delete()
    return out


def g_err(camp):
    """Error de G por canal desde la recta de las posiciones de picos de la fuente de la calibración (índice consecutivo)."""
    pos = {}
    if camp == "may":
        for l in open(PK + "/hd64_may_v8_peaks.csv"):
            f = l.strip().split(",")
            if f[0] != "ch" and f[3] in ("ok", "changed"):
                pos[int(f[0])] = [float(x) for x in f[5].split()]
    else:
        for e in json.load(open(C1))["channels"]:
            if e["ch"] < 32 and isinstance(e.get("G_mV"), (int, float)):
                pos[e["ch"]] = list(e["peak_positions_mV"])
    out = {}
    for c, p in pos.items():
        i = np.arange(len(p))
        if len(p) < 3:
            continue
        A = np.vstack([i, np.ones_like(i)]).T
        coef, rs, _, _ = np.linalg.lstsq(A, np.array(p), rcond=None)
        res = np.array(p) - A @ coef
        s2 = (res ** 2).sum() / (len(p) - 2)
        out[c] = (float(coef[0]), float(math.sqrt(s2 / ((i - i.mean()) ** 2).sum())), float(coef[1]))
    return out


summary = {}
chk_rows = []
for camp in ("may", "oct"):
    R = load(camp)
    prep = json.load(open(f"{OUT}/prep_{camp}.json"))
    G = np.array([float(v) for v in open(f"{W}/G_{camp}.txt").read().split(",")])
    LY = layers(camp)
    lay_of = {c: L for L, cs in LY.items() for c in cs}
    rows, lrows, crow, nrows = [], [], [], []
    sm = {}
    for L, cs in LY.items():
        H = R["H"][cs].sum(0)
        S = R["S"][cs].sum(0)
        Nn = S[:, 0]
        P = H / Nn[:, None]
        for k in range(NN):
            for b in range(NB):
                pl, ph = (EDG[b], EDG[b + 1]) if b < NB - 1 else (60.0, float("inf"))
                rows.append([camp, L, k, f"{pl:.2f}", "inf" if ph == float("inf") else f"{ph:.2f}", int(H[k, b]), f"{P[k, b]:.6e}"])
        mean = S[:, 1] / Nn
        rms = np.sqrt(np.maximum(S[:, 2] / Nn - mean ** 2, 0))
        for k in range(NN):
            lrows.append([camp, L, k, int(Nn[k]), f"{mean[k]:.5f}", f"{mean[k] / k:.5f}" if k else "", f"{rms[k]:.5f}", f"{rms[k] / math.sqrt(Nn[k]):.5f}"])
        # P(A ≥ n + 1 | n)
        pup = {}
        for k in (1, 2, 3, 5, 10, 20, 40):
            j = int(round((k + 1) / 0.25))
            cnt = float(H[k, j:].sum())
            pup[str(k)] = {"p": cnt / Nn[k], "err": math.sqrt(cnt * (1 - cnt / Nn[k])) / Nn[k], "count": cnt, "N": float(Nn[k])}
        # también hacia abajo: P(A < n − 1 | n)
        pdn = {}
        for k in (5, 10, 20, 40):
            j = int(round((k - 1) / 0.25))
            cnt = float(H[k, :j].sum())
            pdn[str(k)] = {"p": cnt / Nn[k], "count": cnt, "N": float(Nn[k])}
        # ruido solo (n = 0): señal y control
        H0, Hc = H[0], R["Hc"][cs].sum(0)
        n0 = {}
        for x in (0.5, 1.0, 1.5, 2.5, 3.0, 5.0):
            j = int(round(x / 0.25))
            n0[str(x)] = {"P_signal_ge": float(H0[j:].sum() / H0.sum()), "P_control_ge": float(Hc[j:].sum() / Hc.sum())}
        for b in range(NB):
            pl, ph = (EDG[b], EDG[b + 1]) if b < NB - 1 else (60.0, float("inf"))
            nrows.append([camp, L, f"{pl:.2f}", "inf" if ph == float("inf") else f"{ph:.2f}", int(H0[b]), f"{H0[b] / H0.sum():.6e}", int(Hc[b]), f"{Hc[b] / Hc.sum():.6e}"])
        # pico de 1 p.e. de la capa (n = 1) y de n = 2
        F1, F2, F0 = R["F"][cs, 1].sum(0), R["F"][cs, 2].sum(0), R["F"][cs, 0].sum(0)
        f1, f2 = gfit(F1), gfit(F2, 1.2, 2.4)
        f0, e1 = gefit(F0), gefit(F1)
        sm[L] = {"channels": cs, "injections_per_n_min": float(Nn.min()), "injections_per_n_max": float(Nn.max()),
                 "mean_over_n": {str(k): float(mean[k] / k) for k in (1, 2, 5, 10, 20, 40, 60)},
                 "mean_over_n_err": {str(k): float(rms[k] / math.sqrt(Nn[k]) / k) for k in (1, 2, 5, 10, 20, 40, 60)},
                 "rms": {str(k): float(rms[k]) for k in (0, 1, 2, 5, 10, 20, 40, 60)},
                 "rms_over_sqrt_n": {str(k): float(rms[k] / math.sqrt(k)) for k in (1, 5, 10, 20, 40, 60)},
                 "P_ge_n_plus_1": pup, "P_lt_n_minus_1": pdn, "noise_n0": n0, "mean_n0": float(mean[0]),
                 "peak_n1_gauss": f1, "peak_n2_gauss": f2, "peak_n1_gauss_expo": e1, "peak_dark_n0_gauss_expo": f0,
                 "slope_mean_20_60": float((mean[60] - mean[20]) / 40.0)}
        print(f"{camp} {L}: canales {cs}, inyecciones por n {Nn.min():.0f}-{Nn.max():.0f}; media/n "
              + " ".join(f"n={k}:{mean[k] / k:.4f}" for k in (1, 5, 10, 20, 40)) + f"; P(≥n+1) n=5 {pup['5']['p']:.4f} n=10 {pup['10']['p']:.4f}; "
              f"pico n=1 (g) {f1}; n=2 (g) {f2}; n=1 (ge) {e1}; oscuro n=0 (ge) {f0}", flush=True)
    # por canal: linealidad y comprobación de 1 p.e.
    ge = g_err(camp)
    for c in range(32):
        if c == DEAD:
            continue
        S = R["S"][c]
        mean = S[:, 1] / S[:, 0]
        rms = np.sqrt(np.maximum(S[:, 2] / S[:, 0] - mean ** 2, 0))
        for k in range(NN):
            crow.append([camp, c, lay_of[c], k, int(S[k, 0]), f"{mean[k]:.5f}", f"{mean[k] / k:.5f}" if k else "", f"{rms[k]:.5f}"])
        e1, e0 = gefit(R["F"][c, 1]), gefit(R["F"][c, 0])
        f1, f2 = gfit(R["F"][c, 1]), gfit(R["F"][c, 2], 1.2, 2.4)
        slope = (mean[60] - mean[20]) / 40.0
        slope_err = math.sqrt(rms[60] ** 2 / S[60, 0] + rms[20] ** 2 / S[20, 0]) / 40.0
        t = prep["templates"][str(c)]
        g = ge.get(c)
        nan4 = (float("nan"),) * 4
        chk_rows.append([camp, c, lay_of[c], f"{G[c]:.5f}", f"{g[1] / G[c]:.5f}" if g else "",
                         *(f"{v:.5f}" for v in (e1 or nan4)), *(f"{v:.5f}" for v in (e0 or nan4)),
                         *(f"{v:.5f}" for v in (f1 or nan4)[:2]), *(f"{v:.5f}" for v in (f2 or nan4)[:2]),
                         f"{(f2[0] - f1[0]) if (f1 and f2) else float('nan'):.5f}", f"{slope:.5f}", f"{slope_err:.5f}",
                         f"{t['peak_est_1pe_pos_over_G']:.5f}", f"{t['peak_est_1pe_sigma_over_G']:.5f}", f"{t['A1_LS_over_G']:.5f}", f"{t['A1_LS_over_G_err']:.5f}"])
    with open(f"{OUT}/response_{camp}{SUF}.csv", "w") as f:
        f.write(f"# hod2018 {camp}: response matrix P(A_meas bin | n) of the tagged-muon amplitude estimator (aS/G in measured p.e.), from n x single-p.e. template injected "
                "into real noise segments of the signal window (no crosstalk); bins of spectra_*.csv (0.25 p.e., 0-60 + overflow); count = injections in the bin\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["campaign", "layer", "n", "pe_lo", "pe_hi", "count", "prob"])
        w.writerows(rows)
    with open(f"{OUT}/linearity_{camp}{SUF}.csv", "w") as f:
        f.write(f"# hod2018 {camp}: mean and RMS of A_meas (measured p.e.) versus injected n per layer; err_mean = RMS/sqrt(N)\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["campaign", "layer", "n", "N_injections", "mean_A", "mean_A_over_n", "rms_A", "err_mean"])
        w.writerows(lrows)
    with open(f"{OUT}/linearity_channels_{camp}{SUF}.csv", "w") as f:
        f.write(f"# hod2018 {camp}: mean and RMS of A_meas (measured p.e.) versus injected n per channel\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["campaign", "ch", "layer", "n", "N_injections", "mean_A", "mean_A_over_n", "rms_A"])
        w.writerows(crow)
    with open(f"{OUT}/noise_n0_{camp}{SUF}.csv", "w") as f:
        f.write(f"# hod2018 {camp}: n = 0 (noise only) distribution of A_meas in the signal window [49,58) and in the control window [99,108) of the same segments\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["campaign", "layer", "pe_lo", "pe_hi", "count_signal", "prob_signal", "count_control", "prob_control"])
        w.writerows(nrows)
    summary[camp] = {"layers": sm, "template_summary": {k: prep[k] for k in ("estimator_check", "events_selected", "muon_time")}}
with open(f"{OUT}/check_1pe{SUF}.csv", "w") as f:
    f.write(f"# hod2018: single-p.e. check, template variant {VAR or 'ls'}. G_cal = calibration_<camp>.csv (peak spacing; the file has no width or error); "
            "G_err_rel = relative slope error of the peak positions of its source; inj1_ge / dark_ge = Gaussian + exponential fit in [0.6,1.5] measured p.e. of the "
            "n = 1 injections / of the real dark pulses in the signal window of the same segments (n = 0); inj1_g, inj2_g = Gaussian in +-0.2 around the maximum; "
            "spacing = inj2_g - inj1_g; slope = (mean(A) at n = 60 - at n = 20)/40; iso = 1 p.e. peak of the peak estimator on isolated pre-trigger dark pulses; "
            "A1_LS = 1 p.e. peak of the least-squares amplitude (template normalization of variant ls)\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "ch", "layer", "G_cal_mV", "G_err_rel", "inj1_ge_pos", "inj1_ge_pos_err", "inj1_ge_sigma", "inj1_ge_sigma_err", "dark_ge_pos", "dark_ge_pos_err",
                "dark_ge_sigma", "dark_ge_sigma_err", "inj1_g_pos", "inj1_g_pos_err", "inj2_g_pos", "inj2_g_pos_err", "spacing", "slope_20_60", "slope_err",
                "iso_pos", "iso_sigma", "A1_LS_over_G", "A1_LS_over_G_err"])
    w.writerows(chk_rows)
# resumen de la comprobación por campaña
for camp in ("may", "oct"):
    rr = [r for r in chk_rows if r[0] == camp]
    def col(i):
        return np.array([float(r[i]) if r[i] not in ("", "nan") else np.nan for r in rr])
    p1, e1, s1, es1 = col(5), col(6), col(7), col(8)
    d, ed, sd, esd = col(9), col(10), col(11), col(12)
    gerr = col(4)
    sp, sl, esl = col(17), col(18), col(19)
    pull_pos = (p1 - d) / np.sqrt(e1 ** 2 + ed ** 2)
    pull_sig = (s1 - sd) / np.sqrt(es1 ** 2 + esd ** 2)
    pull_sl = (sl - 1.0) / np.sqrt(esl ** 2 + gerr ** 2)
    summary[camp]["check_1pe"] = {
        "median_inj1_ge_pos": float(np.nanmedian(p1)), "median_dark_ge_pos": float(np.nanmedian(d)), "median_inj1_ge_sigma": float(np.nanmedian(s1)),
        "median_dark_ge_sigma": float(np.nanmedian(sd)), "median_inj1_g_pos": float(np.nanmedian(col(13))), "median_inj2_g_pos": float(np.nanmedian(col(15))),
        "median_spacing": float(np.nanmedian(sp)), "median_slope_20_60": float(np.nanmedian(sl)), "median_G_err_rel": float(np.nanmedian(gerr)),
        "median_iso_pos": float(np.nanmedian(col(20))), "median_iso_sigma": float(np.nanmedian(col(21))),
        "n_fits_ok_inj1_ge": int(np.isfinite(p1).sum()), "n_fits_ok_dark_ge": int(np.isfinite(d).sum()),
        "pull_pos_inj_vs_dark": [round(float(v), 2) for v in pull_pos], "n_abs_pull_pos_lt2": int((np.abs(pull_pos) < 2).sum()),
        "pull_sigma_inj_vs_dark": [round(float(v), 2) for v in pull_sig], "n_abs_pull_sigma_lt2": int((np.abs(pull_sig) < 2).sum()),
        "pull_slope_vs_G": [round(float(v), 2) for v in pull_sl], "n_abs_pull_slope_lt2": int((np.abs(pull_sl) < 2).sum()),
        "n_channels": len(rr)}
    print(camp, "comprobación 1 p.e.:", {k: v for k, v in summary[camp]["check_1pe"].items() if not k.startswith("pull")}, flush=True)
    for k in ("pull_pos_inj_vs_dark", "pull_sigma_inj_vs_dark", "pull_slope_vs_G"):
        print(" ", k, summary[camp]["check_1pe"][k])
json.dump(summary, open(f"{OUT}/response_matrix{SUF}.json", "w"), indent=1)
