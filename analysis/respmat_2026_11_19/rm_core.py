#!/usr/bin/python3 -I
# rm_core.py — seguimiento 2026-11-19 (matriz de respuesta), núcleo y cola baja de la respuesta, y su origen en la ventana signal.
#  1. Por campaña, variante (ls, gnorm, tcut, gnorm_tcut: plantilla ls o gnorm, tiempos medidos o restringidos a 48-55), capa y n = 5, 10, 20, 40, 60, desde
#     work/<camp>[_<variante>]_resp.bin: gaussiana (bins de 0.25 p.e., verosimilitud de Poisson) en [0.75 n, 1.25 n] → media/n y σ/n del núcleo (sin valor si
#     1.25 n > 60, donde el desborde corta la ventana); P(A < 0.5 n) (cola baja) y la fracción de inyecciones fuera de [0.75 n, 1.25 n].
#  2. Estimador sin ruido sobre la plantilla (mediana de los canales de la campaña, variante ls normalizada a pico 1): para cada muestra entera j del tiempo CFD
#     (43 … 57) y 200 fases uniformes, y_i = T(i − j − φ); pulso = máximo local (y_i ≥ y_{i±1}) con i en [49,58); amplitud = el mayor; 0 si no hay. Media de la
#     amplitud / pico de T y fracción de fases sin pulso, junto a la probabilidad de j de los muones etiquetados (prep_<camp>.json, promedio de las 4 capas).
# Salida: <out>/core_tail.csv, <out>/window_edge.csv, <out>/core_tail.json.
# Uso: rm_core.py <out dir>
import csv, json, math, os, sys

import numpy as np

sys.path.insert(0, os.path.join(os.environ.get("ROOTSYS", os.path.expanduser("~/root6")), "lib"))  # python3 -I ignora PYTHONPATH
import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)
OUT = sys.argv[1]
W = OUT + "/work"
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
ORDER = json.load(open(REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work/stage3.json"))["strip_orders"]
NN, NB = 61, 241
rows, wrows, res = [], [], {}
for camp in ("may", "oct"):
    LY = {}
    for probe in ("X", "Y"):
        for L in ("A", "B"):
            LY[probe + L] = [int(c) for i, c in enumerate(ORDER[camp]["order"][f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != 3]
    for var, suf in (("ls", ""), ("gnorm", "_gnorm"), ("tcut", "_tcut"), ("gnorm_tcut", "_gnorm_tcut")):
        H = np.fromfile(f"{W}/{camp}{suf}_resp.bin", "<f8")[:32 * NN * NB].reshape(32, NN, NB)
        for L, cs in LY.items():
            h = H[cs].sum(0)
            for n in (5, 10, 20, 40, 60):
                tot = h[n].sum()
                th = ROOT.TH1D("hc", "", NB - 1, 0.0, 60.0)
                for b in range(NB - 1):
                    th.SetBinContent(b + 1, float(h[n, b]))
                lo, hi = 0.75 * n, min(1.25 * n, 59.99)
                r = th.Fit("gaus", "QSN0L", "", lo, hi)
                mu, emu, sg, esg = (r.Parameter(1), r.ParError(1), abs(r.Parameter(2)), r.ParError(2)) if int(r) == 0 and 1.25 * n <= 60 else (float("nan"),) * 4
                th.Delete()
                jl, jlo, jhi = int(round(0.5 * n / 0.25)), int(round(lo / 0.25)), int(round(hi / 0.25))
                plow = float(h[n, :jl].sum() / tot)
                pout = float((h[n, :jlo].sum() + h[n, jhi:].sum()) / tot)
                rows.append([camp, var, L, n, int(tot), f"{mu / n:.5f}", f"{emu / n:.5f}", f"{sg / n:.5f}", f"{esg / n:.5f}", f"{plow:.5f}", f"{pout:.5f}"])
                res.setdefault(camp, {}).setdefault(var, {}).setdefault(L, {})[str(n)] = {"core_mean_over_n": mu / n, "core_sigma_over_n": sg / n,
                                                                                           "core_sigma_over_sqrt_n": sg / math.sqrt(n), "P_A_lt_half_n": plow,
                                                                                           "P_outside_core_window": pout}
                print(f"{camp} {var} {L} n={n}: núcleo media/n {mu / n:.4f} σ/n {sg / n:.4f} (σ/√n {sg / math.sqrt(n):.3f}); P(A<n/2) {plow:.4f}; fuera de "
                      f"[0.75n,1.25n] {pout:.4f}", flush=True)
    # 2. estimador sin ruido contra la muestra entera del tiempo
    T = np.fromfile(f"{W}/{camp}_tmpl.bin", "<f8").reshape(32, 333)
    grid = -1.5 + 0.125 * np.arange(333)
    sh = np.median(np.array([T[c] / T[c].max() for c in range(32) if T[c].any()]), axis=0)
    prep = json.load(open(f"{OUT}/prep_{camp}.json"))
    pj = np.zeros(144)
    for v in prep["muon_time"].values():
        for k, p in v["p_integer_sample"].items():
            pj[int(k)] += p / 4.0
    ph = (np.arange(200) + 0.5) / 200
    tot_loss = 0.0
    for j in range(43, 58):
        amp, none = [], 0
        for f in ph:
            i = np.arange(47, 60)
            y = np.interp(i - j - f, grid, sh, left=0.0, right=0.0)
            best = 0.0
            for m in range(2, 11):  # i = 49 … 57
                if y[m] >= y[m - 1] and y[m] >= y[m + 1]:
                    best = max(best, y[m])
            amp.append(best)
            none += best == 0.0
        amp = np.array(amp)
        wrows.append([camp, j, f"{pj[j]:.5f}", f"{amp.mean():.5f}", f"{amp.min():.5f}", f"{amp.max():.5f}", f"{none / len(ph):.4f}"])
        tot_loss += pj[j] * (amp.mean() < 0.5)
        res.setdefault(camp, {}).setdefault("window_edge", {})[str(j)] = {"p_muon": float(pj[j]), "mean_amp_over_peak": float(amp.mean()), "frac_no_pulse": none / len(ph)}
    res[camp]["p_muon_j_with_mean_amp_lt_half"] = float(tot_loss)
    print(f"{camp}: probabilidad de j de los muones con amplitud media < 0.5 del pico (sin ruido): {tot_loss:.4f}", flush=True)
    for r_ in wrows:
        if r_[0] == camp:
            print("   ", r_, flush=True)
with open(f"{OUT}/core_tail.csv", "w") as f:
    f.write("# hod2018: core (Gaussian fit in [0.75n,1.25n], 0.25 p.e. bins) and low tail of the response per campaign, template variant and layer\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "variant", "layer", "n", "N_injections", "core_mean_over_n", "core_mean_err", "core_sigma_over_n", "core_sigma_err", "P_A_lt_half_n",
                "P_outside_0.75n_1.25n"])
    w.writerows(rows)
with open(f"{OUT}/window_edge.csv", "w") as f:
    f.write("# hod2018: noiseless peak estimator (local maximum of y in [49,58)) on the median single-p.e. shape versus the integer sample j of the CFD time, "
            "200 uniform phases; p_muon = probability of j for tagged muons (mean of the 4 layers)\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "j", "p_muon", "mean_amp_over_peak", "min_amp_over_peak", "max_amp_over_peak", "frac_no_pulse"])
    w.writerows(wrows)
json.dump(res, open(f"{OUT}/core_tail.json", "w"), indent=1)
