#!/usr/bin/python3 -I
# iso_xtalk.py — seguimiento 2026-11-18: diafonía de hod2018 con pulsos oscuros aislados (iso_extract.C) y cota sobre juntas de epoxi despegadas.
# Modo prep <out>: lee work/{may_v8,oct_run8}_iso.bin; mayo en tiempo bueno de v8; escribe los pares {uint8 ch; float a} de cada variante en
#   work/pairs_<camp>_<variante>.bin para fit_iso.C. Variantes: pre_w100, pre_w150, pre_w250 = región pre-trigger k ∈ [15, 33], sin otro pulso > 0.5 G
#   del canal a |Δk| ≤ 6, 9, 15 muestras (±96, ±144, ±240 ns: el múltiplo de 16 ns ≤ 100, 150, 250 ns); post_w150 = región posterior k ∈ [73, 128],
#   |Δk| ≤ 9, y el canal sin pulso > 0.5 G en la ventana signal de ese evento.
# Modo summary <out>: junta work/fit_<camp>_<variante>.csv y calcula, por canal y campaña:
#   p, F = 1/(1 − p), P(N ≥ 2) = 1 − e^{−p}, χ²/ndf, fracción de pulsos conservados = aislados / todos (a > 0.5 G, misma región);
#   apilamiento residual esperado = R · 16 ns, R = tasa de pulsos > 0.5 G en la región (todos los pulsos);
#   contaminación de centelleo real: fracción de pulsos aislados con un pulso ≥ 3 G a |Δk| ≤ 2 (32 ns) en una barra vecina del mismo plano (índice ±1,
#   ±2) o en el plano ortogonal, en el mismo evento, menos la misma fracción con el evento siguiente (mezclado), por clase de amplitud.
#   Medianas de la estación, p_oct/p_mayo (por MPPC, barras emparejadas del seguimiento 7, y cociente de medianas), comparación con los muones.
#   Juntas de epoxi (punto 2): razones oct/mayo por canal con la selección 29-54 mm (outputs/hod2018_ratio_reconcile/2026-11-16/r1_C.csv, variante
#   thr3_far); varianza extra τ² por máxima verosimilitud de r_i ~ N(μ, σ_i² + τ²), límite superior al 95 % CL con −2Δln L = 2.71 (una cola);
#   número de canales con una pérdida relativa ≥ 26 % [Geant4 FU18, junta despegada]: modelo de k canales con razón 0.74·μ y el resto μ, los k
#   elegidos para maximizar L, μ libre; límite superior de k al 95 % CL con −2Δln L = 2.71 respecto del mejor k.
# Uso: iso_xtalk.py prep|summary <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

MODE, OUT = sys.argv[1], sys.argv[2]
WK = OUT + "/work"
F10 = REPO + "/outputs/hod2018_mip_tagged/2026-11-01_tagged/work"
R7 = REPO + "/outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv"
RC = REPO + "/outputs/hod2018_ratio_reconcile/2026-11-16/r1_C.csv"
DT = np.dtype([("ev", "<u4"), ("ch", "u1"), ("reg", "u1"), ("k", "<u2"), ("a", "<f4"), ("isod", "u1"), ("flags", "u1"), ("pad", "<u2")])
PRE = {"may": "may_v8", "oct": "oct_run8"}
VARS = {"pre_w100": (0, 6), "pre_w150": (0, 9), "pre_w250": (0, 15), "post_w150": (1, 9)}
NSAMP = {0: 19, 1: 56}
DEAD = 3
LIVE = [c for c in range(32) if c != DEAD]
G = {camp: np.array([float(r["G_mV"]) for r in csv.DictReader(open(f"{F10}/{pre}_dark.csv"))]) for camp, pre in PRE.items()}


def load(camp):
    d = np.fromfile(f"{WK}/{PRE[camp]}_iso.bin", dtype=DT)
    t = np.fromfile(f"{WK}/{PRE[camp]}_iso_t.bin", dtype="<f8")
    m = mc.good_time(t, "v8") if camp == "may" else np.ones(len(t), bool)
    return d[m[d["ev"]]], int(m.sum())


def select(d, camp, v):
    reg, w = VARS[v]
    g = G[camp][d["ch"]]
    base = (d["reg"] == reg) & (d["a"] > 0.5 * g)
    sel = (d["reg"] == reg) & (d["isod"] > w)
    if reg == 1:
        base &= (d["flags"] & 16) == 0
        sel &= (d["flags"] & 16) == 0
    return base, sel


if MODE == "prep":
    info = {}
    for camp in PRE:
        d, nev = load(camp)
        info[camp] = {"events": nev}
        for v in VARS:
            _, sel = select(d, camp, v)
            s = d[sel]
            out = np.empty(len(s), dtype=np.dtype([("ch", "u1"), ("a", "<f4")]))
            out["ch"], out["a"] = s["ch"], s["a"]
            out.tofile(f"{WK}/pairs_{camp}_{v}.bin")
            info[camp][v] = int(len(s))
            print(camp, v, "pulsos aislados (≥ 0.3 G):", len(s), flush=True)
    json.dump(info, open(WK + "/prep.json", "w"), indent=1)
    sys.exit(0)

# ---------------------------------------------------------------- summary
fits = {}
for camp in PRE:
    for v in VARS:
        for r in csv.DictReader(open(f"{WK}/fit_{camp}_{v}.csv")):
            fits[(camp, v, r["variant"], int(r["ch"]))] = r
rows, st = [], {}
cont = {}
for camp in PRE:
    d, nev = load(camp)
    g = G[camp][d["ch"]]
    for v in VARS:
        reg, w = VARS[v]
        base, sel = select(d, camp, v)
        nb = np.bincount(d["ch"][base], minlength=32)
        nk = np.bincount(d["ch"][sel & (d["a"] > 0.5 * g)], minlength=32)
        # tasa de pulsos > 0.5 G en la región (para la región posterior, por evento sin pulso en la ventana signal: aproximación con todos los eventos)
        R = nb / (nev * NSAMP[reg] * 16e-9)
        for fv in ("A", "D"):
            ps, chis, kept, pile = [], [], [], []
            for c in LIVE:
                f = fits[(camp, v, fv, c)]
                if not f["p"]:
                    continue
                p = float(f["p"])
                rows.append([camp, v, fv, c, f["n_fit"], f["p"], f["p_err"], f"{1 / (1 - p):.4f}", f"{1 - math.exp(-p):.4f}", f["chi2_BC"], f["ndf"],
                             f"{float(f['chi2_BC']) / int(f['ndf']):.3f}", f["status"], f"{nk[c] / nb[c]:.4f}", f"{R[c] / 1e3:.1f}", f"{R[c] * 16e-9:.5f}",
                             f["Gf_mV"], f["off_mV"], f["off1_mV"]])
                ps.append(p)
                chis.append(float(f["chi2_BC"]) / int(f["ndf"]))
                kept.append(nk[c] / nb[c])
                pile.append(R[c] * 16e-9)
            ps = np.array(ps)
            st[(camp, v, fv)] = {"n": len(ps), "p_median": float(np.median(ps)), "p_err": float(1.2533 * ps.std() / math.sqrt(len(ps))),
                                 "p16": float(np.percentile(ps, 16)), "p84": float(np.percentile(ps, 84)), "F_of_median": 1 / (1 - float(np.median(ps))),
                                 "P2_of_median": 1 - math.exp(-float(np.median(ps))), "chi2ndf_median": float(np.median(chis)),
                                 "kept_median": float(np.median(kept)), "pileup_median": float(np.median(pile)), "n_p_at_limit": int(((ps >= 0.94) | (ps <= 0.002)).sum())}
        # contaminación de centelleo (solo aislados)
        s = d[sel]
        gs = G[camp][s["ch"]]
        same = (s["flags"] & 3) > 0
        mix = (s["flags"] & 12) > 0
        cls = {"all_gt0.5G": s["a"] > 0.5 * gs, "1pe_0.5-1.5G": (s["a"] > 0.5 * gs) & (s["a"] < 1.5 * gs), "ge2pe_gt1.5G": s["a"] >= 1.5 * gs,
               "ge3pe_gt2.5G": s["a"] >= 2.5 * gs}
        cont[(camp, v)] = {}
        for k, m in cls.items():
            n = int(m.sum())
            fs, fm = float(same[m].mean()), float(mix[m].mean())
            cont[(camp, v)][k] = {"n": n, "frac_same": fs, "frac_mixed": fm, "excess": fs - fm, "excess_err": math.sqrt(fs * (1 - fs) / n + fm * (1 - fm) / n)}
with open(OUT + "/iso_channels.csv", "w") as f:
    wr = csv.writer(f, lineterminator="\n")
    wr.writerow(["campaign", "isolation_variant", "fit_variant", "ch", "n_fit", "p", "p_err", "F", "P_N_ge_2", "chi2_BC", "ndf", "chi2_ndf", "status", "kept_fraction",
                 "rate_gt0.5G_kHz", "residual_pileup_R_x_16ns", "Gf_mV", "off_mV", "off1_mV"])
    wr.writerows(rows)
with open(OUT + "/iso_station.csv", "w") as f:
    wr = csv.writer(f, lineterminator="\n")
    wr.writerow(["campaign", "isolation_variant", "fit_variant", "n_ch", "p_median", "p_err", "p16", "p84", "F_of_median", "P_N_ge_2_of_median", "chi2ndf_median",
                 "kept_median", "residual_pileup_median", "n_p_at_limit"])
    for (camp, v, fv), s in st.items():
        wr.writerow([camp, v, fv, s["n"], f"{s['p_median']:.4f}", f"{s['p_err']:.4f}", f"{s['p16']:.4f}", f"{s['p84']:.4f}", f"{s['F_of_median']:.4f}",
                     f"{s['P2_of_median']:.4f}", f"{s['chi2ndf_median']:.3f}", f"{s['kept_median']:.4f}", f"{s['pileup_median']:.5f}", s["n_p_at_limit"]])
with open(OUT + "/iso_contamination.csv", "w") as f:
    wr = csv.writer(f, lineterminator="\n")
    wr.writerow(["campaign", "isolation_variant", "amplitude_class", "n_isolated", "frac_neighbour_ge3pe_same_event", "frac_mixed", "excess", "err"])
    for (camp, v), cc in cont.items():
        for k, x in cc.items():
            wr.writerow([camp, v, k, x["n"], f"{x['frac_same']:.5f}", f"{x['frac_mixed']:.5f}", f"{x['excess']:.5f}", f"{x['excess_err']:.5f}"])
# p_oct/p_mayo por MPPC
bars = [r for r in csv.DictReader(open(R7)) if r["plane"] in ("X", "Y") and int(r["ch_may"]) < 32]
ratio = {}
for v in VARS:
    for fv in ("A", "D"):
        rr = []
        for b in bars:
            cm, co = int(b["ch_may"]), int(b["ch_oct"])
            if DEAD in (cm, co):
                continue
            pm, po = fits[("may", v, fv, cm)]["p"], fits[("oct", v, fv, co)]["p"]
            if pm and po:
                rr.append(float(po) / float(pm))
        rr = np.array(rr)
        ratio[(v, fv)] = {"n": len(rr), "median_ratio_by_mppc": float(np.median(rr)), "err": float(1.2533 * rr.std() / math.sqrt(len(rr))),
                          "ratio_of_medians": st[("oct", v, fv)]["p_median"] / st[("may", v, fv)]["p_median"]}
# ---------------------------------------------------------------- 2. juntas de epoxi
rc = [r for r in csv.DictReader(open(RC)) if r["variant"] == "thr3_far"]
R_ = np.array([float(r["ratio_oct_over_may"]) for r in rc])
E_ = np.array([float(r["err"]) for r in rc])


def nll_tau(mu, tau2):
    v = E_ ** 2 + tau2
    return 0.5 * float((np.log(v) + (R_ - mu) ** 2 / v).sum())


def prof_tau(tau2):
    v = E_ ** 2 + tau2
    mu = float((R_ / v).sum() / (1 / v).sum())
    return nll_tau(mu, tau2), mu


tgrid = np.linspace(0, 0.03, 30001)
prof = np.array([prof_tau(t)[0] for t in tgrid])
ib = int(np.argmin(prof))
t_best = float(tgrid[ib])
up = float(tgrid[ib:][np.argmax(prof[ib:] - prof[ib] >= 1.355)])  # −2Δln L = 2.71
mu_best = prof_tau(t_best)[1]
W = 1 / E_ ** 2
mw = float((W * R_).sum() / W.sum())
chi = float((W * (R_ - mw) ** 2).sum())


def nll_k(k, mu):
    # cada canal: −ln L normal con la razón μ o 0.74 μ; se eligen los k canales con mayor ganancia
    a = 0.5 * ((R_ - mu) / E_) ** 2
    b = 0.5 * ((R_ - 0.74 * mu) / E_) ** 2
    gain = np.sort(a - b)[::-1]
    return float(a.sum() - gain[:k].sum())


mus = np.linspace(0.4, 1.0, 6001)
kres = []
for k in range(len(R_) + 1):
    v = [nll_k(k, m) for m in mus]
    i = int(np.argmin(v))
    kres.append((k, v[i], float(mus[i])))
kbest = min(kres, key=lambda x: x[1])
kup = max(k for k, v, _ in kres if 2 * (v - kbest[1]) < 2.71 and k >= kbest[0])
epoxy = {"n_channels": len(R_), "weighted_mean": mw, "chi2_const": chi, "ndf": len(R_) - 1, "tau2_mle": t_best, "tau_mle": math.sqrt(t_best),
         "tau2_95CL_upper": up, "tau_95CL_upper": math.sqrt(up), "mu_at_tau_mle": mu_best, "tau_95CL_upper_relative_to_mu": math.sqrt(up) / mu_best,
         "k_best": kbest[0], "mu_at_k_best": kbest[2], "k_95CL_upper": kup,
         "nll_profile_k": [[k, 2 * (v - kbest[1]), m] for k, v, m in kres[:12]], "loss_model": 0.74}
json.dump({"station": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in st.items()}, "p_ratio": {f"{k[0]}|{k[1]}": v for k, v in ratio.items()},
           "contamination": {f"{k[0]}|{k[1]}": v for k, v in cont.items()}, "epoxy": epoxy}, open(OUT + "/iso_xtalk.json", "w"), indent=1)
for (camp, v, fv), s in st.items():
    print(f"{camp} {v} {fv}: p {s['p_median']:.4f} ± {s['p_err']:.4f} (p16 {s['p16']:.4f}, p84 {s['p84']:.4f}), F {s['F_of_median']:.4f}, P(N≥2) {s['P2_of_median']:.4f}, "
          f"χ²/ndf {s['chi2ndf_median']:.3f}, conservados {s['kept_median']:.4f}, apilamiento {s['pileup_median']:.5f}, en el límite {s['n_p_at_limit']}", flush=True)
for (v, fv), r in ratio.items():
    print(f"p_oct/p_mayo {v} {fv}: por MPPC {r['median_ratio_by_mppc']:.4f} ± {r['err']:.4f}, cociente de medianas {r['ratio_of_medians']:.4f}", flush=True)
for (camp, v), cc in cont.items():
    print(camp, v, {k: (x["n"], round(x["frac_same"], 5), round(x["frac_mixed"], 5), round(x["excess"], 5), round(x["excess_err"], 5)) for k, x in cc.items()}, flush=True)
print("epoxi:", json.dumps({k: v for k, v in epoxy.items() if k != "nll_profile_k"}), flush=True)
print("epoxi perfil k:", epoxy["nll_profile_k"], flush=True)
