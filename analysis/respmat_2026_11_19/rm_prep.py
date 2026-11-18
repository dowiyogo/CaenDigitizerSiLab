#!/usr/bin/python3 -I
# rm_prep.py — seguimiento 2026-11-19 (matriz de respuesta), preparación de la inyección para una campaña (mayo v8 en tiempo bueno, octubre run 8).
#  1. Comprobación del estimador: el aS de rm_extract.C (<pref>_ev.bin) contra el aS de hd64_*_hits.bin (hitdef_extract.C, entrada de los espectros), evento por
#     evento y canal por canal (igualdad exacta en float).
#  2. Plantilla de 1 p.e. por canal: pulsos oscuros aislados (±150 ns = |Δk| ≤ 9, región pre-trigger k ∈ [15,33]; mayo solo en tiempo bueno) con CFD válido
#     (0 < k − t_CFD ≤ 2.5) y sin actividad previa (|y| < 0.5 G en k−8 … k−3), con amplitud del estimador de pico en la ventana de 1 p.e.
#     [0.45 G, valle entre 1 y 2 p.e.]. Promedio de y alineado en t_CFD en bins de 1/8 de muestra (τ = muestra − t_CFD ∈ [−7.5, 10)); línea base local =
#     promedio en τ ∈ [−3, −1.5); plantilla = promedio − línea base local para τ ≥ −1.5 y 0 antes; normalizada a pico 1 (máximo del promedio suavizado
#     en 3 bins); cola para τ > 8 = exponencial ajustada en τ ∈ [4, 8] [supuesto], hasta τ = 40.
#     Normalización a 1 p.e.: amplitud por mínimos cuadrados de cada pulso aislado (todas las amplitudes) con la forma fija en su t_CFD y una línea base
#     libre (muestras con τ ∈ [−3, 6]); el pico de 1 p.e. de ese espectro (gaussiana en ±0.2 G alrededor del máximo en [0.6 G, 1.5 G]) = A1 (mV);
#     plantilla en mV = A1 · forma. También el pico de 2 p.e. (separación A2 − A1).
#  3. Tiempos de los muones etiquetados por capa: selección de spectra_export.py (etiqueta = tira ortogonal del intervalo con aS > 10 G, cualquiera de los
#     cuatro intervalos 21-54 mm; sonda > 3 p.e. medidos), t_CFD del pulso de aS; mismo − mezclados (evento siguiente) × N/(N − 1); histograma neto de
#     floor(t_CFD) = distribución de la muestra entera para la inyección (negativos a 0).
#  4. Máscara de segmentos de ruido (orden del archivo): canal c del evento e elegible si (mayo: tiempo bueno) y ninguna tira viva del plano ortogonal con
#     aS > 10 G y ninguna tira del plano de c con índice a ±2 de c (incluida c) con aS > 10 G [supuesto: "sin muón etiquetado" en sentido estricto].
# Salida (<out>/work): <camp>_tmpl.bin (float64 [32][NT], mV), <camp>_tdist.bin (float64 [32][144], probabilidad de la muestra entera),
#                      <camp>_mask.bin (uint8 [nev][32]); <out>/template_<camp>.csv, muon_time_<camp>.csv, prep_<camp>.json.
# Uso: rm_prep.py <out dir> may|oct
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402
sys.path.insert(0, os.path.join(os.environ.get("ROOTSYS", os.path.expanduser("~/root6")), "lib"))  # python3 -I ignora PYTHONPATH
import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)
OUT, CAMP = sys.argv[1], sys.argv[2]
W = OUT + "/work"
PRE = {"may": "may_v8", "oct": "oct_run8"}[CAMP]
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
HITS = PK + {"may": "/hd64_may_v8_hits.bin", "oct": "/hd64_oct8_hits.bin"}[CAMP]
ORDER = json.load(open(PK + "/stage3.json"))["strip_orders"][CAMP]["order"]
SIGN = json.load(open(REPO + "/outputs/hod2018_alongbar_amp/2026-10-31_profiles/layer_from_channels.json"))
SKEY = {"may": "may2026_v8_goodtime", "oct": "oct2026_run8"}[CAMP]
G = np.array([float(v) for v in open(f"{W}/G_{CAMP}.txt").read().split(",")])
DEAD, DEAD_X = 3, 12
D = 37.5 - 15.0 + 2.0 * np.arange(16)
BINS = [(21, 29), (29, 37), (37, 45), (45, 54)]
EVD = np.dtype([("t", "<f8"), ("aS", "<f4", 32), ("tS", "<f4", 32)])
ISD = np.dtype([("ev", "<u4"), ("ch", "u1"), ("isod", "u1"), ("k", "<u2"), ("a", "<f4"), ("t", "<f4"), ("y", "<f4", 18)])
TAU0, DT, NT = -1.5, 0.125, 333  # rejilla de la plantilla: τ = −1.5 … +40
FB0, NFB = -7.5, 140  # rejilla del promedio: τ = −7.5 … +10
res = {"campaign": CAMP}

# 1. comprobación del estimador
ev = np.fromfile(f"{W}/{PRE}_ev.bin", dtype=EVD)
hits = np.fromfile(HITS, dtype=xy.hits_dtype(64))
n = len(ev)
res["estimator_check"] = {"events_rm_extract": n, "events_hits_bin": int(len(hits)),
                          "t_equal": bool(len(hits) == n and np.array_equal(ev["t"], hits["t"])),
                          "aS_mismatch_entries": int((ev["aS"] != hits["aS"][:, :32]).sum()) if len(hits) == n else None,
                          "aS_max_abs_diff_mV": float(np.abs(ev["aS"] - hits["aS"][:, :32]).max()) if len(hits) == n else None}
print("comprobación del estimador:", res["estimator_check"], flush=True)

# orden temporal y tiempo bueno (como spectra_export.py)
idx = np.argsort(ev["t"], kind="stable")
good_sorted = mc.good_time(ev["t"][idx], "v8") if CAMP == "may" else np.ones(n, bool)
good_file = np.zeros(n, bool)
good_file[idx] = good_sorted
es = ev[idx][good_sorted]
N = len(es)
s = N / (N - 1)
res["events_selected"] = int(N)
lay = {}
for probe in ("X", "Y"):
    for i, c in enumerate(ORDER[f"h18_{probe}"]):
        lay[int(c)] = (probe, "A" if i % 2 == 0 else "B", i)

# 3. tiempos de los muones etiquetados por capa
TB = np.arange(40.0, 64.0001, 0.125)
tdist = np.zeros((32, 144))
trows = []
res["muon_time"] = {}
for probe, refc in (("X", "Y"), ("Y", "X")):
    refch = np.asarray(ORDER[f"h18_{refc}"])
    tagS = es["aS"][:, refch] > 10.0 * G[refch]
    tagS[:, refch == DEAD] = False
    for L in ("A", "B"):
        chans = [c for i, c in enumerate(ORDER[f"h18_{probe}"]) if i % 2 == (0 if L == "A" else 1) and c != DEAD]
        T = np.zeros(N, bool)
        for lo, hi in BINS:
            dist = lambda j: D[j] if SIGN[f"{SKEY}/{probe}{L}"]["sign"] < 0 else 75.0 - D[j]  # noqa: E731
            js = [j for j in range(16) if lo <= dist(j) < hi + (1e-9 if hi == 54 else 0)]
            js = [j for j in js if not (probe == "Y" and j == DEAD_X)]
            T |= tagS[:, js].any(1)
        Tn = np.concatenate([T[1:], [False]])
        hs, hm = np.zeros(len(TB) - 1), np.zeros(len(TB) - 1)
        isn, imn = np.zeros(144), np.zeros(144)
        for c in chans:
            sel = es["aS"][:, c] / G[c] > 3.0
            for msk, h, hi_ in ((sel & T, hs, isn), (sel & Tn, hm, imn)):
                v = es["tS"][msk, c]
                h += np.histogram(v, bins=TB)[0]
                hi_ += np.bincount(np.floor(v).astype(int), minlength=144)[:144]
        net = hs - s * hm
        inet = np.clip(isn - s * imn, 0, None)
        p = inet / inet.sum()
        for c in chans:
            tdist[c] = p
        cen = 0.5 * (TB[1:] + TB[:-1])
        mean = float((net * cen).sum() / net.sum())
        rms = float(np.sqrt((net * (cen - mean) ** 2).sum() / net.sum()))
        res["muon_time"][f"{probe}{L}"] = {"n_same": float(hs.sum()), "n_mix_scaled": float(s * hm.sum()), "n_net": float(net.sum()), "mean_tCFD_samples": mean,
                                         "rms_tCFD_samples": rms, "p_integer_sample": {str(j): float(p[j]) for j in range(144) if p[j] > 0}}
        for k in range(len(TB) - 1):
            trows.append([CAMP, f"{probe}{L}", f"{TB[k]:.3f}", f"{TB[k + 1]:.3f}", int(hs[k]), f"{s * hm[k]:.4f}", f"{net[k]:.4f}"])
        print(f"tiempos {probe}{L}: neto {net.sum():.1f}, media {mean:.3f}, RMS {rms:.3f} muestras; enteras",
              {j: round(float(p[j]), 4) for j in range(144) if p[j] > 0.001}, flush=True)
with open(f"{OUT}/muon_time_{CAMP}.csv", "w") as f:
    f.write(f"# hod2018 {CAMP}: CFD time (50 %, samples of the 144-sample span, signal window [49,58)) of the aS pulse of tagged-muon probes (tag > 10 p.e. in any "
            f"position bin 21-54 mm, probe > 3 p.e.); net = same - mixed x {s:.8f}\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "layer", "t_lo", "t_hi", "n_same", "n_mix_scaled", "n_net"])
    w.writerows(trows)
tdist.astype("<f8").tofile(f"{W}/{CAMP}_tdist.bin")

# 4. máscara de ruido (orden temporal → orden del archivo)
big = ev["aS"][idx] > 10.0 * G[None, :]
big[:, DEAD] = False
m_sorted = np.zeros((n, 32), bool)
for c in range(32):
    if c == DEAD:
        continue
    probe, L, i = lay[c]
    orth = [int(x) for x in ORDER[f"h18_{'Y' if probe == 'X' else 'X'}"] if x != DEAD]
    same = [int(ORDER[f"h18_{probe}"][j]) for j in range(max(0, i - 2), min(15, i + 2) + 1) if ORDER[f"h18_{probe}"][j] != DEAD]
    m_sorted[:, c] = good_sorted & ~big[:, orth].any(1) & ~big[:, same].any(1)
mask = np.zeros((n, 32), np.uint8)
mask[idx] = m_sorted
mask.tofile(f"{W}/{CAMP}_mask.bin")
res["noise_segments_per_channel"] = {str(c): int(mask[:, c].sum()) for c in range(32)}
print("segmentos de ruido por canal:", res["noise_segments_per_channel"], flush=True)


# 2. plantillas
def gaus_fit(x, lo, hi, nb):
    h = ROOT.TH1D("h", "", nb, lo, hi)
    for v in x[(x >= lo) & (x < hi)]:
        h.Fill(float(v))
    if h.GetEntries() < 50:
        h.Delete()
        return None
    mx = h.GetBinCenter(h.GetMaximumBin())
    r = h.Fit("gaus", "QSN0", "", mx - 0.2, mx + 0.2)
    out = None
    if int(r) == 0:
        out = (r.Parameter(1), r.ParError(1), abs(r.Parameter(2)), r.ParError(2))
    h.Delete()
    return out


iso = np.fromfile(f"{W}/{PRE}_iso.bin", dtype=ISD)
iso = iso[good_file[iso["ev"]]]
tmpl = np.zeros((32, NT))
trow, tinfo = [], {}
for c in range(32):
    if c == DEAD:
        continue
    q = iso[iso["ch"] == c]
    dk = q["k"].astype(float) - q["t"]
    ok = (dk > 0) & (dk <= 2.5) & (np.abs(q["y"][:, :6]).max(1) < 0.5 * G[c])
    q, dk = q[ok], dk[ok]
    u = q["a"] / G[c]
    hh, ee = np.histogram(u, bins=np.arange(1.0, 1.6001, 0.05))
    hs3 = np.convolve(hh, np.ones(3) / 3, mode="same")
    valley = float(0.5 * (ee[1:] + ee[:-1])[1:-1][np.argmin(hs3[1:-1])])
    one = (u >= 0.45) & (u < valley)
    tau = (np.arange(18) - 8)[None, :] + dk[one][:, None]
    fb = np.floor((tau - FB0) / DT).astype(int)
    okb = (fb >= 0) & (fb < NFB)
    sm = np.bincount(fb[okb], weights=q["y"][one][okb], minlength=NFB)
    cn = np.bincount(fb[okb], minlength=NFB)
    avg = sm / np.maximum(cn, 1)
    fc = FB0 + DT * (np.arange(NFB) + 0.5)
    base = avg[(fc >= -3.0) & (fc < -1.5)].mean()
    sh = np.where(fc >= -1.5, avg - base, 0.0)
    pk = np.convolve(sh, np.ones(3) / 3, mode="same")
    shape = sh / pk.max()
    # cola exponencial para τ > 8
    tl = (fc >= 4) & (fc < 8)
    sl, it = np.polyfit(fc[tl], np.log(np.clip(shape[tl], 1e-6, None)), 1)
    grid = TAU0 + DT * np.arange(NT)
    T = np.interp(grid, fc, shape)
    T[grid > 8.0] = np.exp(it + sl * grid[grid > 8.0])
    # amplitud por mínimos cuadrados (forma fija en t_CFD, línea base libre) de todos los pulsos aislados limpios
    tau_all = (np.arange(18) - 8)[None, :] + dk[:, None]
    fm = np.interp(tau_all, grid, T, left=0.0, right=0.0)
    wsel = (tau_all >= -3) & (tau_all <= 6)
    S1 = wsel.sum(1)
    Sf = (fm * wsel).sum(1)
    Sff = (fm * fm * wsel).sum(1)
    Sy = (q["y"] * wsel).sum(1)
    Sfy = (fm * q["y"] * wsel).sum(1)
    den = S1 * Sff - Sf * Sf
    A = (S1 * Sfy - Sf * Sy) / den
    b = (Sff * Sy - Sf * Sfy) / den
    a1 = gaus_fit(A / G[c], 0.6, 1.5, 45)
    a2 = gaus_fit(A / G[c], a1[0] + 0.6, a1[0] + 1.5, 45) if a1 else None
    # mismo pico con el estimador de pico (amplitud aS-like desde bl) para la comprobación del punto 5
    p1 = gaus_fit(u, 0.45, valley, 40)
    A1 = a1[0] * G[c]
    tmpl[c] = A1 * T
    tinfo[c] = {"n_isolated_clean": int(len(q)), "n_isolated_raw": int((iso["ch"] == c).sum()), "valley_peak_estimator_G": valley, "n_template": int(one.sum()),
                "local_baseline_mV": float(base), "peak_tau_samples": float(fc[np.argmax(pk)]), "tail_decay_samples": float(-1.0 / sl),
                "shape_at_tau_8": float(np.interp(8.0, grid, T)), "A1_LS_over_G": a1[0], "A1_LS_over_G_err": a1[1], "sigma1_LS_over_G": a1[2],
                "A2_LS_over_G": a2[0] if a2 else None, "A2_LS_over_G_err": a2[1] if a2 else None, "A1_mV": A1,
                "LS_local_baseline_median_mV": float(np.median(b)),
                "peak_est_1pe_pos_over_G": p1[0] if p1 else None, "peak_est_1pe_pos_err": p1[1] if p1 else None,
                "peak_est_1pe_sigma_over_G": p1[2] if p1 else None, "peak_est_1pe_sigma_err": p1[3] if p1 else None,
                "frac_kt_hist_8": np.histogram(dk[one] % 1, bins=8, range=(0, 1))[0].tolist()}
    for j in range(NT):
        trow.append([CAMP, c, f"{grid[j]:.3f}", f"{T[j]:.6f}", f"{tmpl[c, j]:.6f}"])
    print(f"ch{c}: limpios {len(q)}, plantilla {one.sum()}, valle {valley:.3f} G, base local {base:.3f} mV, A1/G {a1[0]:.4f} ± {a1[1]:.4f} (σ {a1[2]:.4f}), "
          f"A2/G {a2[0] if a2 else float('nan'):.4f}, cola τ {-1 / sl:.2f}, pico del estimador {p1[0] if p1 else float('nan'):.4f} (σ {p1[2] if p1 else float('nan'):.4f})",
          flush=True)
tmpl.astype("<f8").tofile(f"{W}/{CAMP}_tmpl.bin")
with open(f"{OUT}/template_{CAMP}.csv", "w") as f:
    f.write(f"# hod2018 {CAMP}: single-p.e. template per channel, tau = sample - CFD time (samples of 16 ns); shape normalized to peak 1; template_mV = A1 x shape "
            "(A1 = 1 p.e. peak of the least-squares amplitude of isolated dark pulses); tail for tau > 8 = exponential fitted in [4,8]\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "ch", "tau_samples", "shape", "template_mV"])
    w.writerows(trow)
res["templates"] = {str(c): v for c, v in tinfo.items()}
res["template_grid"] = {"tau0": TAU0, "dtau": DT, "n": NT}
json.dump(res, open(f"{OUT}/prep_{CAMP}.json", "w"), indent=1)
