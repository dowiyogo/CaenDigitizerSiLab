#!/usr/bin/python3 -I
# hitdef_analyse.py — DAQ2-15: fracción X·Y de hod2018 por trigger (accidentales restados por eventos mezclados y por
# ventana fuera de tiempo) y ε(X|Y), ε(Y|X) contra umbral 3-12 p.e., con el mismo código para mayo v8, la noche del
# 2026-10-19 y el run 8, desde work/hd_<dataset>_hits.bin (hitdef_extract.C).
# Umbral: amplitud > T·G (G de la calibración de cada campaña: mayo outputs/pe_calibration_v8.txt, octubre c1); con el
# pedestal común los picos de 1 p.e. caen en múltiplos enteros de G (p0 medido ≈ 0, ver *_peaks.csv).
# Definiciones: "pulso" (nueva) = pulso del buscador con mínimo en la ventana; "DAQ2-13" = mínimo de todo el tramo de 144
# muestras dentro de [m−4, m+5), control [m+46, m+55), con m = moda de la muestra del mínimo de los impactos > 10 p.e.
# "mayo" = definición del análisis de mayo (hod2018_good_time_step2.C): mínimo de la ventana de 10 muestras con pedestal
# local [k−12, k−4], mismo umbral T·G; además el punto publicado "5 p.e." de mayo = p0 + 3.5·G (p0 de pe_calibration_v8.txt).
# Mayo: solo tiempo bueno de runs/good_intervals_2026.csv (fuera EXCLUDE y HV_OFF; hora local UTC−4).
# Uso: hitdef_analyse.py <work> <salida.json>
import calendar, csv, json, math, sys, time

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import xy_estimator as xy  # noqa: E402
import may_common as mc  # noqa: E402

W, OUTF = sys.argv[1], sys.argv[2]
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
DT = xy.HITS_DTYPE
TS = list(range(3, 13))
X, Y = list(range(16)), list(range(16, 32))  # mayo (kCfg_2026_05) y octubre (r4): mismos conjuntos de canales por plano
EXCL = {"sin_ch3": [3], "sin_2_3_9_15": [2, 3, 9, 15]}


# lectores de calibración y tiempo bueno: may_common.py (PROC-1; antes estaban aquí)
g_may, p0_may, g_c1 = mc.g_may, mc.p0_may, mc.g_c1


def good_may(t):
    return mc.good_time(t, "v8")


DS = {"may8": ("mayo v8 (tiempo bueno)", g_may()), "night": ("noche 2026-10-19", g_c1()), "run008": ("run 8 (2026-10-22/23)", g_c1())}


frac = xy.xy_fraction  # estimador compartido (xy_estimator.py)


res = {"definition": {"pedestal": "media recortada (±2 mV de la mediana) de las 40 muestras [s0−49, s0−9)",
                      "signal": "[s0, s0+9) (144 ns); octubre s0 = 49, mayo s0 = 258 (4128 ns)",
                      "off_time": "[s0+50, s0+59) (800 ns después)",
                      "hit": "pulso del buscador del DAQ (mínimo local, separación 4 muestras) con su mínimo dentro de la ventana; "
                             "amplitud = pedestal − mínimo > T·G",
                      "daq2_13": "mínimo de todo el tramo de 144 muestras dentro de [m−4, m+5), amplitud > T·G; control [m+46, m+55)"},
       "datasets": {}}
for k, (lab, G) in DS.items():
    e = np.fromfile(f"{W}/hd_{k}_hits.bin", dtype=DT)
    e = e[np.argsort(e["t"], kind="stable")]
    if k == "may8":
        n_all = len(e)
        e = e[good_may(e["t"])]
    Gf = np.where(np.isfinite(G), G, np.nanmedian(G))
    # moda de la muestra del mínimo (impactos > 10 p.e.) para la definición de DAQ2-13
    big = e["aG"] / Gf > 10
    kk = e["kG"][big].astype(int)
    hk = np.bincount(kk, minlength=144)
    m = 20 + int(hk[20:124].argmax())  # sin los bordes del tramo (colas de pulsos anteriores en la muestra 0)
    d = {"label": lab, "events": int(len(e)), "mode_sample_daq2_13": m, "first_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(e["t"][0])),
         "last_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(e["t"][-1]))}
    if k == "may8":
        d["events_before_good_time"] = n_all
    for ex, chs in EXCL.items():
        xs = [c for c in X if c not in chs]
        ys = [c for c in Y if c not in chs]
        for T in TS:
            thr = T * Gf
            hS, hC = e["aS"] > thr, e["aC"] > thr
            new = frac(hS[:, xs].any(1), hS[:, ys].any(1), hC[:, xs].any(1), hC[:, ys].any(1))
            inw = (e["kG"] >= m - 4) & (e["kG"] < m + 5) & (e["aG"] > thr)
            ctw = (e["kG"] >= m + 46) & (e["kG"] < m + 55) & (e["aG"] > thr)
            old = frac(inw[:, xs].any(1), inw[:, ys].any(1), ctw[:, xs].any(1), ctw[:, ys].any(1))
            hM, hMc = e["aM"] > thr, e["aMc"] > thr
            mayd = frac(hM[:, xs].any(1), hM[:, ys].any(1), hMc[:, xs].any(1), hMc[:, ys].any(1))
            d.setdefault(ex, {})[T] = {"pulso": new, "daq2_13": old, "mayo": mayd}
        if k == "may8":
            thr = p0_may() + 3.5 * Gf
            hM, hMc = e["aM"] > thr, e["aMc"] > thr
            d[ex]["mayo_publicado_5pe"] = frac(hM[:, xs].any(1), hM[:, ys].any(1), hMc[:, xs].any(1), hMc[:, ys].any(1))
    res["datasets"][k] = d
json.dump(res, open(OUTF, "w"), ensure_ascii=False, indent=1)
for k, d in res["datasets"].items():
    print(k, d["label"], d["events"], "moda", d["mode_sample_daq2_13"], d["first_utc"], d["last_utc"])
    for ex in EXCL:
        if "mayo_publicado_5pe" in d[ex]:
            c = d[ex]["mayo_publicado_5pe"]
            print(f"  {ex} mayo p0+3.5G: NXY {c['N_XY']} mix {c['f_mix']:.4f}±{c['f_mix_err']:.4f} ctl {c['f_ctl']:.4f}±{c['f_ctl_err']:.4f} /0.0360 = {c['f_mix'] / 0.036:.3f}")
        for T in TS:
            a, b, c = d[ex][T]["pulso"], d[ex][T]["daq2_13"], d[ex][T]["mayo"]
            print(f"  {ex} T{T:2d} pulso mix {a['f_mix']:.4f}±{a['f_mix_err']:.4f} ctl {a['f_ctl']:.4f}±{a['f_ctl_err']:.4f} "
                  f"X|Y {a['e_X|Y_mix']:.3f}±{a['e_X|Y_mix_err']:.3f} Y|X {a['e_Y|X_mix']:.3f}±{a['e_Y|X_mix_err']:.3f} NXY {a['N_XY']} accmix {a['acc_mix']:.1f} | "
                  f"DAQ2-13 mix {b['f_mix']:.4f}±{b['f_mix_err']:.4f} ctl {b['f_ctl']:.4f} X|Y {b['e_X|Y_mix']:.3f} Y|X {b['e_Y|X_mix']:.3f} | "
                  f"mayo mix {c['f_mix']:.4f}±{c['f_mix_err']:.4f} ctl {c['f_ctl']:.4f} X|Y {c['e_X|Y_mix']:.3f} Y|X {c['e_Y|X_mix']:.3f}")
