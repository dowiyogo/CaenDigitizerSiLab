#!/usr/bin/python3 -I
# xtalk_prep.py — seguimiento 2026-11-15 (diafonía de Borel): entradas de fit_borel.C.
#  - máscara de eventos: mayo = tiempo bueno de v8 (may_common.good_time sobre <prefijo>_dark3_t.bin), octubre = todos;
#  - tasa de pulsos R por canal = pulsos con a > 0.5 G y mínimo en [16,37) / (eventos · 21 · 16 ns) (dark3_rate.csv, todos los eventos del run,
#    sin aislamiento);
#  - longitud del pulso τ = FWHM de la forma media del pulso aislado de 1 p.e. (muestras k−4 … k+12, línea base = media de k−4 y k−3,
#    interpolación lineal);
#  - probabilidad de apilamiento q = R·τ: otro pulso oscuro con el mínimo a menos de ±τ/2 del primero, que el buscador no separa y que suma su
#    amplitud [S: suma completa dentro de la FWHM, nada fuera].
# Salida: <work>/{pre}_mask.bin, <work>/xtalk_prep.csv (campaign, ch, events, rate_Hz, fwhm_ns, q_pileup); imprime las cadenas para fit_borel.C.
# Uso: xtalk_prep.py <work dir>
import csv, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

WK = sys.argv[1]
rows = []
for camp, pre in (("may", "may_v8"), ("oct", "oct_run8")):
    t = np.fromfile(f"{WK}/{pre}_dark3_t.bin", dtype="<f8")
    m = mc.good_time(t, "v8") if camp == "may" else np.ones(len(t), bool)
    m.astype(np.uint8).tofile(f"{WK}/{pre}_mask.bin")
    print(f"{camp}: eventos {len(t)}, en la máscara {int(m.sum())}")
    qs = []
    for r in csv.DictReader(open(f"{WK}/{pre}_dark3_rate.csv")):
        c, ne, npl = int(r["ch"]), int(r["events"]), int(r["n_pulses_r0_gt_half_G"])
        R = npl / (ne * 21 * 16e-9)
        s = np.array([float(r[f"s{j:+d}_mV"]) for j in range(-4, 13)])
        fw = np.nan
        if int(r["n_shape"]) > 100:
            y = s - 0.5 * (s[0] + s[1])
            k = int(np.argmax(y))
            h = 0.5 * y[k]
            i = k
            while i > 0 and y[i - 1] > h:
                i -= 1
            j = k
            while j < len(y) - 1 and y[j + 1] > h:
                j += 1
            left = (i - 1) + (h - y[i - 1]) / (y[i] - y[i - 1]) if i > 0 else float(i)
            right = j + (y[j] - h) / (y[j] - y[j + 1]) if j < len(y) - 1 else float(j)
            fw = (right - left) * 16.0
        q = R * fw * 1e-9 if np.isfinite(fw) else 0.0
        qs.append(q)
        rows.append([camp, c, ne, f"{R:.1f}", f"{fw:.3f}", f"{q:.5f}"])
    print(f"Q_{camp}=" + ",".join(f"{x:.5f}" for x in qs))
with open(f"{WK}/xtalk_prep.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "ch", "events", "rate_Hz", "fwhm_ns", "q_pileup"])
    w.writerows(rows)
for camp in ("may", "oct"):
    rr = [r for r in rows if r[0] == camp and r[1] != 3]
    print(camp, "mediana R %.1f kHz, FWHM %.3f ns, q %.5f" % (np.median([float(r[3]) for r in rr]) / 1e3, np.nanmedian([float(r[4]) for r in rr]),
                                                          np.median([float(r[5]) for r in rr])))
