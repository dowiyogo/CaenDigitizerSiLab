#!/usr/bin/python3 -I
# mip_spectra.py — pedido del proyecto de simulación (2026-10-30): espectro de amplitud de hod2018 por canal, en la ventana signal y en
# la ventana fuera de tiempo del mismo ancho, en unidades de la ganancia de 1 p.e. de cada canal, bins de 1 p.e.
# Fuente: los archivos de impactos de PROC (outputs/proceeding_2026_10/work/hd64_<dataset>_hits.bin, hitdef_extract.C, nch 64):
#   amplitud = la del pulso más alto del buscador del DAQ (daq/core/include/daq/pulse_finder.h) cuyo mínimo cae en la ventana:
#   amplitud de pico (pedestal − mínimo, mV), no carga integrada. Pedestal = media recortada ±2 mV de las 40 muestras previas al
#   tramo; buscador con umbral 0.5 mV y separación 4 muestras. Ventana signal [49,58) del tramo de 144 muestras (9 × 16 ns = 144 ns;
#   mayo: muestras 258-266 del registro = 4128-4256 ns; octubre: 49-57); fuera de tiempo [99,108), mismo ancho, 800 ns después.
#   Sin pulso en la ventana: amplitud 0 (bin 0). Sin corrección de saturación ni de crosstalk.
# G por canal: octubre = c1 (commissioning2_2026_10_21 run 2, 2026-10-21 10:32-10:41 UTC); mayo = G reajustada de PROC-1 etapa 1 sobre
#   v8 (work/hd64_may_v8_peaks.csv, commit 8027437). ch3 sin G (muerto) queda fuera. Canales sin ajuste aceptado: se listan y se usa
#   la mediana (como PROC).
# Conjuntos: octubre = hodoscope_2026_10 run 8 completo (hd64_oct8); mayo = run B (v8), tiempo bueno (may_common.good_time).
# Salida en <out>/<dataset>/: ch00.csv … ch31.csv y sum.csv (pe_lo, pe_hi, n_signal, n_offtime), summary.json.
# Uso: mip_spectra.py <out dir>
import csv, json, math, os, sys, time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
CH = [c for c in range(32) if c != 3]
NB = 300  # bins de 1 p.e. de 0 a 300 (más desborde)


def g_may():
    g, st = np.full(32, np.nan), {}
    for r in csv.DictReader(open(W + "/hd64_may_v8_peaks.csv")):
        c = int(r["ch"])
        st[c] = r["status"]
        if r["status"] in ("ok", "changed"):
            g[c] = float(r["G_fit"])
    return g, st


def g_oct():
    c1 = json.load(open(os.path.join(REPO, "daq/config/calibration/pe_calibration_2026_10_c1.json")))
    g, st = np.full(32, np.nan), {}
    for e in c1["channels"]:
        if e["ch"] < 32:
            st[e["ch"]] = e["status"]
            if isinstance(e.get("G_mV"), (int, float)):
                g[e["ch"]] = e["G_mV"]
    return g, st


summary = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "bins": "1 p.e., [k, k+1) para k = 0…299, más desborde ≥ 300",
           "datasets": {}}
for ds, fn, gfun, gsrc, good in (("oct2026_run8", "hd64_oct8_hits.bin", g_oct, "c1 (commissioning2_2026_10_21 run 2, 2026-10-21)", None),
                                  ("may2026_v8_goodtime", "hd64_may_v8_hits.bin", g_may, "reajuste PROC-1 etapa 1 sobre v8 (hd64_may_v8_peaks.csv, 2026-10-24)", "v8")):
    e = xy.load_hits(os.path.join(W, fn), 64)
    if good:
        e = e[mc.good_time(e["t"], good)]
    G, st = gfun()
    filled = [c for c in CH if not np.isfinite(G[c])]
    med = float(np.nanmedian(G[CH]))
    G = np.where(np.isfinite(G), G, med)
    d = os.path.join(OUT, ds)
    os.makedirs(d, exist_ok=True)
    edges = np.concatenate([np.arange(0, NB + 1, 1.0), [np.inf]])
    tot_s = np.zeros(NB + 1); tot_c = np.zeros(NB + 1)
    chans = {}
    for c in CH:
        s = e["aS"][:, c] / G[c]
        o = e["aC"][:, c] / G[c]
        hs = np.histogram(s, bins=edges)[0]; ho = np.histogram(o, bins=edges)[0]
        tot_s += hs; tot_c += ho
        with open(os.path.join(d, f"ch{c:02d}.csv"), "w") as f:
            f.write(f"# {ds} ch{c:02d}; G = {G[c]:.4f} mV/p.e. ({'mediana, sin ajuste aceptado' if c in filled else gsrc}); N triggers = {len(e)}\n")
            f.write("pe_lo,pe_hi,n_signal,n_offtime\n")
            for k in range(NB + 1):
                f.write(f"{k},{'inf' if k == NB else k + 1},{hs[k]},{ho[k]}\n")
        # resumen: (a) impactos en tiempo > 5 p.e.; (b) espectro sustraído (signal − fuera de tiempo) por encima de 5 p.e.
        hit = s[s > 5]
        sub = (hs - ho).astype(float)[5:]
        cs = np.cumsum(np.clip(sub, 0, None))
        def q(p):
            if cs[-1] <= 0:
                return None
            return float(5 + np.searchsorted(cs, p * cs[-1]) + 0.5)
        chans[c] = {"G_mV": float(G[c]), "G_status": st.get(c), "G_filled_with_median": c in filled,
                    "frac_signal_gt5": float((s > 5).mean()), "frac_offtime_gt5": float((o > 5).mean()),
                    "signal_gt5_median_pe": float(np.median(hit)) if len(hit) else None, "signal_gt5_p90_pe": float(np.percentile(hit, 90)) if len(hit) else None,
                    "subtracted_ge5_median_pe": q(0.5), "subtracted_ge5_p90_pe": q(0.9), "n_signal_ge300": int(hs[NB])}
    with open(os.path.join(d, "sum.csv"), "w") as f:
        f.write(f"# {ds} suma de {len(CH)} canales (sin ch3); N triggers = {len(e)}\n")
        f.write("pe_lo,pe_hi,n_signal,n_offtime\n")
        for k in range(NB + 1):
            f.write(f"{k},{'inf' if k == NB else k + 1},{int(tot_s[k])},{int(tot_c[k])}\n")
    summary["datasets"][ds] = {"source": f"outputs/proceeding_2026_10/work/{fn}", "good_time": good, "n_triggers": int(len(e)), "gain_source": gsrc,
                               "channels_without_gain_used_median": filled, "median_gain_mV": med, "channels": chans}
    print(ds, "N", len(e), "sin G (mediana):", filled)
    for c in CH:
        v = chans[c]
        print(f"  ch{c:02d} G {v['G_mV']:.3f}  >5 p.e. en tiempo {v['frac_signal_gt5']:.4f} fuera {v['frac_offtime_gt5']:.4f}  "
              f"mediana/p90 (en tiempo >5) {v['signal_gt5_median_pe']:.1f}/{v['signal_gt5_p90_pe']:.1f}  "
              f"(sustraído ≥5) {v['subtracted_ge5_median_pe']}/{v['subtracted_ge5_p90_pe']}")
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, ensure_ascii=False)
