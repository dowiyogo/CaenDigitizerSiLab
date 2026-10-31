#!/usr/bin/python3 -I
# h19_sync_quantiles.py — seguimiento 2026-10-31 (6), punto 3: amplitud de los pulsos de hod2019 sincronizados con el trigger en mayo v8 (tiempo
# bueno nuevo), por canal: espectro de aS (ventana signal) > 3 mV menos el de aC (ventana fuera de tiempo, misma anchura) > 3 mV, bins de 0.25 mV,
# con la definición del seguimiento 1 (bins negativos a 0) y, como sensibilidad, sin recortar (acumulada neta). Mediana y percentil 90 del
# espectro neto en mV con error por remuestreo de Poisson (200 réplicas); en p.e. solo con una G aceptada de fit_gain_long.C (archivo
# opcional "ch,G_mV").
# Entradas: work/hd64_may_v8_hits.bin (PROC, 64 canales) y el archivo de intervalos UTC excluidos del seguimiento 4.
# Salida: CSV a stdout. Uso: h19_sync_quantiles.py <may_hits.bin> <excluidos.txt> [<G_hod2019.csv>]
import csv
import sys

import numpy as np

sys.path.insert(0, "daq/tools/characterization")
from xy_estimator import load_hits  # noqa: E402

THR = 3.0
EDG = np.arange(THR, 200.0 + 0.25, 0.25)
CEN = 0.5 * (EDG[1:] + EDG[:-1])
RNG = np.random.default_rng(19)


def quant(net, p, clip):
    h = np.clip(net, 0, None) if clip else net
    cs = np.cumsum(h)
    if cs[-1] <= 0:
        return np.nan
    return float(CEN[int(np.argmax(cs >= p * cs[-1]))])


def main():
    e = load_hits(sys.argv[1], 64)
    ex = np.loadtxt(sys.argv[2], ndmin=2)
    g = np.ones(len(e), bool)
    for a, b in ex:
        g &= ~((e["t"] >= a) & (e["t"] < b))
    e = e[g]
    G = {}
    if len(sys.argv) > 3:
        G = {int(r["ch"]): float(r["G_mV"]) for r in csv.DictReader(open(sys.argv[3]))}
    n = len(e)
    print("ch,n_events,n_signal_gt3mV,n_offtime_gt3mV,net,net_err,median_mV,median_err,p90_mV,p90_err,median_noclip_mV,p90_noclip_mV,G_mV,median_pe,p90_pe")
    for c in range(32, 64):
        s, o = e["aS"][:, c], e["aC"][:, c]
        hs = np.histogram(s[s > THR], bins=EDG)[0].astype(float)
        ho = np.histogram(o[o > THR], bins=EDG)[0].astype(float)
        m, p9 = quant(hs - ho, 0.5, True), quant(hs - ho, 0.9, True)
        bm, bp = [], []
        for _ in range(200):
            d = RNG.poisson(hs) - RNG.poisson(ho)
            bm.append(quant(d, 0.5, True))
            bp.append(quant(d, 0.9, True))
        gc = G.get(c)
        print(f"{c},{n},{int(hs.sum())},{int(ho.sum())},{int(hs.sum() - ho.sum())},{np.sqrt(hs.sum() + ho.sum()):.1f},{m:.3f},{np.nanstd(bm):.3f},{p9:.3f},"
              f"{np.nanstd(bp):.3f},{quant(hs - ho, 0.5, False):.3f},{quant(hs - ho, 0.9, False):.3f},"
              + (f"{gc:.4f},{m / gc:.3f},{p9 / gc:.3f}" if gc else ",,"))


if __name__ == "__main__":
    main()
