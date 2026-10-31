#!/usr/bin/python3 -I
# period_gain.py — seguimiento 2026-10-31 (9), punto 1: G de respaldo por el período del peine de picos de p.e. del espectro de pulsos oscuros,
# para tensiones bajas donde el ajuste de c1 (gain_fit.h) no converge (a −0.5 V: 5 de 31 canales con 12 000 eventos y no mejora con más).
# Espectro: bins de 1 LSB (salida de vbd_gain.C / gain_spectra.C). Ventana [lo, hi]: lo = 1.5 mV, hi = el último bin con ≥ 20 cuentas (≤ 26 mV).
# y = ln(n + 1) menos una tendencia cuadrática ajustada en la ventana; periodograma P(T) = |Σ r_i e^{2πi x_i/T}|² para T entre 0.7 y 1.4 ×
# G esperada (paso 0.005 mV; con 0.6-1.6 a veces elige 1.5 G); G = argmax; un máximo en el borde del intervalo = sin período. Error: desviación de G en réplicas de Poisson del espectro (por omisión 20).
# Uso como módulo (vbd_analysis.py) o: period_gain.py <spec.csv> <G esperada por canal: json {ch: mV}> [réplicas]
import json, math, sys

import numpy as np


def period_gain(h, bin_mv, g_exp, rng=None, replicas=20, lo_mv=1.5, hi_max_mv=26.0, tmin=0.7, tmax=1.4):
    h = np.asarray(h, float)
    x = (np.arange(len(h)) + 0.5) * bin_mv

    def one(hh):
        ok = np.where(hh >= 20)[0]
        if len(ok) == 0:
            return None
        hi = min(hi_max_mv, x[ok.max()])
        m = (x >= lo_mv) & (x <= hi)
        if m.sum() < 8 or hi - lo_mv < 2.5 * g_exp:
            return None
        xx, y = x[m], np.log(hh[m] + 1.0)
        r = y - np.polyval(np.polyfit(xx, y, 2), xx)
        T = np.arange(tmin * g_exp, tmax * g_exp, 0.005)
        P = np.abs(np.exp(2j * np.pi * xx[None, :] / T[:, None]) @ r) ** 2
        k = int(np.argmax(P))
        if k == 0 or k == len(T) - 1:
            return None  # en el borde: sin período claro
        return float(T[k])

    g = one(h)
    if g is None:
        return None, None, 0
    rng = rng or np.random.default_rng(9)
    reps = [one(rng.poisson(h).astype(float)) for _ in range(replicas)]
    reps = [v for v in reps if v is not None]
    sd = float(np.std(reps)) if len(reps) >= 3 else None
    return g, sd, len(reps)


def load_spec(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    bins = np.unique(d[:, 1])
    bin_mv = float(bins[1] - bins[0])
    return {int(c): d[d[:, 0] == c][:, 2] for c in np.unique(d[:, 0])}, bin_mv


if __name__ == "__main__":
    sp, b = load_spec(sys.argv[1])
    gexp = {int(k): v for k, v in json.load(open(sys.argv[2])).items()}
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    print("ch,G_period_mV,G_period_sd_mV,n_replicas")
    for c in sorted(sp):
        if c >= 32 or c not in gexp:
            continue
        g, sd, n = period_gain(sp[c], b, gexp[c], replicas=reps)
        print(f"{c},{'' if g is None else f'{g:.4f}'},{'' if sd is None else f'{sd:.4f}'},{n}")
