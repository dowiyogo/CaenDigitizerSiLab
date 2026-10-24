# xy_estimator.py — DAQ2-16: estimador único de la fracción X·Y de hod2018 por trigger y lector de impactos, compartido
# por charact_analyse.py, hitdef_analyse.py y biasscan.py (no se copia la fórmula en otro lugar).
# Impactos: work/hd_<run>_hits.bin de hitdef_extract.C (definición única de "Comparación con mayo": pulso del buscador
# del DAQ con su mínimo en la ventana signal, pedestal = media recortada de las 40 muestras previas; amplitud > T·G).
# Estimador (ruido independiente por plano; aX, aY = probabilidad de impacto en la ventana fuera de tiempo; rX, rY =
# probabilidad verdadera por plano):
#   mismo evento − mezclados      = (q − rX·rY)(1 − aX)(1 − aY)  ->  q = (N_XY − acc_mix) / (N(1−aX)(1−aY)) + rX·rY
#   mismo evento − fuera de tiempo = q (1 − aX)(1 − aY)            ->  q = (N_XY − acc_ctl) / (N(1−aX)(1−aY))
# (la primera es la fórmula del análisis de mayo, analysis/macros/hod2018_good_time_step2.C). ε(X|Y) = q / rY, ε(Y|X) = q / rX.
import math

import numpy as np


def hits_dtype(nch=32):
    """Registro de hitdef_extract.C con nch canales (32: hod2018; 64: PROC-1, también hod2019)."""
    return np.dtype([("t", "<f8"), ("aS", "<f4", nch), ("aC", "<f4", nch), ("aG", "<f4", nch), ("kG", "u1", nch), ("aM", "<f4", nch),
                     ("aMc", "<f4", nch)])


HITS_DTYPE = hits_dtype(32)


def load_hits(path, nch=32):
    """Registros de hitdef_extract.C ordenados por tiempo."""
    e = np.fromfile(path, dtype=hits_dtype(nch))
    return e[np.argsort(e["t"], kind="stable")]


def plane_hits(e, G, T, chans):
    """(en tiempo, fuera de tiempo): algún canal de chans con amplitud > T·G (definición única)."""
    chans = list(chans)
    thr = T * np.asarray(G, float)[chans]
    return (e["aS"][:, chans] > thr).any(1), (e["aC"][:, chans] > thr).any(1)


def xy_fraction(Xs, Ys, Xc, Yc):
    """q (fracción de triggers con muón en X y en Y) con los dos estimadores, sus errores estadísticos y ε(X|Y), ε(Y|X).
    Xs, Ys: impacto en la ventana signal; Xc, Yc: en la ventana fuera de tiempo (arreglos booleanos por evento, en orden
    de tiempo: los mezclados combinan el evento i con el i+1)."""
    Xs, Ys, Xc, Yc = (np.asarray(a, bool) for a in (Xs, Ys, Xc, Yc))
    N = len(Xs)
    nxy = int((Xs & Ys).sum())
    a1, a2, a3 = int((Xs & Yc).sum()), int((Xc & Ys).sum()), int((Xc & Yc).sum())
    m1, m2 = int((Xs[:-1] & Ys[1:]).sum()), int((Ys[:-1] & Xs[1:]).sum())
    out = {"N": N, "N_XY": nxy, "acc_ctl": a1 + a2 - a3, "acc_mix": 0.5 * (m1 + m2) * N / (N - 1) if N > 1 else float("nan"),
           "N_X": int(Xs.sum()), "N_Y": int(Ys.sum()), "N_Xc": int(Xc.sum()), "N_Yc": int(Yc.sum())}
    if N < 2:
        return out
    aX, aY = out["N_Xc"] / N, out["N_Yc"] / N
    rX = (out["N_X"] / N - aX) / (1 - aX) if aX < 1 else float("nan")
    rY = (out["N_Y"] / N - aY) / (1 - aY) if aY < 1 else float("nan")
    den = (1 - aX) * (1 - aY)
    out.update(aX=aX, aY=aY, rX=rX, rY=rY)
    for tag, acc, var, add in (("mix", out["acc_mix"], (m1 + m2) / 4, rX * rY), ("ctl", out["acc_ctl"], a1 + a2 + a3, 0.0)):
        q = (nxy - acc) / (N * den) + add if den > 0 else float("nan")
        eq = math.sqrt(nxy + var) / (N * den) if den > 0 else float("nan")
        out[f"f_{tag}"], out[f"f_{tag}_err"] = q, eq
        for a, r in (("X|Y", rY), ("Y|X", rX)):
            e = q / r if r > 0 else float("nan")
            out[f"e_{a}_{tag}"] = e
            out[f"e_{a}_{tag}_err"] = math.sqrt((eq / r) ** 2 + max(e * (1 - e), 0) / (r * N)) if r > 0 else float("nan")
    return out
