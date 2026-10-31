# h19_sync.py — seguimiento 2026-10-31 (5), punto 1 (solo lectura): pulsos sincronizados con el trigger por canal, octubre (run 8) contra
# mayo (v8, tiempo bueno nuevo con HV-DIAG-5 excluido). Por canal y umbral T (mV): fracción de eventos con amplitud > T en la ventana signal
# (aS, pulso del buscador con su mínimo en la ventana) y en la ventana fuera de tiempo (aC, misma anchura), y el exceso aS − aC con su error.
# Entrada: work/hd64_*_hits.bin de hitdef_extract.C (PROC-1). Salida: CSV a stdout.
# Uso: python3 h19_sync.py <oct_hits.bin> <may_hits.bin> <may_excluidos_utc.txt>
import sys

import numpy as np

sys.path.insert(0, "daq/tools/characterization")
from xy_estimator import load_hits  # noqa: E402

THR = (1.0, 3.0, 35.0)


def rows(tag, e):
    n = len(e)
    for c in range(64):
        for T in THR:
            s = int((e["aS"][:, c] > T).sum())
            k = int((e["aC"][:, c] > T).sum())
            ex = (s - k) / n
            err = np.sqrt(s + k) / n
            yield f"{tag},{c},{n},{T:g},{s},{k},{ex:.6f},{err:.6f}"


def main():
    oct_e = load_hits(sys.argv[1], 64)
    may_e = load_hits(sys.argv[2], 64)
    ex = np.loadtxt(sys.argv[3], ndmin=2)
    g = np.ones(len(may_e), bool)
    for a, b in ex:
        g &= ~((may_e["t"] >= a) & (may_e["t"] < b))
    may_e = may_e[g]
    print("sample,ch,n_events,thr_mV,n_signal,n_offtime,excess_frac,excess_err")
    for r in rows("oct_run8", oct_e):
        print(r)
    for r in rows("may_v8good", may_e):
        print(r)


if __name__ == "__main__":
    main()
