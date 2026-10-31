#!/usr/bin/python3 -I
# pct_corrected.py — seguimiento 2026-10-31 (5), punto 3: p_ct corregido por los pulsos espurios (colas partidas por el buscador) a partir de
# las salidas de pair_shapes.C (segunda pasada) y de la tabla del seguimiento 4 (dcr_ct_table.py).
#   S05 = segundos pulsos de pares con Δ = 3-6 muestras y E < Ecut·G (Ecut = 0.5 central, 0.75 variante; contados en *_E_hist.csv),
#   S15 = los mismos con A2 > 1.5 p.e. (0 medidos con Ecut = 0.5; con 0.75 se toma 0, ver el informe).
#   p_ct' = (N15 − S15) / (N05 − S05); DCR' = DCR·(1 − S05/N05); con apilamiento (τ = 32 ns, supuesto como en el seguimiento 4):
#   p_pile = 1 − exp(−DCR'·τ), p_ct'' = (p_ct' − p_pile)/(1 − p_pile).
# Para canales sin pares medidos se aplica la mediana de f = S05/N05 de los canales medidos de la misma muestra.
# Uso: pct_corrected.py <tabla_dcr_ct.csv> <prefijo pair_shapes (…_spurious.csv, …_E_hist.csv)> <etiqueta> [<salida.csv>]
import csv
import math
import statistics as st
import sys

TAU = 32e-9


def table(path):
    rows = [l for l in open(path) if not l.startswith("#")]
    return {int(r["ch"]): r for r in csv.DictReader(rows) if r["station"] == "hod2018" and r["p_ct"] not in ("", "nan")}


def spurious(prefix, ecut):
    S = {}
    for r in csv.DictReader(open(prefix + "_spurious.csv")):
        S.setdefault(int(r["ch"]), {"N05": int(r["N_gt0p5"]), "N15": int(r["N_gt1p5"]), "S15": 0, "S05": 0})
        if 3 <= int(r["delta"]) <= 6:
            S[int(r["ch"])]["S15"] += int(r["n_spurious_A2gt1p5"]) if ecut == 0.5 else 0
    for r in csv.DictReader(open(prefix + "_E_hist.csv")):
        if 3 <= int(r["delta"]) <= 6 and float(r["bin_lo_pe"]) < ecut - 1e-9:
            S[int(r["ch"])]["S05"] += int(r["n"])
    return S


def main():
    T = table(sys.argv[1])
    out = []
    for ecut in (0.5, 0.75):
        S = spurious(sys.argv[2], ecut)
        fmed = st.median(s["S05"] / s["N05"] for s in S.values())
        for c, r in sorted(T.items()):
            p = float(r["p_ct"])
            dcr = float(r["DCR_kHz"]) * 1e3
            if c in S:
                s = S[c]
                f = s["S05"] / s["N05"]
                pc = (s["N15"] - s["S15"]) / (s["N05"] - s["S05"])
                src = "medido"
            else:
                f = fmed
                pc = p / (1 - f)
                src = "mediana"
            pp = 1 - math.exp(-dcr * (1 - f) * TAU)
            pcc = (pc - pp) / (1 - pp)
            out.append({"sample": sys.argv[3], "Ecut_pe": ecut, "ch": c, "f_source": src, "f_spurious": f, "p_ct": p, "p_ct_corr": pc,
                        "factor_corr": 1 / (1 - pc), "p_ct_corr_pile": pcc, "factor_corr_pile": 1 / (1 - pcc), "factor_orig": 1 / (1 - p)})
        sel = [o for o in out if o["Ecut_pe"] == ecut]
        meas = [o for o in sel if o["f_source"] == "medido"]
        print(f"{sys.argv[3]} Ecut={ecut}: f medido {', '.join('ch%d %.4f' % (o['ch'], o['f_spurious']) for o in meas)}; mediana f {fmed:.4f}")
        for k in ("p_ct", "p_ct_corr", "factor_orig", "factor_corr", "p_ct_corr_pile", "factor_corr_pile"):
            v = [o[k] for o in sel]
            print(f"   {k}: mediana {st.median(v):.4f}, mín {min(v):.4f}, máx {max(v):.4f} (n = {len(v)})")
        for o in meas:
            print(f"   ch{o['ch']}: p_ct {o['p_ct']:.4f} -> {o['p_ct_corr']:.4f} (factor {o['factor_orig']:.4f} -> {o['factor_corr']:.4f}); "
                  f"con apilamiento {o['p_ct_corr_pile']:.4f} (factor {o['factor_corr_pile']:.4f})")
    if len(sys.argv) > 4:
        with open(sys.argv[4], "w") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]))
            w.writeheader()
            for o in out:
                w.writerow({k: (f"{v:.5f}" if isinstance(v, float) else v) for k, v in o.items()})


if __name__ == "__main__":
    main()
