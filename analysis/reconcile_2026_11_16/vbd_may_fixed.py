#!/usr/bin/python3 -I
# vbd_may_fixed.py — seguimiento 2026-11-16, punto 4: corrección del V_bd de mayo del seguimiento 15.4 (hv_may_oct.py, sección 3; copia sin cambios en
# hv_may_oct_f15_original.py).
# Dónde estaba el problema: fit_may() ajusta I = I0 + a·(V − V_bd)^n a los tres puntos GECO de mayo con a libre. Con a libre, la escala absoluta de
# la corriente se pierde y V_bd queda fijado solo por la pendiente logarítmica d ln(I − I0)/dV = n/(V − V_bd) entre los tres puntos, tomados el 25,
# 27 y 28 de mayo (dos I0Set distintos, temperatura no registrada). La pendiente de mayo entre 65.08 y 65.54 V es mayor que la de octubre, y eso da
# V_bd de mayo mayor que el de octubre: el signo contrario al de la ganancia (G_mayo > G_oct con ΔV ≈ 0 exige V_bd_mayo < V_bd_oct). No es un signo
# invertido en una fórmula: las fórmulas ΔV_bd = −b/k (ganancia) y ΔV_bd = V_bd_mayo − V_bd_oct (I-V) son correctas; el estimador I-V con a libre
# no usa la corriente absoluta, que es la información que fija el signo.
# Corrección: comparar la corriente del mismo grupo en mayo (lectura del 2026-05-28 01:10, la más cercana a run B) y en octubre run 8 (mediana de
# hvmon), con la forma de octubre del grupo (I0, n del ajuste b_free, o b_fixed si b_free está en el borde; seguimiento 8):
#   OV_mayo / OV_oct = ((I_mayo − I0) / (I_oct − I0))^(1/n)   (a se cancela),   OV_oct = VMon_oct − V_bd_oct,
#   ΔV_bd(mayo − oct) = ΔV − OV_oct · (OV_mayo/OV_oct − 1).
# Supuesto [S]: corriente por p.e. y tasa oscura a OV fijo iguales en las dos campañas. Si mayo era más frío, la tasa a OV fijo era menor y el ΔV_bd
# verdadero es más negativo que el de aquí (este valor es una cota superior).
# Intervalo: variando n en su intervalo de octubre y V_bd_oct en su intervalo (cada uno por separado; se da el rango).
# Salida: <out>/vbd_may_fixed.csv, vbd_may_fixed.json. Uso: vbd_may_fixed.py <out dir>
import csv, json, os, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = sys.argv[1]
IV = json.load(open(REPO + "/outputs/hod2018_iv/2026-10-31_iv/iv_fit.json"))["groups"]
HV = {int(r["hv_group"]): r for r in csv.DictReader(l for l in open(REPO + "/outputs/hod2018_hv_may_oct/2026-11-15/hv_delta.csv") if not l.startswith("#"))}
DVBD_DT = 0.060  # V/°C [S], el del seguimiento 15.4
rows, res = [], {}
for g in range(1, 9):
    G = IV[str(g)]
    ref, src = (G["b_free"], "b_free") if not G["b_free"]["at_grid_edge"] else ((G["b_fixed"], "b_fixed") if "b_fixed" in G else (None, "none"))
    h = HV[g]
    vm, vo, im, io = float(h["VMon_may_V"]), float(h["VMon_oct_V"]), float(h["IMon_may_uA"]), float(h["IMon_oct_uA"])
    if ref is None:
        rows.append([g, src, "", "", "", "", "", "", "", "", "no October I-V shape"])
        continue
    i0 = G["I0_floor_uA"]

    def dvbd(n, vb):
        ovo = vo - vb
        r = ((im - i0) / (io - i0)) ** (1.0 / n)
        return (vm - vo) - ovo * (r - 1.0), ovo, r * ovo

    d0, ovo, ovm = dvbd(ref["n"], ref["V_bd_V"])
    alts = [dvbd(ref["n"] + ref["n_err"][0], ref["V_bd_V"])[0], dvbd(ref["n"] + ref["n_err"][1], ref["V_bd_V"])[0],
            dvbd(ref["n"], ref["V_bd_V"] + ref["V_bd_err_V"][0])[0], dvbd(ref["n"], ref["V_bd_V"] + ref["V_bd_err_V"][1])[0]]
    lo, hi = min(alts + [d0]), max(alts + [d0])
    res[g] = {"source": src, "n": ref["n"], "V_bd_oct_V": ref["V_bd_V"], "I0_uA": i0, "VMon_may_V": vm, "VMon_oct_V": vo, "IMon_may_uA": im, "IMon_oct_uA": io,
              "OV_oct_V": ovo, "OV_may_V": ovm, "dVbd_may_minus_oct_V": d0, "range_V": [lo, hi], "V_bd_may_V": ref["V_bd_V"] + d0, "dT_C": d0 / DVBD_DT}
    rows.append([g, src, f"{ref['n']:.3f}", f"{ref['V_bd_V']:.3f}", f"{ovo:.3f}", f"{ovm:.3f}", f"{ref['V_bd_V'] + d0:.3f}", f"{d0:+.3f}", f"{lo:+.3f}", f"{hi:+.3f}",
                 f"{d0 / DVBD_DT:+.2f}"])
    print(f"grupo {g} ({src}): OV oct {ovo:.3f} V, OV mayo {ovm:.3f} V, ΔV_bd {d0:+.3f} V [{lo:+.3f}, {hi:+.3f}], ΔT {d0 / DVBD_DT:+.2f} °C", flush=True)
with open(OUT + "/vbd_may_fixed.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hv_group", "oct_shape", "n", "V_bd_oct_V", "OV_oct_V", "OV_may_V", "V_bd_may_V", "dVbd_may_minus_oct_V", "range_lo_V", "range_hi_V", "dT_may_minus_oct_C"])
    w.writerows(rows)
json.dump(res, open(OUT + "/vbd_may_fixed.json", "w"), indent=1)
