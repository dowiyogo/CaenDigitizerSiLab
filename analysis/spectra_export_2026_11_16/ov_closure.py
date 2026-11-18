#!/usr/bin/python3 -I
# ov_closure.py — seguimiento 2026-11-16 (spectra-export), punto 2: contabilidad del sobrevoltaje por grupo HV de hod2018, octubre − mayo.
# Entradas (solo lectura):
#  - ΔV_bd: outputs/hod2018_ratio_reconcile/2026-11-16/vbd_may_fixed.csv (cociente de corrientes con la forma I-V de octubre; ahí es mayo − octubre,
#    aquí se cambia el signo). Error: el rango de esa tabla (n y V_bd de octubre en sus intervalos), asimétrico; para χ² se usa la mitad del lado
#    hacia el residuo [supuesto].
#  - ΔVMon: outputs/hod2018_hv_may_oct/2026-11-15/hv_delta.csv (VMon de la captura GECO del 2026-05-28 01:10 y mediana de hvmon en run 8); aquí
#    octubre − mayo. Error 0.01 V (resolución de la pantalla GECO) [supuesto].
#  - G_oct/G_mayo por MPPC: G del análisis (mayo = reajuste de PROC, octubre = c1), barras emparejadas por MPPC (seguimiento 7, ratio_by_mppc.csv);
#    fuera ch3 y los MPPC con G de octubre de relleno (ch22, ch31). Por grupo: mediana y error 1.2533·σ/√n (σ = desviación entre MPPC del grupo).
# ΔOV_gain = OV_mayo · (G_oct/G_mayo − 1), OV_mayo = 1.0, 1.2, 1.4 V [supuesto]; cierre: ΔOV_gain = ΔVMon − ΔV_bd; pull = (ΔOV_gain − (ΔVMon − ΔV_bd))/σ,
# σ² = σ_gain² + σ_VMon² + σ_Vbd²; χ² sobre los grupos con ΔV_bd. ΔT = ΔV_bd / 0.060 V/°C [hoja de datos, valor de la tarea]; también el ΔV_bd que
# exige la ganancia, ΔVMon − ΔOV_gain, y su ΔT.
# Salida: <out>/ov_closure.csv, ov_closure.json. Uso: ov_closure.py <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

OUT = sys.argv[1]
PK = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work"
VB = {int(r["hv_group"]): r for r in csv.DictReader(open(REPO + "/outputs/hod2018_ratio_reconcile/2026-11-16/vbd_may_fixed.csv"))}
HV = {int(r["hv_group"]): r for r in csv.DictReader(l for l in open(REPO + "/outputs/hod2018_hv_may_oct/2026-11-15/hv_delta.csv") if not l.startswith("#"))}
BARS = [r for r in csv.DictReader(open(REPO + "/outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv")) if r["plane"] in ("X", "Y")]
DVDT = 0.060
SIG_VMON = 0.01


def g_refit():
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


gm = g_refit()
gm = np.where(np.isfinite(gm), gm, np.nanmedian(gm))
c1 = mc.g_c1()[:32]
grat = {}
for b in BARS:
    cm, co = int(b["ch_may"]), int(b["ch_oct"])
    if cm >= 32 or 3 in (cm, co) or not np.isfinite(c1[co]):
        continue
    grat.setdefault(int(b["hv_group"]), []).append(c1[co] / gm[cm])
rows, res = [], {"groups": {}, "closure": {}}
for g in range(1, 9):
    v = np.array(grat[g])
    gr = float(np.median(v))
    egr = float(1.2533 * np.std(v) / math.sqrt(len(v))) if len(v) > 1 else float("nan")
    dvmon = -float(HV[g]["dVMon_V"])
    r = VB[g]
    if r["dVbd_may_minus_oct_V"]:
        dvbd = -float(r["dVbd_may_minus_oct_V"])
        lo, hi = -float(r["range_hi_V"]), -float(r["range_lo_V"])  # octubre − mayo
    else:
        dvbd, lo, hi = float("nan"), float("nan"), float("nan")
    res["groups"][g] = {"n_mppc": len(v), "G_oct_over_may": gr, "G_err": egr, "dVMon_oct_minus_may_V": dvmon, "dVbd_oct_minus_may_V": dvbd,
                        "dVbd_range_V": [lo, hi], "dT_IV_C": dvbd / DVDT if np.isfinite(dvbd) else None}
    for ov in (1.0, 1.2, 1.4):
        dog = ov * (gr - 1)
        edog = ov * egr
        need = dvmon - dvbd
        resid = dog - need
        # error de ΔV_bd hacia el lado del residuo: residuo > 0 pide ΔV_bd más negativo
        evb = (dvbd - lo) if resid > 0 else (hi - dvbd)
        sig = math.sqrt(edog ** 2 + SIG_VMON ** 2 + evb ** 2) if np.isfinite(dvbd) else float("nan")
        pull = resid / sig if np.isfinite(dvbd) else float("nan")
        dvbd_gain = dvmon - dog
        rows.append([g, len(v), f"{gr:.4f}", f"{egr:.4f}", ov, f"{dog:+.4f}", f"{edog:.4f}", f"{dvmon:+.3f}", f"{dvbd:+.3f}" if np.isfinite(dvbd) else "",
                     f"{lo:+.3f}…{hi:+.3f}" if np.isfinite(lo) else "", f"{need:+.4f}" if np.isfinite(need) else "", f"{resid:+.4f}" if np.isfinite(resid) else "",
                     f"{sig:.4f}" if np.isfinite(sig) else "", f"{pull:+.2f}" if np.isfinite(pull) else "", f"{dvbd / DVDT:+.2f}" if np.isfinite(dvbd) else "",
                     f"{dvbd_gain:+.4f}", f"{dvbd_gain / DVDT:+.2f}"])
        res["closure"].setdefault(str(ov), {})[g] = {"dOV_gain_V": dog, "dOV_gain_err_V": edog, "dVMon_minus_dVbd_V": need, "residual_V": resid, "sigma_V": sig,
                                                     "pull": pull, "dVbd_from_gain_V": dvbd_gain, "dT_from_gain_C": dvbd_gain / DVDT}
for ov in ("1.0", "1.2", "1.4"):
    c = res["closure"][ov]
    pulls = [x["pull"] for x in c.values() if np.isfinite(x["pull"])]
    res["closure"][ov]["chi2"] = float(sum(p * p for p in pulls))
    res["closure"][ov]["ndf"] = len(pulls)
    print(f"OV_mayo {ov} V: χ² del cierre {res['closure'][ov]['chi2']:.2f}/{len(pulls)}; pulls " + " ".join(f"{p:+.2f}" for p in pulls), flush=True)
with open(OUT + "/ov_closure.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["hv_group", "n_mppc", "G_oct_over_may", "err", "OV_may_V_assumed", "dOV_gain_V", "err", "dVMon_oct_minus_may_V", "dVbd_oct_minus_may_V",
                "dVbd_range_V", "dVMon_minus_dVbd_V", "residual_dOVgain_minus_(dVMon-dVbd)_V", "sigma_V", "pull", "dT_from_IV_C", "dVbd_required_by_gain_V",
                "dT_from_gain_C"])
    w.writerows(rows)
json.dump(res, open(OUT + "/ov_closure.json", "w"), indent=1, default=float)
for r in rows:
    print(" | ".join(map(str, r)))
