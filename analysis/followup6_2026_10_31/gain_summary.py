#!/usr/bin/python3 -I
# gain_summary.py — seguimiento 2026-10-31 (6), punto 3: resumen de fit_gain_long.C.
#  - hod2018 (control): G de gain_fit.h con la semilla propia (reajuste de PROC) en la región larga contra la G del reajuste de PROC.
#  - hod2019: por canal y configuración, cuántas semillas de la rejilla dan un ajuste aceptado (ok/changed) y con qué G; "converge" = aceptado con
#    al menos 3 semillas de la rejilla y G de esas semillas dentro de ±5 % de su mediana (ajuste estable frente a la semilla).
#  - Cola del espectro oscuro: razón de cuentas entre bins de 1 mV consecutivos en 2-7 mV (mediana de las 5 razones) y, para un espectro de p.e.
#    con crosstalk geométrico (cuentas ∝ p_ct^n en n·G), la G que daría esa caída con p_ct = 0.1, 0.3 y 0.5 (G = ln p_ct / ln razón, por mV)
#    [inferido: supone que la cola es la serie de p.e.]; y la caída de una cola gaussiana con el RMS de mayo de hod2019 entre 3 y 4 mV.
# Uso: gain_summary.py <prefijo de fit_gain_long.C> <G_may_refit.txt> <salida.json>
import csv, json, math, statistics as st, sys

import numpy as np

pre = sys.argv[1]
gref = [float(x) for x in open(sys.argv[2]).read().strip().split(",")]
fits = list(csv.DictReader(open(pre + "_fits.csv")))
rootf = list(csv.DictReader(open(pre + "_rootfit.csv")))
spec = np.loadtxt(pre + "_spec.csv", delimiter=",", skiprows=1)
out = {"hod2018_control": [], "hod2019": []}
for c in range(32):
    if c == 3:
        continue
    r = next(f for f in fits if int(f["ch"]) == c and f["config"] == "proc" and f["seed_kind"] == "own")
    rf = next((f for f in rootf if int(f["ch"]) == c and f["config"] == "proc"), None)
    out["hod2018_control"].append({"ch": c, "G_proc_refit_mV": gref[c], "G_long_mV": float(r["G_mV"]), "status": r["status"], "peaks": int(r["peaks"]),
                                   "ratio": float(r["G_mV"]) / gref[c],
                                   "rootfit": {k: rf[k] for k in ("G_mV", "G_err", "chi2", "ndf", "fit_status")} if rf else None})
rat = [x["ratio"] for x in out["hod2018_control"] if x["status"] in ("ok", "changed")]
out["hod2018_control_summary"] = {"n_accepted": len(rat), "ratio_median": st.median(rat), "ratio_min": min(rat), "ratio_max": max(rat),
                                  "rootfit_chi2_ndf_range": [min(float(x["rootfit"]["chi2"]) / int(x["rootfit"]["ndf"]) for x in out["hod2018_control"] if x["rootfit"]),
                                                             max(float(x["rootfit"]["chi2"]) / int(x["rootfit"]["ndf"]) for x in out["hod2018_control"] if x["rootfit"])],
                                  "rootfit_status0": sum(1 for x in out["hod2018_control"] if x["rootfit"] and x["rootfit"]["fit_status"] == "0")}
for c in range(32, 64):
    d = {"ch": c}
    for cfg in ("proc", "low"):
        rows = [f for f in fits if int(f["ch"]) == c and f["config"] == cfg and f["seed_kind"] == "grid"]
        acc = [(float(f["seed_mV"]), float(f["G_mV"]), f["status"]) for f in rows if f["status"] in ("ok", "changed")]
        stat = {}
        for f in rows:
            stat[f["status"]] = stat.get(f["status"], 0) + 1
        g = [a[1] for a in acc]
        conv = len(g) >= 3 and all(abs(x / st.median(g) - 1) <= 0.05 for x in g)
        rf = next((f for f in rootf if int(f["ch"]) == c and f["config"] == cfg), None)
        d[cfg] = {"seeds": len(rows), "status_counts": stat, "accepted": acc, "converges": conv,
                  "rootfit": {k: rf[k] for k in ("G_mV", "G_err", "chi2", "ndf", "fit_status")} if rf else None}
    h = spec[spec[:, 0] == c][:, 2]
    mv = h[:400].reshape(100, 4).sum(1)  # bins de 1 mV
    rs = [mv[i + 1] / mv[i] for i in range(2, 7) if mv[i] > 0]
    rr = st.median(rs) if rs else float("nan")
    d["tail"] = {"counts_per_mV_0_10": mv[:10].astype(int).tolist(), "ratio_per_mV_2to7": rr, "n_above_10mV": int(h[40:].sum()), "n_total": int(h.sum()),
                 "G_if_pe_series_mV": {str(p): (math.log(p) / math.log(rr) if 0 < rr < 1 else None) for p in (0.1, 0.3, 0.5)}}
    out["hod2019"].append(d)
out["gaussian_tail_ratio_3_to_4mV_rms0p79"] = math.exp(-(4 ** 2 - 3 ** 2) / (2 * 0.79 ** 2))
out["hod2019_summary"] = {"converging_proc": [d["ch"] for d in out["hod2019"] if d["proc"]["converges"]],
                          "converging_low": [d["ch"] for d in out["hod2019"] if d["low"]["converges"]],
                          "any_accepted_proc": {d["ch"]: d["proc"]["accepted"] for d in out["hod2019"] if d["proc"]["accepted"]},
                          "any_accepted_low": {d["ch"]: d["low"]["accepted"] for d in out["hod2019"] if d["low"]["accepted"]},
                          "tail_ratio_median_over_channels": st.median(d["tail"]["ratio_per_mV_2to7"] for d in out["hod2019"]),
                          "tail_ratio_range": [min(d["tail"]["ratio_per_mV_2to7"] for d in out["hod2019"]), max(d["tail"]["ratio_per_mV_2to7"] for d in out["hod2019"])]}
json.dump(out, open(sys.argv[3], "w"), indent=1)
s = out["hod2018_control_summary"]
print(f"hod2018 control: {s['n_accepted']} aceptados, G_larga/G_PROC mediana {s['ratio_median']:.4f} ({s['ratio_min']:.4f}-{s['ratio_max']:.4f}); "
      f"ROOT χ²/ndf {s['rootfit_chi2_ndf_range'][0]:.1f}-{s['rootfit_chi2_ndf_range'][1]:.1f}, estado 0 en {s['rootfit_status0']}")
for d in out["hod2019"]:
    t = d["tail"]
    print(f"ch{d['ch']}: proc {d['proc']['status_counts']} aceptados {d['proc']['accepted']} converge {d['proc']['converges']} | low {d['low']['status_counts']} "
          f"converge {d['low']['converges']} | cola razón/mV {t['ratio_per_mV_2to7']:.3f}, >10 mV {t['n_above_10mV']} de {t['n_total']}, "
          f"G si serie de p.e. (p_ct 0.1/0.3/0.5) {', '.join('%.2f' % v if v else '-' for v in t['G_if_pe_series_mV'].values())} mV")
print(json.dumps(out["hod2019_summary"]), f"cola gaussiana RMS 0.79 mV, razón 4/3 mV: {out['gaussian_tail_ratio_3_to_4mV_rms0p79']:.4f}")
