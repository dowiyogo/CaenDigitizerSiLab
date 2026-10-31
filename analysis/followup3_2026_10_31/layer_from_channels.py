#!/usr/bin/python3 -I
# layer_from_channels.py — seguimiento 2026-10-31 (3), punto 2: pendiente de cada capa como media ponderada de las pendientes exponenciales de sus
# canales (summary.json de amp_vs_position.py), con el error inflado por √(χ²/ndf) entre canales si es > 1. Más robusto que el ajuste a la suma
# de histogramas de la capa, que mezcla canales de amplitud distinta con pesos que cambian con la posición.
# Salida: <dir>/layer_from_channels.json. Uso: layer_from_channels.py <dir de salidas de amp_vs_position.py>
import json, math, sys

D = sys.argv[1]
S = json.load(open(D + "/summary.json"))
out = {}
for run, R in S["runs"].items():
    for L in ("XA", "XB", "YA", "YB"):
        ks, ws, chs = [], [], []
        for v in R["channels"].values():
            if v["coord"] + v["layer"] != L or not v["fit_mean"]:
                continue
            ex = v["fit_mean"]["exp"]
            ks.append(ex["slope_per_mm"]); ws.append(1 / ex["slope_err"] ** 2); chs.append(v["ch"])
        km = sum(k * w for k, w in zip(ks, ws)) / sum(ws)
        chi = sum(w * (k - km) ** 2 for k, w in zip(ks, ws))
        sc = math.sqrt(max(chi / (len(ks) - 1), 1))
        ek = sc / math.sqrt(sum(ws))
        out[f"{run}/{L}"] = {"slope_per_mm": km, "slope_err": ek, "L_eff_mm": abs(1 / km), "L_eff_err_mm": ek / km ** 2, "sign": int(math.copysign(1, km)),
                             "ratio_54_over_21": math.exp(33 * km), "ratio_err": math.exp(33 * km) * 33 * ek, "chi2_between_channels": chi, "ndf": len(ks) - 1,
                             "n_channels_positive": sum(1 for k in ks if k > 0), "channels": chs}
        o = out[f"{run}/{L}"]
        print(f"{run:20s} {L}: L_eff {o['L_eff_mm']:.0f} ± {o['L_eff_err_mm']:.0f} mm, signo {o['sign']:+d}, razón {o['ratio_54_over_21']:.3f} ± {o['ratio_err']:.3f}, "
              f"χ² {chi:.1f}/{len(ks) - 1}, + {o['n_channels_positive']}/{len(ks)}")
json.dump(out, open(D + "/layer_from_channels.json", "w"), indent=1)
