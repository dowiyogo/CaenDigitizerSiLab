#!/usr/bin/python3 -I
# biasscan_report.py — DAQ2-14: análisis del barrido de tensión a partir de biasscan_results.json.
#  - por canal de hod2018: G (mV/p.e., medida) contra la tensión del canal (VMon en el run SW, que sigue
#    al V0Set con el desvío del canal), ajuste lineal G = k·(V − V_bd), V_bd (estimada, extrapolación) con
#    su error del ajuste y la sobretensión al nominal; canales con menos de tres puntos válidos, marcados;
#  - diafonía P(≥ 2 p.e.) ≈ tasa(> 1.5 p.e.) / tasa(> 0.5 p.e.) por canal y escalón;
#  - ruido sobre N p.e. contra la sobretensión (mediana de canales), eficiencia X·Y contra tensión y
#    umbral, auto-trigger por escalón contra el trigger externo, tasa oscura de hod2019 > 1 mV;
#  - registro de escrituras de HV (biasscan.jsonl) y tiempos.
# Salida: <out>/report_es.md y <out>/biasscan_analysis.json.  Uso: biasscan_report.py <dir de resultados> <jsonl>
import json, math, statistics, sys

out, jl = sys.argv[1], sys.argv[2]
R = json.load(open(out + "/biasscan_results.json"))
steps = R["steps"]
TH = ["5", "8", "10", "12"]


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    k = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    b = my - k * mx
    res = [y - (k * x + b) for x, y in zip(xs, ys)]
    s2 = sum(r * r for r in res) / (n - 2) if n > 2 else float("nan")
    vk = s2 / sxx
    vb = s2 * (1 / n + mx * mx / sxx)
    cov = -mx * s2 / sxx
    vbd = -b / k
    # σ(V_bd) por propagación de (k, b) con su covarianza
    dk, db = b / (k * k), -1 / k
    svbd = math.sqrt(max(0.0, dk * dk * vk + db * db * vb + 2 * dk * db * cov)) if n > 2 else float("nan")
    return k, b, vbd, svbd


ana = {"channels": {}, "steps": []}
lines = ["# DAQ2-14: barrido de tensión de hod2018 y auto-trigger con pares de grupos X·Y", ""]
lines.append(f"Escalones medidos: {[s['delta_V'] for s in steps]}; omitidos: {R.get('skipped')}; aborto: {R.get('aborted')}.")
lines.append(f"Restauración: {R.get('restore', {}).get('result')}; run de física al final: {R.get('physics')}.")
lines.append("")
lines.append("## G contra tensión y tensión de ruptura (G medida; V_bd estimada por extrapolación lineal)")
lines.append("")
lines.append("| canal | puntos válidos | k (mV/p.e./V) | V_bd (V) | σ(V_bd) | sobretensión al nominal (V) |")
lines.append("|---|---|---|---|---|---|")
ch_hv = lambda ch: [0, 1, 2, 3, 4, 5, 6, 7]  # el grupo HV se toma del mapeo en el informe; aquí la tensión del canal HV
mapping = json.load(open("/home/xplus/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r4.json"))
hvch = {c["dig_ch"]: c["hv_ch"] for c in mapping["channels"]}
for ch in range(32):
    pts = []
    for s in steps:
        p = s["sw"]["per_channel"][str(ch)] if str(ch) in s["sw"]["per_channel"] else s["sw"]["per_channel"].get(ch)
        h = hvch.get(ch)
        if p and p.get("G_measured") and h is not None:
            v = s["hv_before_sw"][str(h)]["VMon_V"] if str(h) in s["hv_before_sw"] else s["hv_before_sw"][h]["VMon_V"]
            pts.append((v, p["G_mV"]))
    row = {"points": pts}
    if len(pts) >= 3:
        k, b, vbd, svbd = fit([a for a, _ in pts], [g for _, g in pts])
        vnom = max(a for a, _ in pts)
        row.update(k=k, V_bd=vbd, sigma_V_bd=svbd, overvoltage_nominal=vnom - vbd)
        lines.append(f"| {ch} | {len(pts)} | {k:.3f} | {vbd:.2f} | {svbd:.2f} | {vnom - vbd:.2f} |")
    else:
        lines.append(f"| {ch} | {len(pts)} (menos de 3: sin ajuste) | – | – | – | – |")
    ana["channels"][ch] = row
lines.append("")
lines.append("## Por escalón")
lines.append("")
lines.append("| ΔV | G mediana (mV) | diafonía mediana P(>1.5)/P(>0.5) | ruido mediano > 5 / 8 / 10 / 12 p.e. por canal (Hz) | ε_XY ctl 5/8/10/12 | ε_XY mezclados 5/8/10/12 | hod2019 > 1 mV mediana (Hz) |")
lines.append("|---|---|---|---|---|---|---|")
for s in steps:
    pc = s["sw"]["per_channel"]
    get = lambda ch: pc.get(str(ch), pc.get(ch))
    xt = [get(ch)["rate_1.5pe_Hz"] / get(ch)["rate_0.5pe_Hz"] for ch in range(32) if get(ch).get("rate_0.5pe_Hz")]
    nz = {T: statistics.median([get(ch)[f"rate_{T}pe_Hz"] for ch in range(32) if get(ch).get(f"rate_{T}pe_Hz") is not None]) for T in TH}
    h19 = [get(ch)["rate_1mV_Hz"] for ch in range(32, 64) if get(ch).get("rate_1mV_Hz") is not None]
    e = s["ext"]["all"]
    gete = lambda T, k: e.get(T, e.get(int(T), {})).get(k, float("nan"))
    lines.append(f"| {s['delta_V']:+.1f} | {s['sw']['G_median_mV']:.3f} | {statistics.median(xt):.3f} | " + " / ".join(f"{nz[T]:.4g}" for T in TH)
                 + " | " + "/".join(f"{gete(T, 'eff_ctl'):.3f}" for T in TH) + " | " + "/".join(f"{gete(T, 'eff_mix'):.3f}" for T in TH)
                 + f" | {statistics.median(h19):.4g} |")
    ana["steps"].append({"delta_V": s["delta_V"], "crosstalk_median": statistics.median(xt), "noise_median_Hz": nz, "hod2019_1mV_median_Hz": statistics.median(h19)})
lines.append("")
lines.append("## Auto-trigger con un par de grupos (G_X ∧ G_Y, mayoría 2 en 8 ns, ch2/9/15 fuera)")
lines.append("")
lines.append("| ΔV | T (p.e.) | criterio | par | tasa (Hz) | X·Y en tiempo | accidentales previstos (2W·r_X·r_Y) | fuera de tiempo | razón | muones/h |")
lines.append("|---|---|---|---|---|---|---|---|---|---|")
for s in steps:
    for a in s["self"]:
        lines.append(f"| {s['delta_V']:+.1f} | {s['T_chosen_pe']} | {s['T_choice']} | {a['pair']} | {a['rate_Hz']:.3f} | {a['N_XY_in_time']} | {a['acc_pred']:.1f} | "
                     f"{a['acc_offtime']} | {a['ratio_pred']:.2f} | {a['muons_per_h']:.1f} |")
    lines.append(f"| {s['delta_V']:+.1f} | | | suma | | | | | | {s['self_muons_per_h_sum']:.1f} |")
lines.append("")
lines.append("Referencia, trigger externo: fracción X·Y con accidentales restados 0.0167 (noche del 2026-10-19, eventos mezclados) × tasa del trigger "
             "externo (0.7299 Hz) = 43.9 muones/h en el cuadrado dentro de la aceptancia de las placas.")
lines.append("")
lines.append("## Escrituras de HV (biasscan.jsonl)")
lines.append("")
for l in open(jl):
    j = json.loads(l)
    if j.get("ev") == "v0set" or "bias " in j.get("msg", ""):
        lines.append(f"- {j['t_utc']} {j['msg']}")
open(out + "/report_es.md", "w").write("\n".join(lines) + "\n")
json.dump(ana, open(out + "/biasscan_analysis.json", "w"), indent=1, default=str)
print("\n".join(lines[:60]))
