#!/usr/bin/python3 -I
# biasscan_report.py — DAQ2-14/15: análisis del barrido de tensión a partir de biasscan_results.json.
#  - por canal de hod2018: G (mV/p.e., medida) contra la tensión del canal (VMon en el run SW, que sigue
#    al V0Set con el desvío del canal), ajuste lineal G = k·(V − V_bd), V_bd (estimada, extrapolación) con
#    su error del ajuste y la sobretensión al nominal; canales con menos de tres puntos válidos, marcados;
#  - diafonía P(≥ 2 p.e.) ≈ tasa(> 1.5 p.e.) / tasa(> 0.5 p.e.) por canal y escalón;
#  - ruido sobre N p.e. contra la sobretensión (mediana de canales), eficiencia X·Y contra tensión y
#    umbral, auto-trigger por escalón contra el trigger externo, tasa oscura de hod2019 > 1 mV;
#  - registro de escrituras de HV (biasscan.jsonl) y tiempos.
# Salida: <out>/report_es.md y <out>/biasscan_analysis.json.  Uso: biasscan_report.py <dir de resultados> <jsonl>
import json, math, statistics, sys


def med(v):
    """Mediana o nan si no hay valores (DAQ2-16: bajo la ruptura no hay pulsos y el informe no debe caerse)."""
    v = list(v)
    return statistics.median(v) if v else float("nan")

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
lines = [f"# Barrido de tensión de hod2018 y auto-trigger con pares de grupos X·Y ({out.rstrip('/').split('/')[-1]})", ""]
lines.append(f"Escalones medidos: {[s['delta_V'] for s in steps]}; omitidos: {R.get('skipped')}; aborto: {R.get('aborted')}.")
lines.append(f"Restauración: {R.get('restore', {}).get('result')}; run de física al final: {R.get('physics')}.")
lines.append("")
lines.append("## Hora y temperatura de la placa A7030DP por escalón (única referencia de temperatura; sigue la del rack)")
lines.append("")
lines.append("| ΔV | inicio (UTC) | inicio (Chile) | placa (°C), lecturas | control del auto-trigger |")
lines.append("|---|---|---|---|---|")
for s in steps:
    t = [b["degC"] for b in s.get("board_temp", [])]
    c = s.get("selftrigger_check", {})
    verdict = "funciona" if c.get("working") else (f"no comprobable: {c.get('reason')}" if c.get("working") is None and "rate_Hz" in c else
                                                    f"NO funciona: {c.get('reason')}")
    rg = c.get("registers") or {}
    regs = (f"; releído 0x810C {rg['read'].get('global_trigger_mask_0x810C')}, umbrales ADC {rg['read'].get('group_trigger_threshold_adc')}, "
            f"distintos de lo escrito: {rg.get('mismatch') or 'ninguno'}") if rg.get("read") else ""
    cs = ((f"{c['pair']} a {c['T_pe']} p.e.: {c['rate_Hz']:.3f} Hz contra {c['pred_acc_Hz']:.3f} Hz previstos ({c['frac_of_pred']:.3g}); "
           if "rate_Hz" in c else "") + verdict + regs) if c else "–"
    lines.append(f"| {s['delta_V']:+.1f} | {s.get('t_start_utc')} | {s.get('t_start_chile', '–')} | {t} | {cs} |")
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
_new = any("f_mix" in (v if isinstance(v, dict) else {}) for s in steps for v in s.get("ext", {}).get("all", {}).values())
lines.append("ε_XY = fracción X·Y con el estimador compartido (xy_estimator.py) / 0.0360; impacto = pulso en la ventana signal, sin ch3 (DAQ2-16)." if _new else
             "ε_XY = versión del barrido: impacto = mínimo de la ventana por grupo y estimador sin rX·rY ni (1−aX)(1−aY) (anterior a DAQ2-16); "
             "los valores corregidos están en \"Recálculo offline (DAQ2-16)\".")
lines.append("")
lines.append("| ΔV | G mediana (mV) | diafonía mediana P(>1.5)/P(>0.5) | ruido mediano > 5 / 8 / 10 / 12 p.e. por canal (Hz) | ε_XY ctl 5/8/10/12 | ε_XY mezclados 5/8/10/12 | hod2019 > 1 mV mediana (Hz) |")
lines.append("|---|---|---|---|---|---|---|")
for s in steps:
    pc = s["sw"]["per_channel"]
    get = lambda ch: pc.get(str(ch), pc.get(ch))
    xt = [get(ch)["rate_1.5pe_Hz"] / get(ch)["rate_0.5pe_Hz"] for ch in range(32) if get(ch).get("rate_0.5pe_Hz")]
    nz = {T: med([get(ch)[f"rate_{T}pe_Hz"] for ch in range(32) if get(ch).get(f"rate_{T}pe_Hz") is not None]) for T in TH}
    h19 = [get(ch)["rate_1mV_Hz"] for ch in range(32, 64) if get(ch).get("rate_1mV_Hz") is not None]
    e = s["ext"]["all"]
    gete = lambda T, k: e.get(T, e.get(int(T), {})).get(k, float("nan"))
    lines.append(f"| {s['delta_V']:+.1f} | {s['sw']['G_median_mV']:.3f} | {med(xt):.3f} | " + " / ".join(f"{nz[T]:.4g}" for T in TH)
                 + " | " + "/".join(f"{gete(T, 'eff_ctl'):.3f}" for T in TH) + " | " + "/".join(f"{gete(T, 'eff_mix'):.3f}" for T in TH)
                 + f" | {med(h19):.4g} |")
    ana["steps"].append({"delta_V": s["delta_V"], "crosstalk_median": med(xt), "noise_median_Hz": nz, "hod2019_1mV_median_Hz": med(h19)})
lines.append("")
lines.append("## Auto-trigger con un par de grupos (G_X ∧ G_Y, mayoría 2 en 8 ns, ch2/9/15 fuera)")
lines.append("")
lines.append("| ΔV | T (p.e.) | criterio | par | tasa (Hz) | X·Y en tiempo | accidentales previstos (2W·r_X·r_Y) | fuera de tiempo | razón | muones/h |")
lines.append("|---|---|---|---|---|---|---|---|---|---|")
for s in steps:
    for a in s.get("self", []):
        lines.append(f"| {s['delta_V']:+.1f} | {s['T_chosen_pe']} | {s['T_choice']} | {a['pair']} | {a['rate_Hz']:.3f} | {a['N_XY_in_time']} | {a['acc_pred']:.1f} | "
                     f"{a['acc_offtime']} | {a['ratio_pred']:.2f} | {a['muons_per_h']:.1f} |")
    mu = s.get("self_muons_per_h_sum")
    lines.append(f"| {s['delta_V']:+.1f} | | | suma | | | | | | {f'{mu:.1f}' if mu is not None else 'runs omitidos'} |")
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
