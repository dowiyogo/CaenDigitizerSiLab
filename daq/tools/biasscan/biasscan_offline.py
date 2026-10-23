#!/usr/bin/python3 -I
# biasscan_offline.py — DAQ2-16 (solo lectura de los datos crudos): recálculo offline de un barrido ya hecho con el código
# actual de biasscan.py (se importan sus funciones, no se copian):
#  - G por canal y escalón con la ventana que sigue a la G esperada (analyse_sw: [0.5, 6] × G esperada, piso 2.5 mV,
#    semillas reducidas si no hay G aceptada; la esperada es la del escalón anterior o c1 en ΔV = 0);
#  - eficiencia X·Y con la definición única de impacto y el estimador compartido (analyse_ext, xy_estimator.py), junto a
#    la versión del barrido (impacto = mínimo de la ventana por grupo, estimador sin rX·rY ni (1−aX)(1−aY));
#  - V_bd por canal (estimada): recta G contra VMon del canal HV del canal (mapeo r4) con los escalones con G medida.
# Salida: <dir del barrido>/offline_daq2-16/ (calibraciones, impactos, offline.json) y la sección "Recálculo offline
# (DAQ2-16)" agregada al final de <dir del barrido>/report_es.md.
# Uso: biasscan_offline.py <outputs/daq2-15/<campaña>>
import importlib.util, json, math, os, statistics, sys

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
src = os.path.abspath(sys.argv[1])
campaign = os.path.basename(src.rstrip("/"))
sys.argv = ["biasscan.py", "--campaign", campaign]
spec = importlib.util.spec_from_file_location("bs", REPO + "/daq/tools/biasscan/biasscan.py")
bs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bs)
OFF = src + "/offline_daq2-16"
os.makedirs(OFF, exist_ok=True)
bs.OUT, bs.LOG, bs.JL = OFF, OFF + "/offline.log", OFF + "/offline.jsonl"
R = json.load(open(src + "/biasscan_results.json"))
MAP = json.load(open(bs.MAP))
HVCH = {c["dig_ch"]: c["hv_ch"] for c in MAP["channels"]}
TT = (5, 8, 10, 12)

steps, prev = [], bs.C1
for i, s in enumerate(R["steps"]):
    sw = bs.analyse_sw(s["sw_run"], i, prev)
    prev = sw["calibration"]
    ext = bs.analyse_ext(s["ext_run"], i, prev)
    cal = json.load(open(prev))
    old_cal = json.load(open(f"{src}/pecal_step{i}.json"))
    st = {"delta_V": s["delta_V"], "t_start_utc": s.get("t_start_utc"), "hv": s["hv_before_sw"],
          "G_new": {e["ch"]: {"G_mV": e.get("G_mV"), "status": e.get("status"), "measured": e.get("status") in ("ok", "changed") and "G_mV_fallback" not in e,
                              "seed_mV": e.get("seed_mV"), "seed_used_mV": e.get("seed_used_mV"), "fit_window_mV": e.get("fit_window_mV")}
                    for e in cal["channels"] if e["ch"] < 32},
          "G_old": {e["ch"]: {"G_mV": e.get("G_mV"), "status": e.get("status")} for e in old_cal["channels"] if e["ch"] < 32},
          "G_new_median_mV": sw["G_median_mV"], "n_G_new": sw["n_G_ok"],
          "eff_old": {str(T): {k: s["ext"]["all"].get(str(T), s["ext"]["all"].get(T, {})).get(k) for k in ("N_XY", "acc_mix", "acc_ctl", "eff_mix", "eff_ctl")}
                      for T in TT},
          "eff_new": {str(T): {k: ext["all"][T].get(k) for k in ("N_XY", "acc_mix", "acc_ctl", "f_mix", "f_mix_err", "f_ctl", "f_ctl_err", "eff_mix", "eff_ctl",
                                                                  "e_X|Y_mix", "e_Y|X_mix")} for T in TT},
          "eff_new_masked": {str(T): {k: ext["masked"][T].get(k) for k in ("f_mix", "f_mix_err", "f_ctl", "eff_mix")} for T in TT}}
    steps.append(st)
    json.dump({"campaign": campaign, "steps": steps}, open(OFF + "/offline.json", "w"), indent=1, default=str)


def vbd(ch):
    pts = []
    for s in steps:
        g = s["G_new"][ch]
        h = s["hv"].get(str(HVCH.get(ch))) or s["hv"].get(HVCH.get(ch))
        if g["measured"] and h:
            pts.append((h["VMon_V"], g["G_mV"]))
    if len(pts) < 2:
        return {"points": pts}
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    n = len(pts)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    k = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    b = my - k * mx
    d = {"points": pts, "k_mV_per_V": k, "V_bd_V": -b / k, "overvoltage_nominal_V": max(xs) - (-b / k), "sigma_V_bd_V": None}
    if n > 2:
        s2 = sum((y - (k * x + b)) ** 2 for x, y in zip(xs, ys)) / (n - 2)
        vk, vb, cov = s2 / sxx, s2 * (1 / n + mx * mx / sxx), -mx * s2 / sxx
        dk, db = b / (k * k), -1 / k
        d["sigma_V_bd_V"] = math.sqrt(max(0.0, dk * dk * vk + db * db * vb + 2 * dk * db * cov))
    return d


V = {ch: vbd(ch) for ch in range(32)}
json.dump({"campaign": campaign, "steps": steps, "V_bd": V}, open(OFF + "/offline.json", "w"), indent=1, default=str)

L = ["", "## Recálculo offline (DAQ2-16)", "",
     "Mismos datos crudos, código de DAQ2-16 (`biasscan_offline.py`). \"Barrido\" = lo que calculó el barrido (impacto = mínimo de la ventana por "
     "grupo; estimador sin rX·rY ni (1−aX)(1−aY)); \"corregido\" = definición única de impacto (pulso en la ventana signal, sin ch3) y estimador "
     "compartido (`xy_estimator.py`). ε = fracción / 0.0360.", "",
     "### Eficiencia X·Y (ε con eventos mezclados; entre paréntesis, fuera de tiempo)", "",
     "| ΔV | " + " | ".join(f"{T} p.e. barrido" for T in TT) + " | " + " | ".join(f"{T} p.e. corregido" for T in TT) + " |",
     "|---|" + "---|" * (2 * len(TT))]
for s in steps:
    o = [f"{s['eff_old'][str(T)]['eff_mix']:.3f} ({s['eff_old'][str(T)]['eff_ctl']:.3f})" if s["eff_old"][str(T)]["eff_mix"] is not None else "–" for T in TT]
    nw = [f"{s['eff_new'][str(T)]['eff_mix']:.3f} ± {s['eff_new'][str(T)]['f_mix_err'] / bs.ACC_GEOM:.3f} ({s['eff_new'][str(T)]['eff_ctl']:.3f})" for T in TT]
    L.append(f"| {s['delta_V']:+.1f} | " + " | ".join(o) + " | " + " | ".join(nw) + " |")
L += ["", "### G por canal (mV/p.e.; \"–\" = sin G aceptada; * = changed, a revisar)", "",
      "| canal | " + " | ".join(f"{s['delta_V']:+.1f} V barrido | {s['delta_V']:+.1f} V nuevo" for s in steps) + " | V_bd (V, estimada) | sobretensión al nominal (V) |",
      "|---|" + "---|" * (2 * len(steps) + 2)]
for ch in range(32):
    cells = []
    for s in steps:
        go, gn = s["G_old"][ch], s["G_new"][ch]
        cells.append(f"{go['G_mV']:.3f}" if go.get("G_mV") and go.get("status") in ("ok", "changed") else "–")
        cells.append(f"{gn['G_mV']:.3f}{'*' if gn['status'] == 'changed' else ''}" if gn["measured"] else "–")
    v = V[ch]
    vb = (f"{v['V_bd_V']:.2f}" + (f" ± {v['sigma_V_bd_V']:.2f}" if v.get("sigma_V_bd_V") else " (2 puntos)")) if "V_bd_V" in v else "–"
    ov = f"{v['overvoltage_nominal_V']:.2f}" if "V_bd_V" in v else "–"
    L.append(f"| {ch} | " + " | ".join(cells) + f" | {vb} | {ov} |")
L += ["", "Mediana de G nueva por escalón: " + ", ".join(f"{s['delta_V']:+.1f} V: " + (f"{s['G_new_median_mV']:.3f} mV ({s['n_G_new']} canales)" if s["n_G_new"]
                                                      else "sin G (ningún canal con picos)") for s in steps) + ".",
      "V_bd es estimada: recta G contra VMon con los escalones con G medida; con 2 puntos no hay error del ajuste."]
vb_all = [V[ch]["V_bd_V"] for ch in range(32) if "V_bd_V" in V[ch]]
if vb_all:
    L.append(f"V_bd en {len(vb_all)} canales: mediana {statistics.median(vb_all):.2f} V; sobretensión al nominal mediana "
             f"{statistics.median(V[ch]['overvoltage_nominal_V'] for ch in range(32) if 'V_bd_V' in V[ch]):.2f} V.")
open(src + "/report_es.md", "a").write("\n".join(L) + "\n")
print("\n".join(L))
