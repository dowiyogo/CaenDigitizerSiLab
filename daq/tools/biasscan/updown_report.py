#!/usr/bin/python3 -I
# updown_report.py — DAQ2-17 (solo lectura): informe de la mañana con el barrido de día (DAQ2-15, hacia abajo; G del
# recálculo offline de DAQ2-16) y el de noche (DAQ2-17, hacia arriba).
#  - por canal de hod2018: G contra VMon (día 0 y −0.5 V; noche 0 a +1.0 V; solo escalones con G medida), recta, V_bd ± error
#    (estimada), sobretensión al nominal; V_bd solo de día (2 puntos) y solo de noche, su diferencia, y la diferencia de la
#    temperatura de la placa;
#  - por escalón de la noche: IMon máxima por grupo, tasa oscura (> 0.5 / 1.5 / 3 / 5 / 8 p.e.) y diafonía contra la
#    sobretensión, eficiencia (fracción / 0.03646, ε(X|Y), ε(Y|X)) con errores, espectro en tiempo del plano X en eventos con Y
#    sobre 5 p.e. contra la ventana fuera de tiempo (¿se separa el pico de muones del ruido?);
#  - hod2019: tasa > 1 mV por canal contra ΔV (día y noche): qué canales responden a la tensión;
#  - diagnóstico del trigger (ΔV = 0 de la noche) y tabla de puntos de operación candidatos (sin elegir).
# Uso: updown_report.py <dir del barrido de día> <dir del barrido de noche> <dir de salida>
import json, math, os, statistics, sys

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import xy_estimator as xy  # noqa: E402

DAY, NIGHT, OUT = (os.path.abspath(a) for a in sys.argv[1:4])
os.makedirs(OUT, exist_ok=True)
ACC = 0.03646
MAP = json.load(open(os.path.expanduser("~/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r4.json")))
HVCH = {c["dig_ch"]: c["hv_ch"] for c in MAP["channels"]}
DOFF = json.load(open(DAY + "/offline_daq2-16/offline.json"))
DRES = json.load(open(DAY + "/biasscan_results.json"))
NRES = json.load(open(NIGHT + "/biasscan_results.json"))
NSTEPS = NRES["steps"]
figs = []
L = []


def f(x, d=3):
    return "–" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{d}f}"


def temps(s):
    t = [b["degC"] for b in s.get("board_temp", []) if b.get("degC") is not None]
    return statistics.mean(t) if t else float("nan")


def hv(s, ch, key="hv_before_sw"):
    h = s.get(key, {})
    return h.get(str(HVCH.get(ch))) or h.get(HVCH.get(ch))


def fit(pts):
    n = len(pts)
    if n < 2:
        return None
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx <= 0:
        return None
    k = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    b = my - k * mx
    if k <= 0:  # sin pendiente positiva no hay V_bd
        return None
    d = {"n": n, "k_mV_per_V": k, "V_bd_V": -b / k, "sigma_V_bd_V": None}
    if n > 2:
        s2 = sum((y - (k * x + b)) ** 2 for x, y in zip(xs, ys)) / (n - 2)
        vk, vb, cov = s2 / sxx, s2 * (1 / n + mx * mx / sxx), -mx * s2 / sxx
        dk, db = b / (k * k), -1 / k
        d["sigma_V_bd_V"] = math.sqrt(max(0.0, dk * dk * vk + db * db * vb + 2 * dk * db * cov))
    return d


# ---------------------------------------------------------------- puntos G(VMon) por canal
day_pts, night_pts = {ch: [] for ch in range(32)}, {ch: [] for ch in range(32)}
for s in DOFF["steps"]:
    for ch in range(32):
        g = s["G_new"].get(str(ch)) or s["G_new"].get(ch)
        h = hv(s, ch)
        if g and g.get("measured") and h:
            day_pts[ch].append((h["VMon_V"], g["G_mV"], s["delta_V"]))
for s in NSTEPS:
    pc = s["sw"]["per_channel"]
    for ch in range(32):
        p = pc.get(str(ch)) or pc.get(ch)
        h = hv(s, ch)
        if p and p.get("G_measured") and h:
            night_pts[ch].append((h["VMon_V"], p["G_mV"], s["delta_V"]))
VB = {}
for ch in range(32):
    both = [(v, g) for v, g, _ in day_pts[ch] + night_pts[ch]]
    VB[ch] = {"all": fit(both), "day": fit([(v, g) for v, g, _ in day_pts[ch]]), "night": fit([(v, g) for v, g, _ in night_pts[ch]]),
              "points_day": day_pts[ch], "points_night": night_pts[ch]}
    nom = next((v for v, g, dv in night_pts[ch] if dv == 0), next((v for v, g, dv in day_pts[ch] if dv == 0), None))
    VB[ch]["VMon_nominal_V"] = nom
    if VB[ch]["all"] and nom is not None:
        VB[ch]["overvoltage_nominal_V"] = nom - VB[ch]["all"]["V_bd_V"]
T_day = statistics.mean(temps(s) for s in DRES["steps"][:2] if math.isfinite(temps(s))) if DRES["steps"] else float("nan")
T_night = statistics.mean(t for t in (temps(s) for s in NSTEPS) if math.isfinite(t)) if NSTEPS else float("nan")

L += ["# DAQ2-17: barrido de tensión hacia arriba (noche) junto al de día", "",
      f"Noche: campaña `{os.path.basename(NIGHT)}`, escalones medidos {[s['delta_V'] for s in NSTEPS]}; omitidos {NRES.get('skipped')}; "
      f"omitidos por fallo recuperable {[k['delta_V'] for k in NRES.get('skipped_steps', [])]}; aborto: {NRES.get('aborted')}; "
      f"guarda de corriente: {NRES.get('current_guard', 'no se activó')}; restauración: {NRES.get('restore', {}).get('result')}; "
      f"run de física al final: {NRES.get('physics')}.",
      f"Día: `{os.path.basename(DAY)}` (G del recálculo offline DAQ2-16). Temperatura media de la placa A7030DP: día (0 y −0.5 V) "
      f"{f(T_day, 2)} °C, noche {f(T_night, 2)} °C (resolución 1 °C; sigue la del rack, no la del detector).", "",
      "## G contra VMon y V_bd por canal de hod2018 (V_bd estimada: recta G = k·(VMon − V_bd))", "",
      "| canal | puntos día / noche | k (mV/p.e./V) | V_bd (V) | σ(V_bd) | sobretensión al nominal (V) | V_bd día (2 puntos) | V_bd noche | noche − día (V) |",
      "|---|---|---|---|---|---|---|---|---|"]
diffs = []
for ch in range(32):
    v = VB[ch]
    a, d, n = v["all"], v["day"], v["night"]
    dd = (n["V_bd_V"] - d["V_bd_V"]) if (d and n) else None
    if dd is not None:
        diffs.append(dd)
    L.append(f"| {ch} | {len(v['points_day'])} / {len(v['points_night'])} | {f(a['k_mV_per_V']) if a else '–'} | {f(a['V_bd_V'], 2) if a else '–'} | "
             f"{f(a['sigma_V_bd_V'], 2) if a and a['sigma_V_bd_V'] else '–'} | {f(v.get('overvoltage_nominal_V'), 2)} | {f(d['V_bd_V'], 2) if d else '–'} | "
             f"{f(n['V_bd_V'], 2) if n else '–'} | {f(dd, 2)} |")
vbs = [VB[ch]["all"]["V_bd_V"] for ch in range(32) if VB[ch]["all"]]
ovs = [VB[ch]["overvoltage_nominal_V"] for ch in range(32) if "overvoltage_nominal_V" in VB[ch]]
L += ["", f"V_bd (todos los puntos) en {len(vbs)} canales: mediana {f(statistics.median(vbs), 2) if vbs else '–'} V; sobretensión al nominal mediana "
      f"{f(statistics.median(ovs), 2) if ovs else '–'} V. Noche − día: mediana {f(statistics.median(diffs), 3) if diffs else '–'} V en {len(diffs)} canales "
      f"(media {f(statistics.mean(diffs), 3) if diffs else '–'} ± {f(statistics.stdev(diffs) / math.sqrt(len(diffs)), 3) if len(diffs) > 1 else '–'} V); "
      f"diferencia de temperatura de la placa noche − día {f(T_night - T_day, 2)} °C. Con un coeficiente típico de V_bd de unas decenas de mV/°C "
      "(supuesto, no medido aquí) una diferencia de ~1 °C de la placa no se distingue; la placa no mide el detector.", ""]
for part in (0, 1):
    pads = []
    for ch in range(16 * part, 16 * part + 16):
        pads.append({"title": f"ch{ch}", "x": [p[0] for p in VB[ch]["points_day"] + VB[ch]["points_night"]],
                     "ys": [[p[1] for p in VB[ch]["points_day"] + VB[ch]["points_night"]]]})
    figs.append({"file": f"G_vs_VMon_{part}.png", "title": "hod2018: G contra VMon (día y noche)", "xlabel": "VMon (V)", "ylabel": "G (mV/p.e.)", "kind": "points",
                 "series": [{"name": f"ch{ch} {'día' if k == 0 else 'noche'}", "x": [p[0] for p in pts], "y": [p[1] for p in pts], "yerr": None}
                            for ch in range(16 * part, 16 * part + 16) for k, pts in ((0, VB[ch]["points_day"]), (1, VB[ch]["points_night"])) if pts]})

# ---------------------------------------------------------------- mayo contra octubre con el mismo V0Set (René, DAQ2-17)
G_MAY = {}
for l in open(REPO + "/outputs/pe_calibration_v8.txt"):
    if not l.startswith("#") and l.split() and int(l.split()[0]) < 32:
        G_MAY[int(l.split()[0])] = float(l.split()[3])
L += ["## Mayo contra octubre con el mismo V0Set", "",
      "El V0Set de mayo registrado (captura de GECO2020, handoff línea 21; no se sabe en cuáles de v7-v10 rigió) es igual al nominal de hoy. "
      "Con V0Set igual, la diferencia de G implica un corrimiento de la tensión de ruptura de ΔV_bd ≈ (G_mayo − G_oct) / k, con k la pendiente "
      "del ajuste de esta noche (G de mayo: `outputs/pe_calibration_v8.txt`, con el pedestal común G_ajustada/G_archivo = 0.9946; G de octubre: "
      "escalón ΔV = 0 de la noche). Estimado: supone que la cadena de lectura tiene la misma ganancia en mV por p.e. y que el V0Set de mayo "
      "era el registrado.", "",
      "| canal | G mayo (mV) | G oct. ΔV = 0 (mV) | k noche (mV/V) | ΔV_bd ≈ (G_mayo − G_oct)/k (V) |", "|---|---|---|---|---|"]
dvb = []
for ch in range(32):
    g0 = next((g for v, g, dv in night_pts[ch] if dv == 0), None)
    k = VB[ch]["night"]["k_mV_per_V"] if VB[ch]["night"] else None
    d = (G_MAY[ch] - g0) / k if (ch in G_MAY and g0 and k) else None
    if d is not None:
        dvb.append(d)
    L.append(f"| {ch} | {f(G_MAY.get(ch))} | {f(g0)} | {f(k)} | {f(d, 2)} |")
L += ["", f"Mediana de ΔV_bd ≈ {f(statistics.median(dvb), 2) if dvb else '–'} V en {len(dvb)} canales (G de octubre menor: V_bd más alta en octubre, o "
      "menos ganancia de la cadena). Causas posibles, que con los datos actuales no se distinguen: temperatura del detector (no medida; la placa "
      "A7030DP no la representa), ganancia de la cadena de lectura (cables, amplificadores, DC offset), y un cambio de V0Set no registrado "
      "durante v8 (René cambió V0Set a mano el 2026-05-29, handoff línea 64).", ""]

# ---------------------------------------------------------------- por escalón de la noche
def ov_med(s):
    o = []
    for ch in range(32):
        h = hv(s, ch)
        a = VB[ch]["all"]
        if h and a:
            o.append(h["VMon_V"] - a["V_bd_V"])
    return statistics.median(o) if o else float("nan")


L += ["## Por escalón de la noche", "",
      "| ΔV | sobretensión mediana (V) | IMon máx. por grupo CH00-07 (µA) | G mediana (mV, canales) | tasa oscura mediana > 0.5 / 1.5 / 3 / 5 / 8 p.e. (Hz) | diafonía mediana | placa (°C) |",
      "|---|---|---|---|---|---|---|"]
rows = []
for s in NSTEPS:
    pc = s["sw"]["per_channel"]
    get = lambda ch: pc.get(str(ch), pc.get(ch, {}))
    rate = {T: statistics.median([get(ch)[f"rate_{T}pe_Hz"] for ch in range(32) if ch != 3 and get(ch).get(f"rate_{T}pe_Hz") is not None] or [float("nan")])
            for T in ("0.5", "1.5", "3", "5", "8")}
    xt = [get(ch)["rate_1.5pe_Hz"] / get(ch)["rate_0.5pe_Hz"] for ch in range(32) if get(ch).get("rate_0.5pe_Hz")]
    im = s.get("imon_max_A", {})
    r = {"dv": s["delta_V"], "ov": ov_med(s), "imon": [im.get(str(c), im.get(c)) for c in range(8)], "G": s["sw"]["G_median_mV"], "nG": s["sw"]["n_G_ok"],
         "rate": rate, "xt": statistics.median(xt) if xt else float("nan"), "T": temps(s)}
    rows.append(r)
    L.append(f"| {r['dv']:+.2f} | {f(r['ov'], 2)} | " + ", ".join(f(x * 1e6 if x else None, 2) for x in r["imon"]) + f" | {f(r['G'])} ({r['nG']}) | "
             + " / ".join(f"{rate[T]:.4g}" for T in ("0.5", "1.5", "3", "5", "8")) + f" | {f(r['xt'])} | {f(r['T'], 1)} |")
L += ["", "Eficiencia (definición única de impacto, sin ch3; estimador compartido; fracción / 0.03646, eventos mezclados):", "",
      "| ΔV | eventos | ε 5 p.e. | ε 8 p.e. | ε 10 p.e. | ε 12 p.e. | ε(X\\|Y) / ε(Y\\|X) 5 p.e. | ε(X\\|Y) / ε(Y\\|X) 8 p.e. | ε(X\\|Y) / ε(Y\\|X) 10 p.e. |",
      "|---|---|---|---|---|---|---|---|---|"]
for s in NSTEPS:
    e = s["ext"]["all"]
    g = lambda T: e.get(str(T), e.get(T, {}))
    L.append(f"| {s['delta_V']:+.2f} | {s['ext_run']['events']} | " + " | ".join(f"{g(T)['f_mix'] / ACC:.3f} ± {g(T)['f_mix_err'] / ACC:.3f}" for T in (5, 8, 10, 12))
             + " | " + " | ".join(f"{g(T)['e_X|Y_mix']:.3f} / {g(T)['e_Y|X_mix']:.3f}" for T in (5, 8, 10)) + " |")

# espectro en tiempo del plano X con Y > 5 p.e. contra la ventana fuera de tiempo
L += ["", "Espectro en tiempo del plano X (máximo en p.e.) en eventos con Y sobre 5 p.e., en tiempo contra fuera de tiempo (bins de 1 p.e.):", "",
      "| ΔV | eventos con Y > 5 | X > 10 p.e. en tiempo / fuera | exceso sobre 10 p.e. | máximo del exceso (p.e.) | ¿pico separado del ruido? |", "|---|---|---|---|---|---|"]
spec_series = []
for i, s in enumerate(NSTEPS):
    hf = f"{NIGHT}/hits_step{i}_hits.bin"
    cal = json.load(open(s["sw"]["calibration"]))
    G = np.array([next((e["G_mV"] for e in cal["channels"] if e["ch"] == c and e.get("G_mV")), np.nan) for c in range(32)])
    G = np.where(np.isfinite(G), G, np.nanmedian(G))
    if not os.path.exists(hf):
        L.append(f"| {s['delta_V']:+.2f} | – | – | – | – | sin impactos |")
        continue
    ev = xy.load_hits(hf)
    xs_ = [c for c in range(16) if c != 3]
    ysel = (ev["aS"][:, 16:32] / G[16:32] > 5).any(1)
    xin = (ev["aS"][ysel][:, xs_] / G[xs_]).max(1)
    xct = (ev["aC"][ysel][:, xs_] / G[xs_]).max(1)
    hin, _ = np.histogram(xin, bins=60, range=(0, 60))
    hct, _ = np.histogram(xct, bins=60, range=(0, 60))
    exc = hin - hct
    a10, c10 = int((xin > 10).sum()), int((xct > 10).sum())
    k = int(np.argmax(exc[5:])) + 5
    sig = exc[k] / math.sqrt(max(hin[k] + hct[k], 1))
    sep = "sí" if (k >= 6 and sig > 3 and exc[k] > 2 * max(exc[3:6].min(), 1)) else "no"
    L.append(f"| {s['delta_V']:+.2f} | {int(ysel.sum())} | {a10} / {c10} | {a10 - c10} ± {math.sqrt(a10 + c10):.0f} | {k + 0.5:.1f} ({exc[k]} ± {math.sqrt(hin[k] + hct[k]):.0f}) | {sep} |")
    spec_series += [{"name": f"{s['delta_V']:+.2f} V en tiempo", "x": list(range(60)), "y": hin.tolist(), "yerr": None},
                    {"name": f"{s['delta_V']:+.2f} V fuera de tiempo", "x": list(range(60)), "y": hct.tolist(), "yerr": None}]
figs.append({"file": "spectrum_X_given_Y5.png", "title": "hod2018 X (máximo) con Y > 5 p.e.: en tiempo y fuera de tiempo, por escalón", "xlabel": "p.e.",
             "ylabel": "eventos", "kind": "lines", "logy": True, "series": spec_series})
L += ["", "Criterio de \"pico separado\": el máximo del exceso (en tiempo − fuera) sobre 5 p.e. está en ≥ 6 p.e., con > 3σ y más del doble que el "
      "exceso en 3-5 p.e.", ""]

# ---------------------------------------------------------------- hod2019
def h19(steps, src):
    out = {}
    for s in steps:
        pc = s.get("sw", {}).get("per_channel", {})
        ex = s.get("sw", {}).get("exposure_s")
        for ch in range(32, 64):
            p = pc.get(str(ch), pc.get(ch, {}))
            if p.get("rate_1mV_Hz") is not None:
                out.setdefault(ch, []).append((s["delta_V"], p["rate_1mV_Hz"], ex, src))
    return out


H = {}
for d in (h19(DRES["steps"], "día"), h19(NSTEPS, "noche")):
    for ch, v in d.items():
        H.setdefault(ch, []).extend(v)
L += ["## hod2019: tasa > 1 mV por canal contra ΔV (día y noche)", "",
      "| canal | " + " | ".join(f"{dv:+.2f} V {src}" for dv, _, _, src in sorted(H.get(32, []), key=lambda x: x[0])) + " | responde |",
      "|---|" + "---|" * (len(H.get(32, [])) + 1)]
resp = []
for ch in range(32, 64):
    v = sorted(H.get(ch, []), key=lambda x: x[0])
    cells = [f"{r:.3g}" for _, r, _, _ in v]
    lo = [x for x in v if x[0] <= -1.0]
    hi = [x for x in v if x[0] >= 1.0]
    z = None
    if lo and hi and lo[0][2] and hi[-1][2]:
        nlo, nhi = lo[0][1] * lo[0][2], hi[-1][1] * hi[-1][2]
        z = (nhi - nlo) / math.sqrt(max(nhi + nlo, 1))
    r = "sí" if z is not None and z > 5 else ("no" if z is not None else "–")
    if r == "sí":
        resp.append(ch)
    L.append(f"| {ch} | " + " | ".join(cells) + f" | {r}" + (f" (z = {z:.1f})" if z is not None else "") + " |")
L += ["", f"Criterio: cuentas a ΔV = +1.0 V contra ΔV ≤ −1.0 V con z > 5 (Poisson). Canales que responden: {resp}.", ""]

# ---------------------------------------------------------------- diagnóstico del trigger
diag = next((s.get("trigger_diag") for s in NSTEPS if s.get("trigger_diag")), None)
L += ["## Diagnóstico del trigger (ΔV = 0 de la noche, runs de 30 s)", ""]
if diag:
    L += ["| grupos | mayoría | ventana (relojes) | umbral (p.e.) | predicción (Hz) | eventos | tasa (Hz) | razón medida / prevista | registros releídos = escritos |",
          "|---|---|---|---|---|---|---|---|---|"]
    for d in diag:
        rg = d.get("registers") or {}
        ratio = d["rate_Hz"] / d["pred_Hz"] if d.get("rate_Hz") is not None and d.get("pred_Hz") else float("nan")
        L.append(f"| {'+'.join('G' + str(g) for g in d['groups'])} | {d['majority_groups']} | {d['window_clk']} | {d['T_pe']:g} | {d['pred_Hz']:.3g} | "
                 f"{d.get('events', '–')} | {f(d.get('rate_Hz'))} | {f(ratio, 3)} | "
                 f"{('sí' if not rg.get('mismatch') else 'NO: ' + str(rg.get('mismatch'))) if rg.get('read') else d.get('error', '–')} |")
else:
    L.append("No hay runs de diagnóstico en los resultados de la noche.")
L.append("")

# ---------------------------------------------------------------- candidatos
L += ["## Puntos de operación candidatos (sin elegir; los decide René)", "",
      "| ΔV | sobretensión mediana (V) | ε 5 p.e. | ε 8 p.e. | ε 10 p.e. | ruido mediano > 8 p.e. por canal (Hz) | IMon máx. (µA) |", "|---|---|---|---|---|---|---|"]
for s, r in zip(NSTEPS, rows):
    e = s["ext"]["all"]
    g = lambda T: e.get(str(T), e.get(T, {}))
    im = [x for x in r["imon"] if x]
    L.append(f"| {s['delta_V']:+.2f} | {f(r['ov'], 2)} | " + " | ".join(f"{g(T)['f_mix'] / ACC:.3f} ± {g(T)['f_mix_err'] / ACC:.3f}" for T in (5, 8, 10))
             + f" | {r['rate']['8']:.4g} | {f(max(im) * 1e6 if im else None, 2)} |")
open(OUT + "/report_es.md", "w").write("\n".join(L) + "\n")
json.dump({"V_bd": VB, "steps": rows, "hod2019": {str(k): v for k, v in H.items()}, "trigger_diag": diag, "T_day_degC": T_day, "T_night_degC": T_night},
          open(OUT + "/results.json", "w"), indent=1, default=str)
json.dump(figs, open(OUT + "/figs.json", "w"), default=float)
print("\n".join(L))
