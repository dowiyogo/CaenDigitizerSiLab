#!/usr/bin/python3 -I
# vbd_plan.py — seguimiento 2026-10-31 (9), punto 1: plan del barrido de V_bd (preparado, no armado). Corriente esperada por grupo HV en cada escalón
# con el modelo I-V (b) del seguimiento 8 (I = I_0 + a (V − V_bd)^n, iv_fit.json) y margen hasta I0Set = 15 µA; regla: un grupo con margen < 2 µA a
# +0.15 V va a +0.10 V.
# Modelo por grupo: el ajuste con V_bd libre si no quedó en el borde de la rejilla; si no, el de V_bd fijo (grupos con estimación de dos puntos);
# si tampoco hay, cota por el trip del 2026-10-30 (solo el grupo 4, CH03): forma (b) con I_0 = piso medido, V_bd = mediana de los V_bd libres
# determinados [supuesto] y el menor n compatible con I(VMon 65.80 V) ≥ 15.45 µA (corriente limitada antes del trip): da una COTA INFERIOR de la
# corriente sobre el nominal, así que el margen calculado es una cota superior.
# Escalón ΔV de V0Set = mismo ΔV de VMon desde el VMon nominal del 2026-10-23 [supuesto: desvío VMon − V0Set constante en ±0.4 V].
# Guarda de corriente de biasscan.py: > 12 µA 30 s en algún grupo -> no pasa al escalón siguiente; ≥ 14 µA -> aborto.
# Eventos por escalón: ver vbd_events.C y el informe del seguimiento 9 (el ajuste de c1 no converge a −0.5 V con ningún número de eventos; el
# período sí con 12 000). Duración de una toma SW: eventos / tasa SW (100 Hz) + 10 s; escalón = estabilización + toma + escrituras.
# Salida: JSON y tabla en Markdown a stdout. Uso: vbd_plan.py <iv_fit.json> [eventos por escalón: "0.15:12000,0:12000,-0.25:24000,-0.4:24000"]
import json, math, statistics as st, sys

IV = json.load(open(sys.argv[1]))
EV = {float(k): int(v) for k, v in (x.split(":") for x in (sys.argv[2] if len(sys.argv) > 2 else "0.15:12000,0:12000,-0.25:24000,-0.4:24000").split(","))}
STEPS = [0.15, 0.0, -0.25, -0.4]
I0SET, MARGIN, WARN, STOP = 15.0, 2.0, 12.0, 14.0
SW_HZ, STAB_S, WRITE_S = 100.0, 120.0, 60.0


def model(m, i0, v):
    x = max(0.0, v - m["V_bd_V"])
    return i0 + m["a_uA"] * x ** m["n"]


free_vbd = [g["b_free"]["V_bd_V"] for g in IV["groups"].values() if g.get("b_free") and not g["b_free"]["at_grid_edge"]]
vbd_med = st.median(free_vbd)
plan = {"I0Set_uA": I0SET, "margin_rule_uA": MARGIN, "guards_uA": {"warn_30s": WARN, "stop": STOP}, "V_bd_median_free_V": vbd_med, "groups": {}}
for g, d in sorted(IV["groups"].items(), key=lambda kv: int(kv[0])):
    i0, vnom = d["I0_floor_uA"], d["V_nom_VMon_V"]
    if d.get("b_free") and not d["b_free"]["at_grid_edge"]:
        m, src, bound = d["b_free"], "b_free", False
    elif d.get("b_fixed"):
        m, src, bound = d["b_fixed"], "b_fixed (V_bd de dos puntos)", False
    else:
        # cota por el trip: I(65.80) ≥ 15.45 µA con V_bd = mediana y el nivel nominal medido
        nom = max((p for p in d["points"] if p["scan"] == "scan23"), key=lambda p: p["VMon_V"])
        xr = (65.80 - vbd_med) / (nom["VMon_V"] - vbd_med)
        n = math.log((15.45 - i0) / (nom["IMon_uA"] - i0)) / math.log(xr)
        a = (nom["IMon_uA"] - i0) / (nom["VMon_V"] - vbd_med) ** n
        m, src, bound = {"V_bd_V": vbd_med, "a_uA": a, "n": n}, f"cota por el trip (V_bd {vbd_med:.3f} V [S], n ≥ {n:.2f})", True
    row = {"model": src, "lower_bound": bound, "V_nom_VMon_V": vnom, "I_uA": {}, "margin_uA": {}}
    for dv in STEPS:
        i = model(m, i0, vnom + dv)
        row["I_uA"][str(dv)] = i
        row["margin_uA"][str(dv)] = I0SET - i
    m15 = row["margin_uA"]["0.15"]
    row["step_plus_V"] = 0.10 if (m15 < MARGIN or bound) else 0.15
    row["reason"] = ("margen < 2 µA a +0.15 V" if m15 < MARGIN else ("margen sin cota inferior (solo cota de la corriente por el trip de CH03 del 2026-10-30): "
                                                                     "+0.10 V por precaución" if bound else "margen ≥ 2 µA"))
    if row["step_plus_V"] == 0.10:
        i = model(m, i0, vnom + 0.10)
        row["I_uA"]["0.1"] = i
        row["margin_uA"]["0.1"] = I0SET - i
    row["guard_warn_risk"] = any(v > WARN for k, v in row["I_uA"].items() if float(k) == row["step_plus_V"] or float(k) <= 0)
    plan["groups"][g] = row
tot = 0.0
plan["steps"] = []
for dv in STEPS:
    n = EV[dv]
    sw = n / SW_HZ + 10
    tot += STAB_S + sw + WRITE_S
    plan["steps"].append({"delta_V": dv, "events": n, "sw_s": sw, "step_s": STAB_S + sw + WRITE_S})
plan["total_s_without_restore"] = tot
json.dump(plan, open("/dev/stdout", "w"), indent=1) if "--json" in sys.argv else None
print("| Grupo | Modelo | I nominal | I +0.15 V | I +0.10 V | I −0.25 V | I −0.4 V | Margen a +0.15 V | Escalón positivo | Motivo |")
print("|---|---|---|---|---|---|---|---|---|---|")
for g, r in plan["groups"].items():
    I = r["I_uA"]
    f = lambda k: f"{I[k]:.2f}" + (" (cota inf.)" if r["lower_bound"] and float(k) > 0 else "") if k in I else "—"
    print(f"| {g} | {r['model']} | {I['0.0']:.2f} | {f('0.15')} | {f('0.1')} | {I['-0.25']:.2f} | {I['-0.4']:.2f} | {r['margin_uA']['0.15']:.2f}"
          f"{' (cota sup.)' if r['lower_bound'] else ''} | +{r['step_plus_V']:.2f} V | {r['reason']} |")
print()
for s in plan["steps"]:
    print(f"escalón {s['delta_V']:+.2f} V: {s['events']} eventos SW, toma {s['sw_s']:.0f} s, escalón {s['step_s']:.0f} s")
print(f"total sin restauración ni arranque: {tot:.0f} s ({tot / 60:.1f} min); V_bd mediana de los ajustes libres determinados {vbd_med:.3f} V")
print("grupos con riesgo de la guarda de 12 µA en el escalón positivo:", [g for g, r in plan["groups"].items() if r["guard_warn_risk"]])
json.dump(plan, open(sys.argv[3], "w"), indent=1) if len(sys.argv) > 3 else None
