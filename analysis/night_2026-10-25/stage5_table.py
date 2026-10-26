#!/usr/bin/python3 -I
# stage5_table.py — XPLUS-4b etapa 5: tabla maestra por run (markdown) desde stage5.json. Uso: stage5_table.py <stage5.json> <salida.md>
import datetime, json, sys

d = json.load(open(sys.argv[1]))
f = lambda v, n=4: "—" if v is None else f"{v:.{n}f}"
utc = lambda t: datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime("%m-%d %H:%M")
L = ["| run | inicio UTC | eventos (bueno/total) | horas buenas | G (fuente, canales, mediana mV) | orden X / Y | aX, aY (5 p.e.) | q 5 p.e. mix ± stat (ctl) | q 10 p.e. | ε_XY 5 p.e. | ε_XY 10 p.e. | ocup. X / Y (5 p.e.) | f_2 X / Y (5 p.e.) |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in d["runs"]:
    a, b = r["T5"], r["T10"]
    L.append(f"| {r['run']} | {utc(r['start_utc'])} | {r['events_good']}/{r['events_all']} ({r['good_fraction']:.3f}) | {r['duration_good_h']:.1f} | "
             f"{'run' if r['G_source'] == 'run' else 'v8'}, {r['G_channels_accepted']}, {r['G_median_mV']:.3f} | {r['strip_order']['X']} / {r['strip_order']['Y']} "
             f"{'✔' if r['strip_order_valid'] else '✘'} | {f(a['aX'], 3)}, {f(a['aY'], 3)} | {f(a['q_mix'])} ± {f(a['q_mix_err'])} ({f(a['q_ctl'])}) | "
             f"{f(b['q_mix'])} ± {f(b['q_mix_err'])} | {f(a['eps_XY'], 3)} ± {f(a['eps_XY_err'], 3)} | {f(b['eps_XY'], 3)} ± {f(b['eps_XY_err'], 3)} | "
             f"{f(a.get('occ_X'), 2)} / {f(a.get('occ_Y'), 2)} | {f(a.get('f2_X'), 3)} / {f(a.get('f2_Y'), 3)} |")
open(sys.argv[2], "w").write("\n".join(L) + "\n")
print("\n".join(L))
