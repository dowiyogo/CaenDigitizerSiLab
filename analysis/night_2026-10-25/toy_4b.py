#!/usr/bin/python3 -I
# toy_4b.py — XPLUS-4b etapa 4b.4: propaga al estimador las razones medidas a^oot/a (etapa 4b.2/4b.3) con δ por coordenada.
# Mismo modelo que toy_estimator.py (importa sus funciones params/replica; el generador no cambia) salvo que a^oot = (1+δ_X)·a_X y
# (1+δ_Y)·a_Y por separado. Llama a xy_estimator.xy_fraction. Uso: toy_4b.py <deltas.json> <salida.json> <réplicas> <semilla>
# deltas.json: {"<nombre>": [δ_X, δ_Y] o [δ_X, δ_Y, "v8"|"run8"], ...} (conjunto por defecto: v8)
import json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv_saved = sys.argv
DELTAS, OUT, NREP, SEED = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
# importar el generador de toy_estimator.py sin ejecutar su programa principal
src = open(os.path.join(HERE, "toy_estimator.py")).read()
main_at = src.index("cfs = configs() if GROUP")
ns = {"__file__": os.path.join(HERE, "toy_estimator.py")}
sys.argv = ["toy_estimator.py", "/dev/null", "main", "1", "0"]
exec(compile(src[:main_at], "toy_estimator.py", "exec"), ns)
sys.argv = sys.argv_saved
FG = ns["F18"] * ns["G"]
rng = np.random.default_rng(SEED)
res = {}
for name, spec in json.load(open(DELTAS)).items():
    dX, dY = spec[0], spec[1]
    ds = spec[2] if len(spec) > 2 else "v8"
    cf = dict(ds=ds, q=0.0275, name=name, asym=False, a_scale=1.0, k=1.0, delta=0.0, zero=False)
    p = ns["params"](cf)
    p["aXo"], p["aYo"] = p["aX"] * (1 + dX), p["aY"] * (1 + dY)
    q = {"mix": [], "ctl": []}
    for _ in range(NREP):
        r = ns["replica"](rng, p)
        q["mix"].append(r["f_mix"])
        q["ctl"].append(r["f_ctl"])
    row = {"dataset": ds, "delta_X": dX, "delta_Y": dY, "q_T": p["q_T"]}
    for t in ("mix", "ctl"):
        a = np.array(q[t])
        b, eb = a.mean() - p["q_T"], a.std(ddof=1) / math.sqrt(NREP)
        row[t] = {"bias_q": b, "bias_q_err": eb, "bias_eps": b / FG, "bias_eps_err": eb / FG}
    res[name] = row
    print(f"{name}: δX {dX:+.4f} δY {dY:+.4f}  sesgo mix {row['mix']['bias_q']:+.2e} ± {row['mix']['bias_q_err']:.1e} (ε {row['mix']['bias_eps']:+.4f} ± {row['mix']['bias_eps_err']:.4f})  "
          f"ctl ε {row['ctl']['bias_eps']:+.4f}", flush=True)
json.dump({"replicas": NREP, "seed": SEED, "results": res}, open(OUT, "w"), indent=1)
