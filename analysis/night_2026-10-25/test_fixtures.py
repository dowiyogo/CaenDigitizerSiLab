#!/usr/bin/python3 -I
# test_fixtures.py — XPLUS-4 etapa 0.2: regresión del estimador X·Y (xy_estimator.xy_fraction) contra los fixtures de XPLUS-3.
# Decodifica los cuatro arreglos de cada fixture, llama a xy_fraction y exige igualdad exacta con "expected" y con el q de referencia.
import base64, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import xy_estimator as xy  # noqa: E402

REF = {"xy_estimator_run8_T5.json": (0.02386720834654745, 0.0016840614637531606),
       "xy_estimator_may8good_T5.json": (0.027495083052747792, 0.0022618261790133573)}
ok = True
for name, (q, eq) in REF.items():
    fx = json.load(open(os.path.join(REPO, "docs/split/fixtures", name)))
    n = fx["inputs"]["N"]
    a = fx["inputs"]["arrays_packbits_base64"]
    u = lambda s: np.unpackbits(np.frombuffer(base64.b64decode(s), np.uint8))[:n].astype(bool)
    r = xy.xy_fraction(u(a["Xs"]), u(a["Ys"]), u(a["Xc"]), u(a["Yc"]))
    bad = [k for k in fx["expected"] if r[k] != fx["expected"][k]]
    good = not bad and r["f_mix"] == q and r["f_mix_err"] == eq
    ok &= good
    print(f"{name}: N {n} N_XY {r['N_XY']} q {r['f_mix']!r} ± {r['f_mix_err']!r}  ref {q!r} ± {eq!r}  claves distintas {bad}  {'OK' if good else 'FALLA'}")
sys.exit(0 if ok else 1)
