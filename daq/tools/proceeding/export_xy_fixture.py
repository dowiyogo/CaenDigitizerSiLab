#!/usr/bin/python3 -I
# export_xy_fixture.py — XPLUS-3 (solo lectura de work/): fixture para la prueba de humo de daq_analysis (SPLIT-1, etapa 2.7).
# xy_estimator.xy_fraction() no toma cuentas sino cuatro arreglos booleanos por evento en orden de tiempo (Xs, Ys: impacto en la
# ventana signal; Xc, Yc: fuera de tiempo; los mezclados combinan el evento i con el i+1). Este script los arma igual que
# stage2.py (load + frac: work/hd64_<run>_hits.bin de PROC-1 etapa 1, G de mayo reajustada o c1, tiempo bueno, XCH sin ch3,
# YCH 16-31, T·G), comprueba que xy_fraction() da exactamente los valores guardados en work/stage2.json y los escribe
# empaquetados (np.packbits, base64), con las cuentas suficientes, la salida esperada y la procedencia. No re-extrae nada.
# Uso: export_xy_fixture.py <dir del paquete> <dataset de stage2 (oct8 | may8_good)> <T> <salida.json>
import base64, hashlib, json, os, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT, KEY, T, DEST = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
W = OUT + "/work"
# --- copia literal de stage2.py (DS, XCH, YCH, g_refit, fill, load, frac) ---
XCH = [c for c in range(16) if c != 3]
YCH = list(range(16, 32))
DS = {"may8_good": ("v8", "mayo v8, tiempo bueno", True), "oct8": ("oct8", "octubre run 8", False)}


def g_refit():
    g = np.full(32, np.nan)
    for l in open(W + "/hd64_may_v8_peaks.csv"):
        if l.startswith("ch,"):
            continue
        f = l.strip().split(",")
        if f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return g


G_MAY, G_C1 = g_refit(), mc.g_c1()
fill = lambda g: np.where(np.isfinite(g), g, np.nanmedian(g))
run, label, good = DS[KEY]
hits = f"{W}/hd64_may_{run}_hits.bin" if run.startswith("v") else f"{W}/hd64_{run}_hits.bin"
e = xy.load_hits(hits, 64)
if good:
    e = e[mc.good_time(e["t"], run)]
G = np.concatenate([fill(G_MAY) if run.startswith("v") else fill(G_C1), np.full(32, np.nan)])
Xs, Xc = xy.plane_hits(e, G, T, XCH)
Ys, Yc = xy.plane_hits(e, G, T, YCH)
# --- fin de la copia ---
res = xy.xy_fraction(Xs, Ys, Xc, Yc)
stored = json.load(open(W + "/stage2.json"))
ref = stored["datasets"][KEY]["per_T"][str(T)]
cmp_keys = [k for k in ("N", "N_XY", "acc_mix", "acc_ctl", "aX", "aY", "rX", "rY", "f_mix", "f_mix_err", "f_ctl", "f_ctl_err",
                        "e_X|Y_mix", "e_X|Y_mix_err", "e_Y|X_mix", "e_Y|X_mix_err", "e_X|Y_ctl", "e_Y|X_ctl") if k in ref]
diff = {k: (res[k], ref[k]) for k in cmp_keys if res[k] != ref[k]}
if diff:
    sys.exit(f"no coincide con work/stage2.json: {diff}")


def pack(a):
    return base64.b64encode(np.packbits(np.asarray(a, bool)).tobytes()).decode()


def unpack(s, n):
    return np.unpackbits(np.frombuffer(base64.b64decode(s), np.uint8))[:n].astype(bool)


N = len(Xs)
arrays = {"Xs": pack(Xs), "Ys": pack(Ys), "Xc": pack(Xc), "Yc": pack(Yc)}
back = xy.xy_fraction(*(unpack(arrays[k], N) for k in ("Xs", "Ys", "Xc", "Yc")))
if any(back[k] != res[k] for k in res):
    sys.exit("la ida y vuelta del empaquetado no reproduce xy_fraction")
counts = {"N": N, "N_XY": int((Xs & Ys).sum()), "N_X": int(Xs.sum()), "N_Y": int(Ys.sum()), "N_Xc": int(Xc.sum()), "N_Yc": int(Yc.sum()),
          "XsYc": int((Xs & Yc).sum()), "XcYs": int((Xc & Ys).sum()), "XcYc": int((Xc & Yc).sum()),
          "mix_Xs_i_Ys_i+1": int((Xs[:-1] & Ys[1:]).sum()), "mix_Ys_i_Xs_i+1": int((Ys[:-1] & Xs[1:]).sum())}
sha = hashlib.sha256(open(hits, "rb").read()).hexdigest()
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
fx = {
    "fixture": f"xy_estimator {KEY} T={T}", "dataset": KEY, "label": label, "T_pe": T, "status": "measured (q de PROC-1 etapa 2)",
    "function": "xy_estimator.xy_fraction(Xs, Ys, Xc, Yc)",
    "how_to_decode": "np.unpackbits(np.frombuffer(base64.b64decode(s), np.uint8))[:N].astype(bool), para Xs, Ys, Xc, Yc; orden de tiempo",
    "inputs": {"N": N, "arrays_packbits_base64": arrays},
    "selection": {"hits_file": os.path.relpath(hits, OUT), "hits_sha256": sha, "hits_bytes": os.path.getsize(hits), "events_in_file": None,
                  "good_time": good, "x_channels": XCH, "y_channels": YCH, "threshold": "amplitud > T·G por canal (aS en tiempo, aC fuera de tiempo)",
                  "G_mV_used_ch0_31": [None if not np.isfinite(g) else float(g) for g in G[:32]],
                  "G_source": "G de mayo reajustada (work/hd64_may_v8_peaks.csv; canales sin G = mediana)" if run.startswith("v")
                  else "calibración c1 (daq/config/calibration/pe_calibration_2026_10_c1.json; canales sin G = mediana)"},
    "sufficient_counts": counts,
    "expected": {k: res[k] for k in res},
    "expected_check": f"xy_fraction() sobre estos arreglos = work/stage2.json datasets.{KEY}.per_T.{T} en {len(cmp_keys)} claves, igualdad exacta",
    "provenance": {"hits_extraction": "PROC-1 etapa 1, hitdef_extract.C (nch 64), commit 8027437 (results.json provenance, stage 1)",
                   "q": f"PROC-1 etapa 2, stage2.py, commit {stored['commit']} (work/stage2.json, generado {stored['generated_utc']}); "
                        "results.json v3 lo cita con commit e70e12b",
                   "exported_by": f"daq/tools/proceeding/export_xy_fixture.py, commit {git}, {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}"},
}
fx["selection"]["events_in_file"] = int(len(xy.load_hits(hits, 64))) if good else N
os.makedirs(os.path.dirname(DEST), exist_ok=True)
json.dump(fx, open(DEST, "w"), indent=1, ensure_ascii=False)
print(f"{KEY} T={T}: N {N}, N_XY {counts['N_XY']}, q_mix {res['f_mix']:.5f} ± {res['f_mix_err']:.5f}, q_ctl {res['f_ctl']:.5f} ± {res['f_ctl_err']:.5f}; "
      f"igual a stage2.json en {len(cmp_keys)} claves; {os.path.getsize(DEST)} B -> {DEST}")
