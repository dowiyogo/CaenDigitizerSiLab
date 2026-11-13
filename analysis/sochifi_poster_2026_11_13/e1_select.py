#!/usr/bin/python3 -I
# e1_select.py — SOCHIFI 2026 (póster), revisión 2 (2026-11-13): selección de un muón real en hod2018 para la figura E1 (event display).
# Solo lectura de salidas existentes del proceeding; mismas definiciones que analysis/sochifi_2026_11_12/impact_map.py:
#  - impactos: work/hd64_may_v8_hits.bin (hitdef_extract.C; aS = pulso más alto de la ventana signal, aC = de la ventana fuera de tiempo);
#  - G por canal: reajuste de PROC (work/hd64_may_v8_peaks.csv, estado ok/changed; mediana para el canal sin G), umbral = T·G;
#  - orden de tiras de mayo: work/stage3.json strip_orders["may"]; ch3 (X, posición 12) muerto y fuera;
#  - muestra: estacionario (t ≥ 5.25 h desde el primer disparo) ∩ tiempo bueno de v8 (may_common.good_time): 76 784 disparos.
# Selección (pedido de la tarea, sin umbrales nuevos: 10 p.e. y 5 p.e. son los del análisis):
#  (a) X·Y en tiempo a 10 p.e. por canal;
#  (b) en cada eje, exactamente una tira o dos tiras vecinas (en el orden de mayo) sobre 10 p.e.;
#  (c) ningún otro canal de hod2018 sobre 5 p.e. en la ventana signal;
#  (d) ningún canal de hod2018 sobre 5 p.e. en la ventana fuera de tiempo;
#  (e) tiras encendidas no vecinas de ch3 (posiciones X 11, 12, 13 excluidas).
# Entre los candidatos con una celda de dos tiras en al menos un eje se elige uno a mano con una regla fija (ver elección abajo).
# Los registros de hits.bin están en el orden del archivo crudo (143 814 escritos de 143 814, sin pasos atrás en el tiempo), así que el índice
# en hits.bin = entrada del árbol data_dig; e1_event.C lo comprueba con el timestamp.
# Uso: e1_select.py <out dir>   -> e1_selection.json, e1_candidates.csv
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
PK = os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work")
ORD = json.load(open(PK + "/stage3.json"))["strip_orders"]["may"]["order"]
XO, YO = np.asarray(ORD["h18_X"]), np.asarray(ORD["h18_Y"])
DEAD, T_STAT_H, T_HIT, T_VETO = 3, 5.25, 10, 5
NOISY = {2, 9, 15}
BOOT_REF = 1778886208.8  # hitdef_extract.C: t = BOOT_REF + timestamp·1e-6


def g_refit():  # idéntica a impact_map.py
    g = np.full(32, np.nan)
    for l in open(PK + "/hd64_may_v8_peaks.csv"):
        f = l.strip().split(",")
        if f[0] != "ch" and f[3] in ("ok", "changed"):
            g[int(f[0])] = float(f[2])
    return np.where(np.isfinite(g), g, np.nanmedian(g))


G = g_refit()
raw = np.fromfile(PK + "/hd64_may_v8_hits.bin", dtype=xy.hits_dtype(64))
assert len(raw) == 143814 and (np.diff(raw["t"]) >= 0).all(), "hits.bin no está en orden de tiempo"
e_all = xy.load_hits(PK + "/hd64_may_v8_hits.bin", 64)
assert (e_all["t"] == raw["t"]).all()
t0 = float(e_all["t"].min())
msk = (e_all["t"] >= t0 + T_STAT_H * 3600) & mc.good_time(e_all["t"], "v8")
idx = np.nonzero(msk)[0]
e = e_all[msk]
N = len(e)
aS, aC = e["aS"][:, :32], e["aC"][:, :32]
live = np.ones(32, bool)
live[DEAD] = False
pe_S, pe_C = aS / G, aC / G

# (a) X·Y en tiempo a 10 p.e. (misma definición que xy_estimator.plane_hits)
XCH = [c for c in range(16) if c != DEAD]
YCH = list(range(16, 32))
GG = np.concatenate([G, np.full(32, np.nan)])
Xs, _ = xy.plane_hits(e, GG, T_HIT, XCH)
Ys, _ = xy.plane_hits(e, GG, T_HIT, YCH)
cut = {"triggers_stationary_goodtime": int(N)}
sa = Xs & Ys
cut["a_XY_on_time_10pe"] = int(sa.sum())


def topo(H):
    """H (n, 16) en orden de posición -> (tipo, primera posición): 1 = una tira, 2 = dos vecinas, 0 = otro patrón."""
    n = H.sum(1)
    first = np.argmax(H, 1)
    nxt = H[np.arange(len(H)), np.minimum(first + 1, 15)]
    t = np.zeros(len(H), int)
    t[n == 1] = 1
    t[(n == 2) & (first < 15) & nxt] = 2
    return t, first


HX = pe_S[:, XO] > T_HIT
HX[:, XO == DEAD] = False
HY = pe_S[:, YO] > T_HIT
tX, fX = topo(HX)
tY, fY = topo(HY)
sb = sa & (tX > 0) & (tY > 0)
cut["b_one_or_two_adjacent_each_axis"] = int(sb.sum())
fired = np.zeros((N, 32), bool)
fired[:, XO] = HX
fired[:, YO] = HY
other5 = ((pe_S > T_VETO) & ~fired & live).any(1)
sc = sb & ~other5
cut["c_no_other_channel_gt5pe_signal"] = int(sc.sum())
off5 = ((pe_C > T_VETO) & live).any(1)
sd = sc & ~off5
cut["d_no_channel_gt5pe_offtime"] = int(sd.sum())
posX12 = int(np.nonzero(XO == DEAD)[0][0])
near = HX[:, [p for p in (posX12 - 1, posX12, posX12 + 1) if 0 <= p < 16]].any(1)
se = sd & ~near
cut["e_not_adjacent_to_ch3"] = int(se.sum())
two_any = se & ((tX == 2) | (tY == 2))
two_both = se & (tX == 2) & (tY == 2)
cut["candidates_two_bar_cell_in_at_least_one_axis"] = int(two_any.sum())
cut["candidates_two_bar_cell_in_both_axes"] = int(two_both.sum())
cut["candidates_single_bar_both_axes"] = int((se & (tX == 1) & (tY == 1)).sum())

rows = []
for k in np.nonzero(se)[0]:
    chans = [int(c) for c in np.nonzero(fired[k])[0]]
    rows.append({"entry": int(idx[k]), "t_utc_s": float(e["t"][k]), "type_x": int(tX[k]), "type_y": int(tY[k]),
                 "cell_x": int(2 * fX[k] + (1 if tX[k] == 2 else 0)), "cell_y": int(2 * fY[k] + (1 if tY[k] == 2 else 0)),
                 "channels": chans, "pe": [round(float(pe_S[k, c]), 2) for c in chans], "sum_pe": float(pe_S[k, chans].sum()),
                 "noisy_fired": any(c in NOISY for c in chans), "max_other_pe_signal": float(np.max(np.where(fired[k] | ~live, 0, pe_S[k]))),
                 "max_pe_offtime": float(np.max(np.where(live, pe_C[k], 0)))})
with open(os.path.join(OUT, "e1_candidates.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["entry", "t_utc_s", "type_x", "type_y", "cell_x", "cell_y", "channels", "pe", "sum_pe", "noisy_fired", "max_other_pe_signal", "max_pe_offtime"])
    for r in rows:
        w.writerow([r["entry"], "%.6f" % r["t_utc_s"], r["type_x"], r["type_y"], r["cell_x"], r["cell_y"], " ".join(map(str, r["channels"])),
                    " ".join("%.2f" % p for p in r["pe"]), "%.2f" % r["sum_pe"], int(r["noisy_fired"]), "%.2f" % r["max_other_pe_signal"],
                    "%.2f" % r["max_pe_offtime"]])

# Elección (regla fija, ejemplo escogido a mano): celda de dos tiras en los dos ejes si existe (si no, en al menos uno), sin canales ruidosos
# (2, 9, 15) entre las tiras encendidas, cruce lejos del borde (ninguna tira encendida en las posiciones exteriores 0 y 15, para que la vista
# muestre un cruce dentro del área), y entre esos la suma de p.e. más cercana a la mediana del grupo (un evento típico, no el mayor).
# Nota: la regla del borde se agregó después de ver que la primera elección (sin ella) era la entrada 23190, en la esquina (celdas 29, 29).
pool = [r for r in rows if r["type_x"] == 2 and r["type_y"] == 2] or [r for r in rows if r["type_x"] == 2 or r["type_y"] == 2]
rule = "dos tiras en ambos ejes" if pool and pool[0]["type_x"] == 2 and pool[0]["type_y"] == 2 else "dos tiras en al menos un eje"
quiet = [r for r in pool if not r["noisy_fired"]] or pool
inner = [r for r in quiet if all(int(np.nonzero((XO if c < 16 else YO) == c)[0][0]) not in (0, 15) for c in r["channels"])]
cut["candidates_pool_both_axes_quiet"] = len(quiet)
cut["candidates_pool_both_axes_quiet_inner"] = len(inner)
quiet = inner or quiet
med = float(np.median([r["sum_pe"] for r in quiet]))
ch = min(quiet, key=lambda r: (abs(r["sum_pe"] - med), r["entry"]))
pos = {}
for p, c in enumerate(XO):
    pos[int(c)] = ("X", p)
for p, c in enumerate(YO):
    pos[int(c)] = ("Y", p)
ch["bars"] = [{"channel": c, "bar": "%s%d" % pos[c], "axis": pos[c][0], "position": pos[c][1], "center_mm": -15 + 2 * pos[c][1], "pe": p}
              for c, p in zip(ch["channels"], ch["pe"])]
ch["timestamp_us"] = int(round((ch["t_utc_s"] - BOOT_REF) * 1e6))
ch["t_since_first_trigger_h"] = (ch["t_utc_s"] - t0) / 3600
ch["choice_rule"] = rule + "; sin canales ruidosos 2, 9, 15; sin tiras en las posiciones 0 y 15; suma de p.e. más cercana a la mediana del grupo (%.2f p.e., %d eventos)" % (med, len(quiet))
out = {"run": "v8 (run B, mayo 2026)", "raw_file": "/data/daq/raw/hodoscope_extnim_2026/lectura_vx1740_ext_nim_1M_v8.root",
       "sample": "estacionario (t >= 5.25 h) y tiempo bueno de v8", "thresholds_pe": {"hit": T_HIT, "veto": T_VETO},
       "G_per_channel_mV": [round(float(g), 4) for g in G], "strip_order_may": {"X": XO.tolist(), "Y": YO.tolist()},
       "selection_counts": cut, "chosen": ch}
json.dump(out, open(os.path.join(OUT, "e1_selection.json"), "w"), indent=1, ensure_ascii=False)
for k, v in cut.items():
    print(f"{k}: {v}")
print("elegido:", json.dumps({k: ch[k] for k in ("entry", "timestamp_us", "t_utc_s", "cell_x", "cell_y", "bars", "choice_rule")}, ensure_ascii=False))
