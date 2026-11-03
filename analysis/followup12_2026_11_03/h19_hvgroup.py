#!/usr/bin/python3 -I
# h19_hvgroup.py — seguimiento 2026-11-03 (12), punto A: ¿los pulsos de hod2019 en tiempo van con pulsos de hod2018 del MISMO grupo HV?
# Datos: <work>/{may8,oct8}_ev.bin de ev_dump.C (aS, kS, aC de los 64 canales, = hits de PROC). Mayo: run B (v8) en tiempo bueno
# (may_common.good_time); octubre: run 8 completo (posterior a la pérdida de hod2019 del 2026-10-21: control sin señal de hod2019).
# Definiciones:
#  - impacto de hod2019 en tiempo en el canal c (32-63): aS_c > T, T = 35 mV y 3 mV; impacto de hod2018 en j (0-31 sin ch3): aS_j > 5·G_j con G_j de
#    PROC (G_fit de hd64_*_peaks.csv; sin ajuste aceptado, la semilla). Control fuera de tiempo: lo mismo con aC en los dos.
#  - mapas de grupo HV (grupo g = SY5527 CH0(g−1)):
#      "sheet"  (central, el de la tarea): hv_group_sheet del mapeo r4 = grupo que la planilla «Voltajes Hodoscopio» da a la etiqueta del canal
#               (hereda el supuesto etiqueta ↔ canal del handoff, sección 5);
#      "oct_measured": hv_group del mapeo r4 (prueba de grupos de DAQ2-10, 2026-10-21; canales sin medir fuera) — es el de octubre; para mayo
#               solo vale si los cables de HV y de señal no se movieron (supuesto; el plano Y de hod2018 se recableó);
#      "may_measured_g12": grupos 1 y 2 medidos en mayo con datos (runs/good_intervals_2026.csv: HV_OFF del grupo 2 y el intervalo HV-DIAG-5 de
#               los grupos 1 y 2); solo canales de hod2019 de esos grupos; «otros grupos» = los demás canales de hod2018.
#  - por canal c de hod2019 con grupo g: S_g = canales de hod2018 del grupo g, O_g = canales de hod2018 con grupo conocido ≠ g;
#    n_same, n_other = número de impactos de hod2018 en S_g, O_g en el trigger. Triggers A_c = con impacto en c; triggers B = sin ningún impacto de
#    hod2019 (sobre T). Observado = Σ_c Σ_{A_c} n; esperado = Σ_c N(A_c)·⟨n⟩_B(c). R_same = obs/esp en S, R_other en O; doble razón DR = R_same/R_other
#    (> 1: asociación con el mismo grupo más allá de la asociación general, p. ej. muones que cruzan las dos estaciones). También por grupo.
#    Errores: bootstrap de Poisson sobre triggers (NBOOT réplicas).
#  - amplitud del mayor pulso de hod2018 (max_j aS_j/G_j, p.e. medidos): en triggers con algún impacto de hod2019 contra todos los triggers; y el
#    mayor pulso dentro de S_g (pares trigger-canal de A_c) contra el mismo en B; bins de 0.5 p.e. hasta 100 p.e.
#  - Δk = kS_c − kS_j (muestras de 16 ns) para pares (c en tiempo, j de hod2018 > 5 p.e.) del mismo grupo y de otros grupos.
# Salidas: <out>/h19_hvgroup_ratios.csv, h19_hvgroup_by_group.csv, h19_hvgroup_maxamp.csv, h19_hvgroup_dk.csv, h19_hvgroup.json.
# Uso: h19_hvgroup.py <work dir> <out dir>
import csv, json, math, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402

W, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
NBOOT = 200
TPE = 5.0
THR = (35.0, 3.0)
H18 = np.array([c for c in range(32) if c != 3])
H19 = np.arange(32, 64)
PEAKS = REPO + "/outputs/proceeding_hvdiag5/2026-10-31_full/work/hd64_%s_peaks.csv"
EV_DT = np.dtype([("t", "<f8"), ("aS", "<f4", 64), ("kS", "<u2", 64), ("aC", "<f4", 64)])
EDG = np.arange(0, 100.5, 0.5)
DK = np.arange(-12, 13)
RNG = np.random.default_rng(1103)

MAP = {c["dig_ch"]: c for c in json.load(open(REPO + "/daq/config/mapping/hodo_mapping_2026_10_r4.json"))["channels"]}
maps = {"sheet": {c: MAP[c]["hv_group_sheet"] for c in range(64)}, "oct_measured": {c: MAP[c]["hv_group"] for c in range(64)}}
g12 = {}
for st_, a, b, chs, why in mc.good_intervals("v8"):
    if st_ == "HV_OFF" and "group 2" in why:
        for c in chs.split():
            g12[int(c)] = 2
for st_, a, b, chs, why in mc.good_intervals("v8"):
    if "groups 1 and 2" in why:
        for c in chs.split():
            g12.setdefault(int(c), 1)
maps["may_measured_g12"] = {c: g12.get(c) for c in range(64)}


def g_proc(name, seeds):
    g = {}
    for l in open(PEAKS % name):
        f = l.strip().split(",")
        if f[0] != "ch" and int(f[0]) < 32:
            g[int(f[0])] = float(f[2]) if f[3] in ("ok", "changed") else float(f[1])
    return np.array([g[c] for c in range(32)])


def analyse(e, G, mp, T, win):
    a18 = e["aS" if win == "S" else "aC"][:, :32] / G[None, :]
    h18 = a18[:, H18] > TPE                      # (E, 31)
    a19 = e["aS" if win == "S" else "aC"][:, 32:]
    h19 = a19 > T                                # (E, 32)
    B = ~h19.any(1)
    grp18 = np.array([mp[j] if mp[j] is not None else -1 for j in H18])
    chans = [c for c in H19 if mp[c] is not None and (grp18 == mp[c]).any()]
    if mp is maps["may_measured_g12"]:
        known = np.ones(len(H18), bool)          # grupos 1 y 2 completos: el resto de hod2018 es «otro grupo»
    else:
        known = grp18 > 0
    same = np.stack([(grp18 == mp[c]) for c in chans], 1)        # (31, C)
    other = np.stack([known & (grp18 != mp[c]) for c in chans], 1)
    ns = h18.astype(np.int16) @ same.astype(np.int16)            # (E, C)
    no = h18.astype(np.int16) @ other.astype(np.int16)
    A = h19[:, [c - 32 for c in chans]]                         # (E, C)
    wts = np.vstack([np.ones(len(e))] + [RNG.poisson(1.0, len(e)) for _ in range(NBOOT)]).astype(float)  # (1+NBOOT, E)

    def stats(sel_ch):
        Aw = A[:, sel_ch]
        NA = wts @ Aw                                            # (R, C)
        SAs, SAo = wts @ (Aw * ns[:, sel_ch]), wts @ (Aw * no[:, sel_ch])
        NB = wts @ B
        SBs, SBo = wts @ (B[:, None] * ns[:, sel_ch]), wts @ (B[:, None] * no[:, sel_ch])
        obs_s, obs_o = SAs.sum(1), SAo.sum(1)
        exp_s, exp_o = (NA * SBs).sum(1) / NB, (NA * SBo).sum(1) / NB
        Rs, Ro = obs_s / exp_s, obs_o / exp_o
        DR = Rs / Ro
        r = {"pairs_A": float(NA[0].sum()), "N_B": float(NB[0]), "obs_same": float(obs_s[0]), "exp_same": float(exp_s[0]),
             "obs_other": float(obs_o[0]), "exp_other": float(exp_o[0])}
        for k, v in (("R_same", Rs), ("R_other", Ro), ("DR", DR)):
            r[k] = float(v[0])
            r[k + "_err"] = float(np.nanstd(v[1:], ddof=1))
        # ⟨n⟩ por par, por canal de hod2018 del grupo (para ver tamaños)
        r["mean_same_A"] = float(obs_s[0] / NA[0].sum()) if NA[0].sum() > 0 else np.nan
        r["mean_same_B"] = float(exp_s[0] / NA[0].sum()) if NA[0].sum() > 0 else np.nan
        r["mean_other_A"] = float(obs_o[0] / NA[0].sum()) if NA[0].sum() > 0 else np.nan
        r["mean_other_B"] = float(exp_o[0] / NA[0].sum()) if NA[0].sum() > 0 else np.nan
        return r

    res = {"n_triggers": int(len(e)), "n_with_h19": int((~B).sum()), "h19_channels": [int(c) for c in chans],
           "all": stats(np.arange(len(chans)))}
    res["by_group"] = {}
    for g in sorted(set(mp[c] for c in chans)):
        k = np.array([i for i, c in enumerate(chans) if mp[c] == g])
        res["by_group"][int(g)] = stats(k)
        res["by_group"][int(g)]["n_h18_same"] = int((grp18 == g).sum())
        res["by_group"][int(g)]["h19_channels"] = [int(chans[i]) for i in k]
    # amplitudes y Δk (solo ventana signal)
    if win == "S":
        mx = a18[:, H18].max(1)
        res["maxamp_all"] = np.histogram(np.clip(mx, 0, 99.99), EDG)[0]
        res["maxamp_h19"] = np.histogram(np.clip(mx[~B], 0, 99.99), EDG)[0]
        res["maxamp_q"] = {k: {"median": float(np.median(x)), "p90": float(np.percentile(x, 90)), "f_gt5": float((x > 5).mean()),
                               "f_gt10": float((x > 10).mean()), "n": int(len(x))} for k, x in (("all", mx), ("with_h19", mx[~B]), ("no_h19", mx[B]))}
        sm_A, sm_B = [], []
        aa = a18[:, H18]
        for i, c in enumerate(chans):
            m = np.where(same[:, i][None, :], aa, -np.inf).max(1)
            sm_A.append(m[A[:, i]])
            sm_B.append(m[B])
        sA, sB = np.concatenate(sm_A), np.concatenate(sm_B)
        res["maxamp_same_A"] = np.histogram(np.clip(sA, 0, 99.99), EDG)[0]
        res["maxamp_same_B"] = np.histogram(np.clip(sB, 0, 99.99), EDG)[0]
        res["maxamp_same_q"] = {k: {"median": float(np.median(x)), "p90": float(np.percentile(x, 90)), "f_gt5": float((x > 5).mean()), "n": int(len(x))}
                                for k, x in (("A_pairs", sA), ("B_pairs", sB))}
        k18 = e["kS"][:, H18].astype(int)
        dks, dko = [], []
        for i, c in enumerate(chans):
            ev = np.where(A[:, i])[0]
            if len(ev) == 0:
                continue
            d = e["kS"][ev, c].astype(int)[:, None] - k18[ev]
            hh = h18[ev]
            dks.append(d[hh & same[:, i][None, :]])
            dko.append(d[hh & other[:, i][None, :]])
        dks = np.concatenate(dks) if dks else np.array([], int)
        dko = np.concatenate(dko) if dko else np.array([], int)
        res["dk_same"] = np.histogram(np.clip(dks, DK[0], DK[-1]), np.append(DK, DK[-1] + 1) - 0.5)[0]
        res["dk_other"] = np.histogram(np.clip(dko, DK[0], DK[-1]), np.append(DK, DK[-1] + 1) - 0.5)[0]
        res["dk_q"] = {k: {"n": int(len(x)), "median": float(np.median(x)) if len(x) else np.nan, "mean": float(x.mean()) if len(x) else np.nan,
                           "sd": float(x.std()) if len(x) > 1 else np.nan, "f_abs_le1": float((np.abs(x - np.median(x)) <= 1).mean()) if len(x) else np.nan}
                       for k, x in (("same", dks), ("other", dko))}
    return res


ALL = {}
for camp, tag, peaks in (("may", "may8", "may_v8"), ("oct", "oct8", "oct8")):
    e = np.fromfile(f"{W}/{tag}_ev.bin", dtype=EV_DT)
    e = e[np.argsort(e["t"], kind="stable")]
    if camp == "may":
        e = e[mc.good_time(e["t"], "v8")]
    G = g_proc(peaks, None)
    for mname, mp in maps.items():
        for T in THR:
            for win in ("S", "C"):
                key = f"{camp}|{mname}|{T:g}|{win}"
                ALL[key] = analyse(e, G, mp, T, win)
                r = ALL[key]["all"]
                print(f"{camp} {mname:17s} T {T:>4g} mV {'en tiempo' if win == 'S' else 'fuera   '}: triggers {ALL[key]['n_triggers']}, con hod2019 "
                      f"{ALL[key]['n_with_h19']}, pares {r['pairs_A']:.0f}; R_same {r['R_same']:.3f} ± {r['R_same_err']:.3f}, R_other {r['R_other']:.3f} ± "
                      f"{r['R_other_err']:.3f}, DR {r['DR']:.3f} ± {r['DR_err']:.3f}", flush=True)

with open(f"{OUT}/h19_hvgroup_ratios.csv", "w") as f:
    f.write("# h19_hvgroup.py: hod2019 hit = a > T (mV) in channel c; hod2018 hit = a > 5 p.e. (G of PROC); window S = signal (in time), C = out-of-time "
            "control; R = observed / expected (expected from triggers without any hod2019 hit), DR = R_same / R_other; errors: Poisson bootstrap, "
            f"{NBOOT} replicas; May = run B good time, October = run 8 (after the hod2019 loss)\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "map", "T_mV", "window", "n_triggers", "n_with_h19", "pairs_A", "N_B", "obs_same", "exp_same", "R_same", "R_same_err",
                "obs_other", "exp_other", "R_other", "R_other_err", "DR", "DR_err", "mean_same_A", "mean_same_B", "mean_other_A", "mean_other_B"])
    for key, v in ALL.items():
        r = v["all"]
        w.writerow(key.split("|") + [v["n_triggers"], v["n_with_h19"], f"{r['pairs_A']:.0f}", f"{r['N_B']:.0f}", f"{r['obs_same']:.0f}",
                                     f"{r['exp_same']:.2f}", f"{r['R_same']:.4f}", f"{r['R_same_err']:.4f}", f"{r['obs_other']:.0f}", f"{r['exp_other']:.2f}",
                                     f"{r['R_other']:.4f}", f"{r['R_other_err']:.4f}", f"{r['DR']:.4f}", f"{r['DR_err']:.4f}", f"{r['mean_same_A']:.5f}",
                                     f"{r['mean_same_B']:.5f}", f"{r['mean_other_A']:.5f}", f"{r['mean_other_B']:.5f}"])
with open(f"{OUT}/h19_hvgroup_by_group.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "map", "T_mV", "window", "hv_group", "h19_channels", "n_h18_same", "pairs_A", "obs_same", "exp_same", "R_same", "R_same_err",
                "obs_other", "exp_other", "R_other", "R_other_err", "DR", "DR_err"])
    for key, v in ALL.items():
        for g, r in v["by_group"].items():
            w.writerow(key.split("|") + [g, " ".join(map(str, r["h19_channels"])), r["n_h18_same"], f"{r['pairs_A']:.0f}", f"{r['obs_same']:.0f}",
                                         f"{r['exp_same']:.2f}", f"{r['R_same']:.4f}", f"{r['R_same_err']:.4f}", f"{r['obs_other']:.0f}",
                                         f"{r['exp_other']:.2f}", f"{r['R_other']:.4f}", f"{r['R_other_err']:.4f}", f"{r['DR']:.4f}", f"{r['DR_err']:.4f}"])
with open(f"{OUT}/h19_hvgroup_maxamp.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "map", "T_mV", "bin_lo_pe", "bin_hi_pe", "all_triggers", "triggers_with_h19", "same_group_A_pairs", "same_group_B_pairs"])
    for key, v in ALL.items():
        camp, mname, T, win = key.split("|")
        if win != "S":
            continue
        for i in range(len(EDG) - 1):
            w.writerow([camp, mname, T, EDG[i], EDG[i + 1], int(v["maxamp_all"][i]), int(v["maxamp_h19"][i]), int(v["maxamp_same_A"][i]),
                        int(v["maxamp_same_B"][i])])
with open(f"{OUT}/h19_hvgroup_dk.csv", "w") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["campaign", "map", "T_mV", "dk_samples", "pairs_same_group", "pairs_other_groups"])
    for key, v in ALL.items():
        camp, mname, T, win = key.split("|")
        if win != "S":
            continue
        for i, d in enumerate(DK):
            w.writerow([camp, mname, T, int(d), int(v["dk_same"][i]), int(v["dk_other"][i])])


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, list):
        return [clean(x) for x in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


json.dump({"maps": {k: {str(c): v[c] for c in range(64)} for k, v in maps.items()}, "results": clean(ALL), "nboot": NBOOT, "T_hod2018_pe": TPE},
          open(f"{OUT}/h19_hvgroup.json", "w"), indent=1)
for key, v in ALL.items():
    if key.endswith("|S"):
        q, s, d = v["maxamp_q"], v["maxamp_same_q"], v["dk_q"]
        print(f"{key}: mayor pulso hod2018 (p.e.) mediana todos {q['all']['median']:.2f} / con hod2019 {q['with_h19']['median']:.2f} / sin "
              f"{q['no_h19']['median']:.2f}; f(>5) {q['all']['f_gt5']:.4f} / {q['with_h19']['f_gt5']:.4f} / {q['no_h19']['f_gt5']:.4f}; mismo grupo: "
              f"mediana A {s['A_pairs']['median']:.2f} B {s['B_pairs']['median']:.2f}, f(>5) A {s['A_pairs']['f_gt5']:.4f} B {s['B_pairs']['f_gt5']:.4f}; "
              f"Δk mismo grupo n {d['same']['n']} mediana {d['same']['median']:.1f} sd {d['same']['sd']:.2f}, otros n {d['other']['n']} mediana "
              f"{d['other']['median']:.1f} sd {d['other']['sd']:.2f}")
