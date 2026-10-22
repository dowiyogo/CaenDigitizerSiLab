#!/usr/bin/python3 -I
# biasscan.py — DAQ2-14: barrido de tensión de hod2018 y prueba de auto-trigger con pares de grupos X·Y.
# Corre sin atención (unidad systemd). Secuencia por escalón ΔV = 0, −0.5, −1.0, −1.5 V (CH00-07 juntos):
#   1. V0Set = nominal + ΔV con hvctl bias set (límites y verificación en hvmon), 2 min de estabilización;
#   2. SW trigger 100 Hz, registro completo, 2 min: G (reglas de la fase 5), pedestal, tasas oscuras y de
#      ruido por canal en p.e. de la G del escalón, hod2019 > 1 mV;
#   3. trigger externo 15 min: eficiencia X·Y de hod2018 sobre 5, 8, 10, 12 p.e. (accidentales con
#      eventos mezclados y con ventana fuera de tiempo), espectro en tiempo;
#   4. predicción por par de grupos (G0∧G2, G0∧G3, G1∧G2, G1∧G3; ch2, ch9, ch15 fuera): tasa de trigger
#      accidental 2W·r_X·r_Y, muones por hora y razón; umbral = máximo de muones/h con razón > 3 (si no
#      hay, el de mayor razón, marcado);
#   5. auto-trigger con ese umbral, 4 runs de 10 min (uno por par).
# Después: nominal restaurado (hvctl bias restore), run de física con trigger externo y c1, y un correo
# de resumen (hvctl notify-summary). Cualquier fallo de HV (rechazo, timeout, trip, Status inesperado):
# aborto, restauración del nominal y de Pw = 1, run de física y correo. Si un escalón terminaría después
# de las 10:00Z, se omiten los que faltan.
# Escrituras de HV: solo por hvctl (hvmon), cada una con --authorize. Registro: /data/daq/logs/biasscan.log
# y biasscan.jsonl. Resultados: outputs/daq2-14/.
# Uso: biasscan.py [--sim] [--campaign NOMBRE] [--deadline "AAAA-MM-DD HH:MM:SS"]
#   --sim: hvmon de prueba con FakeHVBackend, daqrun fake, duraciones cortas
import calendar, glob, itertools, json, math, os, statistics, subprocess, sys, time

SIM = "--sim" in sys.argv
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
BIN = os.path.expanduser("~/.local/opt/daqv2/bin")
HVCTL = ([REPO + "/daq/build/apps/hvctl/hvctl", "--socket", "/run/user/1000/daqv2-test/hvmon.sock"] if SIM else [BIN + "/hvctl"])
DAQRUN = (REPO + "/daq/build/apps/daqrun/daqrun") if SIM else BIN + "/daqrun"
def _arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


CAMPAIGN = _arg("--campaign", "biasscan_sim_2026_10_22" if SIM else "biasscan_2026_10_22")
AUTH = "DAQ2-14 bias scan and self-trigger experiment, authorized by René"
NOMINAL = [65.75, 65.75, 65.90, 65.85, 65.75, 65.70, 66.15, 65.90]
STEPS = [0.0, -0.5, -1.0, -1.5]
# Límite: las 10:00Z siguientes al arranque (07:00 Chile, antes de la luz de la mañana), o --deadline "AAAA-MM-DD HH:MM:SS" (UTC).
_now = time.gmtime()
_next10 = calendar.timegm((_now.tm_year, _now.tm_mon, _now.tm_mday, 10, 0, 0, 0, 0, 0))
DEADLINE = calendar.timegm(time.strptime(_arg("--deadline", ""), "%Y-%m-%d %H:%M:%S")) if "--deadline" in sys.argv else (
    _next10 if time.time() < _next10 else _next10 + 86400)
STEP_EST_S = 70 * 60 if not SIM else 60  # duración estimada de un escalón completo
DUR = {"stab": 120, "sw": 120, "ext": 900, "self": 600} if not SIM else {"stab": 5, "sw": 15, "ext": 20, "self": 15}
MASKED = [2, 9, 15]
PAIRS = [(0, 2), (0, 3), (1, 2), (1, 3)]
THRESHOLDS = [4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30]
W_S = 8e-9          # ventana de la mayoría: 1 reloj (0x810C[23:20] = 1)
FLUX_SQ_HZ = 10.89 / 60.0  # modelo: 1 cm^-2 min^-1 (PDG) × cuadrado activo de 10.89 cm²
ACC_GEOM = 0.0360   # fracción de triggers externos que cruza el cuadrado (MC del handoff)
OUT = REPO + "/outputs/daq2-14" + ("-sim" if SIM else "")
LOGD = "/data/daq/sim/logs" if SIM else "/data/daq/logs"
LOG, JL = LOGD + "/biasscan.log", LOGD + "/biasscan.jsonl"
CFGD = ("/data/daq/sim/config/" if SIM else "/data/daq/config/") + CAMPAIGN
RAW = ("/data/daq/sim/raw/" if SIM else "/data/daq/raw/") + CAMPAIGN
PHYS_CFG = "/data/daq/config/hodoscope_2026_10/physics_ext_B.json"
C1 = os.path.expanduser("~/.local/opt/daqv2/share/calibration/pe_calibration_2026_10_c1.json")
MAP = os.path.expanduser("~/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r4.json")
ENV = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser("~/root6/lib"), PATH=os.path.expanduser("~/root6/bin") + ":" + os.environ.get("PATH", ""))
ROOTINC = 'gSystem->AddIncludePath("-I%s/daq/core/include")' % REPO


class Abort(Exception):
    pass


def log(msg, **k):
    t = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    os.makedirs(LOGD, exist_ok=True)
    with open(LOG, "a") as f:
        f.write(f"{t} {msg}\n")
    k.update(t_utc=t, msg=msg)
    with open(JL, "a") as f:
        f.write(json.dumps(k, ensure_ascii=False) + "\n")


def sh(cmd, timeout=3600, check=True, heavy=False):
    if heavy:
        cmd = ["bash", "-c", "ulimit -v 4194304; exec nice -n 19 " + " ".join("'" + c.replace("'", "'\\''") + "'" for c in cmd)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=ENV)
    if check and r.returncode != 0:
        raise Abort(f"falló {' '.join(cmd)[:200]}: rc {r.returncode}: {(r.stdout + r.stderr)[-600:]}")
    return r


def hvctl(*args, timeout=3700):
    r = sh(HVCTL + list(args), timeout=timeout, check=False)
    out = (r.stdout or "").strip().splitlines()
    try:
        res = json.loads(out[-1]) if out else {}
    except json.JSONDecodeError:
        res = {"ok": False, "detail": r.stdout + r.stderr}
    if r.returncode != 0 and res.get("ok", True):
        res["ok"] = False
    return res


def hvstate():
    path = "/dev/shm/daqv2-sim-test/hvmon_state.json" if SIM else "/dev/shm/daqv2/hvmon_state.json"
    s = json.load(open(path))
    return {c["ch"]: c["values"] for c in s["channels"]}


def check_hv(where):
    st = hvstate()
    for ch in range(8):
        v = st[ch]
        if v.get("Pw") != 1 or v.get("Status") not in (1,):
            raise Abort(f"{where}: CH{ch:02d} Pw {v.get('Pw')} Status 0x{v.get('Status', 0):X} (se esperaba encendido y estable)")
    return st


def set_bias(delta):
    """Lleva CH00-07 a nominal + delta en escalones de ≤ 0.5 V (todos los canales por escalón)."""
    st = check_hv("antes de cambiar V0Set")
    cur = [st[ch]["V0Set_V"] for ch in range(8)]
    target = [round(NOMINAL[ch] + delta, 3) for ch in range(8)]
    while any(abs(cur[ch] - target[ch]) > 0.004 for ch in range(8)):
        for ch in range(8):
            if abs(cur[ch] - target[ch]) <= 0.004:
                continue
            step = max(-0.5, min(0.5, target[ch] - cur[ch]))
            v = round(cur[ch] + step, 3)
            res = hvctl("bias", "set", "--ch", str(ch), "--value", f"{v:.3f}", "--authorize", f"{AUTH}: CH{ch:02d} V0Set {v:.2f} V")
            log(f"V0Set CH{ch:02d} -> {v:.3f} V: {res.get('result')} ({res.get('detail', '')})", ev="v0set", ch=ch, value=v, res=res)
            if not res.get("ok") or res.get("result") not in ("verified", "noop"):
                raise Abort(f"V0Set CH{ch:02d} = {v:.3f} V: {res.get('result')}: {res.get('detail')}")
            cur[ch] = v
    st = check_hv("después de cambiar V0Set")
    return {ch: {"V0Set_V": st[ch]["V0Set_V"], "VMon_V": st[ch]["VMon_V"], "IMon_A": st[ch]["IMon_A"]} for ch in range(8)}


def write_cfg(name, base, mutate):
    c = json.load(open(base))
    c["campaign"]["name"] = CAMPAIGN
    c["campaign"]["operator"] = "claude (DAQ2-14, barrido de tensión y auto-trigger, autorizado por René)"
    c["notifier"]["enabled"] = False  # avisos de HV: hvmon; resumen: un correo al final
    c["mapping"]["file"] = MAP
    c["notifier"]["mapping_file"] = MAP
    if SIM:
        c["digitizer"]["backend"] = "fake"
        c["digitizer"]["scenario_file"] = REPO + "/daq/config/fake/dig_synthetic_default.json"
        c["digitizer"]["record_length_samples"] = 144 if c["digitizer"].get("record_length_samples", 144) < 300 else 624
        for k in ("raw_root", "manifests_dir"):
            c["campaign"].pop(k, None)
        c.setdefault("hv", {})["state_file"] = "/dev/shm/daqv2-sim-test/hvmon_state.json"
    mutate(c)
    os.makedirs(CFGD, exist_ok=True)
    p = f"{CFGD}/{name}.json"
    json.dump(c, open(p, "w"), indent=1, ensure_ascii=False)
    return p


BOARD_GAP_S = 10 if not SIM else 1   # mínimo entre cerrar y volver a abrir el VX1740 (DAQ2-14: dos cuelgues)
OPEN_RETRIES = 3
OPEN_RETRY_WAIT_S = 30 if not SIM else 1
_last_close = [0.0]
OPEN_FAIL = ("configuración del digitizador", "CAEN_DGTZ_CommError", "CAEN_DGTZ_Timeout", "OpenDigitizer")


def wait_board_gap(what):
    """Espera a que no quede ningún daqrun y a que pasen BOARD_GAP_S desde el último cierre del digitizador."""
    t0 = time.time()
    while sh(["pgrep", "-f", DAQRUN], check=False).returncode == 0:
        if time.time() - t0 > 330:
            raise Abort(f"{what}: sigue corriendo otro daqrun después de 330 s")
        time.sleep(2)
    wait = _last_close[0] + BOARD_GAP_S - time.time()
    if wait > 0:
        time.sleep(wait)


def daqrun(cfg, seconds):
    """Un run de daqrun; la apertura del digitizador se reintenta OPEN_RETRIES veces si falla por
    comunicación (antes de "digitizador abierto"); siempre BOARD_GAP_S entre cierre y apertura."""
    for attempt in range(1, OPEN_RETRIES + 1):
        try:
            return _daqrun_once(cfg, seconds)
        except OpenFailed as e:
            log(f"  apertura del digitizador falló (intento {attempt} de {OPEN_RETRIES}): {e}")
            if attempt == OPEN_RETRIES:
                raise Abort(f"el digitizador no abre tras {OPEN_RETRIES} intentos: {e}")
            time.sleep(OPEN_RETRY_WAIT_S)


class OpenFailed(Exception):
    pass


def _daqrun_once(cfg, seconds):
    """Un run de daqrun; mientras corre, HV de CH00-07 vigilada cada 5 s: un trip o Status distinto de 0x1
    detiene el run (SIGTERM, cierre limpio) y aborta la secuencia."""
    wait_board_gap(os.path.basename(cfg))
    t0 = time.time()
    logf = f"{OUT}/daqrun_{os.path.basename(cfg)}.log"
    with open(logf, "w") as lf:
        p = subprocess.Popen([DAQRUN, "--config", cfg, "--max-duration", str(seconds), "--no-http"], stdout=lf, stderr=subprocess.STDOUT, env=ENV)
        hv_err = None
        while p.poll() is None:
            time.sleep(5)
            try:
                check_hv("durante el run")
            except Abort as e:
                hv_err = str(e)
                p.terminate()
                p.wait(timeout=300)
                break
            if time.time() - t0 > seconds + 600:
                p.terminate()
                p.wait(timeout=300)
                raise Abort(f"daqrun {os.path.basename(cfg)} no terminó a tiempo")
    _last_close[0] = time.time()
    txt = open(logf).read()
    if p.returncode != 0 and "digitizador abierto" not in txt and any(k in txt for k in OPEN_FAIL):
        raise OpenFailed(next((l for l in txt.splitlines() if "ERROR" in l), txt[-300:])[-300:])
    if hv_err:
        raise Abort(f"{hv_err} (run {os.path.basename(cfg)} detenido)")
    r = subprocess.CompletedProcess(p.args, p.returncode)
    base = None
    for line in txt.splitlines():
        if "daqrun arranca" in line and "(" in line:
            base = line.split("(")[-1].split(")")[0]
        if "daqrun termina" in line:
            log("  " + line[line.find("daqrun termina"):][:300])
    if r.returncode != 0 or not base:
        raise Abort(f"daqrun {os.path.basename(cfg)} falló (rc {r.returncode}): {txt[-500:]}")
    segs = sorted(glob.glob(f"{RAW}/{base}_seg*.root"))
    idx = json.load(open(f"{RAW}/{base}.index.json"))
    return {"base": base, "segments": segs, "events": idx.get("n_events", 0), "wall_s": time.time() - t0,
            "trigger_signature": idx.get("trigger_signature", "")}


def root(macro, *args):
    a = ",".join(json.dumps(x) if isinstance(x, str) else str(x).lower() if isinstance(x, bool) else str(x) for x in args)
    sh(["root", "-l", "-b", "-q", "-e", ROOTINC, f"{REPO}/daq/tools/validation/{macro}.C+({a})"], timeout=3000, heavy=True)


# ---------------------------------------------------------------- análisis de cada escalón

def analyse_sw(run, step, seedsCal):
    pat = RAW + "/" + run["base"] + "_seg*.root"
    seeds = f"{OUT}/seeds_step{step}.json"
    cal0 = json.load(open(seedsCal))
    json.dump([{"ch": c["ch"], "status": "ok", "value": c["G_mV"]} for c in cal0["channels"] if c.get("G_mV")], open(seeds, "w"))
    cal = f"{OUT}/pecal_step{step}.json"
    root("commissioning_pe_calibration", pat, f"daq2-14-step{step}", cal, seeds)
    noise = f"{OUT}/noise_step{step}.csv"
    root("selftrig_noise", pat, "", MAP, noise, 64)
    c = json.load(open(cal))
    G = {e["ch"]: e.get("G_mV") for e in c["channels"]}
    ok = [G[ch] for ch in range(32) if G.get(ch) and c["channels"][ch].get("status") in ("ok", "changed")]
    med = statistics.median(ok) if ok else statistics.median([e["G_mV"] for e in cal0["channels"] if e.get("G_mV")])
    # calibración del escalón para daqrun: canales sin G aceptada -> mediana del escalón (marcados)
    for e in c["channels"]:
        if e["ch"] < 32 and not (e.get("G_mV") and e.get("status") in ("ok", "changed")):
            e["G_mV_fallback"] = med
            e["G_mV"] = med
            e["provisional"] = True
    json.dump(c, open(cal, "w"), indent=1)
    exp, spec = {}, {}
    for line in open(noise):
        if line.startswith("#exposure"):
            _, ch, s = line.split()
            exp[int(ch)] = float(s)
        elif line[0].isdigit():
            ch, b, n = line.strip().split(",")
            spec.setdefault(int(ch), []).append((float(b), int(n)))
    def rate(ch, thr_mv):
        return sum(n for b, n in spec.get(ch, []) if b >= thr_mv) / exp[ch] if exp.get(ch) else float("nan")
    per = {}
    for ch in range(64):
        g = c["channels"][ch].get("G_mV") if ch < 32 else None
        row = {"pedestal_mV": c["channels"][ch].get("pedestal_mV"), "status": c["channels"][ch].get("status"), "G_mV": g,
               "G_measured": ch < 32 and c["channels"][ch].get("status") in ("ok", "changed") and "G_mV_fallback" not in c["channels"][ch]}
        if ch < 32 and g:
            for T in [0.5, 1.5] + THRESHOLDS:
                row[f"rate_{T}pe_Hz"] = rate(ch, T * g)
        if ch >= 32:
            row["rate_1mV_Hz"] = rate(ch, 1.0)
        per[ch] = row
    groups = {}
    for T in THRESHOLDS:
        groups[T] = [sum(per[ch][f"rate_{T}pe_Hz"] for ch in range(8 * g, 8 * g + 8) if ch not in MASKED and per[ch].get(f"rate_{T}pe_Hz") == per[ch].get(f"rate_{T}pe_Hz")) for g in range(4)]
    return {"calibration": cal, "per_channel": per, "group_rates_masked_Hz": groups, "exposure_s": exp.get(0), "G_median_mV": med,
            "n_G_ok": len(ok)}


def analyse_ext(run, step, cal):
    pat = RAW + "/" + run["base"] + "_seg*.root"
    res = {}
    for tag, mask in (("all", ""), ("masked", ",".join(map(str, MASKED)))):
        f = f"{OUT}/muon_step{step}_{tag}.csv"
        root("selftrig_muon", pat, cal, f, mask)
        rows = [{k: float(v) for k, v in r.items()} for r in __import__("csv").DictReader(open(f))]
        N = len(rows)
        out = {"N": N}
        for T in THRESHOLDS:
            sx = [max(r["s0_pe"], r["s1_pe"]) > T for r in rows]
            sy = [max(r["s2_pe"], r["s3_pe"]) > T for r in rows]
            cx = [max(r["c0_pe"], r["c1_pe"]) > T for r in rows]
            cy = [max(r["c2_pe"], r["c3_pe"]) > T for r in rows]
            n = sum(a and b for a, b in zip(sx, sy))
            acc_ctl = sum(a and b for a, b in zip(sx, cy)) + sum(a and b for a, b in zip(cx, sy)) - sum(a and b for a, b in zip(cx, cy))
            acc_mix = 0.5 * (sum(sx[i] and sy[i + 1] for i in range(N - 1)) + sum(sy[i] and sx[i + 1] for i in range(N - 1)))
            d = {"N_XY": n, "acc_ctl": acc_ctl, "acc_mix": acc_mix, "frac_raw": n / N if N else 0,
                 "eff_ctl": (n - acc_ctl) / (ACC_GEOM * N) if N else 0, "eff_mix": (n - acc_mix) / (ACC_GEOM * N) if N else 0}
            pairs = {}
            for gx, gy in PAIRS:
                a = [r[f"s{gx}_pe"] > T for r in rows]; b = [r[f"s{gy}_pe"] > T for r in rows]
                ca = [r[f"c{gx}_pe"] > T for r in rows]; cb = [r[f"c{gy}_pe"] > T for r in rows]
                np_ = sum(x and y for x, y in zip(a, b))
                acc = sum(x and y for x, y in zip(a, cb)) + sum(x and y for x, y in zip(ca, b)) - sum(x and y for x, y in zip(ca, cb))
                pairs[f"G{gx}G{gy}"] = {"N": np_, "acc_ctl": acc, "eff": (np_ - acc) / (ACC_GEOM * N) if N else 0}
            d["pairs"] = pairs
            out[T] = d
        if tag == "all":
            sel = [r for r in rows if max(r["s2_pe"], r["s3_pe"]) > 5]
            xs = sorted(max(r["s0_pe"], r["s1_pe"]) for r in sel)
            out["spectrum_X_given_Y5_quartiles_pe"] = statistics.quantiles(xs, n=4) if len(xs) >= 4 else xs
            out["spectrum_X_given_Y5_hist"] = [sum(1 for x in xs if lo <= x < lo + 2) for lo in range(0, 60, 2)]
        res[tag] = out
    return res


def predict(sw, ext, dur_h_pair=1.0):
    table = []
    for T in THRESHOLDS:
        r = sw["group_rates_masked_Hz"][T]
        tot_mu = tot_acc = 0
        rows = {}
        for gx, gy in PAIRS:
            acc = 2 * W_S * r[gx] * r[gy]
            eff = max(0.0, ext["masked"][T]["pairs"][f"G{gx}G{gy}"]["eff"])
            mu = FLUX_SQ_HZ * eff
            rows[f"G{gx}G{gy}"] = {"acc_Hz": acc, "muon_Hz": mu, "muons_per_h": mu * 3600, "ratio": mu / acc if acc > 0 else float("inf"),
                                   "trigger_Hz": acc + mu}
            tot_mu += mu
            tot_acc += acc
        table.append({"T_pe": T, "pairs": rows, "muons_per_h_sum": tot_mu * 3600, "ratio_sum": tot_mu / tot_acc if tot_acc > 0 else float("inf"),
                      "max_trigger_Hz": max(v["trigger_Hz"] for v in rows.values())})
    ok = [t for t in table if t["ratio_sum"] > 3 and t["max_trigger_Hz"] < 50 and t["muons_per_h_sum"] > 0]
    if ok:
        best = max(ok, key=lambda t: t["muons_per_h_sum"])
        why = "máximo de muones/h con razón > 3"
    else:
        cand = [t for t in table if t["max_trigger_Hz"] < 50 and t["muons_per_h_sum"] > 0] or table
        best = max(cand, key=lambda t: t["ratio_sum"])
        why = "NINGÚN umbral da razón > 3: se usa el de mayor razón"
    return table, best["T_pe"], why


def analyse_self(run, step, T, gx, gy, cal, sw):
    pat = RAW + "/" + run["base"] + "_seg*.root"
    f = f"{OUT}/self_step{step}_G{gx}G{gy}.csv"
    root("selftrig_check", pat, cal, f, 3.0)
    import csv, collections
    rows = [(int(r["ev"]), int(r["ch"]), float(r["amp_pe"]), int(r["kmin"])) for r in csv.DictReader(open(f))]
    inpair = lambda ch: ch // 8 in (gx, gy) and ch not in MASKED
    k = collections.Counter(km for e, ch, a, km in rows if a > T and inpair(ch))
    mode = k.most_common(1)[0][0] if k else 50
    lo, hi = mode - 4, mode + 5
    clo, chi = lo + 50, hi + 50
    inw = collections.defaultdict(set); ctl = collections.defaultdict(set); amps = []
    for e, ch, a, km in rows:
        if a > T and inpair(ch):
            if lo <= km < hi: inw[e].add(ch // 8)
            if clo <= km < chi: ctl[e].add(ch // 8)
    n = run["events"]
    nxy = sum(1 for e in range(n) if {gx, gy} <= inw[e])
    acc_ctl = sum(1 for e in range(n) if gx in inw[e] and gy in ctl[e]) + sum(1 for e in range(n) if gy in inw[e] and gx in ctl[e]) \
        - sum(1 for e in range(n) if {gx, gy} <= ctl[e])
    r = sw["group_rates_masked_Hz"][T]
    dur = DUR["self"]
    acc_pred = 2 * W_S * r[gx] * r[gy] * dur
    for e, ch, a, km in rows:
        if inpair(ch) and lo <= km < hi and ch // 8 == gx and {gx, gy} <= inw[e]:
            amps.append(a)
    return {"pair": f"G{gx}G{gy}", "events": n, "rate_Hz": n / dur, "window": [lo, hi], "N_XY_in_time": nxy, "acc_pred": acc_pred,
            "acc_offtime": acc_ctl, "ratio_pred": nxy / acc_pred if acc_pred > 0 else float("inf"),
            "muons_per_h": (nxy - acc_pred) / dur * 3600, "trigger_signature": run["trigger_signature"],
            "amp_X_quartiles_pe": statistics.quantiles(sorted(amps), n=4) if len(amps) >= 4 else sorted(amps)}


# ---------------------------------------------------------------- secuencia

PHYS_LOG = "/data/daq/logs/daq-physics.log"


def physics_status():
    """Estado del run de física: unidad activa y última línea relevante de su log."""
    if SIM:
        return {"active": False, "last": "(sim)"}
    act = sh(["systemctl", "--user", "is-active", "daq-physics"], check=False).stdout.strip()
    last = ""
    try:
        lines = [l for l in open(PHYS_LOG, errors="replace").read().splitlines()[-400:] if "ERROR" in l or "termina" in l or "arranca" in l or "abierto" in l]
        last = lines[-1][:300] if lines else ""
    except OSError:
        pass
    return {"active": act == "active", "state": act, "last": last}


def physics_stop():
    """Detiene el run de física y espera a que el daqrun termine (cierre del digitizador)."""
    if SIM:
        return {"was_active": False}
    before = physics_status()
    log(f"run de física antes del barrido: {before['state']}; última línea: {before['last']}")
    if before["active"]:
        sh(["systemctl", "--user", "stop", "daq-physics"], timeout=400, check=False)
    t0 = time.time()
    while sh(["pgrep", "-f", DAQRUN], check=False).returncode == 0:
        if time.time() - t0 > 330:
            raise Abort("el daqrun de física no terminó en 330 s")
        time.sleep(2)
    _last_close[0] = time.time()
    after = physics_status()
    log(f"run de física detenido ({after['state']}); última línea: {after['last']}")
    return {"was_active": before["active"], "before_last": before["last"], "after_last": after["last"]}


def physics_start():
    if SIM:
        log("(sim) el run de física se reiniciaría aquí con " + PHYS_CFG)
        return
    cmd = ["systemd-run", "--user", "--unit=daq-physics", "--collect", "--property=RuntimeMaxSec=1209600", "--property=KillSignal=SIGTERM",
           "--property=TimeoutStopSec=300", "--property=StandardOutput=append:/data/daq/logs/daq-physics.log",
           "--property=StandardError=append:/data/daq/logs/daq-physics.log", "--setenv=LD_LIBRARY_PATH=" + os.path.expanduser("~/root6/lib"),
           "--working-directory=" + os.path.expanduser("~"), BIN + "/daqrun", "--config", PHYS_CFG]
    for attempt in range(1, OPEN_RETRIES + 1):
        wait_board_gap("run de física")
        size = os.path.getsize(PHYS_LOG) if os.path.exists(PHYS_LOG) else 0
        r = sh(cmd, timeout=60, check=False)
        time.sleep(15)
        act = sh(["systemctl", "--user", "is-active", "daq-physics"], check=False).stdout.strip()
        new_txt = ""
        try:
            with open(PHYS_LOG, errors="replace") as f:
                f.seek(size)
                new_txt = f.read()
        except OSError:
            pass
        opened = "digitizador abierto" in new_txt
        err = next((l for l in new_txt.splitlines() if "ERROR" in l), "")[:300]
        log(f"run de física (intento {attempt}): unidad {act}, digitizador {'abierto' if opened else 'NO abierto'} {err}")
        if act == "active" and opened:
            return f"activo con trigger externo (intento {attempt})"
        if act != "active":
            _last_close[0] = time.time()
        if attempt < OPEN_RETRIES:
            time.sleep(OPEN_RETRY_WAIT_S)
    return f"NO ARRANCÓ tras {OPEN_RETRIES} intentos: {err or act} (el digitizador puede necesitar reiniciar el crate VME)"


def summary_mail(text):
    p = f"{OUT}/summary_mail.txt"
    open(p, "w").write(text)
    res = hvctl("notify-summary", "--text-file", p, "--kind", "bias_scan_summary", timeout=240)
    log(f"correo de resumen: {res}")


def main():
    os.makedirs(OUT, exist_ok=True)
    log("=== DAQ2-14 barrido de tensión" + (" (SIMULACIÓN)" if SIM else "") + f"; campaña {CAMPAIGN}; límite "
        + time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(DEADLINE)))
    results = {"steps": [], "aborted": None, "skipped": []}
    phys_stopped = False
    try:
        st = hvstate()
        read = [st[ch]["V0Set_V"] for ch in range(8)]
        log(f"V0Set leído: {read}")
        if any(abs(a - b) > 0.005 for a, b in zip(read, NOMINAL)):
            raise Abort(f"V0Set leído {read} distinto del nominal esperado {NOMINAL}: no se empieza")
        check_hv("inicio")
        # Campaña y configuraciones base
        camp = ["daqcampaign"] + (["--sim"] if SIM else ["--authorize", AUTH])
        sh([(REPO + "/daq/build/apps/daqcampaign/" if SIM else BIN + "/") + "daqcampaign"] + camp[1:] + ["open", CAMPAIGN], check=False)
        swcfg = write_cfg("sw", "/data/daq/config/commissioning3_2026_10_21/toggle2_sw100_full.json", lambda c: None)
        extcfg = write_cfg("ext", PHYS_CFG, lambda c: c["run"].__setitem__("max_duration_s", 0))
        results["physics_before"] = physics_stop()
        phys_stopped = True
        res = hvctl("bias", "begin", "--expected", ",".join(f"{v:.2f}" for v in NOMINAL), "--authorize", AUTH)
        log(f"bias begin: {res.get('result')} ({res.get('detail')})", res=res)
        if not res.get("ok"):
            raise Abort("bias begin: " + str(res.get("detail")))
        prevCal = C1
        for i, dv in enumerate(STEPS):
            if time.time() + STEP_EST_S > DEADLINE:
                results["skipped"] = STEPS[i:]
                log(f"escalones {STEPS[i:]} omitidos: terminarían después de las 10:00Z")
                break
            step = {"delta_V": dv, "t_start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
            log(f"--- escalón ΔV = {dv:+.1f} V")
            step["hv_after_set"] = set_bias(dv)
            time.sleep(DUR["stab"])
            step["hv_before_sw"] = {ch: v for ch, v in check_hv("antes del run SW").items() if ch < 8}
            sw = daqrun(swcfg, DUR["sw"])
            step["sw_run"] = sw
            step["sw"] = analyse_sw(sw, i, prevCal)
            prevCal = step["sw"]["calibration"]
            log(f"  SW: {sw['events']} eventos; G mediana {step['sw']['G_median_mV']:.3f} mV ({step['sw']['n_G_ok']} canales con G)")
            check_hv("antes del run externo")
            ext = daqrun(write_cfg(f"ext_step{i}", extcfg, lambda c, cal=prevCal: c["calibration"].__setitem__("file", cal)), DUR["ext"])
            step["ext_run"] = ext
            step["ext"] = analyse_ext(ext, i, prevCal)
            e10 = step["ext"]["all"].get(10, {})
            log(f"  externo: {ext['events']} eventos; ε_XY(10 p.e.) ctl {e10.get('eff_ctl', 0):.3f}, mezclados {e10.get('eff_mix', 0):.3f}")
            table, T, why = predict(step["sw"], step["ext"])
            step["prediction"] = table
            step["T_chosen_pe"] = T
            step["T_choice"] = why
            log(f"  umbral elegido {T} p.e. ({why})")
            step["self"] = []
            for gx, gy in PAIRS:
                check_hv(f"antes del auto-trigger G{gx}∧G{gy}")
                def mut(c, gx=gx, gy=gy, cal=prevCal):
                    c["calibration"]["file"] = cal
                    c["digitizer"]["trigger"] = {"mode": "self_hod2018", "self": {"groups": [gx, gy], "threshold_pe": T, "masked_channels": MASKED,
                                                                                  "majority_groups": 2, "window_clk": 1, "include_external": False,
                                                                                  "target_rate_Hz": 50}}
                run = daqrun(write_cfg(f"self_step{i}_G{gx}G{gy}", extcfg, mut), DUR["self"])
                a = analyse_self(run, i, T, gx, gy, prevCal, step["sw"])
                step["self"].append(a)
                log(f"  auto-trigger G{gx}∧G{gy} T={T}: {run['events']} eventos ({a['rate_Hz']:.2f} Hz), X·Y en tiempo {a['N_XY_in_time']}, "
                    f"accidentales previstos {a['acc_pred']:.1f}, razón {a['ratio_pred']:.2f}, muones/h {a['muons_per_h']:.1f}")
            step["self_muons_per_h_sum"] = sum(a["muons_per_h"] for a in step["self"])
            step["t_end_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            results["steps"].append(step)
            json.dump(results, open(f"{OUT}/biasscan_results.json", "w"), indent=1, ensure_ascii=False, default=str)
    except Exception as e:  # noqa: BLE001 — cualquier fallo: aborto con restauración
        results["aborted"] = str(e)
        log(f"ABORTO: {e}")
    finally:
        res = hvctl("bias", "restore", "--authorize", AUTH + ": restauración del nominal")
        log(f"bias restore: {res.get('result')} ({res.get('detail')})", res=res)
        results["restore"] = res
        try:
            time.sleep(5)  # el archivo de estado de hvmon se publica en el ciclo siguiente
            st = hvstate()
            results["hv_final"] = {ch: {"V0Set_V": st[ch]["V0Set_V"], "VMon_V": st[ch]["VMon_V"], "Pw": st[ch]["Pw"], "Status": st[ch]["Status"]} for ch in range(8)}
            bad = [ch for ch in range(8) if abs(st[ch]["V0Set_V"] - NOMINAL[ch]) > 0.005 or st[ch]["Pw"] != 1]
            log(f"HV final: {'nominal y Pw = 1 en CH00-07' if not bad else 'NO nominal en ' + str(bad)}", hv=results["hv_final"])
        except Exception as e:  # noqa: BLE001
            log(f"no se pudo leer el estado final de HV: {e}")
        # Siempre: el run de física vuelve con trigger externo (también si no corría antes o si se abortó).
        try:
            results["physics"] = physics_start()
        except Exception as e:  # noqa: BLE001
            results["physics"] = f"NO ARRANCÓ: {e}"
        log(f"run de física al final: {results['physics']}")
        json.dump(results, open(f"{OUT}/biasscan_results.json", "w"), indent=1, ensure_ascii=False, default=str)
        lines = [f"DAQ2-14 barrido de tensión de hod2018{' (SIMULACIÓN)' if SIM else ''}: " + ("ABORTADO: " + results["aborted"] if results["aborted"] else "completo"),
                 f"Escalones hechos: {[s['delta_V'] for s in results['steps']]}; omitidos: {results['skipped']}",
                 f"Restauración: {results['restore'].get('result')} {results['restore'].get('detail', '')}",
                 f"Run de física antes del barrido: {'activo' if results.get('physics_before', {}).get('was_active') else 'NO activo'}"
                 f" (última línea del log: {results.get('physics_before', {}).get('before_last', '–')})",
                 f"Run de física al final: {results.get('physics')}"]
        for s in results["steps"]:
            e = s["ext"]["all"]
            lines.append(f"ΔV {s['delta_V']:+.1f} V: G mediana {s['sw']['G_median_mV']:.3f} mV; ε_XY ctl 5/8/10/12 p.e. "
                         + "/".join(f"{e[T]['eff_ctl']:.3f}" for T in (5, 8, 10, 12) if T in e)
                         + f"; umbral {s['T_chosen_pe']} p.e. ({s['T_choice']}); muones/h auto-trigger (suma de pares) {s['self_muons_per_h_sum']:.1f}")
            for a in s["self"]:
                lines.append(f"   {a['pair']}: {a['rate_Hz']:.2f} Hz, X·Y {a['N_XY_in_time']}, accidentales previstos {a['acc_pred']:.1f}, razón {a['ratio_pred']:.2f}")
        try:
            sh(["python3", "-I", REPO + "/daq/tools/biasscan/biasscan_report.py", OUT, JL], timeout=600)
            ana = json.load(open(f"{OUT}/biasscan_analysis.json"))
            vbd = [(ch, r["V_bd"], r["sigma_V_bd"], r["overvoltage_nominal"]) for ch, r in ana["channels"].items() if "V_bd" in r]
            if vbd:
                lines.append(f"V_bd (estimada) en {len(vbd)} canales: mediana {statistics.median(v[1] for v in vbd):.2f} V; sobretensión al nominal "
                             f"mediana {statistics.median(v[3] for v in vbd):.2f} V")
            lines.append(f"Informe: {OUT}/report_es.md")
        except Exception as e:  # noqa: BLE001
            lines.append(f"informe automático falló: {e}")
        lines.append(f"Registro: {LOG}; resultados: {OUT}/biasscan_results.json")
        try:
            summary_mail("\n".join(lines))
        except Exception as e:  # noqa: BLE001
            log(f"correo de resumen falló: {e}")
        log("=== fin")


def _sigterm(signum, frame):
    raise Abort(f"señal {signum} (unidad detenida): se restaura el nominal")


if __name__ == "__main__":
    import signal
    signal.signal(signal.SIGTERM, _sigterm)
    main()
