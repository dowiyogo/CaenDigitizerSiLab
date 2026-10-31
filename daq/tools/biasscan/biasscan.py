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
# Política de fallos (René, 2026-10-23), tres clases:
#   (a) estados esperados (bits de rampa, estabilización): se espera con timeout (6τ, τ = 115 s);
#   (b) fallos recuperables (VMon un poco fuera de tolerancia, error al abrir el digitizador, un run que
#       no arranca, análisis que falla): hasta 3 intentos con 90 s entre ellos; si siguen, el escalón se
#       omite con su motivo, se guarda lo ya tomado y se sigue con el siguiente;
#   (c) riesgos reales (bits de trip o sobrecorriente, Pw que cae a 0 solo, conexión HV perdida, VMon en
#       contra del escalón): aborto, nominal restaurado y run de física rearrancado.
# DAQ2-15 (2026-10-23): barridos de día y de noche. Antes de los runs de auto-trigger de cada escalón, un run
# de control de 2 min con un par a 8 p.e.: si la tasa medida es < 10 % de la tasa accidental prevista, la
# configuración del auto-trigger se marca como que no funciona y se omiten los runs de auto-trigger de ese
# escalón (se sigue con el siguiente). Temperatura de la placa A7030DP y hora del día en cada escalón.
# DAQ2-16 (2026-10-23): eficiencia X·Y con el estimador compartido (daq/tools/characterization/xy_estimator.py) y la
# definición única de impacto (hitdef_extract.C: pulso con mínimo en la ventana signal, pedestal = media recortada de las
# 40 muestras previas); ventana del ajuste de G = [0.5, 6] × G esperada (la del escalón anterior medida, o c1 en ΔV = 0);
# el control del auto-trigger registra los registros escritos y releídos (run_info_json) y es "no comprobable" si la
# predicción da < 3 eventos; --self-check-only: solo el run de control, sin los runs de auto-trigger.
# Escrituras de HV: solo por hvctl (hvmon), cada una con --authorize. Registro: /data/daq/logs/biasscan_<campaña>.log
# y .jsonl. Resultados: outputs/daq2-15/<campaña>/.
# DAQ2-17 (2026-10-23, barrido hacia arriba autorizado por René): --steps "0,0.25,0.5,0.75,1.0", --max-above 1.0 y
# --max-step 0.25 (pasados a hvctl bias begin; el tope de 1.0 V y el escalón de 0.25 V sobre el nominal están en hvmon),
# --authorize "<texto>", --ext-s 2700 (run externo de 45 min), --trigger-diag (en ΔV = 0: runs de 30 s de cada grupo
# solo y de pares con ventanas 1-8 relojes). Siempre: tras cada cambio, |VMon − V0Set − desvío| ≤ 0.05 V; guarda de
# corriente: IMon > 12 µA (80 % de I0Set) en algún grupo durante 30 s -> no se pasa al escalón siguiente; IMon ≥ 14 µA
# -> aborto con restauración inmediata. Al final la restauración se reintenta hasta dejar todo en el nominal.
# DAQ2-17 (revisión de René): la guarda de corriente corre además en un hilo cada 5 s durante todo el barrido (también en
# los análisis entre runs); IMon ≥ 14 µA -> SIGUSR1 al proceso -> Abort en el hilo principal (restauración). El desvío
# VMon − V0Set que se aleje más de 0.15 V del desvío al nominal también aborta (hvmon y check_vmon).
# Uso: biasscan.py [--sim] [--campaign NOMBRE] [--deadline "AAAA-MM-DD HH:MM:SS" | --max-hours H] [--self-check-only]
#        [--steps L] [--max-above V] [--max-step V] [--authorize T] [--ext-s S] [--trigger-diag]
#   --sim: hvmon de prueba con FakeHVBackend, daqrun fake, duraciones cortas
import calendar, glob, itertools, json, math, os, signal, statistics, subprocess, sys, threading, time

SIM = "--sim" in sys.argv
SELF_CHECK_ONLY = "--self-check-only" in sys.argv
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"  # salidas (outputs/daq2-15/…) y binarios de prueba de la simulación
# BIASSCAN_CODE_ROOT (seguimiento 2026-10-31 (4)): el código (este script, xy_estimator.py, macros de daq/tools, cabeceras de daq/core/include)
# se toma del árbol donde está este archivo: el repositorio vivo o una copia fija instalada con daq/tools/biasscan/install_pinned.sh en
# ~/.local/share/daqv2/biasscan/<commit>/ (con BIASSCAN_COMMIT dentro). Las salidas siguen yendo a REPO.
CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
CODE_COMMIT = open(os.path.join(CODE, "BIASSCAN_COMMIT")).read().split()[0] if os.path.exists(os.path.join(CODE, "BIASSCAN_COMMIT")) else "árbol de trabajo"
sys.path.insert(0, CODE + "/daq/tools/characterization")
import xy_estimator as xy  # noqa: E402  estimador X·Y compartido (DAQ2-16)
BIN = os.path.expanduser("~/.local/opt/daqv2/bin")
HVCTL = ([REPO + "/daq/build/apps/hvctl/hvctl", "--socket", "/run/user/1000/daqv2-test/hvmon.sock"] if SIM else [BIN + "/hvctl"])
DAQRUN = (REPO + "/daq/build/apps/daqrun/daqrun") if SIM else BIN + "/daqrun"
def _arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


CAMPAIGN = _arg("--campaign", "biasscan_sim_2026_10_22" if SIM else "biasscan_2026_10_22")
AUTH = _arg("--authorize", "DAQ2-15 day/night bias scans (bounded V0Set and Pw as in DAQ2-14), authorized by René")
NOMINAL = [65.75, 65.75, 65.90, 65.85, 65.75, 65.70, 66.15, 65.90]
STEPS = [float(x) for x in _arg("--steps", "0,-0.5,-1.0,-1.5").split(",")]
MAX_ABOVE_V = float(_arg("--max-above", "0"))   # DAQ2-17: V0Set máximo sobre el nominal en la sesión (hvmon: tope 1.0 V)
MAX_STEP_V = float(_arg("--max-step", "0.5"))   # escalón de cada escritura (hvmon: ≤ 0.25 V sobre el nominal)
TRIGGER_DIAG = "--trigger-diag" in sys.argv
VMON_TOL_V = 0.05          # DAQ2-17: |VMon − V0Set − desvío del último escalón verificado| tras estabilizar
OFFSET_DRIFT_MAX_V = 0.15  # DAQ2-17: |VMon − V0Set − desvío al nominal| máximo (si no: aborto y restauración)
IMON_POLL_S = 5.0          # período del hilo de vigilancia de corriente
IMON_WARN_A, IMON_STOP_A, IMON_WARN_S = 12e-6, 14e-6, 30.0  # guarda de corriente (80 % y 93 % de I0Set = 15 µA)
# Límite: las 10:00Z siguientes al arranque (07:00 Chile, antes de la luz de la mañana), o --deadline "AAAA-MM-DD HH:MM:SS" (UTC).
_now = time.gmtime()
_next10 = calendar.timegm((_now.tm_year, _now.tm_mon, _now.tm_mday, 10, 0, 0, 0, 0, 0))
if "--max-hours" in sys.argv:  # DAQ2-14 (2026-10-22): corte = arranque + N horas
    DEADLINE = time.time() + float(_arg("--max-hours", "6")) * 3600
elif "--deadline" in sys.argv:
    DEADLINE = calendar.timegm(time.strptime(_arg("--deadline", ""), "%Y-%m-%d %H:%M:%S"))
else:
    DEADLINE = _next10 if time.time() < _next10 else _next10 + 86400
DUR = {"stab": 120, "sw": 120, "ext": float(_arg("--ext-s", "900")), "check": 120, "self": 600, "diag": 30} if not SIM else \
    {"stab": 5, "sw": 15, "ext": 20, "check": 10, "self": 15, "diag": 5}
# Seguimiento 2026-10-31 (3): rampa lenta. En el escalón +0.25 V del 2026-10-30, la rampa a RUp = 1 V/s dio saltos de IMon de +2.6 a +7.9 µA
# en segundos en CH00-03, y CH03 (7.5 µA a nominal) pasó I0Set = 15 µA 5 s y se disparó. Si el salto es corriente de carga (I = C·dV/dt),
# baja 10 veces con RUp = 0.1 V/s. hvmon NO escribe RUp (sus únicas escrituras son Pw, ClrAlarm y V0Set de la sesión de barrido): RUp lo fija
# René antes del barrido (GECO o una escritura nueva en hvmon que no existe todavía). El barrido comprueba al empezar que el RUp leído de
# CH00-07 es el esperado (--rup, 0.1 V/s; 1.0 en simulación, que es el del fake) y, tras cada escritura de V0Set verificada, espera
# 3 × (escalón / RUp) antes de la siguiente (0.25 V a 0.1 V/s: 7.5 s), para que la corriente de rampa de un canal no se sume a la del siguiente.
RUP_EXPECTED_V_PER_S = float(_arg("--rup", "1.0" if SIM else "0.1"))
SETTLE_AFTER_WRITE_S = 3.0 * MAX_STEP_V / RUP_EXPECTED_V_PER_S
# Opcional (precaución, no es la causa demostrada del corte del enlace del 2026-10-30): tasa del SW trigger del run de calibración; la duración
# se escala para conservar el número de eventos (100 Hz × 120 s = 12 000).
SW_RATE_HZ = float(_arg("--sw-rate", "100"))
if not SIM:
    DUR["sw"] = 120.0 * 100.0 / SW_RATE_HZ
# duración estimada de un escalón completo: runs + análisis + rampas
STEP_EST_S = (DUR["stab"] + DUR["sw"] + DUR["ext"] + DUR["check"] + (0 if SELF_CHECK_ONLY else 4 * DUR["self"]) + 900) if not SIM else 60
CHECK_T_PE = 8        # run de control del auto-trigger: un par a 8 p.e.
DIAG_SINGLE_PE = 20   # diagnóstico, grupo solo: 20 p.e. (a 8 p.e. un grupo solo dispararía a ~1e5 Hz)
CHECK_MIN_FRAC = 0.10  # tasa medida < 10 % de la accidental prevista -> configuración que no funciona
CHECK_MIN_EXPECTED = 3.0  # menos eventos previstos en el run de control: no comprobable (DAQ2-16: antes 0/0 = "funciona")
GAIN_WINDOW = (0.5, 6.0)  # ventana de picos del ajuste de G en múltiplos de la G esperada (DAQ2-16)
DEAD = [3]                # ch3 muerto: fuera de la eficiencia X·Y
MASKED = [2, 9, 15]
PAIRS = [(0, 2), (0, 3), (1, 2), (1, 3)]
THRESHOLDS = [4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30]
W_S = 8e-9          # ventana de la mayoría: 1 reloj (0x810C[23:20] = 1)
FLUX_SQ_HZ = 10.89 / 60.0  # modelo: 1 cm^-2 min^-1 (PDG) × cuadrado activo de 10.89 cm²
ACC_GEOM = 0.03646  # fracción de triggers externos que cruza el área activa de hod2018 (MC cos²θ con proyección, DAQ2-15)
OUT = REPO + "/outputs/daq2-15" + ("-sim" if SIM else "") + "/" + CAMPAIGN
LOGD = "/data/daq/sim/logs" if SIM else "/data/daq/logs"
LOG, JL = f"{LOGD}/biasscan_{CAMPAIGN}.log", f"{LOGD}/biasscan_{CAMPAIGN}.jsonl"
CFGD = ("/data/daq/sim/config/" if SIM else "/data/daq/config/") + CAMPAIGN
RAW = ("/data/daq/sim/raw/" if SIM else "/data/daq/raw/") + CAMPAIGN
PHYS_CFG = "/data/daq/config/hodoscope_2026_10/physics_ext_B.json"
C1 = os.path.expanduser("~/.local/opt/daqv2/share/calibration/pe_calibration_2026_10_c1.json")
MAP = os.path.expanduser("~/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r4.json")
ENV = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser("~/root6/lib"), PATH=os.path.expanduser("~/root6/bin") + ":" + os.environ.get("PATH", ""))
ROOTINC = 'gSystem->AddIncludePath("-I%s/daq/core/include")' % CODE


class Abort(Exception):
    """Clase (c): riesgo real -> aborto, restauración y run de física."""


class Recoverable(Exception):
    """Clase (b): fallo recuperable agotados los reintentos -> se omite el escalón y se sigue."""


RETRIES = 3
RETRY_WAIT_S = 90 if not SIM else 1
RAMP_TAU_S = 115.0                      # descarga del PCB del detector (DAQ2-10)
STABLE_TIMEOUT_S = 6 * RAMP_TAU_S + 60  # clase (a): espera de fin de rampa en todos los canales
HAZARD_BITS = (1 << 3) | (1 << 6) | (1 << 9)  # sobrecorriente, trip externo, trip interno
RAMP_BITS = (1 << 1) | (1 << 2)


LOG_LOCK = threading.Lock()


def log(msg, **k):
    t = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    os.makedirs(LOGD, exist_ok=True)
    k.update(t_utc=t, msg=msg)
    with LOG_LOCK:  # también escribe el hilo de vigilancia de corriente
        with open(LOG, "a") as f:
            f.write(f"{t} {msg}\n")
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
    if not res.get("ok") and "result" not in res:
        # DAQ2-17: con ok = false hvctl no imprime el JSON; el resultado ("hazard", "rejected", ...) va en stderr como
        # "comando N: <resultado> (<detalle>)". Sin esto un riesgo de hvmon se trataba como un fallo cualquiera.
        import re
        m = re.search(r"comando \d+: (\w+) \((.*)\)\s*$", (r.stderr or "").strip(), re.S)
        if m:
            res.update(result=m.group(1), detail=m.group(2))
        elif r.stderr:
            res.setdefault("detail", r.stderr.strip()[-400:])
    return res


STATE_PATH = "/dev/shm/daqv2-sim-test/hvmon_state.json" if SIM else "/dev/shm/daqv2/hvmon_state.json"


def hvstate_full():
    return json.load(open(STATE_PATH))


def hvstate():
    return {c["ch"]: c["values"] for c in hvstate_full()["channels"]}


def check_hv(where):
    """Clase (c) si hay riesgo: conexión HV perdida o estado viejo, Pw = 0, bit on caído, trip o
    sobrecorriente. Los bits de rampa no son riesgo (clase a). Devuelve el estado por canal."""
    s = hvstate_full()
    conn = s.get("connection", {})
    age = time.time() - s.get("t_utc_ns", 0) / 1e9
    if not conn.get("connected", False) or conn.get("stale", True) or age > 30:
        raise Abort(f"{where}: conexión HV perdida (connected {conn.get('connected')}, stale {conn.get('stale')}, estado de hace {age:.0f} s)")
    st = {c["ch"]: c["values"] for c in s["channels"]}
    for ch in range(8):
        v = st[ch]
        word = int(v.get("Status", 0))
        if v.get("Pw") != 1 or not (word & 1) or (word & HAZARD_BITS):
            raise Abort(f"{where}: riesgo en CH{ch:02d}: Pw {v.get('Pw')} Status 0x{word:X}")
    imon_guard(st, where)
    return st


GUARD = {"over_since": {}, "stop_after_step": None, "imon_max_A": {}, "imon_max_step_A": {}, "finishing": False, "abort_reason": None}
GUARD_LOCK = threading.Lock()


def imon_guard(st, where):
    with GUARD_LOCK:
        _imon_guard(st, where)


def _imon_monitor(stop):
    """Hilo de vigilancia de corriente (DAQ2-17): cada IMON_POLL_S lee el estado de hvmon y aplica la guarda; ante
    IMon ≥ 14 µA manda SIGUSR1 al proceso (el hilo principal levanta Abort y restaura) y termina."""
    while not stop.wait(IMON_POLL_S):
        if GUARD["finishing"]:
            return
        try:
            st = {c["ch"]: c["values"] for c in hvstate_full()["channels"]}
        except Exception:  # noqa: BLE001 — estado ilegible un instante: lo trata check_hv en el hilo principal
            continue
        try:
            imon_guard(st, "vigilancia de corriente (hilo)")
        except Abort as e:
            GUARD["abort_reason"] = str(e)
            log(f"  {e}: aviso al hilo principal")
            os.kill(os.getpid(), signal.SIGUSR1)
            return


def _sigusr1(signum, frame):
    if GUARD["finishing"]:
        return  # ya se está restaurando
    raise Abort(GUARD["abort_reason"] or "guarda de corriente (hilo)")


def _imon_guard(st, where):
    """DAQ2-17: IMon ≥ 14 µA -> Abort (restauración inmediata); > 12 µA sostenido 30 s -> no se pasa al escalón siguiente
    (el escalón actual sigue mientras IMon < 14 µA). Registra el máximo por grupo."""
    now = time.time()
    for ch in range(8):
        i = st[ch].get("IMon_A")
        if i is None:
            continue
        GUARD["imon_max_A"][ch] = max(GUARD["imon_max_A"].get(ch, 0.0), i)
        GUARD["imon_max_step_A"][ch] = max(GUARD["imon_max_step_A"].get(ch, 0.0), i)
        if i >= IMON_STOP_A:
            raise Abort(f"{where}: IMon CH{ch:02d} = {i * 1e6:.2f} µA ≥ {IMON_STOP_A * 1e6:.0f} µA: restauración inmediata")
        if i > IMON_WARN_A:
            t0 = GUARD["over_since"].setdefault(ch, now)
            if now - t0 >= IMON_WARN_S and not GUARD["stop_after_step"]:
                GUARD["stop_after_step"] = f"IMon CH{ch:02d} = {i * 1e6:.2f} µA > {IMON_WARN_A * 1e6:.0f} µA durante {now - t0:.0f} s"
                log(f"  GUARDA DE CORRIENTE: {GUARD['stop_after_step']}: se termina este escalón y no se sigue")
        else:
            GUARD["over_since"].pop(ch, None)


def board_temp(step, where):
    """Temperatura de la placa A7030DP (única referencia de temperatura; sigue la del rack, 1 °C de resolución)."""
    try:
        t = hvstate_full().get("board_temperature_degC")
    except Exception:  # noqa: BLE001
        t = None
    step["board_temp"].append({"t_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "where": where, "degC": t})
    return t


def wait_stable(where, timeout=STABLE_TIMEOUT_S):
    """Clase (a): espera a que ningún canal de CH00-07 esté en rampa (Status 0x1); riesgos -> Abort;
    timeout -> Recoverable."""
    t0 = time.time()
    while True:
        st = check_hv(where)
        ramping = [ch for ch in range(8) if int(st[ch].get("Status", 0)) & RAMP_BITS]
        if not ramping:
            if time.time() - t0 > 5:
                log(f"  {where}: rampas terminadas tras {time.time() - t0:.0f} s")
            return st
        if time.time() - t0 > timeout:
            raise Recoverable(f"{where}: canales {ramping} siguen en rampa tras {timeout:.0f} s")
        time.sleep(2)


def set_bias(delta):
    """Lleva CH00-07 a nominal + delta en escalones de ≤ MAX_STEP_V (todos los canales por escalón)."""
    st = check_hv("antes de cambiar V0Set")
    cur = [st[ch]["V0Set_V"] for ch in range(8)]
    target = [round(NOMINAL[ch] + delta, 3) for ch in range(8)]
    while any(abs(cur[ch] - target[ch]) > 0.004 for ch in range(8)):
        for ch in range(8):
            if abs(cur[ch] - target[ch]) <= 0.004:
                continue
            step = max(-MAX_STEP_V, min(MAX_STEP_V, target[ch] - cur[ch]))
            v = round(cur[ch] + step, 3)
            for attempt in range(1, RETRIES + 1):
                res = hvctl("bias", "set", "--ch", str(ch), "--value", f"{v:.3f}", "--authorize", f"{AUTH}: CH{ch:02d} V0Set {v:.2f} V")
                log(f"V0Set CH{ch:02d} -> {v:.3f} V (intento {attempt}): {res.get('result')} ({res.get('detail', '')})", ev="v0set", ch=ch,
                    value=v, res=res)
                if res.get("ok") and res.get("result") in ("verified", "noop"):
                    if res.get("result") == "verified":
                        log(f"  espera de estabilización tras la escritura: {SETTLE_AFTER_WRITE_S:.1f} s (3 × {abs(step):.2f} V / {RUP_EXPECTED_V_PER_S:g} V/s)")
                        time.sleep(3.0 * abs(step) / RUP_EXPECTED_V_PER_S)
                    break
                if res.get("result") == "hazard":
                    raise Abort(f"V0Set CH{ch:02d} = {v:.3f} V: {res.get('detail')}")
                check_hv(f"tras el fallo de V0Set CH{ch:02d}")  # un riesgo detectado aquí también aborta
                check_offset_drift(f"tras el fallo de V0Set CH{ch:02d}")  # DAQ2-17: un salto del desvío > 0.15 V aborta
                if attempt == RETRIES:
                    raise Recoverable(f"V0Set CH{ch:02d} = {v:.3f} V sin verificar tras {RETRIES} intentos: {res.get('result')}: {res.get('detail')}")
                time.sleep(RETRY_WAIT_S)
            cur[ch] = v
    st = wait_stable("después de cambiar V0Set")
    return {ch: {"V0Set_V": st[ch]["V0Set_V"], "VMon_V": st[ch]["VMon_V"], "IMon_A": st[ch]["IMon_A"]} for ch in range(8)}


def check_offset_drift(where, st=None):
    """DAQ2-17 (René): Abort si |VMon − V0Set − desvío al nominal| > OFFSET_DRIFT_MAX_V en algún canal de CH00-07."""
    sess = hvctl("bias", "status").get("session", {})
    off0 = {int(c): v for c, v in (sess.get("vmon_offset_nominal_V") or sess.get("vmon_offset_V") or {}).items()}
    st = st or check_hv(where)
    drift = {ch: round(st[ch]["VMon_V"] - st[ch]["V0Set_V"] - off0.get(ch, 0.0), 3) for ch in range(8)}
    far = {ch: d for ch, d in drift.items() if abs(d) > OFFSET_DRIFT_MAX_V}
    if far:
        raise Abort(f"{where}: desvío VMon − V0Set a más de {OFFSET_DRIFT_MAX_V} V del desvío al nominal en {far}")
    return drift


def check_vmon(where, tries=3):
    """DAQ2-17: |VMon − V0Set − desvío| ≤ VMON_TOL_V en CH00-07, con el desvío del último escalón verificado por hvmon
    (hvctl bias status). Hasta `tries` lecturas separadas por 60 s; si no, recuperable (se omite el escalón)."""
    for k in range(tries):
        sess = hvctl("bias", "status").get("session", {})
        off = {int(c): v for c, v in (sess.get("vmon_offset_V") or {}).items()}
        st = check_hv(where)
        check_offset_drift(where, st)  # DAQ2-17 (René): riesgo, no recuperable
        dev = {ch: st[ch]["VMon_V"] - st[ch]["V0Set_V"] - off.get(ch, 0.0) for ch in range(8)}
        bad = {ch: round(d, 3) for ch, d in dev.items() if abs(d) > VMON_TOL_V}
        if not bad:
            return {ch: round(d, 4) for ch, d in dev.items()}
        log(f"  {where}: VMon fuera de {VMON_TOL_V} V en {bad} (lectura {k + 1} de {tries})")
        if k + 1 < tries:
            time.sleep(60 if not SIM else 1)
    raise Recoverable(f"{where}: VMon fuera de {VMON_TOL_V} V del V0Set + desvío en {bad}")


def write_cfg(name, base, mutate):
    c = json.load(open(base))
    c["campaign"]["name"] = CAMPAIGN
    c["campaign"]["operator"] = "claude (DAQ2-15, barridos de tensión de día y de noche, autorizado por René)"
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
    for attempt in range(1, RETRIES + 1):
        try:
            return _daqrun_once(cfg, seconds)
        except (OpenFailed, RunFailed) as e:  # clase (b)
            kind = "apertura del digitizador" if isinstance(e, OpenFailed) else "run"
            log(f"  {kind} falló ({os.path.basename(cfg)}, intento {attempt} de {RETRIES}): {e}")
            if attempt == RETRIES:
                raise Recoverable(f"{kind} {os.path.basename(cfg)} falló {RETRIES} veces: {e}")
            time.sleep(RETRY_WAIT_S)


class RunFailed(Exception):
    pass


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
        try:
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
                    _last_close[0] = time.time()
                    raise RunFailed(f"daqrun {os.path.basename(cfg)} no terminó a tiempo")
        except BaseException:  # DAQ2-17: Abort por señal (hilo de corriente, SIGTERM) con el run abierto: cierre limpio del run
            if p.poll() is None:
                p.terminate()
                try:
                    p.wait(timeout=300)
                except subprocess.TimeoutExpired:
                    p.kill()
                    p.wait()
            _last_close[0] = time.time()
            raise
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
        raise RunFailed(f"rc {r.returncode}: {next((l for l in txt.splitlines() if 'ERROR' in l), txt[-300:])[-300:]}")
    segs = sorted(glob.glob(f"{RAW}/{base}_seg*.root"))
    idx = json.load(open(f"{RAW}/{base}.index.json"))
    return {"base": base, "segments": segs, "events": idx.get("n_events", 0), "wall_s": time.time() - t0,
            "trigger_signature": idx.get("trigger_signature", "")}


def root(macro, *args, d="validation"):
    a = ",".join(json.dumps(x) if isinstance(x, str) else str(x).lower() if isinstance(x, bool) else str(x) for x in args)
    sh(["root", "-l", "-b", "-q", "-e", ROOTINC, f"{CODE}/daq/tools/{d}/{macro}.C+({a})"], timeout=3000, heavy=True)


# ---------------------------------------------------------------- análisis de cada escalón

def analyse_sw(run, step, seedsCal):
    pat = RAW + "/" + run["base"] + "_seg*.root"
    seeds = f"{OUT}/seeds_step{step}.json"
    cal0 = json.load(open(seedsCal))
    # G esperada por canal (DAQ2-16): la medida en el escalón anterior (o c1); si no se midió, la esperada que traía.
    expected = {}
    for c in cal0["channels"]:
        if c["ch"] >= 32:
            continue
        measured = c.get("G_mV") and c.get("status") in ("ok", "changed") and "G_mV_fallback" not in c
        g = c["G_mV"] if measured else c.get("G_expected_mV") or c.get("seed_mV") or c.get("G_mV")
        if g:
            expected[c["ch"]] = g
    json.dump([{"ch": ch, "status": "ok", "value": g} for ch, g in expected.items()], open(seeds, "w"))
    cal = f"{OUT}/pecal_step{step}.json"
    root("commissioning_pe_calibration", pat, f"daq2-16-step{step}", cal, seeds, GAIN_WINDOW[0], GAIN_WINDOW[1])
    noise = f"{OUT}/noise_step{step}.csv"
    root("selftrig_noise", pat, "", MAP, noise, 64)
    c = json.load(open(cal))
    G = {e["ch"]: e.get("G_mV") for e in c["channels"]}
    ok = [G[ch] for ch in range(32) if G.get(ch) and c["channels"][ch].get("status") in ("ok", "changed")]
    med = statistics.median(ok) if ok else statistics.median([e["G_mV"] for e in cal0["channels"] if e.get("G_mV")])
    # calibración del escalón para daqrun: canales sin G aceptada -> su G esperada (o la mediana), marcados
    for e in c["channels"]:
        if e["ch"] < 32:
            e["G_expected_mV"] = expected.get(e["ch"], med)
        if e["ch"] < 32 and not (e.get("G_mV") and e.get("status") in ("ok", "changed")):
            e["G_mV_fallback"] = expected.get(e["ch"], med)
            e["G_mV"] = e["G_mV_fallback"]
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
            for T in [0.5, 1.5, 3] + THRESHOLDS:
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
    """Eficiencia X·Y (DAQ2-16): impactos con la definición única (hitdef_extract.C) y el estimador compartido
    xy_estimator.xy_fraction; "all" sin ch3 (muerto), "masked" sin ch2, ch3, ch9, ch15. eff_* = q / ACC_GEOM."""
    pat = RAW + "/" + run["base"] + "_seg*.root"
    G = [0.0] * 32
    for e in json.load(open(cal))["channels"]:
        if e["ch"] < 32 and e.get("G_mV"):
            G[e["ch"]] = e["G_mV"]
    med = statistics.median(g for g in G if g > 0)
    G = [g if g > 0 else med for g in G]
    pref = f"{OUT}/hits_step{step}"
    root("hitdef_extract", pat, "daqv2", ",".join(f"{g:.5f}" for g in G), pref, d="characterization")
    ev = xy.load_hits(pref + "_hits.bin")
    keep = ("N", "N_XY", "acc_ctl", "acc_mix", "N_X", "N_Y", "N_Xc", "N_Yc", "aX", "aY", "rX", "rY", "f_mix", "f_mix_err", "f_ctl", "f_ctl_err",
            "e_X|Y_mix", "e_Y|X_mix")
    res = {}
    for tag, excl in (("all", DEAD), ("masked", sorted(set(MASKED + DEAD)))):
        out = {"N": len(ev), "excluded_channels": excl}
        for T in THRESHOLDS:
            Xs, Xc = xy.plane_hits(ev, G, T, [c for c in range(16) if c not in excl])
            Ys, Yc = xy.plane_hits(ev, G, T, [c for c in range(16, 32) if c not in excl])
            r = xy.xy_fraction(Xs, Ys, Xc, Yc)
            d = {k: r.get(k) for k in keep}
            d.update(frac_raw=r["N_XY"] / r["N"] if r["N"] else 0, eff_mix=r.get("f_mix", 0) / ACC_GEOM, eff_ctl=r.get("f_ctl", 0) / ACC_GEOM)
            pairs = {}
            for gx, gy in PAIRS:
                a, ca = xy.plane_hits(ev, G, T, [c for c in range(8 * gx, 8 * gx + 8) if c not in excl])
                b, cb = xy.plane_hits(ev, G, T, [c for c in range(8 * gy, 8 * gy + 8) if c not in excl])
                p = xy.xy_fraction(a, b, ca, cb)
                pairs[f"G{gx}G{gy}"] = {"N": p["N_XY"], "acc_ctl": p["acc_ctl"], "f_mix": p.get("f_mix"), "eff": p.get("f_mix", 0) / ACC_GEOM}
            d["pairs"] = pairs
            out[T] = d
        if tag == "all":
            Gx = [G[c] for c in range(16) if c not in excl]
            xs_ = [c for c in range(16) if c not in excl]
            ys5 = (ev["aS"][:, 16:32] / [G[c] for c in range(16, 32)] > 5).any(1)
            xs = sorted((ev["aS"][ys5][:, xs_] / Gx).max(1).tolist())
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
    # Muones/h = exceso sobre los accidentales previstos, acotado a ≥ 0; con N = 0, límite superior al 90 % (2.30 cuentas).
    return {"pair": f"G{gx}G{gy}", "events": n, "rate_Hz": n / dur, "window": [lo, hi], "N_XY_in_time": nxy, "acc_pred": acc_pred,
            "acc_offtime": acc_ctl, "ratio_pred": nxy / acc_pred if acc_pred > 0 else float("inf"),
            "muons_per_h": max(0.0, nxy - acc_pred) / dur * 3600, "muons_per_h_raw": (nxy - acc_pred) / dur * 3600,
            "muons_per_h_ul90": (max(nxy, 0) + 2.30) / dur * 3600 if nxy <= 3 else None, "trigger_signature": run["trigger_signature"],
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


def physics_alive():
    """El run de física está vivo: unidad activa, página del run respondiendo, sin errores de lectura y con
    triggers externos en los últimos 10 min (tasa > 0.1 Hz). Devuelve (vivo, detalle)."""
    if SIM:
        return True, "(sim)"
    st = physics_status()
    if not st["active"]:
        return False, f"unidad daq-physics {st['state']}; última línea: {st['last']}"
    try:
        import urllib.request
        r = urllib.request.urlopen("http://127.0.0.1:8080/status/state_json/root.json", timeout=10).read()
        s = json.loads(json.loads(r)["fTitle"])
        d = s["dash"]
        detail = (f"run {d.get('run')}, {d.get('events')} eventos, tasa 10 min {d.get('rate_10min_Hz')} Hz, último trigger hace "
                  f"{d.get('since_last_trigger_s')} s, errores de lectura {s['readout'].get('read_errors')}")
        ok = (s["readout"].get("read_errors", 1) == 0 and not s["readout"].get("finished", False) and (d.get("rate_10min_Hz") or 0) > 0.1)
        return ok, detail
    except Exception as e:  # noqa: BLE001
        return False, f"la página del run no responde: {e}"


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
    if sh(["systemctl", "--user", "is-active", "daq-physics"], check=False).stdout.strip() == "active":
        return "ya estaba activo (no se reinició)"
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
    log("=== DAQ2-15/17 barrido de tensión" + (" (SIMULACIÓN)" if SIM else "") + f"; campaña {CAMPAIGN}; límite "
        + time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(DEADLINE)))
    log(f"código: {CODE} (commit {CODE_COMMIT})", code_root=CODE, code_commit=CODE_COMMIT)
    results = {"steps": [], "aborted": None, "skipped": [], "skipped_steps": [], "current_step": None, "code_root": CODE, "code_commit": CODE_COMMIT}
    phys_stopped = False
    try:
        st = hvstate()
        read = [st[ch]["V0Set_V"] for ch in range(8)]
        log(f"V0Set leído: {read}")
        if any(abs(a - b) > 0.005 for a, b in zip(read, NOMINAL)):
            raise Abort(f"V0Set leído {read} distinto del nominal esperado {NOMINAL}: no se empieza")
        rup = [hvstate()[ch].get("RUp_V_per_s") for ch in range(8)]
        log(f"RUp leído: {rup} V/s (esperado {RUP_EXPECTED_V_PER_S:g} V/s; espera tras cada escritura {SETTLE_AFTER_WRITE_S:.1f} s)")
        if any(r is None or abs(r - RUP_EXPECTED_V_PER_S) > 0.005 for r in rup):
            raise Abort(f"RUp leído {rup} distinto del esperado {RUP_EXPECTED_V_PER_S:g} V/s en CH00-07: no se empieza (hvmon no escribe RUp; "
                        f"fijarlo antes o pasar --rup con el valor real)")
        check_hv("inicio")
        alive, detail = physics_alive()
        log(f"run de física antes de empezar: {'VIVO' if alive else 'NO VIVO'} ({detail})")
        results["physics_check"] = {"alive": alive, "detail": detail}
        if not alive:
            raise Abort("el run de física no está vivo (sin triggers externos o con errores): no se empieza el barrido: " + detail)
        # Campaña y configuraciones base
        camp = ["daqcampaign"] + (["--sim"] if SIM else ["--authorize", AUTH])
        sh([(REPO + "/daq/build/apps/daqcampaign/" if SIM else BIN + "/") + "daqcampaign"] + camp[1:] + ["open", CAMPAIGN], check=False)
        swcfg = write_cfg("sw", "/data/daq/config/commissioning3_2026_10_21/toggle2_sw100_full.json",
                          lambda c: c["digitizer"].__setitem__("sw_trigger_rate_Hz", SW_RATE_HZ))
        extcfg = write_cfg("ext", PHYS_CFG, lambda c: c["run"].__setitem__("max_duration_s", 0))
        results["physics_before"] = physics_stop()
        phys_stopped = True
        res = hvctl("bias", "begin", "--expected", ",".join(f"{v:.2f}" for v in NOMINAL), "--max-above", f"{MAX_ABOVE_V:g}",
                    "--max-step", f"{MAX_STEP_V:g}", "--authorize", AUTH)
        log(f"bias begin: {res.get('result')} ({res.get('detail')})", res=res)
        if not res.get("ok"):
            raise Abort("bias begin: " + str(res.get("detail")))
        signal.signal(signal.SIGUSR1, _sigusr1)
        monitor = threading.Thread(target=_imon_monitor, args=(MONITOR_STOP,), name="imon", daemon=True)
        monitor.start()
        log(f"vigilancia de corriente en un hilo cada {IMON_POLL_S:g} s (> {IMON_WARN_A * 1e6:.0f} µA {IMON_WARN_S:g} s: no se sigue; "
            f"≥ {IMON_STOP_A * 1e6:.0f} µA: restauración inmediata)")
        prevCal = C1
        for i, dv in enumerate(STEPS):
            if time.time() + STEP_EST_S > DEADLINE:
                results["skipped"] = STEPS[i:]
                log(f"escalones {STEPS[i:]} omitidos: terminarían después del límite")
                break
            step = {"delta_V": dv, "t_start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "t_start_chile": time.strftime("%Y-%m-%d %H:%M:%S %z", time.localtime()), "board_temp": []}
            log(f"--- escalón ΔV = {dv:+.2f} V")
            results["current_step"] = dv
            try:
                prevCal = run_step(step, i, dv, prevCal, swcfg, extcfg)
            except Recoverable as e:  # clase (b): se omite el escalón, se guarda lo tomado y se sigue
                step["skipped_reason"] = str(e)
                step["t_end_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                results["skipped_steps"].append({"delta_V": dv, "reason": str(e), "partial": step})
                log(f"ESCALÓN {dv:+.2f} V OMITIDO (recuperable): {e}")
                json.dump(results, open(f"{OUT}/biasscan_results.json", "w"), indent=1, ensure_ascii=False, default=str)
                continue
            results["steps"].append(step)
            json.dump(results, open(f"{OUT}/biasscan_results.json", "w"), indent=1, ensure_ascii=False, default=str)
            if GUARD["stop_after_step"]:
                results["skipped"] = STEPS[i + 1:]
                results["current_guard"] = GUARD["stop_after_step"]
                log(f"escalones {STEPS[i + 1:]} no empezados por la guarda de corriente: {GUARD['stop_after_step']}")
                break
        results["current_step"] = None
    except Exception as e:  # noqa: BLE001 — clase (c) o error inesperado: aborto con restauración
        results["aborted"] = str(e)
        results["aborted_step"] = results.get("current_step")
        log(f"ABORTO: {e}")
    finally:
        finish(results)


def guarded(fn, what):
    """Análisis: un fallo es recuperable (clase b): se omite el resto del escalón con su motivo."""
    try:
        return fn()
    except (Abort, Recoverable):
        raise
    except Exception as e:  # noqa: BLE001
        raise Recoverable(f"{what} falló: {e}")


def run_step(step, i, dv, prevCal, swcfg, extcfg):
    """Un escalón completo; devuelve la calibración del escalón. Recoverable / Abort según la política."""
    GUARD["imon_max_step_A"] = {}
    step["hv_after_set"] = set_bias(dv)
    time.sleep(DUR["stab"])
    step["hv_before_sw"] = {ch: v for ch, v in wait_stable("antes del run SW").items() if ch < 8}
    step["vmon_deviation_V"] = check_vmon("antes del run SW")
    log(f"  temperatura de la placa: {board_temp(step, 'antes del run SW')} °C")
    sw = daqrun(swcfg, DUR["sw"])
    step["sw_run"] = sw
    check_hv("después del run SW")
    step["sw"] = guarded(lambda: analyse_sw(sw, i, prevCal), "análisis SW")
    prevCal = step["sw"]["calibration"]
    log(f"  SW: {sw['events']} eventos; G mediana {step['sw']['G_median_mV']:.3f} mV ({step['sw']['n_G_ok']} canales con G)")
    check_hv("antes del run externo")
    board_temp(step, "antes del run externo")
    ext = daqrun(write_cfg(f"ext_step{i}", extcfg, lambda c, cal=prevCal: c["calibration"].__setitem__("file", cal)), DUR["ext"])
    step["ext_run"] = ext
    step["ext"] = guarded(lambda: analyse_ext(ext, i, prevCal), "análisis externo")
    e10 = step["ext"]["all"].get(10, {})
    log(f"  externo: {ext['events']} eventos; ε_XY(10 p.e.) ctl {e10.get('eff_ctl', 0):.3f}, mezclados {e10.get('eff_mix', 0):.3f}")
    table, T, why = guarded(lambda: predict(step["sw"], step["ext"]), "predicción")
    step["prediction"] = table
    step["T_chosen_pe"] = T
    step["T_choice"] = why
    log(f"  umbral elegido {T} p.e. ({why})")
    step["self"] = []
    board_temp(step, "antes del auto-trigger")
    chk = step["selftrigger_check"] = selftrigger_check(step, i, prevCal, extcfg)
    if TRIGGER_DIAG and dv == 0:
        step["trigger_diag"] = trigger_diag(step, i, prevCal, extcfg)
    step["imon_max_A"] = dict(GUARD["imon_max_step_A"])
    if SELF_CHECK_ONLY or not chk.get("working"):
        log(f"  runs de auto-trigger omitidos: " + ("--self-check-only (DAQ2-16)" if SELF_CHECK_ONLY else
            f"{'NO COMPROBABLE' if chk.get('working') is None else 'NO FUNCIONA'} ({chk.get('reason')})"))
        step["self_muons_per_h_sum"] = None
        step["t_end_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        return prevCal
    for gx, gy in PAIRS:
        check_hv(f"antes del auto-trigger G{gx}∧G{gy}")
        def mut(c, gx=gx, gy=gy, cal=prevCal):
            c["calibration"]["file"] = cal
            c["digitizer"]["trigger"] = {"mode": "self_hod2018", "self": {"groups": [gx, gy], "threshold_pe": T, "masked_channels": MASKED,
                                                                          "majority_groups": 2, "window_clk": 1, "include_external": False,
                                                                          "target_rate_Hz": 50}}
        run = daqrun(write_cfg(f"self_step{i}_G{gx}G{gy}", extcfg, mut), DUR["self"])
        a = guarded(lambda: analyse_self(run, i, T, gx, gy, prevCal, step["sw"]), f"análisis G{gx}∧G{gy}")
        step["self"].append(a)
        log(f"  auto-trigger G{gx}∧G{gy} T={T}: {run['events']} eventos ({a['rate_Hz']:.2f} Hz), X·Y en tiempo {a['N_XY_in_time']}, "
            f"accidentales previstos {a['acc_pred']:.1f}, razón {a['ratio_pred']:.2f}, muones/h {a['muons_per_h']:.1f}"
            + (f" (límite 90 % {a['muons_per_h_ul90']:.1f})" if a.get("muons_per_h_ul90") is not None else ""))
    step["self_muons_per_h_sum"] = sum(a["muons_per_h"] for a in step["self"])
    step["t_end_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    return prevCal


def selftrigger_check(step, i, cal, extcfg):
    """Run de control (DUR["check"]) con el par de mayor tasa accidental prevista a CHECK_T_PE p.e.: tasa medida
    contra la prevista 2W·r_X·r_Y (tasas de ruido del run SW del escalón, ch2/9/15 fuera). Un fallo del run
    (recuperable) también marca el auto-trigger como no comprobado, sin perder el resto del escalón."""
    r = step["sw"]["group_rates_masked_Hz"][CHECK_T_PE]
    gx, gy = max(PAIRS, key=lambda p: r[p[0]] * r[p[1]])
    pred = 2 * W_S * r[gx] * r[gy]
    def mut(c):
        c["calibration"]["file"] = cal
        c["digitizer"]["trigger"] = {"mode": "self_hod2018", "self": {"groups": [gx, gy], "threshold_pe": CHECK_T_PE, "masked_channels": MASKED,
                                                                      "majority_groups": 2, "window_clk": 1, "include_external": False,
                                                                      "target_rate_Hz": max(50.0, 2 * pred)}}
    out = {"pair": f"G{gx}G{gy}", "T_pe": CHECK_T_PE, "pred_acc_Hz": pred, "duration_s": DUR["check"], "expected_events": pred * DUR["check"]}
    check_hv("antes del run de control del auto-trigger")
    try:
        run = daqrun(write_cfg(f"check_step{i}_G{gx}G{gy}", extcfg, mut), DUR["check"])
    except Recoverable as e:
        out.update(working=False, reason=f"run de control falló: {e}")
        return out
    rate = run["events"] / DUR["check"]
    frac = rate / pred if pred > 0 else float("nan")
    out.update(events=run["events"], rate_Hz=rate, frac_of_pred=frac, base=run["base"], trigger_signature=run["trigger_signature"])
    if out["expected_events"] < CHECK_MIN_EXPECTED:
        out["working"] = None
        out["reason"] = f"no comprobable: la predicción da {out['expected_events']:.2f} eventos en {DUR['check']} s (< {CHECK_MIN_EXPECTED:.0f})"
    else:
        out["working"] = frac >= CHECK_MIN_FRAC
        if not out["working"]:
            out["reason"] = f"tasa medida {rate:.3f} Hz = {frac:.3g} de la prevista {pred:.3f} Hz (< {CHECK_MIN_FRAC:.0%})"
    verdict = "FUNCIONA" if out["working"] else "NO COMPROBABLE" if out["working"] is None else "NO FUNCIONA"
    log(f"  control del auto-trigger G{gx}∧G{gy} a {CHECK_T_PE} p.e.: {run['events']} eventos en {DUR['check']} s = {rate:.3f} Hz; "
        f"accidentales previstos {pred:.3f} Hz ({out['expected_events']:.1f} eventos); fracción {frac:.3g}: {verdict}")
    out["registers"] = guarded(lambda: trigger_registers(run, i), "registros del run de control")
    return out


def trigger_diag(step, i, cal, extcfg):
    """DAQ2-17: runs de DUR["diag"] s con auto-trigger: cada grupo de hod2018 solo (mayoría 1, a DIAG_SINGLE_PE p.e.) y pares
    G0∧G1, G2∧G3 (ventana 1 reloj) y G0∧G2 (ventanas 1, 2, 4, 8 relojes) a CHECK_T_PE p.e.; tasa medida, predicción
    (grupo solo: r_g; par: 2·W·r_X·r_Y) y registros releídos de cada run. Un run que falla se anota y se sigue."""
    rates = step["sw"]["group_rates_masked_Hz"]
    plan = [((g,), 1, 1, DIAG_SINGLE_PE) for g in range(4)] + [((0, 1), 2, 1, CHECK_T_PE), ((2, 3), 2, 1, CHECK_T_PE)] + \
        [((0, 2), 2, w, CHECK_T_PE) for w in (1, 2, 4, 8)]
    out = []
    for groups, maj, w, T in plan:
        r = rates.get(T) or rates.get(str(T))
        pred = r[groups[0]] if len(groups) == 1 else 2 * w * W_S * r[groups[0]] * r[groups[1]]
        name = "diag_" + "G".join([""] + [str(g) for g in groups]) + f"_w{w}_T{T:g}"
        def mut(c, groups=groups, maj=maj, w=w, T=T, pred=pred):
            c["calibration"]["file"] = cal
            c["digitizer"]["trigger"] = {"mode": "self_hod2018", "self": {"groups": list(groups), "threshold_pe": T, "masked_channels": MASKED,
                                                                          "majority_groups": maj, "window_clk": w, "include_external": False,
                                                                          "target_rate_Hz": max(50.0, 2 * pred)}}
        d = {"groups": list(groups), "majority_groups": maj, "window_clk": w, "T_pe": T, "pred_Hz": pred, "duration_s": DUR["diag"]}
        check_hv(f"antes del diagnóstico {name}")
        try:
            run = daqrun(write_cfg(f"{name}_step{i}", extcfg, mut), DUR["diag"])
            d.update(events=run["events"], rate_Hz=run["events"] / DUR["diag"], base=run["base"], trigger_signature=run["trigger_signature"])
            d["registers"] = guarded(lambda: trigger_registers(run, i, tag=name), f"registros de {name}")
        except Recoverable as e:
            d["error"] = str(e)
        log(f"  diagnóstico {name}: {d.get('events', '–')} eventos en {DUR['diag']:g} s = {d.get('rate_Hz', float('nan')):.3f} Hz; "
            f"predicción {pred:.3g} Hz" + (f"; ERROR {d['error']}" if "error" in d else ""))
        out.append(d)
    return out


TRIG_KEYS = ("global_trigger_mask_0x810C", "group_trigger_threshold_adc", "channel_group_mask", "group_self_trigger", "trigger_polarity",
             "majority_level", "majority_window_clk", "group_dc_offset", "group_enable_mask", "ext_trigger", "sw_trigger")


def trigger_registers(run, i, tag="check"):
    """Registros de trigger escritos y releídos por daqrun (run_info_json del primer segmento) junto a lo pedido en la
    configuración (umbral en p.e., conversión a ADC, DC offset); cada uno al registro."""
    if not run["segments"]:
        return {"error": "sin segmentos"}
    f = f"{OUT}/{tag}_step{i}_run_info.json"
    sh(["root", "-l", "-b", "-q", f"{CODE}/daq/tools/biasscan/run_info_dump.C({json.dumps(run['segments'][0])}, {json.dumps(f)})"], timeout=300)
    d = json.load(open(f))
    dig = (d.get("run_info") or {}).get("digitizer", {})
    wr, rd = dig.get("config_written", {}), dig.get("config_read", {})
    asked = ((d.get("config") or {}).get("digitizer") or {})
    selfc = asked.get("trigger", {}).get("self", {})
    out = {"written": {k: wr.get(k) for k in TRIG_KEYS}, "read": {k: rd.get(k) for k in TRIG_KEYS},
           "group_registers_read": rd.get("group_registers"), "group_configuration_0x8000_read": rd.get("group_configuration_0x8000"),
           "group_enable_mask_0x8120_read": rd.get("group_enable_mask_0x8120"),
           "asked": {"groups": selfc.get("groups"), "threshold_pe": selfc.get("threshold_pe"), "threshold_adc": selfc.get("threshold_adc"),
                     "majority_groups": selfc.get("majority_groups"), "majority_level_register": selfc.get("majority_level_register"),
                     "window_clk": selfc.get("window_clk"), "channel_masks": selfc.get("channel_masks"), "group_dc_offset": asked.get("group_dc_offset"),
                     "conversion": {g: {k: v for k, v in c.items() if k in ("threshold_adc", "baseline_median_adc", "G_median_mV", "channel_mask")}
                                    for g, c in (selfc.get("conversion") or {}).items() if isinstance(c, dict)}}}
    out["mismatch"] = [k for k in TRIG_KEYS if k not in ("group_self_trigger", "trigger_polarity") and wr.get(k) != rd.get(k)]
    for k in TRIG_KEYS:
        log(f"    registro {k}: escrito {json.dumps(wr.get(k), ensure_ascii=False)}; releído {json.dumps(rd.get(k), ensure_ascii=False)}")
    log(f"    registros por grupo releídos {json.dumps(rd.get('group_registers'))}; 0x8000 {rd.get('group_configuration_0x8000')}; "
        f"0x8120 {rd.get('group_enable_mask_0x8120')}")
    log(f"    pedido: {json.dumps(out['asked'], ensure_ascii=False)}")
    log(f"    distintos entre escrito y releído: {out['mismatch'] or 'ninguno'}")
    return out


MONITOR_STOP = threading.Event()


def finish(results):
    """Siempre: nominal restaurado, run de física rearrancado (verificado) y correo con escalones hechos,
    omitidos (con motivo) y abortados."""
    GUARD["finishing"] = True  # el hilo de corriente ya no interrumpe la restauración
    MONITOR_STOP.set()
    res = hvctl("bias", "restore", "--authorize", AUTH + ": restauración del nominal")
    log(f"bias restore: {res.get('result')} ({res.get('detail')})", res=res)
    for k in range(2, 5):  # DAQ2-17: ningún canal puede quedar fuera del nominal; se reintenta
        if res.get("ok"):
            break
        time.sleep(30 if not SIM else 1)
        res = hvctl("bias", "restore", "--authorize", AUTH + f": restauración del nominal (intento {k})")
        log(f"bias restore (intento {k}): {res.get('result')} ({res.get('detail')})", res=res)
    results["restore"] = res
    results["imon_max_A"] = dict(GUARD["imon_max_A"])
    try:
        time.sleep(5)  # el archivo de estado de hvmon se publica en el ciclo siguiente
        st = hvstate()
        results["hv_final"] = {ch: {"V0Set_V": st[ch]["V0Set_V"], "VMon_V": st[ch]["VMon_V"], "Pw": st[ch]["Pw"], "Status": st[ch]["Status"]} for ch in range(8)}
        bad = [ch for ch in range(8) if abs(st[ch]["V0Set_V"] - NOMINAL[ch]) > 0.005 or st[ch]["Pw"] != 1]
        log(f"HV final: {'nominal y Pw = 1 en CH00-07' if not bad else 'NO nominal en ' + str(bad)}", hv=results["hv_final"])
    except Exception as e:  # noqa: BLE001
        log(f"no se pudo leer el estado final de HV: {e}")
    try:
        results["physics"] = physics_start()
    except Exception as e:  # noqa: BLE001
        results["physics"] = f"NO ARRANCÓ: {e}"
    log(f"run de física al final: {results['physics']}")
    results["physics_ok"] = bool(results["physics"]) and not str(results["physics"]).startswith("NO")
    json.dump(results, open(f"{OUT}/biasscan_results.json", "w"), indent=1, ensure_ascii=False, default=str)
    lines = [f"DAQ2-15 barrido de tensión de hod2018, campaña {CAMPAIGN}{' (SIMULACIÓN)' if SIM else ''}: " + ("ABORTADO" if results["aborted"] else "terminado"),
             "REARRANQUE DEL RUN DE FÍSICA: " + ("OK" if results.get("physics_ok") else "FALLÓ") + f" ({results.get('physics')})",
             f"Escalones completos: {[s['delta_V'] for s in results['steps']]}"]
    for k in results["skipped_steps"]:
        lines.append(f"Escalón OMITIDO {k['delta_V']:+.2f} V: {k['reason']}")
    if results["skipped"]:
        lines.append(f"Escalones no empezados por el límite de tiempo: {results['skipped']}")
    if results["aborted"]:
        lines.append(f"Escalón ABORTADO: {results.get('aborted_step')}: {results['aborted']}")
    lines += [f"Restauración: {results['restore'].get('result')} {results['restore'].get('detail', '')}",
              f"Comprobación previa del run de física: {results.get('physics_check', {}).get('detail', '–')}",
              f"Run de física antes del barrido: {'activo' if results.get('physics_before', {}).get('was_active') else 'NO activo'}"
              f" (última línea del log: {results.get('physics_before', {}).get('before_last', '–')})"]
    for s in results["steps"]:
        e = s["ext"]["all"]
        temps = [b["degC"] for b in s.get("board_temp", []) if b["degC"] is not None]
        mu = s.get("self_muons_per_h_sum")
        lines.append(f"ΔV {s['delta_V']:+.2f} V ({s['t_start_chile']}, placa {min(temps) if temps else '–'}-{max(temps) if temps else '–'} °C): "
                     f"G mediana {s['sw']['G_median_mV']:.3f} mV; ε_XY ctl 5/8/10/12 p.e. "
                     + "/".join(f"{e[T]['eff_ctl']:.3f}" for T in (5, 8, 10, 12) if T in e)
                     + f"; umbral {s['T_chosen_pe']} p.e. ({s['T_choice']}); muones/h auto-trigger (suma de pares) "
                     + (f"{mu:.1f}" if mu is not None else "– (runs omitidos)"))
        c = s.get("selftrigger_check", {})
        lines.append(f"   control {c.get('pair')} a {c.get('T_pe')} p.e.: " + (f"{c.get('rate_Hz', 0):.3f} Hz contra {c.get('pred_acc_Hz', 0):.3f} Hz previstos "
                     f"({c.get('frac_of_pred', 0):.3g}): " if "rate_Hz" in c else "") + ("FUNCIONA" if c.get("working") else
                     f"{'NO COMPROBABLE' if c.get('working') is None and 'rate_Hz' in c else 'NO FUNCIONA'} ({c.get('reason')})"))
        rg = c.get("registers") or {}
        if rg.get("read"):
            lines.append(f"   registros releídos: 0x810C {rg['read'].get('global_trigger_mask_0x810C')}, umbrales {rg['read'].get('group_trigger_threshold_adc')}, "
                         f"máscaras {rg['read'].get('channel_group_mask')}, mayoría {rg['read'].get('majority_level')}, ventana {rg['read'].get('majority_window_clk')}; "
                         f"distintos de lo escrito: {rg.get('mismatch') or 'ninguno'}")
        for a in s["self"]:
            lines.append(f"   {a['pair']}: {a['rate_Hz']:.2f} Hz, X·Y {a['N_XY_in_time']}, accidentales previstos {a['acc_pred']:.1f}, razón {a['ratio_pred']:.2f}")
    try:
        sh(["python3", "-I", CODE + "/daq/tools/biasscan/biasscan_report.py", OUT, JL], timeout=600)
        ana = json.load(open(f"{OUT}/biasscan_analysis.json"))
        vbd = [(ch, r["V_bd"], r["sigma_V_bd"], r["overvoltage_nominal"]) for ch, r in ana["channels"].items() if "V_bd" in r]
        if vbd:
            lines.append(f"V_bd (estimada) en {len(vbd)} canales: mediana {statistics.median(v[1] for v in vbd):.2f} V; sobretensión al nominal "
                         f"mediana {statistics.median(v[3] for v in vbd):.2f} V")
        lines.append(f"Informe: {OUT}/report_es.md")
    except Exception as e:  # noqa: BLE001
        lines.append(f"informe automático falló: {e}")
    lines.append(f"Registro: {LOG}; resultados: {OUT}/biasscan_results.json")
    results["summary_lines"] = lines
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
