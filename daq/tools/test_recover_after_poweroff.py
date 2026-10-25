#!/usr/bin/python3 -I
# test_recover_after_poweroff.py — XPLUS-1 parte B: prueba de punta a punta de recover_after_poweroff.py contra un hvmon de
# prueba con FakeHVBackend (socket /run/user/1000/daqv2-rectest, estado en /data/daq/sim/state-recovery) y un daqrun con el
# digitizador falso (campaña simulada, página en el puerto 8095). Ningún hardware real, ningún correo real (stub que guarda el
# texto). Casos: éxito (V0Set perdido en CH03 vuelve al nominal), sin enlace, IMon demasiado alta, sin triggers, y el rechazo
# sin --authorize.
# Uso: test_recover_after_poweroff.py [caso ...]   (por omisión todos). Sale con 0 si todos pasan.
import json, os, signal, subprocess, sys, time

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
B = REPO + "/daq/build/apps"
CFGD = "/data/daq/sim/config/recovery_test"
STATE = "/data/daq/sim/state-recovery"
SOCK = "/run/user/1000/daqv2-rectest/hvmon.sock"
HVSTATE = "/dev/shm/daqv2-sim-rectest/hvmon_state.json"
ENV = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser("~/root6/lib"))
AUTH = "DAQ2-18 recovery, authorized by René"
NOMINAL = [65.75, 65.75, 65.90, 65.85, 65.75, 65.70, 66.15, 65.90]
CAMPAIGN = "recovery_test_" + time.strftime("%Y%m%d", time.gmtime())
STUB = CFGD + "/mail_stub.sh"
CONET_OK, CONET_NO = CFGD + "/conet_ok.sh", CFGD + "/conet_no.sh"  # stubs: digprobe abre siempre el VX1740 real


def sh(cmd, **k):
    return subprocess.run(cmd, capture_output=True, text=True, env=ENV, **k)


def setup_common():
    os.makedirs(CFGD, exist_ok=True)
    os.makedirs(STATE, exist_ok=True)
    open(STUB, "w").write('#!/bin/sh\n# stub de correo para la prueba: guarda el texto, no envía nada\nmkdir -p "$(dirname "$1")/sent"\ncp "$1" "$(dirname "$1")/sent/"\necho stub-ok\n')
    os.chmod(STUB, 0o755)
    for p, rc in ((CONET_OK, 0), (CONET_NO, 1)):
        open(p, "w").write(f"#!/bin/sh\n# stub del enlace CONET para la prueba (no abre el digitizador)\necho conet-stub {rc}\nexit {rc}\n")
        os.chmod(p, 0o755)
    # nominal registrado por hvmon (estado de sesión persistido y cerrado), como el de 2026-10-23
    json.dump({"schema": "daqv2.hvmon_bias_session/1", "active": False, "reason": "prueba XPLUS-1", "start_utc_ns": 0,
               "nominal_V0Set_V": {str(i): v for i, v in enumerate(NOMINAL)}, "vmon_offset_V": {str(i): -0.2 for i in range(8)},
               "limits": {"max_above_V": 0.0, "max_step_V": 0.5}}, open(STATE + "/hvmon_bias_state.json", "w"))
    os.chmod(STATE + "/hvmon_bias_state.json", 0o600)
    if sh([B + "/daqcampaign/daqcampaign", "--sim", "status", CAMPAIGN]).stdout.find("open") < 0:
        r = sh([B + "/daqcampaign/daqcampaign", "--sim", "open", CAMPAIGN])
        print("campaña:", r.stdout.strip()[-200:], r.stderr.strip()[-200:])


def hv_scenario(name, extra_events):
    # tras el corte: CH00-07 apagados; V0Set de CH03 perdido (0 V)
    ev = [{"t_s": 0.2, "type": "set", "ch": ch, "param": "Pw", "value": 0} for ch in range(8)]
    ev.append({"t_s": 0.2, "type": "set", "ch": 3, "param": "V0Set", "value": 0.0})
    ev += extra_events
    p = f"{CFGD}/fake_hv_{name}.json"
    json.dump({"base": "default", "seed": 25, "events": ev}, open(p, "w"), indent=1)
    return p


def dig_config(name, rate_hz):
    sc = json.load(open(REPO + "/daq/config/fake/dig_synthetic_default.json"))
    sc["trigger_rate_Hz"] = rate_hz
    scp = f"{CFGD}/dig_{name}.json"
    json.dump(sc, open(scp, "w"), indent=1)
    c = json.load(open(REPO + "/daq/config/daqrun.example.json"))
    c["campaign"] = {"name": CAMPAIGN, "operator": "test XPLUS-1"}
    c["digitizer"]["scenario_file"] = scp
    c["run"]["max_duration_s"] = 600
    c["hv"]["state_file"] = HVSTATE
    c["monitoring"]["http"]["port"] = 8095
    c["mapping"]["file"] = REPO + "/daq/config/mapping/hodo_mapping_2026_10_r4.json"
    p = f"{CFGD}/daqrun_{name}.json"
    json.dump(c, open(p, "w"), indent=1)
    return p


def start_hvmon(scenario):
    log = open(f"/data/daq/sim/logs/hvmon_recovery_test_{os.path.basename(scenario)}.out", "w")
    p = subprocess.Popen([B + "/hvmon/hvmon", "--fake", scenario, "--config", CFGD + "/hvmon_test.json"], stdout=log, stderr=subprocess.STDOUT, env=ENV)
    for _ in range(60):
        if os.path.exists(SOCK):
            break
        time.sleep(0.5)
    time.sleep(3)
    return p


def stop(p):
    if p and p.poll() is None:
        p.send_signal(signal.SIGTERM)
        try:
            p.wait(60)
        except subprocess.TimeoutExpired:
            p.kill()


def recovery(name, dig_cfg, extra=(), conet=CONET_OK):
    out = f"/data/daq/sim/logs/recovery_test_{name}_{time.strftime('%H%M%S', time.gmtime())}"
    cmd = ["python3", "-I", "-B", REPO + "/daq/tools/recover_after_poweroff.py", "--authorize", AUTH, "--env", "test",
           "--host", "127.0.0.1", "--hvctl", B + "/hvctl/hvctl", "--socket", SOCK, "--conet-cmd", conet,
           "--daqrun", B + "/daqrun/daqrun", "--phys-config", dig_cfg, "--http-port", "8095", "--notify-cmd", STUB, "--out", out,
           "--poll-s", "2", "--max-wait-h", str(20 / 3600), "--pw-watch-s", "6", "--watch-s", "20", "--watch-period-s", "2",
           "--trig-wait-s", "40", "--fresh-read-s", "65", *extra]
    t0 = time.time()
    r = sh(cmd, cwd=REPO, timeout=3600)
    summ = json.load(open(out + "/recovery_summary.json"))
    sent = sorted(os.listdir(out + "/sent")) if os.path.isdir(out + "/sent") else []
    return r.returncode, summ, sent, round(time.time() - t0, 1), out


def case_ok():
    hv = start_hvmon(hv_scenario("ok", []))
    try:
        rc, s, sent, dt, out = recovery("ok", dig_config("ok", 0.6))
    finally:
        stop(hv)
    v0 = s["steps"].get("3_settings", {}).get("v0set_restored")
    ok = rc == 0 and s["result"] == "ok" and v0 == [3] and len(s["steps"]["4_pw"]["channels"]) == 8 and s["steps"]["6_physics"]["ok"]
    return ok, {"rc": rc, "result": s["result"], "v0set_restored": v0, "emails": [e["kind"] for e in s["emails"]], "sent": sent, "s": dt, "out": out,
                "physics": s["steps"].get("6_physics"), "outside_run8": s["steps"].get("5_watch", {}).get("outside_run8_range")}


def case_no_link():
    hv = start_hvmon(hv_scenario("nolink", [{"t_s": 0.0, "type": "connection_loss", "duration_s": 100000}]))
    try:
        rc, s, sent, dt, out = recovery("nolink", dig_config("ok", 0.6))
    finally:
        stop(hv)
    ok = rc == 2 and s["result"] == "stopped" and "sin enlace" in s.get("stop_reason", "") and [e["kind"] for e in s["emails"]] == ["failure"]
    return ok, {"rc": rc, "result": s["result"], "stop_reason": s.get("stop_reason"), "emails": [e["kind"] for e in s["emails"]], "s": dt, "out": out}


def case_no_conet():
    hv = start_hvmon(hv_scenario("noconet", []))
    try:
        rc, s, sent, dt, out = recovery("noconet", dig_config("ok", 0.6), conet=CONET_NO)
    finally:
        stop(hv)
    ok = rc == 2 and "CONET False" in s.get("stop_reason", "") and [e["kind"] for e in s["emails"]] == ["failure"]
    return ok, {"rc": rc, "result": s["result"], "stop_reason": s.get("stop_reason"), "emails": [e["kind"] for e in s["emails"]], "s": dt, "out": out}


def case_trip_on():
    # +25 µA en CH02 (sobre I0Set = 15 µA): la placa dispara al encender y la verificación de Pw = 1 falla
    hv = start_hvmon(hv_scenario("tripon", [{"t_s": 0.5, "type": "current", "ch": 2, "extra_uA": 25.0, "duration_s": 100000}]))
    try:
        rc, s, sent, dt, out = recovery("tripon", dig_config("ok", 0.6))
    finally:
        stop(hv)
    ok = rc == 2 and "CH02" in s.get("stop_reason", "") and "no se verificó" in s.get("stop_reason", "") and "failure" in [e["kind"] for e in s["emails"]]
    return ok, {"rc": rc, "result": s["result"], "stop_reason": s.get("stop_reason"), "emails": [e["kind"] for e in s["emails"]], "s": dt, "out": out}


def case_imon_high():
    # +10 µA en CH02 (carga ~2 µA: IMon ~12 µA, bajo I0Set): la placa no dispara; la vigilancia de 2 min lo detecta y apaga CH02
    hv = start_hvmon(hv_scenario("imon", [{"t_s": 0.5, "type": "current", "ch": 2, "extra_uA": 10.0, "duration_s": 100000}]))
    try:
        rc, s, sent, dt, out = recovery("imon", dig_config("ok", 0.6))
        pw2 = sh([B + "/hvctl/hvctl", "--socket", SOCK, "read", "--json"]).stdout
        pw2 = {p["name"]: p.get("values", p.get("values_si"))[2] for p in json.loads(pw2)["params"] if p["name"] in ("Pw", "Status", "IMon")}
    finally:
        stop(hv)
    on = [c["ch"] for c in s["steps"].get("4_pw", {}).get("channels", [])] if "4_pw" in s["steps"] else None
    ok = rc == 2 and "CH02" in s.get("stop_reason", "") and "IMon" in s.get("stop_reason", "") and "failure" in [e["kind"] for e in s["emails"]] \
        and pw2.get("Pw") == 0
    return ok, {"rc": rc, "result": s["result"], "stop_reason": s.get("stop_reason"), "emails": [e["kind"] for e in s["emails"]], "s": dt, "out": out,
                "ch2_after": pw2}


def case_no_triggers():
    hv = start_hvmon(hv_scenario("notrig", []))
    try:
        rc, s, sent, dt, out = recovery("notrig", dig_config("notrig", 0.0))
    finally:
        stop(hv)
    ok = rc == 2 and "ningún trigger" in s.get("stop_reason", "") and "failure" in [e["kind"] for e in s["emails"]]
    return ok, {"rc": rc, "result": s["result"], "stop_reason": s.get("stop_reason"), "emails": [e["kind"] for e in s["emails"]], "s": dt, "out": out}


def case_no_authorize():
    r = sh(["python3", "-I", "-B", REPO + "/daq/tools/recover_after_poweroff.py", "--env", "test"])
    r2 = sh(["python3", "-I", "-B", REPO + "/daq/tools/recover_after_poweroff.py", "--authorize", "otra cosa"])
    r3 = sh(["python3", "-I", "-B", REPO + "/daq/tools/recover_after_poweroff.py", "--authorize", AUTH, "--env", "test", "--socket", SOCK,
             "--notify-cmd", STUB, "--conet-cmd", "/nonexistent/digprobe info"])  # rechazado antes de ejecutar nada (ruta inexistente)
    return [r.returncode, r2.returncode, r3.returncode] == [64, 64, 64], {"rc": [r.returncode, r2.returncode, r3.returncode], "stderr": r.stderr.strip()}


CASES = {"no_authorize": case_no_authorize, "ok": case_ok, "no_link": case_no_link, "no_conet": case_no_conet, "imon_high": case_imon_high,
         "trip_on": case_trip_on, "no_triggers": case_no_triggers}
if __name__ == "__main__":
    setup_common()
    names = sys.argv[1:] or list(CASES)
    res = {}
    for n in names:
        ok, info = CASES[n]()
        res[n] = {"pass": ok, **info}
        print(("PASA " if ok else "FALLA ") + n + ": " + json.dumps(info, ensure_ascii=False, default=str)[:900], flush=True)
    json.dump(res, open(CFGD + "/test_results.json", "w"), indent=1, ensure_ascii=False, default=str)
    sys.exit(0 if all(r["pass"] for r in res.values()) else 1)
