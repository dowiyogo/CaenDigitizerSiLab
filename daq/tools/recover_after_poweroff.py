#!/usr/bin/python3 -I
# recover_after_poweroff.py — XPLUS-1 parte B (DAQ2-18a etapa 4): recuperación del montaje tras el corte del 2026-10-23 21:27Z.
# El texto de DAQ2-18 del 2026-10-24 no está en este host; los pasos son los de la tarea XPLUS-1 (René, 2026-10-25):
#  1. esperar a que respondan 192.168.1.118 (ping) y el enlace CONET (digprobe info, solo lectura), revisando cada --poll-s
#     (300 s); además, hvmon conectado a la fuente (hvctl status). Correo una vez cuando todo responde.
#  2. leer y reportar todos los parámetros del SY5527 de CH00-07 (última lectura de hvmon, hvctl read --json; se espera una
#     lectura completa posterior a la reconexión).
#  3. si se perdieron los ajustes: V0Set de un canal apagado distinto del nominal -> hvctl bias nominal (V0Set al nominal que
#     hvmon registró, solo con el canal apagado). Cualquier otro ajuste distinto de la referencia (I0Set, Trip, RUp, RDWn, SVMax,
#     V1Set, I1Set; daq/config/hv_reference_2026_10_23.json) -> se detiene y avisa: hvmon no tiene vía de escritura para ellos.
#  4. Pw = 1 canal por canal (hvctl pw), con 2 min de vigilancia de IMon tras cada uno: si IMon > 10 µA, ese canal se apaga,
#     correo y la recuperación se detiene.
#  5. 30 min de vigilancia de IMon de CH00-07 (> 10 µA: se apaga el canal, correo, se detiene) y comparación de la mediana por
#     canal con el run 8 (3.1-7.3 µA; fuera del intervalo: aviso, no se detiene).
#  6. arranque del run de física (como biasscan.py physics_start: unidad daq-physics, daqrun instalado, configuración de
#     física) y comprobación de triggers en la página del run: si no hay ninguno en 10 min, correo y se detiene.
#  7. correo de resumen.
# Cada paso que falla detiene la recuperación y manda un correo. Todas las escrituras van por hvmon (hvctl) con el texto de
# autorización; el script no habla con la Wrapper ni con la placa salvo digprobe info (lectura de registros).
# Se niega a correr sin --authorize "DAQ2-18 recovery, authorized by René".
# Pruebas: --env test usa un hvmon de prueba con FakeHVBackend (--socket), un daqrun con el digitizador falso
# (--phys-config de prueba, --http-port), --conet-cmd a un stub (digprobe abre SIEMPRE el VX1740 real, sea cual sea el backend de
# la configuración: está prohibido en --env test), --host 127.0.0.1 y --notify-cmd a un stub que guarda los correos en un archivo.
# No correr contra el hardware real salvo el lunes, autorizado por René.
# Salida: <out>/recovery_log.jsonl (un registro por evento) y <out>/recovery_summary.json.
# Uso: recover_after_poweroff.py --authorize "DAQ2-18 recovery, authorized by René" [--env real|test] [opciones; --help]
import argparse, json, os, shlex, statistics, subprocess, sys, time

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
BIN = os.path.expanduser("~/.local/opt/daqv2/bin")
AUTH_TEXT = "DAQ2-18 recovery, authorized by René"
CHANNELS = list(range(8))


class Stop(Exception):
    """Fallo de un paso: se detiene la recuperación (con correo)."""


def parse():
    p = argparse.ArgumentParser(description="Recuperación tras el corte (XPLUS-1 parte B, DAQ2-18a etapa 4)")
    p.add_argument("--authorize", required=False, default="")
    p.add_argument("--env", choices=("real", "test"), default="real")
    p.add_argument("--host", default="192.168.1.118")
    p.add_argument("--hvctl", default=BIN + "/hvctl")
    p.add_argument("--socket", default=None, help="socket de hvmon (por omisión el de hvctl)")
    p.add_argument("--digprobe", default=BIN + "/digprobe")
    p.add_argument("--conet-cmd", default=None, help="comprobación del enlace CONET (por omisión: digprobe info --config <física>; "
                                                     "en --env test es obligatoria y no puede ser digprobe, que abre siempre el VX1740 real)")
    p.add_argument("--daqrun", default=BIN + "/daqrun")
    p.add_argument("--phys-config", default="/data/daq/config/hodoscope_2026_10/physics_ext_B.json")
    p.add_argument("--http-port", type=int, default=8080)
    p.add_argument("--reference", default=REPO + "/daq/config/hv_reference_2026_10_23.json")
    p.add_argument("--notify-cmd", default=None, help="comando que recibe el texto en un archivo (por omisión hvctl notify-summary)")
    p.add_argument("--out", default=None)
    p.add_argument("--poll-s", type=float, default=300.0)
    p.add_argument("--max-wait-h", type=float, default=96.0)
    p.add_argument("--pw-watch-s", type=float, default=120.0)
    p.add_argument("--watch-s", type=float, default=1800.0)
    p.add_argument("--watch-period-s", type=float, default=10.0)
    p.add_argument("--trig-wait-s", type=float, default=600.0)
    p.add_argument("--imon-max-uA", type=float, default=10.0)
    p.add_argument("--run8-uA", default="3.1,7.3")
    p.add_argument("--fresh-read-s", type=float, default=75.0, help="espera tras la reconexión para que hvmon haga una lectura completa (cada 60 s)")
    return p.parse_args()


class Recovery:
    def __init__(self, a):
        self.a = a
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        self.out = a.out or (f"/data/daq/logs/recovery_{stamp}" if a.env == "real" else f"/data/daq/sim/logs/recovery_test_{stamp}")
        os.makedirs(self.out, exist_ok=True)
        self.logf = open(self.out + "/recovery_log.jsonl", "a")
        self.summary = {"env": a.env, "started_utc": self.now(), "steps": {}, "emails": [], "result": None}
        self.hv = [a.hvctl] + (["--socket", a.socket] if a.socket else [])
        self.phys = None  # proceso de daqrun en --env test

    # ---------------------------------------------------------------- utilidades
    @staticmethod
    def now():
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    def log(self, step, msg, **k):
        rec = {"t_utc": self.now(), "step": step, "msg": msg, **k}
        self.logf.write(json.dumps(rec, ensure_ascii=False, default=str) + "\n")
        self.logf.flush()
        print(f"[{rec['t_utc']}] {step}: {msg}", flush=True)

    def run(self, cmd, timeout=600):
        env = dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser("~/root6/lib") + ":" + os.environ.get("LD_LIBRARY_PATH", ""))
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        except subprocess.TimeoutExpired as e:
            return subprocess.CompletedProcess(cmd, 124, e.stdout or "", (e.stderr or "") + f"\ntimeout {timeout} s")

    def hvctl(self, *args, timeout=600):
        return self.run(self.hv + list(args), timeout=timeout)

    def mail(self, kind, subject, body):
        p = f"{self.out}/mail_{len(self.summary['emails']) + 1:02d}_{kind}.txt"
        open(p, "w").write(f"{subject}\n\n{body}\n")
        cmd = shlex.split(self.a.notify_cmd) + [p] if self.a.notify_cmd else self.hv + ["notify-summary", "--text-file", p, "--kind", "recovery_" + kind]
        r = self.run(cmd, timeout=300)
        self.summary["emails"].append({"t_utc": self.now(), "kind": kind, "subject": subject, "file": p, "exit": r.returncode,
                                       "reply": (r.stdout or r.stderr).strip()[:300]})
        self.log("mail", f"{kind}: {subject} (salida {r.returncode})")

    def state(self):
        r = self.hvctl("status", "--json", timeout=60)
        if r.returncode != 0:
            return None
        try:
            return json.loads(r.stdout)
        except ValueError:
            return None

    def reading(self):
        r = self.hvctl("read", "--json", timeout=60)
        if r.returncode != 0:
            raise Stop("hvctl read falló: " + (r.stderr or r.stdout).strip()[:300])
        d = json.loads(r.stdout)
        out = {"last_good_read_utc_ns": d.get("last_good_read_utc_ns"), "channels": {}}
        for prm in d["params"]:
            for ch, v in zip(d["channels"], prm.get("values_si", prm.get("values", []))):  # numéricos en SI; Status/Pw/on-off en "values"
                out["channels"].setdefault(ch, {})[prm["name"]] = v
        return out

    def imon_uA(self):
        rd = self.reading()
        return {ch: rd["channels"][ch]["IMon"] * 1e6 for ch in CHANNELS}, rd

    def pw(self, ch, on):
        r = self.hvctl("pw", "--ch", str(ch), "on" if on else "off", "--authorize", self.a.authorize, timeout=900)
        self.log("pw", f"CH{ch:02d} {'on' if on else 'off'}: salida {r.returncode}", out=(r.stdout or r.stderr).strip()[:400])
        return r.returncode == 0

    # ---------------------------------------------------------------- pasos
    def step1_wait_links(self):
        t0 = time.time()
        tries = 0
        while True:
            tries += 1
            ping = self.run(["ping", "-c", "1", "-W", "3", self.a.host], timeout=20).returncode == 0
            st = self.state()
            hv = bool(st and st.get("connection", {}).get("connected"))
            dg = self.run(shlex.split(self.a.conet_cmd) if self.a.conet_cmd else [self.a.digprobe, "info", "--config", self.a.phys_config], timeout=120)
            conet = dg.returncode == 0
            self.log("1_links", f"intento {tries}: ping {self.a.host} {'sí' if ping else 'no'}, hvmon conectado {'sí' if hv else 'no'}, "
                                f"CONET (digprobe info) {'sí' if conet else 'no'}", digprobe=(dg.stdout if conet else dg.stderr).strip()[-300:])
            if ping and hv and conet:
                self.summary["steps"]["1_links"] = {"ok": True, "tries": tries, "waited_s": round(time.time() - t0, 1)}
                self.mail("links_up", "[DAQ2-18] la fuente HV y el enlace CONET responden",
                          f"A las {self.now()} responden {self.a.host} (ping), hvmon está conectado a la fuente y el digitizador abre por CONET "
                          f"(digprobe info). La recuperación sigue con la lectura de parámetros.")
                return
            if time.time() - t0 > self.a.max_wait_h * 3600:
                raise Stop(f"sin enlace tras {self.a.max_wait_h} h: ping {ping}, hvmon {hv}, CONET {conet}")
            time.sleep(self.a.poll_s)

    def step2_read_params(self):
        # hvmon hace una lectura completa (todos los parámetros) cada full_period_s = 60 s: se espera más que eso desde la
        # reconexión para no reportar ajustes de antes del corte.
        time.sleep(self.a.fresh_read_s)
        st = self.state() or {}
        if not st.get("connection", {}).get("connected"):
            raise Stop("hvmon perdió la conexión antes de la lectura de parámetros")
        rd = self.reading()
        if rd["channels"].get(0, {}).get("SVMax") is None:
            raise Stop("la lectura de hvmon no tiene los parámetros de la lectura completa")
        params = {ch: rd["channels"][ch] for ch in CHANNELS}
        json.dump({"t_utc": self.now(), "board_temperature_degC": st.get("board_temperature_degC"), "channels": params},
                  open(self.out + "/sy5527_params_ch00_07.json", "w"), indent=1)
        lines = []
        for ch in CHANNELS:
            p = params[ch]
            lines.append(f"CH{ch:02d}: " + ", ".join(f"{k} {v:.6g}" for k, v in sorted(p.items()) if isinstance(v, (int, float))))
        self.summary["steps"]["2_params"] = {"ok": True, "file": self.out + "/sy5527_params_ch00_07.json"}
        self.log("2_params", "parámetros de CH00-07 leídos", params=params)
        self.mail("params", "[DAQ2-18] parámetros del SY5527 (CH00-07) tras el corte", "\n".join(lines))
        return params

    def step3_restore_settings(self, params):
        ref = json.load(open(self.a.reference))
        tol = ref["tolerance"]
        names = {"V0Set_V": "V0Set", "I0Set_A": "I0Set", "V1Set_V": "V1Set", "I1Set_A": "I1Set", "RUp_V_per_s": "RUp", "RDWn_V_per_s": "RDWn",
                 "Trip_s": "Trip", "SVMax_V": "SVMax"}
        lost_v0, other = [], []
        for ch in CHANNELS:
            r = ref["channels"][str(ch)]
            for k, n in names.items():
                v = params[ch].get(n)
                if v is None:
                    other.append(f"CH{ch:02d} {n}: sin lectura")
                elif abs(v - r[k]) > tol[k]:
                    (lost_v0 if n == "V0Set" else other).append((ch, n, v, r[k]) if n == "V0Set" else f"CH{ch:02d} {n} = {v:.6g} (referencia {r[k]:.6g})")
        self.log("3_settings", f"V0Set distinto del nominal: {len(lost_v0)}; otros ajustes distintos: {len(other)}", lost_v0=lost_v0, other=other)
        if other:
            raise Stop("ajustes distintos de la referencia sin vía de escritura en hvmon: " + "; ".join(other))
        done = []
        for ch, n, v, nom in lost_v0:
            if params[ch].get("Pw", 0) != 0:
                raise Stop(f"CH{ch:02d}: V0Set {v:.3f} V distinto del nominal {nom:.3f} V con el canal encendido")
            r = self.hvctl("bias", "nominal", "--ch", str(ch), "--authorize", self.a.authorize, timeout=600)
            self.log("3_settings", f"CH{ch:02d}: V0Set {v:.3f} -> nominal: salida {r.returncode}", out=(r.stdout or r.stderr).strip()[:400])
            if r.returncode != 0:
                raise Stop(f"CH{ch:02d}: no se pudo volver al V0Set nominal: {(r.stdout or r.stderr).strip()[:300]}")
            done.append(ch)
        self.summary["steps"]["3_settings"] = {"ok": True, "v0set_restored": done}

    def step4_pw_on(self):
        on = []
        for ch in CHANNELS:
            rd = self.reading()
            if rd["channels"][ch].get("Pw") == 1:
                self.log("4_pw", f"CH{ch:02d} ya encendido")
            elif not self.pw(ch, True):
                raise Stop(f"CH{ch:02d}: Pw = 1 no se verificó (ver recovery_log.jsonl)")
            t_end = time.time() + self.a.pw_watch_s
            mx = 0.0
            while time.time() < t_end:
                im, _ = self.imon_uA()
                mx = max(mx, im[ch])
                if im[ch] > self.a.imon_max_uA:
                    self.pw(ch, False)
                    raise Stop(f"CH{ch:02d}: IMon {im[ch]:.2f} µA > {self.a.imon_max_uA} µA en los {self.a.pw_watch_s:.0f} s tras encender; canal apagado")
                time.sleep(min(5.0, self.a.watch_period_s))
            on.append({"ch": ch, "imon_max_uA": round(mx, 3)})
            self.log("4_pw", f"CH{ch:02d} encendido; IMon máx {mx:.2f} µA en {self.a.pw_watch_s:.0f} s")
        self.summary["steps"]["4_pw"] = {"ok": True, "channels": on}

    def step5_watch(self):
        lo, hi = (float(x) for x in self.a.run8_uA.split(","))
        series = {ch: [] for ch in CHANNELS}
        t_end = time.time() + self.a.watch_s
        while time.time() < t_end:
            im, rd = self.imon_uA()
            for ch in CHANNELS:
                series[ch].append(im[ch])
                if im[ch] > self.a.imon_max_uA:
                    self.pw(ch, False)
                    raise Stop(f"CH{ch:02d}: IMon {im[ch]:.2f} µA > {self.a.imon_max_uA} µA durante la vigilancia de {self.a.watch_s:.0f} s; canal apagado")
                if rd["channels"][ch].get("Pw") != 1:
                    raise Stop(f"CH{ch:02d}: Pw = 0 durante la vigilancia (¿trip?)")
            time.sleep(self.a.watch_period_s)
        res = {ch: {"median_uA": round(statistics.median(v), 3), "max_uA": round(max(v), 3), "n": len(v),
                    "in_run8_range": lo <= statistics.median(v) <= hi} for ch, v in series.items()}
        outside = [ch for ch, r in res.items() if not r["in_run8_range"]]
        self.summary["steps"]["5_watch"] = {"ok": True, "per_channel": res, "outside_run8_range": outside, "run8_range_uA": [lo, hi]}
        self.log("5_watch", f"vigilancia terminada; fuera del intervalo del run 8 ({lo}-{hi} µA): {outside}", per_channel=res)
        if outside:
            self.mail("imon_warning", "[DAQ2-18] IMon fuera del intervalo del run 8 (aviso)",
                      "\n".join(f"CH{ch:02d}: mediana {res[ch]['median_uA']} µA, máx {res[ch]['max_uA']} µA" for ch in outside)
                      + f"\nIntervalo del run 8: {lo}-{hi} µA. La recuperación sigue (no es un fallo).")

    def physics_start(self):
        if self.a.env == "test":
            log = open(self.out + "/daqrun_test.log", "w")
            self.phys = subprocess.Popen([self.a.daqrun, "--config", self.a.phys_config], stdout=log, stderr=subprocess.STDOUT, cwd=REPO,
                                         env=dict(os.environ, LD_LIBRARY_PATH=os.path.expanduser("~/root6/lib")))
            return "daqrun de prueba lanzado (pid %d)" % self.phys.pid
        # igual que biasscan.py physics_start (DAQ2-15): unidad daq-physics con el daqrun instalado y la configuración de física
        act = self.run(["systemctl", "--user", "is-active", "daq-physics"]).stdout.strip()
        if act == "active":
            return "ya estaba activo (no se reinició)"
        cmd = ["systemd-run", "--user", "--unit=daq-physics", "--collect", "--property=RuntimeMaxSec=1209600", "--property=KillSignal=SIGTERM",
               "--property=TimeoutStopSec=300", "--property=StandardOutput=append:/data/daq/logs/daq-physics.log",
               "--property=StandardError=append:/data/daq/logs/daq-physics.log", "--setenv=LD_LIBRARY_PATH=" + os.path.expanduser("~/root6/lib"),
               "--working-directory=" + os.path.expanduser("~"), self.a.daqrun, "--config", self.a.phys_config]
        r = self.run(cmd, timeout=60)
        return f"systemd-run salió con {r.returncode}: {(r.stdout or r.stderr).strip()[:200]}"

    def page(self):
        import urllib.request
        r = urllib.request.urlopen(f"http://127.0.0.1:{self.a.http_port}/status/state_json/root.json", timeout=10).read()
        return json.loads(json.loads(r)["fTitle"])

    def step6_physics(self):
        msg = self.physics_start()
        self.log("6_physics", msg)
        t_end = time.time() + self.a.trig_wait_s
        last = "sin respuesta de la página del run"
        while time.time() < t_end:
            time.sleep(min(15.0, self.a.watch_period_s))
            try:
                s = self.page()
                ev = s["dash"].get("events") or 0
                last = f"run {s['dash'].get('run')}, {ev} eventos, errores de lectura {s['readout'].get('read_errors')}"
                if ev > 0:
                    self.summary["steps"]["6_physics"] = {"ok": True, "start": msg, "detail": last}
                    self.log("6_physics", "hay triggers: " + last)
                    return
            except Exception as e:  # noqa: BLE001
                last = f"la página del run no responde: {e}"
        raise Stop(f"ningún trigger en {self.a.trig_wait_s:.0f} s tras arrancar el run de física ({last})")

    # ---------------------------------------------------------------- orquestación
    def main(self):
        try:
            self.step1_wait_links()
            params = self.step2_read_params()
            self.step3_restore_settings(params)
            self.step4_pw_on()
            self.step5_watch()
            self.step6_physics()
            self.summary["result"] = "ok"
            w = self.summary["steps"].get("5_watch", {})
            self.mail("done", "[DAQ2-18] recuperación terminada",
                      f"CH00-07 encendidos; IMon en 30 min: " + ", ".join(f"CH{ch:02d} {r['median_uA']} µA" for ch, r in w.get("per_channel", {}).items())
                      + f".\nRun de física: {self.summary['steps']['6_physics']['detail']}.")
            return 0
        except Stop as e:
            self.summary["result"] = "stopped"
            self.summary["stop_reason"] = str(e)
            self.log("stop", str(e))
            self.mail("failure", "[DAQ2-18] recuperación DETENIDA", f"{e}\n\nRegistro: {self.out}/recovery_log.jsonl")
            return 2
        except Exception as e:  # noqa: BLE001
            self.summary["result"] = "error"
            self.summary["stop_reason"] = f"error interno: {type(e).__name__}: {e}"
            self.log("stop", self.summary["stop_reason"])
            self.mail("failure", "[DAQ2-18] recuperación DETENIDA (error interno)", f"{self.summary['stop_reason']}\n\nRegistro: {self.out}/recovery_log.jsonl")
            return 3
        finally:
            self.summary["ended_utc"] = self.now()
            json.dump(self.summary, open(self.out + "/recovery_summary.json", "w"), indent=1, ensure_ascii=False, default=str)
            if self.phys is not None and self.phys.poll() is None and self.a.env == "test":
                self.phys.terminate()
                try:
                    self.phys.wait(120)
                except subprocess.TimeoutExpired:
                    self.phys.kill()


if __name__ == "__main__":
    args = parse()
    if args.authorize != AUTH_TEXT:
        print(f'recover_after_poweroff: se necesita --authorize "{AUTH_TEXT}"', file=sys.stderr)
        sys.exit(64)
    if args.env == "test" and (not args.socket or not args.notify_cmd or not args.conet_cmd or "digprobe" in args.conet_cmd):
        # digprobe abre siempre el VX1740 real (CaenDigitizerBackend), sin importar el backend de la configuración
        print("recover_after_poweroff: --env test necesita --socket (hvmon de prueba), --notify-cmd (stub de correo) y --conet-cmd "
              "(stub; no digprobe)", file=sys.stderr)
        sys.exit(64)
    sys.exit(Recovery(args).main())
