#!/usr/bin/python3 -I
# morning_report_waiter.py — XPLUS-1 parte B: espera a que termine el barrido de noche y escribe el informe de la mañana
# (updown_report.py), en una unidad de usuario de systemd (daq/systemd/morning_report.service) en vez de una espera atada a una
# sesión interactiva (DAQ2-17 lo hacía en segundo plano dentro de la sesión).
#  - cada --period-s (60 s) mira `systemctl --user is-active <unidad del barrido>`; cuando deja de estar active/activating y
#    existe <noche>/biasscan_results.json, corre updown_report.py <día> <noche> <salida>;
#  - manda un correo corto con la ruta del informe (hvctl notify-summary, kind morning_report) o, si se pasa la hora límite o el
#    informe falla, un correo de fallo. --notify-cmd reemplaza el correo (pruebas: un stub que guarda el texto).
# Configuración por variables de entorno (EnvironmentFile ~/.config/daqv2/morning_report.env) o argumentos:
#   SCAN_UNIT, DAY_DIR, NIGHT_DIR, OUT_DIR, DEADLINE_UTC ("AAAA-MM-DD HH:MM:SS"), NOTIFY_CMD (opcional)
# Uso: morning_report_waiter.py [--scan-unit U] [--day D] [--night N] [--out O] [--deadline T] [--period-s S] [--notify-cmd C]
import argparse, calendar, os, shlex, subprocess, sys, time

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
HVCTL = os.path.expanduser("~/.local/opt/daqv2/bin/hvctl")


def parse():
    e = os.environ.get
    p = argparse.ArgumentParser()
    p.add_argument("--scan-unit", default=e("SCAN_UNIT", "biasscan_night.service"))
    p.add_argument("--day", default=e("DAY_DIR"))
    p.add_argument("--night", default=e("NIGHT_DIR"))
    p.add_argument("--out", default=e("OUT_DIR"))
    p.add_argument("--deadline", default=e("DEADLINE_UTC"))
    p.add_argument("--period-s", type=float, default=float(e("PERIOD_S", "60")))
    p.add_argument("--notify-cmd", default=e("NOTIFY_CMD") or None)
    a = p.parse_args()
    for k in ("day", "night", "out", "deadline"):
        if not getattr(a, k):
            p.error(f"falta {k} (argumento o variable de entorno)")
    return a


def log(msg):
    print(time.strftime("%Y-%m-%dT%H:%M:%SZ ", time.gmtime()) + msg, flush=True)


def mail(a, kind, subject, body):
    os.makedirs(a.out, exist_ok=True)
    p = f"{a.out}/mail_{kind}.txt"
    open(p, "w").write(f"{subject}\n\n{body}\n")
    cmd = shlex.split(a.notify_cmd) + [p] if a.notify_cmd else [HVCTL, "notify-summary", "--text-file", p, "--kind", "morning_report"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    log(f"correo {kind}: salida {r.returncode} {(r.stdout or r.stderr).strip()[:200]}")


def main():
    a = parse()
    deadline = calendar.timegm(time.strptime(a.deadline, "%Y-%m-%d %H:%M:%S"))
    log(f"esperando el fin de {a.scan_unit} (límite {a.deadline} UTC); noche: {a.night}")
    while True:
        st = subprocess.run(["systemctl", "--user", "is-active", a.scan_unit], capture_output=True, text=True).stdout.strip()
        done = st not in ("active", "activating", "deactivating", "reloading")
        if done and os.path.exists(a.night + "/biasscan_results.json"):
            log(f"{a.scan_unit} terminó ({st}); escribiendo el informe")
            r = subprocess.run(["python3", "-I", "-B", REPO + "/daq/tools/biasscan/updown_report.py", a.day, a.night, a.out],
                               capture_output=True, text=True, timeout=3600, cwd=REPO)
            open(a.out + "/updown_report.log", "w").write(r.stdout + r.stderr)
            if r.returncode == 0:
                mail(a, "done", "[DAQ] informe de la mañana listo", f"Informe de la mañana: {a.out}\nBarrido de noche: {a.night}\nEstado de {a.scan_unit}: {st}")
                return 0
            mail(a, "failure", "[DAQ] el informe de la mañana falló", f"updown_report.py salió con {r.returncode}; registro: {a.out}/updown_report.log")
            return 1
        if done and time.time() > deadline:
            mail(a, "failure", "[DAQ] no hay resultados del barrido de noche", f"{a.scan_unit} está {st} y no existe {a.night}/biasscan_results.json "
                                                                                 f"a la hora límite {a.deadline} UTC.")
            return 2
        if time.time() > deadline + 4 * 3600:
            mail(a, "failure", "[DAQ] el barrido de noche no terminó", f"{a.scan_unit} sigue {st} 4 h después de la hora límite {a.deadline} UTC.")
            return 3
        time.sleep(a.period_s)


if __name__ == "__main__":
    sys.exit(main())
