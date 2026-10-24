# may_common.py — PROC-1: lectores compartidos de calibración y de tiempo bueno (antes dentro de hitdef_analyse.py;
# los usan hitdef_analyse.py y los scripts de daq/tools/proceeding/, sin copiar la lógica).
#  - g_may(), p0_may(): G y p0 por canal de hod2018 de outputs/pe_calibration_v8.txt (convención de mayo);
#  - g_c1(): G por canal de hod2018 de la calibración de octubre 2026_10-c1;
#  - good_time(t, run): máscara de tiempo bueno de runs/good_intervals_2026.csv para el run de mayo `run` (fuera los
#    intervalos EXCLUDE y HV_OFF; horas locales UTC−4; el minuto final incluido), con t en s UTC.
import calendar, csv, json, os, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
MAY_CAL = REPO + "/outputs/pe_calibration_v8.txt"
C1 = os.path.expanduser("~/.local/opt/daqv2/share/calibration/pe_calibration_2026_10_c1.json")
GOOD = REPO + "/runs/good_intervals_2026.csv"


def _may_column(col, fill):
    g = np.full(32, fill)
    with open(MAY_CAL) as fh:
        for l in fh:
            if l.startswith("#") or not l.split():
                continue
            f = l.split()
            if int(f[0]) < 32:
                g[int(f[0])] = float(f[col])
    return g


def g_may():
    return _may_column(3, np.nan)


def p0_may():
    return _may_column(2, 0.0)


def g_c1():
    with open(C1) as fh:
        c = json.load(fh)
    g = np.full(32, np.nan)
    for e in c["channels"]:
        if e["ch"] < 32 and isinstance(e.get("G_mV"), (int, float)):
            g[e["ch"]] = e["G_mV"]
    return g


def good_intervals(run):
    """Intervalos (status, inicio, fin, canales, motivo) del run, con inicio/fin en s 'UTC−4 como UTC'."""
    out = []
    with open(GOOD) as fh:
        for r in csv.DictReader(fh):
            if r["run"] != run or r["status"] not in ("EXCLUDE", "HV_OFF"):
                continue
            a = calendar.timegm(time.strptime(r["start_local"], "%Y-%m-%d %H:%M"))
            b = calendar.timegm(time.strptime(r["end_local"], "%Y-%m-%d %H:%M")) + 60  # el minuto final incluido
            out.append((r["status"], a, b, r["channels"], r["reason"]))
    return out


def good_time(t, run="v8"):
    keep = np.ones(len(t), bool)
    loc = np.asarray(t) - 4 * 3600
    for _, a, b, _, _ in good_intervals(run):
        keep &= ~((loc >= a) & (loc < b))
    return keep
