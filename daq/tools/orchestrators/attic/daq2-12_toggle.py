#!/usr/bin/python3 -I
# DAQ2-12: segunda prueba de grupos HV (un grupo apagado a la vez) con daqrun en SW trigger 100 Hz y
# registro completo, para medir los canales que quedaron sin medir en DAQ2-10. Igual que
# daq2-10_toggle.py, con referencias largas (inicio, medio, final), 240 s de conteo por grupo y un
# corte por luz (IMon o unidad de luz de hod2019 > 1.2 × referencia inicial: se termina el grupo en
# curso, todo queda encendido y se reporta).
# Escrituras solo de Pw por hvctl (hvmon, allow_commands), cada una con --authorize. Cualquier timeout de
# verificación o Status inesperado aborta: no más escrituras salvo encender lo apagado.
import calendar, json, subprocess, time, urllib.request
LOG = '/data/daq/logs/daq2-12_toggle.jsonl'
HV = '/home/xplus/.local/opt/daqv2/bin/hvctl'
STATE = '/dev/shm/daqv2/hvmon_state.json'
MAP = '/home/xplus/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r3.json'
OK_STATUS = {0x0, 0x1, 0x3, 0x5}
REF_S, COUNT_S, SETTLE_MAX_S, AFTER_ON_WAIT_S, AFTER_ON_REF_S = 120, 240, 60, 15, 15
LIGHT_LIMIT = 1.2
DEADLINE_UTC = calendar.timegm(time.strptime('2026-10-21 17:30:00', '%Y-%m-%d %H:%M:%S'))
def log(**k):
    k['t_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); k['t'] = time.time()
    open(LOG, 'a').write(json.dumps(k, ensure_ascii=False) + '\n')
def note(t): subprocess.run([HV, 'note', t], capture_output=True)
def hv():
    s = json.load(open(STATE)); out = {}
    for c in s['channels']:
        v = c['values']; out[c['ch']] = {'pw': v.get('Pw'), 'st': v.get('Status'), 'vmon': v.get('VMon_V'), 'imon': v.get('IMon_A'), 'v0': v.get('V0Set_V')}
    return out
def dash():
    r = urllib.request.urlopen('http://127.0.0.1:8080/status/state_json/root.json', timeout=5).read()
    return json.loads(json.loads(r)['fTitle'])['dash']
class Abort(Exception): pass
def check_status():
    h = hv()
    for ch in range(8):
        if h[ch]['st'] not in OK_STATUS: raise Abort(f'Status inesperado CH{ch:02d} 0x{h[ch]["st"]:04X}')
    return h
def pw(ch, on, why):
    r = subprocess.run([HV, 'pw', '--ch', str(ch), 'on' if on else 'off', '--authorize', why], capture_output=True, text=True, timeout=400)
    out = (r.stdout + r.stderr).strip(); log(ev='pw', ch=ch, on=on, rc=r.returncode, out=out)
    if r.returncode != 0 or 'verified' not in out: raise Abort(f'CH{ch:02d} {"on" if on else "off"}: {out}')
last_period = [0]
def next_period(timeout=20):
    t0 = time.time()
    while time.time() - t0 < timeout:
        d = dash(); p = d.get('last_period_end_utc_ns', 0)
        if p and p != last_period[0]:
            last_period[0] = p; return d
        time.sleep(0.5)
    raise Abort('sin períodos nuevos de daqrun')
def collect(seconds, tag, g=None):
    t0 = time.time(); periods = []
    h0 = check_status()
    while time.time() - t0 < seconds:
        d = next_period(); periods.append(d['channel_counts']); check_status()
    h1 = check_status()
    log(ev='window', tag=tag, group=g, start=t0, end=time.time(), periods=len(periods),
        counts=[sum(p[c] for p in periods) for c in range(64)],
        vmon_start={c: h0[c]['vmon'] for c in range(8)}, vmon_end={c: h1[c]['vmon'] for c in range(8)},
        imon_start={c: h0[c]['imon'] for c in range(8)}, imon_end={c: h1[c]['imon'] for c in range(8)})
    return [sum(p[c] for p in periods) / max(1, len(periods)) for c in range(64)], h1
m = json.load(open(MAP))
expected = {g: [x['dig_ch'] for x in m['channels'] if x['hv_group'] == g + 1] for g in range(8)}
light = sorted(set(sum(m.get('light_units', {}).values(), [])))
light_stop = None
try:
    note('DAQ2-12: inicio de la segunda prueba de grupos HV (SW trigger 100 Hz), autorizada por René')
    time.sleep(10); next_period(60)
    ref, href = collect(REF_S, 'ref_inicio')
    light_ref = sum(ref[c] for c in light)
    for g in range(8):
        if time.time() > DEADLINE_UTC - 600: raise Abort('cerca del límite de las 17:30Z: no se empieza otro grupo')
        why = f'DAQ2-12 grupo {g+1} (CH{g:02d}), authorized by René'
        note(f'DAQ2-12: CH{g:02d} off'); pw(g, False, why + ' off')
        toff = time.time(); watch = set(c for c in expected[g] if ref[c] >= 10)
        while True:
            d = next_period(); cnt = d['channel_counts']; check_status()
            for c in range(64):
                if ref[c] >= 10 and cnt[c] < 0.3 * ref[c]: watch.add(c)
            low = all(cnt[c] <= 0.03 * ref[c] + 2 for c in watch)
            log(ev='wait', group=g, dt=time.time() - toff, watch=sorted(watch), low=low, counts={c: cnt[c] for c in sorted(watch)})
            if low or time.time() - toff > SETTLE_MAX_S: break
        log(ev='settled', group=g, dt=time.time() - toff, reached=low, watch=sorted(watch))
        collect(COUNT_S, 'off', g)
        note(f'DAQ2-12: CH{g:02d} on'); pw(g, True, why + ' on')
        t1 = time.time()
        while True:
            h = check_status()[g]
            if h['st'] == 0x1 and abs(h['vmon'] - h['v0']) < 0.5: break
            if time.time() - t1 > 200: raise Abort(f'CH{g:02d} no llegó a ON en 200 s')
            time.sleep(1)
        log(ev='on_stable', group=g, dt=time.time() - t1)
        time.sleep(AFTER_ON_WAIT_S)
        r2, h2 = collect(AFTER_ON_REF_S, 'ref_entre', g)
        # Corte por luz: unidad de luz de hod2019 (todo encendido) o IMon de los grupos que no se tocaron.
        lr = sum(r2[c] for c in light) / light_ref if light_ref > 0 else 0
        ir = {c: h2[c]['imon'] / href[c]['imon'] for c in range(8) if c != g and href[c]['imon']}
        log(ev='light_check', group=g, light_ratio=lr, imon_ratio=ir)
        if lr > LIGHT_LIMIT or any(v > LIGHT_LIMIT for v in ir.values()):
            light_stop = f'luz: unidad {lr:.2f}, IMon máx {max(ir.values()):.2f} tras el grupo {g+1}'
            note('DAQ2-12: corte por luz, ' + light_stop); log(ev='light_stop', reason=light_stop)
            break
        if g == 3: collect(REF_S, 'ref_medio')
    collect(REF_S, 'ref_final')
    note('DAQ2-12: fin de la prueba de grupos HV' + (f' (cortada: {light_stop})' if light_stop else '')); log(ev='done', light_stop=light_stop)
except Exception as e:
    log(ev='abort', reason=str(e)); note(f'DAQ2-12 ABORTADA: {e}')
finally:
    try:
        h = hv()
        for ch in range(8):
            if h[ch]['pw'] != 1:
                subprocess.run([HV, 'pw', '--ch', str(ch), 'on', '--authorize', f'DAQ2-12 final: dejar CH{ch:02d} ON, authorized by René'], timeout=400)
                log(ev='restore_on', ch=ch)
    except Exception as e2:
        log(ev='restore_failed', reason=str(e2))
    subprocess.run(['systemctl', '--user', 'stop', 'daq2-12-toggle-run'])
    log(ev='daqrun_stopped')
