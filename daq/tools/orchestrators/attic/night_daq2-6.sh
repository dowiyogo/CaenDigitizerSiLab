#!/bin/bash
# DAQ2-6, noche 2026-10-19/20: repetición de 4.1 (trigger externo, ventana B, 30 min) con la sala más
# fresca y 4.2 (SW trigger 100 Hz, registro completo, 2 min). Lo arranca un timer de systemd a las
# 01:00Z. Si la tasa de trigger a los 3 min es < 0.1 Hz, para el run (cierre limpio) y reintenta cada
# 30 min; último intento a las 06:00Z (03:00 en Chile). Después (o si ningún intento sirvió) corre 4.2.
# Log: /data/daq/logs/night_daq2-6.log; marca de fin: /data/daq/logs/night_daq2-6.done
set -u
L=/data/daq/logs/night_daq2-6.log
C=/data/daq/config/commissioning_2026_10
BIN=$HOME/.local/opt/daqv2/bin
log() { echo "$(date -u +%FT%TZ) $*" >> "$L"; }
hv() { log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"; }
startrun() {  # $1 unidad, $2 config, $3 RuntimeMaxSec
  systemd-run --user --unit="$1" --collect --property=RuntimeMaxSec="$3" \
    --property=StandardOutput=append:/data/daq/logs/$1.log --property=StandardError=append:/data/daq/logs/$1.log \
    --setenv=LD_LIBRARY_PATH=$HOME/root6/lib --working-directory=$HOME $BIN/daqrun --config "$2" >> "$L" 2>&1
  systemd-run --user --unit="$1-sampler" --collect --property=RuntimeMaxSec=$(( $3 + 600 )) \
    /data/daq/logs/sampler.sh "$1" 60 /data/daq/logs/$1_ps.log >> "$L" 2>&1
}
events() { curl -s http://127.0.0.1:8080/status/state_json/root.json | python3 -I -c 'import json,sys; print(json.loads(json.load(sys.stdin)["fTitle"])["readout"]["events"])' 2>/dev/null || echo -1; }
log "inicio del orquestador"
hv
ok41=0
attempt=0
deadline=$(date -u -d '2026-10-20 06:00:00 UTC' +%s)
while [ $(date -u +%s) -le $deadline ]; do
  attempt=$((attempt + 1))
  u=daq2-6-n41-$attempt
  log "4.1 intento $attempt: arranca $u"
  startrun $u $C/step4a_night_ext_B_hv_on.json 2400
  sleep 180
  ev=$(events)
  rate=$(python3 -c "print($ev/180.0)")
  log "4.1 intento $attempt: $ev eventos en 180 s (tasa $rate Hz)"
  if python3 -c "import sys; sys.exit(0 if $ev >= 18 else 1)"; then
    log "4.1 intento $attempt: tasa suficiente, el run sigue hasta su fin (30 min)"
    while systemctl --user is-active -q $u; do sleep 20; done
    log "4.1 intento $attempt: terminado; $(grep 'daqrun termina' /data/daq/logs/$u.log | cut -c1-200)"
    ok41=1
    break
  fi
  log "4.1 intento $attempt: tasa < 0.1 Hz, se para el run"
  systemctl --user stop $u
  while systemctl --user is-active -q $u; do sleep 2; done
  next=$(( $(date -u +%s) / 1800 * 1800 + 1800 ))
  [ $next -gt $deadline ] && break
  log "4.1: próximo intento a las $(date -u -d @$next +%FT%TZ)"
  sleep $(( next - $(date -u +%s) ))
done
[ $ok41 = 0 ] && log "4.1: ningún intento con tasa suficiente hasta las 06:00Z"
hv
log "4.2: arranca"
startrun daq2-6-n42 $C/step4b_night_sw_full_hv_on.json 600
sleep 10
while systemctl --user is-active -q daq2-6-n42; do sleep 5; done
log "4.2: terminado; $(grep 'daqrun termina' /data/daq/logs/daq2-6-n42.log | cut -c1-200)"
hv
log "fin del orquestador (4.1 ok = $ok41)"
touch /data/daq/logs/night_daq2-6.done
