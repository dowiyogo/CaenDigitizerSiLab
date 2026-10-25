#!/bin/bash
# DAQ2-6: arranca el run 8 (trigger externo, ventana B, HV encendida, misma config que el run 7 más
# el aviso "sin triggers") cuando termina el run 7 (unidad daq2-6-s4a2); parada planeada 19:30Z.
set -u
L=/data/daq/logs/run8_orchestrator.log
BIN=$HOME/.local/opt/daqv2/bin
CFG=/data/daq/config/commissioning_2026_10/step4a_run8_ext_B_hv_on.json
log() { echo "$(date -u +%FT%TZ) $*" >> "$L"; }
log "orquestador del run 8: espera el fin del run 7 (daq2-6-s4a2)"
while systemctl --user is-active -q daq2-6-s4a2; do sleep 10; done
log "run 7 terminado: $(grep 'daqrun termina' /data/daq/logs/daq2-6-s4a2.log | cut -c1-160)"
stop=$(date -u -d '2026-10-19 19:30:00 UTC' +%s)
dur=$(( stop - $(date -u +%s) - 15 ))
if [ $dur -lt 300 ]; then log "quedan menos de 5 min hasta las 19:30Z: no se arranca el run 8"; exit 0; fi
log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"
log "arranca el run 8 con --max-duration $dur (fin ~19:30Z)"
systemd-run --user --unit=daq2-6-run8 --collect --property=RuntimeMaxSec=$(( dur + 300 )) \
  --property=StandardOutput=append:/data/daq/logs/daq2-6-run8.log --property=StandardError=append:/data/daq/logs/daq2-6-run8.log \
  --setenv=LD_LIBRARY_PATH=$HOME/root6/lib --working-directory=$HOME $BIN/daqrun --config "$CFG" --max-duration $dur >> "$L" 2>&1
systemd-run --user --unit=daq2-6-run8-sampler --collect --property=RuntimeMaxSec=$(( dur + 900 )) \
  /data/daq/logs/sampler.sh daq2-6-run8 60 /data/daq/logs/daq2-6-run8_ps.log >> "$L" 2>&1
sleep 20
while systemctl --user is-active -q daq2-6-run8; do sleep 30; done
log "run 8 terminado: $(grep 'daqrun termina' /data/daq/logs/daq2-6-run8.log | cut -c1-300)"
log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"
touch /data/daq/logs/run8.done
