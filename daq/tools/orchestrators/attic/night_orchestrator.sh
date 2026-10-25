#!/bin/bash
# DAQ2-7: run nocturno en night_2026_10_19 (trigger externo, ventana B, HV encendida, config del
# run 8 con ttt_tick_ns medido 8.000071), parada planeada 13:00Z del 2026-10-20; muestreo de RSS
# (60 s) y de disco (600 s). La HV no se toca.
set -u
L=/data/daq/logs/night_orchestrator.log
BIN=$HOME/.local/opt/daqv2/bin
CFG=/data/daq/config/night_2026_10_19/night_ext_B_hv_on.json
log() { echo "$(date -u +%FT%TZ) $*" >> "$L"; }
stop=$(date -u -d '2026-10-20 13:00:00 UTC' +%s)
dur=$(( stop - $(date -u +%s) - 15 ))
log "orquestador nocturno: --max-duration $dur (fin ~13:00Z)"
log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"
systemd-run --user --unit=daq2-7-night --collect --property=RuntimeMaxSec=$(( dur + 600 )) \
  --property=StandardOutput=append:/data/daq/logs/daq2-7-night.log --property=StandardError=append:/data/daq/logs/daq2-7-night.log \
  --setenv=LD_LIBRARY_PATH=$HOME/root6/lib --working-directory=$HOME $BIN/daqrun --config "$CFG" --max-duration $dur >> "$L" 2>&1
sleep 2
systemd-run --user --unit=daq2-7-night-sampler --collect --property=RuntimeMaxSec=$(( dur + 1200 )) \
  /data/daq/logs/sampler.sh daq2-7-night 60 /data/daq/logs/daq2-7-night_ps.log >> "$L" 2>&1
systemd-run --user --unit=daq2-7-night-disk --collect --property=RuntimeMaxSec=$(( dur + 1200 )) \
  /data/daq/logs/disk_sampler.sh daq2-7-night 600 /data/daq/raw/night_2026_10_19 /data/daq/logs/daq2-7-night_disk.log >> "$L" 2>&1
sleep 20
while systemctl --user is-active -q daq2-7-night; do sleep 30; done
log "run terminado: $(grep 'daqrun termina' /data/daq/logs/daq2-7-night.log | cut -c1-400)"
log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"
log "notificador hvmon: $($BIN/hvctl status | grep notificador)"
touch /data/daq/logs/night.done
