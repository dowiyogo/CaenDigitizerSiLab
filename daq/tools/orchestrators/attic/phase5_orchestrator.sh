#!/bin/bash
# DAQ2-9: run de prueba de la fase 5 (60 min, trigger externo, ventana B, HV encendida) en
# phase5_test_2026_10_20; muestreo de CPU/RSS (60 s) y de disco (600 s). La HV no se toca.
set -u
L=/data/daq/logs/phase5_orchestrator.log
BIN=$HOME/.local/opt/daqv2/bin
CFG=/data/daq/config/phase5_test_2026_10_20/phase5_ext_B_hv_on.json
log() { echo "$(date -u +%FT%TZ) $*" >> "$L"; }
log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"
systemd-run --user --unit=daq2-9-phase5 --collect --property=RuntimeMaxSec=4200 \
  --property=StandardOutput=append:/data/daq/logs/daq2-9-phase5.log --property=StandardError=append:/data/daq/logs/daq2-9-phase5.log \
  --setenv=LD_LIBRARY_PATH=$HOME/root6/lib --working-directory=$HOME $BIN/daqrun --config "$CFG" --max-duration 3600 >> "$L" 2>&1
sleep 2
systemd-run --user --unit=daq2-9-phase5-sampler --collect --property=RuntimeMaxSec=4800 \
  /data/daq/logs/sampler.sh daq2-9-phase5 60 /data/daq/logs/daq2-9-phase5_ps.log >> "$L" 2>&1
systemd-run --user --unit=daq2-9-phase5-disk --collect --property=RuntimeMaxSec=4800 \
  /data/daq/logs/disk_sampler.sh daq2-9-phase5 600 /data/daq/raw/phase5_test_2026_10_20 /data/daq/logs/daq2-9-phase5_disk.log >> "$L" 2>&1
sleep 20
while systemctl --user is-active -q daq2-9-phase5; do sleep 30; done
log "run terminado: $(grep 'daqrun termina' /data/daq/logs/daq2-9-phase5.log | cut -c1-400)"
log "HV: $($BIN/hvctl status | grep -E '^0[0-7]' | awk '{print $1"="$4"V/"$5"A/Pw"$6}' | tr '\n' ' ')"
touch /data/daq/logs/phase5.done
