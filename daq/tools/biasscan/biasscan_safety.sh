#!/bin/bash
# biasscan_safety.sh — DAQ2-14: ExecStopPost de la unidad del barrido. Pase lo que pase con el
# orquestador (fin normal, aborto, señal, RuntimeMaxSec), restaura el V0Set nominal (noop si la sesión
# ya está cerrada) y deja corriendo el run de física con trigger externo.
BIN=$HOME/.local/opt/daqv2/bin
LOG=/data/daq/logs/biasscan.log
echo "$(date -u +%FT%TZ) safety: bias restore" >> $LOG
$BIN/hvctl bias restore --authorize "DAQ2-14 bias scan and self-trigger experiment, authorized by René: restauración de seguridad" >> $LOG 2>&1
$BIN/hvctl bias status >> $LOG 2>&1
if ! systemctl --user is-active --quiet daq-physics; then
  echo "$(date -u +%FT%TZ) safety: run de física no activo: se arranca" >> $LOG
  systemd-run --user --unit=daq-physics --collect --property=RuntimeMaxSec=1209600 --property=KillSignal=SIGTERM --property=TimeoutStopSec=300 \
    --property=StandardOutput=append:/data/daq/logs/daq-physics.log --property=StandardError=append:/data/daq/logs/daq-physics.log \
    --setenv=LD_LIBRARY_PATH=$HOME/root6/lib --working-directory=$HOME $BIN/daqrun --config /data/daq/config/hodoscope_2026_10/physics_ext_B.json >> $LOG 2>&1
fi
