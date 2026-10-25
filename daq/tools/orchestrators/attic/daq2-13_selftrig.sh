#!/bin/bash
# DAQ2-13: tres runs de 5 min con auto-trigger de hod2018 (T = 12, 10, 8 p.e.), uno tras otro.
export LD_LIBRARY_PATH=$HOME/root6/lib
for T in 12 10 08; do
  echo "=== $(date -u +%FT%TZ) inicio T=$T"
  /usr/bin/time -v $HOME/.local/opt/daqv2/bin/daqrun --config /data/daq/config/selftrig_test_2026_10_21/self_T$T.json --max-duration 300
  echo "=== $(date -u +%FT%TZ) fin T=$T rc=$?"
  sleep 5
done
touch /data/daq/logs/daq2-13_selftrig.done
