#!/bin/bash
# Muestrea CPU y RSS de la unidad systemd $1 cada $2 s hasta que termine; salida en $3.
unit=$1; per=$2; out=$3
echo "# unidad $unit" > "$out"
while systemctl --user is-active -q "$unit"; do
  pid=$(systemctl --user show -p MainPID --value "$unit")
  [ "$pid" != "0" ] && echo "$(date -u +%FT%TZ) $(ps -o etimes=,cputimes=,rss= -p "$pid")" >> "$out"
  sleep "$per"
done
