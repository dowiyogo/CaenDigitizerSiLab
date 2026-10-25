#!/bin/bash
# DAQ2-7: uso de disco cada $2 s mientras la unidad $1 esté activa: df de /data y du de $3; salida en $4.
unit=$1; per=$2; dir=$3; out=$4
echo "# t_utc df_used_kB df_avail_kB du_campaign_kB (unidad $unit, $dir)" > "$out"
sleep 30
while systemctl --user is-active -q "$unit"; do
  echo "$(date -u +%FT%TZ) $(df -k --output=used,avail /data | tail -1) $(du -sk "$dir" 2>/dev/null | cut -f1)" >> "$out"
  sleep "$per"
done
echo "$(date -u +%FT%TZ) $(df -k --output=used,avail /data | tail -1) $(du -sk "$dir" 2>/dev/null | cut -f1) fin" >> "$out"
