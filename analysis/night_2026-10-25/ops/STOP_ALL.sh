#!/bin/bash
# STOP_ALL.sh — XPLUS-4: termina los trabajos de esta tarea por PID guardado en pids/<nombre>.pid.
# Uso: STOP_ALL.sh [--dry-run]   (--dry-run lista los PID y si siguen vivos, sin matar nada)
# Solo envía señales a PID cuya línea de comando contiene la marca XPLUS4 (no mata un PID reciclado por otro proceso).
D=/home/xplus/night_2026-10-25/pids
dry=0; [ "${1:-}" = "--dry-run" ] && dry=1
shopt -s nullglob
files=("$D"/*.pid)
[ ${#files[@]} -eq 0 ] && { echo "no hay PID guardados en $D"; exit 0; }
for f in "${files[@]}"; do
  pid=$(cat "$f"); name=$(basename "$f" .pid)
  if ! kill -0 "$pid" 2>/dev/null; then echo "$name: PID $pid ya terminó"; continue; fi
  cmd=$(tr '\0' ' ' < /proc/$pid/cmdline 2>/dev/null)
  case "$cmd" in *XPLUS4*) ;; *) echo "$name: PID $pid vivo pero sin la marca XPLUS4 ($cmd); no se toca"; continue;; esac
  if [ $dry = 1 ]; then echo "$name: PID $pid vivo -> se mataría (grupo de procesos): $cmd"; continue; fi
  pgid=$(ps -o pgid= -p "$pid" | tr -d ' ')
  echo "$name: TERM al grupo $pgid (PID $pid)"; kill -TERM -- -"$pgid" 2>/dev/null || kill -TERM "$pid"
  sleep 5; kill -0 "$pid" 2>/dev/null && { echo "$name: sigue vivo, KILL"; kill -KILL -- -"$pgid" 2>/dev/null || kill -KILL "$pid"; }
done
