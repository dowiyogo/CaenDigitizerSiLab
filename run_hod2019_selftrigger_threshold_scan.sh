#!/usr/bin/env bash

set -u
set -o pipefail

BASE_CONFIG="${1:-config_vx1740_autotrigger_hod2019_ch32_63.conf}"
OUTDIR="outputs/hod2019_threshold_scan"
EXEC="./test_acquire_sw"
SCAN_BIN="./scan_active_channels"

# Keep runs short by default for diagnosis.
NEVENTS="${NEVENTS:-2000}"
TIMEOUT="${TIMEOUT:-45}"
MAX_ROOT_MB="${MAX_ROOT_MB:-512}"

THRESHOLDS=(-80 -40 -20 -10 -5 -2)

stop_requested=0
trap 'echo; echo "[INFO] Ctrl+C detectado. Se detendra el barrido al finalizar la corrida actual."; stop_requested=1' INT

if [[ ! -f "$BASE_CONFIG" ]]; then
  echo "[ERROR] No existe el archivo base: $BASE_CONFIG"
  exit 1
fi

if [[ ! -x "$EXEC" ]]; then
  echo "[ERROR] No se encontro ejecutable: $EXEC"
  echo "        Compila con: make test_acquire_sw"
  exit 1
fi

mkdir -p "$OUTDIR"

echo "[INFO] Base config: $BASE_CONFIG"
echo "[INFO] Outdir: $OUTDIR"
echo "[INFO] Thresholds: ${THRESHOLDS[*]} mV"
echo "[INFO] nevents=$NEVENTS timeout=$TIMEOUT s MAX_ROOT_MB=$MAX_ROOT_MB"

if [[ ! -x "$SCAN_BIN" ]]; then
  echo "[WARN] No se encontro $SCAN_BIN"
  echo "[WARN] Para postprocesar automaticamente, compila scan_active_channels y vuelve a correr el script."
fi

for thr in "${THRESHOLDS[@]}"; do
  if [[ $stop_requested -eq 1 ]]; then
    echo "[INFO] Barrido detenido por usuario."
    break
  fi

  if [[ "$thr" == -* ]]; then
    tag="m${thr#-}"
  else
    tag="p${thr}"
  fi

  cfg_tmp="$OUTDIR/config_tmp_thr_${tag}.conf"
  root_out="$OUTDIR/hod2019_thr_${tag}.root"

  awk -v thr="$thr" -v nevt="$NEVENTS" -v tout="$TIMEOUT" '
    BEGIN { seen_thr=0; seen_nev=0; seen_tout=0 }
    {
      if ($1 == "vthreshold") {
        print "vthreshold " thr
        seen_thr=1
        next
      }
      if ($1 == "nevents") {
        print "nevents " nevt
        seen_nev=1
        next
      }
      if ($1 == "timeout") {
        print "timeout " tout
        seen_tout=1
        next
      }
      print
    }
    END {
      if (!seen_thr) print "vthreshold " thr
      if (!seen_nev) print "nevents " nevt
      if (!seen_tout) print "timeout " tout
    }
  ' "$BASE_CONFIG" > "$cfg_tmp"

  echo
  echo "[RUN] Threshold ${thr} mV"
  echo "[RUN] Config: $cfg_tmp"
  echo "[RUN] Output: $root_out"

  "$EXEC" "$cfg_tmp" "$root_out"
  run_rc=$?
  echo "[RUN] Exit code: $run_rc"

  if [[ -f "$root_out" ]]; then
    root_size_mb=$(du -m "$root_out" | awk '{print $1}')
    echo "[RUN] ROOT size: ${root_size_mb} MB"
    if [[ "$root_size_mb" -ge "$MAX_ROOT_MB" ]]; then
      echo "[WARN] ROOT excede MAX_ROOT_MB=$MAX_ROOT_MB. Se detiene el barrido para evitar corrida larga por ruido."
      stop_requested=1
    fi
  else
    echo "[WARN] No se genero archivo ROOT para threshold ${thr} mV"
  fi

  if [[ -x "$SCAN_BIN" && -f "$root_out" ]]; then
    echo "[POST] Ejecutando scan_active_channels"
    "$SCAN_BIN" "$root_out" -1 4 5

    if [[ -f active_channels_summary.csv ]]; then
      cp active_channels_summary.csv "$OUTDIR/active_channels_summary_thr_${tag}.csv"
      echo "[POST] Copiado: $OUTDIR/active_channels_summary_thr_${tag}.csv"
    else
      echo "[POST] No se encontro active_channels_summary.csv"
    fi

    if [[ -f event_multiplicity.csv ]]; then
      cp event_multiplicity.csv "$OUTDIR/event_multiplicity_thr_${tag}.csv"
      echo "[POST] Copiado: $OUTDIR/event_multiplicity_thr_${tag}.csv"
    else
      echo "[POST] No se encontro event_multiplicity.csv"
    fi
  fi

done

echo
echo "[DONE] Barrido finalizado. Resultados en $OUTDIR"
