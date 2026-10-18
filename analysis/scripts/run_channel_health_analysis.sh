#!/usr/bin/env bash
set -euo pipefail

: "${HODO_RAW_DIR:=/home/xplus/digitizer/data/raw}"

cd "$(dirname "$0")/../.."

OUTDIR="outputs/channel_health_20261016"

for input in \
  lectura_vx1740_ext_nim_1M.root \
  lectura_vx1740_ext_nim_1M_v2.root \
  lectura_vx1740_ext_nim_1M_v3.root \
  lectura_vx1740_ext_nim_1M_v4.root \
  lectura_vx1740_ext_nim_1M_v5.root \
  lectura_vx1740_ext_nim_1M_v6.root \
  lectura_vx1740_ext_nim_1M_v7.root \
  lectura_vx1740_ext_nim_1M_v8.root \
  lectura_vx1740_ext_nim_1M_v9.root; do
  if [[ ! -f "$input" && ! -f "$HODO_RAW_DIR/$input" && ! -f "$HODO_RAW_DIR/hodoscope_extnim_2026/$input" ]]; then
    echo "Missing required input: $input" >&2
    exit 1
  fi
done

make -C analysis bin/channel_health_analysis
analysis/bin/channel_health_analysis "$OUTDIR"
