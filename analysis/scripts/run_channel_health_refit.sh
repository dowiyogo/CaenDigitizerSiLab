#!/usr/bin/env bash
set -euo pipefail

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
  [[ -f "$input" ]] || { echo "Missing required input: $input" >&2; exit 1; }
done

make -C analysis bin/channel_health_refit
analysis/bin/channel_health_refit "$OUTDIR"
cat "$OUTDIR/summary_channel_health.txt"
