#!/usr/bin/env bash
set -euo pipefail

: "${HODO_RAW_DIR:=/home/xplus/digitizer/data/raw}"

cd "$(dirname "$0")/../.."

OUTDIR="outputs/compare_hv_runs"

mkdir -p "${OUTDIR}"

make -C analysis bin/compare_hv_runs

analysis/bin/compare_hv_runs "${OUTDIR}" \
  lectura_vx1740_ext_nim_1M.root \
  lectura_vx1740_ext_nim_1M_v2.root \
  lectura_vx1740_ext_nim_1M_v3.root \
  lectura_vx1740_ext_nim_1M_v4.root

cat outputs/compare_hv_runs/summary.txt
echo "PNGs en outputs/compare_hv_runs/png/"
