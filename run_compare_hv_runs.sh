#!/usr/bin/env bash
set -euo pipefail

OUTDIR="outputs/compare_hv_runs"

mkdir -p "${OUTDIR}"

g++ -O2 -std=c++17 compare_hv_runs.cpp $(root-config --cflags --libs) -o compare_hv_runs

./compare_hv_runs "${OUTDIR}" \
  lectura_vx1740_ext_nim_1M.root \
  lectura_vx1740_ext_nim_1M_v2.root \
  lectura_vx1740_ext_nim_1M_v3.root \
  lectura_vx1740_ext_nim_1M_v4.root

cat outputs/compare_hv_runs/summary.txt
echo "PNGs en outputs/compare_hv_runs/png/"
