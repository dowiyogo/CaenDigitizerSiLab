#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

OUTDIR="outputs/final_hod2018_analysis_v4"

mkdir -p "${OUTDIR}"

make -C analysis bin/final_hodoscope_analysis

analysis/bin/final_hodoscope_analysis lectura_vx1740_ext_nim_1M_v4.root "${OUTDIR}"

echo "Revisar outputs/final_hod2018_analysis_v4/summary.txt"
echo "Revisar outputs/final_hod2018_analysis_v4/png/"
