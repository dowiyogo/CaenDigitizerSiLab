#!/usr/bin/env bash
set -euo pipefail

OUTDIR="outputs/final_hod2018_analysis_v4"

mkdir -p "${OUTDIR}"

g++ -O2 -std=c++17 final_hodoscope_analysis.cpp $(root-config --cflags --libs) -o final_hodoscope_analysis

./final_hodoscope_analysis lectura_vx1740_ext_nim_1M_v4.root "${OUTDIR}"

echo "Revisar outputs/final_hod2018_analysis_v4/summary.txt"
echo "Revisar outputs/final_hod2018_analysis_v4/png/"
