#!/usr/bin/env bash
set -euo pipefail

OUTDIR="outputs/track_reconstruction_interleaved_v7_corrected_mapping"
INPUT="candidate_events_windowed.csv"

if [[ ! -f "${INPUT}" ]]; then
  echo "Primero genera candidate_events_windowed.csv con scan_active_channels_windowed sobre el run v7."
  exit 1
fi

mkdir -p "${OUTDIR}"

g++ -O2 -std=c++17 reconstruct_tracks_interleaved.cpp $(root-config --cflags --libs) -o reconstruct_tracks_interleaved

./reconstruct_tracks_interleaved \
  candidate_events_windowed.csv \
  "${OUTDIR}" \
  33.0 \
  98.4

cat "${OUTDIR}/summary.txt"
ls -lh "${OUTDIR}"/*.png
