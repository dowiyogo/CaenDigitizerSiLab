#!/usr/bin/env bash
set -euo pipefail

OUTDIR="outputs/track_reconstruction_interleaved"

if [[ ! -f candidate_events_windowed.csv ]]; then
  echo "Primero ejecuta scan_active_channels_windowed para generar candidate_events_windowed.csv"
  exit 1
fi

mkdir -p "${OUTDIR}"

g++ -O2 -std=c++17 reconstruct_tracks_interleaved.cpp $(root-config --cflags --libs) -o reconstruct_tracks_interleaved

./reconstruct_tracks_interleaved candidate_events_windowed.csv "${OUTDIR}" 33.0 98.4

cat outputs/track_reconstruction_interleaved/summary.txt
ls -lh outputs/track_reconstruction_interleaved/*.png
