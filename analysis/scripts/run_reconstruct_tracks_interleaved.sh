#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

OUTDIR="outputs/track_reconstruction_interleaved"

if [[ ! -f outputs/v7_windowed/candidate_events_windowed.csv ]]; then
  echo "Primero ejecuta scan_active_channels_windowed para generar candidate_events_windowed.csv"
  exit 1
fi

mkdir -p "${OUTDIR}"

make -C analysis bin/reconstruct_tracks_interleaved

analysis/bin/reconstruct_tracks_interleaved outputs/v7_windowed/candidate_events_windowed.csv "${OUTDIR}" 33.0 98.4

cat outputs/track_reconstruction_interleaved/summary.txt
ls -lh outputs/track_reconstruction_interleaved/*.png
