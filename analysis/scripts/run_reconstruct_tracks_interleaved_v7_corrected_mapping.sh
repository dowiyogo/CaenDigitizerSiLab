#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

OUTDIR="outputs/track_reconstruction_interleaved_v7_corrected_mapping"
INPUT="outputs/v7_windowed/candidate_events_windowed.csv"

if [[ ! -f "${INPUT}" ]]; then
  echo "Primero genera candidate_events_windowed.csv con scan_active_channels_windowed sobre el run v7."
  exit 1
fi

mkdir -p "${OUTDIR}"

make -C analysis bin/reconstruct_tracks_interleaved

analysis/bin/reconstruct_tracks_interleaved \
  "${INPUT}" \
  "${OUTDIR}" \
  33.0 \
  98.4

cat "${OUTDIR}/summary.txt"
ls -lh "${OUTDIR}"/*.png
