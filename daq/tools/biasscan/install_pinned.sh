#!/bin/bash
# install_pinned.sh — seguimiento 2026-10-31 (4): copia fija y versionada del barrido de tensión desde un commit dado, y la unidad que la ejecuta.
# NO instala la unidad en systemd ni arma timers: deja el archivo de unidad junto a la copia y dice cómo instalarlo.
#   1. Resuelve <commit> en el repositorio; exige que su biasscan.py tenga BIASSCAN_CODE_ROOT (si no, la copia importaría el código del árbol vivo).
#   2. git archive de los directorios que usa el barrido (daq/tools/biasscan, daq/tools/characterization, daq/tools/validation, daq/core/include,
#      daq/config) a <prefijo>/<commit corto>/; escribe BIASSCAN_COMMIT (hash completo, fecha, origen). Si el destino ya existe y su
#      BIASSCAN_COMMIT es el mismo, no lo toca; si es otro, se niega.
#   3. Archivos de código de solo lectura (los directorios quedan escribibles: ACLiC compila las macros junto al .C).
#   4. Genera <destino>/biasscan_night.service desde daq/systemd/biasscan_night_pinned.service.in del mismo commit, con @CODE@, @CAMPAIGN@,
#      @DEADLINE@ y @EXTRA@ reemplazados.
# Uso: install_pinned.sh <commit> --campaign NOMBRE --deadline "AAAA-MM-DD HH:MM:SS" [--extra "args del barrido"] [--prefix DIR] [--repo DIR]
set -euo pipefail
REPO=/home/xplus/digitizer/CaenDigitizerSiLab-jai
PREFIX=$HOME/.local/share/daqv2/biasscan
EXTRA="--self-check-only --steps 0,0.25,0.5,0.75,1.0 --max-above 1.0 --max-step 0.25 --ext-s 2700 --trigger-diag --rup 0.1"
COMMIT=${1:?uso: install_pinned.sh <commit> --campaign NOMBRE --deadline "AAAA-MM-DD HH:MM:SS" [...]}
shift
CAMPAIGN="" DEADLINE=""
while [ $# -gt 0 ]; do
  case $1 in
    --campaign) CAMPAIGN=$2; shift 2 ;;
    --deadline) DEADLINE=$2; shift 2 ;;
    --extra) EXTRA=$2; shift 2 ;;
    --prefix) PREFIX=$2; shift 2 ;;
    --repo) REPO=$2; shift 2 ;;
    *) echo "argumento desconocido: $1" >&2; exit 2 ;;
  esac
done
[ -n "$CAMPAIGN" ] && [ -n "$DEADLINE" ] || { echo "faltan --campaign y --deadline" >&2; exit 2; }
FULL=$(git -C "$REPO" rev-parse --verify "$COMMIT^{commit}")
SHORT=$(git -C "$REPO" rev-parse --short "$FULL")
git -C "$REPO" show "$FULL:daq/tools/biasscan/biasscan.py" | grep -q BIASSCAN_CODE_ROOT \
  || { echo "el biasscan.py de $SHORT no tiene BIASSCAN_CODE_ROOT: la copia importaría el código del árbol vivo; usar un commit posterior" >&2; exit 3; }
git -C "$REPO" cat-file -e "$FULL:daq/systemd/biasscan_night_pinned.service.in" \
  || { echo "el commit $SHORT no tiene daq/systemd/biasscan_night_pinned.service.in" >&2; exit 3; }
DEST=$PREFIX/$SHORT
if [ -e "$DEST" ]; then
  if [ "$(head -1 "$DEST/BIASSCAN_COMMIT" 2>/dev/null | cut -d' ' -f1)" = "$FULL" ]; then
    echo "ya instalado: $DEST (commit $FULL); no se toca"
  else
    echo "$DEST existe y no es el commit $FULL: no se sobrescribe" >&2; exit 4
  fi
else
  mkdir -p "$DEST"
  git -C "$REPO" archive "$FULL" daq/tools/biasscan daq/tools/characterization daq/tools/validation daq/core/include daq/config | tar -x -C "$DEST"
  echo "$FULL $(git -C "$REPO" show -s --format=%cI "$FULL") $REPO" > "$DEST/BIASSCAN_COMMIT"
  find "$DEST" -type f \( -name '*.py' -o -name '*.C' -o -name '*.h' -o -name '*.sh' -o -name '*.json' -o -name 'BIASSCAN_COMMIT' \) -exec chmod a-w {} +
  echo "copia instalada: $DEST"
fi
UNIT=$DEST/biasscan_night.service
chmod u+w "$DEST" 2>/dev/null || true
git -C "$REPO" show "$FULL:daq/systemd/biasscan_night_pinned.service.in" \
  | sed -e "s|@CODE@|$DEST|g" -e "s|@CAMPAIGN@|$CAMPAIGN|g" -e "s|@DEADLINE@|$DEADLINE|g" -e "s|@EXTRA@|$EXTRA|g" -e "s|@COMMIT@|$FULL|g" > "$UNIT"
echo "unidad generada (NO instalada): $UNIT"
echo "para instalarla (René): cp $UNIT ~/.config/systemd/user/biasscan_night.service && systemctl --user daemon-reload"
echo "  (el timer ~/.config/systemd/user/biasscan_night.timer necesita su OnCalendar; armarlo es otra decisión: systemctl --user start biasscan_night.timer)"
