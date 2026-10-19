#!/bin/bash
# Instala los binarios de daq/build en ~/.local/opt/daqv2 (ruta estable: recompilar no toca un
# binario en ejecución) y la unidad de usuario de hvmon. No arranca nada.
#   daq/systemd/install.sh            # binarios + mapeo + unidad (y daemon-reload)
# La configuración ~/.config/daqv2/hvmon.json y secrets.env no se tocan si ya existen.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
prefix="$HOME/.local/opt/daqv2"
mkdir -p "$prefix/bin" "$prefix/share/mapping" "$HOME/.config/systemd/user"
for b in apps/hvmon/hvmon apps/hvctl/hvctl apps/sctool/sctool apps/daqrun/daqrun apps/daqcampaign/daqcampaign apps/digprobe/digprobe; do
  src="$repo/daq/build/$b"
  [ -x "$src" ] || { echo "falta $src (compilar daq/build primero)" >&2; exit 1; }
  dst="$prefix/bin/$(basename "$b")"
  install -m 0755 "$src" "$dst.new"   # copia nueva y rename: el proceso en marcha conserva su inodo
  mv -f "$dst.new" "$dst"
done
install -m 0644 "$repo/daq/config/mapping/hodo_mapping_2026_10.json" "$prefix/share/mapping/"
git -C "$repo" rev-parse HEAD > "$prefix/INSTALLED_COMMIT"
install -m 0644 "$repo/daq/systemd/hvmon.service" "$HOME/.config/systemd/user/hvmon.service"
systemctl --user daemon-reload
echo "instalado en $prefix (commit $(cat "$prefix/INSTALLED_COMMIT"))"
