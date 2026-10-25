#!/usr/bin/python3 -I
# write_report_a.py — XPLUS-1 parte A: escribe la parte A del informe /home/xplus/xplus_inventory_2026-10-25_report.md desde
# <dir>/inventory.json y las salidas <dir>/cmd_*.txt de host_inventory.py (cada número va con la salida del comando que lo da).
# Uso: write_report_a.py <dir del inventario> <informe>
import json, os, sys

D, REP = sys.argv[1], sys.argv[2]
inv = json.load(open(D + "/inventory.json"))
cmd = lambda n: open(f"{D}/cmd_{n}.txt").read().rstrip()
GB = lambda b: f"{b / 1e9:.3f} GB"
MB = lambda b: f"{b / 1e6:.1f} MB"
L = []
w = L.append

# ---------------------------------------------------------------- clasificación de datos
data = inv["data_dirs"]
def agg(pred):
    a = {"dirs": 0, "files": 0, "bytes": 0, "first": "9", "last": "0", "idx": set()}
    for h, d in data.items():
        if pred(h, d):
            a["dirs"] += 1; a["files"] += d["files"]; a["bytes"] += d["bytes"]
            a["first"] = min(a["first"], d["first"]); a["last"] = max(a["last"], d["last"]); a["idx"] |= set(d["index_files"])
    return a
raw_campaigns = sorted({h for h in data if h.startswith("/data/daq/raw/")})
sim = agg(lambda h, d: h.startswith("/data/daq/sim/"))
old = agg(lambda h, d: h.startswith("/home/xplus/digitizer/raw.old_sdb_20261018"))
camac = agg(lambda h, d: h.startswith("/home/xplus/datos_camac_USM"))
eic = agg(lambda h, d: h == "/home/xplus")
inrepo = agg(lambda h, d: d["in_repo"])
vendor = agg(lambda h, d: h.startswith(("/home/xplus/Downloads", "/home/xplus/CAENGECO", "/home/xplus/.CAEN")))
raw_total = sum(data[h]["bytes"] for h in raw_campaigns)
repos = {r["path"]: r for r in inv["repos"]}
loose = inv["loose_scripts"]

w("# Inventario del host xplus, scripts de recuperación y análisis menores (XPLUS-1, 2026-10-25)\n")
w(f"Inventario generado por `daq/tools/inventory/host_inventory.py` el {inv['generated_utc']} en `{inv['host']}`. Todo es de solo lectura.")
w("Las salidas completas de los comandos están en `/home/xplus/digitizer/CaenDigitizerSiLab-jai/outputs/xplus_inventory_2026-10-25/cmd_*.txt`, y el resumen en `inventory.json`.")
w("Etiquetas: **medido** = salida de un comando sobre el estado actual; **calculado** = suma o cuenta sobre lo medido; **supuesto** = criterio mío, dicho en el texto.\n")

w("## Resumen (una página)\n")
w("| clase | qué | tamaño | estado | etiqueta |")
w("|---|---|---|---|---|")
j = repos["/home/xplus/digitizer/CaenDigitizerSiLab-jai"]
w(f"| reproducible desde un repo | `CaenDigitizerSiLab-jai` (DAQ v2, análisis, proceeding) | .git 45 MB; directorio {MB(j['size_bytes'])} (outputs/ 1.6 GB y daq/build* 3.8 GB fuera de git) | ramas con commits sin empujar: " + ", ".join(f"{b} {n}" for b, n in j["unpushed_per_branch"].items() if n) + f"; {j['uncommitted']} cambios sin commit | medido |")
p = repos["/home/xplus/digitizer/muographers2026-proceeding"]
w(f"| reproducible desde un repo | `muographers2026-proceeding` (figuras y tablas del artículo) | {MB(p['size_bytes'])} | al día con origin, {p['uncommitted']} cambios | medido |")
for k in ("/home/xplus/camac", "/home/xplus/digitizer/CaenDigitizerSiLab-dev", "/home/xplus/digitizer/CaenDigitizerSiLab-mas"):
    r = repos[k]
    w(f"| reproducible desde un repo | `{k}` | {MB(r['size_bytes'])} | último commit {r['last_commit'][:10]}; sin empujar {sum(r['unpushed_per_branch'].values())}; {r['uncommitted']} cambios ({r['uncommitted_untracked']} sin seguimiento) | medido |")
w(f"| terceros / herramientas | ROOT (`~/root/root_src`, `~/root6`), `.codex`, `.gemini/history` (2 repos sin remoto, 1 commit cada uno) | root_src {MB(repos['/home/xplus/root/root_src']['size_bytes'])} | no son del proyecto | medido |")
nl = sum(len(v) for v in loose.values())
daqlogs = len(loose.get("/data/daq/logs", []))
w(f"| suelto (código) | {nl} scripts fuera de repos en {len(loose)} directorios: {daqlogs} orquestadores del DAQ en `/data/daq/logs` (DAQ2-6 a DAQ2-13), 5 de `datos_camac_USM`, el resto instaladores de CAEN en `Downloads` | < 0.4 MB | sin versionar | medido |")
w(f"| datos crudos del digitizador | `/data/daq/raw` ({len(raw_campaigns)} campañas; mayo `hodoscope_extnim_2026` {GB(data['/data/daq/raw/hodoscope_extnim_2026']['bytes'])}) | {GB(raw_total)} | cada campaña tiene `SHA256SUMS`; las de octubre además `*.index.json` | medido |")
w(f"| copia antigua de datos crudos | `/home/xplus/digitizer/raw.old_sdb_20261018` (mayo, en el disco del sistema sdb) | {GB(old['bytes'])} | `SHA256SUMS` por campaña y `raw_manifest_20261018.sha256` | medido |")
w(f"| datos simulados | `/data/daq/sim` | {GB(sim['bytes'])} | regenerables con las herramientas del repo | medido (tamaño); supuesto (regenerables) |")
w(f"| otros datos | `datos_camac_USM` (CAMAC, no digitizador) {GB(camac['bytes'])}; 3 `.root` de EIC en `/home/xplus` (2022) {GB(eic['bytes'] - 0)} | {GB(camac['bytes'] + eic['bytes'])} | sin índice | medido |")
w(f"| control lento y configuración | `/data/daq/slowcontrol` 105 MB, `/data/daq/config` 492 KB, `/data/daq/logs` 584 KB, `/data/daq/manifests` 64 KB | ~0.11 GB | `.gz` de días cerrados | medido |")
w("")
w("Disco: `/data` está en **sda1** (ST500DM002, 465.8 GB; 28 GB usados, 407 GB libres). El sistema y `/home` están en **sdb2** (ST1000DM010, 931 GB; 58 GB usados, 811 GB libres). Medido; salidas en A.4.\n")

w("## A.1 Repositorios git (medido)\n")
w("Criterio: `find /home/xplus -maxdepth 4 -name .git`. \"Sin empujar\" son los commits de cada rama local que no están en ninguna referencia remota local, contados sin hacer `git fetch`, así que reflejan el último fetch. Todo con `git --no-optional-locks`, que no escribe el índice.\n")
w("| ruta | remoto | rama actual | último commit | sin empujar (por rama) | sin commit | tamaño |")
w("|---|---|---|---|---|---|---|")
for r in inv["repos"]:
    rem = r["remotes"][0].split()[1] if r["remotes"] else "(sin remoto)"
    unp = ", ".join(f"{b}: {n}" for b, n in r["unpushed_per_branch"].items() if n) or "0"
    w(f"| `{r['path']}` | {rem} | {r['branch']} | {r['last_commit'][:25]} | {unp} | {r['uncommitted']} ({r['uncommitted_untracked']} sin seguimiento) | {MB(r['size_bytes'])} |")
w("\nSalida de los comandos:\n\n```\n" + cmd("git_repos") + "\n```\n")
w("Tamaños dentro de `CaenDigitizerSiLab-jai` (lo que no está en git):\n\n```\n" + cmd("du_repo") + "\n```\n")

w("## A.2 Scripts de análisis fuera de los repositorios (medido)\n")
w("Criterio: `.C`, `.cxx`, `.py`, `.ipynb` y `.sh` bajo `/home/xplus` y `/data`, sin los 9 repositorios.")
w("Excluí, por ser de terceros o de herramientas (criterio supuesto): " + ", ".join(f"`{e}`" for e in inv["exclude_for_loose_scripts"]) + ".\n")
for d_, v in loose.items():
    w(f"**`{d_}`**\n")
    w("| archivo | tamaño (B) | mtime | primera línea de comentario |")
    w("|---|---|---|---|")
    for x in v:
        w(f"| {x['file']} | {x['size']} | {x['mtime']} | {x['first_comment'].replace('|', '/') or '–'} |")
    w("")
w("Salida del comando (`find` con `-printf`):\n\n```\n" + cmd("find_scripts") + "\n```\n")

w("## A.3 Datos crudos (medido)\n")
w("Criterio: archivos `.root`, `.dat`, `wave*.txt` y `*caen*` bajo `/data` y `/home/xplus`. Los de `Downloads` son instaladores y manuales de CAEN, no datos, y los de los repositorios son salidas de análisis; van aparte. Fechas: mtime de los archivos.\n")
w("| campaña / directorio | archivos | tamaño | primero | último | índice |")
w("|---|---|---|---|---|---|")
for h in raw_campaigns:
    d = data[h]
    w(f"| `{h}` | {d['files']} | {GB(d['bytes'])} | {d['first']} | {d['last']} | {', '.join(x for x in d['index_files'] if not x.endswith('.index.json')) or '–'}{' + ' + str(sum(1 for x in d['index_files'] if x.endswith('.index.json'))) + ' index.json' if any(x.endswith('.index.json') for x in d['index_files']) else ''} |")
for name, a in (("/data/daq/sim (simulado)", sim), ("/home/xplus/digitizer/raw.old_sdb_20261018 (copia antigua, mayo)", old), ("/home/xplus/datos_camac_USM (CAMAC)", camac),
                ("/home/xplus (.root de EIC, 2022)", eic), ("dentro de repositorios (salidas)", inrepo), ("Downloads y software CAEN (no son datos)", vendor)):
    w(f"| `{name}` | {a['files']} | {GB(a['bytes'])} | {a['first']} | {a['last']} | {', '.join(sorted(x for x in a['idx'] if not x.endswith('.index.json'))[:3]) or '–'} |")
w(f"\nTotal de `/data/daq/raw`: {GB(raw_total)} calculado sobre los archivos que cumplen el criterio; `du -sh /data/daq/raw` da 26G (A.4). `/home/xplus/digitizer/data/raw` es un enlace a `/data/daq/raw`.")
w("La salida completa de `find` está en `cmd_find_data.txt` (una línea por archivo).\n")
w("```\n$ head -5 cmd_find_data.txt\n" + "\n".join(cmd("find_data").splitlines()[:6]) + "\n```\n")

w("## A.4 Discos (medido)\n")
for n in ("df", "lsblk", "mount", "du_data"):
    w("```\n" + cmd(n) + "\n```\n")

w("## A.5 Tareas programadas, sesiones y procesos (medido)\n")
w("- **Unidades de usuario del DAQ:**")
w("  - `hvmon.service`: enabled y activa;")
w("  - `biasscan.service`, `biasscan_day.service`, `biasscan_night.service`: static e inactivas;")
w("  - `biasscan_night.timer`: static, \"elapsed\", sin disparo pendiente.")
w("- **Unidades de sistema:** ninguna del DAQ; solo `snap-caengeco2020-x1.mount`, la de la snap de GECO.")
w("- **Cron:** no hay crontab de xplus, y en `/etc/cron.d` no hay nada del DAQ.")
w("- **screen y tmux:** no están instalados.")
w("- **Procesos del DAQ:** solo `hvmon`, el PID 1221879, desde 2026-10-23 15:41 -03.\n")
for n in ("systemd_user_units", "systemd_user_timers", "user_unit_files", "systemd_system", "cron", "screen_tmux", "processes"):
    w("```\n" + cmd(n) + "\n```\n")
open(REP, "w").write("\n".join(L) + "\n")
print("informe parte A:", REP, len(L), "líneas")
