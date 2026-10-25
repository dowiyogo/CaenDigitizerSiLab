#!/usr/bin/python3 -I
# split_paths.py — XPLUS-2 B.5 (solo lectura): listas de rutas para `git filter-repo --paths-from-file` a partir de
# docs/split_map_2026-10-25.csv. No ejecuta filter-repo; René lo corre en un clon nuevo.
#  - daqv2: clase DAQ, más los archivos de ANALYSIS que el DAQ usa y que se mudan a daqv2 (MOVE_TO_DAQV2), más LEGACY y OTHER
#    que viven bajo daq/ (irán a attic/ y tools/host/ de daqv2);
#  - daq_analysis: clase ANALYSIS (sin MOVE_TO_DAQV2) y DOCS (sin docs/monorepo_handoff/, copia de documentos del monorepo);
#  - RESULTS no entra a ningún git (va a data/results/ con manifest); el LEGACY de la raíz queda en el repositorio antiguo.
# Cada lista incluye los nombres anteriores de cada archivo (renombres de `git log --all -M`), para no perder la historia de
# antes de la reorganización del 2026-10-17. Avisa si un nombre antiguo lo reclaman las dos listas.
# Uso: split_paths.py <repo> <split_map.csv> <dir de salida>
import collections, csv, os, subprocess, sys

REPO, MAP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MOVE_TO_DAQV2 = {
    "daq/tools/characterization/xy_estimator.py",       # importado por biasscan.py (estimador X·Y único, DAQ2-16)
    "daq/tools/characterization/test_xy_estimator.py",  # lo corre daq/tests/CMakeLists.txt
    "daq/tools/characterization/hitdef_extract.C",      # definición única de impacto; biasscan.py la ejecuta
    "daq/tools/biasscan/biasscan_report.py",            # biasscan.py la ejecuta al final del barrido
    "daq/tools/biasscan/run_info_dump.C",               # biasscan.py la ejecuta por run
    "analysis/src/hodo_paths.h",                        # incluido por daq_dig_replay y daq_root_tests
}
rows = list(csv.DictReader(open(MAP)))
target = {}
for r in rows:
    p, c = r["path"], r["class"]
    if p in MOVE_TO_DAQV2 or c == "DAQ" or (p.startswith("daq/") and c in ("LEGACY", "OTHER")):
        target[p] = "daqv2"
    elif (c == "ANALYSIS") or (c == "DOCS" and not p.startswith("docs/monorepo_handoff/")):
        target[p] = "daq_analysis"

# historia de renombres: nombre nuevo -> nombres anteriores (cadena completa)
log = subprocess.run(["git", "-C", REPO, "--no-optional-locks", "log", "--all", "-M", "--name-status", "--format="],
                     capture_output=True, text=True, check=True).stdout
prev = collections.defaultdict(set)
for line in log.splitlines():
    f = line.split("\t")
    if f and f[0].startswith("R") and len(f) == 3:
        prev[f[2]].add(f[1])


def ancestors(p):
    seen, todo = set(), [p]
    while todo:
        for q in prev.get(todo.pop(), ()):
            if q not in seen:
                seen.add(q)
                todo.append(q)
    return seen


lists = collections.defaultdict(set)
claim = collections.defaultdict(set)
for p, t in target.items():
    for q in {p} | ancestors(p):
        lists[t].add(q)
        claim[q].add(t)
os.makedirs(OUT, exist_ok=True)
for t, ps in lists.items():
    with open(os.path.join(OUT, f"{t}.paths"), "w") as f:
        f.write("".join(q + "\n" for q in sorted(ps)))
    n_now = sum(1 for p in target if target[p] == t)
    print(f"{t}: {n_now} rutas actuales, {len(ps) - n_now} nombres anteriores, {len(ps)} líneas -> {OUT}/{t}.paths")
both = sorted(q for q, ts in claim.items() if len(ts) > 1)
print("nombres reclamados por las dos listas:", both if both else "ninguno")
left = collections.Counter(r["class"] for r in rows if r["path"] not in target)
print("rutas actuales que no van a ningún repo nuevo:", dict(left))
