#!/usr/bin/python3 -I
# split_map.py — XPLUS-2 B (solo lectura, no mueve nada): mapa archivo por archivo para dividir el repositorio en daqv2 y
# daq_analysis/ del monorepo. Cada ruta versionada en HEAD recibe una sola clase (DAQ, ANALYSIS, DOCS, RESULTS, LEGACY, OTHER)
# por la primera regla de RULES que coincide, con su motivo.
# Usuarios de cada archivo (quién lo usa), solo desde archivos de código/configuración/unidades (no .md ni datos):
#  - #include "x" / <x> resueltos a rutas versionadas (con '/': por sufijo; sin '/': mismo directorio, si no, nombre único);
#  - import/from de Python resueltos por nombre de módulo a un .py versionado de nombre único;
#  - el nombre del archivo citado como texto fuera de comentarios (rutas en sys.path, configuraciones, ExecStart, llamadas por ruta), solo para
#    nombres distintivos (≥ 8 caracteres, no genéricos);
#  - un encabezado usado arrastra el .cpp del mismo nombre de su biblioteca (p. ej. sc_reader.h → sc_reader.cpp).
# Los usuarios de clase OTHER (herramientas de inventario, que citan nombres como datos) no cuentan como uso.
# Salida: <csv> con path, class, reason, size_bytes, used_by_other_class, other_class_users, n_same_class_users; y un resumen
# por clase en la salida estándar. Uso: split_map.py <repo> <csv>
import collections, csv, os, re, subprocess, sys

REPO, OUT = sys.argv[1], sys.argv[2]

RULES = [
    (r"^daq/tools/orchestrators/", "DAQ", "registro de los orquestadores de DAQ2-6..13 (copias de /data/daq/logs, XPLUS-2 C); va al attic de daqv2"),
    (r"^docs/split_", "OTHER", "plan de división (XPLUS-2 B): mapa y listas para filter-repo; no entra a ningún repo nuevo"),
    (r"^daq/tools/inventory/", "OTHER", "herramienta de inventario y de plan de división del host (XPLUS-1/2), no es DAQ ni análisis"),
    (r"^daq/tools/legacy/(hv_readall\.c|test_sy5527\.c|CMakeLists\.txt)$", "LEGACY", "herramienta HV heredada (2026-10-18), reemplazada por hvctl read/params/map; propuesta: attic/ en la división"),
    (r"^daq/tools/legacy/test_acquire_sw_20260511\.cxx$", "LEGACY", "programa de adquisición heredado, solo referencia (no se compila)"),
    (r"^daq/tools/legacy/diag_sy5527\.sh$", "DAQ", "diagnóstico de red del SY5527, sin reemplazo en hvctl (su paso 9 compila test_sy5527.c)"),
    (r"^daq/tools/characterization/", "ANALYSIS", "caracterización: lee crudos/salidas y produce números (bench MC, hitdef, xy_estimator, ganancia)"),
    (r"^daq/tools/proceeding/", "ANALYSIS", "paquete del proceeding (PROC-1..4, XPLUS-1 C): lee crudos/salidas y produce figuras y números"),
    (r"^daq/tools/validation/(commissioning_|selftrig_)", "ANALYSIS", "puesta en marcha: lee runs y slow control y produce números (pedestal, p.e., tasas)"),
    (r"^daq/tools/validation/validate_", "DAQ", "validación del software DAQ (DAQ2-5) contra salidas conocidas; prueba de aceptación del DAQ"),
    (r"^daq/tools/biasscan/(biasscan\.py|biasscan_safety\.sh|morning_report_waiter\.py)$", "DAQ", "orquestador del barrido de tensión (hvctl bias, daqrun, systemd)"),
    (r"^daq/tools/biasscan/", "ANALYSIS", "análisis del barrido de tensión desde biasscan_results.json o crudos"),
    (r"^daq/tools/(recover_after_poweroff|test_recover_after_poweroff)\.py$", "DAQ", "recuperación tras corte de energía (XPLUS-1 B) y su prueba contra los fakes"),
    (r"^daq/config/mapping/", "DAQ", "mapeo de canales que daqrun/hvmon leen en tiempo de run"),
    (r"^daq/config/calibration/", "DAQ", "calibración p.e. que daqrun/biasscan leen en tiempo de run"),
    (r"^daq/config/", "DAQ", "configuración del DAQ v2 (ejemplos, guiones de los fakes, referencia HV)"),
    (r"^daq/systemd/", "DAQ", "unidades de usuario de systemd e instalador"),
    (r"^daq/tests/", "DAQ", "pruebas GoogleTest del DAQ v2 (contra los fakes)"),
    (r"^daq/core/include/daq/(pulse_finder|gain_fit|trigger_interval_stats|ttt)\.h$", "DAQ", "encabezado compartido del núcleo (header-only, sin ROOT), usado por daqrun y por el análisis"),
    (r"^daq/io/", "DAQ", "E/S ROOT del DAQ (incluye libdaq_screader, lector de slow control)"),
    (r"^daq/", "DAQ", "DAQ v2 (bibliotecas, apps, CMake)"),
    (r"^analysis/baseline_paths/", "RESULTS", "salidas de análisis versionadas (CSV de referencia de rutas, 2026-10-17)"),
    (r"^analysis/", "ANALYSIS", "análisis C++/ROOT y scripts de v1-v10"),
    (r"^runs/", "ANALYSIS", "catálogo de runs y tiempo bueno de mayo, leído por las macros y por stage1"),
    (r"^outputs/", "RESULTS", "resultado versionado (el resto de outputs/ está ignorado)"),
    (r"^docs/daq_v2_design\.md$", "DAQ", "diseño del DAQ v2; acompaña al código en daqv2"),
    (r"^docs/Manual_de_Uso", "LEGACY", "manual de la biblioteca antigua CaenDigitizerSiLab"),
    (r"^docs/", "DOCS", "handoff, notas, cambios de análisis o reorganización"),
    (r"^attic/", "LEGACY", "attic: código o copias retiradas"),
    (r"^CLAUDE\.md$", "DAQ", "reglas del host y del DAQ v2; al dividir, sus secciones de análisis pasan a DOCS de daq_analysis/"),
    (r"^AGENTS\.md$", "DOCS", "registro del DAQ antiguo, de los runs v1-v10 y del análisis"),
    (r"^(\.gitignore|\.directory)$", "OTHER", "metadatos del repositorio / del gestor de archivos"),
    (r"^(config_commissioning/|config[^/]*\.conf$|CaenDigitizerSiLab\.|LinkDef\.h$|Makefile$|configuration\.|example\.cxx$|test_[^/]*\.cxx$|run_hod2019_[^/]*\.sh$|README\.md$)",
     "LEGACY", "DAQ antiguo (biblioteca CaenDigitizerSiLab, .conf, programas de prueba), no usado por el DAQ v2"),
]
RX = [(re.compile(p), c, r) for p, c, r in RULES]
CODE_EXT = (".c", ".cpp", ".cxx", ".h", ".hpp", ".C", ".py", ".sh", ".cmake", ".json", ".conf", ".service", ".timer", ".example")
GENERIC = {"CMakeLists.txt", "README.md", "main.cpp", "Makefile", "__init__.py", "config.h", "config.cpp", "LinkDef.h", "config.conf"}
DISTINCT_EXT = (".py", ".C", ".sh", ".json", ".csv", ".conf", ".cxx", ".c", ".service", ".timer", ".example", ".h", ".cpp")


def git(*a):
    return subprocess.run(["git", "-C", REPO, "--no-optional-locks", *a], capture_output=True, text=True, check=True).stdout


size = {}
for line in git("ls-tree", "-r", "-l", "HEAD").splitlines():
    meta, path = line.split("\t", 1)
    size[path] = int(meta.split()[3])
paths = sorted(size)
cls, why = {}, {}
for p in paths:
    for rx, c, r in RX:
        if rx.search(p):
            cls[p], why[p] = c, r
            break
    else:
        cls[p], why[p] = "OTHER", "sin regla"

by_base = collections.defaultdict(list)
for p in paths:
    by_base[os.path.basename(p)].append(p)
py_mod = {os.path.basename(p)[:-3]: p for p in paths if p.endswith(".py") and len(by_base[os.path.basename(p)]) == 1}
distinct = [b for b, ps in by_base.items() if len(ps) == 1 and len(b) >= 8 and b not in GENERIC and b.endswith(DISTINCT_EXT)]
inc_rx = re.compile(r'#include\s*[<"]([^>"]+)[>"]')
imp_rx = re.compile(r"^\s*(?:from\s+([A-Za-z_]\w*)|import\s+([A-Za-z_][\w ,.]*?))(?:\s+import|\s*#|$)", re.M)

users = collections.defaultdict(set)


def strip_comments(u, text):
    """quita comentarios antes de buscar nombres citados (un comentario que nombra un archivo no es un uso)."""
    if u.endswith((".py", ".sh", ".cmake", ".service", ".timer", ".example")) or os.path.basename(u) == "CMakeLists.txt":
        return re.sub(r"(?m)(^|\s)#.*$", r"\1", text)
    if u.endswith((".c", ".cpp", ".cxx", ".h", ".hpp", ".C")):
        return re.sub(r"//[^\n]*", "", re.sub(r"/\*.*?\*/", "", text, flags=re.S))
    return text


def add(target, user):
    if target != user and cls[user] != "OTHER":
        users[target].add(user)
        if target.endswith(".h"):  # el .cpp de la misma biblioteca
            stem = os.path.basename(target)[:-2]
            for q in by_base.get(stem + ".cpp", []):
                if q.split("/")[:2] == target.split("/")[:2] and q != user:
                    users[q].add(user)


for u in paths:
    if not (u.endswith(CODE_EXT) or os.path.basename(u) == "CMakeLists.txt"):
        continue
    if size[u] > 2_000_000:
        continue
    try:
        text = open(os.path.join(REPO, u), errors="replace").read()
    except OSError:
        continue
    for inc in inc_rx.findall(text):
        if "/" in inc:
            hits = [p for p in paths if p.endswith("/" + inc) or p == inc]
        else:
            same = os.path.join(os.path.dirname(u), inc)
            hits = [same] if same in size else (by_base[inc] if len(by_base.get(inc, [])) == 1 else [])
        for h in hits:
            add(h, u)
    if u.endswith(".py"):
        for a, b in imp_rx.findall(text):
            for name in ([a] if a else [x.strip().split(" as ")[0].split(".")[0] for x in b.split(",")]):
                if name in py_mod:
                    add(py_mod[name], u)
    code = strip_comments(u, text)
    for b in distinct:
        if u.endswith((".py", ".json")) and b.endswith((".h", ".hpp", ".cpp")):
            continue  # un .py o .json que nombra un encabezado lo cita como procedencia; no lo incluye
        if b in code and re.search(r"(?<![\w.-])" + re.escape(b) + r"(?![\w])", code):
            add(by_base[b][0], u)

rows, other_users = [], {}
for p in paths:
    ou = sorted(u for u in users[p] if cls[u] != cls[p])
    other_users[p] = ou
    rows.append({"path": p, "class": cls[p], "reason": why[p], "size_bytes": size[p], "used_by_other_class": "yes" if ou else "no",
                 "other_class_users": "; ".join(f"{cls[u]}:{u}" for u in ou), "n_same_class_users": sum(1 for u in users[p] if cls[u] == cls[p])})
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

tot = collections.Counter()
n = collections.Counter()
for p in paths:
    n[cls[p]] += 1
    tot[cls[p]] += size[p]
print(f"rutas en HEAD: {len(paths)}  bytes: {sum(size.values())}")
for c in ("DAQ", "ANALYSIS", "DOCS", "RESULTS", "LEGACY", "OTHER"):
    print(f"{c:9s} {n[c]:4d} archivos  {tot[c]:>10d} B")
print("\nusados por otra clase:")
for p in paths:
    if other_users[p]:
        ucls = collections.Counter(cls[u] for u in other_users[p])
        print(f"  {cls[p]:8s} {p}  <- {dict(ucls)}")
