#!/usr/bin/python3 -I
# may_raw_manifest.py — XPLUS-1 parte D: manifiesto de los datos crudos de mayo (v1-v10) sin copiarlos ni escribir en /data/daq/raw.
#  - lee cada .root de <raw> en bloques de 8 MiB (nunca el archivo entero en memoria) y calcula SHA256, tamaño y mtime (UTC);
#  - escribe en <out>: MANIFEST.sha256 (columnas sha256, size_bytes, mtime_utc, name; cabecera con '#'), SHA256SUMS.check
#    (formato de sha256sum -c, para verificar desde <raw>), manifest.json (con los triggers por run de PROC-1 etapa 1,
#    work/may_inventory.jsonl), INDEX_line.md (borrador de la línea de INDEX.md del monorepo) y SHA256SUMS de esos archivos;
#  - compara con el SHA256SUMS que ya existe en <raw>.
# Uso: may_raw_manifest.py <raw> <out> <may_inventory.jsonl>
import datetime, hashlib, json, os, socket, subprocess, sys, time

RAW, OUT, INV = sys.argv[1], sys.argv[2], sys.argv[3]
REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
if os.path.realpath(OUT).startswith(os.path.realpath(RAW)) or os.path.realpath(OUT).startswith("/data/daq/raw"):
    sys.exit("la salida no puede estar dentro de /data/daq/raw")
os.makedirs(OUT, exist_ok=True)


def utc(t):
    return datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(8 << 20)
            if not b:
                return h.hexdigest()
            h.update(b)


files = []
t0 = time.monotonic()
for name in sorted(n for n in os.listdir(RAW) if n.endswith(".root")):
    p = os.path.join(RAW, name)
    st = os.stat(p)
    t = time.monotonic()
    d = sha256(p)
    st2 = os.stat(p)
    if (st2.st_size, st2.st_mtime_ns) != (st.st_size, st.st_mtime_ns):
        sys.exit(f"{name} cambió durante la lectura")
    files.append({"name": name, "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns, "mtime_utc": utc(st.st_mtime_ns / 1e9), "sha256": d})
    print(f"{d}  {st.st_size:>12d}  {files[-1]['mtime_utc']}  {name}  ({time.monotonic() - t:.1f} s)", flush=True)

# comparación con el SHA256SUMS existente
old = {}
for line in open(os.path.join(RAW, "SHA256SUMS")):
    if line.strip():
        h, n = line.split(None, 1)
        old[n.strip().lstrip("*")] = h
cmp = {f["name"]: ("igual" if old.get(f["name"]) == f["sha256"] else ("falta en SHA256SUMS" if f["name"] not in old else "DISTINTO")) for f in files}
extra = sorted(set(old) - {f["name"] for f in files})
print("comparación con SHA256SUMS:", {k: v for k, v in sorted(cmp.items())}, "solo en SHA256SUMS:", extra)

# triggers por run (PROC-1 etapa 1)
inv = {}
for line in open(INV):
    r = json.loads(line)
    inv[os.path.basename(r["file"])] = r
runs = []
for f in files:
    r = inv.get(f["name"])
    f["run"] = r["label"] if r else None
    if r:
        f.update({"entries": r["entries"], "utc_first": utc(r["utc_first_s"]), "utc_last": utc(r["utc_last_s"])})
        runs.append(r)
total = sum(f["size_bytes"] for f in files)
order = sorted(runs, key=lambda r: int(r["label"][1:]))
first, last = min(r["utc_first_s"] for r in runs), max(r["utc_last_s"] for r in runs)

with open(os.path.join(OUT, "MANIFEST.sha256"), "w") as m:
    m.write(f"# MANIFEST.sha256 — datos crudos de mayo 2026 (v1-v10), campaña hodoscope_extnim_2026, host {socket.gethostname()}\n")
    m.write(f"# ruta {os.path.realpath(RAW)}; generado {utc(time.time())} por daq/tools/inventory/may_raw_manifest.py; mtime en UTC\n")
    m.write("# sha256  size_bytes  mtime_utc  name\n")
    for f in files:
        m.write(f"{f['sha256']}  {f['size_bytes']}  {f['mtime_utc']}  {f['name']}\n")
with open(os.path.join(OUT, "SHA256SUMS.check"), "w") as c:
    for f in files:
        c.write(f"{f['sha256']}  {f['name']}\n")

trig = ", ".join(f"{r['label']} {r['entries']}" for r in order)
line = (f"| hodoscope_extnim_2026 | {socket.gethostname()} | {os.path.realpath(RAW)}/ | {utc(first)[:10]} – {utc(last)[:10]} (UTC) "
        f"| {trig} (total {sum(r['entries'] for r in runs)}) | {total} B ({total / 1e9:.3f} GB), {len(files)} archivos | MANIFEST.sha256 |")
with open(os.path.join(OUT, "INDEX_line.md"), "w") as o:
    o.write("<!-- borrador XPLUS-1 D: línea de data/raw/INDEX.md (monorepo); no hay un INDEX.md de ejemplo en docs/monorepo_handoff -->\n")
    o.write("| Campaign | Host | Path | Dates | Triggers per run | Size | Manifest |\n|---|---|---|---|---|---|---|\n")
    o.write(line + "\n")
print(line)

git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
json.dump({"stage": "XPLUS-1 D", "commit": git, "generated_utc": utc(time.time()), "host": socket.gethostname(), "raw_dir": os.path.realpath(RAW),
           "hash_seconds": round(time.monotonic() - t0, 1), "total_bytes": total, "n_files": len(files), "files": files,
           "compare_existing_SHA256SUMS": cmp, "only_in_existing_SHA256SUMS": extra,
           "triggers_source": os.path.realpath(INV), "dates_utc": [utc(first), utc(last)], "index_line": line},
          open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
with open(os.path.join(OUT, "SHA256SUMS"), "w") as s:
    for n in ("MANIFEST.sha256", "SHA256SUMS.check", "INDEX_line.md", "manifest.json"):
        s.write(f"{sha256(os.path.join(OUT, n))}  {n}\n")
