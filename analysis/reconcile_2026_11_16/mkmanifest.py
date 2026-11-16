#!/usr/bin/python3 -I
# mkmanifest.py — escribe <out>/manifest.json (comandos, SHA256 de entradas derivadas y crudas, de los archivos del directorio) y <out>/SHA256SUMS.
# Uso: mkmanifest.py <out> <rama> <commit base> <código,…> <entradas derivadas,…> <archivo de SHA256 de crudos o ""> <commands.txt>
import json, hashlib, os, sys, datetime, platform
O, branch, base = sys.argv[1], sys.argv[2], sys.argv[3]
code = sys.argv[4].split(",")
der = sys.argv[5].split(",") if sys.argv[5] else []
raw = {}
if len(sys.argv) > 6 and sys.argv[6]:
    for l in open(sys.argv[6]):
        h, p = l.split(); raw[p + " (READ)"] = h
cmds = open(sys.argv[7]).read().splitlines() if len(sys.argv) > 7 else []
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
files = {}
for d, _, fs in os.walk(O):
    for f in fs:
        if f not in ("manifest.json", "SHA256SUMS"):
            p = os.path.join(d, f); files[os.path.relpath(p, O)] = sha(p)
m = {"created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "branch": branch, "base_commit": base, "code": code,
     "commands": cmds, "machine": platform.node(), "raw_inputs": raw, "derived_inputs": {p: sha(p) for p in der}, "files": files}
json.dump(m, open(O + "/manifest.json", "w"), indent=1)
with open(O + "/SHA256SUMS", "w") as f:
    for k in sorted(files): f.write(f"{files[k]}  {k}\n")
