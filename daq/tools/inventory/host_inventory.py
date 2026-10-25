#!/usr/bin/python3 -I
# host_inventory.py — XPLUS-1 parte A (solo lectura): inventario de este host para decidir qué entra al monorepo muography-thesis.
#  1. Repositorios git bajo /home/xplus (find -name .git, profundidad ≤ 4): remoto, rama, fecha del último commit, commits sin
#     empujar por rama (no alcanzables desde ninguna referencia remota local; no se hace fetch), cambios sin commit, tamaño.
#     Todo con `git --no-optional-locks` (no refresca ni escribe el índice de otros repositorios).
#  2. Scripts de análisis fuera de los repositorios (.C, .cxx, .py, .ipynb, .sh) bajo /home/xplus y /data: ruta, tamaño, mtime y
#     primera línea de comentario, agrupados por directorio. Se excluyen árboles de terceros y de herramientas (lista EXCLUDE).
#  3. Datos crudos del digitizador (.root, .dat, archivos de WaveDump/CAEN) bajo /data y /home/xplus: por directorio, número de archivos,
#     tamaño, fecha del primero y del último, y si hay un índice (SHA256SUMS, *.sha256, MANIFEST*, index.json).
#  4. Discos: df -h, lsblk, mount; 5. unidades systemd del DAQ, timers, cron, screen/tmux, procesos del DAQ.
# Cada comando se guarda con su salida textual en <out>/cmd_*.txt; el resumen en <out>/inventory.json.
# Uso: host_inventory.py <dir de salida>
import json, os, re, subprocess, sys, time
from collections import defaultdict

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
HOME = "/home/xplus"
EXCLUDE = [HOME + p for p in ("/root", "/root6", "/Externals", "/src", "/snap", "/.cache", "/.vscode-server", "/.local/lib", "/.local/share",
                              "/.gemini", "/.codex", "/.claude", "/.config", "/.mozilla", "/.npm", "/.nvm")]
CMDLOG = []


def run(name, cmd, timeout=600):
    """ejecuta y guarda la salida textual (los comandos son de solo lectura)"""
    p = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, timeout=timeout)
    txt = f"$ {cmd if isinstance(cmd, str) else ' '.join(cmd)}\n{p.stdout}{p.stderr}"
    open(f"{OUT}/cmd_{name}.txt", "w").write(txt)
    CMDLOG.append(name)
    return p.stdout


def git(repo, *args):
    p = subprocess.run(["git", "--no-optional-locks", "-C", repo, *args], capture_output=True, text=True, timeout=300)
    return p.stdout.strip()


inv = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": os.uname().nodename, "exclude_for_loose_scripts": EXCLUDE}

# ---------------------------------------------------------------- 1. repositorios
gitdirs = [l for l in run("find_git", f"find {HOME} -maxdepth 4 -name .git -not -path '*/.cache/*'").split() if l]
repos = []
log1 = []
for g in sorted(gitdirs):
    r = os.path.dirname(g)
    branches = git(r, "for-each-ref", "--format=%(refname:short)", "refs/heads").splitlines()
    unp = {b: int(git(r, "rev-list", "--count", b, "--not", "--remotes") or 0) for b in branches}
    st = git(r, "status", "--porcelain")
    du = subprocess.run(["du", "-sb", r], capture_output=True, text=True).stdout.split()
    d = {"path": r, "remotes": git(r, "remote", "-v").splitlines(), "branch": git(r, "rev-parse", "--abbrev-ref", "HEAD"),
         "last_commit": git(r, "log", "-1", "--format=%ci %h %s"), "branches": branches, "unpushed_per_branch": unp,
         "uncommitted": len([l for l in st.splitlines() if l.strip()]), "uncommitted_untracked": len([l for l in st.splitlines() if l.startswith("??")]),
         "size_bytes": int(du[0]) if du else None}
    repos.append(d)
    log1.append(f"== {r}\n$ git remote -v\n" + "\n".join(d["remotes"]) + f"\n$ git rev-parse --abbrev-ref HEAD\n{d['branch']}\n$ git log -1 --format='%ci %h %s'\n{d['last_commit']}\n"
                + "".join(f"$ git rev-list --count {b} --not --remotes\n{n}\n" for b, n in unp.items())
                + f"$ git status --porcelain | wc -l\n{d['uncommitted']}\n$ du -sb\n{d['size_bytes']}\n")
open(f"{OUT}/cmd_git_repos.txt", "w").write("\n".join(log1))
inv["repos"] = repos
repo_paths = [d["path"] for d in repos]

# ---------------------------------------------------------------- 2. scripts sueltos
ext = ("*.C", "*.cxx", "*.py", "*.ipynb", "*.sh")
expr = " -o ".join(f"-name '{e}'" for e in ext)
prune = " -o ".join(f"-path '{p}'" for p in EXCLUDE + repo_paths)
out2 = run("find_scripts", f"find {HOME} /data \\( {prune} \\) -prune -o -type f \\( {expr} \\) -printf '%p\\t%s\\t%TY-%Tm-%Td %TH:%TM\\n'")
loose = defaultdict(list)
for l in out2.splitlines():
    p, s, m = l.split("\t")
    first = ""
    try:
        if p.endswith(".ipynb"):
            j = json.load(open(p))
            cells = j.get("cells", [])
            src = "".join(cells[0].get("source", [])) if cells else ""
            first = src.strip().splitlines()[0][:160] if src.strip() else ""
        else:
            with open(p, errors="replace") as f:
                for line in f:
                    t = line.strip()
                    if t.startswith(("#!",)):
                        continue
                    marks = ("//", "/*", "*") if p.endswith((".C", ".cxx")) else ("#", '"""', "'''")
                    if t.startswith(marks):
                        first = t[:160]
                        break
                    if t:
                        break
    except Exception as ex:  # noqa: BLE001
        first = f"(no legible: {type(ex).__name__})"
    loose[os.path.dirname(p)].append({"file": os.path.basename(p), "size": int(s), "mtime": m, "first_comment": first})
inv["loose_scripts"] = {k: v for k, v in sorted(loose.items())}

# ---------------------------------------------------------------- 3. datos crudos
out3 = run("find_data", f"find /data {HOME} \\( {' -o '.join(f'-path {chr(39)}{p}{chr(39)}' for p in EXCLUDE)} \\) -prune -o -type f "
                        f"\\( -name '*.root' -o -name '*.dat' -o -iname 'wave*.txt' -o -iname '*caen*' \\) -printf '%h\\t%f\\t%s\\t%TY-%Tm-%Td %TH:%TM\\n'")
data = defaultdict(lambda: {"files": 0, "bytes": 0, "first": None, "last": None, "ext": defaultdict(int)})
for l in out3.splitlines():
    h, f, s, m = l.split("\t")
    d = data[h]
    d["files"] += 1
    d["bytes"] += int(s)
    d["first"] = m if d["first"] is None or m < d["first"] else d["first"]
    d["last"] = m if d["last"] is None or m > d["last"] else d["last"]
    d["ext"][os.path.splitext(f)[1] or f] += 1
for h, d in data.items():
    try:
        names = os.listdir(h)
    except OSError:
        names = []
    d["index_files"] = sorted(n for n in names if re.search(r"(SHA256SUMS|\.sha256$|MANIFEST|index\.json$|INDEX)", n))
    d["ext"] = dict(d["ext"])
    d["in_repo"] = any(h.startswith(rp + "/") or h == rp for rp in repo_paths)
inv["data_dirs"] = {k: v for k, v in sorted(data.items())}

# ---------------------------------------------------------------- 4-5. discos, unidades, cron, sesiones, procesos
run("du_repo", f"du -sh --apparent-size {HOME}/digitizer/CaenDigitizerSiLab-jai/.git {HOME}/digitizer/CaenDigitizerSiLab-jai/outputs {HOME}/digitizer/CaenDigitizerSiLab-jai/daq/build*")
run("df", "df -h")
run("lsblk", "lsblk -o NAME,SIZE,TYPE,FSTYPE,MOUNTPOINT,MODEL")
run("mount", "mount | grep -E '^/dev'")
run("du_data", "du -sh /data/daq/* 2>/dev/null")
run("systemd_user_units", "systemctl --user list-unit-files --no-pager; echo; systemctl --user list-units --all --no-pager | grep -i -E 'hv|daq|bias|scan|report|recover' ")
run("systemd_user_timers", "systemctl --user list-timers --all --no-pager")
run("systemd_system", "systemctl list-unit-files --no-pager | grep -i -E 'daq|caen|hv|digit|a4818|conet' ; echo; systemctl list-timers --all --no-pager")
run("user_unit_files", f"ls -la {HOME}/.config/systemd/user {HOME}/.config/systemd/user/default.target.wants 2>&1")
run("cron", "crontab -l 2>&1; echo '--- /etc/crontab'; cat /etc/crontab 2>&1; echo '--- /etc/cron.d'; ls -la /etc/cron.d 2>&1")
run("screen_tmux", "screen -ls 2>&1; tmux ls 2>&1")
run("processes", "ps -eo pid,user,lstart,etime,cmd | grep -i -E 'caen|wavedump|digitizer|daq|acq|hvmon|hvctl|biasscan|root\\.exe|python3' | grep -v -E 'grep|host_inventory'")
inv["commands"] = CMDLOG
json.dump(inv, open(f"{OUT}/inventory.json", "w"), indent=1, default=str)
print("repos", len(repos), "dirs con scripts sueltos", len(loose), "scripts", sum(len(v) for v in loose.values()), "dirs con datos", len(data))
