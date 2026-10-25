#!/usr/bin/python3 -I
# secret_scan.py — XPLUS-2 A.1: busca secretos en toda la historia (git log -p --all), sin reescribir nada.
# Lee el diff en streaming (nunca la historia entera en memoria); por cada línea añadida ('+') aplica la lista PATTERNS y anota
# commit, archivo y patrón. El valor que coincide se enmascara en la salida (primeros 3 caracteres + longitud): el informe no debe
# copiar un secreto si lo hubiera. Además lista los nombres de archivo sospechosos que alguna vez estuvieron en la historia.
# No lee ~/.config/daqv2/secrets.env (regla 10 de CLAUDE.md: solo modo y nombres de variables).
# Uso: secret_scan.py <repo> <salida.json>
import json, re, subprocess, sys

REPO, OUT = sys.argv[1], sys.argv[2]
PATTERNS = {
    "password_assign": r"(?i)\b(pass(word|wd)?|pwd)\b[\"']?\s*[:=]\s*[\"']?([^\s\"',;)]{3,})",
    "daq_secret_env": r"\bDAQ_(HV_PASSWORD|HV_USER|SMTP_[A-Z_]+|NOTIFY_[A-Z_]+)\s*=\s*[\"']?([^\s\"']+)",
    "token_apikey": r"(?i)\b(token|api[_-]?key|secret|auth)\b[\"']?\s*[:=]\s*[\"']?([A-Za-z0-9_\-+/=]{12,})",
    "private_key": r"-----BEGIN [A-Z ]*PRIVATE KEY-----",
    "aws_key": r"\bAKIA[0-9A-Z]{16}\b",
    "github_token": r"\bgh[pousr]_[A-Za-z0-9]{36,}\b",
    "google_key": r"\bAIza[0-9A-Za-z_\-]{35}\b",
    "slack_token": r"\bxox[baprs]-[A-Za-z0-9-]{10,}",
    "url_userinfo": r"[a-z]+://[^\s:/@\"']+:([^\s@/\"']+)@",
    "curl_user": r"curl\b.*\s(-u|--user)\s+[^\s:]+:(\S+)",
    "caenhv_init_literal": r"InitSystem\s*\(.*\"[^\"]+\"\s*,\s*\"[^\"]+\"\s*,\s*\"[^\"]+\"",
    "admin_admin": r"(?i)\badmin\b\W{1,6}\badmin\b",
    "ip_with_cred": r"192\.168\.\d+\.\d+.*(?i:pass|pwd|user)",
}
RX = {k: re.compile(v) for k, v in PATTERNS.items()}
SUSPECT_NAMES = re.compile(r"(^|/)(secrets\.env|[^/]*\.env|hv_credentials\.json|id_(rsa|ed25519|ecdsa)[^/]*|[^/]*\.(pem|key|p12|pfx)|\.netrc|credentials[^/]*)$")


def mask(s):
    return f"{s[:3]}…({len(s)} car.)"


hits = []
p = subprocess.Popen(["git", "-C", REPO, "--no-optional-locks", "log", "--all", "-p", "--no-color", "--no-ext-diff", "--format=@@COMMIT %H"],
                     stdout=subprocess.PIPE, text=True, errors="replace")
commit = path = None
nlines = ncommits = 0
for line in p.stdout:
    if line.startswith("@@COMMIT "):
        commit = line.split()[1]
        ncommits += 1
        continue
    if line.startswith("+++ "):
        path = line[6:].rstrip("\n") if line.startswith("+++ b/") else line[4:].rstrip("\n")
        continue
    if not line.startswith("+") or line.startswith("+++"):
        continue
    nlines += 1
    s = line[1:].rstrip("\n")
    if len(s) > 2000:  # líneas enormes (datos, svg): solo patrones de clave
        s = s[:2000]
    for k, rx in RX.items():
        m = rx.search(s)
        if m:
            val = m.group(m.lastindex) if m.lastindex else m.group(0)
            ctx = s.strip()
            if val and len(val) >= 3:
                ctx = ctx.replace(val, mask(val))
            hits.append({"pattern": k, "commit": commit, "path": path, "value_masked": mask(val or ""), "context_masked": ctx[:200]})
p.wait()

names = subprocess.run(["git", "-C", REPO, "log", "--all", "--format=", "--name-only", "--diff-filter=A"], capture_output=True, text=True).stdout.split("\n")
suspect = sorted({n for n in names if n and SUSPECT_NAMES.search(n)})
json.dump({"patterns": PATTERNS, "commits_scanned": ncommits, "added_lines_scanned": nlines, "hits": hits, "suspect_file_names": suspect},
          open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"commits {ncommits}  líneas añadidas {nlines}  coincidencias {len(hits)}  nombres sospechosos {suspect}")
from collections import Counter
print(Counter(h["pattern"] for h in hits))
