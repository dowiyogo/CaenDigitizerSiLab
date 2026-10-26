#!/usr/bin/python3 -I
# proc_audit.py — XPLUS-4 etapa 7 (PROC-5, solo lectura): auditoría numérica de main_v5.tex (origin/main del proceeding).
#  1. Extrae cada número del abstract, el cuerpo, las leyendas y las tablas (línea, contexto, valor como está escrito).
#  2. Fuentes: outputs/proceeding_2026_10/results.json (v3, aplanado: results, tables, series, maps), tables/data/*.tex de origin/main
#     y outputs/proceeding_2026_10/report_es.md. Estados:
#       - mapa explícito (KEY): igual / redondeo correcto / redondeo incorrecto / discrepancia, con la clave de la fuente;
#       - búsqueda automática: "candidato" si algún valor de las fuentes (en la escala ×1, ×100 o ×0.01) redondea exactamente al
#         número del texto con sus decimales; "sin fuente encontrada" si no; enteros < 20 y constantes de diseño: "no auditado".
#  3. \includegraphics contra el árbol de origin/main; tamaño y dpi (pHYs) de cada PNG de figures/.
#  4. Palabras (cuerpo sin abstract, figuras, tablas ni bibliografía; abstract), palabras clave, citas y bibliografía.
#  5. Marcas \verif. 6. Cantidades repetidas con valores distintos.
# Uso: proc_audit.py <dir de salida> <clon del proceeding>
import csv, json, math, os, re, struct, subprocess, sys, zlib
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT, PROC = sys.argv[1], sys.argv[2]
SRC = os.path.join(OUT, "source")
TEX = open(os.path.join(SRC, "main_v5.tex")).read().split("\n")
R = json.load(open(os.path.join(REPO, "outputs/proceeding_2026_10/results.json")))


def git(*a, binary=False):
    r = subprocess.run(["git", "-C", PROC, "--no-optional-locks", *a], capture_output=True)
    return r.stdout if binary else r.stdout.decode()


# ---------------- fuentes aplanadas
SRCVALS = []  # (valor, ubicación)


def walk(o, p):
    if isinstance(o, dict):
        for k, v in o.items():
            walk(v, f"{p}/{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            walk(v, f"{p}[{i}]")
    elif isinstance(o, (int, float)) and not isinstance(o, bool) and math.isfinite(o):
        SRCVALS.append((float(o), "results.json" + p))


walk(R, "")
NUM = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
for f in sorted(os.listdir(os.path.join(SRC, "tables_data"))):
    for i, l in enumerate(open(os.path.join(SRC, "tables_data", f)), 1):
        for m in NUM.finditer(l.replace("\\,", "")):
            SRCVALS.append((float(m.group()), f"tables/data/{f}:{i}"))
for i, l in enumerate(open(os.path.join(REPO, "outputs/proceeding_2026_10/report_es.md")), 1):
    for m in NUM.finditer(l.replace("\u2009", "").replace("\u202f", "")):
        SRCVALS.append((float(m.group()), f"report_es.md:{i}"))
SRCVALS.sort()
SV = [v for v, _ in SRCVALS]


def lookup(path):
    o = R
    for part in re.findall(r"[^/\[\]]+|\[\d+\]", path):
        if part.startswith("["):
            o = o[int(part[1:-1])]
        elif isinstance(o, list):
            o = next(x for x in o if x.get("id") == part)
        else:
            o = o[part]
    return o


# ---------------- mapa explícito de las cifras principales: (línea, valor en el texto, fuente, función → valor fuente)
def r_(id_, *ks):
    def f():
        o = lookup(f"results/{id_}")
        for k in ks:
            o = o[k]
        return o
    return f, f"results[id={id_}]" + "".join(f".{k}" for k in ks)


def t_(path):
    return (lambda: lookup(path)), f"results.json/{path}"


def calc(expr, desc):
    return (lambda: expr()), desc


qsum = lambda a, b: math.hypot(a, b)
KEY = [
    (144, "99.29", calc(lambda: 100 * 73627 / 74157, "73627/74157 (texto, línea 144)")),
    (45, "99.28", calc(lambda: 100 * 73627 / 74157, "73627/74157 (texto, línea 144)")),
    (395, "99.28", calc(lambda: 100 * 73627 / 74157, "73627/74157 (texto, línea 144)")),
    (234, "3.522", calc(lambda: 100 * lookup("results/f_2018_CAD_cos3")["value"], "results[id=f_2018_CAD_cos3].value ×100")),
    (234, "0.011", calc(lambda: 100 * lookup("results/f_2018_CAD_cos3")["uncertainty"]["stat"], "results[id=f_2018_CAD_cos3].uncertainty.stat ×100")),
    (257, "0.0275", r_("q_may8_good_T5", "value")),
    (257, "0.0023", r_("q_may8_good_T5", "uncertainty", "stat")),
    (259, "0.80", r_("epsXY_may8_good_T5", "value")),
    (259, "0.07", r_("epsXY_may8_good_T5", "uncertainty", "stat")),
    (259, "0.10", r_("epsXY_may8_good_T5", "uncertainty", "syst", "estimator")),
    (259, "0.04", calc(lambda: qsum(lookup("results/epsXY_may8_good_T5")["uncertainty"]["syst"]["acceptance"], lookup("results/epsXY_may8_good_T5")["uncertainty"]["syst"]["g"]),
                       "results[id=epsXY_may8_good_T5]: aceptancia ⊕ g")),
    (45, "0.80", r_("epsXY_may8_good_T5", "value")),
    (409, "0.80", r_("epsXY_may8_good_T5", "value")),
    (262, "0.0208", r_("q_may8_good_T10", "value")),
    (262, "0.0007", r_("q_may8_good_T10", "uncertainty", "stat")),
    (262, "0.60", r_("epsXY_may8_good_T10", "value")),
    (262, "0.03", calc(lambda: qsum(lookup("results/epsXY_may8_good_T10")["uncertainty"]["syst"]["acceptance"], lookup("results/epsXY_may8_good_T10")["uncertainty"]["syst"]["g"]),
                       "results[id=epsXY_may8_good_T10]: aceptancia ⊕ g")),
    (266, "0.69", r_("epsXY_oct8_T5", "value")),
    (266, "0.05", r_("epsXY_oct8_T5", "uncertainty", "stat")),
    (266, "0.06", r_("epsXY_oct8_T5", "uncertainty", "syst", "estimator")),
    (266, "0.35", r_("epsXY_oct8_T10", "value")),
    (264, "0.84", r_("published_convention_v8", "value")),
    (264, "0.12", r_("published_convention_v8", "uncertainty", "stat")),
    (208, "113\\,152", r_("v8_good_time_events", "value")),
    (206, "-0.15", r_("may_p0_refit_median", "value")),
    (282, "3111", r_("map_qN_may8_good", "value")),
    (284, "2.6", r_("occ_maxmin_may8_good_X_T5", "value")),
    (284, "0.8", r_("occ_maxmin_may8_good_X_T5", "uncertainty", "stat")),
    (284, "1.6", r_("occ_maxmin_may8_good_Y_T5", "value")),
    (284, "0.3", r_("occ_maxmin_may8_good_Y_T5", "uncertainty", "stat")),
    (284, "1.9", r_("occ_maxmin_may8_good_X_T10", "value")),
    (284, "2.9", r_("occ_maxmin_may8_good_Y_T10", "value")),
    (334, "0.390", r_("f2_may8_good_X_T5", "value")),
    (334, "0.013", r_("f2_may8_good_X_T5", "uncertainty", "stat")),
    (334, "0.366", r_("f2_may8_good_Y_T5", "value")),
    (334, "0.304", r_("f2_may8_good_X_T10", "value")),
    (334, "0.265", r_("f2_may8_good_Y_T10", "value")),
    (350, "0.77", r_("alongbar_may8_good_XA_T5", "value")),
    (350, "1.27", r_("alongbar_may8_good_XB_T5", "value")),
    (350, "0.69", r_("alongbar_may8_good_XA_T10", "value")),
    (350, "1.31", r_("alongbar_may8_good_XB_T10", "value")),
    (350, "0.85", r_("alongbar_may8_good_YA_T5", "value")),
    (350, "1.03", r_("alongbar_may8_good_YA_T10", "value")),
    (350, "0.56", r_("alongbar_may8_good_YB_T5", "value")),
    (350, "0.57", r_("alongbar_may8_good_YB_T10", "value")),
    (368, "0.0144", r_("h19_rX_35mV", "value")),
    (368, "0.088", r_("h19_rX_3mV", "value")),
    (368, "19", r_("h19_rXrYN_35mV", "value")),
    (368, "950", r_("h19_rXrYN_3mV", "value")),
    (257, "0.0308", t_("tables/stage2/exclude_noisy/5/with/f_ctl")),
    (257, "0.0033", t_("tables/stage2/exclude_noisy/5/with/f_ctl_err")),
    (262, "0.0214", t_("tables/stage2/exclude_noisy/10/with/f_ctl")),
    (262, "0.0008", t_("tables/stage2/exclude_noisy/10/with/f_ctl_err")),
    (262, "1.7", t_("tables/stage2/stability/5/runs_chi2[1]")),
    (262, "2.6", t_("tables/stage2/stability/5/blocks_chi2[1]")),
    (262, "0.0250", t_("tables/stage2/exclude_noisy/5/without_2_9_15/f_mix")),
    (268, "0.37", t_("tables/stage2/exclude_noisy/5/with/e_X|Y_mix")),
    (268, "0.03", t_("tables/stage2/exclude_noisy/5/with/e_X|Y_mix_err")),
    (268, "0.45", t_("tables/stage2/exclude_noisy/5/with/e_Y|X_mix")),
    (268, "0.04", t_("tables/stage2/exclude_noisy/5/with/e_Y|X_mix_err")),
    (264, "2707", t_("tables/stage2/manuscript_8073/run B completo (como el manuscrito)/true_XY_mix")),
    (264, "152", t_("tables/stage2/manuscript_8073/run B completo (como el manuscrito)/err_true_mix")),
    (264, "2515", t_("tables/stage2/manuscript_8073/run B completo (como el manuscrito)/true_XY_ctl")),
    (264, "212", t_("tables/stage2/manuscript_8073/run B completo (como el manuscrito)/err_true_ctl")),
    (321, "2707", t_("tables/stage2/manuscript_8073/run B completo (como el manuscrito)/true_XY_mix")),
    (321, "152", t_("tables/stage2/manuscript_8073/run B completo (como el manuscrito)/err_true_mix")),
    (372, "0.322", r_("block_hod2019_c_X_block_in_overlap", "value")),
    (322, "3985", calc(lambda: lookup("results/f_2018_CAD_cos3")["value"] * 113152, "f_2018 (results[id=f_2018_CAD_cos3]) × 113152")),
    (192, "65536", calc(lambda: 256 ** 2, "256² (pares de píxeles)")),
    (372, "22", r_("sensL_35mV_LX", "value")),
    (372, "30", r_("sensL_3mV_LX", "value")),
    (177, "0.51", r_("G_D150.0", "value")),
    (177, "1.14", r_("G_D98.4", "value")),
]

# ---------------- extracción de números del texto
def clean(line):
    s = line
    s = re.sub(r"\\(label|ref|cite|includegraphics|eqref|graphicspath|bibliography|bibliographystyle|documentclass|journal|vol|auno|footnotemark)(\[[^\]]*\])?\{[^}]*\}", " ", s)
    s = re.sub(r"\\(cite)\[[^\]]*\]\{[^}]*\}", " ", s)
    s = re.sub(r"\[width=[^\]]*\]|\[height=[^\]]*\]", " ", s)
    s = re.sub(r"(\d)\\,(\d)", r"\1\2", s)  # 126\,307 -> 126307
    s = re.sub(r"(\d(?:\.\d+)?)\\times10\^\{?(-?\d+)\}?", r"\1e\2", s)  # 3.5\times10^{-3} -> 3.5e-3
    s = re.sub(r"\\pm|\\mp|\\times|\\approx|\\le|\\ge", " ", s)
    s = s.replace("\\,", " ")
    return s


RE_NUM = re.compile(r"(?<![A-Za-z_\\{\d.])(-?\d+(?:\.\d+)?)(?:e(-?\d+))?(?![\d])")


def decimals(t):
    return len(t.split(".")[1]) if "." in t else 0


def rounds_to(v, t):
    d = decimals(t)
    try:
        x = float(t)
    except ValueError:
        return False
    return abs(round(v, d) - x) < 0.5 * 10 ** (-d) * 1e-6 + 1e-12 or round(v, d) == x


def status_of(v, t):
    x = float(t.replace("\\,", "").replace(",", ""))
    if v == x:
        return "igual"
    d = decimals(t)
    if round(v, d) == x or abs(round(v, d) - x) < 1e-9:
        return "redondeo correcto"
    if abs(v - x) <= 10 ** (-d) * 1.0001:
        return "redondeo incorrecto"
    digits = t.lstrip("-").replace(".", "").lstrip("0")
    sig = len(digits.rstrip("0")) if "." not in t else len(digits)
    if sig and v != 0 and float(f"{v:.{sig}g}") == x:
        return f"redondeo correcto ({sig} cifras significativas)"
    return "discrepancia"


# rango de líneas: del abstract (42) a los agradecimientos (428); se saltan el preámbulo, comentarios y los argumentos estructurales
rows = []
for i in range(41, 429):
    raw = TEX[i]
    if raw.strip().startswith("%"):
        continue
    s = clean(raw)
    for m in RE_NUM.finditer(s):
        t = m.group(1)
        if m.group(2):
            t = f"{t}e{m.group(2)}"
        ctx = s[max(0, m.start() - 45): m.end() + 30].replace("\n", " ")
        rows.append({"line": i + 1, "value_text": t, "context": ctx.strip()})
NOTES = {(137, "126307"): "sin fuente en el repositorio (análisis NA64 de R. Carrera)", (144, "74157"): "sin fuente en el repositorio (análisis NA64)",
         (144, "73627"): "sin fuente en el repositorio (análisis NA64)", (146, "99.29"): "calculado: 73627/74157 = 99.2853 (línea 144); redondeo correcto",
         (148, "1277597"): "sin fuente en el repositorio (análisis NA64)", (148, "1192842"): "sin fuente en el repositorio (análisis NA64)",
         (150, "117700"): "sin fuente en el repositorio (análisis NA64)", (179, "-1.79"): "sin fuente en results.json (ajuste de la ley de potencias, herramientas geométricas)",
         (179, "-1.84"): "sin fuente en results.json (ajuste de la ley de potencias)", (204, "-42.5"): "constante de diseño del banco (CAD; bench_mc.py), no auditado",
         (49, "10.31526"): "no es una cifra (DOI)", (91, "2572"): "no es una cifra (modelo S12572)"}
# estado por fila
keymap = defaultdict(list)
for ln, t, (f, desc) in KEY:
    keymap[(ln, t.replace("\\,", ""))].append((f, desc))
for r in rows:
    t = r["value_text"]
    k = (r["line"], t)
    if k in keymap:
        f, desc = keymap[k][0]
        try:
            v = float(f())
            r.update(source=desc, source_value=repr(v), status=status_of(v, t))
        except Exception as e:  # clave no encontrada
            r.update(source=desc, source_value=f"error: {e}", status="sin fuente (clave no encontrada)")
        continue
    if k in NOTES:
        r.update(source=NOTES[k], source_value="", status=("redondeo correcto" if "redondeo correcto" in NOTES[k] else ("no auditado" if "no " in NOTES[k][:3] or "no auditado" in NOTES[k] else "sin fuente")))
        continue
    try:
        x = float(t)
    except ValueError:
        r.update(source="", source_value="", status="no auditado")
        continue
    if ("e" not in t and "." not in t and abs(x) < 20) or abs(x) in (0,):
        r.update(source="", source_value="", status="no auditado (entero pequeño o estructural)")
        continue
    cands = []
    for scale in (1.0, 0.01, 100.0):
        target = x * scale
        d = decimals(t) + (2 if scale == 0.01 else (-2 if scale == 100.0 else 0))
        tol = 0.5 * 10 ** (-max(d, 0)) if d >= 0 else 0.5 * 10 ** (-d)
        import bisect
        lo, hi = bisect.bisect_left(SV, target - tol), bisect.bisect_right(SV, target + tol)
        for v, loc in SRCVALS[lo:hi]:
            if d >= 0 and round(v, d) == round(target, d):
                cands.append((loc, v, scale))
        if cands:
            break
    if cands:
        loc, v, sc = sorted(cands, key=lambda c: (not c[0].startswith("results.json/results"), not c[0].startswith("tables"), len(c[0])))[0]
        r.update(source=f"{loc}" + ("" if sc == 1 else f" (×{1/sc:g})"), source_value=repr(v), status=f"candidato ({len(cands)} coincidencias)")
    else:
        r.update(source="", source_value="", status="sin fuente encontrada")
with open(os.path.join(OUT, "audit_numbers.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["line", "value_text", "status", "source", "source_value", "context"])
    w.writeheader()
    w.writerows(rows)

# ---------------- figuras
tree = set(git("ls-tree", "-r", "--name-only", "origin/main").split())
figs = []
for i, l in enumerate(TEX, 1):
    for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", l):
        name = m.group(1)
        found = [p for p in ([f"figures/{name}", f"figures/data/{name}", f"figures/sim/{name}"]) if p in tree]
        figs.append({"line": i, "file": name, "in_origin_main": found[0] if found else "NO"})


def png_info(path):
    b = git("show", f"origin/main:{path}", binary=True)
    if b[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    w, h = struct.unpack(">II", b[16:24])
    dpi = None
    pos = 8
    while pos < len(b):
        ln, typ = struct.unpack(">I4s", b[pos:pos + 8])
        if typ == b"pHYs":
            px, py, unit = struct.unpack(">IIB", b[pos + 8:pos + 17])
            dpi = (round(px * 0.0254, 1), round(py * 0.0254, 1)) if unit == 1 else ("sin unidad", px, py)
        if typ == b"IDAT":
            break
        pos += 12 + ln
    return {"width_px": w, "height_px": h, "dpi_pHYs": dpi, "bytes": len(b)}


pngs = {p: png_info(p) for p in sorted(tree) if p.startswith("figures/") and p.endswith(".png")}
# la rama figs/2026-10-25-proc4 contra origin/main
proc4 = {}
for p in ("figures/data/eff_along_bar.png", "figures/data/amp_spectra_hod2018_pe.png"):
    a = git("rev-parse", f"origin/main:{p}").strip()
    b = git("rev-parse", f"origin/figs/2026-10-25-proc4:{p}").strip()
    proc4[p] = {"origin_main_blob": a, "figs_proc4_blob": b, "same": a == b}

# ---------------- palabras, citas
text = "\n".join(TEX)
abs_ = re.search(r"\\begin\{abstract\}(.*?)\\noindent\\textit\{Keywords", text, re.S).group(1)
kw = re.search(r"\\textit\{Keywords:\}(.*?)\\par", text, re.S).group(1)
body = text[text.index("\\maketitle"): text.index("\\section*{Conflicts of Interest}")]
body_nofloat = re.sub(r"\\begin\{(figure\*?|table\*?|equation)\}.*?\\end\{\1\}", " ", body, flags=re.S)


def words(s):
    s = re.sub(r"(?m)(?<!\\)%.*$", " ", s)
    s = re.sub(r"\\verif\{[^}]*\}", " ", s)
    s = re.sub(r"\\(label|ref|cite|section|subsection|footnotetext|renewcommand)\*?(\[[^\]]*\])?\{[^}]*\}", " ", s)
    s = re.sub(r"\$[^$]*\$", " X ", s)
    s = re.sub(r"\\[A-Za-z]+\*?", " ", s)
    s = re.sub(r"[{}~\\]", " ", s)
    return len([w for w in re.split(r"\s+", s) if re.search(r"[A-Za-z0-9]", w)])


caps = re.findall(r"\\caption\{(.*?)\}\n\\label", body, re.S)
cites = re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", text)
cited = sorted({c.strip() for g in cites for c in g.split(",")})
bib = open(os.path.join(SRC, "references.bib")).read()
entries = sorted(re.findall(r"@\w+\{([^,]+),", bib))
counts = {"body_words_no_abstract_figs_tables_refs": words(body_nofloat), "abstract_words": words(abs_), "captions_words": sum(words(c) for c in caps),
          "keywords": [k.strip() for k in kw.split(",")], "n_keywords": len(kw.split(",")), "n_bib_entries": len(entries), "n_cited_keys": len(cited),
          "bib_not_cited": [e for e in entries if e not in cited], "cited_without_entry": [c for c in cited if c not in entries]}
verifs = []
for i, l in enumerate(TEX, 1):
    if l.lstrip().startswith("%"):
        continue  # la línea 16 es el comentario del preámbulo que define \verif, no una marca
    for m in re.finditer(r"\\verif\{((?:[^{}]|\{[^{}]*\})*)\}", l):
        if "#1" not in m.group(1):
            verifs.append({"line": i, "text": m.group(1)})

# ---------------- consistencia interna: mismos nombres con valores distintos
PAT = {"eps_XY (5 p.e.)": r"0\.80\s*\\pm\s*0\.07\s*\\pm\s*0\.10\s*\\pm\s*0\.04", "eps_XY (10 p.e.)": r"0\.60\s*\\pm\s*0\.02\s*\\pm\s*0\.02\s*\\pm\s*0\.03",
       "eficiencia del haz": r"99\.2\d\\?%", "f_2018": r"3\.5\d\d?\s*\\pm|3\.522|3\.52\\%", "f_both": r"1\.65\d|1\.66\\%", "g": r"g=0\.98|0\.980|98\\%",
       "D banco": r"98\.4", "triggers run B": r"143\D{0,3}814", "tiempo bueno": r"113\D{0,3}152", "coincidencias hod2019": r"\b3 (X--Y )?coincidences|3 X--Y", "G a 150 mm": r"0\.51\\,cm", "ratio 60": r"60 times"}
cons = {}
for k, p in PAT.items():
    occ = []
    for i, l in enumerate(TEX, 1):
        for m in re.finditer(p, l):
            occ.append((i, m.group(0)))
    cons[k] = occ
json.dump({"origin_main": open(os.path.join(SRC, "origin_main_commit.txt")).read().strip(), "figures": figs, "pngs": pngs, "proc4_vs_main": proc4,
           "counts": counts, "verif": verifs, "consistency": cons}, open(os.path.join(OUT, "audit.json"), "w"), indent=1, ensure_ascii=False)
from collections import Counter
print("números", len(rows), Counter(r["status"].split(" (")[0] for r in rows))
print("figuras", len(figs), "faltan", [f for f in figs if f["in_origin_main"] == "NO"])
print("proc4", proc4)
print(json.dumps(counts, ensure_ascii=False))
print("verif", len(verifs))
