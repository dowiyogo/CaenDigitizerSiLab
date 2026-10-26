#!/usr/bin/python3 -I
# proc_audit_md.py — XPLUS-4 etapa 7: escribe audit.md desde audit.json y audit_numbers.csv (proc_audit.py). Las discrepancias de la
# sección 1 son las que proc_audit.py marcó (redondeo incorrecto, sin fuente) más las que revisé a mano con los valores calculados que
# se dan en cada fila (f_2019, G de octubre). Uso: proc_audit_md.py <dir de la auditoría>
import csv, json, sys
from collections import Counter

A = sys.argv[1]
J = json.load(open(A + "/audit.json"))
rows = list(csv.DictReader(open(A + "/audit_numbers.csv")))
st = Counter(r["status"].split(" (")[0] for r in rows)
L = []
w = L.append
w("# Auditoría numérica de main_v5.tex (PROC-5, XPLUS-4 etapa 7)\n")
w(f"Fuente: `git show origin/main:main_v5.tex` del clon del proceeding tras `git fetch origin`. `origin/main` = `{J['origin_main']}` "
  "(merge \"Merge PROC-4 figures\"). Copia leída: `source/main_v5.tex` (SHA256 `ebb5d93b…`). No se cambió la rama del clon (sigue en `main`) "
  "ni ningún `.tex`.\n")
w("## 1. Discrepancias y cifras sin fuente (lo primero que hay que mirar)\n")
w("| línea(s) | texto | fuente | valor de la fuente | estado | propuesta |")
w("|---|---|---|---|---|---|")
w("| 45, 395 | 99.28\\% | 73627/74157, de la línea 144 | 99.2853 % | **redondeo incorrecto** (truncado); las líneas 144 y 146 dicen 99.29\\% | 99.29\\% en las cuatro apariciones |")
w("| 45, 259, 409 | ±0.10 (término del estimador de ε_XY a 5 p.e.) | `results[id=epsXY_may8_good_T5].uncertainty.syst.estimator` | 0.0949 | **redondeo incorrecto** "
  "(0.0949 → 0.09; probablemente redondeado dos veces, 0.095 → 0.10) | ±0.09, o el valor que René elija según `deltas_manuscrito.md` |")
w("| 234 | f_2019 = 19.24\\% | `tables.stage2.mc_fractions_1e8.f.h19` | 0.192347 | **redondeo incorrecto** (→ 19.23 %) | 19.23\\% |")
w("| 266 | \"median G of 3.8 mV per p.e.\" | mediana de `tables.calibration[].G_c1_mV` (29 canales) | 3.857 mV | **redondeo incorrecto** (→ 3.9) | 3.9 mV |")
w("| 266, 413 | \"0.87 of the May value\", \"13% lower\" | mediana de G_c1 / mediana de la G de mayo reajustada (4.56 mV, línea 206) = 0.846; "
  "`results[id=G_may_over_c1_median]` = 1.1508, y 1/1.1508 = 0.869, usa la G **del archivo** de mayo y no la de la línea 206 | 0.846 o 0.869 | "
  "**definición inconsistente**: el 0.87 no es el cociente de las dos medianas que da el texto | o \"0.85 of the May value\" y \"15% lower\" (medianas del texto), "
  "o decir que es la mediana por canal de G_c1/G_mayo(archivo) |")
w("| 137, 144, 148, 150 | NA64: 126307, 74157, 73627, 1277597, 1192842, 117700 | — | — | **sin fuente en el repositorio** (análisis NA64 de R. Carrera) | René o R.C. confirman |")
w("| 179 | n entre −1.79 y −1.84 | — | — | **sin fuente en results.json** (ajuste de la ley de potencias, herramientas geométricas) | regenerar o citar la herramienta |\n")
w("No hubo otras discrepancias en las cifras mapeadas explícitamente: f_2018, q, ε_XY, los dos estimadores, el χ², las cuentas de la tabla 2 y de "
  "las 8073 del taller, la ocupación, f_2, la respuesta a lo largo de las barras, hod2019, la longitud sensible y G.\n")
w("## 2. Tabla completa de números\n")
w(f"`audit_numbers.csv`: {len(rows)} números de las líneas 42-428, con línea, valor tal como está, estado, fuente, valor de la fuente y contexto. Estados: "
  + ", ".join(f"{k} {v}" for k, v in st.most_common()) + ".\n")
w("Cómo leer los estados:")
w("- **igual / redondeo correcto / redondeo incorrecto:** cifras mapeadas a mano a una clave de `results.json` v3, o calculadas del propio texto (columna `source`). Son verificaciones.")
w("- **candidato:** existe en `results.json`, en `tables/data/*.tex` (origin/main) o en `report_es.md` un valor que redondea exactamente al del texto "
  "(×1, ×100 o ×0.01). **No es una verificación:** con 2 cifras significativas hay coincidencias por azar. Por ejemplo, χ² = 1.7 coincidía con una razón de "
  "ocupación hasta que se mapeó a mano. Sirve para orientarse.")
w("- **sin fuente:** la búsqueda o el mapa no encontraron fuente (ver 1).")
w("- **no auditado:** enteros pequeños, constantes de diseño (dimensiones de las barras y de las placas) y números estructurales.")
w("- Limitación conocida: el mapa a mano se aplica a todas las apariciones de un valor en su línea. En la línea 268, el segundo \"0.45\" (\"0.45 to 0.48 of the cases\", "
  "factor geométrico) quedó asociado a ε(Y|X); su fuente correcta es el MC de barras.\n")
w("## 3. Figuras\n")
w("| línea | archivo | en origin/main |")
w("|---|---|---|")
for f in J["figures"]:
    w(f"| {f['line']} | `{f['file']}` | {f['in_origin_main']} |")
w(f"\nLos {len(J['figures'])} `\\includegraphics` apuntan a archivos que existen en `origin/main`: "
  f"{'sí' if all(f['in_origin_main'] != 'NO' for f in J['figures']) else 'NO'}.")
p4 = J["proc4_vs_main"]
a, b = p4["figures/data/eff_along_bar.png"], p4["figures/data/amp_spectra_hod2018_pe.png"]
w(f"\n**PROC-4:** el texto de la tarea suponía que `main` todavía tenía las versiones de cuatro paneles. Ya no es así. `eff_along_bar.png` y "
  f"`amp_spectra_hod2018_pe.png` de `origin/main` son el mismo blob que en `origin/figs/2026-10-25-proc4` (`{a['origin_main_blob'][:10]}`, `{b['origin_main_blob'][:10]}`; "
  f"iguales: {a['same']}, {b['same']}). René ya mezcló las de dos paneles (`0158c92`).\n")
w("PNG de `figures/` en origin/main: píxeles y dpi de la cabecera pHYs. No afirmo ningún requisito de JAIS; René compara con las guías.\n")
w("| archivo | ancho px | alto px | dpi (pHYs) | bytes |")
w("|---|---|---|---|---|")
for k, v in J["pngs"].items():
    if v:
        w(f"| `{k}` | {v['width_px']} | {v['height_px']} | {v['dpi_pHYs'] if v['dpi_pHYs'] else 'sin pHYs'} | {v['bytes']} |")
c = J["counts"]
w("\n## 4. Conteos\n")
w(f"- Palabras del cuerpo, sin abstract, figuras, tablas, ecuaciones ni bibliografía (calculado: LaTeX quitado, cada `$…$` cuenta como 1, las marcas "
  f"`\\verif` no cuentan): **{c['body_words_no_abstract_figs_tables_refs']}**. JAIS pide 4000-8000.")
w(f"- Palabras del abstract: **{c['abstract_words']}** (límite 150).")
w(f"- Palabras de las leyendas (aparte): {c['captions_words']}.")
w(f"- Palabras clave: **{c['n_keywords']}** ({', '.join(c['keywords'])}); se piden 7-10.")
w(f"- Bibliografía: {c['n_bib_entries']} entradas, {c['n_cited_keys']} claves citadas. Entradas no citadas: {c['bib_not_cited'] or 'ninguna'}. "
  f"Citas sin entrada: {c['cited_without_entry'] or 'ninguna'}.")
w("- El método de conteo es aproximado (calculado, no el de la revista): otro criterio puede diferir en unas decenas de palabras. El cuerpo queda "
  "a unas 480 palabras del máximo y el abstract a 2 del límite.\n")
w(f"## 5. Marcas \\verif ({len(J['verif'])} en el texto)\n")
w("La tarea habla de 14. Supongo que ese número incluye la línea 16, que es el comentario del preámbulo que define el comando y da su uso: no es una marca.\n")
w("| línea | texto |")
w("|---|---|")
for v in J["verif"]:
    w(f"| {v['line']} | {v['text'].replace('|', '/')} |")
w("\n## 6. Consistencia interna\n")
w("Apariciones de cada cantidad repetida, tal como están escritas:\n")
for k, occ in J["consistency"].items():
    vals = sorted({o[1] for o in occ})
    w(f"- **{k}:** {len(occ)} apariciones; formas: {', '.join('`' + x + '`' for x in vals) if vals else '—'}; líneas {', '.join(str(o[0]) for o in occ)}")
w("\nInconsistencias encontradas (calculado):")
w("- **Eficiencia del haz:** 99.28\\% (líneas 45 y 395) contra 99.29\\% (144 y 146). Ver 1.")
w("- **Ganancia de octubre:** \"0.87 of the May value\" (266) y \"13% lower\" (413) son consistentes entre sí, pero no con la mediana de mayo de la línea 206 (4.56 mV → 0.846). Ver 1.")
w("- **ε_XY a 5 p.e.:** la misma cadena `0.80 ± 0.07 ± 0.10 ± 0.04` en las líneas 45, 259 y 409, así que un cambio del término del estimador hay que hacerlo en las tres.")
w("- **Ángulo cenital medio:** 33.7° (línea 238) y 34° (línea 332); redondeo consistente.")
w("- **f_both:** 1.655 ± 0.008 % (236) y 1.66 % (238, 397, 413); redondeo consistente.")
w("- **Separación D:** 98.4 mm en todas las apariciones.")
w("- **Coincidencias de hod2019:** \"3 coincidences against 19 expected\" (45), \"3 X–Y coincidences, against 19\" (411) y la tabla (3); consistente.")
w("- **Razón de factores geométricos:** \"60 times larger\" (397, 413): 66.4/1.10 = 60.4; consistente.")
w("- **Frecuencia de dos estaciones:** \"9.2 mHz, about 800 per day\" (311): 2380/(71.83 h) = 9.20 mHz y 9.20e-3 × 86400 = 795 por día; consistente.")
open(A + "/audit.md", "w").write("\n".join(L) + "\n")
print(len(L), "líneas")
