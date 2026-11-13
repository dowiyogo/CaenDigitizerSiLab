#!/usr/bin/python3 -I
# check_be9e8de.py — SOCHIFI 2026-11-13: comprueba que los números de P1-P4 son los de la revisión be9e8de del proceeding (rama remota
# text/revision-2026-11-02, leída en un clon de trabajo). Compara (1) el SHA256 de la figura del proceeding con el de la salida de la que se toman
# los datos, (2) los números del texto (main_v5.tex) y de tables/data/doublehit_confirmed.tex con los que dibujan las figuras del póster.
# Uso: check_be9e8de.py <clon de be9e8de> <dir del póster>
import hashlib, json, os, re, sys

CL, OUT = sys.argv[1], sys.argv[2]
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
O = os.path.join(REPO, "outputs")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
tex = open(os.path.join(CL, "main_v5.tex")).read()
nm = json.load(open(os.path.join(OUT, "poster_numbers.json")))["numbers"]
res = {}


def same(a, b):
    return "IGUAL" if a == b else "DISTINTO"


print("HEAD del clon:", os.popen(f"git -C {CL} log -1 --format='%h %s'").read().strip())
# P1
a, b = os.path.join(CL, "figures/data/trigger_interarrival_stationary.pdf"), os.path.join(O, "proceeding_revision/2026-11-02_revision/trigger_interarrival_stationary.pdf")
print("P1 figura be9e8de vs fuente:", same(sha(a), sha(b)), sha(a)[:16])
t1 = ["142\\,544 triggers in 66.58\\,h", "$0.5947\\pm0.0017$\\,Hz", "$p=0.27$"]
p1 = {"triggers": nm["stationary_triggers"]["value"], "horas": nm["stationary_duration"]["value"], "lambda": (nm["trigger_rate"]["value"], nm["trigger_rate"]["error"]),
      "p_KS": round(nm["trigger_ks_p"]["value"], 2)}
print("P1 texto be9e8de contiene:", {s: (s in tex) for s in t1}, "| póster:", p1)
res["P1"] = all(s in tex for s in t1) and p1["triggers"] == 142544 and p1["horas"] == 66.58 and p1["lambda"] == (0.5947, 0.0017) and p1["p_KS"] == 0.27
# P2
a, b = os.path.join(CL, "figures/data/mip_spectra_paper_2p.pdf"), os.path.join(O, "proceeding_revision/2026-11-02_revision/fig_2p/mip_spectra_paper_2p.pdf")
print("P2 figura be9e8de vs fuente:", same(sha(a), sha(b)), sha(a)[:16])
s2 = "has its maximum at 13.2\\,p.e.\\ in the layer shown"
print("P2 texto be9e8de contiene:", {s2: s2 in tex}, "| póster: moda", nm["mode_XB_29_54mm"]["value"], "-> %.1f" % nm["mode_XB_29_54mm"]["value"])
res["P2"] = s2 in tex and round(nm["mode_XB_29_54mm"]["value"], 1) == 13.2
# P3
a, b = os.path.join(CL, "figures/data/eff_vs_threshold.png"), os.path.join(O, "proceeding_revision/2026-11-03_band_legend/eff_vs_threshold.png")
print("P3 figura be9e8de vs fuente:", same(sha(a), sha(b)), sha(a)[:16])
ok3 = True
for T, txt in ((5, "0.81 \\pm 0.09 \\pm 0.07 \\pm 0.04"), (10, "0.63\\pm0.02\\pm0.02\\pm0.03")):
    e = nm[f"eps_XY_{T}pe"]
    r = "%.2f %.2f %.2f %.2f" % (e["value"], e["error"]["stat"], e["error"]["estimator"], e["error"]["acceptance"])
    print(f"P3 T={T}: texto contiene '{txt}': {txt in tex} | póster redondeado {r}")
    ok3 &= (txt in tex) and [float(x) for x in re.findall(r"[\d.]+", txt)] == [float(x) for x in r.split()]
print("P3 banda: texto '4.7\\%' (aceptancia):", "acceptance (4.7\\%" in tex, "| póster", "%.4f" % nm["acceptance_rel_uncertainty"]["value"])
res["P3"] = ok3
# P4
tb = open(os.path.join(CL, "tables/data/doublehit_confirmed.tex")).read()
dh = json.load(open(os.path.join(REPO, "outputs/proceeding_hvdiag5/2026-10-31_full/work/doublehit_confirmed.json")))["data"]["may8_good"]
ok4 = True
for pl in "XY":
    for r in dh[pl]["rows"]:
        m = re.search(rf"mayo v8, tiempo bueno & {pl} & {r['T']} & ([\d.]+) \$\\pm\$ ([\d.]+)", tb)
        if not m:
            continue
        mine = ("%.3f" % r["all"], "%.3f" % r["all_err"])
        ok = mine == (m.group(1), m.group(2))
        ok4 &= ok
        print(f"P4 {pl} T={r['T']}: tabla be9e8de {m.group(1)} ± {m.group(2)} | póster {mine[0]} ± {mine[1]} -> {'IGUAL' if ok else 'DISTINTO'}")
for s in ("$f_2=0.402\\pm0.017\\pm0.010$ for X and $0.366\\pm0.015$ for Y at 5\\,p.e.", "the fractions are $0.320\\pm0.014$ and $0.270\\pm0.013$",
          "0.449 to 0.454", "0.513 to 0.514"):
    print(f"P4 texto contiene '{s}': {s in tex}")
    ok4 &= s in tex
mcx, mcy = nm["two_bar_fraction_X_MC"]["value"], nm["two_bar_fraction_Y_MC"]["value"]
print(f"P4 MC del póster (promedio de los dos órdenes de bloque): X {mcx:.4f} en [0.449, 0.454]: {0.449 <= mcx <= 0.454}; Y {mcy:.4f} en [0.513, 0.514]: {0.513 <= mcy <= 0.5145}")
ok4 &= 0.449 <= mcx <= 0.454 and 0.513 <= mcy <= 0.5145
a, b = os.path.join(CL, "figures/data/double_hit_fraction.pdf"), os.path.join(O, "proceeding_hvdiag5/2026-10-31_full/figs/double_hit_fraction.pdf")
print("P4 figura double_hit_fraction de be9e8de (no está en main_v5.tex) vs la del paquete:", same(sha(a), sha(b)), "| includegraphics en main_v5.tex:", "double_hit_fraction" in tex)
res["P4"] = ok4
print("RESULTADO:", {k: ("coincide" if v else "DIFIERE") for k, v in res.items()})
