#!/usr/bin/python3 -I
# stage6.py — PROC-1 etapa 6 (sin datos de detector): geometría del telescopio y MC del banco.
#  - Factor geométrico de dos cuadrados coaxiales paralelos de lado a = 33 mm a distancia D (Sullivan 1971, ec. G del
#    manuscrito): G = ∫∫ (a−|u|)(a−|v|) D²/(u²+v²+D²)² du dv (isótropo) y G₂ = ∫∫ (a−|u|)(a−|v|) D⁴/(u²+v²+D²)³ du dv
#    (intensidad ∝ cos²θ), por cuadratura de Gauss-Legendre (dos órdenes como control) y con el MC (bench_mc) como control
#    independiente. Ley de potencia entre 30 y 300 mm con tres espaciados, exponente local, efecto de 1 mm a 150 mm, G₂ por
#    par de píxeles (píxel a/16) y tasas en cielo abierto con I_v = 70 m⁻² s⁻¹ sr⁻¹ (el valor del manuscrito).
#  - MC del banco (bench_mc.py, 1e8 rayos, semilla 42): f_2018, f_2019, f_both y ángulo cenital medio para cos³θ y cos²θ,
#    placas OLD y CAD (A5); P(trigger | cruza las dos áreas activas); G₂ del trigger y del telescopio a 98.4 mm; área efectiva.
#  - Sensibilidad: work/geom_sens_<OLD|CAD>_1e8.json (geom_sensitivity.py extendido: desplazamientos, rotaciones, D ± 1 mm,
#    placas z ± 2 mm; f_2018, f_2019, f_both).
# Salida: work/stage6.json, tables/geom_sensitivity.tex.
# Uso: stage6.py <dir del paquete> [N]
import json, math, subprocess, sys, time

import numpy as np

REPO = "/home/xplus/digitizer/CaenDigitizerSiLab-jai"
sys.path.insert(0, REPO + "/daq/tools/characterization")
import bench_mc as bm  # noqa: E402

OUT = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000_000
W = OUT + "/work"
git = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
S = {"stage": 6, "commit": git, "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
A = 33.0
IV = 70.0  # m⁻² s⁻¹ sr⁻¹ (manuscrito, PDG 2022 sec. 30.3.1)


# ---------------------------------------------------------------- factor geométrico
def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


def G(D, m=1, a=A, n=200):
    """m = 1: isótropo; m = 2: pesado con cos²θ. mm² sr."""
    u, wu = gl(n, 0.0, a)
    U, V = np.meshgrid(u, u, indexing="ij")
    f = (a - U) * (a - V) * D ** (2 * m) / (U ** 2 + V ** 2 + D ** 2) ** (m + 1)
    return 4.0 * float((wu[:, None] * wu[None, :] * f).sum())


def G_pair(D, dx, dy, m=2, p=A / 16, n=40):
    """G₂ de un par de píxeles de lado p desplazados (dx, dy) píxeles: ∫∫ tri(u − dx p) tri(v − dy p) cos^(2m)… du dv."""
    u, wu = gl(n, dx * p - p, dx * p + p)
    v, wv = gl(n, dy * p - p, dy * p + p)
    U, V = np.meshgrid(u, v, indexing="ij")
    tri = np.clip(p - np.abs(U - dx * p), 0, None) * np.clip(p - np.abs(V - dy * p), 0, None)
    return float((wu[:, None] * wv[None, :] * tri * D ** (2 * m) / (U ** 2 + V ** 2 + D ** 2) ** (m + 1)).sum())


cm2 = lambda x: x / 100.0
GF = {}
for D in (98.4, 150.0):
    GF[str(D)] = {"G_cm2sr": cm2(G(D)), "G_cm2sr_n400": cm2(G(D, n=400)), "G2_cm2sr": cm2(G(D, 2)), "G2_cm2sr_n400": cm2(G(D, 2, n=400)),
                  "far_field_a4_D2_cm2sr": cm2(A ** 4 / D ** 2), "far_field_over_G": (A ** 4 / D ** 2) / G(D),
                  "angular_pitch_mrad": 2.0 / D * 1e3, "bar_width_over_D_mrad": 3.0 / D * 1e3, "max_angle_deg": math.degrees(math.atan(A / D))}
    GF[str(D)]["rel_precision_n200_n400"] = abs(GF[str(D)]["G_cm2sr"] / GF[str(D)]["G_cm2sr_n400"] - 1)
fits = {}
for name, Ds in (("10 puntos lineales", np.linspace(30, 300, 10)), ("28 puntos cada 10 mm", np.arange(30, 301, 10.0)),
                 ("20 puntos logarítmicos", np.geomspace(30, 300, 20)), ("271 puntos cada 1 mm", np.arange(30, 301, 1.0))):
    g = np.array([G(D, n=120) for D in Ds])
    fits[name] = float(np.polyfit(np.log(Ds), np.log(g), 1)[0])
loc = {}
for D in (30.0, 150.0, 300.0):
    h = 1e-3 * D
    loc[str(D)] = (math.log(G(D + h)) - math.log(G(D - h))) / (math.log(D + h) - math.log(D - h))
d1 = {"G_151_over_150_minus_1": G(151.0) / G(150.0) - 1, "G_149_over_150_minus_1": G(149.0) / G(150.0) - 1, "far_field_2dD_over_D": 2 / 150}
# píxeles y tasas en cielo abierto (D = 150 mm)
gp0 = G_pair(150.0, 0, 0)
gpc = G_pair(150.0, 15, 15)
tot_pairs = sum((16 - abs(i)) * (16 - abs(j)) * G_pair(150.0, i, j) for i in range(-15, 16) for j in range(-15, 16))
rate = lambda g_mm2: IV * g_mm2 * 1e-6  # s⁻¹
S["geometric_factor"] = {"per_D": GF, "power_law_n_30_300": fits, "local_exponent": loc, "one_mm_at_150": d1,
                         "open_sky_D150": {"G2_cm2sr": cm2(G(150.0, 2)), "rate_s": rate(G(150.0, 2)), "per_day": rate(G(150.0, 2)) * 86400,
                                           "directions": 31 * 31, "per_direction_per_day": rate(G(150.0, 2)) * 86400 / 961,
                                           "on_axis_256_pairs_per_day": rate(256 * gp0) * 86400, "corner_1_pair_per_day": rate(gpc) * 86400,
                                           "sum_of_pixel_pairs_over_G2": tot_pairs / G(150.0, 2), "days_for_100_on_axis": 100 / (rate(256 * gp0) * 86400)},
                         "I_v_m2_s_sr": IV, "method": "Gauss-Legendre 200 × 200 (control 400 × 400) en un cuadrante"}
print("G:", {k: (round(v["G_cm2sr"], 4), round(v["G2_cm2sr"], 4)) for k, v in GF.items()}, "n:", fits, "local:", loc, "1 mm:", d1)
print("cielo abierto:", S["geometric_factor"]["open_sky_D150"])

# ---------------------------------------------------------------- MC del banco
AGEN = 400.0 * 400.0  # mm², plano de generación de bench_mc
bench = {}
for geom in ("OLD", "CAD"):
    for model in ("cos3", "cos2"):
        g = bm.GEOMS[geom]
        rng = np.random.default_rng(42)
        c = {"rays": 0, "trig": 0, "h18": 0, "h19": 0, "both": 0, "both_no_trig_req": 0, "th_trig": 0.0, "th_both": 0.0}
        for start in range(0, N, 2_000_000):
            m = min(2_000_000, N - start)
            r = bm.generate(rng, m, model)
            t = bm.trigger(r, g)
            a18, a19 = bm.in_square(r, 0.0), bm.in_square(r, g["H19"])
            th = np.degrees(np.arccos(r.cos_theta()))
            c["rays"] += m
            c["trig"] += int(t.sum())
            c["h18"] += int((a18 & t).sum())
            c["h19"] += int((a19 & t).sum())
            c["both"] += int((a18 & a19 & t).sum())
            c["both_no_trig_req"] += int((a18 & a19).sum())
            c["th_trig"] += float(th[t].sum())
            c["th_both"] += float(th[a18 & a19 & t].sum())
        f = lambda k: (c[k] / c["trig"], math.sqrt(c[k] / c["trig"] * (1 - c[k] / c["trig"]) / c["trig"]))
        d = {"counts": c, "f_2018": f("h18"), "f_2019": f("h19"), "f_both": f("both"), "mean_zenith_trig_deg": c["th_trig"] / c["trig"],
             "mean_zenith_both_deg": c["th_both"] / c["both"], "P_trig_given_both": c["both"] / c["both_no_trig_req"]}
        if model == "cos3":  # G₂ = A_gen · (π/2) · P (direcciones ∝ cos³θ: intensidad cos²θ por el factor de proyección)
            d["G2_trigger_cm2sr"] = cm2(AGEN * math.pi / 2 * c["trig"] / c["rays"])
            d["G2_telescope_cm2sr"] = cm2(AGEN * math.pi / 2 * c["both_no_trig_req"] / c["rays"])
            d["A_eff_trigger_cm2"] = d["G2_trigger_cm2sr"] / (math.pi / 2)
            d["G2_ratio_telescope_over_trigger"] = d["G2_telescope_cm2sr"] / d["G2_trigger_cm2sr"]
        bench[f"{geom}_{model}"] = d
        print(geom, model, "f18 %.5f f19 %.5f fboth %.5f ± %.5f θ %.2f / %.2f P(trig|both) %.6f" % (d["f_2018"][0], d["f_2019"][0], d["f_both"][0], d["f_both"][1],
              d["mean_zenith_trig_deg"], d["mean_zenith_both_deg"], d["P_trig_given_both"]), {k: d[k] for k in d if k.startswith(("G2", "A_eff"))})
S["bench"] = bench
b = bench["OLD_cos3"]
S["derived"] = {"rate_1_per_cm2_min_Hz": b["A_eff_trigger_cm2"] / 60.0, "IvG2_trigger_Hz": IV * b["G2_trigger_cm2sr"] * 1e-4,
                "trigger_over_telescope_G2": b["G2_trigger_cm2sr"] / b["G2_telescope_cm2sr"], "G2_telescope_analytic_98p4": GF["98.4"]["G2_cm2sr"],
                "two_station_rate_mHz_lambda_hat_v8": b["f_both"][0] * 0.5561374216263751 * 1e3,
                "two_station_per_day_lambda_hat_v8": b["f_both"][0] * 0.5561374216263751 * 86400}
print("derivados:", S["derived"])

# ---------------------------------------------------------------- sensibilidad
SENS = {g: json.load(open(f"{W}/geom_sens_{g}_1e8.json")) for g in ("OLD", "CAD")}
S["sensitivity"] = SENS
L = ["% PROC-1 etapa 6: sensibilidad geométrica del MC del banco (bench_mc.py + geom_sensitivity.py, 1e8 rayos, semilla 42, cos^3).",
     "% Cambio relativo respecto de la nominal (mismos rayos); en las variaciones de la estación se mueven las dos estaciones.", f"% Commit {git}.",
     "\\begin{tabular}{lrrr}", "\\hline", "Variación & $\\Delta f_{2018}/f_{2018}$ (\\%) & $\\Delta f_{2019}/f_{2019}$ (\\%) & $\\Delta f_\\mathrm{both}/f_\\mathrm{both}$ (\\%) \\\\",
     "\\hline"]
nom = SENS["OLD"]["cos3"]["nominal"]
L.append(f"Nominal (placas de geom\\_acceptance.py) & {100 * nom['f']:.3f} & {100 * nom['f_2019']:.2f} & {100 * nom['f_both']:.3f} \\\\")
P = lambda v, e: f"{100 * v:+.2f} $\\pm$ {100 * e:.2f}"
for k, v in SENS["OLD"]["cos3"].items():
    if k == "nominal":
        continue
    L.append(f"{k} & {P(v['rel_change'], v['rel_change_err'])} & {P(v['rel_change_2019'], v['rel_change_2019_err'])} & "
             f"{P(v['rel_change_both'], v['rel_change_both_err'])} \\\\")
cad = SENS["CAD"]["cos3"]["nominal"]
L.append(f"Placas del CAD (A5; z +1.5 mm) & {100 * (cad['f'] / nom['f'] - 1):+.2f} & {100 * (cad['f_2019'] / nom['f_2019'] - 1):+.2f} & "
         f"{100 * (cad['f_both'] / nom['f_both'] - 1):+.2f} \\\\")
L += ["\\hline", "\\end{tabular}"]
open(OUT + "/tables/geom_sensitivity.tex", "w").write("\n".join(L) + "\n")
json.dump(S, open(W + "/stage6.json", "w"), indent=1, default=float)
print("hecho")
