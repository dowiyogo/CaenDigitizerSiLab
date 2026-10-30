#!/usr/bin/python3 -I
# mip_gain_vs_mv.py — seguimiento 2026-10-31, punto 4: ¿la diferencia de 30 % entre las medianas MIP de mayo (19.5 p.e.) y octubre (13.5)
# está en la calibración (G) o en la amplitud en mV?
# Mismas entradas y definición que mip_spectra.py (archivos de impactos de PROC, amplitud de pico en la ventana de 144 ns, espectro
# sustraído = signal − fuera de tiempo, bins negativos a 0). Por canal (hod2018 sin ch3):
#   G de mayo (reajuste PROC, hd64_may_v8_peaks.csv) y de octubre (c1), grupo HV de cada canal en octubre (r4, medido) y en mayo
#   (hv_group_previous de r4, medido con apagones de mayo, solo donde existe) y el V0Set de su grupo (octubre: nominal de hvmon;
#   mayo: captura de GECO2020 del handoff, línea 21, igual al nominal);
#   mediana del espectro sustraído en mV con umbral 5 p.e. de su propia G y con umbral fijo de 20 mV, y en p.e.
# Suma de canales en mV y en p.e.; mayo con la G c1 de octubre en los canales con |G_mayo/G_c1 − 1| ≤ 5 %.
# Salida: <out>/mip_gain_vs_mv.json y .csv. Uso: mip_gain_vs_mv.py <out dir>
import csv, json, os, sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "daq/tools/characterization"))
import may_common as mc  # noqa: E402
import xy_estimator as xy  # noqa: E402

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
W = os.path.join(REPO, "outputs/proceeding_2026_10/work")
CH = [c for c in range(32) if c != 3]
V0SET = {1: 65.75, 2: 65.75, 3: 65.90, 4: 65.85, 5: 65.75, 6: 65.70, 7: 66.15, 8: 65.90}  # grupo g = CH0(g−1); mayo = igual (GECO)
gm = {int(r["ch"]): float(r["G_fit"]) for r in csv.DictReader(open(W + "/hd64_may_v8_peaks.csv")) if r["status"] in ("ok", "changed")}
c1 = json.load(open(os.path.join(REPO, "daq/config/calibration/pe_calibration_2026_10_c1.json")))
go = {e["ch"]: e["G_mV"] for e in c1["channels"] if e["ch"] < 32 and isinstance(e.get("G_mV"), (int, float))}
med_o = float(np.median([go[c] for c in CH if c in go]))
mp = {c["dig_ch"]: c for c in json.load(open(os.path.join(REPO, "daq/config/mapping/hodo_mapping_2026_10_r4.json")))["channels"]}
E = {"may": xy.load_hits(W + "/hd64_may_v8_hits.bin", 64), "oct": xy.load_hits(W + "/hd64_oct8_hits.bin", 64)}
E["may"] = E["may"][mc.good_time(E["may"]["t"], "v8")]
EDGES = np.arange(0, 1101, 1.0)  # mV


def sub_hist(e, c, thr_mv):
    s, o = e["aS"][:, c], e["aC"][:, c]
    hs = np.histogram(s[s > thr_mv], bins=EDGES)[0]
    ho = np.histogram(o[o > thr_mv], bins=EDGES)[0]
    return np.clip(hs - ho, 0, None).astype(float), int((s > thr_mv).sum() - (o > thr_mv).sum())


def median_of(h):
    cs = np.cumsum(h)
    return float(EDGES[np.searchsorted(cs, 0.5 * cs[-1])] + 0.5) if cs[-1] > 0 else None


rows, sums = [], {}
for c in CH:
    Gm, Go = gm.get(c), go.get(c, med_o)
    r = {"ch": c, "G_may_mV": Gm, "G_oct_mV": Go, "G_oct_is_median": c not in go, "G_ratio_may_over_oct": Gm / Go if Gm else None,
         "hv_group_oct": mp[c].get("hv_group"), "hv_group_may": (mp[c].get("hv_group_previous") or {}).get("hv_group") if mp[c].get("hv_group_previous") else None}
    r["V0Set_oct_V"] = V0SET.get(r["hv_group_oct"])
    r["V0Set_may_V"] = V0SET.get(r["hv_group_may"]) if r["hv_group_may"] else None
    for k, e, G in (("may", E["may"], Gm), ("oct", E["oct"], Go)):
        h5, n5 = sub_hist(e, c, 5 * G)
        h20, n20 = sub_hist(e, c, 20.0)
        r[f"{k}_median_mV_thr5pe"] = median_of(h5)
        r[f"{k}_median_pe_thr5pe"] = median_of(h5) / G if median_of(h5) else None
        r[f"{k}_median_mV_thr20mV"] = median_of(h20)
        r[f"{k}_net_hits_thr5pe"], r[f"{k}_net_hits_thr20mV"] = n5, n20
        for tag, hh in (("thr5pe", h5), ("thr20mV", h20)):
            sums.setdefault(f"{k}_{tag}_mV", np.zeros(len(EDGES) - 1))
            sums[f"{k}_{tag}_mV"] += hh
    rows.append(r)
# suma en p.e.: cada canal con su G (rebinado a 1 p.e.) — mediana de la suma de espectros en p.e.
def sum_pe(k, gains, chans, thr_pe=5):
    edges_pe = np.arange(0, 301, 1.0)
    tot = np.zeros(300)
    for c in chans:
        e = E[k]
        s, o = e["aS"][:, c] / gains[c], e["aC"][:, c] / gains[c]
        tot += np.clip(np.histogram(s[s > thr_pe], bins=edges_pe)[0] - np.histogram(o[o > thr_pe], bins=edges_pe)[0], 0, None)
    cs = np.cumsum(tot)
    return float(edges_pe[np.searchsorted(cs, 0.5 * cs[-1])] + 0.5), float(cs[-1])


Gm_all = {c: gm[c] for c in CH}
Go_all = {c: go.get(c, med_o) for c in CH}
same = [c for c in CH if c in go and abs(gm[c] / go[c] - 1) <= 0.05]
res = {"channels": rows, "V0Set_by_group_V": V0SET, "V0Set_may_source": "captura de GECO2020 (handoff línea 21; René no recuerda en cuáles de v7-v10)",
       "sum_median_mV": {k: median_of(v) for k, v in sums.items()},
       "sum_median_pe_own_gain": {"may": sum_pe("may", Gm_all, CH), "oct": sum_pe("oct", Go_all, CH)},
       "channels_gain_within_5pct": same,
       "sum_median_pe_same_channels": {"may_with_may_gain": sum_pe("may", Gm_all, same), "may_with_c1_gain": sum_pe("may", Go_all, same),
                                       "oct_with_c1_gain": sum_pe("oct", Go_all, same)},
       "median_G_ratio_may_over_oct": float(np.median([r["G_ratio_may_over_oct"] for r in rows if r["G_ratio_may_over_oct"]]))}
json.dump(res, open(OUT + "/mip_gain_vs_mv.json", "w"), indent=1)
with open(OUT + "/mip_gain_vs_mv.csv", "w") as f:
    keys = list(rows[0].keys())
    f.write(",".join(keys) + "\n")
    for r in rows:
        f.write(",".join("" if r[k] is None else (f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])) for k in keys) + "\n")
for r in rows:
    print(f"ch{r['ch']:02d} G mayo {r['G_may_mV']:.3f} oct {r['G_oct_mV']:.3f}{'*' if r['G_oct_is_median'] else ' '} razón {r['G_ratio_may_over_oct']:.3f} "
          f"grupo oct {r['hv_group_oct']} ({r['V0Set_oct_V']} V) mayo {r['hv_group_may']} | mediana mV (5 p.e.) mayo {r['may_median_mV_thr5pe']} oct {r['oct_median_mV_thr5pe']} "
          f"| (20 mV) mayo {r['may_median_mV_thr20mV']} oct {r['oct_median_mV_thr20mV']} | p.e. mayo {r['may_median_pe_thr5pe']:.1f} oct {r['oct_median_pe_thr5pe']:.1f}")
print(json.dumps({k: v for k, v in res.items() if k != "channels"}, indent=1))
