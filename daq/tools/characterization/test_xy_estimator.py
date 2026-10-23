#!/usr/bin/python3 -I
# test_xy_estimator.py — DAQ2-16: pruebas del estimador X·Y compartido (xy_estimator.py).
#  - sintética: muones con q y rX, rY conocidos más ruido independiente aX, aY: los dos estimadores recuperan q;
#  - regresión con el run 8 (2026-10-22/23): definición única, ch3 fuera, 5 p.e., mezclados = 0.0239 ± 0.0017
#    (outputs/characterization_2026_10_23/report_es.md, "Comparación con mayo"). Se salta si faltan los impactos.
# Uso: python3 -I daq/tools/characterization/test_xy_estimator.py   (también desde ctest)
import json, os, sys, unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import xy_estimator as xy  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
HITS8 = os.path.join(REPO, "outputs/characterization_2026_10_23/work/hd_run008_hits.bin")
C1 = os.path.expanduser("~/.local/opt/daqv2/share/calibration/pe_calibration_2026_10_c1.json")


class Synthetic(unittest.TestCase):
    def test_recovers_q_with_independent_noise(self):
        rng = np.random.default_rng(7)
        N, q, rx_only, ry_only, aX, aY = 400_000, 0.025, 0.035, 0.045, 0.20, 0.15
        u = rng.uniform(size=N)
        mx = u < q + rx_only
        my = (u < q) | ((u >= q + rx_only) & (u < q + rx_only + ry_only))
        nX, nY = rng.uniform(size=N) < aX, rng.uniform(size=N) < aY
        cX, cY = rng.uniform(size=N) < aX, rng.uniform(size=N) < aY
        r = xy.xy_fraction(mx | nX, my | nY, cX, cY)
        for tag in ("mix", "ctl"):
            self.assertAlmostEqual(r[f"f_{tag}"], q, delta=4 * r[f"f_{tag}_err"], msg=tag)
        self.assertAlmostEqual(r["rX"], q + rx_only, delta=0.003)
        self.assertAlmostEqual(r["e_X|Y_mix"], q / (q + ry_only), delta=0.05)

    def test_no_signal_gives_zero(self):
        rng = np.random.default_rng(3)
        N = 200_000
        a = [rng.uniform(size=N) < p for p in (0.2, 0.2, 0.2, 0.2)]
        r = xy.xy_fraction(*a)
        self.assertLess(abs(r["f_ctl"]), 4 * r["f_ctl_err"])
        self.assertLess(abs(r["f_mix"]), 4 * r["f_mix_err"])


@unittest.skipUnless(os.path.exists(HITS8) and os.path.exists(C1), "faltan los impactos del run 8 o la calibración c1")
class Run8Regression(unittest.TestCase):
    def test_run8_5pe_mixed_without_ch3(self):
        G = np.full(32, np.nan)
        with open(C1) as f:
            cal = json.load(f)
        for e in cal["channels"]:
            if e["ch"] < 32 and isinstance(e.get("G_mV"), (int, float)):
                G[e["ch"]] = e["G_mV"]
        G = np.where(np.isfinite(G), G, np.nanmedian(G))
        ev = xy.load_hits(HITS8)
        Xs, Xc = xy.plane_hits(ev, G, 5, [c for c in range(16) if c != 3])
        Ys, Yc = xy.plane_hits(ev, G, 5, range(16, 32))
        r = xy.xy_fraction(Xs, Ys, Xc, Yc)
        self.assertEqual(r["N"], 44176)
        self.assertAlmostEqual(r["f_mix"], 0.0239, delta=0.00005)
        self.assertAlmostEqual(r["f_mix_err"], 0.0017, delta=0.00005)


if __name__ == "__main__":
    unittest.main(verbosity=2)
