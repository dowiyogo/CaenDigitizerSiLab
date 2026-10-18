#!/usr/bin/env python3
"""
analysis/scripts/hv_diag6_shape.py

TASK HV-DIAG-6 — shape of the 2026-05-29 group 1/2 gain loss (existing CSVs only)

INPUTS:
  - outputs/hv_diag5_rates_10min.csv
  - outputs/hv_diag3_gain_vs_time.csv

OUTPUTS:
  - Table 1 (stdout): time r0 r6 r20 r22 r31 rmed_other18
  - Table 2 (stdout): hour G0 G6 G20 G22 G31 Gmed_other18 status_flags
  - Plot: outputs/hv_diag6_shape.png
"""

import sys
import os
import csv
from collections import defaultdict
import statistics
import datetime

# Import ROOT in batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# Constants
HOD2018_OTHER = [
    1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30
]

RATES_CSV = "outputs/hv_diag5_rates_10min.csv"
GAIN_CSV  = "outputs/hv_diag3_gain_vs_time.csv"
OUT_PNG   = "outputs/hv_diag6_shape.png"

def main():
    if not os.path.exists(RATES_CSV):
        sys.exit(f"Error: missing input file {RATES_CSV}")
    if not os.path.exists(GAIN_CSV):
        sys.exit(f"Error: missing input file {GAIN_CSV}")

    # =========================================================================
    # 1. Compute v8 Night Reference Rates R_c from outputs/hv_diag5_rates_10min.csv
    # =========================================================================
    v8_night_rates = defaultdict(list)
    with open(RATES_CSV, mode="r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            if r["run"] == "v8":
                hh = int(r["time"].split(":")[0])
                ch = int(r["ch"])
                rate = float(r["rate"])
                k = int(r["k"])
                n_ev = int(r["n_ev"])
                if 0 <= hh <= 5 and k < n_ev and rate >= 0:
                    v8_night_rates[ch].append(rate)

    Rc = {}
    for ch in range(64):
        rates = sorted(v8_night_rates[ch])
        if rates:
            n = len(rates)
            med = rates[n // 2] if n % 2 == 1 else 0.5 * (rates[n // 2 - 1] + rates[n // 2])
            if med <= 0.001:
                pos = [x for x in rates if x > 0.001]
                if pos:
                    npos = len(pos)
                    med = pos[npos // 2] if npos % 2 == 1 else 0.5 * (pos[npos // 2 - 1] + pos[npos // 2])
            Rc[ch] = med
        else:
            Rc[ch] = 0.0

    # =========================================================================
    # 2. Extract Table 1 data: run v8, 2026-05-29 06:00 to 2026-05-30 02:00
    # =========================================================================
    rates_bins = defaultdict(dict)
    with open(RATES_CSV, mode="r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            if r["run"] == "v8":
                dt = f"{r['date']} {r['time']}"
                if "2026-05-29 06:00" <= dt <= "2026-05-30 02:00":
                    rates_bins[dt][int(r["ch"])] = float(r["rate"])

    sorted_dts = sorted(rates_bins.keys())

    # Build Table 1 rows
    table1_rows = []
    g_r0  = ROOT.TGraph()
    g_r6  = ROOT.TGraph()
    g_r20 = ROOT.TGraph()
    g_r22 = ROOT.TGraph()
    g_r31 = ROOT.TGraph()
    g_rmed = ROOT.TGraph()

    for idx, dt in enumerate(sorted_dts):
        d_str, t_str = dt.split()
        ch_rates = rates_bins[dt]

        r0  = ch_rates[0]  / Rc[0]  if Rc[0]  > 0 else 0.0
        r6  = ch_rates[6]  / Rc[6]  if Rc[6]  > 0 else 0.0
        r20 = ch_rates[20] / Rc[20] if Rc[20] > 0 else 0.0
        r22 = ch_rates[22] / Rc[22] if Rc[22] > 0 else 0.0
        r31 = ch_rates[31] / Rc[31] if Rc[31] > 0 else 0.0

        r_others = [ch_rates[c] / Rc[c] for c in HOD2018_OTHER if Rc[c] > 0 and c in ch_rates]
        rmed = statistics.median(r_others) if r_others else 0.0

        # Avoid negative zero formatting
        r0 = 0.0 if abs(r0) < 1e-9 else r0
        r6 = 0.0 if abs(r6) < 1e-9 else r6
        r20 = 0.0 if abs(r20) < 1e-9 else r20
        r22 = 0.0 if abs(r22) < 1e-9 else r22
        r31 = 0.0 if abs(r31) < 1e-9 else r31
        rmed = 0.0 if abs(rmed) < 1e-9 else rmed

        table1_rows.append((t_str, r0, r6, r20, r22, r31, rmed))

        y, m, d = [int(x) for x in d_str.split("-")]
        hh, mm = [int(x) for x in t_str.split(":")]
        ts = ROOT.TDatime(y, m, d, hh, mm, 0).Convert()

        g_r0.SetPoint(idx, ts, max(0.0, r0))
        g_r6.SetPoint(idx, ts, max(0.0, r6))
        g_r20.SetPoint(idx, ts, max(0.0, r20))
        g_r22.SetPoint(idx, ts, max(0.0, r22))
        g_r31.SetPoint(idx, ts, max(0.0, r31))
        g_rmed.SetPoint(idx, ts, max(0.0, rmed))

    # =========================================================================
    # 3. Extract Table 2 data: run v8, 2026-05-29 05:00 to 23:00
    # =========================================================================
    gain_hours = defaultdict(dict)
    with open(GAIN_CSV, mode="r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            if r["run"] == "v8" and r["date"] == "2026-05-29":
                hh = int(r["hour"].split(":")[0])
                if 5 <= hh <= 23:
                    gain_hours[r["hour"]][int(r["ch"])] = r

    sorted_hours = sorted(gain_hours.keys())

    table2_rows = []
    g_G0  = ROOT.TGraph()
    g_G6  = ROOT.TGraph()
    g_G20 = ROOT.TGraph()
    g_G22 = ROOT.TGraph()
    g_G31 = ROOT.TGraph()
    g_Gmed = ROOT.TGraph()

    for idx, h in enumerate(sorted_hours):
        ch_map = gain_hours[h]
        g0  = float(ch_map[0]["G"])
        g6  = float(ch_map[6]["G"])
        g20 = float(ch_map[20]["G"])
        g22 = float(ch_map[22]["G"])
        g31 = float(ch_map[31]["G"])

        g_others = [float(ch_map[c]["G"]) for c in HOD2018_OTHER if c in ch_map]
        gmed = statistics.median(g_others) if g_others else 0.0

        # Channels with fit status != 0
        bad_chs = [str(c) for c in sorted(ch_map.keys()) if int(ch_map[c]["status"]) != 0]
        status_flags = ",".join(bad_chs) if bad_chs else "none"

        table2_rows.append((h, g0, g6, g20, g22, g31, gmed, status_flags))

        hh = int(h.split(":")[0])
        ts = ROOT.TDatime(2026, 5, 29, hh, 0, 0).Convert()

        g_G0.SetPoint(idx, ts, g0)
        g_G6.SetPoint(idx, ts, g6)
        g_G20.SetPoint(idx, ts, g20)
        g_G22.SetPoint(idx, ts, g22)
        g_G31.SetPoint(idx, ts, g31)
        g_Gmed.SetPoint(idx, ts, gmed)

    # =========================================================================
    # 4. Print Table 1 (stdout)
    # =========================================================================
    print("time  r0  r6  r20  r22  r31  rmed_other18")
    for row in table1_rows:
        t_str, r0, r6, r20, r22, r31, rmed = row
        print(f"{t_str:<5}  {r0:>7.4f}  {r6:>7.4f}  {r20:>7.4f}  {r22:>7.4f}  {r31:>7.4f}  {rmed:>12.4f}")

    print()

    # =========================================================================
    # 5. Print Table 2 (stdout)
    # =========================================================================
    print("hour  G0  G6  G20  G22  G31  Gmed_other18  status_flags")
    for row in table2_rows:
        h, g0, g6, g20, g22, g31, gmed, st_flags = row
        print(f"{h:<5}  {g0:>6.4f}  {g6:>6.4f}  {g20:>6.4f}  {g22:>6.4f}  {g31:>6.4f}  {gmed:>12.4f}  {st_flags}")

    # =========================================================================
    # 6. Generate Plot outputs/hv_diag6_shape.png
    # =========================================================================
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    c_shape = ROOT.TCanvas("c_shape", "HV-DIAG-6 Shape", 1300, 950)
    c_shape.Divide(1, 2, 0.001, 0.001)

    colors  = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kMagenta+2, ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kBlack]
    markers = [20, 21, 22, 23, 33, 24]

    # --- Pad 1: Ratios vs Time ---
    pad1 = c_shape.cd(1)
    pad1.SetGrid()
    pad1.SetLeftMargin(0.10)
    pad1.SetRightMargin(0.06)
    pad1.SetTopMargin(0.10)
    pad1.SetBottomMargin(0.12)

    mg1 = ROOT.TMultiGraph()
    graphs1 = [g_r0, g_r6, g_r20, g_r22, g_r31, g_rmed]
    names1  = ["r0 (ch 0)", "r6 (ch 6)", "r20 (ch 20)", "r22 (ch 22)", "r31 (ch 31)", "rmed_other18"]

    for g, col, m in zip(graphs1, colors, markers):
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetMarkerStyle(m)
        g.SetMarkerSize(0.9 if col != ROOT.kBlack else 1.0)
        g.SetLineWidth(2 if col != ROOT.kBlack else 3)
        if col == ROOT.kBlack:
            g.SetLineStyle(2)
        mg1.Add(g, "LP")

    mg1.Draw("A")
    mg1.SetTitle("Rate Ratios vs Time (2026-05-29 06:00 to 2026-05-30 02:00);Time (Local UTC-4);Rate / R_{c}")
    mg1.GetYaxis().SetRangeUser(0.0, 1.8)
    ax1 = mg1.GetXaxis()
    ax1.SetTimeDisplay(1)
    ax1.SetTimeFormat("%m-%d %H:%M")
    ax1.SetTimeOffset(0)
    ax1.SetLabelSize(0.038)
    ax1.SetTitleSize(0.042)
    mg1.GetYaxis().SetLabelSize(0.038)
    mg1.GetYaxis().SetTitleSize(0.042)

    # Reference dashed line at ratio = 1.0
    ts_min1 = ROOT.TDatime(2026, 5, 29, 6, 0, 0).Convert()
    ts_max1 = ROOT.TDatime(2026, 5, 30, 2, 0, 0).Convert()
    line1 = ROOT.TLine(ts_min1, 1.0, ts_max1, 1.0)
    line1.SetLineColor(ROOT.kGray+2)
    line1.SetLineStyle(3)
    line1.SetLineWidth(2)
    line1.Draw("SAME")

    leg1 = ROOT.TLegend(0.62, 0.65, 0.93, 0.88)
    leg1.SetNColumns(2)
    leg1.SetBorderSize(1)
    leg1.SetFillColor(ROOT.kWhite)
    leg1.SetTextSize(0.036)
    for g, name in zip(graphs1, names1):
        leg1.AddEntry(g, name, "LP")
    leg1.Draw("SAME")

    # --- Pad 2: Gains vs Time ---
    pad2 = c_shape.cd(2)
    pad2.SetGrid()
    pad2.SetLeftMargin(0.10)
    pad2.SetRightMargin(0.06)
    pad2.SetTopMargin(0.10)
    pad2.SetBottomMargin(0.12)

    mg2 = ROOT.TMultiGraph()
    graphs2 = [g_G0, g_G6, g_G20, g_G22, g_G31, g_Gmed]
    names2  = ["G0 (ch 0)", "G6 (ch 6)", "G20 (ch 20)", "G22 (ch 22)", "G31 (ch 31)", "Gmed_other18"]

    for g, col, m in zip(graphs2, colors, markers):
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetMarkerStyle(m)
        g.SetMarkerSize(1.0)
        g.SetLineWidth(2 if col != ROOT.kBlack else 3)
        if col == ROOT.kBlack:
            g.SetLineStyle(2)
        mg2.Add(g, "LP")

    mg2.Draw("A")
    mg2.SetTitle("Hourly Gain G vs Time (2026-05-29 05:00 to 23:00);Time (Local UTC-4);Gain G [mV/p.e.]")
    mg2.GetYaxis().SetRangeUser(3.0, 5.8)
    ax2 = mg2.GetXaxis()
    ax2.SetTimeDisplay(1)
    ax2.SetTimeFormat("%m-%d %H:%M")
    ax2.SetTimeOffset(0)
    ax2.SetLabelSize(0.038)
    ax2.SetTitleSize(0.042)
    mg2.GetYaxis().SetLabelSize(0.038)
    mg2.GetYaxis().SetTitleSize(0.042)

    leg2 = ROOT.TLegend(0.62, 0.65, 0.93, 0.88)
    leg2.SetNColumns(2)
    leg2.SetBorderSize(1)
    leg2.SetFillColor(ROOT.kWhite)
    leg2.SetTextSize(0.036)
    for g, name in zip(graphs2, names2):
        leg2.AddEntry(g, name, "LP")
    leg2.Draw("SAME")

    c_shape.SaveAs(OUT_PNG)

if __name__ == "__main__":
    main()
