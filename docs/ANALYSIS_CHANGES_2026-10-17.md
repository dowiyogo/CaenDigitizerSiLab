# Analysis changes, 2026-10-17

This file records the changes made on branch fix/analysis-units on 2026-10-17, the
definitions now in force, which stored outputs remain usable, and the measurements that
justified each choice. Every number is copied from a program or macro output; the output
is named next to it.

## Commits

| Commit | Change | Effect checked on run v8 |
|---|---|---|
| a164437 | kDt_ns from 1e9/65e6 (15.3846 ns) to 16.0 ns in charge_analysis, multirun_analysis, reconstruct_corrected, reconstruct_cosmic_tracks, tagprobe_analysis. Labels "ADC" to "mV", "ADC·ns" to "mV·ns". Q cuts multiplied by 1.04. | reconstruct_corrected counts unchanged: 54899 events with a hit in the muon window, 8073 hod2018 X-Y, 3 hod2019 X-Y, 486 XZ, 288 YZ, 1 three-dimensional track. The Q window went from 10 to 12 samples summed. |
| d32b58d | Sampling-step check per file, exit on mismatch. Q cuts set to the natural threshold of charge_analysis, 2318.8 mV·ns. | Q_new/Q_old per pulse, first 20000 entries: median 1.099 for 4584 hod2018X pulses, 1.101 for 248 hod2019X pulses. |
| 89858ab | Fixed verdict sentences removed from summaries. Run duration from timestamps. Expected three-dimensional tracks = kFBoth x triggers, kFBoth = 0.0167. | v8 duration printed 71.83 h (was 3600.00 h). Expected tracks 2402 (was 374). |
| ea14ad2 | pulse_shape.C added. | Superseded by later versions of the macro. |
| 25af082 | muon_timing_scan.C added; pulse_shape.C rewritten. | |
| ccc1945 | hodo_consts.h. Muon window [4000, 4400] ns in all programs. Q window per station. Median and percentiles of Q. Landau fit status and limit flag. Column order of the charge_analysis summary fixed. multirun no longer derives its window from a noise peak at 3008 ns. | reconstruct_corrected counts unchanged. kNpost18 = 8 in this commit came from a containment fraction above 1 caused by the global pedestal; replaced in eedcebc. |
| f338b93 | baseline_scan.C added; pulse_shape.C with global and local pedestal. | See Measurements. |
| eedcebc | hod2018 Q with local pedestal and window [-3, +12]; hod2019 with global pedestal and [-3, +22]. ph control cut in reconstruct_cosmic_tracks set to 35 mV. | reconstruct_corrected counts unchanged. Three-dimensional tracks with the ph control cut: 1385 before, 4 after. |
| 50ad284 | Noise of Q from channels without pulse. Two-bin edge guard on the valley search. | v8 noise of Q: 127.6 mV·ns (hod2018X), 137.4 mV·ns (hod2019X). |
| 11252cf | Q histograms to 40000 mV·ns. Valley search in [5 sigma, P10 of signal Q]. kQcut19 used by tagprobe, cosmic_tracks and multirun. | otherwise rule fired (hod2019X v8 valley invalid by edge guard, min_b=64 of bhi=65), kQcut19 = 5 x sigma_Q_noise = 687.0 mV·ns |

## Definitions in force

Sampling period: 16.0 ns, checked against time[1] - time[0] of the first event of each file.

Muon window: [4000, 4400] ns for both stations and all runs v1-v10. In muon_timing_scan.C
(first 3000 entries of each run) the hod2019 pulses above 35 mV peak between 4140 and
4204 ns in every run that has them, and none falls in [3000, 3400] ns.

Pulse height: ph = ped - min(wf), with ped the mean of the samples at t < 2500 ns.
Thresholds 50 mV for hod2018 and 35 mV for hod2019, that is 102 and 72 ADC counts at
2000/4095 mV per count. Before 2026-10-17 the code and figures labelled these values "ADC".

Charge: Q = 16 ns x sum of (ped_Q - wf_i) over samples -3 to +kNpost around the minimum.
For hod2018, ped_Q is the mean of samples -12 to -4 from the minimum and kNpost = 12.
For hod2019, ped_Q = ped and kNpost = 22.

Noise of Q: RMS of Q integrated at t0 = 4200 ns over channels with ph < 5 mV in
[2500, 8000] ns, with the same pedestal and window as the signal.

Q cut for hod2019 in tagprobe_analysis, reconstruct_cosmic_tracks and multirun_analysis:
kQcut19 = 687.0 mV·ns, 5 x noise of Q (hod2019X v8 valley invalid by edge guard: min_b=64 of bhi=65).

## Measurements

Run v8, first 20000 entries.

Baseline, baseline_scan.C: mean over channel-events of (mean of samples in [3700, 3950] ns
minus ped), and its RMS over channel-events: hod2018X +0.016 mV and 4.71 mV; hod2018Y
+0.001 mV and 4.12 mV; hod2019X +0.001 mV and 0.35 mV; hod2019Y -0.001 mV and 0.37 mV.
The cause of the larger spread in hod2018 was not determined.

Noise of Q with the two pedestals, baseline_scan.C step C, channels without pulse:

| Station | Window | RMS, global pedestal [mV·ns] | RMS, local pedestal [mV·ns] |
|---|---|---|---|
| hod2018 | [-3, +8] | 180.03 | 100.68 |
| hod2018 | [-3, +12] | 234.08 | 136.55 |
| hod2019 | [-3, +22] | 129.87 | 155.76 |

Average pulse, pulse_shape.C, outputs/pulse_shape_20261017_164700, 4584 hod2018X pulses
with ph > 50 mV and 248 hod2019X pulses with ph > 35 mV, minimum in [4000, 4400] ns.
hod2018X with local pedestal crosses zero at +19 samples (304 ns); the integral up to the
crossing is 13761.9 mV·ns, of which [-3, +8] contains 0.836 and [-3, +12] 0.950. After the
crossing the average stays negative and reaches -14.4 mV at +40 samples (640 ns).
hod2019X with global pedestal does not cross zero up to +40 samples; the integral over
[-3, +40] is 10951.5 mV·ns, of which [-3, +8] contains 0.709 and [-3, +22] 0.954.

Q with the final definitions, outputs/charge_analysis_20261017_181925:
```
Métrica                                18X-v7        18X-v8        19X-v7        19X-v8
σ_Q noise [mV·ns] (no-pulse channels) 138.8         127.6         125.8         137.4         
N pulsos señal                         3402          33082         219           2087          
Mediana Q [mV·ns]                      9579.2        9489.7        7933.7        8316.8          (Q of pulses above the ph threshold; not a gain estimator)
Percentil 10 Q [mV·ns]                 6479.2        6387.8        5447.7        5461.5        
Percentil 90 Q [mV·ns]                 22769.5       22545.3       17121.6       17612.8       
```

## Validity of stored outputs

| Outputs | Status |
|---|---|
| Directories dated before 2026-10-17 | Counts that depend only on ph and time (hits, coincidences, tracks, timing peaks) are valid. Q values used 15.38 ns and 10 samples and cannot be compared with current Q. "ADC" in labels means mV. |
| multirun_* before 2026-10-17 | Not valid: muon window of v5-v10 taken from a noise peak at 3008 ns; Landau MPV at the 4500 mV·ns limit. |
| charge_analysis_* before 2026-10-17 | Columns 18X-v8 and 19X-v7 swapped in the summary. |
| outputs/_superseded/ | Wrong hod2019 channel mapping. |
| corrected_tracking_20261017_* | Valid. Counts identical in every run of this date. |
| figs_for_beamer_20261017_103956 | Valid. Uses ph and time only. |
| charge_analysis_, multirun_, tagprobe_, cosmic_tracks_ of 2026-10-17 before charge_analysis_20261017_181925, multirun_20261017_182116, tagprobe_20261017_183028, cosmic_tracks_20261017_183153 | Superseded intermediate states. |
| charge_analysis_20261017_181925, multirun_20261017_182116, tagprobe_20261017_183028, cosmic_tracks_20261017_183153 | Current. |
| pulse_shape_20261017_164700 | Current. Other pulse_shape_ directories of this date come from earlier versions of the macro. |

## Not determined

Bias voltages of runs v1-v10. Cause of the baseline spread and of the post-pulse overshoot
in hod2018.

kQcut19 = 687.0 mV·ns is provisional. With it, reconstruct_cosmic_tracks gives 6125 Q+r
three-dimensional tracks in run v8 (outputs/cosmic_tracks_20261017_183153), against 2402 expected
with full efficiency of both stations.
