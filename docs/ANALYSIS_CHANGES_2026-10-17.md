# Analysis Changes and Calibration Record (2026-10-17)

This document records the modifications made to the offline analysis software of the CAEN VX1740 digitizer (hod2018 and hod2019 hodoscopes) on the `fix/analysis-units` branch on 2026-10-17, the final definitions adopted, the validity assessment of stored analysis outputs, and the empirical measurements obtained from diagnostic macros.

---

## 1. Commit History on `fix/analysis-units` (since `9e493ed`)

| Commit | Summary of Changes | Affected Programs | Effect on Analysis Outputs |
|---|---|---|---|
| `a164437` | Corrected sampling period $kDt\_ns$ to 16.0 ns (was 15.3846 ns). Corrected unit labels from ADC counts to mV / mV·ns. Added runtime step check. Clarified threshold conversions in comments (50 mV $\approx$ 102 ADC, 35 mV $\approx$ 72 ADC). | `charge_analysis.cxx`, `multirun_analysis.cxx`, `reconstruct_corrected.cxx`, `reconstruct_cosmic_tracks.cxx`, `tagprobe_analysis.cxx`, `gen_beamer_figs.C` | Track reconstruction counts in `reconstruct_corrected` strictly unchanged ($N_{\text{tracks}}=4443$, $N_{\text{golden}}=1386$, $N_{3D}=1111$, $N_{\text{both}}=1000$) because ph cuts are in mV. Integrated charge Q values rescaled by factor $\sim 1.04$ due to $\Delta t$ and integration window. |
| `d32b58d` | Replaced static step check with per-file check `bool checked_step = false;` placed before event loops with `std::exit(1)` on mismatch. Updated Q cuts in tracking and tag-and-probe from recomputed natural threshold (2372.5 mV·ns for v8 with 10 samples). | `charge_analysis.cxx`, `multirun_analysis.cxx`, `reconstruct_corrected.cxx`, `reconstruct_cosmic_tracks.cxx`, `tagprobe_analysis.cxx` | Runtime verification ensured step is exactly 16.0 ns. Q cuts aligned with natural threshold. |
| `89858ab` | Removed hardcoded verdicts and text interpretations (e.g., "VEREDICTO", "LIMITADO", "Vikuiti", hardcoded SNR table). Replaced hardcoded nominal duration with duration calculated from ROOT timestamps (`ts_last - ts_first`). Replaced hardcoded expected $N_{3D}$ with $kFBoth \times N_{\text{golden}}$. | `analyze_highgain.cxx`, `charge_analysis.cxx`, `infer_mapping.cxx`, `multirun_analysis.cxx`, `poisson_analysis.cxx`, `reconstruct_corrected.cxx`, `reconstruct_cosmic_tracks.cxx`, `reconstruct_v8b.cxx`, `tagprobe_analysis.cxx`, `gen_beamer_figs.C` | Summaries report only computed numerical quantities. Run v8 duration corrected from 70.2 h (nominal) to 65.6 h (timestamp difference: 236208 s). Expected $N_{3D}$ dynamically computed. |
| `ea14ad2` | Added pulse-shape diagnosis macro `analysis/macros/pulse_shape.C`. | `analysis/macros/pulse_shape.C` | Initial diagnostic run for pulse shape inspection. |
| `25af082` | Added `analysis/macros/muon_timing_scan.C` and updated `pulse_shape.C` for reproducible timing and pulse-shape analyses. | `analysis/macros/muon_timing_scan.C`, `analysis/macros/pulse_shape.C` | Identified true muon coincidence peak at $t \approx 4175$ ns in both hodoscopes. Identified electronic pickup in hod2018 at $t \approx 3525$ ns. Showed hod2018 pulse returns to baseline around +8 samples with undershoot, while hod2019 pulse lasts $\sim +22$ samples. |
| `ccc1945` | Set fixed muon window $[4000, 4400]$ ns for both hodoscopes across all analysis programs. Introduced per-station Q integration window ($kNpost18=8$, $kNpost19=22$). Introduced robust Q estimator (median, 10th and 90th percentiles for pulses with $ph >$ threshold). Added Landau fit parameter limit check (`AT_LIMIT`). Fixed swapped column order in `charge_analysis` summary. Fixed noise-derived muon window ($t=3008$ ns) for v5-v10 in `multirun_analysis`. | `charge_analysis.cxx`, `multirun_analysis.cxx`, `reconstruct_corrected.cxx`, `reconstruct_cosmic_tracks.cxx`, `tagprobe_analysis.cxx`, `hodo_consts.h` | `reconstruct_corrected` track counts strictly identical ($N_{\text{tracks}}=4443$, $N_{\text{golden}}=1386$, $N_{3D}=1111$, $N_{\text{both}}=1000$). `charge_analysis` summary columns correctly ordered (18X-v7, 18X-v8, 19X-v7, 19X-v8) with medians 9474.0, 9385.4, 7933.7, 8316.8 mV·ns. `multirun_analysis` peak times for v5-v10 corrected to 4128-4198 ns. |
| `f338b93` | Added macro `analysis/macros/baseline_scan.C` to measure baseline shifts across 9 time windows and evaluate Q noise with global vs local pedestal. Updated `pulse_shape.C` to support dual pedestal definitions. | `analysis/macros/baseline_scan.C`, `analysis/macros/pulse_shape.C` | Measured significant baseline shift and pickup in hod2018 (local ped reduces Q RMS by $44\%$), but negligible shift in hod2019 (global ped preferred). |
| `eedcebc` | Implemented empirical baseline and window decisions: hod2018 uses local pedestal `ped_loc` ($[-12, -4]$ samples from peak) with window $[-3, +12]$ ($kNpost18=12$, 95% containment); hod2019 keeps global pedestal `ped_glob` ($t < 2500$ ns) with window $[-3, +22]$ ($kNpost19=22$, 95% containment). Updated `hodo_consts.h` (`kQLocalPed18=true`, `kQLocalPed19=false`, `kNpost18=12`, `kNpost19=22`). | `hodo_consts.h`, `reconstruct_corrected.cxx`, `reconstruct_cosmic_tracks.cxx`, `charge_analysis.cxx`, `multirun_analysis.cxx`, `AGENTS.md` | `reconstruct_corrected` track counts strictly unchanged ($N_{\text{tracks}}=4443$, $N_{\text{golden}}=1386$, $N_{3D}=1111$, $N_{\text{both}}=1000$). `charge_analysis` median Q for hod2018 changed from 9385.4 to 9489.7 mV·ns (v8) due to $kNpost18=12$. |
| *Current* | Replaced max-Q Gaussian fit $\sigma_Q$ in `charge_analysis` and `multirun_analysis` with true pedestal noise $\sigma_Q$ from no-pulse channels ($ph < 5$ mV in $[2500, 8000]$ ns) integrated over $[t_0-3, t_0+kNpost]$ at $t_0=4200$ ns using each station's pedestal. Added natural threshold boundary guard (`min_b <= blo+2 || min_b >= bhi-2` rejected as "no valley"). Printed `blo, bhi, min_b` next to each natural threshold value. | `charge_analysis.cxx`, `multirun_analysis.cxx`, `docs/ANALYSIS_CHANGES_2026-10-17.md`, `AGENTS.md` | `charge_analysis`: $\sigma_Q$ noise for v8 is 127.6 mV·ns (18X) and 137.4 mV·ns (19X), replacing artificial 1763.3 and 140.5 mV·ns. `multirun_analysis`: $\sigma_Q$ noise for v8 is 137.4 mV·ns (19X) instead of 35.5 mV·ns; edge-artifact valleys in v8 and v10 correctly rejected; global efficiency $\varepsilon_{\text{true, global}} = +0.0053\% \pm 0.0098\%$, upper bound $< 0.02\%$ (95% CL). |

---

## 2. Final Definitions

### 2.1 Sampling Period
- $kDt\_ns = 16.0$ ns, matching the CAEN VX1740 sampling rate of 62.5 MS/s (one sample every 16.0 ns, as recorded in `data_dig.time`).
- Enforced at runtime in all analysis programs with an entry-level check against `time[1] - time[0]`.

### 2.2 Muon Timing Window
- Muon coincidence window: $[4000, 4400]$ ns for both hod2018 and hod2019.
- Peak muon time $t_{\text{peak}} \approx 4175$ ns.
- The window $[3400, 3700]$ ns contains electronic pickup in hod2018 ($t \approx 3525$ ns) and is excluded from signal selection.

### 2.3 Pulse Height (ph)
- Defined as:
  $$ph = ped\_glob - \min_{t \in [t_{\text{lo}}, t_{\text{hi}}]} V(t)$$
  where $ped\_glob$ is the global baseline computed as the mean of waveform samples with $t < 2500$ ns.
- Diagnostic pulse thresholds:
  - **hod2018**: $ph > 50.0$ mV ($\approx 102$ ADC counts for a 2 Vpp 12-bit digitizer, $0.4884$ mV/count).
  - **hod2019**: $ph > 35.0$ mV ($\approx 72$ ADC counts for a 2 Vpp 12-bit digitizer, $0.4884$ mV/count).

### 2.4 Integrated Charge (Q) per Station
Integrated charge is calculated on negative polarity signals as:
$$Q = \sum_{i = i_{\text{peak}} - kNpre}^{i_{\text{peak}} + kNpost} (ped\_Q - V(t_i)) \times kDt\_ns$$
where $kNpre = 3$ (48 ns pre-peak). Station-specific parameters:
- **hod2018**:
  - Pedestal $ped\_Q = ped\_loc$: local baseline computed as the mean of 9 samples in $[-12, -4]$ relative to peak sample ($[i_{\text{peak}} - 12, i_{\text{peak}} - 4]$).
  - Integration window: $[-3, +12]$ samples ($kNpost18 = 12$, corresponding to $+192$ ns post-peak).
  - Containment: $95\%$ containment up to waveform zero crossing at $n_0 = +19$ samples ($+304$ ns).
- **hod2019**:
  - Pedestal $ped\_Q = ped\_glob$: global baseline computed as the mean of samples with $t < 2500$ ns.
  - Integration window: $[-3, +22]$ samples ($kNpost19 = 22$, corresponding to $+352$ ns post-peak).
  - Containment: $95\%$ containment up to waveform zero crossing at $n_0 = +40$ samples ($+640$ ns).

### 2.5 Pedestal Noise of Q ($\sigma_Q$ noise)
- Defined as the sample standard deviation (RMS around the mean) of integrated charge evaluated on channel-events with no pulse:
  $$ph < 5.0\text{ mV in } [2500, 8000]\text{ ns}$$
- Integrated over $[t_0 - 3, t_0 + kNpost]$ at fixed sample $t_0 = 4200$ ns (sample index 262), using each station's baseline definition ($ped\_loc$ for hod2018, $ped\_glob$ for hod2019).
- Replaces earlier Gaussian fits to the max-Q spectrum and uncorrelated $\sigma \times \sqrt{N}$ estimates.

---

## 3. Validity of Stored Analysis Outputs

| Directory Pattern | Status | Reason |
|---|---|---|
| All directories dated before `20261017` | **Invalid** | Hardcoded wrong sampling period (15.3846 ns instead of 16.0 ns); pulse heights and charges labeled as "ADC" when values were in mV / mV·ns; some early runs used uncorrected sequential channel mapping for hod2019. |
| `outputs/timesync_test_*` | **Valid with caveat** | Established coincidence timing peak at $t \approx 4175$ ns and pickup at 3525 ns. However, labels reading "ADC" are actually in mV, and sampling period was 15.3846 ns. |
| `outputs/figs_for_beamer_*` | **Invalid** | Hardcoded interpretations, legacy units ("ADC"), legacy Q integration windows, and nominal durations. |
| `outputs/charge_analysis_*` (before 2026-10-17) | **Invalid** | Swapped columns between 18X-v8 and 19X-v7; legacy sampling period; wrong integration windows. |
| `outputs/charge_analysis_20261017_165730` | **Valid with caveat** | Correct 16.0 ns step and per-station windows ($kNpost18=12$, $kNpost19=22$), but used max-Q Gaussian fit for $\sigma_Q$ ($\sim 1763$ mV·ns) and reported an edge-artifact valley for 18X-v8 at 12377.5 mV·ns. Superseded by `charge_analysis_20261017_174137`. |
| `outputs/charge_analysis_20261017_174137` | **Valid** | Final analysis with per-station windows, local pedestal for hod2018, true $\sigma_Q$ noise from no-pulse channels, and boundary edge guard on valley search. |
| `outputs/multirun_*` (before commit `ccc1945`) | **Invalid** | Peak times for runs v5-v10 derived from noise ($t=3008$ ns); Landau MPV pinned at limit 4500 mV·ns; legacy sampling period. |
| `outputs/multirun_20261017_170204` | **Valid with caveat** | Correct 16.0 ns period and muon window, but $\sigma_Q$ noise was estimated from sample RMS without noise correlation ($\sim 35.5$ mV·ns for v8) and edge valley artifacts in v8/v10 were not filtered. Superseded by `multirun_20261017_174334`. |
| `outputs/multirun_20261017_174334` | **Valid** | Final multi-run analysis with per-station windows, true $\sigma_Q$ noise from no-pulse channels ($137.4$ mV·ns for v8), and valley boundary edge guard. |
| `outputs/tagprobe_*` (before 2026-10-17) | **Invalid** | Legacy sampling period and legacy Q definition. |
| `outputs/tagprobe_20261017_165913` | **Valid** | 16.0 ns period, corrected channel mapping, and natural Q threshold cut (2401.2 mV·ns). |
| `outputs/corrected_tracking_*` (before 2026-10-17) | **Valid with caveat** | Reconstructed tracks, cluster distributions, and angular distributions are geometrically correct because ph cuts in mV are independent of $\Delta t$. However, unit labels state "ADC" instead of mV, and run duration was based on nominal rate. |
| `outputs/corrected_tracking_20261017_165604` | **Valid** | 16.0 ns period, mV units, timestamp-derived duration, and verified track counts. |
| `outputs/cosmic_tracks_*` (before 2026-10-17) | **Valid with caveat** | Reconstructed tracks and angular fits geometrically correct, but unit labels state "ADC" instead of mV. |
| `outputs/cosmic_tracks_20261017_170037` | **Valid** | 16.0 ns period, mV units, and updated Q cuts. |
| `outputs/pulse_shape_20261017_151734` | **Invalid** | Contained pickup and noise selection artifacts. |
| `outputs/pulse_shape_20261017_164700` | **Valid** | Verified pulse shape comparison between `ped_glob` and `ped_loc` on clean muon sample. |

---

## 4. Measured Facts from Diagnostic Macros

### 4.1 Baseline Shift Analysis (`baseline_scan.C`, run v8, first 20,000 entries)
The baseline shift $b\_shift$ is defined as the difference between the mean waveform voltage in the pre-muon interval $[3700, 3950)$ ns and the global pedestal $ped\_glob$ ($t < 2500$ ns).

- **All triggers (no pulse selection)**:
  - **hod2018X**: Mean = $+0.0156$ mV, Event-level RMS = $1.1727$ mV, Channel-event RMS = $4.7128$ mV.
  - **hod2018Y**: Mean = $+0.0007$ mV, Event-level RMS = $1.0327$ mV, Channel-event RMS = $4.1164$ mV.
  - **hod2019X**: Mean = $+0.0008$ mV, Event-level RMS = $0.0949$ mV, Channel-event RMS = $0.3549$ mV.
  - **hod2019Y**: Mean = $-0.0006$ mV, Event-level RMS = $0.0942$ mV, Channel-event RMS = $0.3679$ mV.

- **No-pulse triggers ($ph < 5.0$ mV in $[2500, 8000]$ ns)**:
  - **hod2018X**: Mean = $+0.1447$ mV, Event-level RMS = $0.5843$ mV, Channel-event RMS = $1.0163$ mV.
  - **hod2018Y**: Mean = $+0.4922$ mV, Event-level RMS = $1.3292$ mV, Channel-event RMS = $1.3045$ mV.
  - **hod2019X**: Mean = $+0.0031$ mV, Event-level RMS = $0.0935$ mV, Channel-event RMS = $0.3454$ mV.
  - **hod2019Y**: Mean = $+0.0042$ mV, Event-level RMS = $0.0900$ mV, Channel-event RMS = $0.3458$ mV.

The physical origin of the baseline shift in hod2018 was not recorded; cause not determined.

### 4.2 Pedestal Q Noise Comparison (`baseline_scan.C`, Step C, run v8, first 20,000 entries)
Integrated charge noise evaluated on no-pulse channel-events ($ph < 5.0$ mV in $[2500, 8000]$ ns) at $t_0 = 4200$ ns:

- **hod2018 ($kNpost=8$, $[-3, +8]$ samples)**:
  - $\text{RMS}_Q(ped\_glob) = 180.03$ mV·ns
  - $\text{RMS}_Q(ped\_loc) = 100.68$ mV·ns
  - Ratio $\text{RMS}_{\text{loc}} / \text{RMS}_{\text{glob}} = 0.5592$ (44.1% noise reduction with local pedestal).

- **hod2018 ($kNpost=12$, $[-3, +12]$ samples)**:
  - $\text{RMS}_Q(ped\_glob) = 234.08$ mV·ns
  - $\text{RMS}_Q(ped\_loc) = 120.14$ mV·ns (18X: $136.55$ mV·ns; 18Y: $100.68$ mV·ns)
  - Ratio $\text{RMS}_{\text{loc}} / \text{RMS}_{\text{glob}} = 0.5132$ (48.7% noise reduction with local pedestal).

- **hod2019 ($kNpost=22$, $[-3, +22]$ samples)**:
  - $\text{RMS}_Q(ped\_glob) = 129.87$ mV·ns (19X: $129.87$ mV·ns; 19Y: $130.68$ mV·ns)
  - $\text{RMS}_Q(ped\_loc) = 155.76$ mV·ns (19X: $155.76$ mV·ns; 19Y: $156.40$ mV·ns)
  - Ratio $\text{RMS}_{\text{loc}} / \text{RMS}_{\text{glob}} = 1.1993$ (local pedestal increases noise by 19.9% due to sample variance).

### 4.3 Pulse Shape Containment and Zero-Crossing (`pulse_shape.C`, output `outputs/pulse_shape_20261017_164700`)
Evaluated on run v8 muon candidates ($ph > 50$ mV for hod2018, $ph > 35$ mV for hod2019, $t_{\text{peak}} \in [4000, 4400]$ ns):

- **hod2018X with local pedestal ($ped\_loc$)**:
  - Total integral up to zero crossing: $13761.92$ mV·ns.
  - Zero-crossing sample: $n_0 = +19$ samples ($+304$ ns).
  - 95% charge containment sample: $n_{95} = +12$ samples ($+192$ ns) with integral $13204.64$ mV·ns ($95.95\%$ of total).
  - Adopted window: $[-3, +12]$ samples ($kNpost18 = 12$).

- **hod2019X with global pedestal ($ped\_glob$)**:
  - Total integral up to zero crossing: $10951.52$ mV·ns.
  - Zero-crossing sample: $n_0 = +40$ samples ($+640$ ns).
  - 95% charge containment sample: $n_{95} = +22$ samples ($+352$ ns) with integral $10471.26$ mV·ns ($95.61\%$ of total).
  - Adopted window: $[-3, +22]$ samples ($kNpost19 = 22$).
