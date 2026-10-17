# AGENTS.md — CaenDigitizerSiLab DAQ & Analysis Environment

> **Document Status**: Operational Guide for Autonomous Agents & Developers  
> **Repository Path**: `/home/xplus/digitizer/CaenDigitizerSiLab-jai`  
> **Last Reorganization**: 2026-10-17 (Branch: `snapshot/2026-10-17`)  
> **Host**: Linux / systemd (x86_64), ROOT 6.30/04, GCC 11.4.0  

---

## 1. Scope, Purpose & Architecture

This repository drives a CAEN digitizer system and provides the offline reconstruction and analysis pipeline for a muon hodoscope detector (plastic scintillator bars read out by SiPM arrays, used for detector characterization and cosmic-muon tracking).

The codebase is split strictly into:
1. **DAQ Engine (Root level)**: C++ digitizer driver, CAENComm / CAENDigitizer wrappers, acquisition binaries (`test_acquire_ext_nim`, `test_acquire_sw`, `test_open`), and hardware configuration files.
2. **Analysis Pipeline (`analysis/`)**: Offline event analysis, multi-plane cosmic track reconstruction, pulse charge and time integration, channel health monitoring, and SPE gain calibration.
3. **Run Catalog (`runs/`)**: Production run metadata (`runs/runs_2026.csv`), including calibrated wall-clock intervals, entry counts, and detector status.
4. **Data Isolation (`/home/xplus/digitizer/data/raw/`)**: All raw experimental data files are stored outside the Git repository in immutable, read-only directories, accessible via symlinks in the repository root.

---

## 2. Rules of Engagement for Autonomous Agents

Any autonomous agent operating in this repository must strictly adhere to the following rules:

> [!CAUTION]
> **Data Immutability**: NEVER attempt to delete, rename, or modify files in `/home/xplus/digitizer/data/raw/hodoscope_extnim_2026/` or any other raw data folder. Never open `.root` data files in `UPDATE` or `RECREATE` mode; use `READ` mode only.

1. **Host Execution**: Commands must run directly on the host environment (PID 1 is `systemd`, `pgrep -c .` > 100). Never execute within an isolated container or sandbox where CAEN hardware libraries or external data paths are inaccessible.
2. **Resource Throttling**: Wrap all heavy computation (ROOT macros, C++ compilation, Python scripts) with resource limits:
   ```bash
   ( ulimit -v 4194304; timeout 600 <command> )
   ```
   Before launching any multi-run analysis or compilation phase, check `free -m`. If available RAM is less than 4000 MB, STOP and report.
3. **Zero Deletion Policy**: Superseded code, configurations, or data must be relocated to `attic/`, `config_commissioning/`, or `outputs/_superseded/`, never permanently deleted.
4. **Evidence Tags**: Every statement of fact regarding digitizer configuration, channel mapping, data formats, or run history must carry an evidence tag:
   - `[V: <file>:<line>]` or `[V: cmd <command>]` for directly verified facts.
   - `[I: <reason>]` for deduced/inferred facts.
   - `[U]` for unverified facts. If a parameter is unknown, record it explicitly as `not recorded`.
5. **DAQ Binary Protection**: Do NOT rebuild DAQ binaries (`libCaenDigitizerSiLab.so`, `test_acquire_sw`, `test_acquire_ext_nim`) during routine analysis tasks. Do NOT execute acquisition binaries without explicit physical hardware authorization.

---

## 3. Hardware Architecture & Authoritative Channel Mapping

### Hardware Configuration
- **Digitizer Model**: CAEN VX1740 (VME 6U module, 64 analog channels, 12-bit ADC) [V: `test_open.cxx:16`, `acquire_sw.log:20`].
- **Readout Optical Link**: CAEN A4818 USB 3.0 to CONET2 optical interface (ConnectionType 2 / enum 5, LinkNum 62578, ConetNode 0, VMEBase 0x0) [V: `config_vx1740_ext_nim_1000000.conf:1-5`].
- **Sampling Specification**: Sampling frequency 62.5 MS/s (nominal 16 ns sampling interval, code interval `1e9 / 65e6 = 15.38 ns`) [V: `CaenDigitizerSiLab.cxx:454`, `reconstruct_corrected.cxx:49`].
- **Record Length**: 1024 samples per waveform (1023 samples recorded per event) [V: `config_vx1740_ext_nim_1000000.conf:6`].

### Authoritative Detector Channel Mapping

The detector system consists of two hodoscope units: **Hodoscope 2018 (hod2018)** and **Hodoscope 2019 (hod2019)**, each having two orthogonal planes ($X$ and $Y$) of 16 scintillator bars (2.0 mm pitch, 98.4 mm inter-plane separation) [V: `reconstruct_corrected.cxx:57-58`].

#### 1. Hodoscope 2018 (Direct Sequential Mapping)
Channels 0 through 31 are sequentially mapped directly to physical bars [V: `reconstruct_corrected.cxx:197-198`, `reconstruct_tracks_interleaved.cpp:66-70`, `final_hodoscope_analysis.cpp:78-79`]:
- **Plane X (16 bars, index $0..15$)**: Digitizer Channels `0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15` ($k = c$).
- **Plane Y (16 bars, index $0..15$)**: Digitizer Channels `16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31` ($k = c - 16$).

#### 2. Hodoscope 2019 (Interleaved Cabling Mapping)
Channels 32 through 63 follow the physical cabling recorded in the experimental logbook (Excel "Disposición Correcta") [V: `reconstruct_corrected.cxx:62-63`, `reconstruct_tracks_interleaved.cpp:27-39`]:
- **Plane X (16 bars, index $0..15$)**:
  ```cpp
  // reconstruct_corrected.cxx:62
  static const std::array<int,16> kH19X = {53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56};
  ```
- **Plane Y (16 bars, index $0..15$)**:
  ```cpp
  // reconstruct_corrected.cxx:63
  static const std::array<int,16> kH19Y = {37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40};
  ```
- **Physical Position Conversion**:
  $$\text{coord}_{\text{mm}} = (k + 0.5 - 8.0) \times \text{pitch}_{\text{mm}}$$

> [!IMPORTANT]
> Any analysis output generated using sequential channel assumptions for hod2019 ($X0=\text{Ch48}, Y0=\text{Ch32}$) is invalid. All such runs have been moved to `outputs/_superseded/` [V: `outputs/_superseded/README.md`].

---

## 4. Raw Data Format & ROOT Trees

All production runs (v1 through v10) are stored as ROOT files containing a single `TTree` named `data_dig` [V: `CaenDigitizerSiLab.cxx:585-645`]:

### Tree Schema
| Branch Name | ROOT Type | C++ Type | Meaning | Units / Scale |
| :--- | :--- | :--- | :--- | :--- |
| `Ch0` .. `Ch63` | `vector<float>` | `std::vector<float>` | ADC voltage waveform samples | ADC counts ($0..4095$, baseline $\sim 3600$) |
| `time` | `vector<float>` | `std::vector<float>` | Sample time relative to trigger window | Nanoseconds (ns), step $\approx 15.38$ ns |
| `timestamp` | `Long_t` | `Long_t` | Monotonic host readout clock | Microseconds ($\mu$s) since host boot (`CLOCK_MONOTONIC`) |

### AutoSave & AutoFlush Characteristics
- **Buffer Flushing**: In production runs v1–v5 and v7–v10, ROOT `AutoSave` and `AutoFlush` occurred automatically every **1,092 to 1,548 entries** (corresponding to roughly **35 to 45 minutes** of acquisition at $\approx 0.57$ Hz) [V: C12; `reorg_r0_root_metadata.txt`].
- **Short Run v6**: Run v6 was terminated after only 642 entries and **never reached its first flush** (negative default metadata values present) [V: C12].
- **Unclosed Run v10**: Run v10 was interrupted without calling `TTree::Write()` or `TFile::Close()`. When opening in ROOT, a single key recovery occurs automatically at file offset 5936588098, recovering 242,520 valid entries. Any events acquired after the last AutoSave cycle were lost [V: C11].
- **Timestamp Definition**: `timestamp` represents `std::chrono::steady_clock` (Linux `CLOCK_MONOTONIC`) recorded at host readout [V: `CaenDigitizerSiLab.cxx:592-593`]. It is **not Unix epoch**, **not calendar time**, and **not the digitizer internal trigger time** (CAEN `TriggerTimeTag` is not stored). Calibrated wall-clock intervals are documented in `runs/runs_2026.csv`.

---

## 5. Build System & Compilation

### DAQ Build (`<REPO>/Makefile`)
The root `Makefile` builds only the DAQ library and test binaries. Analysis targets have been detached, and the `clean` target will never delete `.root` or `.txt` data files:
```bash
make all      # Builds libCaenDigitizerSiLab.so.1.0.1, example, test_open, test_acquire_sw, test_acquire_ext_nim
make clean    # Cleans only .o, .d, .pcm, and DAQ test executables (safe against data loss)
```

### Analysis Build (`<REPO>/analysis/Makefile`)
All offline analysis binaries in `analysis/src/` are compiled into `analysis/bin/`:
```bash
# Compile all analysis programs
make -C analysis

# Compile a specific analysis program
make -C analysis bin/scan_active_channels
make -C analysis bin/reconstruct_tracks_interleaved
make -C analysis bin/channel_health_analysis

# Clean analysis binaries only
make -C analysis clean
```
Compilation uses `g++ -std=c++17 -O2 -Wall $(root-config --cflags --ldflags --glibs) -lSpectrum`.

---

## 6. Execution of Analysis Pipelines

All analysis shell scripts reside in `analysis/scripts/` and automatically resolve the repository root portably. They read inputs via the root symlinks and write outputs to `outputs/`:

```bash
# 1. Channel health monitoring across runs v1-v9
./analysis/scripts/run_channel_health_analysis.sh

# 2. High-voltage scan comparison across runs v1-v4
./analysis/scripts/run_compare_hv_runs.sh

# 3. Hodoscope 2018 2D performance and diagnostics (run v4)
./analysis/scripts/run_final_hod2018_analysis_v4.sh

# 4. Candidate track reconstruction with verified interleaved mapping (run v7)
./analysis/scripts/run_reconstruct_tracks_interleaved_v7_corrected_mapping.sh

# 5. Single photo-electron (SPE) gain analysis
./analysis/scripts/run_spe_gain_analysis.sh
```

---

## 7. Known Issues & Operational Pitfalls

1. **Hodoscope 2019 Self-Trigger Scan (Zero Counts)**: The threshold scan script `run_hod2019_selftrigger_threshold_scan.sh` produced 0 entries across all thresholds down to $-2$ mV; the underlying hardware/firmware cause was not determined, and original scan logs are absent [V: C10].
2. **Timing Window in Tracking**: Muon coincidence between hod2018 and hod2019 occurs strictly in the time window $[4000, 4400]$ ns ($\Delta t \approx 0$). Early analysis scripts that looked at $t \approx 3510$ ns sampled background noise rather than synchronized cosmic muons [V: `reconstruct_corrected.cxx:3-4`].
3. **Run v10 Crash / Recovery**: Always open `lectura_vx1740_ext_nim_1M_v10.root` in `READ` mode. Opening corrupted/unclosed ROOT files in `UPDATE` mode can invalidate or overwrite the recovered key header [V: C11].
4. **Bias Voltages**: Operational bias voltages ($V_{\text{bias}}$) for production runs v1 through v10 were not recorded in code, configurations, or logs; they remain marked as `not recorded` [V: R0 audit].

