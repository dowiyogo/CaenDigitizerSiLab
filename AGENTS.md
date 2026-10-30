# AGENTS.md

This file is for agents and people working in this repository. Read it before running
anything. Statements carry [V: source] when verified, [I: reason] when inferred, and
values never recorded are written as "not recorded".

## Purpose

The repository contains the acquisition program for a CAEN VX1740 digitizer and the
offline analysis of the data it records. The digitizer reads the SiPM signals of two
scintillator hodoscopes, hod2018 on channels 0-31 and hod2019 on channels 32-63
[V: config_vx1740_autotrigger_hod2019_ch32_63.conf:2-3]. Each hodoscope has an X plane
and a Y plane of 16 bars.

This computer also controls the digitizer, so an acquisition may be running while you
work here.

## Before writing anything

Check that no acquisition is running and that no file changed in the last hour:

    pgrep -af -i 'caen|wavedump|digitizer|daq|acq|test_acquire'
    find . -type f -mmin -60 -not -path './.git/*'

If either command shows activity, do not write in the repository.

## Rules

1. Do not run test_acquire_ext_nim, test_acquire_sw, test_open, test_caencomm_open,
   example or run_hod2019_selftrigger_threshold_scan.sh, and do not open the digitizer,
   unless the person operating the detector asks for it in the current session.
2. Do not create, modify, rename or delete anything under /home/xplus/digitizer/data/raw/.
   Open data files in READ mode only.
3. Do not run make in the repository root unless asked; use `make -n` to inspect. The
   DAQ binaries in the root are the ones used for data taking.
4. Do not rename configuration keys or change their units.
5. Do not delete files. Superseded material goes to attic/ or outputs/_superseded/, with
   one line in that directory's README stating why it was superseded and the evidence.
6. Do not push and do not rewrite git history.
7. Before running analysis programs check `free -m` and stop if available is below
   4000 MB. Run long ROOT or Python jobs as `( ulimit -v 4194304; timeout 600 <command> )`
   and do not read whole data files into memory. On 2026-10-16 the kernel killed a
   python3 process at 15.4 GB resident memory on this 15 GB machine [V: dmesg].

## Directory map

    CaenDigitizerSiLab.cxx/.h DAQ driver class
    configuration.C/.h       configuration file parser
    test_acquire_ext_nim.cxx acquisition binary for external trigger
    test_acquire_sw.cxx      acquisition binary for software/self trigger
    test_open.cxx            link verification binary
    Makefile                 builds DAQ binaries only; clean does not touch data
    config*.conf             config.conf (2025 DT5730 default), configs used as code
                             defaults, and production configs (config_vx1740_ext_nim_1000000.conf for v1-v10)
    config_commissioning/    test configurations
    analysis/src/            C++ analysis sources (compiled to analysis/bin/)
    analysis/scripts/        analysis pipeline scripts (run from repository root)
    analysis/macros/         ROOT macros
    analysis/python/         Python analysis scripts
    outputs/                 analysis outputs (not in git; results stay on disk)
    outputs/_superseded/     outputs produced with wrong channel mapping
    runs/runs_2026.csv       run catalog
    logs/                    commissioning logs (not in git)
    docs/                    library manual and reorganization record
    attic/                   legacy DAQ code and pre-reorganization binaries

Raw data are never stored inside the repository.

## Environment

Linux Mint 22 (Ubuntu 24.04 base), GCC 13.3.0, ROOT 6.32.10 loaded by ~/.bashrc:118
(`source ~/root6/bin/thisroot.sh`) [V]. CAEN libraries in /usr/lib: CAENDigitizer 2.17.3,
CAENComm 1.7.0, CAENVMELib 4.0.2 [V: headers]. The VX1740 is reached through a CAEN A4818
USB3-to-CONET2 adapter, USB ID 21e1:0016 [V: lsusb].

## Acquisition program

The program is built from the repository root. This command was verified on 2026-10-17
in a copy of the sources (exit code 0):

    source ~/root6/bin/thisroot.sh
    make libCaenDigitizerSiLab.so.1.0.1 link test_open test_acquire_sw test_acquire_ext_nim

The executables link CaenDigitizerSiLab.o and the dictionary object directly, so
libCaenDigitizerSiLab.so is not loaded at run time [V: readelf -d]. The library file in
the repository is the 2025 committed build; a newer local build was overwritten on
2026-10-17. The copy in /usr/local/lib dates from May 2025 and differs from both.

Runs v1-v10 were started with [V: ~/.bash_history:1338-1353]:

    ./test_acquire_ext_nim config_vx1740_ext_nim_1000000.conf lectura_vx1740_ext_nim_1M_vN.root

The second argument is the output path, relative to the current directory. New runs go
to a directory under /home/xplus/digitizer/data/raw/ (see Data).

### Configuration keys

configuration.C reads the keys below [V: configuration.C:67-96]. The last column is the
value in config_vx1740_ext_nim_1000000.conf, used for v1-v10.

| key | meaning | unit or values | v1-v10 |
|---|---|---|---|
| model | digitizer model | 1740 = VX1740 | 1740 |
| connectionType | link type | 0 USB direct (DT57xx), 1 A2818/A3818 optical, 2 A4818 | 2 |
| linkNum | link identifier; USB PID for the A4818 | integer | 62578 |
| conetNode | position in the CONET chain | integer | 0 |
| vmeBaseAddress | VME base address | hex | 0x00000000 |
| polarization | centring of the 2 Vpp input range | 0, 1, -1 | 0 |
| vthreshold | channel self-trigger threshold; not used with external trigger | mV | -100 |
| acqsamples | samples per record | samples | 625 |
| nevents | events per run | events | 1000000 |
| bunches | bunches per measurement (config comment) | integer | 1 |
| ptriggersize | post-trigger size | % of the record | 60 |
| triggerpolaritymode | channel trigger edge | 0 falling, 1 rising | 0 |
| timeout | maximum run duration | s | 1209600 (14 d) |
| triggerSource | trigger source | 0-7 channel, 8 OR of x740 groups, 9 external TRG-IN | 9 |
| triggerGroupMask | groups enabled for x740 self-trigger | hex | not set (default 0xFF) |

The configuration is not stored in the output file. Record it in runs/runs_2026.csv and
copy it next to the data.

### Output file

Each trigger produces one entry of the TTree data_dig [V: CaenDigitizerSiLab.cxx:181-205]:

| branch | type | content |
|---|---|---|
| Ch0 ... Ch63 | vector<float> | waveform in mV, converted by adc2mV [V: CaenDigitizerSiLab.h:103, CaenDigitizerSiLab.cxx:636] |
| time | vector<float> | sample time in ns, step 16.0 ns [V: CaenDigitizerSiLab.h:163, CaenDigitizerSiLab.cxx:642, B2 root check] |
| event | int | sequential event index |
| timestamp | int64 | µs since host boot (std::chrono::steady_clock), read once per readout |

Records configured with acqsamples 625 contain 624 samples per channel, 0 to 9968 ns [V: ROOT check on v6, entry 0].

timestamp is neither Unix time nor the trigger time. The digitizer TriggerTimeTag is
printed for the first three events only and is not stored
[V: CaenDigitizerSiLab.cxx:592-593, 608, 645]. For v1-v10, taken during a single host
boot, calendar time = boot_ref + timestamp with boot_ref = 2026-05-15T19:03:28.8-04:00,
uncertainty about 1 min; the method is in docs/REORG_2026-10-17.md.

The tree is written with TTree::Write when the run ends. During v1-v5 and v7-v10, ROOT
AutoSave stored it every 1092-1548 entries, about 35-45 min at 0.57 Hz [V: GetAutoSave].
A run that stops abnormally loses the entries after the last AutoSave; this happened to v10.

## Data

Raw data live in /home/xplus/digitizer/data/raw/, one directory per campaign, each with a
SHA256SUMS file. Files and directories are read-only.

    hodoscope_extnim_2026/            v1-v10, external NIM trigger, 16.03 GB
    commissioning_2026-05-22/         test runs of 2026-05-22 and 2026-05-26
    hod2019_threshold_scan_2026-05-26/ self-trigger scan (0 entries) and its configs

To store a new series of runs, create a new campaign directory, write the runs there,
copy the config file used, then:

    cd /home/xplus/digitizer/data/raw/<campaign>
    sha256sum * > SHA256SUMS
    chmod a-w * .

and add one row per run to runs/runs_2026.csv, including the command line and the bias
voltage.

The ten lectura_vx1740_ext_nim_1M*.root files in the repository root are symlinks into
hodoscope_extnim_2026/. They exist because these sources open the files by bare name
relative to the repository root: analysis/src/infer_mapping.cxx:61,
analysis/src/reconstruct_cosmic_tracks.cxx:210-211,
analysis/src/reconstruct_corrected.cxx:149, analysis/src/tagprobe_analysis.cxx:324.

## Analysis

Build all programs with `make -C analysis`; binaries go to analysis/bin/. Run everything
from the repository root, because programs write to relative paths under outputs/. The
scripts in analysis/scripts/ change to the root themselves:

    analysis/scripts/run_channel_health_analysis.sh       v1-v9
    analysis/scripts/run_channel_health_refit.sh          v1-v9
    analysis/scripts/run_compare_hv_runs.sh               v1-v4
    analysis/scripts/run_final_hod2018_analysis_v4.sh     v4
    analysis/scripts/run_reconstruct_tracks_interleaved_v7_corrected_mapping.sh   v7
    analysis/scripts/run_spe_gain_analysis.sh             v1-v9

### Channel mapping

hod2018 uses consecutive channels: X bar k = channel 0-15 with k = channel, Y bar
k = channel - 16 for channels 16-31 [V: reconstruct_corrected.cxx:197-198]. hod2019
follows the cabling record "Disposición Correcta" [V: reconstruct_corrected.cxx:62-63]:

    X bars 0-15: 53 61 54 62 55 63 52 60 51 59 50 58 49 57 48 56
    Y bars 0-15: 37 45 38 46 39 47 36 44 35 43 34 42 33 41 32 40

Results made with the sequential hod2019 assumption (X0 = Ch48, Y0 = Ch32) are wrong and
are kept in outputs/_superseded/.

The track scripts pass active_size_mm = 33.0 and z_sep_mm = 98.4 to
reconstruct_tracks_interleaved (where z_sep_mm is the vertical distance along Z between
hod2018 at z=0 and hod2019 at z=z_sep_mm [V: analysis/src/reconstruct_tracks_interleaved.cpp:265,269,358]).
These are code parameters, not measurements recorded in this repository.

Muon signals of both hodoscopes appear at t ≈ 4175 ns in the time branch, window
[4000, 4400] ns; hod2018 shows pickup at t ≈ 3525 ns
[V: outputs/timesync_test_20260531_174521/summary_timesync.txt]. Analyses that used
t ≈ 3510 ns selected that pickup, not muons.

## Known issues

1. The hod2019 self-trigger scan of 2026-05-26 recorded 0 entries at every threshold down
   to -2 mV. The cause was not determined and the scan logs are absent.
2. The acquisition loop polls ReadData without waiting and prints the event counter on
   every iteration; ext_nim.log reached 7.9 MB for a 1000-event run [V: logs/].
3. The configuration and the bias voltage are not stored with the data.
4. Bias voltages for v1-v10: not recorded.
5. outputs/ is not in git; analysis results exist only on this disk.
6. Hardcoded sampling period (15.38 ns vs 16.0 ns) and unit labels (ADC vs mV): fixed; see docs/ANALYSIS_CHANGES_2026-10-17.md.
7. Thresholds labeled as ADC counts applied directly to mV waveforms: fixed; see docs/ANALYSIS_CHANGES_2026-10-17.md.
8. Swapped columns 18X-v8 and 19X-v7 in charge_analysis before commit ccc1945: fixed; see docs/ANALYSIS_CHANGES_2026-10-17.md.
9. Noise-derived muon window (3008 ns) in multirun_analysis before commit ccc1945: fixed; see docs/ANALYSIS_CHANGES_2026-10-17.md.
10. Q integration windows and pedestals (hod2018 ped_loc [-3,+12], hod2019 ped_glob [-3,+22]): see docs/ANALYSIS_CHANGES_2026-10-17.md.
11. Baseline shift and Q noise measurements (mean, RMS, ped_loc vs ped_glob): see docs/ANALYSIS_CHANGES_2026-10-17.md.

## Open questions

- What v<N> in the run names means. compare_hv_runs assumes v1-v4 differ in bias voltage [I].
- What changed in the trigger setup between v1 (1.39 Hz), v2-v4 (0.27-0.44 Hz) and
  v5-v10 (0.56-0.64 Hz) [V: runs/runs_2026.csv].
- Whether a copy of the raw data exists outside this machine.
- Branch reorg/2026-10-17 is not pushed and not merged into jairo_branch.

## hod2018 construction and orientation (René, 2026-10-30)

Vikuiti held with Kapton, no adhesive (air gap); each bar end glued with optical epoxy to a PMMA light guide painted on the
outside, glued to the MPPC. Origin (0,0) at the lower left corner; X bars parallel to the Y axis, Y bars parallel to the X axis.
Upper Y layer read out on the right, lower Y on the left; upper X read out at the bottom, lower X at the top. Details, and the
status of the absolute sense of the strip order (still open): `docs/hod2018_construction_2026-10-30.md`.
