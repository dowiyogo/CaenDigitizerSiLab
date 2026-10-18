# CLAUDE.md

Guidance for Claude Code sessions in this repository. Read `AGENTS.md` as well: it holds
the detailed record of the old DAQ, runs v1-v10 and the analysis. Where the two disagree,
the rules in this file are the stricter ones and apply; the factual differences found on
2026-10-18 are listed at the end.

## Project

Acquisition and analysis for two scintillator hodoscopes read by MPPCs:

- hod2018 on digitizer channels 0-31, hod2019 on channels 32-63. Each station has an X and
  a Y plane of 16 bars (75 x 3 x 1 mm), one MPPC per bar.
- Trigger: coincidence of two external scintillator plates (N845 at -100 mV, N405 AND) into
  TRG IN of the digitizer, about 0.56-0.62 Hz.
- Digitizer: CAEN VX1740 (64 ch, 8 groups of 8, 12 bit, 2 Vpp, 16 ns/sample), reached
  through a CAEN A4818 USB3-to-CONET2 adapter (`CAEN_DGTZ_USB_A4818`, linkNum = USB PID 62578).
- HV: CAEN SY5527LC mainframe at 192.168.1.118 (private network on `enp3s0`), board A7030DP
  in slot 3. CH00-07 bias HV groups 1-8 (8 MPPCs each, groups mix both stations);
  CH08-11 unused. The mainframe clock is about 81 min ahead: never use its timestamps.
- Work in progress: a new DAQ with HV slow control (`hvmon` daemon + `daqrun` per run).
  Design: `docs/daq_v2_design.md` (Spanish). Implementation goes under `daq/`.
- Inputs: `docs/handoff_analysis_2026-10-18.md` (analysis of v1-v10 and requirements; the
  project-side version is the reference). Items for the analysis chat go to
  `docs/notes_for_analysis.md`; do not fix analysis code from DAQ tasks.
- October 2026 campaign: channel mapping `kCfg_2026_10`; HV group g = SY5527 CH0(g-1).
  `outputs/pe_calibration_v8.txt` is from May and must not be used for October: until the
  commissioning calibration exists the DAQ runs in the "no calibration" state.

## Paths

    repository            /home/xplus/digitizer/CaenDigitizerSiLab-jai
    raw data (read-only)  /data/daq/raw/<campaign>/   (also /home/xplus/digitizer/data/raw -> /data/daq/raw)
    slow control (new)    /data/daq/slowcontrol/      (does not exist yet)
    raw manifests         /data/daq/manifests/
    DAQ v2 user config    ~/.config/daqv2/ (secrets.env mode 600, hvmon.json, hvctl.json)
    DAQ v2 runtime        $XDG_RUNTIME_DIR/daqv2/ (hvmon.lock, hvmon.sock), /dev/shm/daqv2/
    analysis outputs      outputs/  (not in git)
    ROOT                  ~/root6 (6.32.10, built with C++17)
    CAEN HV Wrapper 6.6   /usr/lib64/libcaenhvwrapper.so, /usr/include/CAENHVWrapper.h (needs -DUNIX)
    CAEN digitizer libs   /usr/lib: CAENDigitizer 2.17.3, CAENComm 1.7.0, CAENVMELib 4.0.2
                          (identical unversioned copies in /usr/local/lib; daq/ links the
                          versioned /usr/lib files by full path)

## Hard rules

1. Never acquire data and never touch the digitizer: no CAENDigitizer calls, no running
   test_acquire_*, test_open, test_caencomm_open, example or
   run_hod2019_selftrigger_threshold_scan.sh.
2. HV power supply: reading parameters is allowed only in tasks that say so. Any write
   (`CAENHV_SetChParam`, `CAENHV_SetBdParam`, any `Set*`, `ExecComm`) requires explicit
   authorization written in the task prompt.
3. `/data/daq/raw/` is read-only: never write, move, chmod or delete there. Open ROOT
   files there in READ mode only. The only exception is a new campaign directory created
   by the new DAQ itself in an authorized acquisition task.
4. Never delete files: move them to `attic/` (code) or `outputs/_superseded/` (results)
   with one line in that directory's README. Never `git push`, never rewrite history.
5. Before anything heavy, check that `free -m` shows more than 4000 MB available, and run
   it as `( ulimit -v 4194304; timeout 3600 <command> )`. Do not read whole data files into
   memory (a python3 process was OOM-killed at 15.4 GB on this 15 GB machine).
   Exception (accepted in DAQ2-3): ASan/UBSan/TSan builds and test runs reserve terabytes of
   virtual memory and cannot run under `ulimit -v`; run them as
   `timeout 3600 systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 <command>`, and
   TSan binaries under `setarch x86_64 -R` (the test CMake does this automatically). Every
   other command keeps `ulimit -v 4194304`.
10. Credentials: never read shell history or other personal files to obtain passwords, and
   never create, copy, move or print credential files. Secrets live only in
   `~/.config/daqv2/secrets.env` (created by René, mode 600); programs read them themselves.
   You may check its mode and variable names (`stat`, `grep -o '^[A-Z_]*'`), never its values.
   If secrets are missing, stop and ask René. Secrets never go into argv, logs, state files,
   HTTP pages, test output or commits.
6. CAENGECO2020 (HV GUI, a snap) may be running. Do not close or touch it.
7. Only one task runs on this host at a time; the analysis chat also sends tasks. Do not
   start background jobs that outlive the task. At the start of every task run `git status`,
   `pgrep -af -i 'caen|wavedump|digitizer|daq|acq|test_acquire'` (stop if an acquisition
   runs) and check the modification times of the files you will edit. If a file you are
   about to edit was changed in the last hour by someone other than you and the task prompt
   does not mention it, stop and ask. Files the task prompt names explicitly are
   authorized. This replaces the AGENTS.md rule "do not write if any file changed in the
   last hour".
8. Copy every reported number verbatim from command output. If something was not
   measured, say so. Do not run `make` in the repository root (DAQ binaries used for data
   taking live there); use `make -n`.
9. Do not rename keys of the old `.conf` files or change their units.

## Build and test

ROOT environment: `source ~/root6/bin/thisroot.sh` (from `~/.bashrc`). In the Claude Code
shell this script may fail outside the repository directory; then use
`export PATH=$HOME/root6/bin:$PATH LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH`.

Analysis (Makefile, C++17, run everything from the repository root):

    make -C analysis                      # binaries in analysis/bin/
    analysis/scripts/run_*.sh             # pipeline scripts, cd to the root themselves

New DAQ (CMake, `daq/`; phase 0 exists: core, read-only HV library, `hvctl`). Tests use
GoogleTest from apt (`libgtest-dev` 1.14) through `gtest_discover_tests`. If GoogleTest is
missing, stop and say so; do not vendor a test framework. Details for René: `daq/README.md`.

    # normal build (from the repository root)
    cmake -S daq -B daq/build -DCMAKE_BUILD_TYPE=RelWithDebInfo
    ( ulimit -v 4194304; timeout 3600 cmake --build daq/build -j4 )
    ( ulimit -v 4194304; timeout 3600 ctest --test-dir daq/build --output-on-failure -j4 )

    # sanitizers: ASan+UBSan together, TSan alone. They reserve terabytes of virtual memory and
    # do not run under ulimit -v; cap real memory with a cgroup instead. TSan tests run under
    # setarch -R automatically (GCC 13 TSan vs this kernel's ASLR).
    cmake -S daq -B daq/build-address-undefined -DDAQ_SANITIZE=address,undefined
    cmake -S daq -B daq/build-thread -DDAQ_SANITIZE=thread
    timeout 3600 systemd-run --user --scope --quiet -p MemoryMax=4G -p MemorySwapMax=0 \
        cmake --build daq/build-thread -j4
    timeout 3600 systemd-run --user --scope --quiet -p MemoryMax=4G -p MemorySwapMax=0 \
        ctest --test-dir daq/build-thread --output-on-failure -j4

ROOT must be findable at configure time (`ROOTSYS` or `~/root6`); no phase 0 target links it.
The CAEN HV Wrapper is the IMPORTED target `caen::hvwrapper` (adds `-DUNIX`).

Tests must run against the fakes only (`FakeHVBackend`, fake digitizer, fake temperature
source). No test may open the digitizer or connect to the SY5527.

HV slow control (connects to the mainframe: only when the task authorizes HV reads). `hvmon`
is the only owner of the Wrapper session (lock `$XDG_RUNTIME_DIR/daqv2/hvmon.lock`); while it
runs, `hvctl` talks to it through `$XDG_RUNTIME_DIR/daqv2/hvmon.sock` and never opens a session.

    export LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH          # hvmon links ROOT
    ( ulimit -v 4194304; daq/build/apps/hvmon/hvmon --config daq/config/hvmon.example.json )
    daq/build/apps/hvmon/hvmon --fake default --config <cfg with scratch slowcontrol_dir>  # no hardware
    daq/build/apps/hvctl/hvctl status | read | changes --n 20 | note "<text>"
    daq/build/apps/hvctl/hvctl measure-offsets --seconds 60          # only channels with Pw = 0
    daq/build/apps/hvctl/hvctl pw --ch N on|off --authorize "<reason>"  # needs allow_commands=true
    # without hvmon (direct, read-only):
    daq/build/apps/hvctl/hvctl map | params --slot 3 | read --slot 3 | bench --slot 3 --n 100
    daq/build/apps/hvctl/hvctl --fake default read --slot 3

Output: `/data/daq/slowcontrol/YYYY/` (daily CSV/JSONL with `.schema.json`), state
`/dev/shm/daqv2/hvmon_state.json`, read-only page `http://127.0.0.1:8081/` (remote: SSH tunnel
over Tailscale). Unit file `daq/systemd/hvmon.service` (not installed; René enables it).
Writes: the only write in the code base is `HVBackend::setChannelPw()`; `hvmon` refuses it unless
its config has `allow_commands: true`, which needs explicit authorization in the task prompt.
Secrets: `DAQ_HV_HOST`, `DAQ_HV_USER`, `DAQ_HV_PASSWORD` (and `DAQ_SMTP_*`, `DAQ_NOTIFY_*`) from
the environment or `~/.config/daqv2/secrets.env`; template `daq/config/secrets.env.example`.
`~/.config/daqv2/hv_credentials.json` is a deprecated fallback (warning); René removes it.

## Coding conventions for daq/

- Identifiers in English, comments in Spanish. Reports and design docs in Spanish.
- C++17 everywhere in daq/ (ROOT 6.32.10 was built with C++17; C++20 compiles but ROOT
  warns). Use the small in-house `daq::JThread`/`StopToken`, not `std::jthread`. Standard
  library concurrency only.
- SI units internally (V, A, s). Convert at the hardware boundary (the A7030 reports
  current in uA). Put units in branch titles, CSV headers and JSON keys (`_V`, `_A`, `_s`).
- Times: UTC nanoseconds since the Unix epoch (`int64`), from `CLOCK_REALTIME`; intervals
  from `CLOCK_MONOTONIC`. Never local time in data files. The host is in America/Santiago:
  UTC-3 in October, UTC-4 in May-June (runs v1-v10 and the local times in
  `runs/good_intervals_2026.csv` are UTC-4).
- HV writes go only through `hvmon` (auto-retry, HV-group test, `hvctl` commands), each with
  the authorization text logged. HTTP pages are read-only and bind to 127.0.0.1 (or the
  Tailscale address), never to the university interface.
- Every ROOT object has exactly one owning thread. Every hardware dependency sits behind an
  interface with a fake implementation; time comes from an injectable clock.

## Known differences between AGENTS.md and the system (2026-10-18)

- AGENTS.md gives the raw-data location as /home/xplus/digitizer/data/raw; that path is now
  a symlink to /data/daq/raw.
- AGENTS.md says hod2018 X bar k = channel k; `analysis/src/hodo_mapping.h` (kCfg_2026_05,
  used for v1-v10) has a different, measured hod2018 order.
- AGENTS.md rule 7 uses `timeout 600`; the current rule is `timeout 3600`.
- AGENTS.md says not to write if any file changed in the last hour; rule 7 above replaces it.
- Legacy `daq/tools/legacy/hv_readall` takes the password as an argument (visible in `ps`); do
  not use it when `hvctl` can do the job.
- AGENTS.md places hod2018 pickup at t ≈ 3525 ns; a read-only check on v8 (2026-10-18) did
  not find it (`docs/notes_for_analysis.md`, item 3).
