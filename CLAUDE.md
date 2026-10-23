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
    slow control          /data/daq/slowcontrol/YYYY/ (hvmon; closed days gzipped)
    raw manifests         /data/daq/manifests/
    simulated output      /data/daq/sim/ (raw/, manifests/, slowcontrol/, logs/, config/, offline/): every
                          run or hvmon with a fake/replay digitizer, FakeHVBackend or test_mode
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
   The Antigravity editor extension (`agy`, Gemini) may be running: it starts with the editor
   and René keeps it idle on this repository. Its presence is not a reason to stop. Instead,
   before EVERY commit run `git status` and compare it with the files you changed in the task;
   if any other file changed (besides the known ones, e.g. `analysis/scripts/hv_diag6_shape.py`),
   stop and report instead of committing.
8. Copy every reported number verbatim from command output. If something was not
   measured, say so. Do not run `make` in the repository root (DAQ binaries used for data
   taking live there); use `make -n`.
9. Do not rename keys of the old `.conf` files or change their units.
11. Unstable network (René's connection and the internet drop in the afternoon, Chile time):
   long hardware operations (acquisition runs, validation runs, anything over a few minutes)
   are started detached under systemd, never tied to an interactive terminal:
   `systemd-run --user --unit=<name> --property=RuntimeMaxSec=<s> --collect ...` (or the
   `hvmon` user service), with a planned stop time and logs to files under `/data/daq/...` or
   the journal; their status is checked afterwards (`systemctl --user status <name>`,
   `journalctl --user -u <name>`). Nothing depends on a terminal staying open.
12. `hvmon` runs as the user service `hvmon.service` (since DAQ2-6, binaries in
   `~/.local/opt/daqv2/bin/`, config `~/.config/daqv2/hvmon.json`, `allow_commands: false`).
   Tasks use it through `hvctl` (socket) and must never start another `hvmon` against the
   real mainframe; `systemctl --user stop|restart hvmon` only when the task says so. Rebuilding
   `daq/build` does not touch the running binary; `daq/systemd/install.sh` installs new
   binaries (then `systemctl --user restart hvmon`, a gap of ~2.4 s in the slow control).

## Build and test

ROOT environment: `source ~/root6/bin/thisroot.sh` (from `~/.bashrc`). In the Claude Code
shell this script may fail outside the repository directory; then use
`export PATH=$HOME/root6/bin:$PATH LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH`.

Analysis (Makefile, C++17, run everything from the repository root):

    make -C analysis                      # binaries in analysis/bin/
    analysis/scripts/run_*.sh             # pipeline scripts, cd to the root themselves

New DAQ (CMake, `daq/`; phases 0-5 exist: core, HV library, `hvctl`, `hvmon` with auto-retry,
notifier and compression, `sctool`, `daqrun` and `daqcampaign` with fake, replay or the real
VX1740 (`digitizer.backend = "caen"`, phase 4), `digprobe`; phase 5: data-driven flags 6/7/8/11,
IMon warning, in-time panels, hourly gain, per-run TTT clock, `daqoffline`).
Tests (`daq_tests`, and `daq_root_tests` for the parts that need ROOT) use
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

ROOT must be findable at configure time (`ROOTSYS` or `~/root6`); it is linked by `hvmon`,
`daqrun`, `daq_io`, `daq_dig_replay` and `libdaq_screader.so`, not by the other libraries.
The CAEN HV Wrapper is the IMPORTED target `caen::hvwrapper` (adds `-DUNIX`).

Tests must run against the fakes only (`FakeHVBackend`, fake digitizer, fake temperature
source). No test may open the digitizer or connect to the SY5527.

HV slow control (connects to the mainframe: only when the task authorizes HV reads). `hvmon`
is the only owner of the Wrapper session (lock `$XDG_RUNTIME_DIR/daqv2/hvmon.lock`); while it
runs, `hvctl` talks to it through `$XDG_RUNTIME_DIR/daqv2/hvmon.sock` and never opens a session.

    systemctl --user status hvmon ; journalctl --user -u hvmon -n 50   # the real hvmon (service)
    daq/systemd/install.sh && systemctl --user restart hvmon            # deploy new binaries
    export LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH          # hvmon links ROOT
    daq/build/apps/hvmon/hvmon --fake default --config <cfg without output keys>  # no hardware, /data/daq/sim
    daq/build/apps/hvctl/hvctl status | read | changes --n 20 | note "<text>"
    daq/build/apps/hvctl/hvctl measure-offsets --seconds 60          # only channels with Pw = 0
    daq/build/apps/hvctl/hvctl pw --ch N on|off --authorize "<reason>"  # needs allow_commands=true (true since DAQ2-10)
    daq/build/apps/hvctl/hvctl notify-test [--ch N]   # sends ONE real e-mail: only if the task authorizes it
    daq/build/apps/sctool/sctool compress [--keep] FILE... | cat FILE | digest FILE   # plain or .gz
    # without hvmon (direct, read-only):
    daq/build/apps/hvctl/hvctl map | params --slot 3 | read --slot 3 | bench --slot 3 --n 100
    daq/build/apps/hvctl/hvctl --fake default read --slot 3

Output: `/data/daq/slowcontrol/YYYY/` (daily CSV/JSONL with `.schema.json`), state
`/dev/shm/daqv2/hvmon_state.json`, read-only page `http://127.0.0.1:8081/` (remote: SSH tunnel
over Tailscale). Unit file `daq/systemd/hvmon.service` (not installed; René enables it).
Since DAQ2-10 `allow_commands` is true in `~/.config/daqv2/hvmon.json` so René can switch channels
from the console; every write still needs `--authorize`, and a task may write only if its prompt
authorizes it (rule 2). Pw off is verified by the Status `on` bit, not by VMon (the detector PCB
discharges with tau ~110-150 s). Current mapping `hodo_mapping_2026_10_r2.json` (groups measured),
October calibration `pe_calibration_2026_10_c1.json` (design 15.8).
Writes: the only writes in the code base are `HVBackend::setChannelPw()`, `clearBoardAlarm()`
(`ClrAlarm`) and, since DAQ2-14, `setChannelV0Set()`; `hvmon` refuses all of them unless its config has
`allow_commands: true`, which needs explicit authorization in the task prompt. V0Set is a bounded write
available only inside a `hvmon` scan session (`daq/hvmon/include/daq/hvmon/bias_scan.h`;
`hvctl bias begin|set|restore|status`): the session starts only if the V0Set read back matches the
expected nominal; each write is <= nominal (default), >= nominal − 2.0 V, a step <= 0.5 V, on a channel that is on
and stable; V0Set above nominal is allowed only inside a session opened with `--max-above` (hard cap nominal + 1.0 V
in code, steps <= 0.25 V above nominal) and only when the task explicitly authorizes it; the ramp is verified (V0Set read back, VMon − V0Set within 0.10 V of the offset at the previous
verified step, VMon moved with the step by at least half of it; any Status bit other than on/ramp aborts);
`bias restore` returns every channel to nominal in <= 0.5 V steps and Pw = 1; an active session is
persisted (`~/.local/state/daqv2/hvmon_bias_state.json`) and restored automatically when `hvmon` starts.
Like any HV write, using it needs explicit authorization in the task prompt (rule 2). Auto-retry (`auto_retry.enabled_channels`, design 8)
requires `allow_commands: true`; `test_mode` with CH00-07 also needs
`test_mode_allow_detector_channels`. Retry state persists in
`$XDG_STATE_HOME/daqv2/hvmon_retry_state.json` (default `~/.local/state/daqv2/`). Notifier
(`notifier.enabled`, design 8.4) sends real e-mail through `curl --config -`; keep it disabled
unless the task authorizes e-mail and count the messages. Closed days in `/data/daq/slowcontrol`
are gzipped by `hvmon` after verification (design 8.5); readers must accept `.gz`. A test
`hvmon` beside a running one needs its own `runtime_dir` (socket path ≤ 107 bytes, e.g.
`/run/user/1000/daqv2-test`), `state_dir`, `slowcontrol_dir`, `state_file` and HTTP port.
`hvmon --fake` (or `test_mode`) refuses `slowcontrol_dir` inside `/data/daq/slowcontrol` and, for
fakes, the real state dir; with the keys absent it uses `/data/daq/sim/slowcontrol`,
`/dev/shm/daqv2-sim`, `$XDG_RUNTIME_DIR/daqv2-sim`. E-mails from a fake or test_mode start with
`[SIMULADO]`, notify-test with `[PRUEBA]`.
Secrets: `DAQ_HV_HOST`, `DAQ_HV_USER`, `DAQ_HV_PASSWORD` (and `DAQ_SMTP_*`, `DAQ_NOTIFY_*`) from
the environment or `~/.config/daqv2/secrets.env`; template `daq/config/secrets.env.example`.
`~/.config/daqv2/hv_credentials.json` is a deprecated fallback (warning); René removes it.

DAQ runs. A simulated run (`fake`/`replay`) refuses `raw_root`/`manifests_dir` inside
`/data/daq/raw`/`/data/daq/manifests` and defaults to `/data/daq/sim/...`. A run needs an open
campaign. Real runs (`caen`, VX1740 via A4818) only in tasks that authorize the digitizer (rule 1);
check first that nothing else holds the A4818 (`fuser /dev/bus/usb/002/002`, no old DAQ process).
Real configs live in `/data/daq/config/<campaign>/`, logs in `/data/daq/logs/`; runs are started
detached from the installed binaries (`~/.local/opt/daqv2/bin/daqrun`, rule 11). Board facts
(design 15.5): record 625 -> 624, 0x8114 read back, memory full = bit 4 of 0x8104 (1024 events;
lost triggers do not advance the counter), TTT 8 ns / 31 bit (bit 31 = rollover flag), the
external trigger depends on the NIM8302 crate (trips on temperature) and its cables.

    daq/build/apps/daqcampaign/daqcampaign --sim open|close|status <name>   # real: --authorize "<reason>"
    ~/.local/opt/daqv2/bin/digprobe info | memfull --config <cfg>   # board info/registers (needs authorization)
    daq/build/apps/daqrun/daqrun --config <cfg.json> [--campaign N] [--max-duration S] [--no-http]
    # example config daq/config/daqrun.example.json, synthetic scenario daq/config/fake/
    # page http://127.0.0.1:8080/ ; SIGINT/SIGTERM closes the segment and the index cleanly
    root -l -b -q 'analysis/macros/plot_hv_imon_day.C("20261019")'   # slow-control reader example
    # phase 5: monitor + flags + panels + TTT clock over written runs (read-only), or the fake with a
    # manual clock; output only under /data/daq/sim/offline/ (design 15.7)
    daq/build/apps/daqoffline/daqoffline --config <cfg.json> --out /data/daq/sim/offline/<name> RUN_seg*.root
    daq/build/apps/daqoffline/daqoffline --config <cfg.json> --out /data/daq/sim/offline/<name> --fake <scenario.json> --duration S

Run files: `<raw_root>/<campaign>/<campaign>_<UTC>_runNNN_segMMM.root` (read-only after close),
`..._runNNN.index.json`, `..._runNNN.config.json`, `SHA256SUMS`, `CAMPAIGN_OPEN|CAMPAIGN_CLOSED`;
design 4.1-4.4 and 15.4. Shared with the analysis (header-only, no ROOT):
`daq/core/include/daq/pulse_finder.h`, `daq/core/include/daq/trigger_interval_stats.h`,
`daq/core/include/daq/ttt.h`, `daq/core/include/daq/gain_fit.h`. Data-driven flags and e-mail kinds:
`data_flags` config key (design 15.7); a task can restrict daqrun e-mails with `notifier.only_kinds`
and `notifier.max_messages`.

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
