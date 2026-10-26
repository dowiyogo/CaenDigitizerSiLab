#!/bin/bash
# run_macro.sh — XPLUS-4b: corre una macro de analysis/night_2026-10-25/ compilada (ACLiC) en primer plano, nice 10, ionice -c3,
# ulimit -v 4 GB, timeout 3600. Uso: run_macro.sh '<macro>.C+("args")'
cd /home/xplus/night_2026-10-25/build
cp -u /home/xplus/digitizer/CaenDigitizerSiLab-jai/analysis/night_2026-10-25/*.C .
export PATH=$HOME/root6/bin:$PATH LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH
ulimit -v 4194304
exec nice -n 10 ionice -c3 timeout 3600 root -l -b -q -e 'gSystem->AddIncludePath("-I/home/xplus/digitizer/CaenDigitizerSiLab-jai/daq/core/include")' -e 'int XPLUS4_tag = 4;' "$1" 2>&1 | grep -v -E '^Info in <TUnixSystem::ACLiC>|^$'
