#!/bin/bash
# run_ntuple.sh — XPLUS-4b etapa 1: lanza ntuple_extract.C sobre un archivo o patrón, desacoplado (setsid nohup), con nice 10, ionice -c3,
# ulimit -v 4 GB y timeout 3600; log en logs/ntuple_<nombre>.log y PID en pids/ntuple_<nombre>.pid. La marca XPLUS4 va en la línea de comando.
# Uso: run_ntuple.sh <nombre> <patrón de archivos> <may|daqv2|sw>
N=/home/xplus/night_2026-10-25; D=/data/daq/derived/night_2026-10-25/ntuples
name=$1; pat=$2; fmt=$3
cd $N/build
export PATH=$HOME/root6/bin:$PATH LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH
setsid nohup bash -c "ulimit -v 4194304; exec nice -n 10 ionice -c3 timeout 3600 root -l -b -q -e 'gSystem->AddIncludePath(\"-I/home/xplus/digitizer/CaenDigitizerSiLab-jai/daq/core/include\")' -e 'int XPLUS4_tag = 4;' 'ntuple_extract.C+(\"$pat\", \"$fmt\", \"$D/$name\", 4)'" > $N/logs/ntuple_$name.log 2>&1 < /dev/null &
sleep 1
pid=$(pgrep -f -- "^timeout 3600 root .*ntuples/$name\"" | head -1)
echo "$pid" > $N/pids/ntuple_$name.pid
echo "$name PID $pid"
