#!/bin/bash
# Uso: analysis/proc_hvdiag5_2026_10_31/run_full.sh (desde cualquier lugar). El paquete de salida debe existir con work/ preparado: *.bin enlazados
# a outputs/proceeding_2026_10/work/ y el resto copiado (ver el README del paquete). outputs/INDEX.md: guardar una copia antes de stage8.
# Pipeline PROC (código de cálculo = 055fe38; gráficos PROC-4) sobre outputs/proceeding_hvdiag5/2026-10-31_full con el good_intervals nuevo.
cd /home/xplus/digitizer/CaenDigitizerSiLab-jai || exit 1
P=outputs/proceeding_hvdiag5/2026-10-31_full
L=$P/logs; mkdir -p $L
export PATH=$HOME/root6/bin:$PATH LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH
run() { local name=$1; shift; local t0=$(date +%s); ( ulimit -v 4194304; nice -n 19 timeout 3600 "$@" ) > $L/$name.log 2>&1; local rc=$?; echo "$(date -u +%FT%TZ) $name rc=$rc $(( $(date +%s) - t0 )) s" | tee -a $L/run.log; return $rc; }
T=daq/tools/proceeding
run stage1 python3 -I -B $T/stage1.py $P
run stage2 python3 -I -B $T/stage2.py $P
run stage3 python3 -I -B $T/stage3.py $P
run stage4 python3 -I -B $T/stage4.py $P
run stage5 python3 -I -B $T/stage5.py $P 30000000
run stage6 python3 -I -B $T/stage6.py $P 100000000
run stage7 python3 -I -B $T/stage7.py $P
run headline python3 -I -B $T/headline.py $P
run along_bar python3 -I -B $T/along_bar.py $P 30000000
run h19_xy python3 -I -B $T/h19_xy.py $P
run doublehit2 python3 -I -B $T/doublehit2.py $P 10000000
run trigger_diurnal python3 -I -B $T/trigger_diurnal.py $P
run block_fractions python3 -I -B $T/block_fractions.py $P 100000000
run sensitive_length python3 -I -B $T/sensitive_length.py $P 100000000
run doublehit_confirmed python3 -I -B $T/doublehit_confirmed.py $P 30000000
run stage8 python3 -I -B $T/stage8.py $P
S="root -l -b -q $T/proc_style.C"
run fig_sens $S "$T/proc_plot_sensitive_length.C(\"$P/work/fig_sensitive_length.json\",\"$P/figs/sensitive_length_hod2019\")"
run fig_eff $S "$T/proc_plot_eff.C(\"$P/work/fig_eff.json\",\"$P/figs/eff_vs_threshold\")"
run fig_stage3 $S "$T/proc_plot_stage3.C(\"$P/work/fig_stage3.json\",\"$P/figs\")"
run fig_stage4 $S "$T/proc_plot_stage4.C(\"$P/work/fig_stage4.json\",\"$P/figs/amp_spectra_hod2018_pe\")"
run fig_stage5 $S "$T/proc_plot_stage5.C(\"$P/work/fig_stage5.json\",\"$P/figs/double_hit_fraction\")"
run fig_stage7 $S "$T/proc_plot_stage7.C(\"$P/work/fig_stage7.json\",\"$P/figs/hod2019_rates\")"
run fig_along_bar $S "$T/proc_plot_along_bar.C(\"$P/work/fig_along_bar.json\",\"$P/figs/eff_along_bar\")"
run fig_diurnal $S "$T/proc_plot_diurnal.C(\"$P/work/fig_trigger_diurnal.json\",\"$P/figs/trigger_rate_vs_time\")"
# PROC-4: las de cuatro paneles pasan a *_4p y se generan las de dos paneles (paneles 1,3; canales 13,2)
for b in amp_spectra_hod2018_pe eff_along_bar; do for e in pdf png; do mv $P/figs/$b.$e $P/figs/${b}_4p.$e; done; done
run fig_along_bar_2p $S "$T/proc_plot_along_bar.C(\"$P/work/fig_along_bar.json\",\"$P/figs/eff_along_bar\",\"1,3\")"
run fig_stage4_2p $S "$T/proc_plot_stage4.C(\"$P/work/fig_stage4.json\",\"$P/figs/amp_spectra_hod2018_pe\",\"13,2\")"
# stage8 otra vez, con figs/ completo (reescribe la línea proceeding_2026_10/ de outputs/INDEX.md: restaurarla a mano)
run stage8_rerun python3 -I -B $T/stage8.py $P
echo "$(date -u +%FT%TZ) fin" | tee -a $L/run.log
