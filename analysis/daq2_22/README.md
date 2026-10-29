# analysis/daq2_22 — caracterización del trigger nuevo (tejas START) con run 10

Todo es solo lectura sobre los crudos (`/data/daq/raw/hodoscope_2026_10`, segmentos cerrados según el `index.json`); las salidas van a
`/data/daq/derived/daq2_22/` (`out/`, `ntuples/`, `build/`, `logs/`). Informe: `/home/xplus/daq2_22_trigger_remote_2026-10-29_report.md`.

| Script | Etapa | Qué hace |
|---|---|---|
| `tile_scenarios.py` | 1 | MC del banco (cos³, hod2018 CAD): f_2018, rX/rX₈ y R/R₈ por escenario de tejas; a y ε del par con lo medido en DAQ2-21 |
| `temp_corr.py` | 2 | por segmento cerrado de run 10 (incremental) y por hora de run 8 (control): R, aX, aY, rX, IMon; regresiones ln R y ln rX contra ln aX |
| `timing_interp.C` | 3 | tiempo interpolado (parábola y CFD 50 %) de los impactos > 20 p.e. de hod2018, con `ttt_raw` |
| `timing_fold.py` | 3 | tiempo plegado módulo 32 ns: una gaussiana, dos gaussianas, pico + meseta; dependencia con la tira |
| `imon_steps.py` | 4 | escalones de IMon (> 0.3 µA en < 1 min) en el slow control de 1 s; IMon contra aX por bloque |

Relanzamiento de las etapas 2 y 4 con los segmentos nuevos (desde la raíz del repositorio):

    ( ulimit -v 4194304; nice -n 10 timeout 3600 python3 -I analysis/daq2_22/temp_corr.py [--since 2026-10-29T12:00] )
    ( ulimit -v 4194304; nice -n 10 timeout 1800 python3 -I analysis/daq2_22/imon_steps.py 1793222461 /data/daq/derived/daq2_22/out )

`temp_corr.py` guarda los bloques ya procesados en `out/temp_corr_blocks_run010.json` y solo extrae los segmentos nuevos;
`--since` restringe la regresión. Usa la ganancia de todo run 10 de DAQ2-21 (`/data/daq/derived/daq2_21/out/run10_gain.csv`),
fija, porque el ajuste por segmento (~1000 eventos) no converge.
