# Notas para el chat de análisis

Escritas desde el trabajo del DAQ v2 (tareas DAQ2-0 y DAQ2-1, 2026-10-18). No se modificó
código de análisis; cada punto queda para que el chat de análisis decida.

## 1. `ch_to_plane()` y `ch_to_strip()` solo valen para `kCfg_2026_05`

`analysis/src/hodo_mapping.h:114-129` fija los rangos de 2026-05: 19X = ch 48-63 y
19Y = ch 32-47. Con `kCfg_2026_10` (la configuración de la campaña nueva, X de hod2019 =
ch 32-47 e Y = ch 48-63) esas funciones devuelven plano y tira equivocados para todo
hod2019. Lo mismo pasa con `hod2019_x_index()`/`hod2019_y_index()` (líneas 134-135), y
todas usan las inversas de `kH19X`/`kH19Y`, que son alias de `kCfg_2026_05` (líneas 70-73).
Hay que parametrizarlas por configuración antes de analizar datos de octubre.

## 2. AGENTS.md contradice `hodo_mapping.h` para hod2018

AGENTS.md (sección "Channel mapping") dice "hod2018 uses consecutive channels: X bar k =
channel 0-15 with k = channel". `kCfg_2026_05` (`hodo_mapping.h:44-50`), que usan v1-v10,
tiene X = {15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14} e
Y = {16, 24, 17, 25, ...}, medidos con correlaciones de vecinos. El handoff
(`docs/handoff_analysis_2026-10-18.md`, sección 5) confirma que el orden secuencial está mal
para mayo. AGENTS.md debería corregirse o marcar esa frase como superada.

## 3. El pickup de hod2018 en 3510-3525 ns no aparece en v8

Estudio de solo lectura sobre v8 (9057 eventos, 20 bloques de 700 entradas, fuera de los
intervalos excluidos), `outputs/daq2_window_study_20261018/`:

- Posición de los pulsos > 3 p.e. de hod2018 por muestra (`window_study_pulsepos_hod2018.csv`):
  plana entre 3280 y 4112 ns, |z| ≤ 2.4 respecto del nivel base de 1507.7 pulsos por muestra
  en [30,200). El muón aparece en las muestras 260-264 (4160-4224 ns), z hasta +14.0.
- Forma de onda media bajo la línea base en las muestras 210-229 (3360-3664 ns), mediana
  sobre canales de hod2018: entre 1.312 y 1.408 mV, sin estructura (es el efecto medio de
  los pulsos oscuros).
- El "t_peak = 3525 ns" de Ch2 y Ch9 en
  `outputs/timesync_test_20260531_174521/summary_timesync.txt` cae justo al inicio de la
  ventana de búsqueda [3500, 7000] ns de ese análisis, lo que sugiere un efecto de borde y no
  una señal [I].

Si el análisis tiene otra evidencia del pickup (otro run, otra amplitud), conviene
indicarla, porque la restricción de evitar 3510-3525 ns condicionaba la ventana del DAQ.

## 4. Cuantización de la línea base

El paso del ADC es 2000/4095 = 0.4884 mV. En hod2019 (ruido ~0.35 mV) la mediana de pocas
decenas de muestras cae casi siempre en el mismo código: en v8 la MAD de la mediana de 27 o
40 muestras dio 0 para casi todos los canales de hod2019. Para pedestales de hod2019 conviene
una media recortada (en el estudio: promedio de las muestras a ±2 mV de la mediana), que dio
σ = 0.2495 mV (mediana de canales, 40 muestras).

## 5. El pedestal posterior a la señal está desplazado en canales con muón

En v8, para canales de hod2018 con pulso en tiempo > 10 p.e., la mediana de las muestras
311-365 (≥ 624 ns después del fin de la ventana de señal) queda +3.907 mV sobre la línea
base de referencia (mediana sobre canales), por el sobreimpulso. Cualquier pedestal que
mezcle muestras posteriores a la señal está sesgado en esos eventos.

## 6. Control lento: `hv_full_20261019.csv` mezcla dos formatos (DAQ2-4, 2026-10-19)

`/data/daq/slowcontrol/2026/hv_full_20261019.csv` tiene el encabezado de 24 columnas (corrida de
DAQ2-3, filas 2-1345, última `t_req_utc_ns` = 1792374668848097448), pero desde la fila 1346
(2026-10-19T08:32:37Z, `hvmon` con las columnas de la placa) las filas tienen 28 columnas: al final
`bd_BdStatus,bd_HVMax_V,bd_HIMax_A,bd_Temp_degC`. El `.schema.json` describe solo las 24. Para
leerlo, separar por número de campos. Desde DAQ2-4 el escritor ya no mezcla: abre
`hv_full_YYYYMMDD.2.csv` (diseño 4.3). El archivo no se reescribió.
