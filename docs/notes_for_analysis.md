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

## 6. Control lento: `hv_full_20261019` reparado en dos partes (DAQ2-4, 2026-10-19)

`/data/daq/slowcontrol/2026/hv_full_20261019.csv` mezclaba dos formatos: encabezado de 24
columnas (corrida de DAQ2-3, filas 2-1345, última `t_req_utc_ns` = 1792374668848097448) y, desde
la fila 1346 (2026-10-19T08:32:37Z, `hvmon` con las columnas de la placa), filas de 28 columnas
(al final `bd_BdStatus,bd_HVMax_V,bd_HIMax_A,bd_Temp_degC`). Tras parar `hvmon` se separó:

| archivo | filas de datos | columnas | esquema | bytes | crc32 |
|---|---|---|---|---|---|
| `hv_full_20261019.csv` | 1344 | 24 | `.schema.json` existente (`daqv2.hv_full/1`) | 160328 | e961c2c5 |
| `hv_full_20261019.2.csv` | 1092 | 28 | nuevo `.2.csv.schema.json` (`daqv2.hv_full/2`) | 150829 | ddd5d531 |

1344 + 1092 = 2436 = filas de datos del original (2437 líneas con el encabezado); las filas son
idénticas y en el mismo orden. Original (310910 bytes, crc32 2f92343c) y su esquema guardados en
`/data/daq/slowcontrol/repair_backup/`. Es el mismo esquema de nombres que usa el escritor
corregido: los lectores deben tomar todas las partes del día (`hv_full_YYYYMMDD*.csv*`) y ordenar
por `t_resp_utc_ns`.

## 7. Código compartido con el DAQ (DAQ2-5, 2026-10-19)

- `daq/core/include/daq/pulse_finder.h` (C++17, header-only, sin ROOT): búsqueda de pulsos del DAQ
  (mínimo local, separación mínima, media recortada de la línea base; funciona con muestras ADC o
  en mV). `daq/core/include/daq/trigger_interval_stats.h`: tasa por máxima verosimilitud, KS y χ²
  de exponencialidad, huecos y ráfagas, con resolución declarada para v7-v10. Se usan con
  `-I daq/core/include`; no hace falta el build de CMake.
- Lector de control lento a TTree (`libdaq_screader.so`, planos y `.gz`, partes `.N`, esquemas):
  ejemplo `analysis/macros/plot_hv_imon_day.C`.
- Dato del replay de v8 (entradas 30000-34999): tasa fuera de tiempo de hod2018 > 3 p.e. (ganancias
  de v8, 91 muestras de la ventana B) mediana 377.75 kHz con el buscador del DAQ y 378.85 kHz con
  el método de DAQ2-1 sobre los mismos eventos; hod2019 > 3 mV 29.67 / 29.26 kHz. El estudio de
  DAQ2-1 (otros 9057 eventos repartidos en el run) dio 331.4 y 24.27 kHz: las tasas cambian a lo
  largo de v8. No se analizó por qué.

## 8. Tiempos muertos de trigger en v7-v10: posible apagado térmico del crate NIM (DAQ2-6, 2026-10-19)

El 2026-10-19 el trigger externo no llegaba al VX1740: la fuente del crate NIM (CAEN) se había
apagado por sobretemperatura (la sala cerrada se calienta de noche o con el sol), dejando sin
alimentación el discriminador N845 y el AND N405. Al volver a encenderla, los triggers volvieron.
Los tiempos muertos de trigger de v7-v10 se atribuyeron al cable; pueden venir también (o en
cambio) de estos apagados. Sugerencia: listar los huecos de trigger de v7-v10 (Δt mayores que el
umbral de `trigger_interval_stats.h`, ojo con que `timestamp` es por lectura y en hora UTC-4) y ver
si se agrupan en las horas más calurosas del día o con sol directo, y si su duración se parece a la
de un enfriamiento de la fuente (no a un corte de cable intermitente). desde DAQ2-6 daqrun avisa por
correo tras 5 min sin triggers en runs con trigger externo, con la temperatura de la placa HV.
