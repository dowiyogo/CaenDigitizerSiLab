# Resumen de los resultados vigentes, seguimientos 2026-10-31 (1) a (8)

Resumen único, sin números nuevos. Cada valor está copiado del informe que lo dio, con su sección. Los informes están en `/home/xplus/`:

| Seguimiento | Informe |
|---|---|
| (1) | `xplus_followup_2026-10-31_report.md` |
| (2) | `xplus_proc_hvdiag5_2026-10-31_report.md` |
| (3) a (8) | `xplus_followup<N>_2026-10-31_report.md` |

Etiquetas: [M] medido, [I] inferido, [S] supuesto, como en cada informe. Escrito el 2026-10-31 (seguimiento 9).

## Reglas que valen para todo lo que sigue

- **Comparaciones de Y entre mayo y octubre: por barra (MPPC), no por canal.** Entre mayo y octubre se recableó el plano Y de hod2018 (X no se tocó).
  - Cada MPPC de Y quedó en otro canal del digitalizador con su misma línea de HV.
  - El emparejamiento por barra (tira i del orden de mayo ↔ tira i del orden de octubre) coincide con el grupo HV medido en las dos campañas en 5 de 5 canales, y en 0 de 5 con la orientación invertida [M] (seguimiento 7, §1.1).
  - Tabla: `outputs/hod2018_recabling/2026-10-31_groups/ratio_by_mppc.csv`.
- **hod2019 no tiene señal desde el 2026-10-21, entre las 17:13 y las 20:51 UTC (dentro del run 4)** [M] (seguimiento 5, §1.3).
  - En octubre los 32 canales tienen el RMS de un canal sin entrada (0.186-0.257 mV, como ch3: 0.225 mV), sin pulsos oscuros ni de muones [M] (seguimiento 5, §1.1-1.2).
  - La HV estaba nominal y hod2018, en los mismos grupos, siguió normal. La causa más probable es la cadena de lectura o la distribución de bias de hod2019 [I] (seguimiento 5, §1).
  - **Todo resultado de octubre es solo de hod2018.**
  - El paso e.0 del checklist del lunes (`docs/lab_checklist_2026-11-02.md`) dice qué revisar (seguimiento 5, §1.5).
- **El sobrevoltaje nominal de hod2018 es de ~1.2 V** [I], con una estimación de dos puntos en 5 canales: OV de 1.02 a 1.47 V, mediana 1.24 V; k mediana 3.26 mV/V (seguimiento 8, §1.2).
  - A −1.0 y −1.5 V los pulsos oscuros desaparecen; G solo es medible a nominal y −0.5 V [M] (seguimiento 8, §1.2).
  - Con un OV así, ±0.1 V es ±8 % del OV [I] (seguimiento 8, §1.2).

## Resultados vigentes

### HV, trip de CH03 y barridos

| Resultado | Valor vigente | Informe, sección |
|---|---|---|
| I0Set, RUp, RDWn y Trip de CH00-07 | 15.0 µA, 1 V/s, 50 V/s, 5 s [M] | (1) §1b |
| CH03 a nominal antes del trip | 7.552 µA (margen 7.448 µA) [M] | (1) §1b |
| Trip de CH03, 2026-10-30 | corriente **estable** sobre I0Set: limitada en 15.45 µA durante 5 s hasta el trip a las 02:08:30 UTC; CH00 y CH02 llegan a una meseta en 6-8 s y no decaen en 30 s [M en el registro, I en la lectura] | (8) §2.1 |
| Rampa lenta (RUp 0.1 V/s) para evitar el trip | no lo evitaría en CH03 a +0.25 V [I] | (8) §2.1 |
| Corriente contra voltaje, modelo con piso de fuga (grupos 1, 2, 3, 5, 7 y 8) | V_bd (VMon) 63.3-64.4 V, OV nominal 1.2-2.2 V, n 3.1-6.5; esperado a +0.5 V: 10-20 µA; +0.75 V: 17-30 µA; +1.0 V: 28-45 µA; grupos 4 y 6 no determinados [I] | (8) §2 |
| RUp mínimo que informa la A7030DP | 1.0 V/s (rango 1-500 V/s, `hvctl params` del 2026-10-18) [M de esa lectura] | (7) §2 |
| Escritura de RUp en hvmon | rama `daq/hvmon-rup-2026-10-31`, sin desplegar; 0.05-5 V/s y rango de la placa: 0.1 V/s es rechazado [M en el simulado] | (7) §2 |
| Fallos del enlace del digitalizador (2026-10-23 y 2026-10-30) | primer error `ReadData: CAEN_DGTZ_CommError (−1)` a mitad de una adquisición; los datos no apuntan a la toma SW ni a la rampa de HV como causa [M] | (3) §3 |

### Eficiencia y proceeding (tiempo bueno nuevo, HV-DIAG-5 excluido)

| Resultado | Valor vigente | Informe, sección |
|---|---|---|
| Tiempo bueno de mayo v8 | 77 558 triggers (antes 113 152) [M] | (2) §1 |
| ε_XY a 5 p.e. | 0.79 ± 0.08 ± 0.07 ± 0.04 [M] | (2) Resumen; (1) §3 |
| ε_XY a 10 p.e. | 0.61 ± 0.02 ± 0.02 ± 0.03 [M] | (2) Resumen; (1) §3 |
| «Signos opuestos» de las capas X | a 10 p.e., B − A a 2.8σ en mayo y 3.0σ en octubre (3.9σ juntas); a 5 p.e. no se sostiene (1.0-1.2σ) [M] | (2) §4 |
| Eficiencia absoluta por capa y posición, sin dividir por ε_MC | octubre: 0.20-0.33 (5 p.e.) y 0.11-0.24 (10 p.e.); mayo: 0.23-0.53 (5 p.e.) y 0.17-0.30 (10 p.e.); ε_MC 0.31-0.36 [M] | (6) §2 |
| Afirmaciones de `main_v5.tex` sobre octubre | ningún número de hod2019 ni de dos estaciones viene de octubre; las líneas 201 y 268 dan a entender que funcionaban las dos estaciones [M] | (5) §1.4 |

### Amplitud, crosstalk y luz

| Resultado | Valor vigente | Informe, sección |
|---|---|---|
| L_eff por capa (pendientes de los canales) | octubre: XA 136 ± 42, XB 131 ± 44, YA 105 ± 18, YB 175 ± 43 mm; capas de cada coordenada con signos opuestos [M][I] | (3) §2 |
| Espectro de coincidencias X–Y por posición (octubre) | mediana XA 17.5 → 10.5 p.e. y XB 14.5 → 9.5 al alejarse del extremo de lectura; Y plana o débil [M] | (4) §2 |
| DCR > 0.5 p.e. de hod2018 | mediana 2209 ± 57 kHz en octubre, 2641 ± 79 kHz en mayo [M] | (4) §1 |
| Colas partidas por el buscador (Δ = 3-6 muestras) | 0.523 de los pares; f = 0.159-0.173 en la región pre-trigger de octubre [M] | (5) §3.1-3.2 |
| 1/(1 − p_ct), canal mediano de octubre | **1.38 a 1.52, central 1.46** [I] (corregido por colas partidas) | (5) §3.3 |
| Longitud de camino en la barra sonda | mediana 1.015 mm, media 0.981 mm (MC de PROC) [M] | (5) §2 |
| Razón mayo/octubre en p.e. primarios, suma de canales | 1.276 ± 0.032 [M]; el 1.44 en p.e. se reparte en ~1/3 de crosstalk (1.13) y ~2/3 de luz o PDE (1.28) [I] | (6) §1 |
| Razón mayo/octubre por MPPC, por plano | X 1.208 ± 0.057, Y 1.207 ± 0.051; X − Y = +0.001 ± 0.077 [M]: no sigue al recableado | (7) §1.2 |
| Razón por grupo HV | grupos 1-3: 1.36-1.46; grupos 4-6: 1.09-1.16; p = 0.012 entre grupos (≈ 0.08 corregido por 7 agrupaciones) [M][I] | (7) §1.2 |
| Forma del pulso de 1 p.e. | cola más larga en octubre en los 31 MPPC (ΔT(6) +0.060 a +0.131), igual en X e Y; depende del grupo HV (p = 1.1·10⁻⁹) y sigue al MPPC, no al canal [M][I] | (7) §1.3 |
| Razón de p.e. primarios contra ΔOV | Pearson +0.29 a +0.35 (t = 1.46-1.78), Spearman +0.11 a +0.19, n 25-29: no significativa; ΔOV mediana +0.19 a +0.20 V (OV más alto en mayo) [M] | (8) §3 |

### hod2019 en mayo

| Resultado | Valor vigente | Informe, sección |
|---|---|---|
| G de hod2019 en mayo | el ajuste de c1/PROC no converge en ningún canal [M]; la cola oscura sugiere G ≈ 0.5-2.5 mV/p.e. [I] | (6) §3 |
| Pulsos de hod2019 sincronizados con el trigger, mayo | medianas de 4.9 a 12.4 mV, percentiles 90 de 10 a 90 mV [M] | (6) §3 |

## Erratas

Valores que un informe posterior corrigió:

| Valor anterior | Valor vigente | Motivo | Informe que corrige |
|---|---|---|---|
| `good_intervals` con 4 filas | 6 filas: faltaban dos HV_OFF de v10 en junio | un `grep` filtraba por fechas de mayo | (1), errata del informe del 2026-10-30 |
| Tiempo bueno de mayo v8: 113 152 triggers | 77 558 | episodios de HV-DIAG-5 excluidos (intervalo completo) | (2) §1 |
| ε_XY a 5 p.e.: 0.797 ± 0.066 ± 0.095 ± 0.037 | 0.785 ± 0.084 ± 0.070 ± 0.037 (0.79 ± 0.08 ± 0.07 ± 0.04) | ídem | (1) §3; (2) |
| ε_XY a 10 p.e.: 0.602 | 0.614 ± 0.024 ± 0.015 ± 0.029 (0.61 ± 0.02 ± 0.02 ± 0.03) | ídem | (1) §3; (2) |
| `main_v5.tex` l. 286: 15 tiras X vivas sobre 3σ a 5 p.e. | 14 | ídem | (2) §3 |
| `main_v5.tex` l. 352: referencia ch29 «about 0.38 for both X layers» | 0.56 ± 0.27 (capa A) y 0.37 ± 0.23 (capa B) | ídem | (2) §3 |
| Razones de la capa X B a lo largo de la barra: 1.27 y 1.31; pendiente +0.7 %/mm | 1.08 y 1.15; +0.2 %/mm | ídem | (2) §3 |
| Trip de CH03: «el salto es de la rampa misma» (corriente de carga, ver también (3) §3) | corriente estable sobre I0Set; la rampa lenta no lo evita | registro de 1 s: CH00 y CH02 en meseta sin decaer | (8) §2.1 |
| 1/(1 − p_ct) de octubre: 1.360 ± 0.019 (mediana) | 1.38-1.52, central 1.46 | 16-28 % de los pulsos > 0.5 p.e. son colas partidas de un pulso de ~1 p.e. | (5) §3.3 |
| ch11 y ch26: canales con ganancia ≤ 5 % distinta entre campañas | solo ch11 (X04); ch26 es otro MPPC en cada campaña: G_mayo/G_oct = 1.178 por MPPC | recableado de Y | (7) §1.4 |
| Razón en p.e. primarios con G igual: ch11 + ch26 1.26 ± 0.13 | ch11 solo: 1.26 ± 0.20 (1.3σ de 1); se debilita «la pérdida no es sobrevoltaje» | ídem | (7) §1.4 |
| ch29 (octubre) y ch30 (mayo): canales anómalos | son el mismo MPPC (Y13), con razón normal 1.129 ± 0.102 | ídem | (7) §1.1 |
| ch17, ch20 y ch22 «cambiaron de grupo HV» entre campañas | no cambiaron los MPPC; cambió qué MPPC lee cada canal | ídem | (7) §1.1 |
| Mediana de p.e. primarios con bins de 0.25 p.e. y recorte (mayo 17.4 en p.e. medidos) | definición del seguimiento 1 (bins de 1 p.e.); 19.5 en mayo | el recorte por bin depende del ancho del bin | (6), Notas |
| RUp = 0.1 V/s para el barrido lento (seguimientos 3 y 7) | la placa informa mínimo 1.0 V/s; hvmon lo rechazaría | rango del parámetro en la placa | (7) §2 |
| Tasa de run B en tiempo vivo (huecos restados): 0.590 Hz | tasa estacionaria con t ≥ 5.25 h: 0.5947 ± 0.0017 Hz (142 544 triggers en 66.58 h, sin huecos) | los 15 huecos y la subida de la tasa están antes de 5.25 h | «Run B start», más abajo |
| Modelo de longitud sensible de hod2019 (PROC-3, `sensitive_length.py`): L ≈ 21 mm reproduce r_X a 35 mV; 86 X·Y predichas | z ≈ 3.2 mm (X) y 2.1 mm (Y); 0 X·Y predichas (< 0.16 al 95 %) | el histograma 2D descartaba las trayectorias sin cruce en la otra coordenada: su P_X era P(X ∧ alguna Y) | seguimiento 12, punto C |

## Pendientes que siguen abiertos

- **V_bd por canal con 3 o más voltajes:** necesita una toma nueva. El barrido está preparado, sin armar, en `daq/biasscan-vbd-2026-10-31` (`docs/biasscan_vbd_2026-10-31.md` en esa rama). Se arma solo con una decisión de René (seguimiento 9).
- **Causa de la pérdida de hod2019:** la inspección del 2026-11-02 (checklist e.0).
- **Sentido absoluto del orden de tiras de hod2018:** sin fijar. La prueba de orientación está preparada y no ejecutada (seguimiento 1, §5).

## Run B start

Comprobación del 2026-11-01 sobre run B (v8, 143 814 triggers), solo lectura. t se mide desde el primer trigger de v8 (2026-05-28 14:11:19 UTC). Todos los números son [M].
- Ítem 1: `may_common.good_time` sobre `outputs/proceeding_hvdiag5/2026-10-31_full/work/hd64_may8_hits.bin`.
- Ítems 2 a 4: rama `timestamp` de `lectura_vx1740_ext_nim_1M_v8.root`; el ajuste y el KS son los de `fit_exp_ks()` de `analysis/src/poisson_analysis.cxx`.

**1. Tiempo bueno antes de 5.25 h.**
- Hay 774 triggers con t < 5.25 h dentro del tiempo bueno (77 558 triggers en total).
- El primero está en t = 4.8121 h; el intervalo EXCLUDE «start of run» termina en 4.8114 h.
- Con t < 5.25 h hay 1270 triggers en total.

**2. Huecos de DAQ (ΔT > 50 T/(N−1) = 89.906 s).**

| # | inicio (h) | fin (h) | duración (s) |
|---|---|---|---|
| 1 | 0.0738 | 0.9399 | 3118.0 |
| 2 | 0.9634 | 2.5442 | 5691.1 |
| 3 | 2.6808 | 2.7101 | 105.6 |
| 4 | 2.7583 | 2.7857 | 98.6 |
| 5 | 2.7884 | 2.8515 | 227.2 |
| 6 | 2.8545 | 2.8843 | 107.4 |
| 7 | 2.9404 | 2.9916 | 184.3 |
| 8 | 2.9916 | 3.0250 | 120.4 |
| 9 | 3.0250 | 3.1413 | 418.8 |
| 10 | 3.1413 | 3.1671 | 92.8 |
| 11 | 3.1767 | 3.2127 | 129.7 |
| 12 | 3.2413 | 3.3084 | 241.4 |
| 13 | 3.3231 | 3.3546 | 113.1 |
| 14 | 3.3572 | 3.4471 | 323.4 |
| 15 | 3.4471 | 4.5112 | 3830.8 |

En total suman 14 802.6 s. Los 15 huecos, con sus 14 802.6 s, caen dentro de las primeras 5.25 h.

**3. Tasa de triggers en bines de 15 min, primeras 6 h.** La tasa es el número de triggers dividido por 900 s, sin corregir por huecos. La última columna es el tiempo de huecos dentro del bin.

| bin (h) | triggers | tasa (Hz) | huecos en el bin (s) |
|---|---|---|---|
| 0.00-0.25 | 8 | 0.0089 | 634.4 |
| 0.25-0.50 | 0 | 0.0000 | 900.0 |
| 0.50-0.75 | 0 | 0.0000 | 900.0 |
| 0.75-1.00 | 3 | 0.0033 | 815.4 |
| 1.00-1.25 | 0 | 0.0000 | 900.0 |
| 1.25-1.50 | 0 | 0.0000 | 900.0 |
| 1.50-1.75 | 0 | 0.0000 | 900.0 |
| 1.75-2.00 | 0 | 0.0000 | 900.0 |
| 2.00-2.25 | 0 | 0.0000 | 900.0 |
| 2.25-2.50 | 0 | 0.0000 | 900.0 |
| 2.50-2.75 | 30 | 0.0333 | 264.9 |
| 2.75-3.00 | 11 | 0.0122 | 647.8 |
| 3.00-3.25 | 9 | 0.0100 | 762.5 |
| 3.25-3.50 | 7 | 0.0078 | 837.3 |
| 3.50-3.75 | 0 | 0.0000 | 900.0 |
| 3.75-4.00 | 0 | 0.0000 | 900.0 |
| 4.00-4.25 | 0 | 0.0000 | 900.0 |
| 4.25-4.50 | 0 | 0.0000 | 900.0 |
| 4.50-4.75 | 320 | 0.3556 | 40.2 |
| 4.75-5.00 | 422 | 0.4689 | 0.0 |
| 5.00-5.25 | 460 | 0.5111 | 0.0 |
| 5.25-5.50 | 543 | 0.6033 | 0.0 |
| 5.50-5.75 | 484 | 0.5378 | 0.0 |
| 5.75-6.00 | 497 | 0.5522 | 0.0 |

**4. Solo t ≥ 5.25 h.**
- Hay 142 544 triggers y 142 543 intervalos.
- Ningún intervalo supera el umbral de huecos, ni el del run completo (89.906 s) ni el propio de este conjunto (50 T/(N−1) = 84.077 s).

| | todo el run | t ≥ 5.25 h |
|---|---|---|
| 1/⟨ΔT⟩ | 0.5898 Hz | 0.5947 Hz (⟨ΔT⟩ = 1.681544 s) |
| λ_fit | 0.5943 ± 0.0017 Hz | 0.5947 ± 0.0017 Hz |
| KS | D = 0.00300, p = 0.1502 | D = 0.00265, p = 0.2688 |

## Objetivos por barra para Geant4 (2026-11-01)

`analysis/followup10_2026_11_01/bar_targets.py`; tabla en `outputs/hod2018_bar_targets/2026-11-01_g4/bar_targets.csv` (SHA256 `1419f9a7…6339b`). Todos los números son [M].
- **Filas:** 32 barras × 2 campañas × 4 intervalos de distancia al extremo de lectura (21-29, 29-37, 37-45, 45-54 mm) = 256 filas. Las 8 de X12 (ch3, muerto) van sin datos.
- **Selección:** alguna tira ortogonal del intervalo > 10 p.e. en el mismo evento, menos mezclados. Amplitud de la sonda > 5 p.e. medidos. Mediana en bins de 1 p.e. medido, interpolada dentro del bin. p.e. primarios = medidos / F, con F = 1.75 (mayo) y 1.46 (octubre).
- **Mediana fiable si n_net ≥ 15:** desde ahí el error relativo por remuestreo (200) queda bajo ~0.3 (una fila 0.34); debajo llega a 0.9-5.9.
- **ε:** q/r_ref con mezclados, a 5 y 10 p.e. ε_MC (barra cruzada dada la etiqueta) va de 0.0388 a 0.0475.

| | mayo v8 (77 558 triggers) | octubre run 8 (44 176) |
|---|---|---|
| filas con datos | 124 | 124 |
| filas fiables (n_net ≥ 15) | 115 | 79 |
| n_net por fila: mediana (máx.) | 42.5 (81.0) | 18.5 (35.0) |
| mediana en p.e. primarios, filas fiables (mediana) | 9.81 | 8.65 |
| ε/ε_MC a 5 p.e.: mediana (error mediano) | 0.950 (0.258) | 0.855 (0.300) |
| ε/ε_MC a 10 p.e.: mediana (error mediano) | 0.822 (0.201) | 0.557 (0.217) |
| L_eff determinado (de 31 barras) | 25 | 14 |
| L_eff: mediana (rango) | 88 mm (38-1435) | 82 mm (34-606) |
| L_eff con error < L_eff/2 | 6 | 1 |
| no determinado: < 3 intervalos fiables / pendiente ≥ 0 | 2 / 4 | 13 / 4 |

**L_eff común por campaña** (media ponderada de L_eff con 1/σ², χ² contra el valor común) [M]:

| | barras | ⟨L_eff⟩ | χ²/ndf | comprobación en 1/L_eff: 1/⟨1/L⟩ | χ²/ndf |
|---|---|---|---|---|---|
| mayo, error < L_eff/2 | 6 (X06, X08, X11, X15, Y09, Y14) | 43.6 ± 5.7 mm | 2.07/5 = 0.41 | 49.9 +7.6 −5.8 mm | 3.28/5 |
| mayo, todas las determinadas | 25 | 48.6 ± 5.3 mm | 9.37/24 = 0.39 | 82.4 +11.1 −8.7 mm | 22.02/24 |
| octubre, error < L_eff/2 | 1 (X02) | 33.6 ± 12.7 mm | — | 33.6 +20.5 −9.2 mm | — |
| octubre, todas las determinadas | 14 | 43.6 ± 10.6 mm | 3.17/13 = 0.24 | 89.0 +31.5 −18.4 mm | 5.33/13 |

El ajuste da la pendiente 1/L_eff con error gaussiano; σ_L = σ_k · L_eff², así que el peso 1/σ_L² y el corte «error < L_eff/2» favorecen los L_eff cortos.

## Espectros MIP etiquetados, amplitud y carga (2026-11-01)

`analysis/followup10_2026_11_01/{charge_extract.C, mip_tagged.py, mip_fit.C}`; salida en `outputs/hod2018_mip_tagged/2026-11-01_tagged/` (figura `mip_tagged_spectra`). Todos los números son [M].
- **Selección:** etiqueta = alguna tira ortogonal del intervalo > 10 p.e.; sonda = cada tira viva de la capa > 3 p.e. medidos en su variable; mismo evento − mezclados; todas las tiras de la capa sumadas. Mayo v8 en tiempo bueno, octubre run 8 completo. p.e. primarios con F = 1.75 / 1.46.
- **Carga:** Σ (p_loc − x) · 16 ns en [k−3, k+12] muestras alrededor del mínimo del pulso de la ventana signal, pedestal local p_loc = media de [k−12, k−4]. 1 p.e. en carga = mediana de pulsos oscuros aislados (sin otro pulso > 0.5 G en ±16 muestras) de amplitud 0.5-1.5 G, fuera de la ventana signal: 491-658 mV·ns en mayo y 395-575 mV·ns en octubre (Q1/G ≈ 115-130 ns).
- **Ajuste:** Landau ⊗ Gauss (langaus) desde el primer bin sobre el umbral (2 p.e. primarios en mayo, 3 en octubre) hasta 60. «n.d.»: ajuste en el límite de un parámetro o error de MPV mayor que el MPV.

MPV de la Landau (p.e. primarios) y razón MPV / mediana de la selección del seguimiento 6 (por capa, F de campaña; mayo XA 10.57, XB 11.71, YA 12.29, YB 10.00; octubre XA 9.25, XB 9.93, YA 8.56, YB 8.56):

| Campaña | Capa | Intervalo | Amplitud: MPV | MPV/mediana | Carga: MPV | MPV/mediana |
|---|---|---|---|---|---|---|
| mayo | XA | 21-29 | 7.40 ± 1.60 | 0.70 | 7.85 ± 1.38 | 0.74 |
| mayo | XA | 29-54 | 5.84 ± 1.37 | 0.55 | 6.19 ± 1.40 | 0.59 |
| mayo | XB | 21-29 | 8.12 ± 1.55 | 0.69 | 8.11 ± 1.72 | 0.69 |
| mayo | XB | 29-54 | 7.18 ± 0.69 | 0.61 | 7.72 ± 0.93 | 0.66 |
| mayo | YA | 21-29 | 7.57 ± 1.61 | 0.62 | 8.84 ± 1.00 | 0.72 |
| mayo | YA | 29-54 | 6.21 ± 1.46 | 0.51 | 4.95 ± 3.35 | 0.40 |
| mayo | YB | 21-29 | 4.83 ± 3.48 | 0.48 | 5.63 ± 1.51 | 0.56 |
| mayo | YB | 29-54 | 5.79 ± 1.00 | 0.58 | n.d. | — |
| octubre | XA | 21-29 | 9.19 ± 1.79 | 0.99 | 8.65 ± 1.92 | 0.94 |
| octubre | XA | 29-54 | 3.59 ± 2.45 | 0.39 | n.d. | — |
| octubre | XB | 21-29 | n.d. | — | n.d. | — |
| octubre | XB | 29-54 | 5.05 ± 1.40 | 0.51 | 6.01 ± 0.55 | 0.61 |
| octubre | YA | 21-29 | n.d. | — | n.d. | — |
| octubre | YA | 29-54 | n.d. | — | 5.66 ± 1.27 | 0.66 |
| octubre | YB | 21-29 | 4.99 ± 1.84 | 0.58 | n.d. | — |
| octubre | YB | 29-54 | n.d. | — | 4.48 ± 1.41 | 0.52 |

Pico visible a ojo (figura completa `mip_tagged_spectra_all`): mayo XA, XB e YA en los dos intervalos (~8-10 p.e.), YB en ~5-6 p.e., cerca del umbral; octubre solo XA 21-29 (~10 p.e.) y XB 29-54 (~6 p.e.); en el resto de octubre los espectros suben hasta el umbral. Los parámetros de ajuste completos (anchura, σ, χ²/ndf, máximo de la función) están en `fits.csv`.

**Suma de las cuatro capas de mayo** (`mip_sum.C`; netos sumados, errores en cuadratura; mismo ajuste desde 2 hasta 60 p.e. primarios; `fits_may_sum.csv`) [M]:

| Intervalo | Variable | n_net | MPV | Anchura de Landau | σ de la Gauss | χ²/ndf | Máximo de la función | f(> 5 p.e. medidos) | f(> 10 p.e. medidos) |
|---|---|---|---|---|---|---|---|---|---|
| 21-29 mm | amplitud | 1370 | 5.51 ± 1.96 | 0.90 ± 0.31 | 7.70 ± 1.72 | 70.2/54 | 7.36 | 0.881 | 0.737 |
| 21-29 mm | carga | 1280 | 6.50 ± 1.59 | 1.35 ± 0.34 | 6.65 ± 1.72 | 55.8/54 | 8.60 | 0.906 | 0.778 |
| 29-54 mm | amplitud | 3490 | 1.92 ± 1.46 | 0.65 ± 0.16 | 8.37 ± 1.07 | 76.7/54 | 3.21 | 0.822 | 0.639 |
| 29-54 mm | carga | 3501 | 5.14 ± 1.49 | 1.26 ± 0.26 | 5.74 ± 1.56 | 51.4/54 | 7.00 | 0.880 | 0.721 |

- f(> T) = fracción de la función ajustada sobre T/1.75 p.e. primarios (2.86 y 5.71), normalizada a su integral en [0, 150]: por debajo de 2 p.e. primarios la función está extrapolada.
- Referencia: ε_XY = 0.79 (5 p.e.) y 0.61 (10 p.e.). La comparación es aproximada: ε_XY exige las dos coordenadas sobre el umbral.
- En 29-54 mm la amplitud queda con una Landau angosta y una Gauss ancha (σ 8.4 p.e.), con el máximo a 3.21 p.e., cerca del inicio del ajuste.

**Versión para presentación** (`mip_talk.C`; figura `mip_tagged_spectra_talk.{png,pdf}`: mayo, 29-54 mm, carga, cuatro capas, 16 × 10 cm). Por capa e intervalo, moda de la función ajustada (no el parámetro MPV), fracción de la función ajustada sobre 5 y 10 p.e. medidos (2.86 y 5.71 primarios; normalizada a [0, 150], extrapolada bajo 2 p.e.) y su cuadrado, junto a ε_XY = 0.79 (5 p.e.) y 0.61 (10 p.e.); `talk_table.csv` [M]:

| Capa | Intervalo | Variable | Moda del ajuste | f(>5) | f(>5)² | f(>10) | f(>10)² |
|---|---|---|---|---|---|---|---|
| XA | 21-29 mm | carga | 9.08 | 0.964 | 0.929 | 0.862 | 0.743 |
| XA | 21-29 mm | amplitud | 8.81 | 0.950 | 0.902 | 0.834 | 0.695 |
| XA | 29-54 mm | carga | 7.74 | 0.907 | 0.823 | 0.749 | 0.560 |
| XA | 29-54 mm | amplitud | 7.32 | 0.900 | 0.809 | 0.728 | 0.529 |
| XB | 21-29 mm | carga | 10.18 | 0.940 | 0.884 | 0.841 | 0.707 |
| XB | 21-29 mm | amplitud | 10.06 | 0.936 | 0.875 | 0.829 | 0.688 |
| XB | 29-54 mm | carga | 8.68 | 0.972 | 0.945 | 0.868 | 0.754 |
| XB | 29-54 mm | amplitud | 8.02 | 0.974 | 0.948 | 0.850 | 0.723 |
| YA | 21-29 mm | carga | 10.06 | 0.980 | 0.961 | 0.906 | 0.820 |
| YA | 21-29 mm | amplitud | 9.31 | 0.925 | 0.856 | 0.802 | 0.643 |
| YA | 29-54 mm | carga | 6.67 | 0.868 | 0.753 | 0.704 | 0.495 |
| YA | 29-54 mm | amplitud | 7.68 | 0.899 | 0.808 | 0.735 | 0.540 |
| YB | 21-29 mm | carga | 5.63 | 0.944 | 0.892 | 0.743 | 0.552 |
| YB | 21-29 mm | amplitud | 5.21 | 0.898 | 0.806 | 0.689 | 0.475 |
| YB | 29-54 mm | carga | pico bajo 3 p.e. (2.34) | 0.819 | 0.671 | 0.642 | 0.412 |
| YB | 29-54 mm | amplitud | 6.03 | 0.958 | 0.917 | 0.756 | 0.572 |

f² es una aproximación de exigir las dos coordenadas sobre el umbral con la misma distribución; ε_XY = 0.79 y 0.61 son las de mayo.

## Revisión del proceeding (2026-11-02)

Código en `analysis/followup11_2026_11_02/`; salidas en `outputs/proceeding_revision/2026-11-02_revision/`. Todos los números son [M] salvo que se indique.

**1. Espectros MIP para el artículo: mayo, 29-54 mm, amplitud en p.e. medidos (sin dividir por F).** Selección de `mip_tagged.py`, bins de 1 p.e. medido (los espectros divididos por 1.75 reproducen `spectra.csv`, máx |Δ| 5e-5). El ajuste Landau ⊗ Gauss es el de la figura de la charla con el rango en p.e. medidos: de 4 a 105. Antes empezaba en 2 primarios = 3.5 medidos, y 4 es el borde de bin siguiente. Empezar en 3 incluye el bin [3, 4), que en YB tiene 109 ± 43 cuentas de exceso en el umbral, y el ajuste llega al límite (MPV 0.10). La moda es el máximo de la función ajustada. f(>T) es la fracción de la función por encima de T sobre [0, 262.5], extrapolada por debajo de 4.

| Capa | n_net (4-105) | Moda (p.e.) | χ²/ndf | f(>5) | f(>5)² | f(>10) | f(>10)² |
|---|---|---|---|---|---|---|---|
| XA | 787 | 13.44 | 75.1/74 | 0.916 | 0.839 | 0.747 | 0.557 |
| XB | 865 | 13.23 | 60.2/80 | 0.991 | 0.982 | 0.877 | 0.770 |
| YA | 892 | 13.85 | 81.1/77 | 0.949 | 0.901 | 0.803 | 0.644 |
| YB | 779 | 10.11 | 70.2/86 | 0.964 | 0.929 | 0.755 | 0.570 |

Referencia: ε_XY = 0.79 (5 p.e.) y 0.61 (10 p.e.). En YB el σ de la Gauss queda en 0.01 (Landau sola).

**2. Intervalos entre triggers con t ≥ 5.25 h.** 142 544 triggers, 142 543 intervalos, ninguno sobre 50 T/(N−1); sin corte de huecos. ⟨ΔT⟩ = 1.681544 s (1/⟨ΔT⟩ = 0.5947 Hz), λ_fit = 0.5947 ± 0.0017 Hz, KS D = 0.00265, p = 0.2688. Pendiente de ln(densidad) en [0.5, 6]: −1.0027 ± 0.0036 (χ²/ndf 47.8/53).

**3. Incertidumbre de aceptancia de ε_XY (sin cambios, solo origen).**
- **0.037 (5 p.e.) y 0.029 (10 p.e.):** son la suma en cuadratura de dos términos de `stage2.py`, ambos relativos a ε_XY.
  - Aceptancia (f_2018): 3.50 % = estadística del MC (0.32 %) ⊕ placas z ± 2 mm (2.7 %) ⊕ cota de posición de 10 mm (2.2 %). Da 0.0275 a 5 p.e. y 0.0215 a 10 p.e.
  - g: 0.03 / 0.98 = 3.06 % (diferencia con Meiga, 0.951). Da 0.0240 a 5 p.e. y 0.0188 a 10 p.e.
  - Total: 4.65 % relativo. A 5 p.e. la suma da 0.0365 [I: cálculo]. El 0.037 viene del valor antes de HV-DIAG-5 (0.0279 ⊕ 0.0244 = 0.0371) y no se actualizó; redondeado, el 0.04 del texto no cambia. A 10 p.e. da 0.0285, de ahí el 0.029.
- **Lo que NO entra:**
  - el canal 3: f_2018 y g se calculan con las 16 barras X vivas;
  - el orden de los bloques: g es el promedio de X arriba 0.9794 y Y arriba 0.9804, y su semidiferencia (0.0005) no se suma.
- **Banda gris de `eff_vs_threshold.png`:** se dibujó con ε × 0.034977, es decir, solo el término de f_2018 (3.5 %), sin g. Lo hizo `proc_plot_eff.C` con `fig_eff_mayonly.json` (`f18_rel_syst` = 0.034977). El PNG publicado coincide en SHA256 con esa salida.

**4. Figura de dos paneles `mip_spectra_paper_2p`.**
- (a) ch10, el canal de XB con G más cerca de la mediana de la estación: G = 4.564 mV contra 4.5616 mV. Mayo v8 tiempo bueno, como amp_spectra_hod2018_pe.
- (b) el panel (b) de `mip_tagged_spectra_paper` (moda 13.23 p.e., n_net 865).
- Código: `amp_spectrum_xb.py` y `mip_paper_2p.C`.

**5. g con el canal 3 muerto (`g_dead_ch3.py`).** MC del banco de stage5 (cos³θ, placas del CAD, semilla 42, primeros 3e7 rayos: 791 485 triggers, 28 109 cruzan el cuadrado de 33 mm). g = P(X ∧ Y | cuadrado) con la tira X12 (ch3, centro 9 mm, capa A) apagada [M]:
- X arriba: g16 0.9802 ± 0.0008, g sin X12 0.9545 ± 0.0012; Y arriba: 0.9805 ± 0.0008 y 0.9558 ± 0.0012. Razón a 0.98: 0.9739 y 0.9753 (promedio 0.9746).
- Con 1e8 rayos (la muestra de `bar_mc_g.json`, cuyas cuentas de 16 barras reproduce): 0.9529 y 0.9549.
- Con g = 0.9551 (promedio de los dos órdenes) y aceptancia f_2018 3.50 % ⊕ g 3.06 % = 4.65 %: ε_XY(5) = 0.8056 ± 0.0858 ± 0.0717 ± 0.0374; ε_XY(10) = 0.6296 ± 0.0248 ± 0.0155 ± 0.0293 (stat, estimador, aceptancia). Con g = 0.98: 0.7851 y 0.6136.

**6. Tiempo bueno de run B.** 71.83 h, 143 814 triggers; quedan 77 558. Fuera [M]: arranque con todo bajo el 10 % de la mediana (4.81 h, 496), dos intervalos de bias inestable de HV-DIAG-5 (4.68 h, 9725, y 12.35 h, 25 869: 17.03 h, 35 594 triggers), grupo 2 apagado (14.00 h, 30 123) y la última hora parcial (0.02 h, 43); 35.86 h en total. Magnitud vigilada: tasa de pulsos pre-trigger de cada canal (> p0 + 7.5 G, corregida por ocupación, bins de 10 min) dividida por su mediana de 00:00-05:59 locales. No hubo registro de VMon/IMon en mayo. Los bordes de los dos intervalos inestables se fijaron a mano desde el inventario de episodios, sin un umbral único: el 2026-05-29 los grupos 1 y 2 siguen en 0.53-0.93 de su tasa nocturna hasta las 22:10, después del borde de las 20:50.

**7. Capas A y B.** Índice de tira 0-15 = orden por centro (−15 … 15 mm, fig:schematic). Capa A = índices pares (centros −15 … 13), capa B = impares (−13 … 15), en X y en Y. Mayo: X A = ch15, 13, 11, 9, 7, 5, 3, 1 y X B = ch0, 2, …, 14; Y A = ch16-23 y Y B = ch24-31. El MC y el esquema ponen la capa A arriba dentro de cada bloque [S]: no hay registro de cuál está arriba, ni de qué bloque (X o Y) está arriba. El sentido absoluto del orden de tiras de hod2018 tampoco está fijado: invertirlo cambiaría A por B en los nombres, no qué canales comparten capa.

**8. Leyenda de `eff_vs_threshold`.** Regenerada con «run B» y «run B, 50 mV threshold» (`outputs/.../eff_band/runB/`), pero no copiada al proceeding. La estrella es la convención p0 + 3.5 G de `pe_calibration_v8.txt` (15.0-24.0 mV por canal, mediana 19.7 mV, T equivalente 4.43 p.e.), no 50 mV. Los 50 mV son el umbral de las 8073 coincidencias del manuscrito (`manuscript_xy.C`, ≈ 12 p.e.); el pie de `main_v5.tex` l. 223 tiene la misma confusión. Panel (a) de `mip_spectra_paper_2p`: ch10 (su histograma coincide bin a bin con `fig_amp_xb.json`; el de ch13 no). G de ch13 = 4.56155 mV, igual a la mediana, pero ch13 está en XA.

**9. Objetivos por capa de mayo para Geant4 (`may_layer_targets.py`, `layer_fit.C`).** `outputs/hod2018_may_targets/2026-11-02/may_layer_targets.csv`: 4 capas × 4 intervalos de distancia, p.e. medidos, F = 1.75 (p = 0.429) en la cabecera. La fila de control 29-54 reproduce las modas y n de `paper_table.csv`. L_eff de la mediana [M]: XA 114.5 ± 43.0 mm (χ² 14.20/2), XB 79.8 ± 28.4 (1.82/2), YA 79.1 ± 20.4 (1.85/2), YB 97.5 ± 42.6 (2.17/2). XA 37-45: moda en el borde del ajuste (4.000), no determinada.

**10. `eff_vs_threshold` con g = 0.9551.** ε, errores y estrella × 0.98/0.95512; banda 4.648 %; leyendas «run B» y «run B, workshop threshold (p0 + 3.5 G)». ε_XY(5) = 0.8056, ε_XY(10) = 0.6296; estrella 0.930 ± 0.162 en T = 4.43 p.e. Copiada a `figs/2026-11-02-revision`.

**11. Leyenda de la banda de `eff_vs_threshold` (2026-11-03).** «acceptance uncertainty (May v8)» → «acceptance uncertainty (run B)»; nada más. Clave opcional `band_name` en `proc_plot_eff.C`; con el JSON anterior el PNG sale idéntico byte a byte al publicado (956714bb…). Salida `outputs/proceeding_revision/2026-11-03_band_legend/`; commit a49d72b en `figs/2026-11-02-revision` (PNG f6b3c530…, PDF e68fa94b…).

## hod2019 contra grupo HV y estabilidad diurna de hod2018 (2026-11-03)

Seguimiento 12, solo lectura de las ntuplas XPLUS-4b de v8 y run 8. Código en `analysis/followup12_2026_11_03/`; tablas, figuras de inspección y README en
`outputs/hod_stability/2026-11-03_followup12/`. Mayo = run B; octubre = run 8 (2026-10-22 20:20 a 10-23 13:36 UTC, después de la pérdida de hod2019).

**A. hod2019 en tiempo contra el grupo HV de hod2018** (`h19_hvgroup.py`).
- **Definiciones:**
  - Impacto de hod2019: aS > 35 mV, o > 3 mV.
  - Impacto de hod2018: aS > 5·G (G de PROC).
  - R = impactos de hod2018 en triggers con hod2019 / esperado de los triggers sin hod2019. Se calcula para el mismo grupo HV y para los otros.
  - DR = R_same/R_other; bootstrap de 200 réplicas.
- **Mapas:** planilla «Voltajes Hodoscopio» (central), grupos medidos en DAQ2-10 (r4) y grupos 1-2 medidos en mayo.
- **Resultado en mayo, tiempo bueno (77 558 triggers) [M]:**

  | | DR, planilla | DR, DAQ2-10 | DR, grupos 1-2 de mayo |
  |---|---|---|---|
  | 35 mV (2086 triggers) | 1.065 ± 0.082 | 0.894 ± 0.090 | 1.211 ± 0.138 |
  | 3 mV | 0.960 ± 0.017 | 1.009 ± 0.019 | 1.055 ± 0.034 (fuera de tiempo 1.083 ± 0.046) |

  - Δk (hod2019 − hod2018) igual en el mismo grupo y en los otros: con 35 mV, mediana −1 y sd 2.88 contra 2.81 muestras.
  - Mayor pulso de hod2018: f(> 5 p.e.) 0.515 con hod2019 > 35 mV, contra 0.535 en todos los triggers.
- **Conclusión:** los impactos en tiempo de hod2019 de mayo no están asociados a pulsos de hod2018 del mismo grupo HV [M].
- **Pulsos > 35 mV:** los 2086 están en la ventana signal y ninguno fuera de tiempo, pero no van con más actividad de hod2018 (R_other 0.94 ± 0.03) [M].
- **Octubre (hod2019 muerto):**
  - 19 triggers con > 35 mV y 112 con > 3 mV de 44 176. Son triggers con mucha actividad en hod2018 (mediana del mayor pulso 10.04 contra 3.87 p.e.).
  - Sin dependencia del grupo con el mapa de la planilla: DR 0.83 ± 0.13 y 1.03 ± 0.05 [M].

**B. Estabilidad diurna de hod2018** (`fit_gain_hourly.C`, `diurnal.py`).
- **Bins:** 1 h UTC. Mayo: 66 bins desde t = 5.81 h, 34 en tiempo bueno. Octubre: 18 bins, de 17:20 a 10:36 hora local, menos de 24 h.
- **G por canal y hora:** ajuste de picos de 2-5 p.e. con p0 fijo al valor del run, sobre el espectro `spec` (método de PROC).
  - gain_fit.h por hora se descartó: acepta 777 de 2336 canal-hora.
- **DCR > 0.5 p.e. y p_ct = N(> 1.5)/N(> 0.5):** de la región pre-trigger, con el umbral de la G de cada hora.
- **Niveles de la estación (mediana de los bins buenos) [M]:**
  - mayo: G 4.712 mV, DCR 2734 kHz, p_ct 0.326 (0.264 corregido por apilamiento);
  - octubre: G 3.869 mV, DCR 2343 kHz, p_ct 0.259 (0.199 corregido).
- **Modulación de G/G_ref** (sinusoide de 24 h más deriva, errores escalados por √(χ²/ndf)) [M]:

  | | pico a pico | máximo (hora local) | χ²/ndf | ΔT equivalente [I] |
  |---|---|---|---|---|
  | mayo, tiempo bueno | 1.60 ± 0.78 % | 7.1 ± 1.5 h | 216.4/30 | 0.38 ± 0.19 °C (plegado 0.55) |
  | octubre, sin deriva | 3.06 ± 0.57 % | 11.4 ± 0.6 h (borde de los datos) | 87.4/15 | 0.74 ± 0.14 °C (rango 1.41) |

  - Variantes de mayo, pico a pico: región pre-trigger 1.28 ± 0.51 %; p0 libre 1.08 ± 0.80 %.
  - Mayo: el cambio mayor es la primera noche, del 28 al 29: G sube de 0.97 a 1.01 entre las 21 h y las 06 h, justo después del arranque inestable. En las noches 2 y 3, G queda dentro de ±1 %.
  - Octubre: G sube casi monótonamente, de 0.99 (21-01 h) a 1.04 (10 h).
  - En las dos campañas G es mínimo al anochecer y sube hacia la mañana.
  - El DCR de mayo deriva −2.18 ± 0.12 %/día a lo largo de run B; G no deriva (−0.03 ± 0.18 %/día).
- **Correlación entre bins (tiempo bueno) [M]:**
  - G–DCR: Pearson −0.554 (mayo) y −0.730 (octubre).
  - G–p_ct: +0.731 y +0.870.
  - d ln p_ct / d ln G: 0.95 ± 0.16 y 1.26 ± 0.18.
  - Por canal: mediana de Pearson G–DCR −0.51 y −0.49.
  - [I] El signo es el de un OV que cambia por V_bd(T) con V_bias fijo. Un cambio de V_bias daría G y DCR con el mismo signo.
- **ΔT equivalente:** ΔT = (ΔG/G)·OV/(dV_bd/dT) con OV = 1.2 V y 50 mV/°C [S]: derivada bajo ese supuesto. No hay registro de temperatura para comprobarla.

**C. Prueba geométrica: ¿hod2019 detecta solo cerca del extremo con MPPC?** (`h19_end_mc.py`, `h19_end_test.py`; mismo directorio, `C_h19_end/`).
- **MC del banco de stage5:** 1e8 rayos, 2 639 979 triggers. hod2019 en z = 98.4 mm, eficiencia 1 a menos de z del extremo. Extremos de lectura por capa según kCfg_2026_05: X par = Delantera, impar = Trasera; Y par = Opuesta, impar = Frontal. Se prueban las dos asignaciones de signo.
- **Probabilidades [estimado]:**

  | | z = 2 mm | 3 mm | 5 mm | 8 mm |
  |---|---|---|---|---|
  | P_X | 0.0086 | 0.0133 | 0.0235 | 0.0409 |
  | P_Y | 0.0087 | 0.0134 | 0.0236 | 0.0410 |

- **z que reproduce lo medido:** z_X = 3.2 mm y z_Y = 2.1 mm a 35 mV; 14.7 y 13.1 mm a 3 mV.
- **Coincidencias predichas:** en las z ajustadas, P(Y | X) = 0 y P(hod2018 | X) = 0 en todas las asignaciones. Las primeras X·Y aparecen con z ≥ 16.75 mm y las primeras trayectorias que cruzan hod2018 con z ≥ 8.25 mm. Predichas: 0 X·Y y 0 de dos estaciones (< 0.16 al 95 %), contra 3 y 1 observadas.
- **Datos [M]:**
  - Las 3 X·Y > 35 mV tienen X e Y en la misma muestra (kS 261, 261 y 259); el evento de dos estaciones es una de ellas (X8 ch51, Y15 ch40).
  - El Y encendido es la barra Y con el extremo más cercano al de la barra X solo en 4 de 16 orientaciones, y solo en ese evento; en los otros dos tiene rango 3-14.
  - Por canal: Delantera 1696 contra Trasera 398 en X, Opuesta 588 contra Frontal 767 en Y; ch49, ch39 y ch43 sin impactos. El MC predice capas iguales y barras casi uniformes (χ²/15 de 974 a 1924).
- **Conclusión [I]:** la hipótesis del extremo reproduce r_X y r_Y solo porque z queda libre. El déficit de X·Y lo sobre-explica: predice ninguna coincidencia de un mismo muón, mientras que hay 3 X·Y sincronizadas y 1 de dos estaciones. Tampoco reproduce la distribución por canal.

## Resolución temporal y posición en el solapamiento de hod2018, run B (2026-11-03)

Seguimiento 13. Solo lectura del crudo v8, tiempo bueno (77 558 triggers). Código en `analysis/followup13_2026_11_03/`; tablas, figuras de inspección y README en `outputs/hod2018_timing/2026-11-03_followup13/`. Nada para el artículo.

**A. Resolución temporal** (`time_extract.C`, `time_res.py`).
- **Selección:** muón etiquetado = tira X y tira Y más altas, las dos > 10 p.e. medidos (ch3 fuera). 2459 eventos; mezclados (Y del evento siguiente) 1243. Los histogramas son mismo − mezclados.
- **Tiempos:**
  - CFD 0.5: interpolación lineal en el frente de subida, con el pedestal local [k−12, k−4].
  - El frente de subida medio dura unas 2 muestras: 0.23-0.27 en k−2, 0.80-0.82 en k−1.
  - kmin = muestra del mínimo × 16 ns.
- **Correcciones:**
  - Desplazamiento por canal: rango −0.82 a +1.95 ns, RMS 0.72 ns.
  - Propagación, s·(d_X − d_Y): d se toma de la tira ortogonal y del extremo de lectura por capa. s ajustada = 0.0155 ± 0.0094 ns/mm, no significativa; σ(d_X − d_Y) = 13.8 mm. Cambia σ en menos de 0.06 ns.
- **CFD 0.5, con desplazamientos y propagación [M]:**

  | menor amplitud (p.e.) | n (mismo) | σ(Δt) ns | RMS ±3σ ns | σ por tira = σ/√2 ns |
  |---|---|---|---|---|
  | 10-15 | 1105 | 2.59 ± 0.12 | 2.68 | 1.83 ± 0.09 |
  | 15-20 | 573 | 2.17 ± 0.10 | 2.32 | 1.53 ± 0.07 |
  | 20-30 | 472 | 1.90 ± 0.10 | 2.15 | 1.34 ± 0.07 |
  | > 30 | 182 | 1.91 ± 0.15 | 1.95 | 1.35 ± 0.11 |
  | todas | 2332 | 2.27 ± 0.06 | 2.42 | 1.60 ± 0.04 |

  - Sin correcciones: σ = 2.65, 2.45, 2.23 y 1.78 ns. Con el pedestal del tramo (cfd50_bl): σ = 2.44, 2.51, 1.96 y 1.96 ns.
- **Mínimo de muestra [M]:**
  - σ gaussiana: 9.22, 8.95, 7.94 y 7.61 ns. RMS ±3σ: 8.77, 8.51, 7.33 y 6.81 ns. Por tira (RMS/√2): 6.20 a 4.81 ns.
  - La gaussiana no describe bien la distribución, que es discreta (χ² 25.1/2 en 10-15).
  - Fracción neta en Δk = 0: 0.58, 0.66, 0.77 y 0.92.
  - El CFD mejora la resolución por tira en un factor 3.4-3.9 (RMS/√2 del mínimo de muestra contra σ/√2 del CFD).
- [I] σ por tira de 1.3-1.8 ns con muestreo de 16 ns, mejor a mayor amplitud hasta unos 20 p.e. y plano por encima. El supuesto X = Y no está comprobado.

**B. Posición dentro del solapamiento** (`overlap_mc.py`, `overlap_ratio.py`).
- **Selección:** etiqueta ortogonal > 10 p.e.; exactamente dos tiras de la coordenada > 5 p.e., vecinas en el orden de posición; R = a1/(a1 + a2). X: 653 pares, mezclados 234.0. Y: 637 pares, mezclados 232.0.
- **MC del banco** (stage5, 3e7 rayos, 791 485 triggers; misma selección por camino > 0; X arriba e Y arriba):
  - L1 = L2 en el 77.5-79 % de los pares, porque las dos capas son planos separados 3 mm en z y la trayectoria atraviesa el espesor completo de las dos barras. R_path = 0.5 en esos pares; el resto son cruces por el borde.
  - Posición del cruce en el plano medio: RMS 0.40 mm; 77-79 % dentro de ±0.5 mm.
  - dR_path/du = −0.138 a −0.144 /mm, con correlación −0.50 a −0.51, causada por los cruces por el borde.
- **Anchos de R [M datos, estimado MC]:**

  | | X | Y |
  |---|---|---|
  | datos netos: media, RMS, σ_IQR | 0.502, 0.173, 0.204 | 0.544, 0.206, 0.222 |
  | datos dentro del par: RMS, σ_IQR | 0.171, 0.182 | 0.162, 0.151 |
  | MC camino: RMS, σ_IQR | 0.115-0.111, 0 | 0.110-0.113, 0 |
  | MC + fotoestadística (13 p.e., F = 1.75): RMS, σ_IQR | 0.128, 0.139-0.140 | 0.128, 0.136-0.142 |
  | MC + espectro medido: RMS, σ_IQR | 0.193, 0.198-0.199 | 0.191, 0.188-0.189 |

  - Medianas por par de R: de 0.43 a 0.60 en X y de 0.27 a 0.75 en Y. La de Y llega a 0.75 con ch30, de mediana 43 p.e. Son diferencias de respuesta entre tiras, que ensanchan la distribución global de Y.
  - La corrección de atenuación con el L_eff por capa no cambia el ancho: X 0.175, Y 0.202.
- **Resolución en posición [estimado]:**
  - σ_R de la fotoestadística sola = 0.128. Con |dR/du| ≈ 0.14 /mm, σ_u = 0.89-0.93 mm.
  - Eso es más que el RMS de la posición dentro del solapamiento (0.40 mm): el cociente no aporta posición dentro del solapamiento.
  - La resolución la da el solapamiento mismo: 1 mm/√12 = 0.29 mm si fuera uniforme.
- **Conclusión [I]:** R no sigue el cociente de caminos. El MC lo predice concentrado en 0.5, y el ancho medido (0.16-0.17 dentro del par) supera el límite fotoestadístico (0.128) y queda cerca del modelo con el espectro medido, que incluye Landau y la dispersión de ganancias. La amplitud no da posición dentro del solapamiento de 1 mm.

