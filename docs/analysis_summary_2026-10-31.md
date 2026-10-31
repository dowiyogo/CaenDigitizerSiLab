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

## Pendientes que siguen abiertos

- **V_bd por canal con 3 o más voltajes:** necesita una toma nueva. El barrido está preparado, sin armar, en `daq/biasscan-vbd-2026-10-31` (`docs/biasscan_vbd_2026-10-31.md` en esa rama). Se arma solo con una decisión de René (seguimiento 9).
- **Causa de la pérdida de hod2019:** la inspección del 2026-11-02 (checklist e.0).
- **Sentido absoluto del orden de tiras de hod2018:** sin fijar. La prueba de orientación está preparada y no ejecutada (seguimiento 1, §5).
