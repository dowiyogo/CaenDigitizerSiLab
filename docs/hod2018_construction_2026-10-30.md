# hod2018: construcción y orientación (datos de René, 2026-10-30)

## Datos de construcción confirmados por René

- **Barras:** 16 X y 16 Y de 75 × 3 × 1 mm, dos capas por coordenada (intercaladas), un MPPC por barra.
- **Reflector:** Vikuiti (ESR) sujeto con Kapton, **sin adhesivo**: hay una capa de aire entre la barra y el Vikuiti.
- **Acoplamiento:** cada extremo de barra va pegado con **epoxi óptico** a una **guía de luz de PMMA pintada por fuera**, y la guía va pegada al MPPC.
- **Sistema de coordenadas:** origen (0,0) en la esquina inferior izquierda. Las barras X son paralelas al eje Y (miden x); las barras Y son paralelas al eje X (miden y).
- **Lado de lectura de cada capa:**
  - capa Y superior: por la derecha (x alto);
  - capa Y inferior: por la izquierda (x bajo);
  - capa X superior: por abajo (y bajo);
  - capa X inferior: por arriba (y alto).

## Orden de tiras y sentido absoluto (estado al 2026-10-30)

**El orden relativo está medido.**
- Coincidencias entre vecinas en la noche del 2026-10-19 (`outputs/daq2-10/neighbours.txt`) y en mayo v8.
- Mapeo r4 (`daq/config/mapping/hodo_mapping_2026_10_r4.json`):
  - X: índice i = 0…15 → ch 15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14 (igual a `kCfg_2026_05`);
  - Y: índice j → ch 16 + j (igual a la hoja «Octubre 2026»).
- El sentido absoluto (i ↔ 15 − i, j ↔ 15 − j) no está determinado.

**Capas.** Las tiras alternan capa: la capa A tiene los índices pares (centros −15, −11, … mm en el sistema de índices) y la B los impares.
- X: A = canales impares (15, 13, …, 1) y B = canales pares (0, 2, …, 14).
- Y: A = ch 16, 18, …, 30 y B = ch 17, 19, …, 31.

**Lo que fija el lado de lectura.** Combino el dato de René con la eficiencia a lo largo de las barras por capa de PROC-2 etapa 2 (`daq/tools/proceeding/along_bar.py`, `outputs/proceeding_2026_10/work/along_bar.json`, ε/ε_MC, 10 p.e.):
- **Barras X** (posición = tira Y de referencia j):
  - la capa A baja con j: −0.0078 ± 0.0027 /mm en mayo v8, −0.0072 ± 0.0038 /mm en octubre;
  - la capa B sube con j: +0.0060 ± 0.0025 y +0.0094 ± 0.0038 /mm.
  - Signos opuestos: las dos capas X se leen por extremos opuestos, como dice René.
  - La capa X superior se lee por abajo, así que su eficiencia baja con y físico. Por eso: **si j crece con y (sentido «+» de Y), la capa X superior es la A (canales impares); si j decrece con y, es la B (canales pares).**
- **Barras Y** (posición = tira X de referencia i):
  - capa A plana: +0.0009 ± 0.0028 (mayo), +0.0010 ± 0.0039 (octubre);
  - capa B baja con i: −0.0110 ± 0.0025 y −0.0091 ± 0.0034.
  - No hay signos opuestos, así que la prueba es débil.
  - Si se toma la tendencia de B: la capa Y superior se lee por la derecha y su eficiencia sube con x físico. Por eso: **si i crece con x, la capa B de Y (ch 17, 19, …, 31) es la inferior; si i decrece con x, es la superior.**

**Lo que falta para fijar el sentido.** Hace falta un dato que una canales concretos con un lado físico. Por ejemplo:
- en qué borde (abajo o arriba) están los MPPC de los canales impares de X;
- o qué canales salen por el conector de la derecha en Y.

La hoja «Octubre 2026» de «Disposicion Correcta Hodoscopios» no está en este host; solo su resumen (tira k → ch k) está en `analysis/src/hodo_mapping.h`. Si la hoja anota el lado de cada canal, ese es el dato.

Saber qué capa (superior o inferior) tiene la barra más cercana al origen no alcanza por sí solo. Junto con lo anterior solo da una relación entre los dos sentidos: queda libre el giro de 180° alrededor de z, que cambia a la vez el sentido de X, el de Y y qué capa es la superior.

**Coherencia con el cableado (inferido).** En X, los canales impares van en orden decreciente y los pares en orden creciente a lo largo de i. Es lo que se espera si cada capa X se cablea desde su propio extremo con la misma numeración de conector, y es coherente con la lectura por extremos opuestos. En Y las dos capas van en el mismo sentido (16…31), así que el cableado de Y sigue otra convención. Esto no fija el sentido.
