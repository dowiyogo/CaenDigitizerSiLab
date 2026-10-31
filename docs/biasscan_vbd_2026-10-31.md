# Barrido de bias para medir V_bd por canal de hod2018 (preparado, no armado)

Seguimiento 2026-10-31 (9), punto 1. Rama `daq/biasscan-vbd-2026-10-31`, encima de `daq/hvmon-rup-2026-10-31`.

**Se arma solo con una decisión explícita de René.** Este documento no autoriza nada: lanzarlo escribe V0Set en la SY5527 y usa el digitalizador, así que la tarea que lo lance tiene que autorizar las dos cosas (CLAUDE.md, reglas 1 y 2).

## Por qué

Con las tomas existentes (2026-10-23), G de hod2018 solo se mide en 2 voltajes, nominal y −0.5 V, con el método de c1; a −1.0 V los pulsos oscuros desaparecen.
- La estimación de dos puntos da un OV nominal de ~1.24 V (seguimiento 8).
- Este barrido agrega escalones donde hay pulsos: +0.15, 0, −0.25 y −0.4 V. Con eso hay 3-4 puntos por canal para G = k (V − V_bd).

## Escalones

Todos los grupos van con los mismos ΔV de V0Set respecto del nominal (65.75, 65.75, 65.90, 65.85, 65.75, 65.70, 66.15, 65.90 V), salvo el positivo del grupo 4.
- RUp sin cambiar: 1 V/s, el actual de CH00-07 según el slow control. No se usa `--set-rup`.
- I0Set sin cambiar: 15 µA.

**Corriente esperada** con el modelo I-V (b) del seguimiento 8: I = I_0 + a (V − V_bd)^n con piso de fuga, en `outputs/hod2018_iv/2026-10-31_iv/iv_fit.json`. Salida de `vbd/vbd_plan.py`:

| Grupo | Modelo | I nominal | I +0.15 V | I +0.10 V | I −0.25 V | I −0.4 V | Margen a +0.15 V | Escalón positivo |
|---|---|---|---|---|---|---|---|---|
| 1 | b, V_bd libre | 4.38 | 6.32 | — | 2.33 | 1.61 | 8.68 | +0.15 V |
| 2 | b, V_bd libre | 3.17 | 4.46 | — | 1.86 | 1.42 | 10.54 | +0.15 V |
| 3 | b, V_bd libre | 6.44 | 9.24 | — | 3.29 | 2.11 | 5.76 | +0.15 V |
| 4 | cota por el trip de CH03 | 7.42 | ≥ 11.00 | ≥ 9.69 | 3.57 | 2.21 | ≤ 4.00 | **+0.10 V** |
| 5 | b, V_bd libre | 6.98 | 9.87 | — | 3.59 | 2.27 | 5.13 | +0.15 V |
| 6 | b, V_bd fijo (dos puntos) | 4.59 | 6.34 | — | 2.59 | 1.83 | 8.66 | +0.15 V |
| 7 | b, V_bd libre | 5.17 | 7.25 | — | 2.92 | 2.10 | 7.75 | +0.15 V |
| 8 | b, V_bd libre | 6.73 | 9.34 | — | 3.62 | 2.38 | 5.66 | +0.15 V |

Corrientes en µA.

**Grupo 4 (CH03):** el ajuste I-V no quedó determinado.
- La única información sobre su corriente encima del nominal es el trip del 2026-10-30: a V0Set +0.25 V superó I0Set y se quedó limitado en 15.45 µA.
- Eso da una cota inferior de la corriente (≥ 11.0 µA a +0.15 V, con V_bd = 64.096 V, la mediana de los ajustes libres, y n ≥ 4.09); el margen a +0.15 V es como mucho 4.0 µA, sin cota inferior.
- Por precaución va a +0.10 V (`--group-pos-cap 4:0.10`). Su punto positivo queda en +0.10 V.

**Guarda de corriente de biasscan:** > 12 µA durante 30 s en algún grupo impide seguir al escalón siguiente, y ≥ 14 µA aborta. Ningún grupo la alcanza según el modelo; el grupo 4 a +0.10 V tiene una cota inferior de 9.7 µA.

## Eventos y duración

**Ajuste de c1 a −0.4 V:** no converge con ningún número de eventos.
- Prueba hecha con el espectro medido a −0.5 V, 12 000 eventos y 6.7 pulsos > 1.5 mV por evento: con 1 a 20 veces los eventos, el ajuste acepta 2-8 de 31 canales, sin tendencia (`outputs/hod2018_vbd/2026-10-31_scanplan/events_*.csv`, `vbd/vbd_events.C`).
- Con más eventos los fallos pasan de «insufficient» a «rejected»: el peine de picos con G ≈ 2.6 mV no pasa la prueba de espaciado del 20 %.
- Con otras ventanas y tolerancias (2-14 mV o 4-20 mV, tolerancia 0.25 o 0.3) tampoco mejora.

**Método de respaldo, por período:** `vbd/period_gain.py` mide el período del peine de picos con un periodograma del espectro.
- Converge a −0.5 V en 31 de 31 canales con los mismos 12 000 eventos, con un error relativo mediano de 2.5 %.
- A nominal lee 3.5-3.9 % más que c1. Un factor común no cambia V_bd si se usa el mismo método en todos los escalones.
- Como −0.4 V tiene más OV que −0.5 V, ahí debería ir igual o mejor [inferido].

**Eventos por escalón:** 12 000 a +0.15 y 0 V (donde c1 converge en 28-31 canales), 24 000 a −0.25 y −0.4 V (para bajar el error del período a ~1.8 %).

| Escalón | Eventos SW | Toma | Escalón completo |
|---|---|---|---|
| +0.15 V (grupo 4 +0.10 V) | 12 000 | 130 s | 310 s |
| 0 | 12 000 | 130 s | 310 s |
| −0.25 V | 24 000 | 250 s | 430 s |
| −0.4 V | 24 000 | 250 s | 430 s |

- Escalón completo = 120 s de estabilización + toma + ~60 s de escrituras.
- Total de escalones: 1480 s (24.7 min).
- Con la detención del run de física, `bias begin`, la restauración desde −0.4 V y el rearranque de la física: **~35-40 min** [inferido].
- Datos crudos: ~150 MB por cada 12 000 eventos SW (como las tomas del 2026-10-23), unos 900 MB en total.

## Qué revisar antes

1. **Decisión de René por escrito.** La tarea que lo lance autoriza el digitalizador y las escrituras de V0Set.
2. **Que no corra otra adquisición** (`pgrep -af -i 'caen|daq|acq'`). `biasscan.py` detiene él mismo el run de física (`daq-physics`) y lo rearranca al final; comprobar después que quedó vivo.
3. **hvmon activo y conectado:** `hvctl status`, `allow_commands: true` (lo está desde DAQ2-10), sin sesión de barrido abierta (`hvctl bias status`).
4. **V0Set en el nominal y RUp = 1 V/s** en CH00-07. El barrido lo comprueba con `--rup 1.0`; sin `--set-rup` no escribe RUp. Las escrituras de RUp de la rama `daq/hvmon-rup-2026-10-31` no se necesitan aquí.
5. **hod2019 (resuelto o no):** los grupos HV mezclan las dos estaciones, así que el barrido también cambia la tensión de los MPPC de hod2019. Sin la lectura de hod2019 (perdida desde el 2026-10-21), el barrido solo mide hod2018. Si la lectura se arregló el 2026-11-02, las mismas tomas dan también los espectros de hod2019, pero su ajuste de G no converge con el método de c1 (seguimiento 6).
6. **Temperatura de la placa A7030** estable: anotarla. La resolución es de 1 °C; con 54 mV/°C [supuesto], 1 °C mueve V_bd ~0.05 V.
7. **CAENGECO2020** puede estar abierto: no tocarlo.
8. **Copia fija:** para no depender del árbol de trabajo, instalar la copia fija de este commit con `daq/tools/biasscan/install_pinned.sh` (como el barrido del 2026-10-30) y lanzar desde ahí.

## Cómo lanzarlo (cuando René lo decida)

Separado de la terminal (regla 11), con la restauración de seguridad en ExecStopPost:

    CODE=<copia fija o raíz del repositorio>
    systemd-run --user --unit=biasscan-vbd --collect --property=RuntimeMaxSec=7200 --property=TimeoutStopSec=900 \
      --property=KillSignal=SIGTERM --property=ExecStopPost=$CODE/daq/tools/biasscan/biasscan_safety.sh \
      /usr/bin/python3 -I $CODE/daq/tools/biasscan/biasscan.py --campaign biasscan_vbd_AAAA_MM_DD --max-hours 2 \
        --steps 0.15,0,-0.25,-0.4 --max-above 0.15 --max-step 0.25 --rup 1.0 \
        --sw-only --sw-events 0.15:12000,0:12000,-0.25:24000,-0.4:24000 --group-pos-cap 4:0.10 \
        --authorize "<texto de la autorización de René>"

- Seguimiento: `systemctl --user status biasscan-vbd`, `journalctl --user -u biasscan-vbd`, `/data/daq/logs/biasscan_<campaña>.log`.
- Resultados: `outputs/daq2-15/<campaña>/biasscan_results.json`.
- `--max-above 0.15` y `--max-step 0.25`: hvmon exige escalones ≤ 0.25 V sobre el nominal.
- El orden +0.15 → 0 → −0.25 → −0.4 V baja siempre al final; la restauración sube de −0.4 V al nominal en escalones de 0.25 V.

## Análisis después del barrido

    python3 -I daq/tools/biasscan/vbd/vbd_analysis.py --out outputs/hod2018_vbd/<fecha>_vbd_scan \
      --biasscan-results outputs/daq2-15/<campaña>/biasscan_results.json

**Por toma:**
- espectro de pulsos oscuros con `vbd_gain.C` (reglas de c1);
- G con c1, donde converge, y con el período;
- VMon por grupo del slow control de 1 s.

**Por canal y método:** ajuste G = k (V − V_bd) ponderado, con 3-4 escalones; V_bd, OV nominal y k con errores (escalados por √(χ²/ndf) si es > 1). Además, resumen por grupo y por estación.

**Prueba hecha** sobre las tomas del 2026-10-23 (tres nominales y −0.5 V): con `--unweighted` reproduce la estimación de dos puntos del seguimiento 8 en los 5 canales (`outputs/hod2018_vbd/2026-10-31_scanplan/test_2026-10-23/`).

**Herramientas** en `daq/tools/biasscan/vbd/`:
- `vbd_plan.py`: corriente esperada y escalones;
- `vbd_events.C`: eventos contra convergencia;
- `period_gain.py`: G por período;
- `vbd_gain.C`: espectro y ajuste de c1;
- `vbd_analysis.py`: análisis.
