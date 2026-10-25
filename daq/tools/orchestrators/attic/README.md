# Orquestadores de DAQ2-6 a DAQ2-13 (attic)

Copias (`cp -p`, mtime conservado; SHA256 idénticos) de los scripts de `/data/daq/logs/`, hechas en XPLUS-2 C (2026-10-25). Los originales siguen allí. Se les quitó el bit de ejecución: son registro, no herramientas. Tienen rutas fijas a configuraciones de campañas cerradas y dos de ellos escriben HV. No correrlos.

- `sampler.sh` — DAQ2-6, 2026-10-19: muestrea CPU y RSS del MainPID de una unidad de systemd cada N s mientras está activa (lo usan los orquestadores de abajo).
- `night_daq2-6.sh` — DAQ2-6, noche 2026-10-19/20 (timer a las 01:00Z): paso 4.1. Run de daqrun con trigger externo, ventana B, 30 min, HV encendida; si la tasa a los 3 min es < 0.1 Hz para el run y reintenta cada 30 min hasta las 06:00Z. Después, paso 4.2: SW trigger 100 Hz con registro completo, 2 min. HV solo leída (`hvctl status`).
- `run8_orchestrator.sh` — DAQ2-6, 2026-10-19: espera el fin del run 7 (unidad `daq2-6-s4a2`) y arranca el run 8 (trigger externo, ventana B, HV encendida) hasta las 19:30Z, con `sampler.sh`.
- `disk_sampler.sh` — DAQ2-7, 2026-10-19: `df` de /data y `du` del directorio de la campaña cada N s mientras una unidad está activa.
- `night_orchestrator.sh` — DAQ2-7, noche 2026-10-19/20: run nocturno de daqrun en `night_2026_10_19` (config del run 8 con ttt_tick_ns 8.000071) hasta las 13:00Z del 2026-10-20, con muestreo de RSS (60 s) y de disco (600 s); la HV no se toca.
- `phase5_orchestrator.sh` — DAQ2-9, 2026-10-20: run de prueba de la fase 5 (60 min, trigger externo, ventana B, HV encendida) en `phase5_test_2026_10_20`, con muestreo de CPU/RSS y de disco.
- `daq2-10_toggle.py` — DAQ2-10 parte C, 2026-10-21: prueba de grupos HV, un grupo apagado a la vez (escrituras de Pw por `hvctl pw --authorize`), con daqrun en SW trigger 100 Hz; mapeo `hodo_mapping_2026_10.json`; si aborta, vuelve a encender todo.
- `daq2-12_toggle.py` — DAQ2-12, 2026-10-21: segunda prueba de grupos HV con referencias largas, 240 s por grupo, corte por luz (> 1.2 × referencia) y límite a las 17:30Z; mapeo `hodo_mapping_2026_10_r3.json`.
- `daq2-13_selftrig.sh` — DAQ2-13, 2026-10-21: tres runs de daqrun de 5 min con auto-trigger de hod2018 (T = 12, 10 y 8 p.e.), uno tras otro.
