# DAQ v2

Código nuevo de adquisición y control lento del hodoscopio. Diseño completo en
`docs/daq_v2_design.md`. Estado: **fase 1** (`hvmon`: registro continuo de la HV, registro de
cambios, estado publicado, página web, ruta de comandos). Nada de esto toca el digitizador. La
única escritura posible a la fuente (Pw de un canal) está deshabilitada por configuración
(`allow_commands: false`).

## Requisitos (ya instalados en el PC del DAQ)

GCC 13, CMake ≥ 3.20, ROOT 6.32 en `~/root6` (solo para configurar; la fase 0 no lo usa),
nlohmann-json3-dev, libgtest-dev 1.14, CAEN HV Wrapper 6.6 en `/usr/lib64`.

## Compilar y probar

Desde la raíz del repositorio:

```bash
source ~/root6/bin/thisroot.sh            # o: export ROOTSYS=$HOME/root6
cmake -S daq -B daq/build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build daq/build -j4
ctest --test-dir daq/build --output-on-failure
```

Con sanitizers (cada uno en su propio directorio de build):

```bash
cmake -S daq -B daq/build-address-undefined -DDAQ_SANITIZE=address,undefined
cmake -S daq -B daq/build-thread -DDAQ_SANITIZE=thread
cmake --build daq/build-thread -j4 && ctest --test-dir daq/build-thread --output-on-failure
```

TSan corre los tests con `setarch x86_64 -R` (ASLR desactivado solo para ese proceso), porque
el TSan de GCC 13 falla con la configuración de ASLR de este kernel. Los sanitizers no funcionan
con `ulimit -v`; para limitar memoria usar
`systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 <comando>`.

## Secretos

Todos los secretos van en `~/.config/daqv2/secrets.env` (modo 600, directorio 700), nunca en el
repositorio. Plantilla: `daq/config/secrets.env.example` (`DAQ_HV_HOST`, `DAQ_HV_USER`,
`DAQ_HV_PASSWORD`, `DAQ_SMTP_*`, `DAQ_NOTIFY_*`). Los programas lo leen solos (o del entorno, p. ej.
`EnvironmentFile=` de systemd) y se niegan si el modo no es 600. El antiguo
`~/.config/daqv2/hv_credentials.json` solo se usa como respaldo, con un aviso; se puede borrar.

## hvmon (control lento continuo)

```bash
export LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH
daq/build/apps/hvmon/hvmon --config daq/config/hvmon.example.json      # Ctrl+C para parar
```

- Lee la fuente cada 1 s (V0Set, I0Set, VMon, IMon, Status, Pw) y todo cada 60 s.
- Escribe `/data/daq/slowcontrol/AAAA/hv_AAAAMMDD.csv`, `hv_full_...csv`, `hv_changes_...jsonl`,
  `notes_...jsonl` (cada uno con su `.schema.json`) y `offsets.jsonl`.
- Estado cada segundo en `/dev/shm/daqv2/hvmon_state.json`.
- Página de solo lectura en `http://127.0.0.1:8081/` (gráficos de VMon e IMon de CH00-07 de las
  últimas 24 h y tabla de estado). Desde otro equipo: `ssh -L 8081:127.0.0.1:8081
  xplus@100.114.138.29` y abrir `http://localhost:8081/`.
- Un solo `hvmon` a la vez (el segundo sale con un mensaje).
- Servicio de usuario: `daq/systemd/hvmon.service` (instrucciones dentro del archivo; requiere
  `sudo loginctl enable-linger xplus`).

## hvctl

Configuración opcional en `~/.config/daqv2/hvctl.json` (ejemplo: `daq/config/hvctl.example.json`).

Con `hvmon` corriendo, `hvctl` le pregunta a `hvmon` (no abre otra sesión con la fuente):

```bash
H=daq/build/apps/hvctl/hvctl
$H status                               # conexión, contadores, valores por canal
$H read                                 # tabla igual a hv_readall (última lectura de hvmon)
$H changes --n 20                       # últimos cambios con su origen (trip, hvmon_command, external)
$H note "abrí GECO para revisar CH03"   # nota en notes_AAAAMMDD.jsonl
$H measure-offsets --seconds 60         # offsets de IMon/VMon de los canales con Pw = 0
$H pw --ch 3 on --authorize "motivo"    # rechazado mientras allow_commands = false
```

Sin `hvmon`, modo directo de solo lectura:

```bash
H=daq/build/apps/hvctl/hvctl
$H map                                  # mapa del crate (solo sin hvmon)
$H params --slot 3 --ch 0               # parámetros: tipo, modo, unidad, rango
$H read --slot 3                        # tabla igual a hv_readall
$H read --slot 3 --decode               # además, bits de Status con nombre
$H read --slot 3 --json                 # JSON con UTC ns y valores en SI (V, A, V/s, s)
$H bench --slot 3 --n 100               # latencia de 100 lecturas: min, mediana, p95, max
$H --fake default read --slot 3         # lo mismo contra la fuente simulada
$H --fake daq/config/fake/trip_ch2.json read --slot 3 --decode
```

En el modo directo ningún subcomando escribe en la fuente. Las únicas escrituras del código son Pw
de un canal y `ClrAlarm` de la placa; las hace `hvmon` (por `hvctl pw` o el auto-retry) solo con
`allow_commands: true`.

## Auto-retry, correos y compresión (fase 2)

Detalle en `docs/daq_v2_design.md`, secciones 8.1-8.5. Todo viene deshabilitado en
`daq/config/hvmon.example.json`.

- Auto-retry: `"auto_retry": {"enabled_channels": [..]}` (exige `allow_commands: true`). Tras un
  trip espera 5, 30 y 60 min y escribe solo `Pw = 1`; si vuelve a saltar tras el tercer intento
  queda en `GAVE_UP`. Cualquier cambio de `Pw` o `V0Set` hecho a mano (GECO) pasa el canal a
  `MANUAL` y el auto-retry no lo toca hasta que esté estable 10 min. Estado en `hvctl status`,
  en la tabla "Auto-retry" de la página y en `~/.local/state/daqv2/hvmon_retry_state.json`.
  `test_mode` (esperas de 1/2/3 min) solo para canales de prueba (CH08-11).
- Correos: `"notifier": {"enabled": true}` con `DAQ_SMTP_URL`, `DAQ_SMTP_USER`,
  `DAQ_SMTP_PASSWORD`, `DAQ_NOTIFY_FROM`, `DAQ_NOTIFY_TO` (varios separados por comas) en
  `secrets.env`. Por omisión solo avisa `GAVE_UP`, errores de escritura y pérdida de conexión
  (> 5 min). Prueba: `hvctl notify-test` (un correo marcado PRUEBA).
- Compresión: al cambiar de día `hvmon` comprime el día anterior a `.gz` y borra el plano solo
  si el `.gz` releído coincide (bytes, líneas, CRC32). A mano:

```bash
S=daq/build/apps/sctool/sctool
$S compress --keep /data/daq/slowcontrol/2026/hv_20261018.csv   # --keep conserva el plano
$S cat /data/daq/slowcontrol/2026/hv_20261018.csv.gz | head       # lee plano o .gz
$S digest /data/daq/slowcontrol/2026/hv_20261018.csv.gz           # bytes, líneas, crc32
```

## daqrun (runs de adquisición; fase 3: solo digitizador simulado)

Un run necesita una campaña abierta. Con el digitizador simulado (`fake` o `replay`) todo va a
`/data/daq/sim/` (nunca a `/data/daq/raw`):

```bash
C=daq/build/apps/daqcampaign/daqcampaign
$C --sim open prueba_1          # crea /data/daq/sim/raw/prueba_1/ y su manifiesto
$C --sim status prueba_1
export LD_LIBRARY_PATH=$HOME/root6/lib:$LD_LIBRARY_PATH
daq/build/apps/daqrun/daqrun --config daq/config/daqrun.example.json --campaign prueba_1 --max-duration 600
$C --sim close prueba_1         # después del último run: la campaña queda de solo lectura
```

- Página de solo lectura en `http://127.0.0.1:8080/` (tasa de trigger, Δt, tasas por canal,
  amplitudes en mV con umbrales PROVISIONALES, colas, flags; aviso si hvmon no corre).
- Ctrl+C (o SIGTERM) cierra el segmento y el índice de forma limpia.
- Archivos: `<campaña>_<UTC>_runNNN_segMMM.root` (sin permiso de escritura al cerrarse),
  `..._runNNN.index.json`, `..._runNNN.config.json` y `SHA256SUMS` (`sha256sum -c SHA256SUMS`).
- El escenario sintético está en `daq/config/fake/dig_synthetic_default.json` (tasa, ruido,
  pulsos oscuros, guion de grupo HV apagado, fuga de luz, huecos, ráfagas). Para reproducir v8:
  `"digitizer": {"backend": "replay", "scenario": {"file": "lectura_vx1740_ext_nim_1M_v8.root",
  "first_entry": 30000, "max_entries": 5000, "speed": 1}}`.
- Para leer los archivos de control lento desde ROOT: `analysis/macros/plot_hv_imon_day.C`.


```
daq/
  core/          JThread/StopToken, BoundedQueue, Snapshot, Clock, Logger, unidades, secretos,
                 archivos, lock de instancia, socket Unix
  hv/            HVBackend, CaenHVBackend (Wrapper), FakeHVBackend, Status, credenciales, JSON
  temp/          interfaz de temperatura y backend fake
  hvmon/         lógica de hvmon: sondeo, cambios, escritor, comandos, estado, auto-retry,
                 notificador, compresor, mapeo
  dig/           digitizador: interfaz, FakeDigitizer (sintético), replay de v8, bucle de lectura
  daqrun/        configuración, campañas, cliente HV y flags, monitor, notas
  io/            escritor ROOT de los runs; lector de control lento (libdaq_screader.so)
  apps/hvmon/    daemon hvmon (único que enlaza ROOT: THttpServer)
  apps/hvctl/    herramienta de línea de comandos
  apps/sctool/   compresión y lectura de archivos de control lento (plano o .gz)
  apps/daqrun/   un run de adquisición (THttpServer :8080)
  apps/daqcampaign/  abrir, cerrar y ver campañas
  tools/validation/  macros de validación de DAQ2-5 (solo lectura)
  systemd/       unidad de usuario de hvmon (no instalada)
  tools/legacy/  hv_readall.c, test_sy5527.c, diag_sy5527.sh (heredados)
  tests/         GoogleTest, solo con fakes
  config/        ejemplos de configuración y guiones del fake
```
