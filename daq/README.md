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

En el modo directo ningún subcomando escribe en la fuente. La única escritura del código es Pw
de un canal y la hace `hvmon` (por `hvctl pw`) solo con `allow_commands: true`.

## Estructura

```
daq/
  core/          JThread/StopToken, BoundedQueue, Snapshot, Clock, Logger, unidades, secretos,
                 archivos, lock de instancia, socket Unix
  hv/            HVBackend, CaenHVBackend (Wrapper), FakeHVBackend, Status, credenciales, JSON
  temp/          interfaz de temperatura y backend fake
  hvmon/         lógica de hvmon: sondeo, cambios, escritor, comandos, estado
  apps/hvmon/    daemon hvmon (único que enlaza ROOT: THttpServer)
  apps/hvctl/    herramienta de línea de comandos
  systemd/       unidad de usuario de hvmon (no instalada)
  tools/legacy/  hv_readall.c, test_sy5527.c, diag_sy5527.sh (heredados)
  tests/         GoogleTest, solo con fakes
  config/        ejemplos de configuración y guiones del fake
```
