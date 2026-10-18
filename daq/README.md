# DAQ v2

Código nuevo de adquisición y control lento del hodoscopio. Diseño completo en
`docs/daq_v2_design.md`. Estado: **fase 0** (primitivas de concurrencia, biblioteca HV de solo
lectura, `hvctl`). Nada de esto toca el digitizador ni escribe en la fuente de HV.

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

## hvctl (solo lectura)

Credenciales en `~/.config/daqv2/hv_credentials.json`, con modo 600 (si no, `hvctl` se niega):

```json
{"user": "admin", "password": "..."}
```

Configuración opcional en `~/.config/daqv2/hvctl.json` (ejemplo: `daq/config/hvctl.example.json`).

```bash
H=daq/build/apps/hvctl/hvctl
$H map                                  # mapa del crate
$H params --slot 3 --ch 0               # parámetros: tipo, modo, unidad, rango
$H read --slot 3                        # tabla igual a hv_readall
$H read --slot 3 --decode               # además, bits de Status con nombre
$H read --slot 3 --json                 # JSON con UTC ns y valores en SI (V, A, V/s, s)
$H bench --slot 3 --n 100               # latencia de 100 lecturas: min, mediana, p95, max
$H --fake default read --slot 3         # lo mismo contra la fuente simulada
$H --fake daq/config/fake/trip_ch2.json read --slot 3 --decode
```

Ningún subcomando escribe en la fuente; el binario no importa `CAENHV_Set*` ni
`CAENHV_ExecComm` (`nm -D --undefined-only daq/build/apps/hvctl/hvctl`).

## Estructura

```
daq/
  core/          JThread/StopToken, BoundedQueue, Snapshot, Clock (real y manual), Logger, unidades
  hv/            HVBackend, CaenHVBackend (Wrapper), FakeHVBackend, Status, credenciales
  apps/hvctl/    herramienta de línea de comandos
  tools/legacy/  hv_readall.c, test_sy5527.c, diag_sy5527.sh (heredados)
  tests/         GoogleTest, solo con fakes
  config/        ejemplos de configuración y guiones del fake
```
