# DAQ v2 con control lento de HV: documento de diseño

Estado: borrador para revisión (tarea DAQ2-0, 2026-10-18). No hay código implementado.
Las afirmaciones llevan [V: fuente] si se verificaron en este host, [I: razón] si son
inferidas, y "no medido" si nadie las midió.

## 1. Punto de partida

### 1.1 DAQ actual (CaenDigitizerSiLab)

- Conexión: `CAEN_DGTZ_OpenDigitizer(CAEN_DGTZ_USB_A4818, 62578, 0, 0x0, &handle)`. La
  clave `connectionType 2` del `.conf` se traduce a `CAEN_DGTZ_USB_A4818` (= 5 en el enum de
  CAENDigitizer 2.17.3) [V: CaenDigitizerSiLab.h:211-213, CAENDigitizerType.h:281,
  config_vx1740_ext_nim_1000000.conf]. El enlace es A4818 USB3 → CONET2 óptico → VX1740.
- Configuración (CaenDigitizerSiLab.cxx:64-160, CaenDigitizerSiLab.h:328-394):
  `Reset`, `SetRecordLength`, `SetGroupEnableMask(0xFF)`, por grupo `SetGroupDCOffset`
  (0x7FFF con polarization 0), `SetGroupTriggerThreshold`, `SetTriggerPolarity`;
  `SetAcquisitionMode(SW_CONTROLLED)`, `SetMaxNumEventsBLT(1)`. El post-trigger del x740 no
  se fija con `SetPostTriggerSize` sino escribiendo el registro 0x8114 con
  `ptriggersize*N/100 - 10` (CaenDigitizerSiLab.h:372-376). Con N = 624 y 60 % el valor
  calculado es 364 [I: aritmética entera del código; el valor escrito no está en los logs].
- Disparo externo (CaenDigitizerSiLab.cxx:519-546): self-trigger de grupos deshabilitado,
  `SetExtTriggerInputMode(ACQ_ONLY)`, `SetIOLevel(NIM)`, SW trigger deshabilitado. `init()`
  escribe antes 0x811C = 1 (TTL), pero `SetIOLevel(NIM)` lo sobrescribe después.
- Lectura (CaenDigitizerSiLab.cxx:573-656): un solo hilo; bucle que llama `ReadData`
  (MBLT) sin espera ni sleep, imprime el contador en cada vuelta, decodifica y llena el
  TTree. No hay hilos, ni cola, ni manejo de errores de `ReadData`. El `TriggerTimeTag`
  solo se imprime para los 3 primeros eventos y no se guarda. El `EventCounter` del
  hardware no se guarda.
- Salida: TTree `data_dig` con `Ch0..Ch63` (`vector<float>`, mV), `time` (`vector<float>`,
  ns), `event` (int), `timestamp` (int64, µs de `steady_clock`, una vez por lectura)
  [V: CaenDigitizerSiLab.cxx:181-205]. Se escribe en el directorio actual con el nombre del
  segundo argumento. El TTree se escribe con `Write` al final; si el proceso muere se
  pierde lo posterior al último AutoSave (caso v10).
- v8: 143814 entradas, 3692974698 bytes en disco, factor de compresión 6.35, AutoSave
  cada 1176 entradas [V: ROOT, lectura READ].
- Bibliotecas: `test_acquire_ext_nim` enlaza `libCAENDigitizer.so` (que trae CAENComm y
  CAENVME), libCore, libRIO, libTree; no carga `libCaenDigitizerSiLab.so` en ejecución
  [V: readelf -d, ldd]. La copia de `/usr/local/lib/libCaenDigitizerSiLab.so*` (mayo 2025)
  son tres archivos regulares iguales, no symlinks [V: ls -la].
- Placa: VX1740, ROC FPGA 04.29 build 8716, AMC FPGA 00.13 build 8A15
  [V: logs/commissioning_2026-05/*.log].

### 1.2 Ventana reducida: ¿es posible en el VX1740?

- Granularidad de la longitud de registro: pedidos de 1024 y 625 muestras dieron 1023 y
  624 [V: logs/commissioning_2026-05, entrada 0 de v6]. Ambos son el múltiplo de 3 inferior;
  la regla "múltiplo de 3" es [I] y debe confirmarse en la fase 4 con `GetRecordLength`.
  Valores candidatos entre 128 y 160: 129, 132, 135, 138, 141, 144, 147, 150, 153, 156, 159.
- Post-trigger: con 0x8114 = 364 el pulso de muón cae en 4128-4272 ns = muestras
  258-267, y 624 - 364 = 260 muestras de pre-trigger = 4160 ns. El registro se comporta
  como "muestras post-trigger" con paso de 1 muestra [I: coincidencia numérica; el manual
  V1740 no está en el repositorio]. El "-10" del código es una corrección empírica.
- Conclusión: una ventana de ~144 muestras es configurable (record length 144, 0x8114
  ≈ 106 para el diseño A de abajo). Debe verificarse en la fase 4 midiendo la posición del
  pico de muón (aceptación: ±1 muestra respecto a lo esperado).
- Conflicto a resolver: el pickup de hod2018 en 3510-3525 ns está ~38 muestras (603-618 ns)
  antes del inicio del pulso de muón. "~40 muestras de pedestal + ≥144 ns de control antes
  de la señal" (= 49 muestras) mete el pickup dentro del pedestal. Dos diseños:
  - A (propuesto por defecto): N = 144, inicio del registro en 3552 ns del registro
    actual. Pre-señal 36 muestras = pedestal 27 [0,27) + control 9 [27,36); señal 9
    [36,45); cola 44 [45,89) = 704 ns; pedestal posterior 55 [89,144). Pedestal total 82
    muestras, sin pickup.
  - B: N = 144 con 49 muestras de pre-señal (inicio 3344 ns); el pickup cae en las
    muestras ~10-11 y se excluye del pedestal en el análisis.
  La disposición de la ventana se declara en la configuración (rangos de muestras) y se
  valida al arrancar; el operador elige (ver Preguntas abiertas).
- Volumen con N = 144 en ADC crudo de 16 bit: 64 × 144 × 2 = 18432 B por evento; a
  0.6 Hz, 0.955 GB/día sin comprimir [I: cálculo]. `/data` tiene 419 GB libres [V: df].

### 1.3 Lo que se reutiliza del análisis

- `analysis/src/hodo_mapping.h`: `kCfg_2026_05` (v1-v10; hod2018 por correlaciones de
  vecinos, hod2019 por la hoja "Mayo 2026") y `kCfg_2026_10` (hoja "Octubre 2026",
  secuencial). Las caras de hod2019 por canal de digitizador están solo en un comentario
  (líneas 51-53): Delantera 48-55, Trasera 56-63, Frontal 40-47, Opuesta 32-39. Problema:
  `ch_to_plane()` y `ch_to_strip()` fijan 19X = 48-63 y 19Y = 32-47, válido solo para
  `kCfg_2026_05`; con `kCfg_2026_10` (19X = 32-47) dan resultados erróneos.
- Grupos de HV conocidos por los datos, no por un registro de cableado
  [V: analysis/macros/hv_diag3_pretrip.C:92-94, runs/good_intervals_2026.csv]:
  grupo 1 = {0, 22, 31, 43, 47, 54, 58, 63}, grupo 2 = {6, 20, 46, 50, 56, 59, 60, 62},
  grupo 6 = {1, 16, 17, 34, 53, 55, 57} (7 de 8). Grupos 3, 4, 5, 7, 8: no determinados.
- `analysis/src/hodo_consts.h`: `kDt_ns = 16.0`, ventana de señal 4128-4272 ns, ventana de
  muón 4000-4400 ns, umbrales de diagnóstico 50 mV (hod2018) y 35 mV (hod2019). En la ventana
  reducida estas constantes en ns absolutos dejan de valer: el DAQ guarda en la cabecera
  la posición de la ventana en "ns del registro antiguo" para poder traducirlas.
- `outputs/pe_calibration_v8.txt`: columnas `channel plane p0 G s0 s1 chi2 ndf status`,
  63 filas (falta ch3, muerto). hod2018: 31 filas, status 0 en todas, G entre 2.98524 y
  5.50853 mV/p.e.; hod2019: 13 filas con status 0 (varias con valores en los límites,
  p. ej. p0 = 0.5) y 19 con status 1, 2 o 4. Sirve como calibración inicial versionada solo
  para hod2018, como exigen los requisitos.
- `analysis/src/hodo_paths.h`: `resolveRawFile()` busca en `<raw>/<campaña>/` y devuelve ""
  si el nombre existe en dos campañas. Los nombres nuevos son únicos entre campañas, así
  que funciona sin cambios con `HODO_RAW_DIR=/data/daq/raw`.

### 1.4 Toolchain verificado

GCC 13.3.0, CMake 3.28.3, ROOT 6.32.10 compilado con C++17 (`root-config --cxxstandard` =
17), `libRHTTP.so` presente (THttpServer), nlohmann-json3-dev 3.11.3, curl 8.5.0. No hay
GoogleTest, Catch2 ni doctest instalados. Un TU C++20 que incluye `TFile.h` y
`THttpServer.h` y usa `std::jthread` compila y enlaza, con el aviso de `RConfigure.h:30`
("The C++ standard in this build does not match ROOT configuration (201703L)").

## 2. Arquitectura

Dos ejecutables (`hvmon`, `daqrun`) y una herramienta de línea de comandos (`hvctl`),
sobre una biblioteca común `libdaqcore` (CMake, directorio `daq/`).

    daq/
      CMakeLists.txt
      core/        tiempo (Clock), cola acotada, StopSource/JThread, config JSON, mapeo, flags
      hv/          IHVBackend, CaenHVBackend (Wrapper), FakeHVBackend, máquina de auto-retry
      dig/         IDigitizer, CaenDigitizer (CAENDigitizer), FakeDigitizer (replay/sintético)
      temp/        ITemperatureSource, Ds18b20Sysfs, SerialSource, FakeTemperature
      io/          escritor ROOT, CSV/JSONL de slow control, SHA256, campañas
      apps/        hvmon, daqrun, hvctl, daqnote, daqcampaign
      tests/
      config/      ejemplos de configuración y mapeos versionados

### 2.1 Procesos e hilos

| Proceso / hilo | Posee | Consume | Produce |
|---|---|---|---|
| `hvmon` hilo HV | la sesión del HV Wrapper (único en el host, con `flock`) | Wrapper (1 Hz rápido, 1/min lento) | `HvSnapshot` con número de secuencia; eventos de cambio |
| `hvmon` hilo lógica | máquina de auto-retry, detector de cambios, estimador de corriente | snapshots | comandos `Pw=1` (cola hacia el hilo HV), eventos, alertas |
| `hvmon` hilo escritor | archivos de `/data/daq/slowcontrol/` | snapshots, eventos, temperaturas | CSV/JSONL diarios con fsync |
| `hvmon` hilo temperatura | dispositivo de temperatura | sensor (10-60 s) | lecturas de temperatura |
| `hvmon` hilo notificador | procesos `curl` hijos | alertas | correos (un proceso por correo, con timeout) |
| `hvmon` hilo principal | THttpServer :8081, histogramas HV, `/dev/shm/hvmon_state.json` | snapshots | página web de solo lectura, estado publicado, señales |
| `daqrun` hilo lectura | handle CAENDigitizer (único llamador) | digitizador | bloques de eventos → cola de escritura y copia → cola de monitor |
| `daqrun` hilo escritor | `TFile`/`TTree` del run | cola de eventos, HV, notas, resúmenes | segmentos ROOT, SHA256SUMS |
| `daqrun` hilo monitor | datos planos (sin ROOT) | copia de eventos (cola con descarte) | resúmenes por minuto, entradas de ajustes de ganancia, flags de datos (bits 6, 7) |
| `daqrun` hilo cliente HV | último estado HV recibido | `/dev/shm/hvmon_state.json` cada 1 s | snapshots HV con edad; marca "stale" |
| `daqrun` hilo principal | THttpServer :8080 y todos sus histogramas/gráficos | resúmenes del monitor | página web de solo lectura; `ProcessRequests()` a 10-20 Hz; señales |
| `hvctl` | nada mientras `hvmon` corre | estado de `hvmon` (lectura) o Wrapper directo si `hvmon` no corre | volcado legible / JSON |

`hvmon` es un proceso aparte porque los trips ocurren entre runs y de noche; reiniciar o
caerse `daqrun` no debe parar el registro de HV, el auto-retry ni las alertas.

### 2.2 Primitivas de concurrencia

- C++17 en todo el código (ver Diferencias, punto 1): `std::mutex` + `std::scoped_lock`,
  `std::condition_variable`, `std::atomic`. `daq::StopSource`/`StopToken`/`JThread` propios
  (~60 líneas, misma interfaz que C++20 para migrar con un `using`).
- `BoundedQueue<T>`: capacidad fija, `push` bloqueante (cola de escritura: nunca se
  descartan eventos; si se llena, la lectura espera y el digitizador acumula en su memoria
  interna) o `try_push` con descarte del más antiguo y contador (cola del monitor),
  `close()` que despierta a los consumidores y hace que `pop` devuelva vacío al vaciarse.
- `Latest<T>`: valor más reciente protegido por mutex con número de secuencia.
- Regla: cada objeto ROOT tiene un solo hilo dueño. `ROOT::EnableThreadSafety()` al
  arrancar y `TH1::AddDirectory(kFALSE)`. THttpServer con `SetTimer(0, kTRUE)` y
  `SetReadOnly(kTRUE)`: solo atiende peticiones cuando el hilo principal llama
  `ProcessRequests()`, de modo que nunca lee un histograma mientras otro hilo lo modifica.
- Sin TBB, sin colas lock-free, sin ZeroMQ: ~0.6 Hz de eventos de 18 kB y 1 Hz de HV.
- Sin diccionarios ROOT propios: los TTree usan tipos básicos, arreglos fijos y
  `std::string`/`std::vector` de tipos básicos. No hace falta `rootcling`.

### 2.3 Tiempo

- `IClock` con `utc_ns()` (`CLOCK_REALTIME`), `mono_ns()` (`CLOCK_MONOTONIC`) y
  `sleep_until()`; `FakeClock` avanza a mano para tests acelerados.
- Toda marca en archivos es UTC ns (`int64`). Hora local nunca (Chile alterna UTC-4/UTC-3;
  v1-v10 fueron a UTC-4). Las marcas del SY5527 se ignoran (reloj ~81 min adelantado).
- Cada lectura HV guarda `t_req_utc_ns` (antes de la llamada) y `t_resp_utc_ns` (después):
  el valor leído corresponde a algún instante de ese intervalo.
- Digitizador: el `TriggerTimeTag` se reinicia al arrancar la adquisición. Se guarda el par
  (UTC antes de `SWStartAcquisition`, UTC después) como origen. Según el manual del V1740
  el TTT es de 31 bit con paso de 8 ns, vuelta cada ~17.2 s [I: no verificado en este host;
  `ttt_tick_ns` va en la configuración y se mide en la fase 4 ajustando ΔTTT contra ΔUTC].
  La vuelta se desenrolla con la hora del PC de cada lectura, que es inequívoca mientras
  el error del PC sea menor que medio período.

## 3. Flujo de datos y apagado

### 3.1 Flujo

    SY5527 ──Wrapper──> hvmon[HV] ──snapshot──> hvmon[lógica] ──> hvmon[escritor] ──> /data/daq/slowcontrol/
                             ^                    │  └──> notificador ──> curl SMTP
                             └──── Pw=1 ──────────┘
                                                  └──> hvmon[principal] ──> /dev/shm/hvmon_state.json, :8081

    VX1740 ──A4818──> daqrun[lectura] ──bloques──> cola escritura ──> daqrun[escritor] ──> /data/daq/raw/<campaña>/
                                     └──copia──> cola monitor ──> daqrun[monitor] ──resúmenes──> escritor y principal (:8080)
    /dev/shm/hvmon_state.json ──> daqrun[cliente HV] ──> escritor (árbol hv, flags) y principal
    daqnote ──O_APPEND──> /data/daq/slowcontrol/notes_YYYYMMDD.jsonl ──> hvmon y daqrun (leen nuevas líneas)

- Estado publicado por `hvmon`: JSON reemplazado atómicamente (escribir a temporal en el
  mismo directorio + `rename`) cada segundo, con `seq`, `t_resp_utc_ns`, `alive`, valores
  por canal, estado de auto-retry, y un anillo con los últimos 256 eventos de cambio, cada
  uno con id creciente. `daqrun` deduplica por id: si se pierde una lectura no pierde
  cambios. Para los comandos de `hvctl` (fase 2+) se usa un socket Unix
  `/run/user/1000/hvmon.sock` con modo 0600.
- Línea de retardo en el escritor de `daqrun`: los eventos esperan
  `flag_delay_s` (por defecto 5 s, ~3 eventos) antes de escribirse, para que los flags que
  dependen de lecturas HV posteriores al evento (bit 5) se apliquen sin reescribir el
  archivo. Los flags se pueden recalcular offline con los archivos de slow control.
- Lectura del digitizador: `ReadData` en bucle con espera de 5 ms cuando no hay datos
  (el actual no espera nunca), `SetMaxNumEventsBLT` > 1, y lectura del registro de estado
  de adquisición para detectar memoria llena. El contador de eventos del hardware se guarda
  para detectar eventos perdidos.
- Segmentos: un run se divide en archivos de como máximo `segment_max_s` (3600 s) o
  `segment_max_bytes`. Cada segmento se cierra, se le quita permiso de escritura y se añade
  su SHA256 a `SHA256SUMS` de la campaña (con `flock`). La adquisición no se detiene al
  cambiar de segmento. AutoSave cada `autosave_s` (300 s). Un fallo pierde como máximo
  5 min de un segmento, no el run entero.

### 3.2 Secuencia de apagado de `daqrun`

1. SIGINT/SIGTERM: el manejador solo pone un `std::atomic<bool>` (async-signal-safe).
2. El hilo principal lo ve y pide parada al hilo de lectura.
3. Lectura: `SWStopAcquisition`, vacía lo que quede con `ReadData` hasta 0 bytes, empuja
   los últimos bloques, `CloseDigitizer`, cierra la cola de escritura y la del monitor.
4. Escritor: vacía la cola, vacía la línea de retardo, escribe HV/notas/resúmenes
   pendientes y la cabecera final (fin UTC, eventos, tiempo vivo), `Write` + `Close`,
   `fsync`, `chmod a-w`, SHA256 a `SHA256SUMS`.
5. Monitor y cliente HV terminan al cerrarse sus colas o al pedirles parada.
6. Principal: borra THttpServer, une todos los hilos, código de salida 0 si el cierre fue
   limpio.
7. Una segunda señal no salta el cierre del archivo; solo acorta los tiempos de espera.
   Si el escritor no termina en 120 s se registra el error y se sale con código ≠ 0.

`hvmon` ante SIGTERM: termina la llamada Wrapper en curso, no inicia comandos nuevos,
escribe y sincroniza los archivos, guarda el estado de la máquina de auto-retry en
`/data/daq/slowcontrol/hvmon_state_persist.json`, publica `alive: false` y sale.

## 4. Formato de archivos

### 4.1 Archivo ROOT de `daqrun`

Nombre: `<campaña>_<YYYYMMDDThhmmssZ>_run<NNN>_seg<MMM>.root` en
`/data/daq/raw/<campaña>/`. `NNN` y `MMM` empiezan en 001.

TTree `events` (una entrada por trigger):

| rama | tipo | unidad / contenido |
|---|---|---|
| `evt_seq` | `ULong64_t` | índice secuencial dentro del run |
| `evt_counter_hw` | `UInt_t` | contador de eventos de la placa (24 bit) |
| `ttt_raw` | `UInt_t` | TriggerTimeTag tal como llega |
| `ttt_ns` | `ULong64_t` | TTT desenrollado, ns desde el arranque de adquisición |
| `t_utc_ns` | `Long64_t` | UTC estimado del trigger = origen + `ttt_ns` |
| `t_read_utc_ns` | `Long64_t` | UTC del PC al leer el bloque |
| `t_read_mono_ns` | `Long64_t` | `CLOCK_MONOTONIC` al leer el bloque |
| `dt_prev_ns` | `Long64_t` | `ttt_ns` menos el del evento anterior |
| `adc` | `UShort_t[64][N]` | muestras ADC crudas (N fijo por run, en la cabecera) |
| `ch_flags` | `UShort_t[64]` | flags por canal (sección 6) |
| `evt_flags` | `UInt_t` | flags de evento (sección 6) |
| `hv_seq` | `ULong64_t` | `seq` del snapshot HV usado para los flags |

La conversión a mV (`adc2mV`, DC offset por grupo) va en la cabecera. Lector de
compatibilidad en `analysis/` para producir `Ch*`/`time` si un programa viejo lo necesita.

TTree `hv` (1 Hz): `seq`, `t_req_utc_ns`, `t_resp_utc_ns`, `ok`, y por canal HV (12):
`v0set_V`, `i0set_A`, `vmon_V`, `imon_A`, `status`, `pw`, `retry_state`; por canal de
digitizador (64): `i_mppc_est_A` (título "ESTIMADO"), `i_est_mode`.

TTree `hv_changes`: `t_detect_utc_ns`, `t_prev_read_utc_ns`, `hv_ch`, `param`, `old_value`,
`new_value`, `origin` (`daq_auto`, `trip`, `external`), `change_id`.

TTree `notes`: `t_utc_ns`, `user`, `host`, `text`.

TTree `monitor` (una entrada por minuto): tasa por canal fuera de tiempo > 5 p.e. y en
tiempo (Hz), canales vivos por plano, tasa de trigger, tasa por cara física, ocupación de
colas, eventos descartados por el monitor.

TTree `gain` (una entrada por canal y hora, solo hod2018): `t_start_utc_ns`,
`t_end_utc_ns`, `ch`, `p0_mV`, `G_mV_per_pe`, `G_err`, `n_pulses`, `chi2`, `ndf`, `status`.

TTree `temperature`: `t_utc_ns`, `sensor_id`, `T_degC` (la única magnitud fuera de SI, con
la unidad en el nombre).

Cabecera (objetos `TObjString`/`TNamed` en el directorio raíz): `config_json` (configuración
completa ya resuelta), `mapping_json` y `mapping_version`, `calibration_json` y
`calibration_version`, `run_info_json` (campaña, run, segmento, host, commit git del DAQ y
si estaba sucio, versiones de CAENDigitizer/Wrapper/ROOT, modelo, número de serie y
firmware de la placa, par de origen TTT/UTC, inicio/fin UTC, eventos, tiempo vivo,
`record_length`, registro de post-trigger escrito y leído, disposición de la ventana).

### 4.2 Slow control (`/data/daq/slowcontrol/`)

Archivos por día UTC, texto plano con `fsync` cada 60 s y al cerrar. La primera línea de
cada CSV es `# schema=<n> units=...`.

| archivo | contenido |
|---|---|
| `hv_fast_YYYYMMDD.csv` | 1 Hz: `seq,t_req_utc_ns,t_resp_utc_ns,ok` + por canal 0-11 `V0Set_V,I0Set_A,VMon_V,IMon_A,Status,Pw` |
| `hv_slow_YYYYMMDD.csv` | 1/min: todos los parámetros de canal (lista de `GetChParamInfo`) y de la placa |
| `hv_changes_YYYYMMDD.jsonl` | un cambio por línea: id, tiempos, canal, parámetro, viejo, nuevo, origen |
| `hvmon_events_YYYYMMDD.jsonl` | transiciones de auto-retry, comandos enviados y su resultado, pérdidas de conexión, GECO detectado, correos |
| `hv_offsets.jsonl` | offsets de IMon/VMon medidos con Pw = 0 |
| `temperature_YYYYMMDD.csv` | `t_utc_ns,sensor_id,T_degC` |
| `notes_YYYYMMDD.jsonl` | notas de `daqnote` |

Tamaño estimado de `hv_fast`: 72 columnas × ~10 caracteres × 86400 s ≈ 60 MB/día [I].

## 5. Configuración JSON (borrador)

```json
{
  "schema_version": 1,
  "campaign": {
    "name": "hodo_2026_11",
    "raw_root": "/data/daq/raw",
    "slowcontrol_dir": "/data/daq/slowcontrol",
    "operator": "",
    "bias_note": ""
  },
  "digitizer": {
    "model": "VX1740",
    "link_type": "USB_A4818",
    "link_num": 62578,
    "conet_node": 0,
    "vme_base_address": "0x00000000",
    "record_length_samples": 144,
    "post_trigger_register": 106,
    "group_enable_mask": "0xFF",
    "group_dc_offset": [32767, 32767, 32767, 32767, 32767, 32767, 32767, 32767],
    "trigger": {"source": "ext_trg_in", "io_level": "NIM"},
    "max_events_blt": 16,
    "poll_sleep_s": 0.005,
    "ttt_tick_ns": 8.0,
    "window_layout": {
      "ref_old_record_start_ns": 3552,
      "pedestal": [[0, 27], [89, 144]],
      "control": [27, 36],
      "signal": [36, 45],
      "tail": [45, 89]
    }
  },
  "run": {
    "max_events": 0,
    "max_duration_s": 0,
    "segment_max_s": 3600,
    "segment_max_bytes": 2000000000,
    "autosave_s": 300,
    "flag_delay_s": 5.0
  },
  "hv": {
    "system": "SY5527",
    "address": "192.168.1.118",
    "credentials_file": "/home/xplus/.config/daq/hv_credentials.json",
    "slot": 3,
    "channels": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    "expected_on": [0, 1, 2, 3, 4, 5, 6, 7],
    "fast_params": ["V0Set", "I0Set", "VMon", "IMon", "Status", "Pw"],
    "fast_period_s": 1.0,
    "slow_period_s": 60.0,
    "call_timeout_s": 5.0,
    "stale_after_s": 3.0,
    "geco_policy": "warn"
  },
  "auto_retry": {
    "enabled": false,
    "channels": [],
    "wait_s": [300, 1800, 3600],
    "ramp_timeout_extra_s": 60,
    "stability_window_s": 600,
    "stability_imon_max_A": 1.0e-5,
    "stability_max_s": 1800,
    "reset_after_stable_s": 7200,
    "min_spacing_between_channels_s": 60,
    "inhibit_if_geco_running": true
  },
  "flags": {
    "ovc_k_over_baseline": 2.0,
    "imon_baseline_window_s": 3600,
    "vmon_tolerance_V": 0.5,
    "trip_backdate_s": 5.0,
    "blackout_window_events": 30,
    "lightleak_window_s": 120,
    "lightleak_k_over_baseline": 2.0,
    "lightleak_baseline_window_s": 3600,
    "gap_probability": 1.0e-6,
    "burst_probability": 1.0e-6,
    "pulse_threshold_pe": 5.0
  },
  "mapping": {"file": "config/mapping/hodo_mapping_2026_05.json"},
  "calibration": {"file": "config/calibration/pe_calibration_v8.json"},
  "current_estimate": {
    "mode": "identical",
    "passive_resistance_ohm": null,
    "series_resistance_ohm": null,
    "reference_currents_file": null
  },
  "temperature": {"source": "none", "period_s": 30, "device": ""},
  "notifications": {
    "enabled": false,
    "smtp_credentials_file": "/home/xplus/.config/daq/smtp.json",
    "to": [],
    "on_gave_up": true,
    "on_trip": false,
    "curl_timeout_s": 30
  },
  "monitoring": {
    "daqrun_http": "http:127.0.0.1:8080?loopback",
    "hvmon_http": "http:127.0.0.1:8081?loopback",
    "process_requests_hz": 15,
    "monitor_queue_capacity": 64
  }
}
```

Los archivos de credenciales están fuera del repositorio con modo 600; el programa se
niega a arrancar si el modo es más abierto. Al arrancar se valida todo (rangos de la
ventana dentro de N, N múltiplo de 3, canales de `auto_retry.channels` ⊂ `expected_on`,
etc.) y la configuración resuelta se copia junto a los datos y dentro de cada segmento.

### 5.1 Mapeo versionado

Archivo JSON propio (`daq/config/mapping/<version>.json`), fuente única para el DAQ:

```json
{
  "mapping_version": "2026_05-r1",
  "channels": [
    {"dig_ch": 0, "station": "hod2018", "plane": "X", "strip": 1,
     "hv_group": 1, "hv_ch": 0, "face": null,
     "confirmed": {"strip": "measurement", "hv_group": "measurement", "face": "unknown"},
     "source": "hodo_mapping.h kCfg_2026_05; hv_diag3_pretrip.C:92"}
  ]
}
```

Cada campo lleva su procedencia (`cabling_record`, `measurement`, `assumed`, `unknown`).
Fuente única compartida con el análisis: un test de `ctest` compara el JSON con
`kCfg_2026_05`/`kCfg_2026_10` de `hodo_mapping.h` (canal → plano y tira) y falla si
difieren. No se genera `hodo_mapping.h` desde el JSON: ese archivo es del chat de análisis
y tiene historia propia (ver Diferencias, punto 9). La relación grupo HV → canal HV
(grupo g → CH0(g-1)) es [I] y debe confirmarse.

## 6. Flags

### 6.1 Palabra por canal (`ch_flags[64]`, 16 bit)

Para el canal c del digitizador, g = su grupo HV, h = su canal HV, f = su cara física, y s
el snapshot HV más reciente con `t_resp_utc_ns ≤ t_evento` (si hay uno posterior ya
recibido gracias a la línea de retardo, se usa también para el bit 5).

| bit | nombre | cálculo |
|---|---|---|
| 0 | `HV_OFF` | h en `expected_on` y (`Pw` = 0 o bit ON de `Status` = 0 o bits de trip activos) |
| 1 | `HV_OVC` | bit OVC de `Status` o `IMon - offset > k × línea base`; línea base = mediana de IMon en NORMAL en la última hora |
| 2 | `HV_RAMP` | bits RUP o RDWN de `Status` |
| 3 | `HV_VDEV` | `Pw` = 1, sin rampa y \|VMon − V0Set\| > `vmon_tolerance_V` |
| 4 | `HV_UNKNOWN` | sin snapshot, snapshot con edad > `stale_after_s`, lectura fallida, `hvmon` no responde, o grupo HV del canal desconocido |
| 5 | `HV_AMBIGUOUS` | el evento cae en un intervalo en que el estado de h cambió sin saber cuándo exactamente (ver 6.3) |
| 6 | `DATA_BLACKOUT` | todos los canales del grupo g sin pulsos > umbral en las últimas `blackout_window_events` mientras los demás grupos tienen actividad normal |
| 7 | `LIGHT_LEAK` | tasa fuera de tiempo de la cara f en `lightleak_window_s` > `k` × su línea base de la última hora (excluyendo periodos con flags) |
| 8 | `RETRY_STAB` | h está en la ventana de estabilización de un auto-retry |
| 9 | `MAP_UNCONFIRMED` | grupo HV o cara del canal sin `confirmed = measurement/cabling_record` (añadido, ver Diferencias) |
| 10-15 | reservados | 0 |

Los bits de `Status` (ON = bit 0, RUP = 1, RDWN = 2, OVC = 3, OVV = 4, UNV = 5, EXTTRIP = 6,
MAXV = 7, EXT_DIS = 8, INTTRIP = 9, ...) son [I: convención habitual de CAEN para SY5527].
Se definen en un solo archivo y se verifican contra el manual del A7030 y con la carga
resistiva en la fase 2 antes de usarlos para decisiones.

### 6.2 Palabra de evento (`evt_flags`, 32 bit)

| bit | nombre | cálculo |
|---|---|---|
| 0 | `GAP` | `dt_prev` > T con exp(−r·T) = `gap_probability`; r = tasa media de la última hora (a 0.6 Hz, T ≈ 23 s) |
| 1 | `BURST` | número de triggers en una ventana deslizante con probabilidad de Poisson < `burst_probability` |
| 2 | `HW_COUNTER_JUMP` | `evt_counter_hw` salta más de 1 (eventos perdidos) |
| 3 | `BOARD_FULL` | el registro de estado indicó memoria llena desde la lectura anterior |
| 4 | `TTT_PC_MISMATCH` | \|`t_utc_ns` − `t_read_utc_ns`\| incompatible con la latencia de lectura |
| 5 | `HV_ANY_FLAG` | OR de los bits 0-5 de cualquier canal (para cortes rápidos) |
| 6 | `MONITOR_PROVISIONAL` | los bits 6/7 de este evento se calcularon con estadística incompleta (inicio del run) |

### 6.3 Ambigüedad temporal (bit 5)

Una lectura HV mide un valor en algún instante de [`t_req`, `t_resp`]. Si la lectura k−1
mostró el estado viejo y la lectura k el nuevo, el cambio ocurrió en
(`t_req(k−1)`, `t_resp(k)`]. Todo evento en ese intervalo lleva el bit 5, y los bits 0-3 se
fijan con el estado más pesimista de los dos. Para un trip, el intervalo se extiende hacia
atrás `trip_backdate_s` (= parámetro Trip, 5 s), porque la placa limita la corriente hasta
5 s antes de apagar. A 1 Hz el intervalo típico es ~1-2 s, más la latencia de llamada (no
medida). La línea de retardo de `flag_delay_s` permite poner el bit sin reescribir; si la
lectura k llega más tarde que el retardo, el evento se escribe con bit 4 y el recálculo
offline lo corrige.

Los bits 6 y 7 dependen de ventanas de varios eventos: online se ponen desde que se
detecta la condición y el monitor registra el inicio estimado del episodio en `monitor`;
el recálculo offline los extiende hacia atrás.

## 7. Auto-retry

Por canal HV, deshabilitado por defecto (`auto_retry.enabled = false` y lista de canales
vacía). Mientras está deshabilitado `hvmon` solo registra.

### 7.1 Estados

`DISABLED`, `NORMAL`, `MANUAL_OFF`, `WAIT(k)`, `RETRY_RAMP(k)`, `RETRY_STAB(k)`,
`GAVE_UP`, `UNKNOWN` (sin conexión). k = intento (1..len(wait_s)).

### 7.2 Transiciones

| desde | condición | hacia | acción |
|---|---|---|---|
| `NORMAL` | Pw 1→0 con bits de trip, o OVC/IMon alta en los `Trip` s previos | `WAIT(k+1)` | evento `trip`; correo si `on_trip`; temporizador `wait_s[k]` |
| `NORMAL` | Pw 1→0 sin indicio de trip | `MANUAL_OFF` | evento `external`; sin acción |
| `NORMAL` | estable `reset_after_stable_s` (2 h) | `NORMAL` con k = 0 | — |
| `WAIT(k)` | vence el temporizador, conexión fresca (≥ 2 lecturas buenas), sin GECO si `inhibit_if_geco_running`, ningún otro canal en rampa en los últimos `min_spacing_between_channels_s` | `RETRY_RAMP(k)` | un solo `SetChParam(Pw=1)`; nunca V0Set, I0Set ni Trip |
| `WAIT(k)` | Pw = 1 hecho por otro (detectado) | `NORMAL` (k se conserva) | evento `external` |
| `WAIT(k)` | cambio externo de V0Set, I0Set o Trip | `MANUAL_OFF` | evento `external`; automatización suspendida en el canal |
| `RETRY_RAMP(k)` | Pw leído ≠ 1 tras 5 s | `WAIT(k+1)` o `GAVE_UP` | comando fallido |
| `RETRY_RAMP(k)` | sin RUP y \|VMon − V0Set\| < tolerancia | `RETRY_STAB(k)` | inicia ventana de estabilidad |
| `RETRY_RAMP(k)` | supera V0Set/RUp + `ramp_timeout_extra_s` (~126 s con 65.9 V a 1 V/s) | `WAIT(k+1)` o `GAVE_UP` | — |
| `RETRY_STAB(k)` | `stability_window_s` seguidos con IMon < `stability_imon_max_A` y sin trip | `NORMAL` | éxito |
| `RETRY_STAB(k)` | IMon ≥ umbral | `RETRY_STAB(k)` | la ventana vuelve a empezar (el canal sigue encendido; la protección es I0Set/Trip de la placa) |
| `RETRY_STAB(k)` | `stability_max_s` sin completar la ventana | `GAVE_UP` | correo "inestable, encendido" |
| `RETRY_*(k)` | trip | `WAIT(k+1)`, o `GAVE_UP` si k = len(wait_s) | — |
| `RETRY_*` | acción manual (Pw=0 externo, cambio de V0Set/I0Set/Trip) | `MANUAL_OFF` | la acción manual manda |
| `GAVE_UP` | — | — | un correo por canal y episodio |
| `GAVE_UP` | Pw = 1 externo | `NORMAL` (k se conserva hasta 2 h estables) | evento `external` |
| `MANUAL_OFF` | Pw = 1 externo | `NORMAL` | — |
| cualquiera | pérdida de conexión | `UNKNOWN` (se recuerda el estado anterior y sus temporizadores) | bit 4 en `daqrun` |
| `UNKNOWN` | reconexión y ≥ 2 lecturas completas | estado anterior, reconciliado | si Pw/V0Set/I0Set cambiaron durante el corte, el cambio es `external` y se aplican las reglas de acción manual; no se envía ningún comando antes de la reconciliación |

### 7.3 Otras reglas

- Origen de cada cambio: `daq_auto` si coincide con un comando enviado por `hvmon` en los
  10 s previos; `trip` según 7.2; `external` en cualquier otro caso (GECO u otro).
- Arranque de `hvmon`: un canal encontrado apagado queda en `MANUAL_OFF`, salvo que el
  estado persistido diga `WAIT` con un trip registrado hace menos de `max(wait_s)`; nunca se
  actúa en el primer minuto tras arrancar.
- Correo: proceso `curl --max-time` aparte; el fallo de un correo se registra y no bloquea
  nada. Credenciales en archivo con modo 600 fuera del repositorio.
- Rehabilitar tras un trip podría requerir más que `Pw = 1` en el SY5527 (por ejemplo
  borrar la alarma); no se sabe. Se verifica en la fase 2 con la carga resistiva.

## 8. Corriente estimada por MPPC

Siempre marcada "ESTIMADA" en ramas, CSV y web.

- `identical`: I_est = (IMon − offset − V/R_pasiva) / N_g, con N_g = número de MPPC del
  grupo según el mapeo (8; el grupo 6 tiene 7 conocidos).
- `reference_weighted`: I_c = (IMon − offset − V/R_pasiva) × I_ref,c / Σ I_ref del grupo,
  con I_ref de "Voltajes Hodoscopio.xlsx" (198-278 nA). Esa planilla no está en este host
  (búsqueda en / sin resultado).
- Corrección opcional V_MPPC,est = VMon − I × R_s cuando se conozcan R_s y R_pasiva
  (desconocidas). Mientras sean `null`, el término pasivo es 0 y la cabecera lo dice.
- Offsets: `hvmon` registra IMon/VMon de todo canal con Pw = 0 y VMon < 2 V durante más de
  60 s; `hvctl offsets` hace la medida explícita de inicio/fin de campaña y se niega si
  algún canal tiene Pw = 1 (no apaga nada).

## 9. Coexistencia y seguridad

- `hvmon` toma `flock` exclusivo de `/run/user/1000/hv_wrapper.lock`; `hvctl` en modo
  directo también, y si no lo obtiene lee a través de `hvmon`. `daqrun` toma
  `/run/user/1000/digitizer.lock`.
- Al arrancar `hvmon`: buscar GECO2020 en ejecución (patrón configurable; `pgrep -f caengeco`
  no se pudo verificar porque GECO no estaba corriendo). `geco_policy` = `warn` o `refuse`.
  Con GECO detectado se inhibe el auto-retry si `inhibit_if_geco_running`. Se re-chequea
  cada minuto y se registra.
- HTTP de solo lectura (`SetReadOnly(kTRUE)`), enlazado a 127.0.0.1 por defecto; acceso
  remoto por túnel SSH o Tailscale. Ningún comando HV por HTTP.
- `daqrun` crea un directorio de campaña nuevo con un marcador `CAMPAIGN_OPEN`; nunca
  escribe en un directorio sin marcador (todas las campañas existentes) ni en uno sin
  permiso de escritura. `daqcampaign close` borra el marcador moviéndolo a un registro dentro
  del manifiesto y aplica `chmod a-w` al directorio. `daqcampaign` nunca toca campañas
  existentes.
- Las llamadas al Wrapper pueden bloquearse (TCP). Corren en un hilo propio con watchdog;
  si una llamada supera `call_timeout_s`, el estado pasa a `UNKNOWN`. Si no vuelve en
  5 min, `hvmon` sale con error para que systemd lo reinicie.

## 10. Temperatura

`ITemperatureSource::read() → vector<{sensor_id, T_degC, t_utc_ns}>`. Implementaciones:
`Ds18b20Sysfs` (adaptador USB 1-Wire con el driver w1 del kernel: `/sys/bus/w1/devices/28-*`;
hoy `/sys/bus/w1` no existe), `SerialLineSource` (microcontrolador por USB serie con una
línea de texto por lectura; el usuario `xplus` no está en el grupo `dialout`) y
`FakeTemperature`. La lee `hvmon` (proceso continuo), no `daqrun`; `daqrun` copia las
lecturas del periodo del run en su árbol `temperature`.

## 11. Pruebas

### 11.1 Fakes

- `FakeHVBackend`: canales con V0Set/I0Set/RUp/Trip, rampa a RUp, IMon con offset y ruido,
  y un guion temporal: trip en t, sobrecorriente creciente, cambio manual (simula GECO),
  pérdida de conexión de duración d, llamadas lentas, `Pw = 1` que no prende,
  re-trip durante la estabilización.
- `FakeDigitizer`: (a) reproducción de eventos de v8/v10 abiertos en READ desde
  `/data/daq/raw` recortados a la ventana reducida, con tiempos regenerados como Poisson o
  copiados de `timestamp`; (b) sintético con pedestal, ruido, pulsos de 1 p.e. y muones,
  canales muertos, apagones de grupo y huecos de trigger programables.
- `FakeTemperature` y `FakeClock` (tiempo acelerado: la máquina de auto-retry recorre 2 h
  en milisegundos).

### 11.2 Framework

No hay framework instalado. Recomendado: GoogleTest del paquete `libgtest-dev` de Ubuntu
24.04 (requiere `sudo apt install`, decisión del operador). Alternativa sin sudo: doctest
(un solo header MIT) copiado en `daq/third_party/`. Ambos se ejecutan con `ctest`.

### 11.3 Criterios de aceptación por fase

| fase | aceptación |
|---|---|
| 0 HV lib + `FakeHVBackend` + `hvctl` | tests unitarios pasan; `hvctl --fake` y `hvctl` real (solo lectura, tarea autorizada) dan la misma tabla que `hv_readall`; latencia por llamada medida y registrada; ninguna llamada `Set*` en el binario de fase 0 (verificado con `nm`) |
| 1 `hvmon` registro | 24 h con fake acelerado sin pérdida de líneas; cambios del guion detectados con el origen correcto en el 100 % de los casos; reconexión tras corte; 8 h reales de solo lectura con GECO abierto, cada cambio hecho en GECO aparece como `external` |
| 2 auto-retry y notificador | todas las transiciones de 7.2 cubiertas por tests con tiempo acelerado; acciones manuales siempre mandan; un correo por episodio; luego una vez en hardware real con carga resistiva e I0Set bajo, con autorización explícita: trip provocado → WAIT → RETRY → NORMAL, y GAVE_UP con el correo |
| 3 `daqrun` con digitizador fake | 24 h aceleradas sin perder eventos (contador = entradas); apagado limpio por señal en cualquier momento; segmentos con SHA256 correctos; flags del guion correctos; memoria estable (RSS sin deriva > 10 %) |
| 4 digitizador real | autorizado por el operador; `GetRecordLength` = N pedido; pico de muón en la muestra esperada ±1; paso del TTT medido; 1 h de datos comparables con v8 (tasa, ocupación por canal) |
| 5 monitoreo y flags | todos los paneles de la sección 3 de requisitos en :8080 y :8081; flags coherentes con `hv_changes` y con los apagones de datos en una reproducción de v8 |
| 6 resistencia 24 h | 24 h reales con HV y digitizador: 0 eventos perdidos sin flag, `hvmon` continuo a través de un reinicio de `daqrun`, CPU < 25 % de un núcleo, memoria estable |

## 12. Plan por fases

0. Biblioteca HV (`IHVBackend`, `CaenHVBackend` solo lectura, `FakeHVBackend`), `IClock`,
   configuración, `hvctl` de solo lectura (port de `hv_readall`), framework de tests.
1. `hvmon`: sondeo 1 Hz / 1 min, archivos de slow control, detector de cambios y origen,
   estado en `/dev/shm`, página :8081, offsets, temperatura con fake.
2. Auto-retry y notificador: primero con fakes y tiempo acelerado; después una vez en
   hardware real con carga resistiva e I0Set bajo, con autorización explícita para `Pw = 1`.
3. `daqrun` con `FakeDigitizer`: hilos, colas, escritor con segmentos y SHA256, cliente HV,
   flags, notas, campañas.
4. Digitizador real (autorización del operador): configuración de ventana reducida, TTT,
   contador de hardware, memoria llena.
5. Paneles de monitoreo y flags de punta a punta (incluye ganancia horaria de hod2018).
6. Prueba de 24 h.

## 13. Diferencias con la propuesta

1. **C++17 en vez de C++20.** ROOT 6.32.10 está compilado con C++17. Un TU C++20 con ROOT
   compila y enlaza, pero `RConfigure.h:30` avisa de la discrepancia y no se probó la
   interacción con cling. Se usa C++17 en todo `daq/` con un `JThread`/`StopToken` propio de
   interfaz idéntica, para migrar cuando ROOT se recompile con C++20.
2. **Framework de tests.** GoogleTest y Catch2 no están instalados. Se propone
   `libgtest-dev` (requiere sudo) o doctest copiado en el repositorio.
3. **Slow control en texto, no en ROOT.** Un daemon que escribe un TFile durante días
   pierde datos si muere (como v10). CSV/JSONL con `fsync` periódico son a prueba de cortes y
   legibles con `grep`. Las mismas lecturas se copian al árbol `hv` de cada run.
4. **Runs segmentados en archivos de 1 h y AutoSave cada 5 min.** Limita la pérdida ante un
   fallo y permite cerrar con SHA256 durante una campaña larga. El nombre lleva `_seg<MMM>`.
5. **ADC crudo `UShort_t[64][N]` en vez de `vector<float>` en mV.** Sin pérdida, la mitad de
   bytes, tamaño fijo por run; la conversión va en la cabecera. Requiere un lector de
   compatibilidad para los programas de análisis existentes.
6. **Estado de `hvmon` con secuencia y anillo de cambios.** Solo el último valor perdería
   cambios si `daqrun` se salta una lectura. Además `flock` garantiza un único dueño del
   Wrapper entre nuestras herramientas.
7. **Línea de retardo de flags en el escritor** (5 s) para poner el bit 5 retroactivamente,
   y flags recalculables offline. El cálculo puramente online no puede conocer el intervalo
   de ambigüedad.
8. **Temperatura en `hvmon`**, porque es continua y relevante entre runs.
9. **Mapeo: JSON propio + test de consistencia contra `hodo_mapping.h`**, no generación de
   uno desde el otro. `hodo_mapping.h` es del chat de análisis, tiene dos configuraciones y
   funciones (`ch_to_plane`) que solo valen para una; además grupos HV y caras no están en
   él. El JSON añade procedencia por campo.
10. **Bit 9 `MAP_UNCONFIRMED` y flags de evento extra** (salto del contador de hardware,
    memoria llena, discrepancia TTT/PC), porque hoy 5 de 8 grupos HV son desconocidos y el
    DAQ actual no detecta eventos perdidos.
11. **Lectura por sondeo con espera de 5 ms y varios eventos por BLT**, no interrupciones,
    en la primera versión. Las IRQ por CONET quedan como opción si la medida lo justifica.
12. **THttpServer con `SetTimer(0)` y enlazado a 127.0.0.1** para que solo el hilo principal
    procese peticiones y la página no quede abierta en la red de la universidad.
13. **Notas por archivo append-only** (`daqnote` escribe con `O_APPEND`), sin IPC; funcionan
    con o sin `daqrun` corriendo.
14. **Salvaguardas extra del auto-retry:** inhibido con GECO abierto (configurable), un solo
    canal en rampa a la vez con separación mínima, ninguna acción al arrancar sobre canales
    encontrados apagados, estado persistido entre reinicios.
15. **Herramienta `daqcampaign`** para abrir y cerrar campañas con marcador, porque "nunca
    escribir en una campaña existente" y "añadir el SHA256 al cerrar cada run" exigen
    distinguir la campaña abierta propia de las ya cerradas.

## 14. Preguntas abiertas

1. ¿Qué mapeo vale para la campaña nueva: `kCfg_2026_05` (medido en v8) o `kCfg_2026_10`
   (hoja "Octubre 2026")? ¿Se recableó entre junio y octubre? AGENTS.md dice que hod2018 es
   secuencial, en contra de `hodo_mapping.h`.
2. Composición de los grupos HV 3, 4, 5, 7, 8 y del octavo canal del grupo 6; ¿existe un
   registro de la PCB de distribución? ¿El grupo g corresponde siempre a CH0(g-1)?
3. Caras físicas de hod2018 (para el aviso de fuga de luz): no hay registro.
4. Valores de las resistencias de la PCB (término pasivo y R_s): desconocidos.
5. "Voltajes Hodoscopio.xlsx" no está en este host.
6. Ventana reducida: ¿diseño A (sin pickup, pedestal en dos tramos) o B (pickup dentro del
   pedestal)? ¿N = 144 está bien?
7. Con la ventana reducida la región sin señal cae de ~250 a ~80 muestras por evento:
   ¿alcanza la estadística de pulsos de 1 p.e. para un ajuste horario de G_c? No medido.
   Opción: triggers de software intercalados (p. ej. 2 Hz), si se puede distinguir su origen
   en el evento.
8. TTT del V1740 con ROC 04.29: ancho, paso (8 o 16 ns) y si hay TTT extendido. Se mide en
   la fase 4.
9. Bits de `Status` del A7030 y qué hace falta para rehabilitar un canal tras un trip
   (`Pw = 1` solo, o borrar alarma). Se verifica en la fase 2.
10. Latencia de `CAENHV_GetChParam` contra este SY5527 (firmware antiguo): no medida. Define si
    caben 6 llamadas por segundo. ¿Conviene el modo suscripción (`SubscribeChannelParams`)?
    El log de GECO (`CAENGECO2020SubscrFail_SY5527.log`) muestra fallos de suscripción de dos
    parámetros de sistema.
11. Nombre exacto del proceso de GECO2020 (no estaba corriendo; no se pudo verificar
    `pgrep -f caengeco`).
12. Servicio: `loginctl` da `Linger=no`, así que un servicio systemd de usuario se detiene al
    cerrar sesión. ¿`enable-linger` o servicio de sistema? Ambos requieren al operador.
13. ¿Quién debe ver las páginas web y desde dónde (localhost, Tailscale, red de la
    universidad)? ¿Autenticación?
14. Proveedor SMTP, remitente, destinatarios.
15. Hardware de temperatura (adaptador 1-Wire USB o microcontrolador); permisos `dialout`.
16. Instalación del framework de tests con sudo.
17. ¿Uso de `/data/daq/manifests/` por el DAQ nuevo (hoy contiene
    `raw_manifest_20261018.sha256`)?
18. Hay un proceso `agy` (Gemini/Antigravity) con este repositorio añadido como directorio.
    ¿Escribe aquí? Afecta la regla de una tarea a la vez.
