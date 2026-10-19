# DAQ v2 con control lento de HV: documento de diseño

Estado: borrador revisado (DAQ2-0 y DAQ2-1, 2026-10-18); implementadas las fases 0-3 (DAQ2-2 a
DAQ2-5, hasta 2026-10-19): `hvmon`, `hvctl`, auto-retry, notificador y `daqrun` con digitizador
simulado. Las decisiones fijadas por la implementación están en 15.1-15.4.
Las afirmaciones llevan [V: fuente] si se verificaron en este host, [I: razón] si son
inferidas, [H: sección] si vienen de `docs/handoff_analysis_2026-10-18.md` (versión de
referencia: la del proyecto), y "no medido" si nadie las midió.

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
  (MBLT) sin espera, imprime el contador en cada vuelta, decodifica y llena el TTree. Sin
  hilos, sin cola, sin manejo de errores de `ReadData`. El `TriggerTimeTag` solo se imprime
  para los 3 primeros eventos. El `EventCounter` del hardware no se guarda.
- Salida: TTree `data_dig` con `Ch0..Ch63` (`vector<float>`, mV), `time` (`vector<float>`,
  ns), `event` (int), `timestamp` (`Long_t`, µs de `steady_clock`, una vez por lectura)
  [V: CaenDigitizerSiLab.cxx:181-205; tipo leído con ROOT]. Se escribe en el directorio
  actual. El TTree se escribe con `Write` al final; si el proceso muere se pierde lo
  posterior al último AutoSave (caso v10).
- v8: 143814 entradas, 3692974698 bytes en disco, factor de compresión 6.35, AutoSave
  cada 1176 entradas [V: ROOT, lectura READ].
- Placa: VX1740, ROC FPGA 04.29 build 8716, AMC FPGA 00.13 build 8A15
  [V: logs/commissioning_2026-05/*.log].

### 1.2 Trigger por software en el código actual

- La biblioteca sabe enviar triggers de software: `readEvents`/`readEventsTree` llaman
  `CAEN_DGTZ_SendSWtrigger` en cada vuelta del bucle si `automatic == true`
  (CaenDigitizerSiLab.cxx:323-324, 454-455, 581-582), e `init()` deja
  `SetSWTriggerMode(ACQ_ONLY)` (línea 149).
- `test_acquire_sw.cxx` (versión del 2026-05-15) **no** usa trigger por software: llama
  `readEventsTree(nevents, false, 0, timeout, triggerSource)` (test_acquire_sw.cxx:47), es
  decir `automatic = false`. Con `triggerSource 8` (config_vx1740_autotrigger_or.conf:25)
  configura el OR de self-trigger de grupos y deshabilita el SW trigger
  (CaenDigitizerSiLab.cxx:1020); con `triggerSource 9` (config_vx1740_swtrigger.conf) usa
  el trigger externo y también deshabilita el SW trigger (líneas 545-546). Imprime
  "Starting VX1740 self-trigger OR acquisition..." (línea 46).
- Sí hubo una versión anterior con trigger por software en este VX1740: el log del
  2026-05-11 (`logs/commissioning_2026-05/acquire_sw.log`) muestra "Starting
  software-trigger acquisition...", máscara de fuentes 0x810C = 0xC0000000 (SW y externo
  habilitados), 10 eventos, "<rate>: 294.10 Hz" con 1023 muestras. Esa fuente no está en
  `attic/` y las dos versiones de `test_acquire_sw.cxx` en git (`885795e` y `6acf5ea`, ambas
  del 2026-10-16) son la del 15 de mayo, sin SW trigger. Conclusión: el SW trigger funciona en esta placa con
  CAENDigitizer 2.17.3, pero hay que escribirlo de nuevo.
- Fuente recuperada (DAQ2-3, solo lectura, nada restaurado al árbol): en el árbol colgante
  `eccefd4ead41bd8b0fe493ade20467a28a645183` de `.git` está el blob
  `31394dbdc72a803e98071620766a3b4ac7c9b1f0` (1690 bytes, ruta `test_acquire_sw.cxx`) que imprime
  "Starting software-trigger acquisition...", lee `config_vx1740_swtrigger.conf` y llama
  `readEventsTree(nevents, true, 0, timeout, triggerSource)`: el "programa de SW trigger" era
  el mismo `main` con `automatic = true`, que hace que la biblioteca llame
  `CAEN_DGTZ_SendSWtrigger` en cada vuelta del bucle (CaenDigitizerSiLab.cxx:581-582). Se ve con
  `git cat-file -p 31394dbd`.

### 1.3 Ventana reducida: viabilidad y elección (estudio con v8)

Granularidad y post-trigger:

- Pedidos de 1024 y 625 muestras dieron 1023 y 624 [V: logs, entrada 0 de v6]. Ambos son
  el múltiplo de 3 inferior; la regla "la longitud se redondea hacia abajo a un múltiplo de
  3" es una inferencia [I] que se confirma en la fase 4 con `GetRecordLength`. N se elige
  siempre múltiplo de 3.
- Con 0x8114 = 364 el pulso de muón cae en las muestras 258-267 (4128-4272 ns) y
  624 − 364 = 260 muestras de pre-trigger = 4160 ns. El registro se comporta como
  "muestras post-trigger" con paso de 1 muestra [I].

Estudio (DAQ2-1). Macro `outputs/daq2_window_study_20261018/window_study.C`, solo lectura
sobre v8: 20 bloques de 700 entradas consecutivas repartidos en el run, 9057 eventos
usados y 4943 excluidos (intervalos de `good_intervals_2026.csv` y de HV-DIAG-5:
2026-05-28 10:00-20:40, 2026-05-29 08:30-11:50, 2026-05-30 11:00 a 05-31 00:59 y la última
hora). Índices en el registro de 624 muestras:

| región | diseño A (inicio 3552 ns, N = 144) | diseño B (inicio 3344 ns, N = 144) |
|---|---|---|
| pedestal pre | [222,249) = 27 muestras | [209,249) = 40 muestras |
| control | [249,258) = 9 | [249,258) = 9 |
| señal | [258,267) = 9 | [258,267) = 9 |
| cola | [267,311) = 44 (704 ns) | [267,311) = 44 |
| pedestal post | [311,366) = 55 | [311,353) = 42 |

Referencia por evento y canal: línea base en [20,200), fuera de ambos diseños.
Estimadores: mediana (cuantizada en pasos de 0.4884 mV = 2000/4095) y media recortada
(promedio de las muestras a ±2 mV de la mediana). Precisión = desviación estándar de
(estimador − referencia), con |x − mediana| < 10 mV, restando en cuadratura el ruido de la
referencia medido con la diferencia de sus dos mitades. Resultados (mediana sobre canales;
hod2018 sin ch3, que está muerto; entre paréntesis mínimo-máximo):

| magnitud | A | B |
|---|---|---|
| (a) σ pedestal pre, media recortada, hod2018 [mV] | 2.272 (1.481-3.217) | 1.927 (1.147-2.945) |
| (a) σ pedestal pre, media recortada, hod2019 [mV] | 0.2916 (0.0656-0.4685) | 0.2495 (0.0559-0.3881) |
| (a) eventos con \|Δ\| > 10 mV, hod2018 | 0.01137 | 0.00397 |
| (a) sesgo del pedestal pre+post en canales con muón (> 10 p.e. en tiempo), hod2018 [mV, mediana] | +3.907 | +2.93 |
| (b) registros con pulso > 3 p.e. en el pedestal pre, hod2018 | 0.1463 (0.07331-0.3214) | 0.2027 (0.1009-0.4094) |
| (b) registros con pulso > 3 mV en el pedestal pre, hod2019 | 0.01071 (0-0.08215) | 0.0159 (0-0.1164) |
| (c) tiempo sin señal por registro (pre + control + post) | 91 muestras = 1456 ns | 91 muestras = 1456 ns |
| (c) tasa de pulsos > 3 p.e. en esa zona, hod2018 [kHz] | 330.6 (150.3-749.3) | 331.4 (151.5-736.1) |
| (c) tasa de pulsos > 0.5 p.e., hod2018 [kHz] | 2340 (1547-3595) | 2310 (1529-3545) |
| (c) tasa de pulsos > 3 mV, hod2019 [kHz] | 23.81 (0-214.8) | 24.27 (0-213.7) |

Para comparar: el registro completo de v8 tiene 531 muestras sin señal (8496 ns) y da
332 kHz (> 3 p.e., hod2018) y 22.46 kHz (> 3 mV, hod2019). La ganancia usada para los
umbrales en p.e. es la de `pe_calibration_v8.txt`, válida para v8 (no para octubre).

Pickup en 3510-3525 ns: **no se observa en v8**. La posición de los pulsos > 3 p.e. de
hod2018 por muestra es plana entre 3280 y 4112 ns (|z| ≤ 2.4 respecto del nivel base de
1507.7 pulsos por muestra en [30,200)), mientras el muón aparece en las muestras 260-264
(z hasta +14.0). La forma de onda media de hod2018 en las muestras 210-229 tampoco
muestra una estructura: varía entre 1.312 y 1.408 mV bajo la referencia (es el efecto
medio de los pulsos oscuros). El "t_peak = 3525 ns" de
`outputs/timesync_test_20260531_174521/summary_timesync.txt` corresponde al inicio de la
ventana de búsqueda [3500, 7000] ns de ese análisis [I: artefacto de borde]. Queda anotado
para el chat de análisis (`docs/notes_for_analysis.md`).

Decisión: **diseño B**, N = 144, registro desde 3344 ns del registro antiguo.

- B da mejor precisión de pedestal pre que A (1.927 contra 2.272 mV en hod2018, 0.2495
  contra 0.2916 mV en hod2019) y menos eventos atípicos (0.00397 contra 0.01137), a costa
  de más registros con algún pulso dentro del pedestal (0.2027 contra 0.1463), que la media
  recortada tolera.
- El pedestal posterior no sirve para la línea base de un canal con muón: el sobreimpulso
  de hod2018 lo desplaza +2.93 mV (B) o +3.907 mV (A). El pedestal del evento se calcula
  solo con la región pre; la región post se usa para conteo de pulsos oscuros y, en canales
  sin pulso en tiempo, como verificación de la línea base.
- El tiempo sin señal es el mismo en los dos diseños (1456 ns) y las tasas medidas en él
  coinciden con las del registro completo, así que el conteo de pulsos no favorece a ninguno.
- No se excluye la zona 3510-3525 ns porque no hay pickup en v8. El DAQ registra la
  posición de pulsos por muestra en el monitor; si el pickup aparece en octubre, se excluye
  en el análisis sin cambiar la ventana.
- Registro de post-trigger esperado: trigger en la muestra 260 del registro antiguo = 51
  del nuevo; post = 144 − 51 = 93 → 0x8114 ≈ 93 [I]. Aceptación en la fase 4: pico de muón
  en las muestras 49-58 del registro nuevo ±1.
- Volumen: 64 × 144 × 2 B = 18432 B por evento; a 0.6 Hz, 0.955 GB/día sin comprimir [I].

Estadística por hora (0.6 Hz → 2160 registros/h): hod2018 ~3.4 pulsos > 0.5 p.e. por
registro y canal (~7400/h) en la zona sin señal, suficiente para un ajuste de ganancia
horario [I: HV-DIAG-3 ajustó con ~2.7 veces más tiempo por registro]; hod2019 ~0.035
pulsos > 3 mV por registro y canal (~76/h): las medidas de hod2019 necesitan agregación o
el modo de calibración.

### 1.4 Lo que se reutiliza del análisis

- `analysis/src/hodo_mapping.h`: `kCfg_2026_10` es el mapeo de la campaña nueva [H: 5].
  Defecto conocido: `ch_to_plane()` y `ch_to_strip()` fijan 19X = 48-63 y 19Y = 32-47
  (válido solo para `kCfg_2026_05`) y dan resultados erróneos con `kCfg_2026_10`
  (19X = 32-47). Además AGENTS.md dice "hod2018 X bar k = channel k", lo que contradice
  `kCfg_2026_05`. No se corrige aquí: ambos van en `docs/notes_for_analysis.md`.
- `analysis/src/hodo_consts.h`: ventanas en ns absolutos del registro de 624 muestras. La
  cabecera del run guarda `ref_old_record_start_ns` para traducirlas.
- `outputs/pe_calibration_v8.txt`: **no se usa para octubre** (otros voltajes) [H: 8.8]. Sirve
  solo para el estudio de la sección 1.3 y como formato de referencia.
- `analysis/src/hodo_paths.h`: `resolveRawFile()` funciona sin cambios con nombres únicos.

### 1.5 Toolchain y bibliotecas verificados

- GCC 13.3.0, CMake 3.28.3, ROOT 6.32.10 compilado con C++17 (`root-config --cxxstandard`
  = 17), `libRHTTP.so` presente, nlohmann-json3-dev 3.11.3, curl 8.5.0. GoogleTest no está
  instalado todavía (René lo instala antes de la fase 0).
- Un TU C++20 con `TFile.h`, `THttpServer.h` y `std::jthread` compila y enlaza, con el aviso
  de `RConfigure.h:30` ("The C++ standard in this build does not match ROOT configuration
  (201703L)").
- Bibliotecas CAEN duplicadas [V: sha256sum, readelf, ldconfig -p, ldd]:

  | biblioteca | /usr/lib | /usr/local/lib |
  |---|---|---|
  | CAENDigitizer | `libCAENDigitizer.so.v2.17.3` (+ symlink `.so`), 220272 B, 2025-05-06 | `libCAENDigitizer.so` archivo regular, 220272 B, 2024-11-11 |
  | CAENComm | `libCAENComm.so.v1.7.0` (+ symlink), 23216 B, cadena "1.7.0" | `libCAENComm.so` regular, 23216 B, "1.7.0" |
  | CAENVME | `libCAENVME.so.v4.0.2` (+ symlink), 85160 B | `libCAENVME.so` regular, 85160 B |

  Cada par es idéntico byte a byte (mismo SHA256). Ninguna tiene SONAME. `ldconfig -p`
  lista primero `/usr/lib` (`/etc/ld.so.conf.d` incluye `/usr/lib` en varios archivos que
  se leen antes que `libc.conf`, que trae `/usr/local/lib`), y `test_acquire_ext_nim`
  carga las tres desde `/usr/lib`.
  Decisión: el DAQ nuevo enlaza contra los archivos versionados de `/usr/lib` con target
  IMPORTED de CMake y ruta completa (`/usr/lib/libCAENDigitizer.so.v2.17.3`, etc.). Como no
  tienen SONAME, `DT_NEEDED` queda con la ruta absoluta y no depende del orden de búsqueda;
  se verifica con `readelf -d` en la fase 4. CAENDigitizer resuelve CAENComm y CAENVME por
  nombre, pero las copias son idénticas. Al arrancar se registra la versión de la biblioteca
  (`CAEN_DGTZ_SWRelease`). HV Wrapper: solo `/usr/lib64/libcaenhvwrapper.so.6.6`.
- `libCaenDigitizerSiLab`: ningún binario del repositorio la tiene en `DT_NEEDED`
  (test_acquire_ext_nim, test_acquire_sw, test_caencomm_open, analysis/bin, attic/bin), el
  Makefile la pasa con `-L. -lCaenDigitizerSiLab` (la copia del repositorio, no la de
  `/usr/local/lib`) pero el enlazador la descarta porque los objetos ya traen los símbolos,
  ninguna macro la carga con `gSystem->Load`, y no hay `.rootrc`/`rootlogon.C` ni
  referencia en `~/.bashrc`. La copia de `/usr/local/lib` difiere de la del repositorio y
  no la usa nada del repositorio. El DAQ nuevo no la usa.

## 2. Arquitectura

Dos ejecutables (`hvmon`, `daqrun`) y herramientas de línea de comandos (`hvctl`,
`daqnote`, `daqcampaign`), sobre una biblioteca común `libdaqcore` (CMake, directorio
`daq/`).

    daq/
      CMakeLists.txt
      core/        tiempo (Clock), cola acotada, StopSource/JThread, logger, unidades (fase 0);
                   config JSON, mapeo, flags (fases siguientes). Headers en core/include/daq/core/
      core/include/daq/trigger_interval_stats.h   estadística de Δt (header-only, también para analysis/)
      hv/          HVBackend, CaenHVBackend (Wrapper), FakeHVBackend, Status, credenciales (fase 0);
                   máquina de auto-retry (fase 2). Headers en hv/include/daq/hv/
      dig/         IDigitizer, CaenDigitizer (CAENDigitizer), FakeDigitizer (replay/sintético)
      temp/        TempBackend, Ds18b20Sysfs, SerialLineSource, FakeTemperature
      pulse/       búsqueda de pulsos, histogramas de amplitud, ajuste de ganancia
      io/          escritor ROOT, CSV/JSONL de slow control y su lector a TTree/RDataFrame,
                   SHA256, índice de run, campañas
      apps/        hvmon, daqrun, hvctl, daqnote, daqcampaign
      tools/legacy/ hv_readall.c, test_sy5527.c, diag_sy5527.sh (heredados, compilados tal cual)
      tests/
      config/      ejemplos de configuración, guiones del fake y mapeos versionados

### 2.1 Procesos e hilos

| Proceso / hilo | Posee | Consume | Produce |
|---|---|---|---|
| `hvmon` hilo HV | la sesión del HV Wrapper (único en el host, con `flock`) | Wrapper (1 Hz rápido, 1/min lento), cola de comandos | `HvSnapshot` con número de secuencia; eventos de cambio; resultado de comandos |
| `hvmon` hilo lógica | máquina de auto-retry, detector de cambios, estimador de corriente, secuenciador de la prueba de grupos | snapshots | comandos `Pw` (cola hacia el hilo HV), eventos, alertas |
| `hvmon` hilo comandos | socket Unix `/run/user/1000/hvmon.sock` (0600) | peticiones de `hvctl` | peticiones validadas hacia el hilo lógica; respuestas |
| `hvmon` hilo escritor | archivos de `/data/daq/slowcontrol/` | snapshots, eventos, temperaturas | CSV/JSONL diarios, `fsync` cada `fsync_period_s` |
| `hvmon` hilo temperatura | dispositivo de temperatura | sensor (10-60 s) | lecturas |
| `hvmon` hilo notificador | procesos `curl` hijos | alertas | correos (un proceso por correo, con timeout) |
| `hvmon` hilo principal | THttpServer :8081, histogramas HV, `/dev/shm/hvmon_state.json` | snapshots | página de solo lectura, estado publicado, señales |
| `daqrun` hilo lectura | handle CAENDigitizer (único llamador), temporizador de SW trigger en modo calibración | digitizador | bloques de eventos → cola de escritura y copia → cola del monitor |
| `daqrun` hilo escritor | `TFile`/`TTree` del run, índice del run | cola de eventos, HV, notas, resúmenes | segmentos ROOT, `SHA256SUMS`, `run<NNN>.index.json` |
| `daqrun` hilo monitor | datos planos (sin ROOT) | copia de eventos (cola con descarte) | resúmenes por minuto, histogramas de amplitud, entradas de ajuste de ganancia, flags de datos (bits 6, 7) |
| `daqrun` hilo cliente HV | último estado HV recibido | `/dev/shm/hvmon_state.json` cada 1 s | snapshots HV con edad; marca "stale" |
| `daqrun` hilo principal | THttpServer :8080 y todos sus histogramas/gráficos | resúmenes del monitor | página de solo lectura; `ProcessRequests()` a 10-20 Hz; señales |
| `hvctl` | nada mientras `hvmon` corre | `hvmon` por socket (lectura y comandos) o Wrapper directo de solo lectura si `hvmon` no corre | volcado legible / JSON, tiempos de lectura |

`hvmon` es un proceso aparte porque los trips ocurren entre runs y de noche; reiniciar o
caerse `daqrun` no debe parar el registro de HV, el auto-retry ni las alertas. Toda
escritura al SY5527 pasa por `hvmon` (auto-retry, prueba de grupos, comandos de `hvctl`).

### 2.2 Primitivas de concurrencia

- C++17 en todo `daq/` (decisión 1 de la sección 13): `std::mutex` + `std::scoped_lock`,
  `std::condition_variable`, `std::atomic`. `daq::StopSource`/`StopToken`/`JThread` propios,
  pequeños (~60 líneas), con la misma interfaz que C++20 y con tests propios (parada
  cooperativa, join en el destructor, `request_stop` idempotente, callback de parada).
- `BoundedQueue<T>`: capacidad fija; `push` bloqueante (cola de escritura: nunca se
  descartan eventos; si se llena, la lectura espera y el digitizador acumula en su memoria)
  o `try_push` con descarte del más antiguo y contador (cola del monitor); `close()` despierta
  a los consumidores y `pop` devuelve vacío al vaciarse.
- `Latest<T>`: valor más reciente protegido por mutex con número de secuencia.
- Cada objeto ROOT tiene un solo hilo dueño. `ROOT::EnableThreadSafety()` al arrancar y
  `TH1::AddDirectory(kFALSE)`. THttpServer con `SetTimer(0, kTRUE)` y `SetReadOnly(kTRUE)`:
  solo atiende peticiones cuando el hilo principal llama `ProcessRequests()`.
- Sin TBB, sin colas lock-free, sin ZeroMQ. Sin diccionarios ROOT propios (ramas de tipos
  básicos, arreglos fijos, `std::string`, `std::vector` de tipos básicos).

### 2.3 Tiempo

- `IClock` con `utc_ns()` (`CLOCK_REALTIME`), `mono_ns()` (`CLOCK_MONOTONIC`) y
  `sleep_until()`; `FakeClock` avanza a mano para tests acelerados.
- Toda marca en archivos es UTC ns (`int64`). Nunca hora local: el host está en
  America/Santiago, −03 en octubre y −04 en mayo-junio (v1-v10 y `good_intervals_2026.csv`
  están en −04). Las marcas del SY5527 se ignoran (reloj ~81 min adelantado).
- Cada lectura HV guarda `t_req_utc_ns` (antes de la llamada) y `t_resp_utc_ns` (después).
- TTT: no hay manual del V1740 en el host (búsqueda en `/usr/share`, `/usr/local/share`,
  `~/`, `/opt`: solo `/usr/share/CAEN/HVWrapper/CAENHVWrapper.pdf` y XML de DPP x740). El
  ancho y el paso se miden en la fase 4: ajuste lineal de TTT desenrollado contra UTC del PC
  sobre un run; la pendiente da el paso y los saltos de los residuos el período de vuelta.
  Indicio, no medida: en el log de SW trigger del 2026-05-11 los TTT de los tres primeros
  eventos son 128233, 537049 y 831531, con una tasa media impresa de 294.10 Hz; un paso de
  8 ns da intervalos de 3.27 y 2.36 ms, compatible con ~3.4 ms de media, y 16 ns daría el
  doble. `ttt_tick_ns` va en la configuración (inicial 8.0) y el desenrollado usa la hora del
  PC de cada lectura.
- Origen del run: par (UTC antes de `SWStartAcquisition`, UTC después); el TTT se reinicia
  al arrancar la adquisición [I: comportamiento habitual del x740, se verifica en la fase 4].

## 3. Flujo de datos y apagado

### 3.1 Flujo

    SY5527 ──Wrapper──> hvmon[HV] ──snapshot──> hvmon[lógica] ──> hvmon[escritor] ──> /data/daq/slowcontrol/
                             ^                    │  └──> notificador ──> curl SMTP
                             └──── Pw ────────────┘<── hvmon[comandos] <── hvctl (socket Unix)
                                                  └──> hvmon[principal] ──> /dev/shm/hvmon_state.json, :8081

    VX1740 ──A4818──> daqrun[lectura] ──bloques──> cola escritura ──> daqrun[escritor] ──> /data/daq/raw/<campaña>/
                                     └──copia──> cola monitor ──> daqrun[monitor] ──resúmenes──> escritor y principal (:8080)
    /dev/shm/hvmon_state.json ──> daqrun[cliente HV] ──> escritor (árbol hv, flags) y principal
    daqnote ──O_APPEND──> /data/daq/slowcontrol/notes_YYYYMMDD.jsonl ──> hvmon y daqrun

- Estado publicado por `hvmon`: JSON reemplazado atómicamente (temporal en el mismo
  directorio + `rename`) cada segundo, con `seq`, `t_resp_utc_ns`, `alive`, valores por
  canal, estado de auto-retry y de la prueba de grupos, y un anillo con los últimos 256
  eventos de cambio con id creciente. `daqrun` deduplica por id.
- Línea de retardo en el escritor de `daqrun`: `flag_delay_s` (5 s) antes de escribir un
  evento, para fijar el bit 5 sin reescribir. Los flags se pueden recalcular offline con los
  archivos de slow control.
- Lectura del digitizador: `ReadData` en bucle con espera de 5 ms cuando no hay datos,
  `SetMaxNumEventsBLT` > 1, y lectura del registro de estado de adquisición para detectar
  memoria llena. El contador de eventos del hardware se guarda. Modo IRQ
  (`CAEN_DGTZ_IRQWait`) solo si el sondeo muestra un problema medido.
- Segmentos: máximo `segment_max_s` (3600 s) o `segment_max_bytes` por archivo. Al cerrar
  cada segmento: `fsync`, `chmod a-w`, SHA256 añadido a `SHA256SUMS` (con `flock`) y entrada
  añadida al índice del run (4.2). AutoSave cada `autosave_s` (300 s).

### 3.2 Apagado de `daqrun`

1. SIGINT/SIGTERM: el manejador solo pone un `std::atomic<bool>`.
2. El hilo principal pide parada al hilo de lectura.
3. Lectura: `SWStopAcquisition`, vacía con `ReadData` hasta 0 bytes, empuja los últimos
   bloques, `CloseDigitizer`, cierra la cola de escritura y la del monitor.
4. Escritor: vacía la cola y la línea de retardo, escribe HV/notas/resúmenes pendientes y la
   cabecera final, `Write` + `Close`, `fsync`, `chmod a-w`, SHA256, índice con
   `"complete": true`.
5. Monitor y cliente HV terminan al cerrarse sus colas o al pedirles parada.
6. Principal: borra THttpServer, une los hilos, código 0 si el cierre fue limpio.
7. Una segunda señal no salta el cierre del archivo. Si el escritor no termina en 120 s se
   registra el error y se sale con código ≠ 0; el índice queda con `"complete": false`.

`hvmon` ante SIGTERM: termina la llamada Wrapper en curso, no inicia comandos nuevos (una
prueba de grupos en curso se aborta dejando todos los canales como estaban al inicio de la
prueba, ver 9.3), sincroniza los archivos, persiste el estado de auto-retry, publica
`alive: false` y sale.

## 4. Formato de archivos

### 4.1 Archivo ROOT de `daqrun`

Implementado en la fase 3 (DAQ2-5); lo que difiere de esta sección está en 15.4 (decisiones 54-56).

Nombre: `<campaña>_<YYYYMMDDThhmmssZ>_run<NNN>_seg<MMM>.root` en `/data/daq/raw/<campaña>/`.

TTree `events` (una entrada por trigger):

| rama | tipo | unidad / contenido |
|---|---|---|
| `evt_seq` | `ULong64_t` | índice secuencial dentro del run |
| `evt_counter_hw` | `UInt_t` | contador de eventos de la placa |
| `ttt_raw` | `UInt_t` | TriggerTimeTag tal como llega |
| `ttt_ns` | `ULong64_t` | TTT desenrollado, ns desde el arranque de adquisición |
| `t_utc_ns` | `Long64_t` | UTC estimado del trigger = origen + `ttt_ns` |
| `t_read_utc_ns` | `Long64_t` | UTC del PC al leer el bloque |
| `t_read_mono_ns` | `Long64_t` | `CLOCK_MONOTONIC` al leer el bloque |
| `dt_prev_ns` | `Long64_t` | `ttt_ns` menos el del evento anterior |
| `trigger_kind` | `UChar_t` | 0 externo, 1 software (modo calibración) |
| `adc` | `UShort_t[64][N]` | muestras ADC crudas (N fijo por run, en la cabecera) |
| `ch_flags` | `UShort_t[64]` | flags por canal (sección 6) |
| `evt_flags` | `UInt_t` | flags de evento (sección 6) |
| `hv_seq` | `ULong64_t` | `seq` del snapshot HV usado para los flags |

Conversión a mV en la cabecera (`adc_conversion`: fórmula, `NBit` = 12, `MaxVpp` = 2.0 V, DC
offset por grupo escrito y leído), para no depender de `adc2mV` del código antiguo.

TTree `hv` (1 Hz): `seq`, `t_req_utc_ns`, `t_resp_utc_ns`, `ok`, y por canal HV (12):
`v0set_V`, `i0set_A`, `vmon_V`, `imon_A`, `status`, `pw`, `retry_state`; por canal de
digitizador (64): `i_mppc_est_A` (título "ESTIMADO"), `i_est_mode`.

TTree `hv_changes`: `t_detect_utc_ns`, `t_prev_read_utc_ns`, `hv_ch`, `param`, `old_value`,
`new_value`, `origin` (`daq_auto`, `daq_command`, `trip`, `external`), `change_id`.

TTree `notes`: `t_utc_ns`, `user`, `host`, `text`.

TTree `monitor` (por minuto): por canal, tasa de pulsos en la zona sin señal y en tiempo
(Hz) sobre los umbrales activos (sección 7.1), canales vivos por plano, tasa de trigger,
tasa por cara física y por grupo HV, ocupación de colas, eventos descartados por el monitor,
estado de calibración.

TTree `amp_hist` (por canal y hora): histograma de amplitud de pulsos fuera de tiempo con
paso de 1 LSB (0.4884 mV), base de los ajustes de ganancia y del modo sin calibración.

TTree `gain` (por canal y hora, hod2018): `t_start_utc_ns`, `t_end_utc_ns`, `ch`, `p0_mV`,
`G_mV_per_pe`, `G_err`, `n_pulses`, `chi2`, `ndf`, `status`, `ambiguous`, `seed_G`.
Para hod2019: `dark_rate_3mV_Hz` y su error (indicador de ganancia).

TTree `temperature`: `t_utc_ns`, `sensor_id`, `T_degC`.

Cabecera (`TObjString`/`TNamed`): `config_json`, `mapping_json` y `mapping_version`,
`calibration_json` y `calibration_version` (o `"none"`), `run_info_json` (campaña, run,
segmento, modo `physics`/`calibration`, host, commit git del DAQ y si estaba sucio,
versiones de CAENDigitizer/Wrapper/ROOT, modelo, serie y firmware de la placa, par de
origen TTT/UTC, inicio/fin UTC, eventos, tiempo vivo, `record_length`, registro 0x8114
escrito y leído, disposición de la ventana, `ref_old_record_start_ns`).

### 4.2 Índice del run

`<campaña>_<YYYYMMDDThhmmssZ>_run<NNN>.index.json`, junto a los segmentos, reemplazado
atómicamente tras cerrar cada segmento:

```json
{
  "campaign": "hodo_2026_10", "run": 1, "mode": "physics",
  "start_utc_ns": 0, "end_utc_ns": 0, "complete": false,
  "config_sha256": "...", "mapping_version": "2026_10-r1", "calibration_version": "none",
  "segments": [
    {"seg": 1, "file": "..._run001_seg001.root", "sha256": "...",
     "first_utc_ns": 0, "last_utc_ns": 0, "n_events": 0,
     "first_evt_seq": 0, "last_evt_seq": 0, "closed_cleanly": true}
  ]
}
```

El análisis encadena los segmentos en el orden de `segments`. `first_evt_seq` y
`last_evt_seq` permiten comprobar que no faltan eventos entre segmentos. Al cerrar el run el
índice queda sin permiso de escritura y su SHA256 entra en `SHA256SUMS`.

### 4.3 Slow control (`/data/daq/slowcontrol/`), implementado en la fase 1

Archivos diarios bajo `/data/daq/slowcontrol/YYYY/`, rotados a medianoche UTC según la marca de
tiempo de cada dato. Cada línea se escribe con una llamada `write()` (O_APPEND) en cuanto llega;
`fsync` cada `fsync_period_s` (10 s) sobre los archivos con datos nuevos, al rotar y al
cerrar (no por línea: `/data` es un HDD compartido con los datos crudos). Cada archivo tiene un
hermano `<archivo>.schema.json` (`schema`, `format`, `columns`/`fields` con `name`, `type`,
`unit`, `description`); los CSV tienen una única línea de cabecera, así que se leen directo con
`ROOT::RDF::FromCSV`.

| archivo | esquema | contenido |
|---|---|---|
| `YYYY/hv_YYYYMMDD.csv` | `daqv2.hv_fast/1` | una fila por canal por lectura rápida (1 Hz): `t_req_utc_ns,t_resp_utc_ns,read_seq,slot,ch,V0Set_V,I0Set_A,VMon_V,IMon_A,Status,Pw,read_duration_ns` |
| `YYYY/hv_full_YYYYMMDD.csv` | `daqv2.hv_full/2` (`/1` sin placa) | una fila por canal por lectura completa (60 s): los 18 parámetros, numéricos en SI (`RUp_V_per_s`, `Trip_s`, ...), y desde `d9bcfda` los de la placa `bd_BdStatus,bd_HVMax_V,bd_HIMax_A,bd_Temp_degC` |
| `YYYY/hv_changes_YYYYMMDD.jsonl` | `daqv2.hv_changes/1` | `type: change` (id, t_utc_ns, t_prev_read_utc_ns, ch, param, unit, old, new, origin, command_id, during_disconnect, read_kind) y `type: command` (command_id, phase, value, reason, requester, detail) |
| `YYYY/notes_YYYYMMDD.jsonl` | `daqv2.notes/1` | notas del operador: t_utc_ns, user, host, text |
| `YYYY/temp_YYYYMMDD.csv` | `daqv2.temp/1` | `t_utc_ns,sensor_id,T_degC` (solo si hay backend de temperatura) |
| `offsets.jsonl` | — | una línea por medida de offsets (fsync inmediato) |

**Partes** (DAQ2-4): un CSV nunca mezcla columnas. Si al abrir el archivo del día ya existe con
otro encabezado (otra versión de `hvmon`), o si cambian las columnas durante la corrida (falla la
lectura de la placa en un ciclo y la fila sale sin `bd_*`), las filas van a
`<tipo>_YYYYMMDD.2.csv` (`.3`, ...), cada una con su `.schema.json`; se vuelve a la primera parte
cuyo encabezado coincide. Los lectores deben tomar todas las partes del día (`<tipo>_YYYYMMDD*.csv*`)
y ordenar por `t_resp_utc_ns`. `hv_full_20261019.csv`, escrito antes de la corrección, mezclaba 24
y 28 columnas; se separó en `hv_full_20261019.csv` (1344 filas) y `.2.csv` (1092), con el original
en `/data/daq/slowcontrol/repair_backup/` (`docs/notes_for_analysis.md`, ítem 6).

Los eventos de conexión, GECO y comandos rechazados quedan además en el log de `hvmon`
(stderr / journal de systemd). Notas: `daqrun` copiará a su archivo las notas del intervalo de su
run (fase 3). Lector a TTree/RDataFrame: pendiente para la fase 3, que es la primera que lo
necesita; el esquema ya lo permite. Tamaño medido en la validación: sección 15.2.

### 4.5 Estado publicado por `hvmon`

`/dev/shm/daqv2/hvmon_state.json` (directorio 700, archivo 600), reescrito cada segundo con
temporal + `rename` (un lector nunca ve un archivo parcial; probado con un lector en bucle
contra 3000 reescrituras). Es el mismo JSON que devuelve el comando `status` del socket.
Esquema `daqv2.hvmon_state/1`:

| campo | contenido |
|---|---|
| `schema`, `seq`, `t_utc_ns`, `t_utc`, `final`, `pid`, `started_utc_ns` | identificación; `seq` crece con cada publicación; `final: true` en el último estado al apagar |
| `system`, `slot`, `allow_commands` | configuración relevante |
| `connection` | `connected`, `ever_connected`, `stale` (sin conexión o última lectura buena hace más de `stale_after_s` = 3 s), `last_good_read_utc_ns`, `age_s`, `last_error`, `last_error_utc_ns`, `next_connect_in_s` |
| `read_seq`, `poller_seq` | número de lectura y de publicación del hilo de sondeo |
| `stats` | `fast_reads`, `full_reads`, `read_errors`, `connect_attempts`, `connect_failures`, `reconnections`, `fast_read_mean_ms`, `fast_read_max_ms`, `full_read_mean_ms`, `full_read_max_ms` |
| `channels[]` | por canal: `ch`, `values` (numéricos en SI con la unidad en la clave: `V0Set_V`, `IMon_A`, ...; enteros: `Status`, `Pw`, ...), `status_flags` (bits decodificados), `auto_retry` (objeto de 8.3 desde la fase 2; en la fase 1 era siempre `"NORMAL"`); además `notifier` y `compressor` en la raíz |
| `changes[]` | anillo de los últimos 50 cambios (mismo formato que el JSONL) |
| `geco_running` | GECO2020 detectado (sondeo cada 10 s) |
| `writer` | `ok`, `last_error`, `last_error_utc_ns`, `errors`, `fsyncs`, `bytes`, `lines`, `open_files`, `dropped` |
| `active_command`, `active_offsets`, `last_offsets` | comando de Pw en verificación, medida de offsets en curso, última medida |
| `temperature` | `null` sin backend; si hay: `backend` y `sensors[]` (`id`, `T_degC`, `t_utc_ns`) o `error` |

Ningún campo contiene secretos.

### 4.4 Manifiestos de campaña

Implementado en la fase 3 con marcadores `CAMPAIGN_OPEN` / `CAMPAIGN_CLOSED` y la regla de
simulación (`/data/daq/sim/manifests`); ver 15.4 (decisiones 56-57).

`/data/daq/manifests/<campaña>.json`: creación y cierre (UTC), operador, configuración
inicial y su SHA256, lista de runs con su índice, SHA256 del `SHA256SUMS` final, y el
estado del mapeo y de la calibración usados. Lo escribe `daqcampaign`; el directorio ya
existe y contiene `raw_manifest_20261018.sha256`, que no se toca.

## 5. Configuración JSON (borrador)

```json
{
  "schema_version": 1,
  "campaign": {
    "name": "hodo_2026_10",
    "raw_root": "/data/daq/raw",
    "slowcontrol_dir": "/data/daq/slowcontrol",
    "manifests_dir": "/data/daq/manifests",
    "operator": "",
    "bias_note": ""
  },
  "digitizer": {
    "model": "VX1740",
    "link_type": "USB_A4818",
    "link_num": 62578,
    "conet_node": 0,
    "vme_base_address": "0x00000000",
    "library": "/usr/lib/libCAENDigitizer.so.v2.17.3",
    "record_length_samples": 144,
    "post_trigger_register": 93,
    "group_enable_mask": "0xFF",
    "group_dc_offset": [32767, 32767, 32767, 32767, 32767, 32767, 32767, 32767],
    "trigger": {"source": "ext_trg_in", "io_level": "NIM"},
    "max_events_blt": 16,
    "poll_sleep_s": 0.005,
    "ttt_tick_ns": 8.0,
    "window_layout": {
      "ref_old_record_start_ns": 3344,
      "pedestal": [0, 40],
      "control": [40, 49],
      "signal": [49, 58],
      "tail": [58, 102],
      "post": [102, 144],
      "excluded": []
    }
  },
  "run": {
    "mode": "physics",
    "max_events": 0,
    "max_duration_s": 0,
    "segment_max_s": 3600,
    "segment_max_bytes": 2000000000,
    "autosave_s": 300,
    "flag_delay_s": 5.0
  },
  "calibration_mode": {
    "sw_trigger_rate_Hz": 100,
    "record_length_samples": 624,
    "post_trigger_register": 364,
    "external_trigger": false,
    "store_waveform_fraction": 0.01,
    "pulse_thresholds_mV": [1.5, 3.0, 5.0],
    "store_pulse_lists": false,
    "count_window_s": 90
  },
  "hv": {
    "system": "SY5527",
    "address": "192.168.1.118",
    "credentials_file": "~/.config/daqv2/hv_credentials.json",
    "slot": 3,
    "channels": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    "expected_on": [0, 1, 2, 3, 4, 5, 6, 7],
    "fast_params": ["V0Set", "I0Set", "VMon", "IMon", "Status", "Pw"],
    "fast_period_s": 1.0,
    "slow_period_s": 60.0,
    "call_timeout_s": 5.0,
    "stale_after_s": 3.0,
    "geco_process_pattern": "caengeco",
    "geco_policy": "warn",
    "status_bits": {
      "on": 0, "ramp_up": 1, "ramp_down": 2, "overcurrent": 3, "overvoltage": 4,
      "undervoltage": 5, "external_trip": 6, "max_v": 7, "external_disable": 8,
      "internal_trip": 9, "calibration_error": 10, "unplugged": 11
    }
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
    "inhibit_when_geco_running": false
  },
  "flags": {
    "ovc_k_over_baseline": 2.0,
    "imon_baseline_window_s": 3600,
    "vmon_tolerance_V": 0.5,
    "trip_backdate_s": 5.0,
    "blackout_window_s": 60,
    "blackout_aggregation": "hv_group",
    "lightleak_window_s": 120,
    "lightleak_k_over_baseline": 2.0,
    "lightleak_aggregation": "face",
    "rate_baseline_window_s": 3600,
    "significance_p": 1.0e-4,
    "min_expected_counts": 10,
    "gap_probability": 1.0e-6,
    "burst_probability": 1.0e-6
  },
  "pulse_finder": {
    "min_separation_samples": 4,
    "provisional_thresholds_mV": {"hod2018": 20.0, "hod2019": 3.0},
    "thresholds_pe": {"hod2018": [0.5, 3.0, 5.0, 10.0]},
    "thresholds_mV": {"hod2019": [3.0, 5.0]}
  },
  "gain_fit": {
    "period_s": 3600,
    "max_peaks": 5,
    "spacing_tolerance": 0.2,
    "ambiguity_chi2_ratio": 1.2,
    "max_relative_change": 0.3
  },
  "mapping": {"file": "config/mapping/hodo_mapping_2026_10.json"},
  "calibration": {"file": null},
  "current_estimate": {
    "mode": "identical",
    "common_resistance_ohm": null,
    "series_resistance_ohm": null,
    "reference_currents_file": null
  },
  "temperature": {"backend": "none", "period_s": 30, "device": ""},
  "notifications": {
    "enabled": false,
    "config_file": "~/.config/daqv2/notify.json",
    "on_gave_up": true,
    "on_trip": false,
    "curl_timeout_s": 30
  },
  "monitoring": {
    "daqrun_http": {"bind": "127.0.0.1", "port": 8080},
    "hvmon_http": {"bind": "127.0.0.1", "port": 8081},
    "process_requests_hz": 15,
    "monitor_queue_capacity": 64
  },
  "slowcontrol": {"fsync_period_s": 10}
}
```

- Borrador: los bloques `auto_retry` y `notifications` de arriba se reemplazaron en la fase 2
  por `auto_retry` y `notifier` de `hvmon.json` (8.1, 8.4); `notify.json` no existe: el SMTP
  va en `secrets.env` (`DAQ_SMTP_*`, `DAQ_NOTIFY_*`), modo 600.
- `monitoring.*.bind`: `"127.0.0.1"` (por defecto) o `"tailscale"` (se resuelve a la
  dirección de `tailscale0` al arrancar; hoy 100.114.138.29 [V: ip -br addr]). Cualquier otro
  valor, `0.0.0.0` o la dirección de `enp0s31f6` (red de la universidad) se rechazan al
  validar.
- Validación al arrancar: rangos de la ventana dentro de N y sin solaparse, N múltiplo de 3,
  `auto_retry.channels` ⊂ `expected_on`, bits de `status_bits` distintos, etc. La
  configuración resuelta se copia junto a los datos y dentro de cada segmento.

### 5.1 Mapeo versionado (contenido inicial: `hodo_mapping_2026_10.json`)

Fuente única para el DAQ. Cada canal de digitizador tiene `station`, `plane`, `strip`,
`hv_group`, `hv_ch` y `face`, y para cada campo `confirmed` (`measurement`,
`cabling_record`, `derived`, `assumed`, `unknown`) y `source`. Un test de `ctest` compara
canal → plano y tira con `kCfg_2026_10` de `hodo_mapping.h`.

Grupos HV: grupo g = canal SY5527 CH0(g−1), confirmado por René [H: 3]. Tabla de la planilla
en etiquetas (se guarda en el JSON como `hv_groups_by_label`, fuente para la prueba de grupos):

| Grupo (canal SY5527) | hod2018 X (Canal) | hod2018 Y (Canal) | hod2019 (cara canal) |
|---|---|---|---|
| 1 (CH00) | 13 | 11, 15 | Delantera 7; Frontal 4, 6; Trasera 3, 6 |
| 2 (CH01) | 7 | 7 | Delantera 3; Frontal 7; Trasera 1, 4, 5, 7 |
| 3 (CH02) | 3, 12, 14 | 5, 10 | Delantera 5; Frontal 8; Trasera 8 |
| 4 (CH03) | 1 | 2, 4, 8, 12, 13 | Frontal 2; Opuesta 7 |
| 5 (CH04) | 4, 5, 6, 8, 9 | 6 | Delantera 4; Opuesta 6 |
| 6 (CH05) | 0, 2 | 1 | Delantera 6, 8; Opuesta 3, 8; Trasera 2 |
| 7 (CH06) | 10, 11 | — | Delantera 1, 2; Frontal 1, 5; Opuesta 1, 5 |
| 8 (CH07) | 15 | 0, 3, 9, 14 | Frontal 3; Opuesta 2, 4 |

hod2019 con el cableado de octubre (`kCfg_2026_10`: X = ch 32-47, Y = ch 48-63, coordenada
creciente). Regla [H: 5]: X tira k par = Delantera canal 8 − k/2, impar = Trasera canal
8 − (k−1)/2; Y tira k par = Opuesta canal 8 − k/2, impar = Frontal canal 8 − (k−1)/2. Caras:
Delantera = ch 32, 34, …, 46; Trasera = 33, 35, …, 47; Opuesta = 48, 50, …, 62;
Frontal = 49, 51, …, 63. Grupos derivados (`confirmed = derived`):

| Grupo | canales hod2019 (octubre) |
|---|---|
| 1 (CH00) | 34, 37, 43, 53, 57 |
| 2 (CH01) | 35, 39, 41, 42, 47, 51 |
| 3 (CH02) | 33, 38, 49 |
| 4 (CH03) | 50, 61 |
| 5 (CH04) | 40, 52 |
| 6 (CH05) | 32, 36, 45, 48, 58 |
| 7 (CH06) | 44, 46, 54, 55, 62, 63 |
| 8 (CH07) | 56, 59, 60 |

Excepción: ch48 (Opuesta 8, grupo 6) queda `confirmed = measurement` (HV-DIAG-5 [H: 5]).
Comprobación de la regla con la tabla de etiquetas: grupo 1 = Delantera 7 → ch34,
Trasera 3 → ch43, Trasera 6 → ch37, Frontal 4 → ch57, Frontal 6 → ch53, que coincide con la
fila 1 [I: cálculo en este documento].

hod2018 con el cableado de octubre: `strip` según la hoja (X tira k = ch k, Y tira k = ch
16 + k), `confirmed = cabling_record` pero marcado sin confirmar, porque en mayo esa misma
hoja no se cumplió [H: 5]. `hv_group`: solo las relaciones medidas en mayo, válidas si los
cables de hod2018 no se tocaron, todas con `confirmed = assumed`:

| canal | grupo | etiqueta medida en mayo |
|---|---|---|
| 0 | 1 | X Canal 13 |
| 6 | 2 | X Canal 7 |
| 20 | 2 | Y Canal 7 |
| 22, 31 | 1 | Y Canales {11, 15} |
| 1, 3 | 6 | X Canales {0, 2} |
| 16, 17 | 6 | uno es Y Canal 1 |

Los otros 23 canales de hod2018 (32 − 9; el borrador decía 22) empiezan con `hv_group = null` (`unknown`): bit 4 por
mapeo desconocido (sección 6) hasta la prueba de grupos de la fase 6.

Discrepancia conocida: en mayo el grupo 6 contenía ch1, ch16 y ch17 de hod2018 (más ch3,
muerto, deducido), con dos MPPC del plano Y, mientras la planilla le asigna un solo MPPC Y
(Canal 1). Con las 5 entradas de hod2019 de la planilla, el grupo 6 tendría 9 MPPC, uno más
que los 8 por grupo. El JSON lo guarda en `known_discrepancies` y la prueba de grupos lo
resuelve.

## 6. Flags

### 6.1 Palabra por canal (`ch_flags[64]`, 16 bit)

Fase 3: implementados los bits 0-5 y 9 de `ch_flags` y los bits 0-5 y 7 de `evt_flags` (15.4,
decisiones 50-53); los demás quedan en 0 hasta las fases 5 y 6.

Para el canal c, g = su grupo HV, h = su canal HV, f = su cara física, s = snapshot HV más
reciente con `t_resp_utc_ns ≤ t_evento` (y el siguiente, si ya llegó, para el bit 5).

| bit | nombre | cálculo |
|---|---|---|
| 0 | `HV_OFF` | h en `expected_on` y (`Pw` = 0 o bit `on` de `Status` = 0 o bits de trip) |
| 1 | `HV_OVC` | bit `overcurrent` o `IMon − offset > k × línea base` (mediana de IMon en NORMAL en la última hora) |
| 2 | `HV_RAMP` | bits `ramp_up` o `ramp_down` |
| 3 | `HV_VDEV` | `Pw` = 1, sin rampa, \|VMon − V0Set\| > `vmon_tolerance_V` |
| 4 | `HV_UNKNOWN` | sin snapshot, edad > `stale_after_s`, lectura fallida, `hvmon` no responde, o `hv_group` del canal desconocido |
| 5 | `HV_AMBIGUOUS` | el evento cae en el intervalo de ambigüedad de un cambio de h (6.3) |
| 6 | `DATA_BLACKOUT` | prueba de Poisson de déficit (6.4) en la agregación configurada |
| 7 | `LIGHT_LEAK` | prueba de Poisson de exceso (6.4) en la agregación configurada |
| 8 | `RETRY_STAB` | h en la ventana de estabilización de un auto-retry |
| 9 | `MAP_UNCONFIRMED` | `hv_group` o `face` del canal sin `confirmed = measurement` |
| 10 | `PE_PROVISIONAL` | el canal no tiene calibración de p.e. válida (sección 7.4) |
| 11 | `STAT_INSUFFICIENT` | la prueba de los bits 6/7 no se pudo hacer (cuentas esperadas < mínimo) |
| 12-15 | reservados | 0 |

Posiciones de `Status`: en la configuración (`hv.status_bits`). Valores iniciales tomados
del manual del HV Wrapper (Rev. 13, 2016, tipo `PARAM_TYPE_CHSTATUS`): bit 0 on, 1 ramp up,
2 ramp down, 3 overcurrent, 4 overvoltage, 5 undervoltage, 6 external trip, 7 max V,
8 external disable, 9 internal trip, 10 calibration error, 11 unplugged
[V: /usr/share/CAEN/HVWrapper/CAENHVWrapper.pdf]. Ese manual describe la convención
genérica del Wrapper; el significado para el A7030DP se confirma con el manual de la placa
(René lo busca) y con la carga resistiva en la fase 2.

### 6.2 Palabra de evento (`evt_flags`, 32 bit)

| bit | nombre | cálculo |
|---|---|---|
| 0 | `GAP` | `dt_prev` > T con exp(−r·T) = `gap_probability`; r = tasa de la última hora (a 0.6 Hz, T ≈ 23 s) |
| 1 | `BURST` | número de triggers en ventana deslizante con probabilidad de Poisson < `burst_probability` |
| 2 | `HW_COUNTER_JUMP` | `evt_counter_hw` salta más de 1 |
| 3 | `BOARD_FULL` | memoria llena desde la lectura anterior |
| 4 | `TTT_PC_MISMATCH` | \|`t_utc_ns` − `t_read_utc_ns`\| incompatible con la latencia de lectura |
| 5 | `HV_ANY_FLAG` | OR de los bits 0-5 de cualquier canal |
| 6 | `MONITOR_PROVISIONAL` | bits 6/7 calculados con estadística incompleta (inicio de run) |
| 7 | `NO_CALIBRATION` | el run está en estado sin calibración (sección 7.4) |
| 8 | `HV_GROUP_TEST` | evento tomado durante un paso de la prueba de grupos |

Los bits 0 y 1 usan `daq/core/include/daq/trigger_interval_stats.h` (sección 7.5).

### 6.3 Ambigüedad temporal (bit 5)

Una lectura HV mide en algún instante de [`t_req`, `t_resp`]. Si la lectura k−1 mostró el
estado viejo y la k el nuevo, el cambio ocurrió en (`t_req(k−1)`, `t_resp(k)`]. Los eventos
en ese intervalo llevan el bit 5 y los bits 0-3 con el estado más pesimista. Para un trip el
intervalo se extiende hacia atrás `trip_backdate_s` (= Trip, 5 s). Si la lectura k llega
después de `flag_delay_s`, el evento se escribe con bit 4 y el recálculo offline lo corrige.

### 6.4 Pruebas de Poisson para los bits 6 y 7

- Unidad de agregación configurable: `channel`, `hv_group` o `face`. Por defecto `hv_group`
  para apagones (bit 6) y `face` para fuga de luz (bit 7), porque un trip afecta a un grupo y
  la luz a una cara [H: 4].
- Tasa base r por canal: tasa de pulsos sobre el umbral activo en la zona sin señal, en
  NORMAL y sin flags, durante `rate_baseline_window_s`. Esperado en la ventana W:
  μ = Σ_c r_c · t_sin_señal · n_eventos(W).
- Bit 6: se marca si μ ≥ `min_expected_counts` y P(N ≤ n_obs | μ) < `significance_p`. Bit 7:
  si μ ≥ mínimo, n_obs/μ > `lightleak_k_over_baseline` y P(N ≥ n_obs | μ) < `significance_p`.
  Si μ < mínimo: bit 11 y no se decide.
- Escala con la ventana de 1456 ns sin señal y 0.6 Hz: una cara de hod2019 (8 canales,
  ~24 kHz cada uno con los datos de v8) da μ ≈ 0.28 cuentas por evento, ≈ 20 en 2 min; un
  canal solo de hod2019, μ ≈ 2.5 en 2 min (insuficiente). Un grupo HV tiene entre 2 y 6
  canales de hod2018 según la planilla a ~331 kHz (> 3 p.e.): μ ≈ 0.48 por evento y canal [I: cálculo con
  las tasas de 1.3]. Por eso la agregación por canal no sirve para hod2019.
- Estos números dependen de la ganancia de octubre; el monitor publica μ por unidad para
  ajustar umbrales y ventanas con datos reales.

## 7. Monitoreo

Fase 3: el hilo monitor calcula solo tasa de trigger, Δt, tasa de pulsos fuera de tiempo y espectro
de amplitud por canal (mV, umbrales provisionales), niveles de cola y resumen de flags (15.4,
decisión 58). El resto de los paneles y los bits 6-8 son de la fase 5.

### 7.1 Búsqueda de pulsos

Todos los pulsos de la zona considerada, no uno por registro [H: 8.9]: mínimo local
(muestra menor o igual que sus dos vecinas) con amplitud `línea_base − v` sobre el umbral,
y separación mínima `min_separation_samples` (4 = 64 ns); dentro de esa separación se
queda el más grande. Línea base del evento: media recortada (±2 mV de la mediana) de la
región `pedestal`. Zonas: sin señal = `pedestal` + `control` + `post`, más `tail` en canales
sin pulso en tiempo; en tiempo = `signal`. Si se usa ocupación por registro en algún panel,
se convierte con tasa = −ln(1 − occ)/T. Se publica también la posición de los pulsos por
muestra (para detectar pickup).

### 7.2 Paneles (:8080 `daqrun`, :8081 `hvmon`)

Por canal: tasa fuera de tiempo y en tiempo sobre los umbrales activos; canales vivos por
plano; tasa de trigger; distribución de Δt con ajuste exponencial, huecos y ráfagas; IMon y
VMon por grupo HV; tasa y μ por cara y por grupo con su p-valor; ocupación de colas;
panel de flags; estado de calibración; histograma de amplitud por canal; G_c(t) de hod2018
y tasa oscura > 3 mV de hod2019.

### 7.3 Ajuste de ganancia por hora (hod2018)

- Entrada: histograma de amplitud de pulsos fuera de tiempo de la hora, paso de 1 LSB.
- Semilla: G y p0 de la hora anterior del mismo canal (o de la calibración de puesta en
  marcha para la primera hora).
- Ajuste multigaussiano de hasta `max_peaks` picos con posiciones p0 + n·G. Además se
  ajusta con semillas G/2 y 2G; si la mejor alternativa tiene χ²/ndf a menos de
  `ambiguity_chi2_ratio` de la elegida, el ajuste se marca `ambiguous` y no se acepta.
- Prueba de espaciado: las diferencias entre picos consecutivos ajustados libres deben
  coincidir entre sí y con G dentro de `spacing_tolerance` (20 %); si no, se rechaza (caso
  ch9 en v8: 5.58 o 2.79 mV/p.e. según la semilla [H: 8.11]).
- Si G cambia más de `max_relative_change` respecto de la hora anterior sin que haya un
  cambio de HV registrado, se marca para revisión pero se guarda.
- hod2019 no tiene picos resueltos: su indicador es la tasa oscura > 3 mV por canal y hora.
- Selecciones por amplitud (tira o cluster dominante, umbrales) en p.e.:
  (ph − p0)/G_c [H: 8.6].

### 7.4 Estado "sin calibración"

Mientras no exista la calibración de puesta en marcha de octubre (`calibration.file =
null`), el sistema funciona así y lo dice en la cabecera (`calibration_version = "none"`),
en `evt_flags` (bit 7) y en el panel:

- Todo el monitoreo funciona en mV con los umbrales `provisional_thresholds_mV` (hod2018
  20 mV ≈ 5 p.e. con las ganancias de mayo, hod2019 3 mV), con la etiqueta "PROVISIONAL".
- Umbrales en p.e., selecciones en p.e. y el panel de p.e. están deshabilitados. Los bits 6
  y 7 funcionan con los umbrales en mV (sus tasas base se aprenden igual) y llevan el bit 10.
- El ajuste de ganancia horario corre igual y publica G_c(t), pero sin semilla de
  calibración: semilla = máximo del espectro de autocorrelación del histograma; sus
  resultados no se usan para umbrales hasta que la calibración exista.
- Al cargar una calibración, el cambio queda en la cabecera del siguiente segmento; no se
  recalculan segmentos ya cerrados.

### 7.5 Estadística de Δt compartida con el análisis

`daq/core/include/daq/trigger_interval_stats.h`: C++17 header-only, sin ROOT, para que
`analysis/` lo use con `-I../daq/core/include` sin depender del build de CMake. Entrada:
tiempos de trigger en ns. Funciones: tasa por máxima verosimilitud (1/media) con error,
prueba de exponencialidad (Kolmogorov-Smirnov y χ² binado), huecos (Δt con probabilidad
< p a la tasa local), ráfagas (estadística de barrido en ventana deslizante), y resumen por
intervalo. Tests en `ctest` con secuencias sintéticas. Para v7-v10 el análisis solo tiene
`timestamp` (µs, una vez por lectura, sin TTT): la resolución es la del bucle de lectura y
eventos leídos juntos comparten timestamp; el header acepta una resolución declarada.

## 8. Auto-retry (implementado en la fase 2)

Por canal HV, deshabilitado por defecto (`auto_retry.enabled_channels` vacío). Mientras está
deshabilitado `hvmon` solo registra. La lógica es pura (`daq/hvmon/src/retry.cpp`,
`RetryMachine`): recibe la lectura rápida, los cambios detectados y la hora UTC; devuelve
transiciones, eventos para el notificador y acciones. Las acciones (`Pw = 1`, y opcionalmente
`ClrAlarm` antes) las ejecuta el hilo de sondeo como solicitudes internas por la misma ruta
que los comandos de `hvctl` (lectura previa, compuerta `allow_commands`, verificación).

### 8.1 Estados

`DISABLED` (canal fuera de `enabled_channels`), `NORMAL`, `MANUAL`, `WAIT(k)`, `RETRY(k)`
(Pw = 1 escrito, rampa en curso), `STABILIZING(k)`, `GAVE_UP`. k es el número de intento
(1..len(wait_min)). `UNKNOWN` y `GROUP_TEST` del borrador no existen: sin lectura fresca la
máquina no avanza y no escribe (bloqueo "sin lectura fresca"); `GROUP_TEST` llega en la fase 6.

Configuración (`auto_retry` en `hvmon.json`): `enabled_channels` [], `wait_min` [5, 30, 60],
`stable_min` 10, `reset_after_h` 2, `imon_stable_max_uA` 6.0, `max_ramping_channels` 1,
`inhibit_when_geco_running` false, `clear_alarm_before_retry` false, `state_file`
(por omisión `$XDG_STATE_HOME/daqv2/hvmon_retry_state.json`, o
`~/.local/state/daqv2/hvmon_retry_state.json`), `test_mode` {`enabled`, `channels`,
`wait_min` [1, 2, 3], `stable_min` 1}, `test_mode_allow_detector_channels` false. Validación:
`enabled_channels` ⊂ `channels` y exige `allow_commands = true`; los canales de `test_mode`
deben estar en `enabled_channels`, y CH00-07 (detector) se rechazan en `test_mode` salvo
`test_mode_allow_detector_channels = true`.

"Estable" = Pw = 1, bit `on`, sin bits de trip ni de rampa, e IMon − offset < `imon_stable_max_uA`.

### 8.2 Transiciones (tal como están en el código)

| desde | condición | hacia | efecto |
|---|---|---|---|
| (arranque) | canal fuera de `enabled_channels` | `DISABLED` | sin registro de cambio |
| (arranque) | persistido `WAIT(k)`, Pw = 0, y ahora < trip + wait_min[k−1] + 10 min | `WAIT(k)` | reanuda episodio, k, hora del trip y línea de tiempo |
| (arranque) | persistido `WAIT(k)` vencido hace más de 10 min | `MANUAL` | sin acción automática |
| (arranque) | persistido `GAVE_UP` y Pw = 0 | `GAVE_UP` | — |
| (arranque) | Pw = 1 | `NORMAL` | k se conserva si lo persistido era NORMAL/STABILIZING |
| (arranque) | Pw = 0 | `MANUAL` | "apagado al arrancar" |
| cualquiera salvo `MANUAL`, `DISABLED` | cambio de `Pw` o `V0Set` con origen `external` | `MANUAL` | la acción manual manda |
| `NORMAL` | trip (bit de trip, o Pw 1→0 con origen `trip`) | `WAIT(k+1)` | episodio nuevo si k = 0; evento `trip` |
| `NORMAL` | trip y k = len(wait_min) | `GAVE_UP` | eventos `trip` y `gave_up` |
| `NORMAL` | k > 0 y `reset_after_h` en NORMAL | `NORMAL`, k = 0 | registro del reinicio |
| `WAIT(k)` | Pw leído = 1 | `MANUAL` | — |
| `WAIT(k)` | `Status` (sin bits de rampa) o `V0Set` distinto del de la entrada a WAIT | `MANUAL` | — |
| `WAIT(k)` | vence wait_min[k−1] desde el trip y no hay bloqueo | `RETRY(k)` | acción `Pw = 1` (con `ClrAlarm` antes si está configurado) |
| `WAIT(k)` | vence pero hay bloqueo: sin lectura fresca; GECO corre e `inhibit_when_geco_running`; otros canales en rampa o en RETRY ≥ `max_ramping_channels` | `WAIT(k)` | `block_reason` publicado; se reintenta en cada lectura |
| `RETRY(k)` | resultado `verified` o `noop` | `STABILIZING(k)` | — |
| `RETRY(k)` | resultado `busy` (otro comando en curso) o `write_failed` | `WAIT(k)` | pospuesto, se reintenta |
| `RETRY(k)` | resultado `aborted` (lectura previa distinta) o `rejected` | `MANUAL` | — |
| `RETRY(k)` | resultado `verify_failed` (rampa no terminó) | `WAIT(k+1)` o `GAVE_UP` | intento fallido, sin evento `trip` |
| `RETRY(k)`, `STABILIZING(k)` | trip nuevo: Pw 1→0 con origen `trip`, o bit de trip con Pw = 1 | `WAIT(k+1)`, o `GAVE_UP` si k = len(wait_min) | evento `trip` (y `gave_up`) |
| `STABILIZING(k)` | estable durante `stable_min` | `NORMAL` | evento `recovery` |
| `STABILIZING(k)` | deja de estar estable | `STABILIZING(k)` | la ventana vuelve a empezar |
| `MANUAL` | estable durante `stable_min` | `NORMAL`, k = 0 | "estable tras acción manual" |
| `GAVE_UP` | — | — | solo sale por acción manual (Pw o V0Set externo → `MANUAL` → `NORMAL`) |

Un bit de trip retenido con el canal apagado no cuenta como trip nuevo en RETRY/STABILIZING
(la placa lo mantiene hasta que la relectura muestra Pw = 1). Un trip en `MANUAL` no inicia un
episodio: el canal debe haber pasado antes a `NORMAL` (prueba de la parte C, 15.3).

### 8.3 Publicación, registro y persistencia

- Cada transición es un `ChangeRecord` con `param = "auto_retry"`, `origin = "hvmon_retry"`,
  `old`/`new` como etiquetas (`"WAIT(1)"`) y `detail` {`k`, `reason`}; va a
  `hv_changes_YYYYMMDD.jsonl` y al log. Las escrituras son `CommandRecord` con autor
  `hvmon (auto-retry)` y texto de autorización `auto-retry intento k`.
- Estado (`channels[i].auto_retry`): `state`, `label`, `k`, `episode`, `enabled`, `test_mode`,
  `since_utc_ns`, `trip_utc(_ns)` (trip que abrió la espera actual), `next_retry_in_s`,
  `stable_for_s`, `block_reason`, `last_outcome`. La página tiene la tabla "Auto-retry" y
  `hvctl status` una columna y líneas por canal habilitado.
- Persistencia: `state_file` (esquema `daqv2.hvmon_retry_state/1`) por escritura atómica en
  cada transición; guarda estado, k, episodio, hora del trip y línea de tiempo por canal.
- GECO abierto no inhibe por defecto; al reintentar con GECO abierto se registra un aviso por
  intento.
- Si rehabilitar tras un trip requiere más que `Pw = 1` en el A7030DP real no se sabe
  (`ClrAlarm` existe como opción): se verifica en la prueba con carga resistiva (DAQ2-5).

### 8.4 Notificador

Hilo propio en `hvmon` (`daq/hvmon/src/notifier.cpp`). Eventos: `trip`, `recovery`,
`gave_up` (de la máquina), `writer_error` (escritor con errores), `connection_lost`
(desconectado más de `connection_lost_min`). Configuración `notifier`: `enabled` false,
`notify_each_trip` false, `notify_recovery` false, `notify_writer_errors` true,
`notify_connection_lost` true, `connection_lost_min` 5, `problem_interval_min` 30,
`retry_delays_min` [1, 5, 15], `curl` "curl", `max_time_s` 30, `mapping_file`.

- Envío: `curl --config -` lanzado con `posix_spawn`; URL, usuario, contraseña, `mail-from`
  y `mail-rcpt` van por stdin (nunca en argv); el mensaje va en un temporal 0600 en el
  directorio de runtime que se borra al terminar. `ssl-reqd`, `max-time`. Timeout del hijo con
  SIGTERM y luego SIGKILL.
- Mensaje: cabeceras `Date` (RFC 5322, UTC), `Message-ID`, `From`, `To` (varios
  destinatarios separados por comas en `DAQ_NOTIFY_TO`), `Subject` codificado RFC 2047 que
  `[DAQ hodoscopio] [PRUEBA ]<tipo> CHnn (grupo g)`; cuerpo en español, CRLF: evento, hora UTC
  del evento y del envío, host, canal y grupo HV, intento k, canales del digitizador del
  grupo según el mapeo (los no confirmados por medición marcados con `*`), motivo, línea de
  tiempo del episodio, último VMon/IMon/Status, temperatura de la placa y si GECO corre.
- Límites: un correo por canal, tipo y episodio para `trip`/`recovery`/`gave_up`; uno cada
  `problem_interval_min` por tipo para `writer_error` y `connection_lost`. Un envío fallido
  se reintenta a 1, 5 y 15 min y luego se descarta con registro. Los descartes se cuentan
  (`suppressed`).
- `hvctl notify-test [--ch N]`: correo marcado PRUEBA, sin límites ni reintentos; `hvctl`
  espera la respuesta (hasta 180 s) e imprime el código de salida de `curl`.
- Secretos: `DAQ_SMTP_URL`, `DAQ_SMTP_USER`, `DAQ_SMTP_PASSWORD`, `DAQ_NOTIFY_FROM`,
  `DAQ_NOTIFY_TO` (entorno o `secrets.env`). Con `enabled = true` y sin ellos `hvmon` no
  arranca; con `enabled = false`, aviso. Reemplaza al `notify.json` del borrador.

### 8.5 Compresión de días cerrados

- `compress_closed_days` (true por omisión). Al rotar, el escritor pasa la ruta del día
  cerrado a un hilo compresor (nice 19); al arrancar se buscan días cerrados sin comprimir
  (`hv_*`, `hv_full_*`, `hv_changes_*`, `notes_*`, incluidas las partes `.N`, con fecha < hoy;
  los `.schema.json` y `offsets.jsonl` quedan planos). Al rotar se encolan todas las partes del día.
- `gzipFileVerified` (zlib, nivel 6): escribe `.gz`, lo relee y compara bytes, líneas y CRC32
  con el original; solo entonces borra el plano. Si falla, el plano queda y el `.gz` pasa a
  `.gz.failed`.
- Los lectores (`LineReader`, `sctool cat`) aceptan plano y `.gz`. Herramienta
  `daq/build/apps/sctool/sctool compress [--keep] | cat | digest`.

## 9. Coexistencia, seguridad y comandos HV

### 9.1 Locks y GECO (implementado en la fase 1)

- `hvmon` toma `flock` exclusivo de `$XDG_RUNTIME_DIR/daqv2/hvmon.lock` (directorio 700, archivo
  600): un segundo `hvmon` sale con código 3. `hvctl` en modo directo toma el mismo lock mientras
  usa la Wrapper; si lo tiene `hvmon`, `hvctl read/status/changes/note/measure-offsets/pw` van
  por el socket y `map/params/bench` se rechazan. `daqrun` tomará
  `$XDG_RUNTIME_DIR/daqv2/digitizer.lock`.
- GECO2020: cada 10 s se busca el patrón `bin/CAENGECO2020` (ejecutable ELF del snap) en
  `/proc/<pid>/exe` o en `argv[0]` de cada proceso; el estado publica `geco_running` y el log
  cada transición. GECO abierto no inhibe nada en la fase 1. (Hasta `6eb49b2` se buscaba en toda
  la línea de comandos, lo que dio un falso positivo: sección 15.2.)
- Coexistencia medida (2026-10-18, sección 15.2): con `hvmon` leyendo cada 1 s, GECO conectado
  al SY5527 unos 17 min no produjo errores ni reconexiones en `hvmon`, ni cambios en el registro;
  GECO mostró CONNECTED y los mismos valores (René). Al conectar, GECO escribe
  `~/CAENGECO2020SubscrFail_SY5527.log` (240 B, 23:32:13Z) con dos parámetros de sistema que no
  pudo suscribir (`CMDExecMode`, `DHCPDStatus`); el contenido es idéntico al de la copia de la raíz
  del repositorio (11:11 −03, antes de que existiera `hvmon`) y `~/CAENGECO2020SubscrFail_System.log`
  es de 2026-05-11: el aviso no tiene relación con `hvmon`.

### 9.2 Ruta de comandos `hvctl` → `hvmon` (implementada en la fase 1)

- Socket Unix `$XDG_RUNTIME_DIR/daqv2/hvmon.sock`, modo 600. Una petición JSON por línea (máximo
  64 KiB, 5 s para llegar completa), una respuesta JSON por línea, una petición por conexión. Un
  cliente que corta a mitad de la petición o que no lee la respuesta no afecta al servidor.
  Las peticiones largas (`pw`, `measure-offsets`) se atienden en un hilo propio.

| petición | respuesta |
|---|---|
| `{"cmd":"ping"}` | `{"ok":true,"pid":N}` |
| `{"cmd":"status"}` | `{"ok":true,"state":{...}}` (esquema de 4.5) |
| `{"cmd":"read"}` | `{"ok":true,"reading":{...}}` (formato de `readingToJson`: `params` es un arreglo ordenado) |
| `{"cmd":"changes","n":20}` | `{"ok":true,"changes":[...]}` (últimos n del anillo de 50) |
| `{"cmd":"note","text":"..."}` | `{"ok":true,"t_utc_ns":...}`; va a `notes_YYYYMMDD.jsonl` |
| `{"cmd":"measure-offsets","seconds":N,"channels":[...]}` | `{"ok":true,"offsets":{...}}`; rechaza los canales con Pw = 1 (`refused_pw_on`) |
| `{"cmd":"pw","ch":N,"on":true,"authorize":"motivo"}` | `{"ok":...,"command_id":...,"result":"rejected|aborted|noop|write_failed|verified|verify_failed","detail":...}` |

  Error: `{"ok":false,"error":"..."}`. El usuario del cliente sale de `SO_PEERCRED`.
- Única escritura del código base: `HVBackend::setChannelPw()` (Pw de un canal a 0 o 1, en la
  Wrapper `CAENHV_SetChParam("Pw")`). No hay escritura de V0Set, I0Set, Trip ni otro parámetro.
- Secuencia en el hilo de sondeo: (1) rechazo si `allow_commands = false`, canal fuera de la
  configuración, sin conexión, otro comando en curso, o petición interna con un canal en rampa;
  (2) relectura inmediata del canal y aborto si Pw, Status o V0Set difieren de la última lectura
  publicada (carrera con GECO o con un trip; la ventana es una llamada de la Wrapper, ~2-3 ms);
  (3) `noop` si ya está en el estado pedido; (4) escritura; (5) verificación en las lecturas
  rápidas siguientes: encendido = Pw 1, bit `on`, sin `ramp_up` y \|VMon − V0Set\| ≤ 2 V;
  apagado = Pw 0, sin `ramp_down`, VMon ≤ 2 V (con HV apagada el A7030DP lee VMon de 0 a
  1.47 V, por eso el margen); timeout V0Set/RUp (o V0Set/RDWn) + 15 s.
- Cada fase queda en `hv_changes` como `type: command` (`request`, `rejected`, `aborted`,
  `noop`, `written`, `write_failed`, `verified`, `verify_failed`) con el motivo (`--authorize`) y
  el usuario; los cambios de Pw y Status de ese canal desde la escritura hasta 5 s después de
  terminar llevan `origin: hvmon_command` y `command_id`.
- HTTP sigue siendo de solo lectura: ningún comando HV por HTTP.

### 9.3 Prueba de grupos (fase 6)

Secuencia ejecutada por `hvmon` (no por `daqrun`), con `daqrun` en modo calibración en
paralelo:

1. Estado inicial: todos los grupos en Pw = 1, estables 10 min. Conteo de referencia
   `count_window_s` (90 s).
2. Para g = 1..8: `Pw = 0` en CH0(g−1); esperar VMon < 1 V (RDWn = 50 V/s, ~1.3 s desde
   65.9 V) más 10 s; contar `count_window_s`; `Pw = 1`; esperar fin de rampa (RUp = 1 V/s,
   ~66 s) y \|VMon − V0Set\| < tolerancia, más 60 s de estabilización; contar
   `count_window_s` (referencia posterior). Cada paso deja en `hvmon_events` sus tiempos y
   lo publica para que `daqrun` marque los eventos (bit 8 de evento) con el paso.
3. Aborto (señal, pérdida de conexión, trip en otro grupo, acción externa): se devuelve
   `Pw = 1` al grupo en prueba solo si estaba en 1 al inicio, y se termina.
4. Resultado: para cada canal de digitizador, el grupo en el que su tasa cae
   significativamente (Poisson, como 6.4) frente a las dos referencias. Llena `hv_group` con
   `confirmed = measurement`.

Duración estimada: ~4.5 min por grupo, ~40 min en total [I]. Escribe `Pw`: requiere la
autorización explícita de la sección 9.2.

## 10. Modo de calibración (trigger por software)

Requisito [H: 8.8]: runs cortos con `CAEN_DGTZ_SendSWtrigger` a tasa configurable, sin
autotrigger por canal, para la prueba de grupos y la calibración de p.e. de puesta en
marcha.

- `daqrun --mode calibration`: trigger externo deshabilitado (`external_trigger = false`), el
  hilo de lectura envía SW triggers con un temporizador a `sw_trigger_rate_Hz` (100 Hz por
  defecto) y lee entre envíos. `trigger_kind = 1` en cada evento.
- Registro largo en este modo (624 muestras por defecto) porque todo el registro es zona sin
  señal: a 100 Hz, un canal de hod2019 a ~20 kHz sobre 3 mV da ~1800 cuentas en 90 s
  [H: 8.8], y hod2018 a ~2.3 MHz sobre 0.5 p.e. da ~2.3 × 10³ pulsos/s por canal [I].
- Volumen con registros completos: 64 × 624 × 2 B × 100 Hz ≈ 8.0 MB/s [I]. Decisión de
  almacenamiento:
  1. Por evento, conteos de pulsos por canal sobre cada umbral de `pulse_thresholds_mV`
     (`UShort_t[64][n_umbrales]`, 384 B con 3 umbrales) y la línea base: ~38 kB/s.
  2. Histogramas de amplitud por canal con paso de 1 LSB, acumulados por paso de la prueba
     de grupos o por minuto: son la entrada de la calibración de p.e.
  3. Formas de onda completas de una fracción `store_waveform_fraction` (1 %) de los
     eventos, elegidos por contador (no al azar) para auditoría y reanálisis: ~80 kB/s.
  4. Listas de pulsos (tiempo, amplitud, carga) opcionales (`store_pulse_lists`), apagadas
     por defecto: con hod2018 serían ~7 × 10⁴ pulsos/s [I].
  Total por defecto ≈ 0.12 MB/s; un run de 1 h ≈ 0.43 GB [I].
- Calibración de p.e. de puesta en marcha: ajuste de los histogramas de 2. con el método de
  7.3, resultado versionado (`pe_calibration_2026_10_<n>.json`: p0, G, error, estado,
  ambigüedad, método, run de origen, voltajes de cada grupo en ese momento). Ese archivo es
  el que sale del estado "sin calibración".

## 11. Corriente estimada por MPPC

Siempre marcada "ESTIMADA".

- `identical`: I_est = (IMon − offset − I_pasiva) / N_g, con N_g = número de MPPC del grupo
  según el mapeo.
- `reference_weighted`: I_c = (IMon − offset − I_pasiva) × I_ref,c / Σ I_ref del grupo, con
  I_ref de "Voltajes Hodoscopio.xlsx" (198-278 nA por SiPM [H: 3]); la planilla no está en
  este host y se exporta a JSON cuando René la entregue.
- `common_resistance_ohm` y `series_resistance_ohm` opcionales (`null` = desconocidos, hoy).
  Con resistencia común R_c: V_MPPC,est = VMon − IMon × R_c para todo el grupo; con
  resistencia por MPPC R_s: V_c,est = VMon − IMon × R_c − I_c × R_s. El episodio del
  2026-05-29 muestra que importa distinguirlas [H: 4]; VMon e IMon a 1 Hz lo permiten.
- Offsets: `hvmon` registra IMon/VMon de todo canal con Pw = 0 y VMon < 2 V durante más de
  60 s; `hvctl offsets` hace la medida explícita de inicio/fin de campaña y se niega si
  algún canal tiene Pw = 1.
- Corriente típica de un grupo: ~2 µA, con I0Set = 15 µA (margen ~7 veces) [H: 3].

## 12. Temperatura

`TempBackend::read() → vector<{sensor_id, T_degC, t_utc_ns}>`. Implementaciones:
`Ds18b20Sysfs` (adaptador USB 1-Wire con el driver w1: `/sys/bus/w1/devices/28-*`; hoy
`/sys/bus/w1` no existe), `SerialLineSource` (microcontrolador por USB serie; René agrega
`xplus` a `dialout` cuando llegue el hardware) y `FakeTemperature`. La lee `hvmon`. Sensores
previstos: junto a cada estación y junto a las placas de trigger [H: 7.7].

## 13. Pruebas

### 13.1 Fakes

- `FakeHVBackend`: canales con V0Set/I0Set/RUp/RDWn/Trip, rampas, IMon con offset y ruido,
  guion temporal (trip, sobrecorriente creciente, cambio manual tipo GECO, pérdida de
  conexión, llamadas lentas, `Pw = 1` que no prende, re-trip en estabilización), y un modelo
  de qué canales de digitizador alimenta cada grupo para la prueba de grupos.
- `FakeDigitizer`: (a) reproducción de eventos de v8/v10 abiertos en READ desde
  `/data/daq/raw` recortados a la ventana reducida, con tiempos Poisson o copiados de
  `timestamp`; (b) sintético con pedestal, ruido, pulsos de 1 p.e. (incluido un canal con
  espectro ambiguo tipo ch9), muones, canales muertos, apagones por grupo según el mapeo,
  fuga de luz por cara, huecos y ráfagas de trigger; (c) modo SW trigger.
- `FakeTemperature` y `FakeClock` (tiempo acelerado).

### 13.2 Framework

GoogleTest de `libgtest-dev` (apt), integrado con `ctest` mediante `gtest_discover_tests`.
René lo instala antes de la fase 0; si no está instalado al empezar la fase 0, la tarea se
detiene y lo informa (no se copia un framework al repositorio).

### 13.3 Criterios de aceptación por fase

| fase | aceptación |
|---|---|
| 0 HV lib + `FakeHVBackend` + `hvctl` | tests pasan (normal, ASan+UBSan, TSan); `hvctl --fake` y `hvctl` real (solo lectura, tarea autorizada) dan la misma tabla que `hv_readall`; latencia de una lectura completa de 12 canales medida y registrada; el binario de fase 0 no importa `CAENHV_Set*` ni `CAENHV_ExecComm` (verificado con `nm -D --undefined-only`). Resultado DAQ2-2: cumplido salvo la prueba con GECO abierto (GECO estaba cerrado) |
| 1 `hvmon` + ruta de comandos | Resultado DAQ2-3: 2 h 30 min reales de solo lectura sin errores ni pérdidas (15.2), tests de todos los puntos con fake y tiempo acelerado; detección de GECO por ejecutable confirmada en DAQ2-4 (15.3); pendiente la evaluación de la suscripción. Criterios originales: 24 h con fake acelerado sin pérdida de líneas; cambios del guion con el origen correcto en el 100 % de los casos; reconexión; lector a TTree/RDataFrame validado; 8 h reales de solo lectura con GECO abierto, cada cambio hecho en GECO aparece como `external`; patrón de GECO confirmado; suscripción evaluada; comandos `pw` probados solo contra el fake (con hardware real recién en fase 2 o 6, con autorización) |
| 2 auto-retry y notificador | Resultado DAQ2-4: transiciones de 8.2 con fake y reloj acelerado (`test_hvmon_retry`), notificador con `FakeSpawner`, compresión verificada; secuencia completa NORMAL → WAIT(1..3) → GAVE_UP con fake en tiempo real y correo real (15.3). Pendiente (DAQ2-5): una vez en hardware real con carga resistiva. Criterio original: todas las transiciones de 8.2 cubiertas con tiempo acelerado; acción manual siempre manda; un correo por episodio; luego una vez en hardware real con carga resistiva e I0Set bajo, con autorización explícita: trip provocado → WAIT → RETRY → NORMAL, GAVE_UP con correo, y verificación de los bits de `Status` |
| 3 `daqrun` con digitizador fake | Resultado DAQ2-5 (15.4): 4 corridas de validación con fake y replay, sin pérdidas en la cola de escritura, segmentos e índice con SHA256 correctos, apagado limpio por señal, bits 0 y 9 del guion correctos, estado sin calibración en cabecera y `evt_flags`. Pendientes: 24 h aceleradas, pruebas de Poisson (bits 6-7, fase 5), modo calibración con almacenamiento propio (sección 10). Criterio original: 24 h aceleradas sin perder eventos; apagado limpio por señal en cualquier momento; segmentos e índice con SHA256 correctos; flags del guion correctos (incluidas pruebas de Poisson con agregación); estado sin calibración visible; modo calibración con SW trigger fake; RSS sin deriva > 10 % |
| 4 digitizador real | autorizado por el operador; `GetRecordLength` = 144; pico de muón en la muestra esperada ±1; paso del TTT medido; SW trigger a 100 Hz sin pérdidas; `DT_NEEDED` con la ruta de `/usr/lib`; 1 h de datos comparables con v8 |
| 5 monitoreo y flags | todos los paneles de 7.2; flags coherentes con `hv_changes` y con los apagones en una reproducción de v8/v10 |
| 6 puesta en marcha de octubre | con autorización explícita de escritura HV: prueba de grupos completa; prueba de coincidencias entre vecinos en la primera hora; mapeo con `hv_group` y `strip` confirmados por medición en los 64 canales (o lista de los que no se pudieron determinar); calibración de p.e. versionada |
| 7 resistencia 24 h | 24 h reales: 0 eventos perdidos sin flag, `hvmon` continuo a través de un reinicio de `daqrun`, CPU < 25 % de un núcleo, memoria estable |

## 14. Plan por fases

0. Biblioteca HV (`IHVBackend`, `CaenHVBackend` solo lectura, `FakeHVBackend`), `IClock`,
   `JThread` propio, configuración, `hvctl` de solo lectura (port de `hv_readall`) con
   medida de latencia, GoogleTest.
1. `hvmon`: sondeo 1 Hz / 1 min, slow control con su lector, detector de cambios y origen,
   estado en `/dev/shm`, página :8081, offsets, temperatura con fake, ruta de comandos
   `hvctl` → `hvmon` (probada contra el fake), servicio de usuario con
   `loginctl enable-linger xplus` (René, con sudo), evaluación de suscripción.
2. Auto-retry y notificador: con fakes y tiempo acelerado; después una vez en hardware real
   con carga resistiva e I0Set bajo, con autorización explícita para `Pw`.
3. `daqrun` con `FakeDigitizer`: hilos, colas, escritor, segmentos, índice, SHA256, cliente
   HV, flags, notas, campañas y manifiestos, modo calibración, estadística de Δt.
4. Digitizador real (autorización del operador): ventana B, TTT, contador de hardware,
   memoria llena, SW trigger.
5. Paneles de monitoreo y flags de punta a punta, ajustes de ganancia.
6. Puesta en marcha de octubre: prueba de grupos (escribe `Pw`, autorización explícita),
   prueba de vecinos, confirmación del mapeo, calibración de p.e.
7. Prueba de 24 h.

## 15. Diferencias con la propuesta y decisiones de la revisión

1. **C++17 con `JThread` propio.** ROOT 6.32.10 está compilado con C++17 y avisa con C++20
   (`RConfigure.h:30`). Decisión: aceptado; el reemplazo es pequeño, con equivalente de
   `stop_token` y cubierto por tests.
2. **Framework de tests.** No había ninguno instalado. Decisión: GoogleTest de apt con
   `gtest_discover_tests`; René lo instala antes de la fase 0; si falta, se detiene la tarea.
   No se copia un framework.
3. **Slow control en texto, no en ROOT.** Un daemon que escribe un TFile durante días pierde
   datos si muere. Decisión: aceptado, con archivos diarios y `fsync` con período
   configurable (10 s) porque `/data` es un HDD compartido; lector en la biblioteca común a
   TTree/RDataFrame. Se añade el `.schema.json` hermano para que el CSV sea legible por
   `RDF::FromCSV`.
4. **Segmentos de 1 h, AutoSave cada 5 min, SHA256 por segmento.** Limita la pérdida ante
   un fallo. Decisión: aceptado, con índice JSON por run (4.2).
5. **ADC crudo `UShort_t[64][N]`.** Sin pérdida y la mitad de bytes. Decisión: aceptado; la
   conversión ADC → mV y los DC offsets van en la cabecera del run.
6. **Secuencia, anillo de cambios y `flock`.** Decisión: aceptado.
7. **Línea de retardo de 5 s para el bit 5 y flags recalculables offline.** Decisión:
   aceptado.
8. **Temperatura en `hvmon`.** Decisión: aceptado.
9. **Mapeo JSON propio + test de consistencia con `hodo_mapping.h`.** Decisión: aceptado.
   `ch_to_plane()`/`ch_to_strip()` fijan los rangos de 2026-05 y fallan con `kCfg_2026_10`, y
   AGENTS.md contradice `hodo_mapping.h` para hod2018: ambos en `docs/notes_for_analysis.md`,
   sin tocar el código de análisis.
10. **Bit 9 `MAP_UNCONFIRMED` y flags de evento extra.** Decisión: aceptado. En esta revisión
    se añaden los bits de canal 10 (`PE_PROVISIONAL`) y 11 (`STAT_INSUFFICIENT`) y los de
    evento 7 (`NO_CALIBRATION`) y 8 (`HV_GROUP_TEST`).
11. **Sondeo con 5 ms de espera y varios eventos por BLT.** Decisión: aceptado; IRQ
    (`CAEN_DGTZ_IRQWait`) solo si el sondeo muestra un problema.
12. **THttpServer con `SetTimer(0)` en 127.0.0.1.** Decisión: aceptado como valor por defecto.
    Acceso remoto: túnel SSH por Tailscale,
    `ssh -L 8080:127.0.0.1:8080 -L 8081:127.0.0.1:8081 xplus@100.114.138.29` y abrir
    `http://localhost:8080`. Opción `bind = "tailscale"`. Nunca enlace sin autenticación a la
    interfaz de la universidad ni a `0.0.0.0` (rechazado en la validación).
13. **`daqnote` con `O_APPEND`.** Decisión: aceptado.
14. **Salvaguardas del auto-retry.** Decisión: aceptado salvo GECO: abierto no inhibe por
    defecto; aviso en el registro y detección de cambios externos (cualquier cambio de
    `Pw`/`V0Set` no hecho por `hvmon` cancela el reintento y pasa a `MANUAL`); opción
    `inhibit_when_geco_running` (por defecto false).
15. **`daqcampaign` con marcador `CAMPAIGN_OPEN`.** Decisión: aceptado; manifiestos en
    `/data/daq/manifests/` (4.4).
16. **Nuevo en DAQ2-1: diseño B de ventana.** Se descarta evitar 3510-3525 ns porque el
    pickup no aparece en v8 (sección 1.3).
17. **Nuevo en DAQ2-1: escrituras HV centralizadas en `hvmon`** (auto-retry, prueba de
    grupos, `hvctl`), con `hv.allow_commands` y autorización registrada por comando.
18. **Nuevo en DAQ2-1: estadística de Δt header-only**, para que el análisis (Makefile) la
    use sin depender de CMake.

### 15.1 Cambios fijados en la fase 0 (DAQ2-2)

19. **Credenciales y configuración en `~/.config/daqv2/`** (`hv_credentials.json`,
    `notify.json`, `hvctl.json`), no en `~/.config/daq/`. Modo 600 obligatorio, archivo del
    usuario actual.
20. **Estructura de `daq/`**: cada módulo con `include/daq/<módulo>/` y `src/`; la Wrapper se
    expone como target IMPORTED `caen::hvwrapper` (ruta `/usr/lib64/libcaenhvwrapper.so`,
    `UNIX` como definición INTERFACE). La biblioteca no tiene SONAME, así que `DT_NEEDED` de
    `hvctl` es la ruta absoluta [V: readelf -d]. `daq_hv_common` (tipos, Status, credenciales,
    fake) no depende de la Wrapper; `daq_hv_caen` sí. ROOT se busca en CMake pero ningún target
    de la fase 0 lo enlaza.
21. **Nombres**: las primitivas de `daq/core` imitan la biblioteca estándar (`request_stop`,
    `try_push`, `now_utc_ns`, `sleep_until`); las clases de dominio usan camelCase
    (`readChannels`, `crateMap`), como pidió la tarea.
22. **El header de la Wrapper define macros `uchar`, `ushort`, `ulong`** (y `ulong` es
    `unsigned int`, 4 bytes, no `unsigned long`). Solo se incluye en `caen_hv_backend.cpp` y
    las macros se anulan después. La propiedad `Type`/`Mode` de `GetChParamProp` se lee como
    `unsigned int`.
23. **Unidades de la placa**: el A7030DP informa `Unit` y `Exp` para todos los parámetros
    numéricos (IMon, I0Set, I1Set en µA, Exp = −6; V0Set, V1Set, SVMax, VMon en V; RUp, RDWn en
    V/s; Trip en s) [V: hvctl params, canales 0 y 8]; la tabla de respaldo por nombre no se usó.
    La conversión a SI divide por 10^−Exp (15 µA → 1.5e-05 A exacto) y el float de la Wrapper se
    pasa a double por su forma decimal de 7 cifras (0.048 y no 0.0480000004).
24. **Crate**: 4 slots (0-3), solo el 3 ocupado: A7030DP, descripción
    `" 12 Ch Pos 3KV 1mA (1.5W) "`, SN 58015, FW 1.15. ON/OFF: Pw On/Off, POn En/Dis, PDwn
    Ramp/Kill, ZCDetect On/Off, ZCAdjust En/Dis. Con PDwn = 0 la placa está en **Kill**; el fake
    interpreta Kill como caída inmediata al apagar y Ramp como bajada a RDWn [I: a verificar en
    fase 2]. La secuencia de la prueba de grupos (9.3) no depende de esto: espera VMon < 1 V.
25. **CH08-11** no tienen los mismos ajustes que CH00-07: I0Set 5 µA, I1Set 101 µA, RUp 5 V/s,
    Trip 10 s [V: hv_readall 2026-10-18T21:56:18Z].
26. **Latencia de la Wrapper** (hvctl bench, 100 lecturas de 6 parámetros × 12 canales,
    2026-10-18T22:10:42Z, HV apagada, GECO cerrado): mediana 1.954 ms, p95 2.772 ms, máximo
    3.145 ms; repetido con el binario final a las 22:21:38Z: mediana 2.063 ms, p95 2.860 ms,
    máximo 6.208 ms. Un sondeo de 1 Hz usa ~0.2 % del tiempo; el modo suscripción no hace falta por
    latencia (se evalúa en la fase 1 solo por detección de cambios).
27. **Sanitizers**: ASan+UBSan combinados (`DAQ_SANITIZE=address,undefined`) y TSan
    (`thread`) por separado. TSan de GCC 13 aborta con "unexpected memory mapping" con la
    entropía ASLR de este kernel; los tests TSan corren con `setarch x86_64 -R` (ASLR apagado
    solo para ese proceso) mediante `CROSSCOMPILING_EMULATOR`. Los sanitizers reservan
    terabytes de espacio virtual y no funcionan con `ulimit -v 4194304`: se corren dentro de
    `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0`, que limita la memoria real
    a 4 GiB.

### 15.2 Fase 1 (DAQ2-3): decisiones fijadas y validación

Decisiones fijadas por la implementación:

28. **Secretos** en `~/.config/daqv2/secrets.env` (o el entorno; `EnvironmentFile=` en la unidad
    de systemd); el JSON heredado `hv_credentials.json` solo como respaldo con aviso. Ningún
    secreto en argv, logs, estado, HTTP ni salida de tests (test de redacción de `*PASSWORD*`).
    `hvmon` y `hvctl` ponen `RLIMIT_CORE = 0` (sin volcados de memoria); la unidad, `LimitCORE=0`.
29. **Rutas**: lock y socket en `$XDG_RUNTIME_DIR/daqv2/` (`hvmon.lock`, `hvmon.sock`), estado en
    `/dev/shm/daqv2/hvmon_state.json`; el anillo de cambios es de 50 (antes se proponían 256).
30. **Archivos de control lento** con subdirectorio por año y nombres `hv_`, `hv_full_`,
    `hv_changes_`, `notes_`, `temp_` (sección 4.3). Los eventos de conexión y GECO van al log de
    `hvmon`, no a un `hvmon_events` separado (se agregará con el auto-retry de la fase 2).
31. **Orígenes** `trip`, `hvmon_command`, `external` (antes `daq_auto`/`daq_command`). `trip` solo
    para cambios de Pw o Status con bit de trip en el Status nuevo; los cambios de Status por
    sobrecorriente o rampa sin comando quedan como `external`. En las lecturas completas se siguen
    todos los parámetros no monitores (superconjunto de RUp, RDWn, Trip, I1Set, V1Set).
32. **Verificación de Pw** por Status (bit `on`, sin bits de rampa) con VMon como control grueso
    (2 V): con HV apagada la placa lee VMon de 0 a 1.47 V. La tolerancia de 0.5 V del borrador
    nunca habría verificado con el offset de CH03 (0.78 V); lo mostró el test con el fake.
33. **JSON de lecturas** con `params` como arreglo ordenado (en la fase 0 era un objeto, que en
    nlohmann pierde el orden); `hvctl read --json` cambia en consecuencia.
34. **Parámetros de placa** (pedido de René): en cada lectura completa, `GetBdParamInfo` +
    `GetBdParamProp` + `GetBdParam` (solo lectura). La placa informa 5: `BdStatus` (bdstatus),
    `HVMax` 3127 V, `HIMax` 1111 µA (Exp −6), `Temp` 25.000 degC (unidad 4) y `ClrAlarm` (on/off,
    lectura-escritura) [V: hvctl board, 2026-10-19T01:53:50Z]. La lista de nombres no trae cantidad:
    bloques de `MAX_PARAM_NAME` terminados en un nombre vacío (comprobado: 5 nombres). `ClrAlarm` es
    el comando de borrado de alarmas: queda excluido de las lecturas periódicas
    (`board_params_exclude`) y es relevante para la pregunta de cómo rehabilitar tras un trip (fase
    2). Columnas `bd_BdStatus, bd_HVMax_V, bd_HIMax_A, bd_Temp_degC` en `hv_full` (esquema
    `daqv2.hv_full/2`); `board` y `board_temperature_degC` en el estado, con la nota de que es la
    temperatura de la placa A7030DP (sigue la del rack), no la del detector; gráfico de 24 h en
    `/board/BoardTemp_degC` y línea en la tabla de la página; `hvctl board` en modo directo.
35. **Versión de la Wrapper**: `CAENHVLibSwRel()` devuelve "6.5" aunque el archivo es
    `libcaenhvwrapper.so.6.6` [V: log de hvmon].
36. **Fuente del SW trigger** recuperada de un árbol colgante de git (sección 1.2).
37. **Lector a TTree/RDataFrame** de los archivos de control lento: pospuesto a la fase 3.

Validación contra el SY5527 real (solo lectura, `allow_commands: false`, HV con Pw = 0 en todos
los canales; `hvmon` del commit `163abcd`, binario normal; números copiados de su estado final y
del muestreo cada 60 s):

| magnitud | valor |
|---|---|
| duración | 8999.96 s (2026-10-18T23:22:08.49Z a 2026-10-19T01:52:08.46Z, SIGTERM por `timeout`) |
| lecturas rápidas / completas | 9000 / 150 |
| errores de lectura / intentos de conexión / reconexiones | 0 / 1 / 0 |
| duración de lectura rápida (12 canales × 6 parámetros) | media 3.080 ms, máximo 7.198 ms |
| duración de lectura completa (12 × 18) | media 8.591 ms, máximo 9.362 ms |
| CPU (`ps -o %cpu`, promedio desde el arranque) | 1.1 % a los 44 s; 0.7 % al final |
| RSS | 249176 kB al inicio, 256088 kB a los 15 min, 256728 kB al final (+640 kB en las últimas 2 h, en pasos de 128 kB) |
| escritor | 109807 líneas, 9581534 bytes, 1045 fsync, 0 errores, 0 descartes |
| archivos | `hv_20261018.csv` 2355361 B (27265 líneas), `hv_20261019.csv` 7009223 B (80737), `hv_full_20261018.csv` 54313 B (457), `hv_full_20261019.csv` 160328 B (1345), `notes_20261018.jsonl` 373 B (2), `offsets.jsonl` 1936 B (1) |
| filas | 108002 = 9000 × 12 + 2 cabeceras; 1802 = 150 × 12 + 2 |
| rotación | a 00:00:00Z: el primer dato del 19 fue a `hv_20261019.csv` |
| cambios | 0 (no se creó `hv_changes_*.jsonl`) |
| estado | `stale` = 0 en las 140 muestras; estado final con `final: true` |
| segundo `hvmon` | rechazado por el lock, código de salida 3 |
| `hvctl read` / `status` por socket | misma tabla que `hv_readall`; `map` rechazado mientras corre `hvmon` |
| offsets (60 s, 62 lecturas, ningún canal rechazado) | IMon CH00-11 = −0.0480, −0.0360, −0.0460, −0.0360, −0.0440, −0.0460, −0.0460, −0.0740, −0.1180, −0.0940, −0.0740, −0.0840 µA, desviación 0.0000; igual a la lectura de la fase 0. VMon 0.014 a 1.476 V |
| HTTP | 16 `TGraph` (VMon, IMon de CH00-07), `/status/table`, `/status/state_json`; solo 127.0.0.1:8081; `exe.json` responde 404 |

GECO (René abrió GECO, notas 23:29:14Z "GECO abierto" y 23:49:12Z "GECO cerrado"; conectado al
SY5527 desde ~23:32Z, sin operar nada):

- `hvmon`: "GECO2020 detectado en ejecución" 23:30:59Z, "ya no está en ejecución" 23:49:10Z; entre
  23:27 y 23:52 ninguna muestra con errores, reconexiones ni `stale`; la duración máxima de lectura
  pasó de 3.339 a 3.423 ms. Ningún cambio en el registro.
- Falso positivo: una segunda detección de 23:55:50Z a 23:57:01Z con GECO cerrado. Causa [I]: el
  detector de entonces buscaba el patrón en toda la línea de comandos, y un proceso `bash -c` del
  agente, cuyo texto contenía `bin/CAENGECO2020`, compiló y corrió los tests de `build-dev` en esos
  ~71 s (la compilación siguiente empezó a las 23:57:12Z); GECO no reescribió su log de suscripción
  (sigue con 23:32:13Z). El bucle de muestreo, cuya línea contiene `geco_running` pero no el
  patrón, no coincidió (el indicador volvió a 0 con el bucle corriendo). Corregido en `6eb49b2`:
  solo `/proc/<pid>/exe` o `argv[0]`. Detección nueva confirmada con GECO abierto en DAQ2-4 (15.3).
- Escribir por HTTP o con GECO no se probó (fuera del alcance de la fase 1).

### 15.3 Fase 2 (DAQ2-4): decisiones fijadas y validación

Decisiones fijadas por la implementación:

38. **Estados**: `RETRY_RAMP`/`RETRY_STAB` se llaman `RETRY`/`STABILIZING`; `UNKNOWN` se
    reemplaza por "sin lectura fresca → no avanza ni escribe"; `GROUP_TEST` queda para la
    fase 6. Tiempos en minutos (`wait_min`, `stable_min`) y horas (`reset_after_h`).
39. **`DISABLED` silencioso**: los canales fuera de `enabled_channels` no generan registro al
    arrancar (una versión previa registraba `UNINIT -> DISABLED` por canal).
40. **Trip nuevo** en RETRY/STABILIZING: Pw 1→0 con origen `trip`, o bit de trip con Pw = 1.
    Un bit de trip retenido con Pw = 0 no cuenta.
41. **Fase `busy`**: si la ruta de comandos está ocupada, la solicitud interna del auto-retry
    vuelve a `WAIT(k)` sin contar intento. `write_failed` igual; `verify_failed` cuenta como
    intento fallido; `aborted` y `rejected` pasan a `MANUAL`.
42. **Segunda escritura** del código: `clearBoardAlarm()` (`CAENHV_SetBdParam("ClrAlarm")`),
    solo por el poller, con la misma compuerta `allow_commands` y su propio `CommandRecord`;
    se usa solo con `clear_alarm_before_retry = true` (false por omisión). Con la de `Pw` son
    las dos únicas escrituras.
43. **Persistencia del auto-retry** en `$XDG_STATE_HOME/daqv2/hvmon_retry_state.json` (no en
    `/dev/shm`, que se pierde al reiniciar), escrita atómicamente en cada transición; una
    espera persistida se reanuda solo si no venció hace más de 10 min.
44. **Notificador** sin `notify.json`: SMTP en `secrets.env`, `curl --config -` por stdin,
    cuerpo en un temporal 0600 (8.4). Los correos de prueba saltan límites y reintentos.
45. **Compresión** en `hvmon` (hilo con nice 19) y herramienta `sctool` con el mismo código;
    verificación por bytes, líneas y CRC32 antes de borrar el plano (8.5).
46. **Mapeo** `daq/config/mapping/hodo_mapping_2026_10.json` (`daqv2.mapping/1`) para el texto
    de los correos; los grupos de hod2018 no medidos van marcados como supuestos o
    desconocidos.

Validación (DAQ2-4):

- Tests: 132/132 en las tres builds (normal, ASan+UBSan, TSan), sin reportes de sanitizers
  (tabla en el informe de DAQ2-4).
- Defecto encontrado en la corrida real y corregido: el escritor agregaba filas de 28 columnas a un
  `hv_full` del mismo día con encabezado de 24 (4.3, "Partes"); tests
  `HvmonWriter.ExistingFileWithOtherHeaderStartsNewPart` y `MissingBoardReadingNeverMixesColumns`.
- Compresión de copias de `/data/daq/slowcontrol/2026` con `sctool`: total 11628312 → 789041
  bytes (0.0679); `hv_*.csv` 0.0686/0.0680, `hv_full_*.csv` 0.0598/0.0552, `notes` 0.6059.
  Luego los originales del día cerrado 20261018 con la misma herramienta: 2410047 → 165164
  bytes (0.0685), digests (bytes, líneas, CRC32) idénticos antes y después.
- Correo real: `hvctl notify-test` → "curl salió con 0 (intentos 1)". Secuencia con fake
  (CH08 a 10 V, I0Set 5 µA, Trip 10 s, excursión de +20 µA continua, `test_mode` 1/2/3 min,
  `stable_min` 1): MANUAL → NORMAL (09:19:52Z) → trip 09:20:37Z → WAIT(1) → RETRY(1)
  09:21:38Z → STABILIZING(1) → trip → WAIT(2) → RETRY(2) 09:23:50Z → … → RETRY(3) 09:27:02Z →
  STABILIZING(3) → trip → GAVE_UP 09:27:13Z; correo `gave_up CH8` a las 09:27:17Z, curl 0;
  4 eventos `trip` descartados por `notify_each_trip = false`. Dos correos en total.
- Un primer intento con la excursión a t = 20 s no inició el episodio: el canal estaba en
  `MANUAL` (apagado al arrancar) y el trip llegó antes de `stable_min`; es el comportamiento
  de 8.2.
- GECO con el detector nuevo (`/proc/<pid>/exe` o `argv[0]`), corrida de solo lectura desde
  08:32:37Z: nota "GECO abierto" 09:48:36.798Z → "GECO2020 detectado en ejecución" 09:48:55.626Z;
  nota "GECO cerrado" 10:01:35.983Z → "ya no está en ejecución" 10:01:36.321Z. Una sola
  detección, sin falsos positivos. Al final: 5400 lecturas rápidas y 90 completas, 0 errores,
  0 reconexiones, 0 cambios, escritor 0 errores. `board_temperature_degC` = 25.0 en el estado,
  en la página (`/board/BoardTemp_degC`) y en `bd_Temp_degC` de `hv_full`.

### 15.4 Fase 3 (DAQ2-5): decisiones fijadas y validación

Parte A (correcciones de DAQ2-4):

47. **Marca de simulación en los correos**: `[SIMULADO]` al inicio del asunto y primera línea del
    cuerpo "*** SIMULADO *** Evento simulado (FakeHVBackend | test_mode del canal CHnn): ningún
    canal real estuvo involucrado." para hvmon con fuente fake o canal en `test_mode`;
    `[PRUEBA]` para `notify-test`. Cuerpo: IMon en µA con 3 decimales, VMon en V con 2, hora
    `America/Santiago` (`notifier.local_timezone`, TZ del proceso fijada antes de crear hilos)
    junto a cada hora UTC, espera aplicada en cada `retry` de la secuencia (`wait_min`), y línea
    `ssh -N -L <puerto>:127.0.0.1:<puerto> <page_ssh_user>@<page_ssh_host>` con el puerto de
    `http.port`.
48. **Salida simulada separada**: hvmon con `--fake` o `test_mode` rechaza `slowcontrol_dir`
    dentro de `/data/daq/slowcontrol` (ruta resuelta, con enlaces); con fake también el
    `state_dir` real. Valores por omisión: `/data/daq/sim/slowcontrol`, `/dev/shm/daqv2-sim`,
    `$XDG_RUNTIME_DIR/daqv2-sim`. Las corridas fake de DAQ2-4 escribieron en el scratchpad de la
    sesión (copias en `/data/daq/sim/archive_daq2-4/`), nunca en `/data/daq/slowcontrol`.
49. **`.gz` con el modo y dueño del original** (gzopen creaba 0666 & umask = 664).
    hvmon publica además `last_read` (`read_seq`, `kind`, `t_req_utc_ns`, `t_resp_utc_ns`).

Fase 3:

50. **Bibliotecas**: `daq_dig` (interfaz, `FakeDigitizer`, `ReadoutLoop`; sin ROOT),
    `daq_dig_replay` (ROOT, READ), `daq_daqrun` (configuración, campañas, cliente HV y flags,
    monitor, notas; sin ROOT), `daq_io` (escritor ROOT), `libdaq_screader.so` (lector de control
    lento para macros). Compartidos con el análisis, header-only y sin ROOT:
    `daq/core/include/daq/pulse_finder.h` y `daq/core/include/daq/trigger_interval_stats.h`.
    Tests con ROOT en un binario aparte (`daq_root_tests`); LeakSanitizer suprime solo
    `libCling.so` en ese binario (cling pierde 227 asignaciones al salir también en un programa
    mínimo con un TFile).
51. **Bloques compartidos**: `shared_ptr<const EventBlock>` a la cola de escritura (push
    bloqueante, contado) y a la del monitor (descarte del más antiguo); el diseño hablaba de copia.
52. **Línea de retardo**: un evento se escribe cuando la hora del escritor pasa t_utc +
    `flag_delay_s`; al cerrar el run se espera `flag_delay_s` + 1 s con lecturas HV nuevas y sin
    rotar segmento.
53. **Bits HV** (FlagEngine, hilo escritor): lecturas deduplicadas por `read_seq`; bit 5 en
    (t_req(a) − backdate, t_resp(b)] de cada par consecutivo visto por daqrun (si daqrun saltea
    una lectura el intervalo crece: conservador); bit 4 sin lectura posterior al liberar el
    evento, con hvmon no vivo/stale, lectura previa con más de `stale_after_s` o grupo
    desconocido; bit 1 = bit `overcurrent` o IMon > `ovc_k_over_baseline` × mediana de IMon en
    NORMAL de la última hora (≥ 60 lecturas; sin offset todavía); bit 9 estático del mapeo (con
    el JSON actual, todos los canales salvo ch48). `evt_flags`: 0 y 1 con la tasa de la última hora
    (≥ 20 triggers; ráfaga con μ = 5 en la ventana), 2 (contador de 24 bit), 3, 4 (|t_utc − t_read| >
    `ttt_pc_tolerance_s` = 1 s), 5 y 7. El resto en 0 (fases 5-6). hvmon ausente: bit 4 en todos
    los canales y aviso en la página.
54. **TTT**: t_utc = UTC antes del start + ttt_ns, desenrollado con la hora de lectura
    (`ttt_wrap_bits` = 31 por omisión: 17.18 s con 8 ns) y sin retroceder. Replay: TTT al ritmo de
    entrega (intervalos originales / speed; speed 0 = hora de lectura): con speed 0 los eventos de
    un bloque comparten tiempo y los bits 0/1 no significan nada (usar speed 1 para tiempos).
55. **Archivo**: ZSTD nivel 5 (`compression` = 505); `adc[64][N]/s`; árbol `hv` con
    received/state_seq/read_seq/t_req/t_resp/file_ok/alive/stale/connected y arreglos [12]
    (`retry_state` como `vector<string>`; la estimación de corriente es de la fase 5);
    `hv_changes` con old/new en texto JSON; `monitor` cada `summary_period_s` (60 s) con tasas,
    colas y eventos con cada bit en el período. Cabecera por segmento (TNamed): `config_json`,
    `config_sha256`, `mapping_json`, `mapping_version`, `mapping_sha256`, `calibration_json` y
    `calibration_version` ("none"), `run_info_json` (con `digitizer` al arrancar y
    `digitizer_final` al terminar la lectura), `hvmon_state_at_start`.
56. **Junto a los datos**: `<base>.config.json` (configuración resuelta, 0444), segmentos 0444,
    `<base>.index.json` (reemplazo atómico; 0444 al terminar), `SHA256SUMS` (flock), marcadores
    `CAMPAIGN_OPEN` → `CAMPAIGN_CLOSED` (rename) y directorio 0555 al cerrar; lock `.run.lock` (un
    run por campaña); número de run = 1 + el mayor encontrado. Manifiesto por campaña
    `<manifests>/<campaña>.json` con los runs y el SHA256 de `SHA256SUMS`. Campaña real:
    `daqcampaign` exige `--authorize`.
57. **Simulación**: digitizador fake o replay → `raw_root` y `manifests_dir` nunca dentro de
    `/data/daq/raw` ni `/data/daq/manifests`; por omisión `/data/daq/sim/raw` y
    `/data/daq/sim/manifests`. `hodo::resolveRawFile` encuentra los segmentos simulados solo con
    `HODO_RAW_DIR=/data/daq/sim/raw`.
58. **Monitor de la fase 3**: umbrales provisionales hod2018 20 mV y hod2019 3 mV; zonas fuera de
    tiempo pedestal + control + post (+ tail si el canal no tiene pulso en tiempo); espectro desde
    2 LSB hasta 100 mV con paso de 1 LSB; Δt en log10(Δt/s) de −4 a 3.
59. **Commit del build** en cada build (`daq_build_info.h`); `git_dirty` mira solo `daq/`.
60. El mapeo JSON tiene 23 canales de hod2018 sin grupo (5.1 decía 22).

Validación (DAQ2-5; números copiados de la salida):

- Builds: compilación limpia (`--clean-first`) de las tres variantes, con el código de la parte A
  (0 warnings, 51 unidades compiladas en cada una) y otra vez con el código final de la fase 3 (0
  warnings, 70 unidades en cada una). Tests con el código final: 155/155 en normal, ASan+UBSan
  (0 reportes) y TSan (0 reportes).
- Validación 1 (sintético 0.6 Hz, 1 h, hvmon real de solo lectura): 2208 eventos en 2204 bloques
  (span 3599.304712608 s), 1 segmento de 12215430 bytes, SHA256 correctos, índice completo; colas:
  escritura máx. 1, 0 bloqueos, monitor 0 descartes; CPU 53 s en 3602 s; RSS 263672 kB al
  arrancar, 296056 kB a los 5 min (primer AutoSave), 300148 kB al final (+4092 kB en 55 min).
  Árbol hv: 343 filas, todas con hvmon vivo, las 343 con su read_seq en `hv_`/`hv_full_20261019`
  y 0 distintas (Pw, Status, V0Set, VMon, IMon de los 12 canales; diferencia máxima 0 V,
  1.32e-23 A).
- Defecto encontrado en la validación 1 y corregido: el lector de control lento dejaba las ramas
  apuntando a variables locales; un macro que ponía dirección solo a algunas ramas escribía en
  memoria liberada (violación de segmento). Ahora `ResetBranchAddresses()` antes de devolver el
  árbol, con test.
- Correo: `hvctl notify-test` desde el hvmon real (solo lectura) → "curl salió con 0 (intentos 1)".
- Validación 2 (sintético 200 Hz, 10 min, segmentos de 2 min, hvmon real): 120217 eventos en
  79707 bloques (máx. 16 eventos por bloque), 5 segmentos contiguos de 120 s (22870, 24552,
  24278, 24449, 24068 eventos), 621476986 bytes, SHA256 correctos, índice completo; cola de
  escritura: 0 bloqueos, ocupación máxima 228 de 256; monitor: 0 descartes, máx. 2; CPU 346 s en
  605 s; RSS 267416 kB al arrancar, 396416 kB al final. Entrada ADC 2215839744 bytes (3.69 MB/s),
  escrito 1.04 MB/s (span de eventos 600.267 s).
- Validación 3 (replay de las entradas 30000-34999 de v8, speed 0, run 2 con el TTT corregido):
  5000 eventos en 313 bloques, 22101516 bytes, RSS máximo 661300 kB. 100 eventos al azar:
  921600 muestras, 0 fuera de LSB/2, diferencia máxima 0.000014 mV, contador = `event` en todos;
  residuo máximo de redondeo 6.242741153528186e-05 ADC. Tasas fuera de tiempo en esos 5000
  eventos (91 muestras): hod2018 > 3 p.e. (ganancias de v8) mediana 377.75 kHz con el buscador
  nuevo y 378.85 kHz con el método de DAQ2-1 (31 canales); hod2019 > 3 mV 29.67 y 29.26 kHz (32
  canales). DAQ2-1 midió 331.4 y 24.27 kHz con otros 9057 eventos repartidos en el run.
- Validación 4 (sintético 5 Hz, hvmon simulado con CH02 apagado 299.0 s según el árbol hv):
  2743 eventos; 1457 eventos con CH02 apagado lejos de los cambios, bit 0 exactamente en los
  canales 33, 38 y 49 (grupo 3 del mapeo) en todos; 1119 con todo encendido, ningún bit 0; bit 9
  distinto de lo esperado en 0 casos (63 canales sin confirmar); bit 4 distinto de lo esperado en
  0 casos (23 canales sin grupo); bit 5 en 24 eventos, todos a ±8 s de un cambio y solo en 33, 38,
  49.
- Apagado por SIGINT (dos señales) a mitad de un run de 20 Hz: salida 0, índice
  `complete: true`, `closed_cleanly: true`, SHA256 correctos.
- Pendiente respecto del criterio de 13.3: 24 h aceleradas, pruebas de Poisson (bits 6-7, fase
  5) y el almacenamiento del modo calibración (sección 10; el fake y el bucle ya envían SW
  triggers).

## 16. Preguntas resueltas por el handoff y la revisión

| pregunta (DAQ2-0) | respuesta |
|---|---|
| mapeo de la campaña nueva | `kCfg_2026_10` [H: 5] |
| grupos HV y relación con CH00-07 | grupo g = CH0(g−1), confirmado; tabla de etiquetas y tabla derivada de hod2019 en 5.1 [H: 3, 5] |
| caras de hod2019 | Delantera, Trasera, Opuesta, Frontal por regla de 5.1 [H: 5]; hod2018 sin caras para fuga de luz (la luz llegó también a hod2018 en mayo [H: 4], se agrega por grupo o canal) |
| planilla "Voltajes Hodoscopio.xlsx" | contenido relevante en [H: 3]; la planilla misma no está en el host |
| ventana A o B | B, con datos (1.3) |
| estadística de 1 p.e. con ventana reducida | ~7400 pulsos > 0.5 p.e. por canal y hora en hod2018 (1.3); hod2019 con agregación o modo calibración |
| TTT | no hay manual en el host; medida en la fase 4 (2.3) |
| bits de `Status` | posiciones en configuración; valores iniciales del manual del Wrapper; verificación en fase 2 |
| latencia y suscripción | latencia en fase 0, suscripción en fase 1 |
| nombre de GECO | snap `caengeco2020`, ejecutable `bin/CAENGECO2020`; patrón confirmado en fase 1 |
| servicio | servicio de usuario con `enable-linger` (fase 1) |
| acceso web | localhost + túnel SSH por Tailscale; opción Tailscale |
| SMTP | `secrets.env` (600, fuera del repo), `DAQ_SMTP_*` y `DAQ_NOTIFY_*`; `curl --config -` por stdin (8.4) |
| temperatura y `dialout` | `TempBackend` + fake; `dialout` cuando llegue el hardware |
| resistencias de la PCB | desconocidas; opcionales en configuración |
| framework de tests | GoogleTest de apt |
| manifiestos | `/data/daq/manifests/` |
| proceso `agy` | René lo cierra o quita el repositorio de su espacio de trabajo |

## 17. Preguntas abiertas

1. Si los cables de hod2018 se tocaron entre mayo y octubre (define si valen las relaciones
   de grupo medidas en mayo) [H: 10].
2. Orden de las etiquetas "Canal" del eje Y de hod2018 y orientación absoluta de los ejes
   [H: 10]; los resuelven la prueba de grupos y la de vecinos solo en parte.
3. Modelo de MPPC y voltaje de ruptura (afecta la interpretación de G_c(t) y una eventual
   compensación de bias por temperatura) [H: 10].
4. Manual del A7030DP / SY5527: bits de `Status` específicos y procedimiento para rehabilitar
   tras un trip (René lo busca).
5. Canales de hod2019 con tasa oscura ~0 (en v8 el mínimo fue 0 kHz sobre 3 mV): la prueba de
   grupos no puede asignarles grupo por tasa oscura. ¿Se acepta asignarlos por la planilla o
   se usa otra señal (p. ej. un LED o la línea base)?
6. Si el pickup reaparece en octubre con otro cableado; el monitor lo vigila (7.1).
7. Proveedor SMTP y destinatarios (antes de la fase 2).
