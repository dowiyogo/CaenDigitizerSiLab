# Checklist de laboratorio, lunes 2026-11-02

Para la persona en el laboratorio del hodoscopio (y para René o el agente, que corren los comandos de comprobación desde el host
`xplus`). Escrito el 2026-10-31.

## Contexto

- **Física detenida:** no se toman datos desde el 2026-10-30 a las 02:12 UTC.
- **El enlace con el digitalizador está colgado:** el VX1740, conectado por fibra CONET al adaptador USB A4818, no responde desde las 02:03:53 UTC del 30 (`CAEN_DGTZ_CommError` y luego `CAEN_DGTZ_Timeout`).
  - El A4818 sigue visible por USB.
  - Un corte igual del 2026-10-23 solo se arregló con un apagado.
- **HV de los MPPC:** encendida en valores nominales (SY5527, CH00-07). `hvmon` la vigila.

## Lo que NO se hace en ningún paso

- **No tocar V0Set ni Pw** de ningún canal de la SY5527: ni desde GECO, ni desde la fuente, ni desde la línea de comandos. Solo se toca RUp en el paso d, y solo si René lo decide.
- **No abrir las cajas de los hodoscopios** ni retirar tapas o cinta. Con HV encendida, la luz ambiente en un MPPC dispara la corriente y el trip, y puede dañarlo.
- **No apagar la SY5527** ni el PC `xplus`.
- **No desconectar cables de señal** de los hodoscopios ni los cables de las tejas de trigger.
- **No mover las tejas de trigger** ni los hodoscopios: su posición es la referencia de los datos. Si hay que moverlos para la prueba de orientación, avisar a René antes.
- **No cerrar CAENGECO2020** en el PC (puede estar abierto) y no cambiar en él nada fuera del paso d.

---

## a) Antes de tocar nada

1. **Fotos del montaje** desde arriba y de lado: hodoscopios, tejas T1 y T2, crate VME con el VX1740, crate NIM, A4818 y fibra CONET. Que se vean las etiquetas y los cables. Enviarlas a René.
2. **Temperaturas:**
   - Crate VME: si la fuente o el panel muestran temperatura, anotarla. Si no, un termómetro de ambiente junto a la entrada de aire.
   - VX1740: si el panel frontal tiene indicador de temperatura o LED de alarma, anotar su estado.
   - Anotar la hora (hora de Chile) de cada lectura.
3. **LEDs del frente del VX1740 y del A4818:** cuáles están encendidos, apagados o parpadeando, con foto.

Riesgo: ninguno; solo observación.

## b) Fibra CONET y reinicio del enlace

### b.1 Revisión visual (sin desconectar)

Revisar la fibra entre el A4818 y el VX1740:
- que no tenga curvas cerradas ni esté aplastada;
- que los conectores estén bien asentados en las dos puntas;
- que no estén sueltos ni tirantes.

Foto de los dos conectores. Si se ve un conector suelto, **no** limpiar ni reinsertar todavía: avisar a René.

Riesgo: bajo. Desconectar la fibra no daña nada, pero se pierde la información de cómo estaba.

### b.2 Comprobación antes de reiniciar

René o el agente, en el host:

    lsusb | grep 21e1
    fuser -v /dev/bus/usb/002/* 2>&1 | head         # que ningún proceso tenga tomado el A4818
    pgrep -af -i 'daqrun|digprobe|wavedump'          # que no corra ningún DAQ
    ~/.local/opt/daqv2/bin/digprobe info --config /data/daq/config/hodoscope_2026_10/physics_ext_B.json

Hoy se espera:
- `Bus 002 Device 003: ID 21e1:0016 CAEN SpA A4818 USB3/CONET2 ADAPTER`;
- que `digprobe info` falle con `CommError` o `Timeout`.

`digprobe info` abre el VX1740 (solo lectura de registros): lo corre René o el agente con una tarea que lo autorice.

### b.3 Paso 1: desconectar y reconectar el A4818 (USB)

1. Confirmar con René que no corre ningún DAQ (b.2).
2. Desconectar el cable USB del A4818 en el lado del PC, esperar 10 s y reconectarlo.
3. Comprobación en el host:

        lsusb | grep 21e1                      # debe reaparecer (el número de Device cambia)
        journalctl -k --since "-5 min" | grep -i -E 'usb|21e1|caen'
        ~/.local/opt/daqv2/bin/digprobe info --config /data/daq/config/hodoscope_2026_10/physics_ext_B.json

4. **Resultado bueno:** `digprobe info` muestra `VX1740`, serie 21, ROC 04.29 y los registros, sin `ERROR`. Si es así, saltar al paso e.

Riesgo: bajo. Solo si el DAQ estuviera corriendo se perdería el run, y no está corriendo.

### b.4 Paso 2, solo si b.3 no basta: apagar y prender el crate VME

1. Confirmar con René. En el crate VME está **solo** el VX1740. La HV de los MPPC va en la SY5527, que **no** se apaga.
2. Apagar el crate VME con su interruptor, esperar 30 s y prenderlo. Esperar 1 min a que el VX1740 termine de arrancar (LEDs).
3. Si el A4818 se desconectó, reconectarlo como en b.3.
4. La misma comprobación de b.3: `lsusb`, `journalctl -k` y `digprobe info`.
5. Si `digprobe info` sigue fallando, **parar** y avisar a René con fotos de los LEDs. No probar otros cables ni otras ranuras.

Riesgos:
- Bajo para el detector, porque el crate VME no alimenta los MPPC.
- Apagar el **crate NIM** en vez del VME detendría el trigger: comprobar la etiqueta del crate antes de apagar.

## c) Prueba de orientación con fuente β (preferida) o LED en la esquina (0,0)

Objetivo: saber qué extremo del orden de tiras de hod2018 corresponde al origen (0,0) de René. La macro ya existe: `analysis/followup_2026_10_31/orientation_corner_test.C`. Se hace después de b, con el enlace funcionando y **antes** de volver a arrancar la física (paso e).

### Convención de René

Visto desde arriba:
- el origen (0,0) es la esquina inferior izquierda;
- las barras X son paralelas al eje Y (verticales en el dibujo) y miden x;
- las barras Y son paralelas al eje X (horizontales) y miden y.

Lado de lectura de cada capa:
- capa X superior: abajo;
- capa X inferior: arriba;
- capa Y superior: derecha;
- capa Y inferior: izquierda.

          y ↑                MPPC de la capa X inferior (arriba)
            |      ┌──────────────────────────────────────┐
            |      │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │      │
    MPPC    |      │─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─│      │   MPPC de la
    capa Y  |      │ │ │ │  área activa 33 × 33 mm  │      │   capa Y
    inferior|      │─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─┼─│      │   superior
    (izq.)  |      │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │      │   (der.)
            |      │ ★ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │      │
            |      └──────────────────────────────────────┘
            |      ★ = fuente, sobre la ESQUINA del área activa
            +——————— MPPC de la capa X superior (abajo) ——————→ x
          (0,0)

### Cómo colocar la fuente

- **Dónde:** sobre la tapa de la caja de hod2018, **sin abrirla**, encima de la esquina (0,0) del área activa, a unos 2-3 mm de los dos bordes del área activa (no del borde de la caja). Si no se ve el área activa desde fuera, medir desde las marcas de la caja o pedir a René la posición.
- **Colimación:** fuente β (p. ej. Sr-90) con colimador, apuntando hacia abajo. Fijarla con cinta para que no se mueva.
- **LED:** solo si existe un puerto óptico ya previsto. **No** abrir la caja para meter un LED.
- **Seguridad radiológica:** manipular la fuente según las reglas del laboratorio, solo personal autorizado, y anotar fuente, actividad y hora.

### Runs

Los lanza René o el agente, con una tarea que autorice el digitalizador. La física sigue detenida. Cada run dura **300 s** con SW trigger a 100 Hz (unos 30 000 eventos):

    # campaña de la prueba
    ~/.local/opt/daqv2/bin/daqcampaign --authorize "orientación hod2018, autorizado por René" open orientation_2026_11_02
    # 1) referencia SIN fuente (la fuente lejos, a > 1 m)
    systemd-run --user --unit=orient-ref --collect --property=RuntimeMaxSec=900 \
      --setenv=LD_LIBRARY_PATH=$HOME/root6/lib --working-directory=$HOME \
      ~/.local/opt/daqv2/bin/daqrun --config /data/daq/config/biasscan_2026_10_30_night/sw.json \
      --campaign orientation_2026_11_02 --max-duration 300 --no-http
    # 2) con la fuente en la esquina (0,0): el mismo comando con --unit=orient-src
    # comprobar cada run: systemctl --user status orient-ref ; journalctl --user -u orient-ref | tail

Análisis (solo lectura de los archivos de la campaña):

    root -l -b -q 'analysis/followup_2026_10_31/orientation_corner_test.C+("<segmentos del run con fuente, separados por comas>", "<segmentos del run de referencia>", 10.0)'

### Qué debe disparar, según la orientación

| Coordenada | Orientación «+» (el índice crece desde el origen) | Orientación «−» |
|---|---|---|
| X | ch15, ch0, ch13 | ch14, ch1, ch12 |
| Y | ch16, ch17, ch18 | ch31, ch30, ch29 |

La macro decide cada coordenada cuando la diferencia entre los dos grupos supera 5σ. Además compara la capa superior con la inferior en la barra de la esquina: con una β desde arriba recibe más la capa superior.

Riesgos:
- **Fuente:** radiológico (seguir las reglas del laboratorio).
- **Equipo:** ninguno con la caja cerrada.
- **No apoyar peso** sobre la caja de hod2018.

## d) RUp = 0.1 V/s en CH00-07 desde GECO

Solo si René decide hacer el barrido de bias esa noche.

1. En CAENGECO2020, en la SY5527, slot 3, canales CH00 a CH07, **anotar el RUp actual** de cada canal (se espera 1 V/s) y sacar una foto de la pantalla.
2. Cambiar **solo RUp** a **0.1** V/s en CH00-07. No tocar V0Set, I0Set, Trip, RDWn, Pw ni otro parámetro, ni CH08-11.
3. Comprobación en el host, leída por `hvmon` (solo lectura):

        python3 -I -c "import json;s=json.load(open('/dev/shm/daqv2/hvmon_state.json'));print([c['values']['RUp_V_per_s'] for c in s['channels'] if c['ch']<8])"

   Debe dar `[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]`.
4. Después del barrido, **volver a RUp = 1 V/s**: el barrido no lo restaura.

Riesgos:
- RUp solo cambia la velocidad de la rampa. Con 0.1 V/s una subida de 0 a 65 V tarda 11 min en vez de 1, lo que importa si un canal se enciende.
- Un error de tecleo en otro parámetro (V0Set, I0Set) sí es peligroso: revisar antes de confirmar.

## e) Volver a arrancar la física

Después de b (enlace bueno) y de c si se hizo, René o el agente:

    systemctl --user is-active daq-physics            # debe dar inactive
    systemd-run --user --unit=daq-physics --collect --property=RuntimeMaxSec=1209600 --property=KillSignal=SIGTERM \
      --property=TimeoutStopSec=300 --property=StandardOutput=append:/data/daq/logs/daq-physics.log \
      --property=StandardError=append:/data/daq/logs/daq-physics.log --setenv=LD_LIBRARY_PATH=$HOME/root6/lib \
      --working-directory=$HOME ~/.local/opt/daqv2/bin/daqrun --config /data/daq/config/hodoscope_2026_10/physics_ext_B.json
    tail -f /data/daq/logs/daq-physics.log

(Es el mismo comando que usan `biasscan.py` y `biasscan_safety.sh`. La campaña `hodoscope_2026_10` sigue abierta: está `CAMPAIGN_OPEN`.)

**Qué debe verse en el log**, en el primer minuto:
- `INFO [main] daqrun arranca (pid …): campaña hodoscope_2026_10, run 11 (…_run011) …`
- `INFO [main] digitizador abierto: VX1740 serie 21, ROC 04.29 - Build 8716, AMC 00.13 - Build 8A15 …`
- `INFO [writer] segmento abierto: /data/daq/raw/hodoscope_2026_10/…_run011_seg001.root`
- **Ningún** `ERROR`.

**En los primeros 10 min:**
- **No** debe aparecer `WARN [monitor] sin triggers`.
- La página `http://127.0.0.1:8080/` (o el log de vigilancia) debe dar una tasa de ~0.26-0.28 Hz, la del trigger nuevo con las tejas START.

**Si aparece `CommError` o `Timeout`:** detener la unidad (`systemctl --user stop daq-physics`), volver a b.4 y avisar a René.

Riesgos: ninguno para el hardware.
