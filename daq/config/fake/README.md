# Guiones de `FakeHVBackend`

Se usan con `hvctl --fake <archivo>` (o `--fake default`) y en los tests. Los tiempos `t_s`
son segundos desde que se crea el backend; con `hvctl` corren en tiempo real, en los tests con
un reloj manual (tiempo acelerado).

- Sin `"base"`: el archivo describe todo (`slots`, `events`, `seed`, `crate_slots`); ver
  `FakeHVBackend::defaultScenario()` en `daq/hv/src/fake_hv_backend.cpp`.
- `"base": "default"`: usa el A7030DP del slot 3 con los valores leídos el 2026-10-18 con HV
  apagada, y el archivo solo aporta `events` y `seed`.

Valores en unidades de la placa (V, µA, V/s, s), como los muestra `hv_readall`.

Eventos:

| tipo | campos | efecto |
|---|---|---|
| `set` | `ch`, `param`, `value` | cambio externo de un parámetro (como si se hiciera en GECO); `Pw = 1` borra el trip |
| `current` | `ch`, `extra_uA`, `duration_s` | excursión de corriente (p. ej. luz); escala con V/V0Set |
| `connection_loss` | `duration_s` | todas las llamadas fallan; al terminar hay que reconectar |

Archivos:

- `trip_ch2.json`: enciende CH02 en t = 1 s (rampa de 65.9 s a 1 V/s); en t = 100 s una
  excursión de +20 µA supera I0Set = 15 µA y, tras Trip = 5 s, el canal se apaga y queda
  apagado con `internal_trip`.
- `connection_loss.json`: el mainframe deja de responder entre t = 10 s y t = 30 s.
