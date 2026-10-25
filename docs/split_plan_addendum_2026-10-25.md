# Adenda a B.5 del plan de división (XPLUS-3, 2026-10-25)

Complementa las secciones B.4 y B.5 de `/home/xplus/xplus_split_plan_2026-10-25_report.md` (XPLUS-2), en los dos puntos que encontró la revisión de SPLIT-1. Lo que esta adenda no menciona sigue igual.

## 1. `hodo_paths.h`: usuarios en el análisis y cambio del Makefile

**Corrección de un número.** En la revisión de SPLIT-1 dije "28 archivos que lo incluyen". Ese 28 era `git grep -l hodo_paths` sobre `analysis/` y `daq/tools/`, es decir, archivos que *mencionan* el nombre: incluía el propio `hodo_paths.h`, el `Makefile`, `validate_replay_v8.C` (DAQ) y `split_paths.py`. Los números correctos (medido, `git grep -n -E '#include\s*[<"][^>"]*hodo_paths\.h[>"]'` en `eafd3b7`):
- **24 archivos del análisis lo incluyen**: 19 en `analysis/src`, 5 en `analysis/macros`;
- **`analysis/Makefile`** lo usa como prerrequisito;
- **3 archivos del DAQ** lo incluyen (`daq/dig/src/replay_digitizer.cpp`, `daq/tests/test_root_io.cpp`, `daq/tools/validation/validate_replay_v8.C`). Estos ya están en la lista de ajuste de daqv2 de B.5.

`hodo_paths.h` se muda a daqv2 como `core/include/daq/hodo_paths.h` (D3). En daq_analysis se lee desde el submódulo, `external/daqv2/core/include/daq/hodo_paths.h`. Rutas después de la división (`daq_analysis/src`, `daq_analysis/macros`); entre paréntesis, la ruta de hoy:

| archivo en daq_analysis | línea | hoy | después |
|---|---|---|---|
| `daq_analysis/macros/hod2018_amplitude.C` (`analysis/macros/hod2018_amplitude.C`) | 38 | `#include "../src/hodo_paths.h"` | `R__ADD_INCLUDE_PATH(external/daqv2/core/include) + #include "daq/hodo_paths.h"` |
| `daq_analysis/macros/hv_diag3_gain_vs_time.C` (`analysis/macros/hv_diag3_gain_vs_time.C`) | 28 | `#include "../src/hodo_paths.h"` | `R__ADD_INCLUDE_PATH(external/daqv2/core/include) + #include "daq/hodo_paths.h"` |
| `daq_analysis/macros/hv_diag3_pretrip.C` (`analysis/macros/hv_diag3_pretrip.C`) | 23 | `#include "../src/hodo_paths.h"` | `R__ADD_INCLUDE_PATH(external/daqv2/core/include) + #include "daq/hodo_paths.h"` |
| `daq_analysis/macros/hv_diag4.C` (`analysis/macros/hv_diag4.C`) | 26 | `#include "../src/hodo_paths.h"` | `R__ADD_INCLUDE_PATH(external/daqv2/core/include) + #include "daq/hodo_paths.h"` |
| `daq_analysis/macros/hv_diag5_episodes.C` (`analysis/macros/hv_diag5_episodes.C`) | 22 | `#include "../src/hodo_paths.h"` | `R__ADD_INCLUDE_PATH(external/daqv2/core/include) + #include "daq/hodo_paths.h"` |
| `daq_analysis/src/analyze_highgain.cxx` (`analysis/src/analyze_highgain.cxx`) | 44 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/analyze_hod2018.cpp` (`analysis/src/analyze_hod2018.cpp`) | 16 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/analyze_independent.cxx` (`analysis/src/analyze_independent.cxx`) | 45 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/channel_health_analysis.cpp` (`analysis/src/channel_health_analysis.cpp`) | 31 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/channel_health_refit.cpp` (`analysis/src/channel_health_refit.cpp`) | 22 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/charge_analysis.cxx` (`analysis/src/charge_analysis.cxx`) | 47 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/compare_hv_runs.cpp` (`analysis/src/compare_hv_runs.cpp`) | 26 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/final_hodoscope_analysis.cpp` (`analysis/src/final_hodoscope_analysis.cpp`) | 25 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/infer_mapping.cxx` (`analysis/src/infer_mapping.cxx`) | 54 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/multirun_analysis.cxx` (`analysis/src/multirun_analysis.cxx`) | 40 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/poisson_analysis.cxx` (`analysis/src/poisson_analysis.cxx`) | 38 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/reconstruct_corrected.cxx` (`analysis/src/reconstruct_corrected.cxx`) | 43 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/reconstruct_cosmic_tracks.cxx` (`analysis/src/reconstruct_cosmic_tracks.cxx`) | 45 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/reconstruct_v8.cxx` (`analysis/src/reconstruct_v8.cxx`) | 63 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/reconstruct_v8b.cxx` (`analysis/src/reconstruct_v8b.cxx`) | 52 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/scan_active_channels.cpp` (`analysis/src/scan_active_channels.cpp`) | 15 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/scan_active_channels_windowed.cpp` (`analysis/src/scan_active_channels_windowed.cpp`) | 19 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/spe_gain_analysis.cpp` (`analysis/src/spe_gain_analysis.cpp`) | 24 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |
| `daq_analysis/src/tagprobe_analysis.cxx` (`analysis/src/tagprobe_analysis.cxx`) | 43 | `#include "hodo_paths.h"` | `#include "daq/hodo_paths.h"` |

En las 5 macros se usa la misma forma que ya tiene `analysis/macros/plot_hv_imon_day.C` (`R__ADD_INCLUDE_PATH(daq/io/include)`): `R__ADD_INCLUDE_PATH(external/daqv2/core/include)` antes del `#include`. Como hoy, las macros se corren desde la raíz del repositorio, que pasa a ser la raíz del monorepo.

**`daq_analysis/Makefile`** (hoy `analysis/Makefile`). Se corre con `make -C daq_analysis`, así que el submódulo queda en `../external/daqv2`:

```make
# antes
CXXFLAGS := -O2 -Wall -std=c++17 -Isrc $(shell root-config --cflags)
...
$(addprefix $(BIN_DIR)/, $(PROGRAMS_PATHS)): $(SRC_DIR)/hodo_paths.h

# después
DAQV2 ?= ../external/daqv2
CXXFLAGS := -O2 -Wall -std=c++17 -Isrc -I$(DAQV2)/core/include $(shell root-config --cflags)
...
$(addprefix $(BIN_DIR)/, $(PROGRAMS_PATHS)): $(DAQV2)/core/include/daq/hodo_paths.h
```

`PROGRAMS_PATHS` ya lista exactamente los 19 programas de `src/` que incluyen el encabezado (medido: los mismos 19 nombres). No cambia nada más: `hodo_mapping.h` y `hodo_consts.h` siguen en `src/`. Si el submódulo falta, `make` debe fallar con un mensaje claro; propongo `$(if $(wildcard $(DAQV2)/core/include/daq/hodo_paths.h),,$(error falta el submódulo external/daqv2: git submodule update --init))`.

**Punto abierto, sin decidir aquí.** `hodo_paths.h` tiene como valor por defecto `/home/xplus/digitizer/data/raw`, que hoy es un enlace a `/data/daq/raw`, y se puede reemplazar con `HODO_RAW_DIR`. Es una ruta de datos, no relativa al archivo. D5 la deja en daqv2 tal cual. Para el análisis, propongo que el README del módulo diga `HODO_RAW_DIR=$MUOGRAPHY_DATA/raw`, sin tocar el encabezado. El significado del encabezado no cambia en la división.

## 2. Secciones de análisis de CLAUDE.md: fuera de daqv2 desde el principio (opción a)

Esto reemplaza el ítem de la etapa 1.5 de SPLIT-1 ("move them to a file `docs/analysis_sections_from_CLAUDE.md` … then delete that file from daqv2 in the same commit"). **Ese archivo nunca entra a daqv2.** En el commit de ajuste de daqv2:
- `CLAUDE.md` queda sin las líneas de abajo;
- el texto quitado se guarda, fuera de todo repositorio, en `$HOME/split_2026-10/analysis_sections_from_CLAUDE.md`;
- en la etapa 2 ese archivo es la fuente de `daq_analysis/README.md` (y de las reglas de `daq_analysis/AGENTS.md`), y queda en `$HOME/split_2026-10/` hasta que René cierre la migración.

Líneas del `CLAUDE.md` actual (`eafd3b7`, 259 líneas):

| líneas | contenido | en daqv2 | a `analysis_sections_from_CLAUDE.md` |
|---|---|---|---|
| 3-6 | "Read `AGENTS.md` as well: it holds the detailed record of the old DAQ, runs v1-v10 and the analysis…" | se reescribe sin AGENTS.md, que va a daq_analysis: "Guidance for Claude Code sessions in the daqv2 repository." | sí, completas |
| 23-25 | Inputs: `docs/handoff_analysis_2026-10-18.md`, `docs/notes_for_analysis.md`, "do not fix analysis code from DAQ tasks" | se quita (esos docs van a daq_analysis). Queda una línea: "Analysis code lives in the monorepo module `daq_analysis/`; do not fix it from DAQ tasks." | sí |
| 26-28 | calibración de octubre / `pe_calibration_v8.txt` de mayo, estado "no calibration" | **se queda** (es estado del DAQ) | también (copia) |
| 40 | `analysis outputs outputs/ (not in git)` | se quita | sí, como "results → `data/results/daq_analysis/<fecha>_<tag>/`" |
| 112-115 | bloque "Analysis (Makefile…)": `make -C analysis`, `analysis/scripts/run_*.sh` | se quita | sí, con las rutas nuevas (`make -C daq_analysis`) |
| 215 | ejemplo `root … analysis/macros/plot_hv_imon_day.C(...)` | se quita; el lector de slow control sigue documentado por `libdaq_screader` (líneas 223-227) | sí |
| 223-227 | "Shared with the analysis (header-only, no ROOT): pulse_finder.h, trigger_interval_stats.h, ttt.h, gain_fit.h" | **se queda** y se añade `hodo_paths.h`, `shared/xy_estimator.py` y `shared/hitdef_extract.C` (D3) | no |
| 237-240 | convención de tiempos; "runs v1-v10 and `runs/good_intervals_2026.csv` are UTC-4" | se queda la convención; se quita la mención de `runs/good_intervals_2026.csv`, que va a daq_analysis | sí, la frase de UTC-4 de v1-v10 |
| 247-254, 257-258 | "Known differences between AGENTS.md and the system": ruta de crudos, `hodo_mapping.h` kCfg_2026_05, `timeout 600`, regla de la última hora, pickup en 3525 ns | se quita la sección | sí, completas |
| 255-256 | `hv_readall` toma el password como argumento | se reescribe en la sección de rutas o reglas: "`attic/legacy/hv_readall` and `test_sy5527`: superseded by `hvctl`; do not use." | no |

Las reglas de host (1-12) **se quedan todas en daqv2**, porque gobiernan el host y el hardware. Hay que ajustar dos menciones que dejan de existir en daqv2:
- regla 4, "`outputs/_superseded/` (results)": queda como regla de daq_analysis y en daqv2 se dice solo `attic/`;
- regla 7, el ejemplo "`analysis/scripts/hv_diag6_shape.py`": se quita el ejemplo y queda la regla.

Ninguna regla cambia de contenido.

## 3. Fixtures para la prueba de humo de la etapa 2.7

`docs/split/fixtures/xy_estimator_run8_T5.json` y `docs/split/fixtures/xy_estimator_may8good_T5.json`, de `daq/tools/proceeding/export_xy_fixture.py` (XPLUS-3). Ver el informe XPLUS-3, sección 3.
