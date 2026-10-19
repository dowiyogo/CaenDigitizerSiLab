# Genera daq_build_info.h con el commit de git y si el árbol estaba sucio, en cada build (solo
# reescribe el archivo si cambió). Uso: cmake -DSRC=<repo> -DOUT=<archivo> -P git_info.cmake
execute_process(COMMAND git -C "${SRC}" rev-parse HEAD OUTPUT_VARIABLE commit OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE r1 ERROR_QUIET)
execute_process(COMMAND git -C "${SRC}" status --porcelain --untracked-files=no -- daq OUTPUT_VARIABLE dirty
                OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE r2 ERROR_QUIET)
if(NOT r1 EQUAL 0)
  set(commit "unknown")
endif()
if(dirty STREQUAL "")
  set(isdirty "false")
else()
  set(isdirty "true")
endif()
set(content "// Generado por daq/cmake/git_info.cmake en cada build.\n#pragma once\n#define DAQ_GIT_COMMIT \"${commit}\"\n#define DAQ_GIT_DIRTY ${isdirty}\n")
if(EXISTS "${OUT}")
  file(READ "${OUT}" old)
else()
  set(old "")
endif()
if(NOT old STREQUAL content)
  file(WRITE "${OUT}" "${content}")
endif()
