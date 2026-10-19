# daq_hvmon: lógica de hvmon (sondeo, cambios, escritor, comandos, estado). Sin ROOT: el
# servidor HTTP vive en apps/hvmon, que es el único que enlaza ROOT.
add_library(daq_hvmon
  src/config.cpp
  src/state.cpp
  src/change_detector.cpp
  src/writer.cpp
  src/poller.cpp
  src/commands.cpp
  src/system_probe.cpp
  src/retry.cpp
  src/notifier.cpp
  src/mapping.cpp
  src/compressor.cpp)
target_include_directories(daq_hvmon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_hvmon PUBLIC daq_hv_common daq_temp daq_core PRIVATE daq_warnings)
