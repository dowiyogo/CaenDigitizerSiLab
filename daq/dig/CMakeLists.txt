# daq_dig: interfaz del digitizador, FakeDigitizer sintético y bucle de lectura. Sin ROOT.
add_library(daq_dig src/digitizer.cpp src/fake_digitizer.cpp src/readout.cpp)
target_include_directories(daq_dig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_dig PUBLIC daq_core PRIVATE daq_warnings)

# daq_dig_replay: replay de los archivos v8 (lee ROOT en modo READ; hodo_paths.h de analysis/).
add_library(daq_dig_replay src/replay_digitizer.cpp)
target_include_directories(daq_dig_replay PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                           PRIVATE ${PROJECT_SOURCE_DIR}/../analysis/src)
target_link_libraries(daq_dig_replay PUBLIC daq_dig ROOT::Core ROOT::RIO ROOT::Tree PRIVATE daq_warnings)
