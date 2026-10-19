# digprobe: información, registros y prueba de memoria llena del VX1740 (fase 4, con autorización).
add_executable(digprobe main.cpp)
target_link_libraries(digprobe PRIVATE daq_dig_caen daq_dig daq_core daq_warnings)
