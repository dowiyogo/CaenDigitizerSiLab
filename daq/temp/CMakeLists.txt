# daq_temp: interfaz de temperatura y backend fake (fase 1). Sin hardware.
add_library(daq_temp src/temp_backend.cpp)
target_include_directories(daq_temp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_temp PUBLIC daq_core PRIVATE daq_warnings)
