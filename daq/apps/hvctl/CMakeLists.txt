add_executable(hvctl main.cpp)
target_link_libraries(hvctl PRIVATE daq_hvmon daq_hv_caen daq_hv_common daq_core daq_warnings)
