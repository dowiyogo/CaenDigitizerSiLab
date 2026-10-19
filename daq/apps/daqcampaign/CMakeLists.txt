# daqcampaign: abre, cierra y muestra campañas.
add_executable(daqcampaign main.cpp)
target_link_libraries(daqcampaign PRIVATE daq_daqrun daq_core daq_warnings)
