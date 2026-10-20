# daq_daqrun: configuración, campañas, cliente HV y flags, monitor y notas de daqrun. Sin ROOT.
add_library(daq_daqrun src/config.cpp src/campaign.cpp src/hv_client.cpp src/monitor.cpp src/notes.cpp src/dash.cpp src/data_flags.cpp src/panels.cpp)
target_include_directories(daq_daqrun PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_daqrun PUBLIC daq_dig daq_hvmon daq_hv_common daq_core PRIVATE daq_warnings)
