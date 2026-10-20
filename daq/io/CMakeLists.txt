# daq_io: escritor ROOT de daqrun (segmentos, índice, SHA256).
add_library(daq_io src/run_writer.cpp src/panels_root.cpp)
target_include_directories(daq_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_io PUBLIC daq_daqrun ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::RHTTP PRIVATE daq_warnings)

# libdaq_screader.so: lector de control lento a TTree, para macros (R__LOAD_LIBRARY).
add_library(daq_screader SHARED src/sc_reader.cpp)
target_include_directories(daq_screader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_screader PUBLIC ROOT::Core ROOT::Tree PRIVATE daq_core daq_warnings)
