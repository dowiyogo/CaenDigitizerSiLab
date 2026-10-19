# daqrun: un run de adquisición. Fase 3: solo digitizador fake o replay (ningún enlace a CAENDigitizer).
add_executable(daqrun main.cpp http_page.cpp)
add_dependencies(daqrun daq_build_info)
target_include_directories(daqrun PRIVATE ${DAQ_GENERATED_DIR})
target_link_libraries(daqrun PRIVATE daq_io daq_dig_replay daq_dig daq_daqrun daq_hvmon daq_core
                      ROOT::Core ROOT::Hist ROOT::RHTTP daq_warnings)
