# hvmon: único ejecutable de la fase 1 que enlaza ROOT (THttpServer, TGraph).
add_executable(hvmon main.cpp http_page.cpp)
target_link_libraries(hvmon PRIVATE daq_hvmon daq_hv_caen daq_temp daq_core ROOT::Core ROOT::Hist ROOT::RHTTP daq_warnings)
