# daqoffline: monitor, bits de datos, paneles y reloj del TTT de daqrun sobre runs ya escritos (o
# sobre el FakeDigitizer con reloj manual). Solo lectura de los datos; salida en /data/daq/sim/offline.
add_executable(daqoffline main.cpp)
target_link_libraries(daqoffline PRIVATE daq_io daq_dig daq_daqrun daq_hvmon daq_core ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
                      daq_warnings)
