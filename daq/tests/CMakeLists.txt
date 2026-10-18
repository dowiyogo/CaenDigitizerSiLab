# Tests de la fase 0 (GoogleTest). Solo usan fakes: ninguno abre el digitizador ni conecta al SY5527.
add_executable(daq_tests
  test_bounded_queue.cpp
  test_jthread.cpp
  test_snapshot.cpp
  test_clock_time.cpp
  test_logger.cpp
  test_units_status.cpp
  test_credentials.cpp
  test_fake_hv.cpp)
target_link_libraries(daq_tests PRIVATE daq_hv_common daq_core GTest::gtest GTest::gtest_main daq_warnings)
# TSan de GCC 13 falla con "unexpected memory mapping" con la entropía ASLR de este kernel; se
# corre el binario con ASLR deshabilitado solo para ese proceso (setarch -R), sin tocar sysctl.
if(DAQ_SANITIZE MATCHES "thread")
  set_property(TARGET daq_tests PROPERTY CROSSCOMPILING_EMULATOR setarch ${CMAKE_HOST_SYSTEM_PROCESSOR} -R)
endif()
gtest_discover_tests(daq_tests DISCOVERY_TIMEOUT 60 PROPERTIES TIMEOUT 300)
