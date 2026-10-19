# Tests (GoogleTest). Solo usan fakes: ninguno abre el digitizador ni conecta al SY5527, y ninguno
# lee el secrets.env real (DAQV2_SECRETS apunta a temporales).
add_executable(daq_tests
  test_bounded_queue.cpp
  test_jthread.cpp
  test_snapshot.cpp
  test_clock_time.cpp
  test_logger.cpp
  test_units_status.cpp
  test_credentials.cpp
  test_fake_hv.cpp
  test_secrets.cpp
  test_infra.cpp
  test_hvmon_poller.cpp
  test_hvmon_writer.cpp
  test_hvmon_commands.cpp
  test_hvmon_misc.cpp
  test_hvmon_retry.cpp
  test_notifier.cpp
  test_gzip.cpp
  test_analysis_shared.cpp
  test_dig.cpp
  test_daqrun.cpp
  test_x740.cpp)
target_compile_definitions(daq_tests PRIVATE DAQ_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(daq_tests PRIVATE daq_daqrun daq_dig daq_hvmon daq_temp daq_hv_common daq_core GTest::gtest
                      GTest::gtest_main daq_warnings)

# Tests con ROOT (escritor, replay, lector de control lento): binario aparte para que daq_tests no
# dependa de ROOT. Solo archivos temporales; el replay usa un archivo pequeño con el formato de v8.
add_executable(daq_root_tests test_root_io.cpp)
target_compile_definitions(daq_root_tests PRIVATE DAQ_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_include_directories(daq_root_tests PRIVATE ${PROJECT_SOURCE_DIR}/../analysis/src)
target_link_libraries(daq_root_tests PRIVATE daq_io daq_dig_replay daq_screader daq_daqrun daq_dig daq_hvmon daq_core
                      ROOT::Core ROOT::RIO ROOT::Tree GTest::gtest GTest::gtest_main daq_warnings)
# TSan de GCC 13 falla con "unexpected memory mapping" con la entropía ASLR de este kernel; se
# corre el binario con ASLR deshabilitado solo para ese proceso (setarch -R), sin tocar sysctl.
if(DAQ_SANITIZE MATCHES "thread")
  set_property(TARGET daq_tests PROPERTY CROSSCOMPILING_EMULATOR setarch ${CMAKE_HOST_SYSTEM_PROCESSOR} -R)
  set_property(TARGET daq_root_tests PROPERTY CROSSCOMPILING_EMULATOR setarch ${CMAKE_HOST_SYSTEM_PROCESSOR} -R)
endif()
gtest_discover_tests(daq_tests DISCOVERY_TIMEOUT 60 PROPERTIES TIMEOUT 300 ENVIRONMENT "DAQV2_SECRETS=/nonexistent/daq_tests_secrets.env")
gtest_discover_tests(daq_root_tests DISCOVERY_TIMEOUT 120 PROPERTIES TIMEOUT 600
                     ENVIRONMENT "DAQV2_SECRETS=/nonexistent/daq_tests_secrets.env")
