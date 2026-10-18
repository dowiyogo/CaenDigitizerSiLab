# daq_core: primitivas de concurrencia y utilidades. Sin hardware, sin ROOT.
add_library(daq_core
  src/time.cpp
  src/clock.cpp
  src/logger.cpp
  src/stop_token.cpp
  src/secrets.cpp
  src/file_system.cpp
  src/instance_lock.cpp
  src/unix_socket.cpp)
target_include_directories(daq_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_core PUBLIC Threads::Threads nlohmann_json::nlohmann_json PRIVATE daq_warnings)
