# daq_hv: biblioteca HV de solo lectura (fase 0). Sin ROOT.
# daq_hv_common: tipos, Status, credenciales y backend fake (sin la Wrapper).
add_library(daq_hv_common
  src/hv_types.cpp
  src/status_bits.cpp
  src/credentials.cpp
  src/fake_hv_backend.cpp
  src/reading_json.cpp)
target_include_directories(daq_hv_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(daq_hv_common PUBLIC daq_core nlohmann_json::nlohmann_json PRIVATE daq_warnings)

# daq_hv_caen: backend real sobre CAEN HV Wrapper 6.6.
add_library(daq_hv_caen src/caen_hv_backend.cpp)
target_link_libraries(daq_hv_caen PUBLIC daq_hv_common PRIVATE caen::hvwrapper daq_warnings)
