# daqrun: un run de adquisición. Fase 3: solo digitizador fake o replay (ningún enlace a CAENDigitizer).
set(DASH_H ${CMAKE_CURRENT_BINARY_DIR}/dash_html.h)
add_custom_command(OUTPUT ${DASH_H}
  COMMAND ${CMAKE_COMMAND} -DIN=${CMAKE_CURRENT_SOURCE_DIR}/dash.html -DOUT=${DASH_H} -DNAME=kDashHtml
          -P ${PROJECT_SOURCE_DIR}/cmake/embed_text.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/dash.html ${PROJECT_SOURCE_DIR}/cmake/embed_text.cmake
  COMMENT "dash.html -> dash_html.h")
add_executable(daqrun main.cpp http_page.cpp ${DASH_H})
add_dependencies(daqrun daq_build_info)
target_include_directories(daqrun PRIVATE ${DAQ_GENERATED_DIR} ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(daqrun PRIVATE daq_io daq_dig_replay daq_dig_caen daq_dig daq_daqrun daq_hvmon daq_hv_common daq_core
                      ROOT::Core ROOT::Hist ROOT::RHTTP daq_warnings)
