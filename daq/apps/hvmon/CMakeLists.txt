# hvmon: único ejecutable de la fase 1 que enlaza ROOT (THttpServer, TGraph).
set(DASH_H ${CMAKE_CURRENT_BINARY_DIR}/dash_html.h)
add_custom_command(OUTPUT ${DASH_H}
  COMMAND ${CMAKE_COMMAND} -DIN=${CMAKE_CURRENT_SOURCE_DIR}/dash.html -DOUT=${DASH_H} -DNAME=kDashHtml
          -P ${PROJECT_SOURCE_DIR}/cmake/embed_text.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/dash.html ${PROJECT_SOURCE_DIR}/cmake/embed_text.cmake
  COMMENT "dash.html -> dash_html.h")
add_executable(hvmon main.cpp http_page.cpp ${DASH_H})
target_include_directories(hvmon PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/io/include)
target_link_libraries(hvmon PRIVATE daq_hvmon daq_hv_caen daq_temp daq_core ROOT::Core ROOT::Hist ROOT::RHTTP daq_warnings)
