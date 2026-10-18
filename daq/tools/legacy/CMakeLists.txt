# Herramientas HV heredadas del chat de HV (2026-10-18), compiladas tal cual desde sus fuentes,
# sin las advertencias estrictas del resto de daq/. hv_readall y test_sy5527 solo leen.
add_executable(hv_readall hv_readall.c)
target_link_libraries(hv_readall PRIVATE caen::hvwrapper)
target_compile_options(hv_readall PRIVATE -Wall)

add_executable(test_sy5527 test_sy5527.c)
target_link_libraries(test_sy5527 PRIVATE caen::hvwrapper)
target_compile_options(test_sy5527 PRIVATE -Wall)
