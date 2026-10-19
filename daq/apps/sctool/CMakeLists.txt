# sctool: utilidades de archivos de control lento (comprimir y verificar, leer planos o .gz).
add_executable(sctool main.cpp)
target_link_libraries(sctool PRIVATE daq_core daq_warnings)
