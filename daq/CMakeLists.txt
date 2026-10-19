# DAQ v2: adquisición del VX1740 y control lento del SY5527.
# Diseño: docs/daq_v2_design.md. Construir desde la raíz del repositorio:
#   cmake -S daq -B daq/build && cmake --build daq/build -j4 && ctest --test-dir daq/build
cmake_minimum_required(VERSION 3.20)
project(daqv2 VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# libdaq_screader.so (macros de análisis) enlaza las bibliotecas estáticas: todo con -fPIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Tipo de build" FORCE)
endif()

option(DAQ_BUILD_TESTS "Compilar los tests (GoogleTest)" ON)
set(DAQ_SANITIZE "OFF" CACHE STRING "Sanitizer: OFF, address, undefined, address,undefined, thread")
set_property(CACHE DAQ_SANITIZE PROPERTY STRINGS OFF address undefined "address,undefined" thread)

if(NOT DAQ_SANITIZE STREQUAL "OFF")
  if(DAQ_SANITIZE MATCHES "thread" AND DAQ_SANITIZE MATCHES "address")
    message(FATAL_ERROR "DAQ_SANITIZE: thread no se combina con address")
  endif()
  add_compile_options(-fsanitize=${DAQ_SANITIZE} -fno-omit-frame-pointer -fno-sanitize-recover=all)
  add_link_options(-fsanitize=${DAQ_SANITIZE})
endif()

# Advertencias estrictas solo para nuestro código (los fuentes legacy en C no las usan).
add_library(daq_warnings INTERFACE)
target_compile_options(daq_warnings INTERFACE
  $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow -Wconversion>)

find_package(Threads REQUIRED)
find_package(nlohmann_json 3.11 REQUIRED)

# ROOT: lo enlazan apps/hvmon, apps/daqrun, daq_io (escritor y lector de control lento) y
# daq_dig_replay; daq_core, daq_hv*, daq_hvmon, daq_dig y daq_daqrun no dependen de ROOT.
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}" "$ENV{HOME}/root6")
find_package(ROOT 6.32 REQUIRED COMPONENTS Core RIO Tree Hist RHTTP)

# CAEN HV Wrapper 6.6: ruta explícita. La biblioteca no tiene SONAME, así que DT_NEEDED
# queda con esta ruta. UNIX es obligatorio para el header en Linux.
add_library(caen::hvwrapper SHARED IMPORTED GLOBAL)
set_target_properties(caen::hvwrapper PROPERTIES
  IMPORTED_LOCATION "/usr/lib64/libcaenhvwrapper.so"
  INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
  INTERFACE_COMPILE_DEFINITIONS "UNIX")

# Commit de git del build (cabecera de cada run de daqrun).
set(DAQ_GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
add_custom_target(daq_build_info ALL
  COMMAND ${CMAKE_COMMAND} -DSRC=${PROJECT_SOURCE_DIR}/.. -DOUT=${DAQ_GENERATED_DIR}/daq_build_info.h
          -P ${PROJECT_SOURCE_DIR}/cmake/git_info.cmake
  BYPRODUCTS ${DAQ_GENERATED_DIR}/daq_build_info.h
  COMMENT "commit de git del build")

# CAENDigitizer 2.17.3 (VX1740, fase 4): copias versionadas de /usr/lib por ruta completa (sin
# SONAME, así DT_NEEDED queda con la ruta). CAENDigitizer resuelve CAENComm por nombre; las copias de
# /usr/lib y /usr/local/lib son idénticas (DAQ2-1).
add_library(caen::digitizer SHARED IMPORTED GLOBAL)
set_target_properties(caen::digitizer PROPERTIES IMPORTED_LOCATION "/usr/lib/libCAENDigitizer.so.v2.17.3"
                      INTERFACE_INCLUDE_DIRECTORIES "/usr/include" INTERFACE_COMPILE_DEFINITIONS "LINUX;UNIX")
add_library(caen::comm SHARED IMPORTED GLOBAL)
set_target_properties(caen::comm PROPERTIES IMPORTED_LOCATION "/usr/lib/libCAENComm.so.v1.7.0")
add_library(caen::vme SHARED IMPORTED GLOBAL)
set_target_properties(caen::vme PROPERTIES IMPORTED_LOCATION "/usr/lib/libCAENVME.so.v4.0.2")

add_subdirectory(core)
add_subdirectory(hv)
add_subdirectory(temp)
add_subdirectory(hvmon)
add_subdirectory(dig)
add_subdirectory(daqrun)
add_subdirectory(io)
add_subdirectory(apps/hvctl)
add_subdirectory(apps/hvmon)
add_subdirectory(apps/sctool)
add_subdirectory(apps/daqrun)
add_subdirectory(apps/daqcampaign)
add_subdirectory(apps/digprobe)
add_subdirectory(tools/legacy)

if(DAQ_BUILD_TESTS)
  enable_testing()
  find_package(GTest 1.14 REQUIRED)
  include(GoogleTest)
  add_subdirectory(tests)
endif()
