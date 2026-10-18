# DAQ v2: adquisición del VX1740 y control lento del SY5527.
# Diseño: docs/daq_v2_design.md. Construir desde la raíz del repositorio:
#   cmake -S daq -B daq/build && cmake --build daq/build -j4 && ctest --test-dir daq/build
cmake_minimum_required(VERSION 3.20)
project(daqv2 VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
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

# ROOT: se busca para las fases siguientes (daqrun, hvmon). Ningún target de la fase 0 lo usa.
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}" "$ENV{HOME}/root6")
find_package(ROOT 6.32 REQUIRED COMPONENTS Core RIO Tree RHTTP)

# CAEN HV Wrapper 6.6: ruta explícita. La biblioteca no tiene SONAME, así que DT_NEEDED
# queda con esta ruta. UNIX es obligatorio para el header en Linux.
add_library(caen::hvwrapper SHARED IMPORTED GLOBAL)
set_target_properties(caen::hvwrapper PROPERTIES
  IMPORTED_LOCATION "/usr/lib64/libcaenhvwrapper.so"
  INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
  INTERFACE_COMPILE_DEFINITIONS "UNIX")

add_subdirectory(core)
add_subdirectory(hv)
add_subdirectory(apps/hvctl)
add_subdirectory(tools/legacy)

if(DAQ_BUILD_TESTS)
  enable_testing()
  find_package(GTest 1.14 REQUIRED)
  include(GoogleTest)
  add_subdirectory(tests)
endif()
