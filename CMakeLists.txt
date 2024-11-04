cmake_minimum_required(VERSION 3.12)

# Fortran needed
find_program(CMAKE_Fortran_COMPILER $ENV{FC} gfortran PATH)

# Python needed
find_program(PYTHON python3 python PATH)
message(STATUS "Found Python: ${PYTHON}")
add_custom_target(python_interface ALL)

project(tiktaalik)
set(tiktaalik_VERSION_MAJOR 0)
set(tiktaalik_VERSION_MINOR 1)

#IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#  SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/binaries/" CACHE PATH "Potato" FORCE)
#ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

enable_language(Fortran)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(CMAKE_Fortran_FLAGS "-Ofast -march=native -fopenmp \
  -Wall -Wno-trampolines -Wno-unused-dummy-argument -Wno-conversion \
  -mcmodel=large")
else()
  set (CMAKE_Fortran_FLAGS "-O2")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  CACHE PATH "Fortran library directory")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
  CACHE PATH "Fortran module directory")

add_subdirectory(f90src)
add_subdirectory(python)

################################################################################

# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
