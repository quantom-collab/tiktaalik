cmake_minimum_required(VERSION 3.18)

# Fortran needed
find_program(CMAKE_Fortran_COMPILER $ENV{FC} gfortran PATH)
enable_language(Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_program(CMAKE_Fortran_COMPILER $ENV{FC} gfortran PATH)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(CMAKE_Fortran_FLAGS "-Ofast -march=native -fopenmp \
  -ftrampoline-impl=heap \
  -Wall -Wno-trampolines -Wno-unused-dummy-argument -Wno-conversion \
  -mcmodel=large")
else()
  set (CMAKE_Fortran_FLAGS "-O2 -ftrampoline-impl=heap")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  CACHE PATH "Fortran library directory")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
  CACHE PATH "Fortran module directory")

# Add all subdirectories
add_subdirectory(evolution)
add_subdirectory(helpers)
add_subdirectory(kernels)
add_subdirectory(model)
add_subdirectory(qcd)

set(tiktaalik_libs "tiktaalik")
add_library(${tiktaalik_libs} SHARED ${sources})

# Need to pass this to the parent scope so Python knows it's a dependency
set(tiktaalik_libs "${tiktaalik_libs}" PARENT_SCOPE)

# Install a big library containing everything
install(TARGETS ${tiktaalik_libs} DESTINATION lib)

# Install module files
install(FILES ${modules} DESTINATION include)
