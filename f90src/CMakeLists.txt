# Fortran needed
find_program(CMAKE_Fortran_COMPILER $ENV{FC} gfortran PATH)
enable_language(Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_program(CMAKE_Fortran_COMPILER $ENV{FC} gfortran PATH)

# TODO: clean up/sanitize Fortran flags

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(CMAKE_Fortran_FLAGS "-Ofast -march=native -fopenmp \
  -Wall -Wno-trampolines -Wno-unused-dummy-argument -Wno-conversion \
  -mcmodel=large")
else()
  set (CMAKE_Fortran_FLAGS "-O2")
endif()

# See if the nested functions can be implemented on the heap
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-ftrampoline-impl=heap" can_use_heap)
if(can_use_heap)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ftrampoline-impl=heap")
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
