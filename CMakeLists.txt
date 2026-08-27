cmake_minimum_required(VERSION 3.18)

project(tiktaalik C Fortran)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.27")
	find_package(Python REQUIRED COMPONENTS Interpreter Development)
	include_directories(${Python_INCLUDE_DIRS})
else()
  find_package(PythonExtensions REQUIRED)
endif()
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)

add_subdirectory(f90src)
add_subdirectory(tiktaalik)
