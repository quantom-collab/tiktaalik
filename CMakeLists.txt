cmake_minimum_required(VERSION 3.18)

project(tiktaalik C Fortran)

find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)

add_subdirectory(f90src)
add_subdirectory(tiktaalik)
