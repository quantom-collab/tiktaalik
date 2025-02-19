cmake_minimum_required(VERSION 3.12)

set(localfiles
  "gridspace"
  "integration"
  "pixelation"
  "specfun"
  )

set(localsources ${localfiles})
list(TRANSFORM localsources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
list(TRANSFORM localsources APPEND ".f90")
list(APPEND sources ${localsources})
set(sources "${sources}" PARENT_SCOPE)

set(localmodules ${localfiles})
list(TRANSFORM localmodules PREPEND "${CMAKE_Fortran_MODULE_DIRECTORY}/")
list(TRANSFORM localmodules APPEND ".mod")
list(APPEND modules ${localmodules})
set(modules "${modules}" PARENT_SCOPE)
