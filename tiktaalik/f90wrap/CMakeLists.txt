# Function to generate an f2py module and add it to a list of modules to install
function(gen_f2pymod f2py_module_name module_list)
  set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/${f2py_module_name}.f90")
  set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
  # Custom command to execture f2py, since cmake's built-in methods don't seem to work.
  add_custom_command(
    OUTPUT ${generated_module_file}
    COMMAND ${F2PY_EXECUTABLE} -m ${f2py_module_name} -c ${fortran_src_file} -I${CMAKE_BINARY_DIR}/include -L${CMAKE_BINARY_DIR}/lib -ltiktaalik
    DEPENDS ${fortran_src_file}
    VERBATIM
  )
  add_custom_target(${f2py_module_name}
    ALL DEPENDS ${generated_module_file} ${tiktaalik_libs}
  )
  # Add the generated module to a list
  list(APPEND module_list "${generated_module_file}")
  set(module_list "${module_list}" PARENT_SCOPE)
endfunction()

set(module_list "")

# List of f2py modules
gen_f2pymod(matrices "${module_list}")
gen_f2pymod(model "${module_list}")
gen_f2pymod(qcd "${module_list}")

install(FILES ${module_list} DESTINATION tiktaalik/f90wrap)
