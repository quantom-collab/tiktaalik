cmake_minimum_required(VERSION 3.12)

set(extlibs "" CACHE INTERNAL "") # is this needed...?

# Add all subdirectories
add_subdirectory(evolution)
add_subdirectory(helpers)
add_subdirectory(kernels)
add_subdirectory(model)
add_subdirectory(qcd)

set(tiktaalik_libs "tiktaalik")
add_library(${tiktaalik_libs} SHARED ${sources})

# Install a big library containing everything
install(TARGETS ${tiktaalik_libs} DESTINATION lib)

# Install module files
install(FILES ${modules} DESTINATION include)
