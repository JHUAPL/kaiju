# Glob all Fortran source files
file(GLOB volt_mpi_srcs *.F90)

# Define the static library target
add_library(voltmpilib ${volt_mpi_srcs})

# Link the library with its dependencies
target_link_libraries(voltmpilib voltlib rcmlib remixlib chimplib gammpilib gamlib basempilib baselib)

# Add a custom target to run ranlib on the static library
add_custom_target(run_ranlib_voltmpilib ALL
    COMMAND ${CMAKE_RANLIB} $<TARGET_FILE:voltmpilib>
    DEPENDS voltmpilib
    COMMENT "Running ranlib on libvoltmpilib.a"
)
