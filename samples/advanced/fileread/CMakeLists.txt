cmake_minimum_required(VERSION 3.15)
project(AMGExamples Fortran)

# Installation directories (passed as CMake variables)
set(AMG_INSTALL_DIR "" CACHE PATH "Path to AMG installation")
set(PSBLAS_INSTALL_DIR "" CACHE PATH "Path to PSBLAS installation")

# Check if installation directories are set
if(NOT AMG_INSTALL_DIR)
    message(FATAL_ERROR "AMG_INSTALL_DIR must be set.  Use -DAMG_INSTALL_DIR=/path/to/amg when running CMake.")
endif()

if(NOT PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "PSBLAS_INSTALL_DIR must be set.  Use -DPSBLAS_INSTALL_DIR=/path/to/psblas when running CMake.")
endif()

# Include directories
include_directories(
    "."
    "${AMG_INSTALL_DIR}/include"
    "${AMG_INSTALL_DIR}/modules"
    "${PSBLAS_INSTALL_DIR}/include" # Assuming PSBLAS include directory is here
)

# Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(FMFLAG "-I") # Adjust if needed

# Library directories
link_directories(
    "${AMG_INSTALL_DIR}/lib"
    "${PSBLAS_INSTALL_DIR}/lib"
)

# Libraries
set(AMG_LIBS psb_linsolve amg_prec psb_prec)
set(PSBLAS_LIBS psb_util psb_linsolve psb_prec psb_base) # Example, adjust as needed

# Executable directory
set(EXEDIR "${CMAKE_CURRENT_BINARY_DIR}/runs")
file(MAKE_DIRECTORY "${EXEDIR}")

# Source files
set(DFSOBJS amg_df_sample.f90 data_input.f90)
set(SFSOBJS amg_sf_sample.f90 data_input.f90)
set(CFSOBJS amg_cf_sample.f90 data_input.f90)
set(ZFSOBJS amg_zf_sample.f90 data_input.f90)

# Function to create executable
macro(create_amg_executable target sources)
    add_executable(${target} ${sources})
    target_link_libraries(${target} PUBLIC
        ${AMG_LIBS}
        ${PSBLAS_LIBS}
        #${LDLIBS} # Assuming this variable is defined elsewhere if needed
    )

    # Move executable to EXEDIR (post-build)
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E move $<TARGET_FILE:${target}> ${EXEDIR}
    )
endmacro()

# Create executables
create_amg_executable(amg_df_sample ${DFSOBJS})
create_amg_executable(amg_sf_sample ${SFSOBJS})
create_amg_executable(amg_cf_sample ${CFSOBJS})
create_amg_executable(amg_zf_sample ${ZFSOBJS})

# Create "runs" directory (if it doesn't exist)
add_custom_target(create_runs_dir
    COMMAND ${CMAKE_COMMAND} -E make_directory "${EXEDIR}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating runs directory"
)
add_dependencies(amg_df_sample create_runs_dir)
add_dependencies(amg_sf_sample create_runs_dir)
add_dependencies(amg_cf_sample create_runs_dir)
add_dependencies(amg_zf_sample create_runs_dir)



# lib target
add_custom_target(lib
    COMMAND ${CMAKE_COMMAND} -E chdir "../../" ${CMAKE_COMMAND} --build . --target library
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building library"
)

# verycleanlib target
add_custom_target(verycleanlib
    COMMAND ${CMAKE_COMMAND} -E chdir "../../" ${CMAKE_COMMAND} --build . --target veryclean
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Cleaning library"
)

# Install target (optional)
install(TARGETS
    amg_sf_sample
    amg_df_sample
    amg_cf_sample
    amg_zf_sample
    DESTINATION bin
)
