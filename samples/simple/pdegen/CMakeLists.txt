cmake_minimum_required(VERSION 3.15)
project(AMGExamples Fortran)

# Installation directories (passed as CMake variables)
set(AMG4PSBLAS_INSTALL_DIR "" CACHE PATH "Path to AMG installation")
set(PSBLAS_INSTALL_DIR "" CACHE PATH "Path to PSBLAS installation")

# Check if installation directories are set
if(NOT AMG4PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "AMG_INSTALL_DIR must be set.  Use -DAMG_INSTALL_DIR=/path/to/amg when running CMake.")
endif()

if(NOT PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "PSBLAS_INSTALL_DIR must be set.  Use -DPSBLAS_INSTALL_DIR=/path/to/psblas when running CMake.")
endif()


find_package(psblas REQUIRED PATHS ${PSBLAS_INSTALL_DIR})
find_package(amg4psblas REQUIRED PATHS ${AMG4PSBLAS_INSTALL_DIR})


# Include directories
include_directories(
    "."
    "${AMG4PSBLAS_INSTALL_DIR}/include"
    "${AMG4PSBLAS_INSTALL_DIR}/modules"
    "${PSBLAS_INSTALL_DIR}/include" 
)

# Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(FMFLAG "-I") 

# Library directories
link_directories(
    "${AMG4PSBLAS_INSTALL_DIR}/lib"
    "${PSBLAS_INSTALL_DIR}/lib"
)

# Libraries
set(AMG_LIBS amg4psblas::amgprec)
set(PSBLAS_LIBS  psblas::util psblas::linsolve psblas::prec psblas::base) 

# Executable directory
set(EXEDIR "${CMAKE_CURRENT_BINARY_DIR}/runs")
file(MAKE_DIRECTORY "${EXEDIR}")

set(COMMON_SOURCE data_input.f90)


# Source files
set(DMOBJS ${COMMON_SOURCE} amg_dexample_ml.f90 amg_dpde_mod.f90)
set(D1OBJS ${COMMON_SOURCE} amg_dexample_1lev.f90 amg_dpde_mod.f90)
set(SMOBJS ${COMMON_SOURCE} amg_sexample_ml.f90 amg_spde_mod.f90)
set(S1OBJS ${COMMON_SOURCE} amg_sexample_1lev.f90 amg_spde_mod.f90)


# Function to create executable
macro(create_amg_executable target sources)
    add_executable(${target} ${sources})
    target_link_libraries(${target} PUBLIC
        ${AMG_LIBS}
        ${PSBLAS_LIBS}

    )

    # Move executable to EXEDIR (post-build)
    #add_custom_command(TARGET ${target} POST_BUILD
    #    COMMAND ${CMAKE_COMMAND} -E move $<TARGET_FILE:${target}> ${EXEDIR}
    #)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
        )
endmacro()
# Create executables
create_amg_executable(amg_dexample_ml "${DMOBJS}")
create_amg_executable(amg_dexample_1lev "${D1OBJS}")
create_amg_executable(amg_sexample_ml "${SMOBJS}")
create_amg_executable(amg_sexample_1lev "${S1OBJS}")

# Create "runs" directory (if it doesn't exist)
add_custom_target(create_runs_dir
    COMMAND ${CMAKE_COMMAND} -E make_directory "${EXEDIR}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating runs directory"
)
add_dependencies(amg_dexample_ml create_runs_dir)
add_dependencies(amg_dexample_1lev create_runs_dir)
add_dependencies(amg_sexample_ml create_runs_dir)
add_dependencies(amg_sexample_1lev create_runs_dir)


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
    amg_dexample_ml
    amg_dexample_1lev
    amg_sexample_ml
    amg_sexample_1lev
    DESTINATION bin
)
