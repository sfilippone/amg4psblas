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
    "${PSBLAS_INSTALL_DIR}/include" 
)

# Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(FMFLAG "-I")

# Library directories
link_directories(
    "${AMG_INSTALL_DIR}/lib"
    "${PSBLAS_INSTALL_DIR}/lib"
)

# Libraries
set(AMG_LIBS psb_linsolve amg_prec psb_prec)
set(PSBLAS_LIBS psb_util psb_linsolve psb_prec psb_base) 

# Executable directory
set(EXEDIR "${CMAKE_CURRENT_BINARY_DIR}/runs")
file(MAKE_DIRECTORY "${EXEDIR}")

# Source files
set(DGEN2D
    amg_d_pde2d_poisson_mod.f90
    amg_d_pde2d_exp_mod.f90
    amg_d_pde2d_gauss_mod.f90
    amg_d_pde2d_box_mod.f90
)
set(DGEN3D
    amg_d_pde3d_poisson_mod.f90
    amg_d_pde3d_exp_mod.f90
    amg_d_pde3d_gauss_mod.f90
    amg_d_pde3d_box_mod.f90
)
set(SGEN2D
    amg_s_pde2d_poisson_mod.f90
    amg_s_pde2d_exp_mod.f90
    amg_s_pde2d_gauss_mod.f90
    amg_s_pde2d_box_mod.f90
)
set(SGEN3D
    amg_s_pde3d_poisson_mod.f90
    amg_s_pde3d_exp_mod.f90
    amg_s_pde3d_gauss_mod.f90
    amg_s_pde3d_box_mod.f90
)

# Define executables and their sources
set(amg_d_pde3d_SOURCES
    amg_d_pde3d.f90
    amg_d_genpde_mod.f90
    ${DGEN3D}
    data_input.f90
)
set(amg_s_pde3d_SOURCES
    amg_s_pde3d.f90
    amg_s_genpde_mod.f90
    ${SGEN3D}
    data_input.f90
)
set(amg_d_pde2d_SOURCES
    amg_d_pde2d.f90
    amg_d_genpde_mod.f90
    ${DGEN2D}
    data_input.f90
)
set(amg_s_pde2d_SOURCES
    amg_s_pde2d.f90
    amg_s_genpde_mod.f90
    ${SGEN2D}
    data_input.f90
)

# Function to create executable
macro(create_amg_executable target sources)
    add_executable(${target} ${sources})
    target_link_libraries(${target} PUBLIC
        ${AMG_LIBS}
        ${PSBLAS_LIBS}

    )

    # Move executable to EXEDIR (post-build)
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E move $<TARGET_FILE:${target}> ${EXEDIR}
    )
endmacro()

# Create executables
create_amg_executable(amg_d_pde3d ${amg_d_pde3d_SOURCES})
create_amg_executable(amg_s_pde3d ${amg_s_pde3d_SOURCES})
create_amg_executable(amg_d_pde2d ${amg_d_pde2d_SOURCES})
create_amg_executable(amg_s_pde2d ${amg_s_pde2d_SOURCES})

# Create "runs" directory (if it doesn't exist)
add_custom_target(create_runs_dir
    COMMAND ${CMAKE_COMMAND} -E make_directory "${EXEDIR}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating runs directory"
)
add_dependencies(amg_d_pde3d create_runs_dir)
add_dependencies(amg_s_pde3d create_runs_dir)
add_dependencies(amg_d_pde2d create_runs_dir)
add_dependencies(amg_s_pde2d create_runs_dir)


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

# Check target (simulated)
add_custom_target(check
    COMMAND ${CMAKE_COMMAND} -E chdir "${EXEDIR}" "${EXEDIR}/amg_d_pde2d" < "amg_pde2d.inp"
    COMMAND ${CMAKE_COMMAND} -E chdir "${EXEDIR}" "${EXEDIR}/amg_s_pde2d" < "amg_pde2d.inp"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS amg_d_pde2d amg_s_pde2d
    COMMENT "Running check"
)

# Install target (optional)
install(TARGETS
    amg_d_pde3d
    amg_s_pde3d
    amg_d_pde2d
    amg_s_pde2d
    DESTINATION bin
)
