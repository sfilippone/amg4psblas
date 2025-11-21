cmake_minimum_required(VERSION 3.15)
project(AMGExamples Fortran CXX)

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



find_package(MPI REQUIRED Fortran CXX )

if(MPI_FOUND)


  if( (MPI_CXX_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
    message ( WARNING
      "The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is
known to cause segmentation faults for some Fortran codes. See, e.g.,
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729 or
https://github.com/sourceryinstitute/OpenCoarrays/issues/317.

`noexecstack` is being replaced with `execstack`"
      )
    string(REPLACE "noexecstack"
      "execstack" MPI_CXX_LINK_FLAGS_FIXED ${MPI_CXX_LINK_FLAGS})
    string(REPLACE "noexecstack"
      "execstack" MPI_C_LINK_FLAGS_FIXED ${MPI_C_LINK_FLAGS})
    string(REPLACE "noexecstack"
      "execstack" MPI_Fortran_LINK_FLAGS_FIXED ${MPI_Fortran_LINK_FLAGS})
    set(MPI_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS_FIXED}" CACHE STRING
      "MPI CXX linking flags" FORCE)
    set(MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS_FIXED}" CACHE STRING
      "MPI C linking flags" FORCE)
    set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS_FIXED}" CACHE STRING
      "MPI Fortran linking flags" FORCE)
      
  endif()

    message(STATUS "Found MPI: ${MPI_C_LIBRARIES} - ${MPI_CXX_LIBRARIES} - ${MPI_Fortran_LIBRARIES}")
    
    
  #----------------
  # Setup MPI compiler
  #----------------  
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C++ compiler" FORCE)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE FILEPATH "Fortran compiler" FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Add -g to the Fortran compiler flags.
    # We use STRING(APPEND) to ensure we don't overwrite other important flags.
    string(APPEND CMAKE_Fortran_FLAGS " -g")
    message(STATUS "Fortran debug flags added: -g")
endif()

# Include directories
include_directories(
    "."
    "${AMG4PSBLAS_INSTALL_DIR}/include"
    "${AMG4PSBLAS_INSTALL_DIR}/modules"
    "${PSBLAS_INSTALL_DIR}/include"
    "${PSBLAS_INSTALL_DIR}/modules"
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
set(DFSOBJS ${COMMON_SOURCE} amg_df_sample.f90)
set(SFSOBJS ${COMMON_SOURCE} amg_sf_sample.f90)
set(CFSOBJS ${COMMON_SOURCE} amg_cf_sample.f90)
set(ZFSOBJS ${COMMON_SOURCE} amg_zf_sample.f90)


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
    ${COMMON_SOURCE} 
    amg_d_pde3d.F90
    amg_d_genpde_mod.F90
    ${DGEN3D}
)
set(amg_s_pde3d_SOURCES
    ${COMMON_SOURCE} 
    amg_s_pde3d.F90
    amg_s_genpde_mod.F90
    ${SGEN3D}
)
set(amg_d_pde2d_SOURCES
    ${COMMON_SOURCE} 
    amg_d_pde2d.F90
    amg_d_genpde_mod.F90
    ${DGEN2D}
)
set(amg_s_pde2d_SOURCES
    ${COMMON_SOURCE} 
    amg_s_pde2d.F90
    amg_s_genpde_mod.F90
    ${SGEN2D}
)




# Function to create executable
function(create_amg_executable target sources)
    add_executable(${target} ${sources})
    target_link_libraries(${target} PUBLIC
        ${AMG_LIBS}
        ${PSBLAS_LIBS}
        ${MPI_LIBRARIES}
    )

    # Move executable to EXEDIR (post-build)
    #add_custom_command(TARGET ${target} POST_BUILD
    #    COMMAND ${CMAKE_COMMAND} -E move $<TARGET_FILE:${target}> ${EXEDIR}
    #)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
        LINKER_LANGUAGE Fortran
        )
endfunction()


# Create executables
create_amg_executable(amg_d_pde3d "${amg_d_pde3d_SOURCES}")
create_amg_executable(amg_s_pde3d "${amg_s_pde3d_SOURCES}")
create_amg_executable(amg_d_pde2d "${amg_d_pde2d_SOURCES}")
create_amg_executable(amg_s_pde2d "${amg_s_pde2d_SOURCES}")

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
